using System.Diagnostics;
using System.Windows.Input;
using Avalonia;
using Avalonia.Controls;
using Avalonia.Input;
using Avalonia.Interactivity;
using Avalonia.Media.Imaging;
using Avalonia.Platform;
using Avalonia.Platform.Storage;
using Avalonia.Threading;
using Avalonia.VisualTree;
using CapsStudio.ViewModels;

namespace CapsStudio.Views;

public sealed class RelayCommand(Func<Task> run) : ICommand
{
    public event EventHandler? CanExecuteChanged { add { } remove { } }
    public bool CanExecute(object? parameter) => true;
    public async void Execute(object? parameter) => await run();
}

public partial class MainWindow : Window
{
    private readonly MainViewModel _vm = new();
    private readonly string? _samples = FindSamples();

    // Render loop state: latest request wins; one render in flight at a time.
    private int _requested, _rendered;
    private bool _busy;
    private int _pixW, _pixH;
    private double _scaling = 1;

    private Point _press, _last;
    private Control? _hostField;
    private Control _host => _hostField ?? ViewHost;
    private bool _dragging, _moved;
    // lasso and move tools: the lasso's points, the atoms a move carries and the projected atoms that fix the view plane
    private List<Point>? _lassoPts;
    private int[]? _moveAtoms;
    private List<(double X, double Y, double Z, double Sx, double Sy)>? _moveFit;
    private Interop.CapsCamera _lastCam;
    private Interop.CapsRenderOpts _lastOpt;
    private bool _haveLast;
    private bool _pan, _addPick;

    private readonly DispatcherTimer _playTimer = new() { Interval = TimeSpan.FromMilliseconds(125) };

    public MainViewModel ViewModel => _vm;
    public void SelectAnalysisTab(int k) => AnalysisTabs.SelectedIndex = k;
    public ICommand OpenCommand { get; }
    public ICommand SaveCommand { get; }
    public bool HasSamples => _samples != null;

    private static bool _namesInstalled;

    public MainWindow()
    {
        if (!_namesInstalled) { AccessibleNames.Install(); _namesInstalled = true; }
        OpenCommand = new RelayCommand(OpenDialog);
        InitializeComponent();
        DataContext = _vm;
        // Tab order follows the accessibility map: rail, tools, project, 3D view, inspector, dock
        var order = new Control?[] { Rail, Toolbar, ProjectPanel, ViewHost, InspectorPanel, AnalysisDock };
        for (var k = 0; k < order.Length; ++k) if (order[k] is { } c) KeyboardNavigation.SetTabIndex(c, k + 1);
        _vm.ScaleChanged += k => ScaleRoot.LayoutTransform = Math.Abs(k - 1) < 1e-9 ? null : new Avalonia.Media.ScaleTransform(k, k);
        _vm.LoadSettings();
        _vm.HookJobs();
        _vm.LoadRecent();
        _vm.LoadLastSession();
        Closing += (_, _) => _vm.SaveSession();   // the project tree comes back from Start next time
        AddWindowCommands();
        _vm.InitProtocol();
        _vm.LoadReactionSet();
        if (Paths.Python is { } py) CapsStudio.Interop.Native.SetPython(py, Environment.GetEnvironmentVariable("CAPS_PYTHON"));
        KeyBindings.Add(new KeyBinding { Gesture = new KeyGesture(Key.O, KeyModifiers.Meta), Command = OpenCommand });
        KeyBindings.Add(new KeyBinding { Gesture = new KeyGesture(Key.O, KeyModifiers.Control), Command = OpenCommand });
        var recentCommand = new RelayCommand(() => { if (_vm.Idle) OpenMostRecent(); return Task.CompletedTask; });
        KeyBindings.Add(new KeyBinding { Gesture = new KeyGesture(Key.O, KeyModifiers.Meta | KeyModifiers.Shift), Command = recentCommand });
        KeyBindings.Add(new KeyBinding { Gesture = new KeyGesture(Key.O, KeyModifiers.Control | KeyModifiers.Shift), Command = recentCommand });
        SaveCommand = new RelayCommand(() => _vm.HasDocument && _vm.Idle ? SaveAs("data", "LAMMPS data") : Task.CompletedTask);
        KeyBindings.Add(new KeyBinding { Gesture = new KeyGesture(Key.S, KeyModifiers.Meta), Command = SaveCommand });
        KeyBindings.Add(new KeyBinding { Gesture = new KeyGesture(Key.S, KeyModifiers.Control), Command = SaveCommand });
        _vm.RenderRequested += RequestRender;
        _vm.CompareStatesChanged += () =>
        {
            var sh = _vm.CompareShifts;
            this.FindControl<LinePlot>("ShiftPlot")?.SetData(sh.Select((y, i) => ((double)(i + 1), y)).ToArray());
        };
        _vm.ViewRequested += RequestViewRender;
        _vm.PropertyChanged += (_, e) =>
        {
            if (e.PropertyName == nameof(ViewModels.MainViewModel.RunLiveDoc) && _vm.RunLiveDoc != null && _vm.Busy) RequestRender();
            else if (e.PropertyName == nameof(ViewModels.MainViewModel.Busy)) RequestRender();
        };
        ViewGl.FitStale += RequestViewRender;
        ViewGl.ReadyChanged += () => Dispatcher.UIThread.Post(() =>
        {
            _vm.GpuStatus = ViewGl.Ready ? "In use · " + ViewGl.Status : "Not in use: " + ViewGl.Status;
            if (!ViewGl.Ready) { ViewGl.IsVisible = false; ViewImage.IsVisible = true; }
            RequestRender();
        });
        RenderGuide.Vm = _vm;
        PipeTablePlot.Brushable = true;
        PipeTablePlot.Brushed += (x0, x1, y0, y1) => _vm.ApplyBrush(x0, x1, y0, y1);
        _vm.FullViewRequested += live => OpenFullView(live);
        _vm.PipeTableChanged += () =>
        {
            PipeTablePlot.XLabel = _vm.PipeTableXLabel;
            PipeTablePlot.YLabel = _vm.PipeTableYLabel;
            if (_vm.PipeTableHeat is { } heat) { PipeTablePlot.ZLabel = _vm.PipeTableZLabel; PipeTablePlot.SetHeat(heat.X, heat.Y, heat.Z); return; }
            var pts = _vm.PipeTableX.Zip(_vm.PipeTableY).Where(p => double.IsFinite(p.Second)).ToArray();
            if (_vm.PipeTableScatter) { PipeTablePlot.RefY = null; PipeTablePlot.SetData(pts, []); }
            else { PipeTablePlot.RefY = 1.0; PipeTablePlot.Markers = false; PipeTablePlot.SetData(pts); }
        };
        _vm.RenderOverlayChanged += () => RenderGuide.InvalidateVisual();
        _vm.PropertyChanged += (_, e) =>
        {
            if (e.PropertyName == nameof(MainViewModel.HasFocusAtom) && _vm.HasFocusAtom && !_focusTabShown)
            {
                _focusTabShown = true;   // the first walk shows the neighbour table; later the user's tab choice stays
                AnalysisTabs.SelectedItem = FocusTab;
            }
            if (e.PropertyName == nameof(MainViewModel.ViewIsLight))
            {
                ViewHost.Background = _vm.ViewIsLight ? Avalonia.Media.Brushes.White : (Avalonia.Media.IBrush)this.FindResource("Bg0B")!;
                HintText.Foreground = _vm.ViewIsLight ? new Avalonia.Media.SolidColorBrush(Avalonia.Media.Color.Parse("#5A6168")) : (Avalonia.Media.IBrush)this.FindResource("DimB")!;
            }
            if (e.PropertyName == nameof(MainViewModel.RdfCurve)) RdfPlot.SetData(_vm.RdfCurve.Select(p => (p.R, p.G)).ToArray());
            if (e.PropertyName == nameof(MainViewModel.Relaxing) && _vm.Relaxing)
            {
                EnergyPlot.RefY = null;
                ForcePlot.RefY = Math.Log10((double)(_vm.RelaxFtolD ?? 0.5m));
                _vm.ShowRdf = true;
                AnalysisTabs.SelectedIndex = 1;
            }
            if (e.PropertyName == nameof(MainViewModel.MdRunning) && _vm.MdRunning)
            {
                TempPlot.RefY = (double)(_vm.MdTempD ?? 300m);
                PressPlot.RefY = _vm.MdHasBarostat ? (double)(_vm.MdPressureD ?? 1m) : null;
                DensPlot.RefY = null;
                _vm.ShowRdf = true;
                AnalysisTabs.SelectedIndex = 2;
            }
            if (e.PropertyName == nameof(MainViewModel.EqRunning) && _vm.EqRunning)
            {
                TempPlot.RefY = null;
                PressPlot.RefY = null;
                DensPlot.RefY = null;
                _vm.ShowRdf = true;
                AnalysisTabs.SelectedIndex = 2;
            }
            if (e.PropertyName == nameof(MainViewModel.ChainCurve))
            {
                if (_vm.RisCurve.Length > 0) ChainPlot.SetCompare(_vm.ChainCurve, _vm.RisCurve);
                else ChainPlot.SetData(_vm.ChainCurve);
            }
            if (e.PropertyName == nameof(MainViewModel.Reacting) && _vm.Reacting)
            {
                ConvPlot.RefY = null;
                GelPlot.RefY = null;
                _vm.ShowRdf = true;
                AnalysisTabs.SelectedIndex = 4;
            }
            if (e.PropertyName == nameof(MainViewModel.Busy) && !_vm.Busy) RequestRender();
        };
        // one 3D view: it moves into the page that shows a live view (Grow's live cell), and back to the Studio
        _viewHome = (Panel)ViewHost.Parent!;
        _vm.PropertyChanged += (_, e) => { if (e.PropertyName == nameof(MainViewModel.Module)) PlaceViewport(); };
        _vm.Analyze.PropertyChanged += (_, e) => { if (e.PropertyName == nameof(AnalyzeViewModel.Curve)) ShowCurve(_vm.Analyze.Curve); };
        _vm.ThermoChanged += () =>
        {
            var t = _vm.Thermo;
            TempPlot.SetData(t.Select(r => (r.TimePs, r.Temperature)).ToArray());
            PressPlot.SetData(t.Select(r => (r.TimePs, r.Pressure)).ToArray());
            DensPlot.SetData(t.Select(r => (r.TimePs, r.Density)).ToArray());
        };
        _vm.ReactChanged += () =>
        {
            var rows = _vm.ReactRows;
            ConvPlot.SetData(rows.Select(r => ((double)r.Cycle, r.Conversion)).ToArray());
            GelPlot.SetData(rows.Select(r => (r.Conversion, r.LargestFraction)).ToArray());
        };
        _vm.RelaxCurvesChanged += () =>
        {
            EnergyPlot.SetData(_vm.RelaxEnergyCurve);
            ForcePlot.SetData(_vm.RelaxForceCurve);
        };

        foreach (var host in new[] { ViewHost, FieldViewHost })
        {
            host.SizeChanged += (_, _) => RequestRender();
            host.PointerPressed += OnPointerPressed;
            host.PointerMoved += OnPointerMoved;
            host.PointerReleased += OnPointerReleased;
            host.PointerWheelChanged += OnWheel;
        }
        _vm.Field.PropertyChanged += (_, e) => { if (e.PropertyName == nameof(FieldViewModel.SelectedRow)) RequestRender(); };
        _playTimer.Tick += (_, _) => { if (!_vm.PlayStep()) _playTimer.Stop(); };
        _vm.PlaybackChanged += () => _playTimer.Interval = TimeSpan.FromMilliseconds(1000.0 / _vm.PlaybackFps);
        Timeline.Vm = _vm;
        _vm.TimelineChanged += () => Timeline.InvalidateVisual();
        _vm.PropertyChanged += (_, e) => { if (e.PropertyName == nameof(MainViewModel.Frame)) Timeline.InvalidateVisual(); };
        AddHandler(KeyDownEvent, OnKey, RoutingStrategies.Tunnel);
        AddHandler(KeyUpEvent, (_, e) => { if (e.Key == Key.L) _vm.LensHold = false; }, RoutingStrategies.Tunnel);
        DragDrop.SetAllowDrop(ViewHost, true);
        ViewHost.AddHandler(DragDrop.DropEvent, OnDrop);
        NativeMenu.SetMenu(this, BuildMenu());   // the menu bar: the palette's commands and the pages
    }

    public void RunPackExampleForTest() { if (_samples != null) _vm.AddPackExample(_samples); }

    public void PickForTest(params int[] atoms)
    {
        for (var k = 0; k < atoms.Length; k++) _vm.Pick(atoms[k], k > 0);
        RequestRender();
    }

    public void OpenOnStart(string path, string? topology) => Dispatcher.UIThread.Post(() => TryOpen(path, topology));

    private static string? FindSamples() => Paths.Samples;

    private void TryOpen(string path, string? topology = null)
    {
        try
        {
            // A dump opened on its own picks up a data file with the same stem, for types, masses and bonds.
            if (topology == null && (path.EndsWith(".lammpstrj") || path.EndsWith(".dump")))
            {
                var data = Path.ChangeExtension(path, ".data");
                if (File.Exists(data)) topology = data;
            }
            _vm.Open(path, topology);
        }
        catch (Exception e)
        {
            _vm.Status = $"Could not open {Path.GetFileName(path)}: {e.Message}";
        }
    }

    /// <summary>File picker, then the Open preview page (format, columns, types, topology) before reading.</summary>
    public async Task OpenWithPreview()
    {
        var files = await StorageProvider.OpenFilePickerAsync(new FilePickerOpenOptions
        {
            Title = "Open a structure or trajectory", AllowMultiple = false,
            FileTypeFilter =
            [
                new FilePickerFileType("Structures and trajectories") { Patterns = ["*.data", "*.lmp", "*.lammpstrj", "*.dump", "*.dcd", "*.gro", "*.xtc", "*.trr", "*.prmtop", "*.parm7", "*.inpcrd", "*.rst7", "*.restrt", "*.ncrst", "*.nc", "*.mdcrd", "*.pdb", "*.ent", "*.xyz", "*.extxyz", "*.mol2", "*.sdf", "*.mol", "*.cif", "*.car", "*.vasp", "POSCAR*", "CONTCAR*", "*.gz"] },
                new FilePickerFileType("All files") { Patterns = ["*"] },
            ],
        });
        if (files.Count > 0 && files[0].TryGetLocalPath() is { } p) _vm.PreviewOpen(p);
    }

    public async Task OpenDialog()
    {
        var files = await StorageProvider.OpenFilePickerAsync(new FilePickerOpenOptions
        {
            Title = "Open structure or trajectory",
            AllowMultiple = true,
            FileTypeFilter =
            [
                new FilePickerFileType("Structures and trajectories") { Patterns = ["*.data", "*.lmp", "*.lammpstrj", "*.dump", "*.dcd", "*.gro", "*.xtc", "*.trr", "*.prmtop", "*.parm7", "*.inpcrd", "*.rst7", "*.restrt", "*.ncrst", "*.nc", "*.mdcrd", "*.pdb", "*.ent", "*.xyz", "*.extxyz", "*.mol2", "*.sdf", "*.mol", "*.cif", "*.car", "*.vasp", "POSCAR*", "CONTCAR*", "*.gz"] },
                new FilePickerFileType("All files") { Patterns = ["*"] },
            ],
        });
        OpenMany(files.Select(f => f.TryGetLocalPath()).OfType<string>().ToList());
    }

    public void OpenMany(List<string> paths)
    {
        if (paths.Count == 0) return;
        // Dropped together: a dump plus a data file pair up (dump is the trajectory, data is the topology).
        var dump = paths.FirstOrDefault(p => p.EndsWith(".lammpstrj") || p.EndsWith(".dump"));
        var data = paths.FirstOrDefault(p => p.EndsWith(".data") || p.EndsWith(".lmp"));
        if (dump != null) TryOpen(dump, data);
        else TryOpen(paths[0]);
    }

    private void OnDrop(object? sender, DragEventArgs e)
    {
        if (!_vm.Idle) return;
        var files = e.Data.GetFiles();
        if (files != null) OpenMany(files.Select(f => f.TryGetLocalPath()).OfType<string>().ToList());
    }

    private bool _focusTabShown;

    public void OnOpenSample(object? sender, RoutedEventArgs e) => OpenSample("ps");

    /// <summary>A sample made by CAPS: "ps" (dump with its data file), "gro" (the same cell as .gro) or "water" (PDB).</summary>
    public void OpenSample(string which)
    {
        if (_samples == null) return;
        if (which == "gro") TryOpen(Path.Combine(_samples, "ps_melt.gro"));
        else if (which == "water") TryOpen(Path.Combine(_samples, "water.pdb"));
        else TryOpen(Path.Combine(_samples, "ps_melt.lammpstrj"), Path.Combine(_samples, "ps_melt.data"));
    }

    /// <summary>⌘⇧O: the most recent file, with its topology.</summary>
    public void OpenMostRecent()
    {
        var r = RecentFiles.Load().FirstOrDefault(x => File.Exists(x.Path));
        if (r == null) { _vm.Status = "No recent file to open"; return; }
        OpenMany(r.Topology != null && File.Exists(r.Topology) ? [r.Path, r.Topology] : [r.Path]);
    }

    private void OnPlay(object? s, RoutedEventArgs e) => TogglePlay();
    private void OnTrajectoryPlayer(object? s, RoutedEventArgs e) => _vm.OpenTrajectory();
    private void OnTorsionScan(object? s, RoutedEventArgs e) => _vm.OpenTorsion();
    private void OnSplit(object? s, RoutedEventArgs e) => _vm.OpenSplit();
    private void OnFragmentLibrary(object? s, RoutedEventArgs e) => _vm.OpenFragments();
    private async void OnQuickFragment(object? s, RoutedEventArgs e) { if ((s as Control)?.Tag is ViewModels.FragmentItem f) await _vm.UseFragment(f, false); }
    private void OnFragmentSearchKey(object? s, KeyEventArgs e) { if (e.Key == Key.Enter) _vm.OpenFragments(_vm.FragmentQuery); }
    private void OnLoadCancel(object? s, RoutedEventArgs e) => _vm.CancelLoad();
    private void OnLoadBackground(object? s, RoutedEventArgs e) => _vm.LoadToBackground();

    private void OnModuleGrow(object? s, RoutedEventArgs e) => _vm.SetModule(0);
    private void OnModuleAnalyze(object? s, RoutedEventArgs e) => _vm.SetModule(1);
    private void OnModuleRelax(object? s, RoutedEventArgs e) => _vm.SetModule(2);
    private void OnCheckField(object? s, RoutedEventArgs e) => _vm.CheckField();
    private async void OnRelax(object? s, RoutedEventArgs e) => await _vm.Relax();
    private void OnRelaxCancel(object? s, RoutedEventArgs e) => _vm.CancelRelax();
    private void OnModuleDynamics(object? s, RoutedEventArgs e) => _vm.SetModule(3);
    private async void OnMdRun(object? s, RoutedEventArgs e) => await _vm.RunMd();
    private void OnMdCancel(object? s, RoutedEventArgs e) => _vm.CancelMd();
    private void OnModuleEquilibrate(object? s, RoutedEventArgs e) => _vm.SetModule(4);
    private void OnModulePack(object? s, RoutedEventArgs e) => _vm.SetModule(5);
    private void OnModuleReact(object? s, RoutedEventArgs e) => _vm.SetModule(6);
    private void OnModuleField(object? s, RoutedEventArgs e) => _vm.SetModule(7);
    private void OnModuleStudio(object? s, RoutedEventArgs e) => _vm.SetModule(8);
    private void OnModuleBuild(object? s, RoutedEventArgs e) => _vm.SetModule(13);
    private void OnExportCenter(object? s, RoutedEventArgs e) => _vm.OpenExportCenter();
    private void OnPipelineNext(object? s, RoutedEventArgs e) => _vm.GoPipelineNext();
    private void OnPipelineStep(object? s, RoutedEventArgs e) { if ((s as Control)?.DataContext is ViewModels.PipelineStep st) _vm.GoPipelineStep(st); }

    private readonly Panel _viewHome;
    private void PlaceViewport()
    {
        Decorator? slot = _vm.IsGrow ? GrowPageView.Slot : _vm.IsPack ? PackPageView.Slot : null;
        if (slot != null && ViewHost.Parent != slot)
        {
            _viewHome.Children.Remove(ViewHost);
            slot.Child = ViewHost;
        }
        else if (slot == null && ViewHost.Parent is Decorator d)
        {
            d.Child = null;
            _viewHome.Children.Add(ViewHost);
        }
        RequestRender();
    }
    // Studio toolbar
    private void OnToolSelect(object? s, RoutedEventArgs e) => _vm.MeasureTool = false;
    private void OnViewBgDark(object? s, RoutedEventArgs e) => _vm.ViewBackground = 0;
    private void OnViewBgLight(object? s, RoutedEventArgs e) => _vm.ViewBackground = 1;
    private void OnToolMeasure(object? s, RoutedEventArgs e) => _vm.MeasureTool = true;
    private void OnStyleItem(object? s, RoutedEventArgs e) { if (s is MenuItem { Tag: string t }) _vm.StyleIndex = int.Parse(t); }
    private void OnColourItem(object? s, RoutedEventArgs e) { if (s is MenuItem { Tag: string t }) _vm.ColourIndex = int.Parse(t); }
    private void OnPerspectiveOn(object? s, RoutedEventArgs e) => _vm.Perspective = true;
    private void OnPerspectiveOff(object? s, RoutedEventArgs e) => _vm.Perspective = false;
    public void ShowSettings() => _vm.SetModule(10);
    private void OnSettingsRail(object? s, RoutedEventArgs e) => _vm.SetModule(10);
    private void OnModuleJobs(object? s, RoutedEventArgs e) => _vm.SetModule(11);
    private void OnModuleBench(object? s, RoutedEventArgs e) => _vm.SetModule(12);
    private void OnCloseDocument(object? s, RoutedEventArgs e) { e.Handled = true; _vm.CloseDocument(); }
    private void OnProjectTab(object? s, RoutedEventArgs e) { if ((s as Control)?.Tag is ViewModels.ProjectItem it) _vm.Activate(it); }
    private void OnCloseProjectTab(object? s, RoutedEventArgs e)
    {
        e.Handled = true;
        if ((s as Control)?.Tag is not ViewModels.ProjectItem it) return;
        _vm.Activate(it);
        if (_vm.ActiveItem == it) _vm.CloseDocument();
    }
    private void OnDuplicateStructure(object? s, RoutedEventArgs e) => _vm.DuplicateStructure();
    private void OnThemeDark(object? s, RoutedEventArgs e) { Tokens.Use(false); RequestRender(); }
    private void OnThemeLight(object? s, RoutedEventArgs e) { Tokens.Use(true); RequestRender(); }
    private void OnCommandPalette(object? s, RoutedEventArgs e) => TogglePalette();

    public void TogglePalette()
    {
        _vm.PaletteOpen = !_vm.PaletteOpen;
        if (_vm.PaletteOpen) Dispatcher.UIThread.Post(() => PaletteBox.Focus(), DispatcherPriority.Input);
    }

    private void OnPaletteBackdrop(object? s, PointerPressedEventArgs e) => _vm.PaletteOpen = false;

    private void OnPaletteKey(object? s, KeyEventArgs e)
    {
        switch (e.Key)
        {
            case Key.Down: _vm.PaletteMove(1); PaletteList.ScrollIntoView(_vm.PaletteIndex); e.Handled = true; break;
            case Key.Up: _vm.PaletteMove(-1); PaletteList.ScrollIntoView(_vm.PaletteIndex); e.Handled = true; break;
            case Key.Enter: _vm.PaletteRun(); e.Handled = true; break;
            case Key.Escape: _vm.PaletteOpen = false; e.Handled = true; break;
        }
    }

    private void OnPaletteTapped(object? s, TappedEventArgs e)
    {
        if ((e.Source as Control)?.DataContext is PaletteRow r && r.IsCommand) _vm.PaletteRun(r);
    }

    /// <summary>Palette commands that need the window (file dialogs, the view's size).</summary>
    private void AddWindowCommands()
    {
        _vm.AddCommand(new PaletteCommand { Title = "Open with preview…", Id = "document.open preview", Icon = "folder", Section = "File",
            Keywords = "inspect columns types topology format first lines", Run = () => _ = OpenWithPreview() });
        _vm.AddCommand(new PaletteCommand { Title = "Open a structure or trajectory…", Id = "document.open", Icon = "folder", Shortcut = "⌘O", Section = "File",
            Keywords = "load file lammps gromacs pdb xyz mol2", Run = () => _ = OpenDialog() });
        _vm.AddCommand(new PaletteCommand { Title = "Save as LAMMPS data…", Id = "document.save data", Icon = "save", Shortcut = "⌘S", Section = "File",
            Enabled = () => _vm.HasDocument && _vm.Idle, Run = () => _ = SaveAs("data", "LAMMPS data") });
        _vm.AddCommand(new PaletteCommand { Title = "Save as PDB…", Id = "document.save pdb", Icon = "save", Section = "File",
            Enabled = () => _vm.HasDocument && _vm.Idle, Run = () => _ = SaveAs("pdb", "PDB") });
        _vm.AddCommand(new PaletteCommand { Title = "Save as mol2…", Id = "document.save mol2", Icon = "save", Section = "File",
            Enabled = () => _vm.HasDocument && _vm.Idle, Run = () => _ = SaveAs("mol2", "Tripos mol2") });
        _vm.AddCommand(new PaletteCommand { Title = "Save as Materials Studio (.car / .mdf)…", Id = "document.save car", Icon = "save", Section = "File",
            Keywords = "materials studio biovia discover car mdf msi2lmp interface",
            Enabled = () => _vm.HasDocument && _vm.Idle, Run = () => _ = SaveAs("car", "Materials Studio .car (with its .mdf)") });
        _vm.AddCommand(new PaletteCommand { Title = "Save the trajectory (LAMMPS dump)…", Id = "trajectory.save", Icon = "save", Section = "File",
            Enabled = () => _vm.HasDocument && _vm.HasFrames && _vm.Idle, Run = () => _ = SaveTrajectoryAsync() });
        _vm.AddCommand(new PaletteCommand { Title = "Export figure (PNG)…", Id = "export.png", Icon = "download", Section = "File", Keywords = "image picture render",
            Enabled = () => _vm.HasDocument, Run = () => _ = Export("png") });
        _vm.AddCommand(new PaletteCommand { Title = "Export figure (SVG)…", Id = "export.svg", Icon = "download", Section = "File", Keywords = "vector image",
            Enabled = () => _vm.HasDocument, Run = () => _ = Export("svg") });
        _vm.CompactChanged += ApplyCompact;
        SizeChanged += (_, e) => { _vm.Compact = e.NewSize.Width < 1440; ToolbarRight.Classes.Set("narrow", e.NewSize.Width < 1700); };
        // a folded dock opens when one of its tabs is chosen
        AnalysisTabs.AddHandler(PointerReleasedEvent, (_, _) => { if (_vm.Compact && !_vm.DockOpen) _vm.DockOpen = true; }, RoutingStrategies.Tunnel | RoutingStrategies.Bubble, true);
        _vm.TourChanged += UpdateTour;
        _vm.PropertyChanged += (_, e) => { if (e.PropertyName == nameof(MainViewModel.Document)) Dispatcher.UIThread.Post(_vm.MaybeStartTour, DispatcherPriority.Background); };
        _vm.OpenRequested += what =>
        {
            var parts = what.Split('\n');
            OpenMany(parts.Where(File.Exists).ToList());
        };
    }

    // ---- Analyze › Properties
    private async void OnAnalyzeRun(object? s, RoutedEventArgs e) { if (_vm.Idle) await _vm.Analyze.Run(); }
    private void OnAnalyzeCancel(object? s, RoutedEventArgs e) => _vm.Analyze.Cancel();

    private void OnAnalyzeExport(object? s, RoutedEventArgs e) => ExportAnalysis();

    /// <summary>Asks for a folder and writes the Analyze results there (also from the focused Analyze pages).</summary>
    public async void ExportAnalysis()
    {
        var dirs = await StorageProvider.OpenFolderPickerAsync(new FolderPickerOpenOptions { Title = "Folder for results.csv, results.tex and the curves", AllowMultiple = false });
        var dir = dirs.Count > 0 ? dirs[0].TryGetLocalPath() : null;
        if (dir == null) return;
        try
        {
            var n = _vm.Analyze.Export(dir);
            _vm.Status = $"Wrote {n} files to {dir}: results.csv, results.tex, results.json and one CSV per curve";
        }
        catch (Exception ex) { _vm.Status = "Could not export: " + ex.Message; }
    }

    private void ShowCurve(SeriesItem? c)
    {
        if (c == null) { PropPlot.SetData([]); return; }
        PropPlot.Markers = false;
        var pts = Enumerable.Range(0, c.X.Length).Select(k => (X: c.X[k], Y: c.Y[k])).Where(p => double.IsFinite(p.X) && double.IsFinite(p.Y));
        if (c.LogLog)
        {
            pts = pts.Where(p => p.X > 0 && p.Y > 0).Select(p => (Math.Log10(p.X), Math.Log10(p.Y)));
            PropPlot.XLabel = "log₁₀ " + c.XLabel;
            PropPlot.YLabel = "log₁₀ " + c.YLabel;
            PropPlot.RefY = null;
        }
        else
        {
            PropPlot.XLabel = c.XLabel;
            PropPlot.YLabel = c.YLabel;
            PropPlot.RefY = c.RefY;
        }
        if (c.Markers && c.OverlayX != null && c.OverlayY != null)
            PropPlot.SetData(pts.ToArray(), Enumerable.Range(0, c.OverlayX.Length).Select(k => (c.OverlayX[k], c.OverlayY[k])).ToArray());
        else
            PropPlot.SetData(pts.ToArray());
    }

    // ---------------------------------------------------------------- Field
    private async void OnFieldAssign(object? s, RoutedEventArgs e) { if (!_vm.Busy) await _vm.Field.Assign(); }
    private async void OnFieldClear(object? s, RoutedEventArgs e) { if (!_vm.Busy) await _vm.Field.Clear(); }
    private async void OnFieldUseFile(object? s, RoutedEventArgs e) { if (!_vm.Busy) await _vm.Field.UseFileForceField(); }
    private void OnApplyCoordination(object? s, RoutedEventArgs e) => _vm.ApplyCoordination();
    private async void OnFieldAssignGroups(object? s, RoutedEventArgs e) { if (!_vm.Busy) await _vm.Field.AssignGroups(); }
    private void OnFieldGroupAdd(object? s, RoutedEventArgs e) => _vm.Field.AddGroup();
    private void OnFieldGroupSuggest(object? s, RoutedEventArgs e) => _vm.Field.SuggestGroups();
    private void OnFieldGroupRemove(object? s, RoutedEventArgs e) { if ((s as Control)?.Tag is FieldGroupRow g) _vm.Field.RemoveGroup(g); }

    /// <summary>MEAM's alloy parameter file for a group.</summary>
    private async void OnFieldGroupBrowse2(object? s, RoutedEventArgs e)
    {
        if ((s as Control)?.Tag is not FieldGroupRow g) return;
        var files = await StorageProvider.OpenFilePickerAsync(new FilePickerOpenOptions
        {
            Title = "The MEAM parameter file", AllowMultiple = false,
            FileTypeFilter = [new FilePickerFileType("MEAM parameter files") { Patterns = ["*.meam", "*"] }],
        });
        if (files.Count > 0 && files[0].TryGetLocalPath() is { } p) g.File2 = p;
    }

    /// <summary>A literature potential's file for a group (Tersoff, EAM …).</summary>
    private async void OnFieldGroupBrowse(object? s, RoutedEventArgs e)
    {
        if ((s as Control)?.Tag is not FieldGroupRow g) return;
        var files = await StorageProvider.OpenFilePickerAsync(new FilePickerOpenOptions
        {
            Title = "The potential file (" + FieldViewModel.PotentialStyles[Math.Clamp(g.Style, 0, FieldViewModel.PotentialStyles.Length - 1)] + ")", AllowMultiple = false,
            FileTypeFilter =
            [
                new FilePickerFileType("LAMMPS potential files") { Patterns = ["*.tersoff", "*.tersoff.*", "*.sw", "*.vashishta", "*.gw", "*.eam.alloy", "*.eam.fs", "*.setfl", "*.airebo", "*.airebo-m", "*.rebo", "*.meam"] },
                new FilePickerFileType("All files") { Patterns = ["*"] },
            ],
        });
        if (files.Count > 0 && files[0].TryGetLocalPath() is { } p)
        {
            g.File = p;
            // the style from the file's name when it says (Si.tersoff, CuNi.eam.alloy)
            var n = Path.GetFileName(p).ToLowerInvariant();
            var k = n == "library.meam" || n.StartsWith("library") && n.EndsWith(".meam") ? "meam" : n.EndsWith(".airebo") ? "airebo" : n.EndsWith(".airebo-m") ? "airebo/morse" : n.EndsWith(".rebo") ? "rebo" : n.EndsWith(".eam.alloy") ? "eam/alloy" : n.EndsWith(".eam.fs") ? "eam/fs" : n.EndsWith(".sw") ? "sw" : n.EndsWith(".vashishta") ? "vashishta"
                  : n.EndsWith(".gw") ? "gw" : n.Contains(".tersoff") ? (n.Contains("zbl") ? "tersoff/zbl" : n.Contains(".mod.c") ? "tersoff/mod/c" : n.Contains(".mod") ? "tersoff/mod" : "tersoff") : "";
            var i = Array.IndexOf(FieldViewModel.PotentialStyles, k);
            if (i >= 0) g.Style = i;
        }
    }
    private async void OnFieldOverride(object? s, RoutedEventArgs e) { if (!_vm.Busy) await _vm.Field.ApplyOverride(); }
    private async void OnFieldResetOverride(object? s, RoutedEventArgs e) { if (!_vm.Busy) await _vm.Field.ResetOverride(); }
    private async void OnFieldAddRule(object? s, RoutedEventArgs e) { if (!_vm.Busy) await _vm.Field.AddRule(); }
    private async void OnFieldRemoveRules(object? s, RoutedEventArgs e) { if (!_vm.Busy) await _vm.Field.RemoveRules(); }

    private async void OnFieldCompanionCharges(object? s, RoutedEventArgs e)
    {
        if (!_vm.Busy) await _vm.Field.UseCompanionCharges();
    }

    private void OnSetMeasured(object? s, RoutedEventArgs e) => _vm.SetMeasured();
    private void OnFullView(object? s, RoutedEventArgs e) => OpenFullView(false);
    /// <summary>The structure (or the running job's live snapshots) alone in a full-screen window.</summary>
    public void OpenFullView(bool live)
    {
        if (_vm.Document == null && _vm.RunLiveDoc == null) return;
        new VisualizationWindow(_vm, live).Show(this);
    }

    private async void OnExportInspectorCsv(object? s, RoutedEventArgs e)
    {
        var name = _vm.InspectorTab switch { 0 => "particles.csv", 1 => "bonds.csv", 2 => "attributes.csv", _ => (_vm.PipeTableName ?? "table") + ".csv" };
        var file = await StorageProvider.SaveFilePickerAsync(new FilePickerSaveOptions
        {
            Title = "Export as CSV",
            SuggestedFileName = string.Concat(name.Select(c => char.IsLetterOrDigit(c) || c is '.' or '_' or '-' ? c : '_')),
            FileTypeChoices = [new FilePickerFileType("CSV") { Patterns = ["*.csv"] }],
        });
        if (file?.TryGetLocalPath() is { } path) _vm.ExportInspectorCsv(path);
    }

    private async void OnExportGrid(object? s, RoutedEventArgs e)
    {
        var file = await StorageProvider.SaveFilePickerAsync(new FilePickerSaveOptions
        {
            Title = "Export the grid",
            SuggestedFileName = "density.cube",
            FileTypeChoices = [new FilePickerFileType("Gaussian cube") { Patterns = ["*.cube"] }, new FilePickerFileType("VTK structured grid") { Patterns = ["*.vtk"] },
                               new FilePickerFileType("NumPy array") { Patterns = ["*.npy"] }],
        });
        if (file?.TryGetLocalPath() is { } path) _vm.ExportGrid(path);
    }

    private async void OnFieldFillSuggested(object? s, RoutedEventArgs e)
    {
        if (!_vm.Busy) await _vm.Field.FillSuggested();
    }

    private async void OnFieldFill(object? s, RoutedEventArgs e)
    {
        if (_vm.Busy) return;
        var files = await StorageProvider.OpenFilePickerAsync(new FilePickerOpenOptions
        {
            Title = "Fill missing terms from another force field",
            AllowMultiple = false,
            SuggestedStartLocation = Paths.ForceFields is { } lib ? await StorageProvider.TryGetFolderFromPathAsync(lib) : null,
            FileTypeFilter = [new FilePickerFileType("Force fields") { Patterns = ["*.json", "*.lt", "*.frcmod", "frcmod*", "*.itp", "*.top"] }],
        });
        var path = files.Count > 0 ? files[0].TryGetLocalPath() : null;
        if (path != null) await _vm.Field.FillGaps(path);
    }

    private async void OnFieldImport(object? s, RoutedEventArgs e)
    {
        if (_vm.Busy) return;
        var files = await StorageProvider.OpenFilePickerAsync(new FilePickerOpenOptions
        {
            Title = "Import force-field parameters",
            AllowMultiple = false,
            FileTypeFilter =
            [
                new FilePickerFileType("Parameters (CAPS .json, moltemplate .lt, AMBER frcmod, GROMACS .itp/.top)") { Patterns = ["*.json", "*.lt", "*.frcmod", "frcmod*", "*.dat", "*.itp", "*.top"] },
                new FilePickerFileType("All files") { Patterns = ["*"] },
            ],
        });
        var path = files.Count > 0 ? files[0].TryGetLocalPath() : null;
        if (path != null) await _vm.Field.Import(path);
    }

    private async void OnFieldSaveTypes(object? s, RoutedEventArgs e)
    {
        if (_vm.Document == null) return;
        var file = await StorageProvider.SaveFilePickerAsync(new FilePickerSaveOptions
        {
            Title = "Save atom types",
            SuggestedFileName = Path.GetFileNameWithoutExtension(_vm.Document.Path) + "_types.txt",
            DefaultExtension = "txt",
        });
        var path = file?.TryGetLocalPath();
        if (path == null) return;
        try { _vm.Field.SaveTypes(path); }
        catch (Exception ex) { _vm.Status = "Could not save the types: " + ex.Message; }
    }
    private async void OnReactRun(object? s, RoutedEventArgs e) => await _vm.RunReact();
    private void OnReactCancel(object? s, RoutedEventArgs e) => _vm.CancelReact();
    private void OnPackNew(object? s, RoutedEventArgs e) => _vm.NewPackInput();
    private void OnPackExample(object? s, RoutedEventArgs e) => PackExample();
    public void PackExample() { if (_samples != null) _vm.AddPackExample(_samples); }
    public Task PackAddAsync() => PackAdd();
    public Task PackOpenAsync() => PackOpen();

    /// <summary>Start › From a recipe: a CAPS recipe (run, then opened) or a Packmol input (into Pack).</summary>
    public async Task RecipeOpenAsync()
    {
        var files = await StorageProvider.OpenFilePickerAsync(new FilePickerOpenOptions
        {
            Title = "Import a recipe",
            FileTypeFilter = [new FilePickerFileType("CAPS recipe or Packmol input") { Patterns = ["*.yaml", "*.yml", "*.json", "*.inp"] }, new FilePickerFileType("All files") { Patterns = ["*"] }],
        });
        if (files.Count == 0 || files[0].TryGetLocalPath() is not string p) return;
        if (p.EndsWith(".inp", StringComparison.OrdinalIgnoreCase) || p.EndsWith(".txt", StringComparison.OrdinalIgnoreCase))
        {
            _vm.SetModule(5);
            try { _vm.LoadPackInput(p); }
            catch (Exception ex) { _vm.Status = "Could not open: " + ex.Message; }
            return;
        }
        await _vm.RunRecipeFile(p);
    }
    private async void OnPackRun(object? s, RoutedEventArgs e) => await _vm.RunPack();
    private void OnPackCancel(object? s, RoutedEventArgs e) => _vm.CancelPack();

    private async void OnPackAdd(object? s, RoutedEventArgs e) => await PackAdd();
    private async Task PackAdd()
    {
        var files = await StorageProvider.OpenFilePickerAsync(new FilePickerOpenOptions
        {
            Title = "Add a molecule",
            AllowMultiple = true,
            FileTypeFilter = [new FilePickerFileType("Molecules") { Patterns = ["*.pdb", "*.xyz", "*.data", "*.gro", "*.mol2"] }, new FilePickerFileType("All files") { Patterns = ["*"] }],
        });
        foreach (var f in files) if (f.TryGetLocalPath() is string p) _vm.AddPackStructure(p);
    }

    private async void OnPackOpen(object? s, RoutedEventArgs e) => await PackOpen();
    private async Task PackOpen()
    {
        var files = await StorageProvider.OpenFilePickerAsync(new FilePickerOpenOptions
        {
            Title = "Open packmol input",
            FileTypeFilter = [new FilePickerFileType("packmol input") { Patterns = ["*.inp", "*.txt"] }, new FilePickerFileType("All files") { Patterns = ["*"] }],
        });
        if (files.Count > 0 && files[0].TryGetLocalPath() is string p)
        {
            try { _vm.LoadPackInput(p); }
            catch (Exception ex) { _vm.Status = "Could not open: " + ex.Message; }
        }
    }
    private async void OnEqRun(object? s, RoutedEventArgs e) => await _vm.RunEquilibrate();
    private void OnEqCancel(object? s, RoutedEventArgs e) => _vm.CancelEquilibrate();

    private async void OnSaveTrajectory(object? s, RoutedEventArgs e) => await SaveTrajectoryAsync();
    public async Task SaveTrajectoryAsync()
    {
        if (_vm.Document == null) return;
        var file = await StorageProvider.SaveFilePickerAsync(new FilePickerSaveOptions
        {
            Title = "Save trajectory (LAMMPS dump)",
            SuggestedFileName = $"{Path.GetFileNameWithoutExtension(_vm.Document.Path)}.lammpstrj",
            DefaultExtension = "lammpstrj",
        });
        var path = file?.TryGetLocalPath();
        if (path == null) return;
        try { _vm.SaveTrajectory(path); }
        catch (Exception ex) { _vm.Status = "Save failed: " + ex.Message; }
    }
    private void OnGrowDensityMode(object? s, RoutedEventArgs e) => _vm.GrowUseBox = false;
    private async void OnGrowBuild(object? s, RoutedEventArgs e) => await _vm.Grow();
    private void OnGrowCancel(object? s, RoutedEventArgs e) => _vm.CancelGrow();
    private async void OnSaveData(object? s, RoutedEventArgs e) => await SaveAs("data", "LAMMPS data");
    private async void OnSavePdb(object? s, RoutedEventArgs e) => await SaveAs("pdb", "PDB");
    private async void OnSaveXyz(object? s, RoutedEventArgs e) => await SaveAs("xyz", "Extended XYZ");
    private async void OnSaveCar(object? s, RoutedEventArgs e) => await SaveAs("car", "Materials Studio .car (with its .mdf)");

    /// <summary>Save dialog for the open document (pages call this).</summary>
    public Task SaveAsAsync(string ext, string label) => SaveAs(ext, label);

    private async Task SaveAs(string ext, string label)
    {
        if (_vm.Document == null) return;
        var file = await StorageProvider.SaveFilePickerAsync(new FilePickerSaveOptions
        {
            Title = $"Save as {label}",
            SuggestedFileName = $"{Path.GetFileNameWithoutExtension(_vm.Document.Path)}.{ext}",
            DefaultExtension = ext,
        });
        var path = file?.TryGetLocalPath();
        if (path == null) return;
        try { _vm.SaveDocument(path); }
        catch (Exception ex) { _vm.Status = "Save failed: " + ex.Message; }
    }

    /// <summary>Play or pause the frames (the Studio's frame bar and the trajectory player).</summary>
    public void TogglePlayback() => TogglePlay();

    private void TogglePlay()
    {
        if (_vm.Document == null || !_vm.HasFrames) return;
        _vm.IsPlaying = !_vm.IsPlaying;
        if (_vm.IsPlaying) _playTimer.Start(); else _playTimer.Stop();
    }

    // Viewer keys; ignored while typing in a text box or choosing in a list.
    private void OnKey(object? sender, KeyEventArgs e)
    {
        if (_vm.RecordingShortcut && _vm.RecordShortcutKey(e.Key, e.KeyModifiers)) { e.Handled = true; return; }   // Settings › Your shortcuts
        if (e.Key == Key.K && e.KeyModifiers is KeyModifiers.Meta or KeyModifiers.Control) { TogglePalette(); e.Handled = true; return; }
        if (_vm.PaletteOpen) return;
        if (_vm.TryUserShortcut(e.Key, e.KeyModifiers, FocusManager?.GetFocusedElement() is TextBox)) { e.Handled = true; return; }
        if (e.KeyModifiers is KeyModifiers.Meta or KeyModifiers.Control)
        {
            if (e.Key is Key.OemPlus or Key.Add) { _vm.StepScale(1); e.Handled = true; return; }
            if (e.Key is Key.OemMinus or Key.Subtract) { _vm.StepScale(-1); e.Handled = true; return; }
        }
        if (e.Key == Key.W && e.KeyModifiers is KeyModifiers.Meta or KeyModifiers.Control && _vm.Document != null) { _vm.CloseDocument(); e.Handled = true; return; }
        if (_vm.Document == null || _vm.Busy || FocusManager?.GetFocusedElement() is TextBox or ComboBox) return;
        // editing keys in the Studio: ⌘Z / ⇧⌘Z undo and redo, ⇧E the element picker, Delete the picked atoms
        if (_vm.IsStudio)
        {
            var cmd = OperatingSystem.IsMacOS() ? KeyModifiers.Meta : KeyModifiers.Control;
            if (e.Key == Key.Z && e.KeyModifiers == cmd) { _vm.UndoEdit(false); e.Handled = true; return; }
            if (e.Key == Key.Z && e.KeyModifiers == (cmd | KeyModifiers.Shift)) { _vm.UndoEdit(true); e.Handled = true; return; }
            if (e.Key == Key.E && e.KeyModifiers == KeyModifiers.Shift) { OpenElementPicker(); e.Handled = true; return; }
            if (e.Key == Key.F && e.KeyModifiers == cmd && _vm.HasDocument) { OpenQuery(); e.Handled = true; return; }
            if (e.Key is Key.Delete or Key.Back && e.KeyModifiers == KeyModifiers.None && _vm.HasPicked) { _vm.DeletePicked(); e.Handled = true; return; }
        }
        // the keyboard map (design/boards/InteractionMap)
        {
            var cmd = OperatingSystem.IsMacOS() ? KeyModifiers.Meta : KeyModifiers.Control;
            if (e.KeyModifiers == cmd && e.Key >= Key.D1 && e.Key <= Key.D9) { _vm.GoRailPage(e.Key - Key.D0); e.Handled = true; return; }
            if (e.Key == Key.Tab && e.KeyModifiers == KeyModifiers.Control) { _vm.NextStructure(); e.Handled = true; return; }
            if (e.Key == Key.F1) { _vm.GoModule(41); e.Handled = true; return; }
            if (_vm.IsStudio)
            {
                if (e.Key == Key.I && e.KeyModifiers == cmd) { _vm.InvertSelectionKey(); e.Handled = true; return; }
                if (e.Key == Key.H && e.KeyModifiers == cmd) { _vm.AddHydrogensAll(); e.Handled = true; return; }
                if (e.Key == Key.C && e.KeyModifiers == (cmd | KeyModifiers.Shift)) { _ = _vm.AutoClean(); e.Handled = true; return; }
                if (e.Key == Key.I && e.KeyModifiers == KeyModifiers.Alt) { _vm.InvertPicked(); e.Handled = true; return; }
                if (e.KeyModifiers == KeyModifiers.None)
                    switch (e.Key)
                    {
                        case Key.D1: _vm.SetView(0, 0); e.Handled = true; return;               // front
                        case Key.D2: _vm.SetView(0, Math.PI / 2); e.Handled = true; return;     // top
                        case Key.D3: _vm.SetView(Math.PI / 2, 0); e.Handled = true; return;     // side
                        case Key.D5: _vm.Perspective = !_vm.Perspective; e.Handled = true; return;
                        case Key.OemCloseBrackets when !_vm.HasFocusAtom: _vm.GrowSelectionKey(); e.Handled = true; return;
                    }
            }
        }
        if (_vm.QueryOpen && e.Key == Key.Escape) { _vm.QueryOpen = false; ViewHost.Focus(); e.Handled = true; return; }
        if (_vm.UpdateOpen && e.Key == Key.Escape) { _vm.CloseUpdate(); e.Handled = true; return; }
        if (_vm.ImportOpen && e.Key == Key.Escape) { _vm.CloseImport(); e.Handled = true; return; }
        if (_vm.ExportDialogOpen && e.Key == Key.Escape) { _vm.CancelExport(); _vm.ExportDialogOpen = false; e.Handled = true; return; }
        if (e.Key == Key.E && e.KeyModifiers == (OperatingSystem.IsMacOS() ? KeyModifiers.Meta : KeyModifiers.Control) && _vm.HasDocument) { _vm.OpenExportDialog(); e.Handled = true; return; }
        // Keyboard walk (design/boards/VisAccess): in the 3D view, or anywhere once an atom has the focus ring
        var focused = FocusManager?.GetFocusedElement();
        var walk = _vm.IsStudio && (focused == ViewHost || (_vm.HasFocusAtom && focused is not (ListBox or Slider or TreeView or TabItem)));
        if (e.Key == Key.A && e.KeyModifiers == ((OperatingSystem.IsMacOS() ? KeyModifiers.Meta : KeyModifiers.Control) | KeyModifiers.Shift))
        {
            _vm.AnnounceSelection(); e.Handled = true; return;
        }
        if (walk && e.KeyModifiers is KeyModifiers.None or KeyModifiers.Shift)
        {
            switch (e.Key)
            {
                case Key.Up: _vm.FocusStep(-1); e.Handled = true; return;
                case Key.Down: _vm.FocusStep(1); e.Handled = true; return;
                case Key.OemOpenBrackets: _vm.FocusMolecule(-1); e.Handled = true; return;
                case Key.OemCloseBrackets: _vm.FocusMolecule(1); e.Handled = true; return;
                case Key.B: _vm.FocusBond(); e.Handled = true; return;
                case Key.M when _vm.HasFocusAtom: _vm.FocusMeasure(); e.Handled = true; return;
                case Key.Space when _vm.HasFocusAtom: _vm.FocusSelect(); e.Handled = true; return;
                case Key.Escape when _vm.HasFocusAtom: _vm.ClearFocus(); e.Handled = true; return;
            }
        }
        switch (e.Key)
        {
            case Key.Left: _vm.StepFrame(-1); e.Handled = true; break;
            case Key.Right: _vm.StepFrame(1); e.Handled = true; break;
            case Key.Space: TogglePlay(); e.Handled = true; break;
            case Key.R when e.KeyModifiers == KeyModifiers.None: _vm.ResetView(); e.Handled = true; break;
            case Key.A when e.KeyModifiers == KeyModifiers.None: _vm.ToggleAutoClean(); e.Handled = true; break;
            case Key.F when e.KeyModifiers == KeyModifiers.None: _vm.FrameSelection(); e.Handled = true; break;
            case Key.L when e.KeyModifiers == KeyModifiers.None && _vm.LensOn: _vm.LensHold = true; e.Handled = true; break;
            case Key.Escape: _vm.ClearSelection(); e.Handled = true; break;
        }
    }

    private void OnResetView(object? s, RoutedEventArgs e) => _vm.ResetView();
    private void OnUseExpression(object? s, RoutedEventArgs e) { if (s is Control { Tag: string x }) _vm.UseExpression(x); }
    private void OnSaveExpression(object? s, RoutedEventArgs e) => _vm.SaveExpression();
    private void OnMakeReal(object? s, RoutedEventArgs e) => _vm.MakeReal();
    private void OnFrameSelection(object? s, RoutedEventArgs e) => _vm.FrameSelection();
    private async void OnCopyMolTable(object? s, RoutedEventArgs e)
    {
        if (Clipboard != null) await Clipboard.SetTextAsync(_vm.MolInfoTable);
        _vm.Status = "Copied the molecule's table";
    }
    private void OnAddMolLibrary(object? s, RoutedEventArgs e) => _vm.AddMoleculeToLibrary();
    private void OnHistoryStep(object? s, RoutedEventArgs e) { if ((s as Control)?.Tag is HistoryRow r) _vm.JumpToStep(r.Step); }
    private void OnTakeSnapshot(object? s, RoutedEventArgs e) => _vm.TakeSnapshot();
    private void OnSnapshotCompare(object? s, RoutedEventArgs e) { if ((s as Control)?.Tag is SnapshotRow r) _vm.CompareWithSnapshot(r); }
    private void OnStatesClose(object? s, RoutedEventArgs e) => _vm.StatesOpen = false;
    private void OnStatesSplit(object? s, RoutedEventArgs e) => _vm.OpenStatesInSplit();
    // the project tree: a job folder opens its outputs (twice: its report), an output opens what it is
    private void OnJobFolder(object? s, RoutedEventArgs e) { if ((s as Control)?.Tag is Job j) j.Expanded = !j.Expanded; }
    private void OnJobOutput(object? s, RoutedEventArgs e) { if ((s as Control)?.Tag is JobOutput o) o.Open(); }
    private async void OnStatesCsv(object? s, RoutedEventArgs e)
    {
        var f = await StorageProvider.SaveFilePickerAsync(new Avalonia.Platform.Storage.FilePickerSaveOptions
        {
            Title = "Shift of each atom", SuggestedFileName = "shifts.csv",
            FileTypeChoices = [new Avalonia.Platform.Storage.FilePickerFileType("CSV") { Patterns = ["*.csv"] }],
        });
        if (f?.TryGetLocalPath() is { } p) _vm.ExportShiftsCsv(p);
    }
    private void OnSnapshotMenu(object? s, RoutedEventArgs e) { }
    private void OnSnapshotRestore(object? s, RoutedEventArgs e) { if ((s as Control)?.Tag is SnapshotRow r) _vm.RestoreSnapshot(r); }
    private void OnSnapshotDelete(object? s, RoutedEventArgs e) { if ((s as Control)?.Tag is SnapshotRow r) _vm.DeleteSnapshot(r); }
    private void OnBranchRestore(object? s, RoutedEventArgs e) { if ((s as Control)?.Tag is BranchRow r) _vm.RestoreBranch(r); }
    private void OnBranchDelete(object? s, RoutedEventArgs e) { if ((s as Control)?.Tag is BranchRow r) _vm.DeleteBranch(r); }
    /// <summary>⌘F: the query bar over the view, focused.</summary>
    public void OpenQuery()
    {
        _vm.QueryOpen = true;
        Dispatcher.UIThread.Post(() => { QueryBox.Focus(); QueryBox.SelectAll(); }, DispatcherPriority.Input);
    }
    private void OnQueryKey(object? s, KeyEventArgs e)
    {
        if (e.Key != Key.Enter) return;
        _vm.ApplyQuery(e.KeyModifiers.HasFlag(KeyModifiers.Shift) ? "add" : "replace");
        e.Handled = true;
    }
    private void OnQuerySelect(object? s, RoutedEventArgs e) => _vm.ApplyQuery();
    private void OnQueryGroup(object? s, RoutedEventArgs e) => _vm.QueryAsGroup();
    private void OnQuerySave(object? s, RoutedEventArgs e) => _vm.SaveQuery(QueryName.Text);
    private void OnQueryClose(object? s, RoutedEventArgs e) => _vm.QueryOpen = false;
    private void OnSavedQuery(object? s, RoutedEventArgs e) { if ((s as Control)?.Tag is SavedQueryRow r) _vm.UseSavedQuery(r); }
    private void OnSavedQueryDelete(object? s, RoutedEventArgs e) { if ((s as Control)?.Tag is SavedQueryRow r) _vm.DeleteSavedQuery(r); }
    private void OnNoticeClose(object? s, RoutedEventArgs e) { if ((s as Control)?.Tag is Notice n) _vm.DismissNotice(n); }
    private void OnNoticePrimary(object? s, RoutedEventArgs e) { if ((s as Control)?.Tag is Notice n) _vm.NoticePrimary(n); }
    private void OnNoticeSecondary(object? s, RoutedEventArgs e) { if ((s as Control)?.Tag is Notice n) _vm.NoticeSecondary(n); }
    private void OnViewTop(object? s, RoutedEventArgs e) => _vm.SetView(0, Math.PI / 2);
    private void OnViewFront(object? s, RoutedEventArgs e) => _vm.SetView(0, 0);
    private void OnViewSide(object? s, RoutedEventArgs e) => _vm.SetView(Math.PI / 2, 0);

    // ---------------------------------------------------------------- rendering

    private void RequestRender()
    {
        _sceneDirty = true;   // content may have changed: the GPU view rebuilds its scene
        RequestViewRender();
    }

    /// <summary>Only the camera moved (drag, wheel, flight): the GPU view turns the scene it has.</summary>
    private void RequestViewRender()
    {
        if (_vm.IsRender) RenderGuide.InvalidateVisual();
        _requested++;
        if (!_busy) _ = RenderLoop();
    }

    // ---- GPU view (GlMolView): the scene is uploaded on content changes; turning only sends the camera. Picking and the
    // label overlay use the CPU renderer's id buffer, refreshed off screen once the view stops moving, or before a click.
    private bool _sceneDirty = true, _gpuCpuOnly, _cpuStale, _gpuTried;
    private object? _sceneDoc;
    private DispatcherTimer? _idle;
    private byte[] _pickBuf = [];

    private bool GpuPath()
    {
        if (!_vm.GpuView || _vm.IsRender) return false;
        if (!ViewGl.Ready)
        {
            if (!_gpuTried) { _gpuTried = true; ViewGl.IsVisible = true; }   // shown once: OpenGL starts, ReadyChanged follows
            return false;
        }
        return true;
    }

    private void RestartIdle()
    {
        if (_idle == null)
        {
            _idle = new DispatcherTimer { Interval = TimeSpan.FromMilliseconds(250) };
            _idle.Tick += async (_, _) =>
            {
                _idle!.Stop();
                if (_cpuStale && !_busy) await RefreshPickBuffer();
            };
        }
        _idle.Stop();
        _idle.Start();
    }

    /// <summary>The CPU id buffer for the GPU view's last camera (off the UI thread), then the overlays from it.</summary>
    private async Task RefreshPickBuffer()
    {
        if (_vm.Document is not { } doc || _vm.Busy || !_haveLast) return;
        var cam = _lastCam;
        var opt = _lastOpt;
        var n = _pixW * _pixH * 4;
        if (_pickBuf.Length != n) _pickBuf = new byte[n];
        var buf = _pickBuf;
        try { await Task.Run(() => doc.Render(cam, opt, buf)); } catch { return; }
        if (!ReferenceEquals(doc, _vm.Document) || !SameCamera(cam, _lastCam)) return;   // moved meanwhile: the next idle does it
        _cpuStale = false;
        UpdateOverlays(cam, opt);
    }

    /// <summary>Before a pick: the id buffer for what the GPU view shows now.</summary>
    private void EnsurePickBuffer()
    {
        if (!_cpuStale || _vm.Document is not { } doc || !_haveLast || _vm.Busy) return;
        var n = _pixW * _pixH * 4;
        if (_pickBuf.Length != n) _pickBuf = new byte[n];
        try { doc.Render(_lastCam, _lastOpt, _pickBuf); _cpuStale = false; } catch { }
    }

    private static bool SameCamera(in CapsStudio.Interop.CapsCamera a, in CapsStudio.Interop.CapsCamera b) =>
        a.Yaw == b.Yaw && a.Pitch == b.Pitch && a.Zoom == b.Zoom && a.PanX == b.PanX && a.PanY == b.PanY && a.Perspective == b.Perspective;

    private CapsStudio.Interop.CapsSceneData? _gpuScene;
    private double _panSpan = 30;

    private void ClearOverlays()
    {
        Labels.SetLabels(new List<ViewModels.ViewLabel>());
        Labels.SetLens(null);
        Labels.SetMonitors(new List<ViewModels.MonitorMark>());
    }

    /// <summary>Atom labels, the lens and pinned monitors over the view for this camera.</summary>
    private void UpdateOverlays(CapsStudio.Interop.CapsCamera cam, CapsStudio.Interop.CapsRenderOpts opt)
    {
        try { Labels.SetLabels(_vm.AnyLabels && !_vm.IsVisualize ? _vm.ViewLabels(cam, opt, _scaling) : new List<ViewModels.ViewLabel>()); }
        catch { Labels.SetLabels(new List<ViewModels.ViewLabel>()); }
        Labels.SetLens(_vm.LensCircle(cam, opt, _scaling));
        try { Labels.SetMonitors(_vm.MonitorMarks(cam, opt, _scaling)); } catch { Labels.SetMonitors(new List<ViewModels.MonitorMark>()); }
    }

    private byte[] _frameBuf = [];
    private readonly WriteableBitmap?[] _frames = new WriteableBitmap?[3];   // main view A/B, Field view
    private bool _frameFlip;

    private async Task RenderLoop()
    {
        _busy = true;
        try
        {
            while (_rendered != _requested)
            {
                var ticket = _requested;
                var doc = _vm.Document;
                // the Field page has its own view (coloured by force-field type); otherwise the main viewport
                var field = _vm.IsField;
                var host = field ? FieldViewHost : ViewHost;
                var image = field ? FieldImage : ViewImage;
                if (field && !_vm.Field.Assigned) { FieldImage.Source = null; _rendered = ticket; break; }
                var w = (int)Math.Max(16, host.Bounds.Width);
                var h = (int)Math.Max(16, host.Bounds.Height);
                _scaling = VisualRoot?.RenderScaling ?? 1;
                if (doc == null) { image.Source = null; _rendered = ticket; break; }
                // While a run holds the document: its live snapshots (MD, equilibration) are drawn instead, and the GPU
                // view turns the scene it has; nothing here waits on the document.
                var busy = _vm.Busy;
                var src = busy ? _vm.RunLiveDoc : doc;
                var gpu = !field && GpuPath();
                if (busy && (field || (src == null && !(gpu && _gpuScene != null && !_gpuCpuOnly)))) { _rendered = ticket; break; }
                var pw = (int)(w * _scaling);
                var ph = (int)(h * _scaling);
                var cam = busy ? _vm.Camera : _vm.ViewCamera(w, h);
                if (gpu)
                {
                    var gopt = _vm.ViewOptions(pw, ph, 1);
                    if (src != null && (_sceneDirty || !ReferenceEquals(_sceneDoc, src)))
                    {
                        _sceneDirty = false;   // a change during the build marks it again
                        CapsStudio.Interop.CapsSceneData? sc;
                        try { sc = await Task.Run(() => src.RenderScene(gopt)); } catch { sc = null; }   // a snapshot may be replaced meanwhile
                        if (sc != null)
                        {
                            _gpuCpuOnly = sc.CpuOnly;   // surfaces, polyhedra, colour-vision preview: the CPU draws them
                            if (!_gpuCpuOnly) { ViewGl.SetScene(sc); _gpuScene = sc; _sceneDoc = src; }
                            else { _gpuScene = null; _sceneDoc = null; }
                        }
                        else if (!busy) _sceneDirty = true;
                    }
                    if (!_gpuCpuOnly && _gpuScene != null && ReferenceEquals(doc, _vm.Document))
                    {
                        var gsw = Stopwatch.StartNew();
                        var fit = _gpuScene.Fit(cam, pw, ph);   // from the scene: no call into the document
                        ViewGl.SetView(fit);
                        if (!ViewGl.IsVisible) ViewGl.IsVisible = true;
                        if (ViewImage.IsVisible) ViewImage.IsVisible = false;
                        if (!busy)
                        {
                            var cpuOpt = _vm.ViewOptions(pw, ph, _scaling >= 1.5 ? 1 : 2);   // what picks and overlays are made with
                            _pixW = pw; _pixH = ph;
                            _lastCam = cam; _lastOpt = cpuOpt; _haveLast = true;
                            _cpuStale = true;
                            UpdateOverlays(cam, cpuOpt);
                            RestartIdle();
                        }
                        else if (_haveLast) { _cpuStale = true; ClearOverlays(); }   // picks and labels wait for the run to end
                        RenderStat.Text = $"{pw}×{ph} px · GPU{(busy ? src != null ? " · live" : " · run" : "")} · {gsw.Elapsed.TotalMilliseconds + ViewGl.LastFrameMs:0} ms";
                        _rendered = ticket;
                        continue;
                    }
                    if (busy && src == null) { _rendered = ticket; break; }
                }
                if (!field && ViewGl.IsVisible && ViewGl.Ready) { ViewGl.IsVisible = false; }
                if (!field && !ViewImage.IsVisible) ViewImage.IsVisible = true;
                // Retina already gives 2 samples per point; supersample only on 1× screens.
                var opt = field ? _vm.FieldViewOptions(pw, ph, _scaling >= 1.5 ? 1 : 2) : _vm.ViewOptions(pw, ph, _scaling >= 1.5 ? 1 : 2);
                // the pixel buffer and two bitmaps are kept while the size holds (a drag renders tens of frames a second;
                // new 6 MB buffers each time would keep the garbage collector busy and the rotation uneven)
                if (_frameBuf.Length != pw * ph * 4) _frameBuf = new byte[pw * ph * 4];
                var buf = _frameBuf;
                var sw = Stopwatch.StartNew();
                var rdoc = src!;   // the document, or during a run its live snapshot
                try { await Task.Run(() => rdoc.Render(cam, opt, buf)); }
                catch (ObjectDisposedException) when (busy) { _rendered = ticket; continue; }   // a snapshot replaced meanwhile
                sw.Stop();
                var slot = field ? 2 : _frameFlip ? 1 : 0;
                _frameFlip = !_frameFlip;
                var bmp = _frames[slot];
                if (bmp == null || bmp.PixelSize.Width != pw || bmp.PixelSize.Height != ph || Math.Abs(bmp.Dpi.X - 96 * _scaling) > 1e-6)
                {
                    if (bmp != null && !ReferenceEquals(image.Source, bmp)) bmp.Dispose();
                    bmp = _frames[slot] = new WriteableBitmap(new PixelSize(pw, ph), new Vector(96 * _scaling, 96 * _scaling), PixelFormat.Rgba8888, AlphaFormat.Unpremul);
                }
                using (var fb = bmp.Lock())
                {
                    if (fb.RowBytes == pw * 4) System.Runtime.InteropServices.Marshal.Copy(buf, 0, fb.Address, pw * ph * 4);
                    else
                        for (var y = 0; y < ph; y++)
                            System.Runtime.InteropServices.Marshal.Copy(buf, y * pw * 4, fb.Address + y * fb.RowBytes, pw * 4);
                }
                if (ReferenceEquals(image.Source, bmp)) image.InvalidateVisual();
                else image.Source = bmp;
                _pixW = pw; _pixH = ph;
                RenderStat.Text = $"{pw}×{ph} px · {sw.ElapsedMilliseconds} ms";
                if (!field && !busy) _vm.ReportFrame(sw.Elapsed.TotalMilliseconds);
                if (busy) { ClearOverlays(); _rendered = ticket; continue; }   // a snapshot: no labels or picks from it
                // atom labels (Appearance): the visible atoms' screen positions from this render
                if (!field)
                {
                    try { Labels.SetLabels(_vm.AnyLabels && !_vm.IsVisualize ? _vm.ViewLabels(cam, opt, _scaling) : new List<ViewModels.ViewLabel>()); }
                    catch { Labels.SetLabels(new List<ViewModels.ViewLabel>()); }
                    Labels.SetLens(_vm.LensCircle(cam, opt, _scaling));
                    try { Labels.SetMonitors(_vm.MonitorMarks(cam, opt, _scaling)); } catch { Labels.SetMonitors(new List<ViewModels.MonitorMark>()); }
                    _lastCam = cam; _lastOpt = opt; _haveLast = true;
                }
                _rendered = ticket;
            }
        }
        catch (Exception e)
        {
            _vm.Status = "Render failed: " + e.Message;
            _rendered = _requested;
        }
        finally
        {
            _busy = false;
        }
    }

    // ---------------------------------------------------------------- builder tools and the Element picker

    private void OnToolPlace(object? s, RoutedEventArgs e) => _vm.EditTool = _vm.EditTool == 1 ? 0 : 1;
    private void OnToolBond(object? s, RoutedEventArgs e) => _vm.EditTool = _vm.EditTool == 2 ? 0 : 2;
    private void OnToolDelete(object? s, RoutedEventArgs e) { if (_vm.HasPicked && _vm.EditTool != 3) _vm.DeletePicked(); else _vm.EditTool = _vm.EditTool == 3 ? 0 : 3; }
    private void OnAddHydrogens(object? s, RoutedEventArgs e) => _vm.AddHydrogensAll();
    private void OnInvert(object? s, RoutedEventArgs e) => _vm.InvertPicked();
    private void OnFuseRing(object? s, RoutedEventArgs e) => _vm.FuseRingPicked();
    private void OnToolLasso(object? s, RoutedEventArgs e) => _vm.EditTool = _vm.EditTool == 4 ? 0 : 4;
    private void OnToolMove(object? s, RoutedEventArgs e) => _vm.EditTool = _vm.EditTool == 5 ? 0 : 5;
    private void OnPinMonitor(object? s, RoutedEventArgs e) { _vm.PinMeasurement(); RequestRender(); }
    private void OnUnpin(object? s, RoutedEventArgs e) { if ((s as Control)?.Tag is ViewModels.MonitorRow m) { _vm.UnpinMonitor(m); RequestRender(); } }

    /// <summary>Even–odd rule: is (x, y) inside the polygon?</summary>
    private static bool InPolygon(List<Point> poly, double x, double y)
    {
        var inside = false;
        for (int i = 0, j = poly.Count - 1; i < poly.Count; j = i++)
            if ((poly[i].Y > y) != (poly[j].Y > y) && x < (poly[j].X - poly[i].X) * (y - poly[i].Y) / (poly[j].Y - poly[i].Y) + poly[i].X) inside = !inside;
        return inside;
    }
    private void OnFragmentTool(object? s, RoutedEventArgs e) => _vm.OpenFragments();
    private async void OnAutoClean(object? s, RoutedEventArgs e) => await _vm.AutoClean();
    private void OnUndo(object? s, RoutedEventArgs e) => _vm.UndoEdit(false);
    private void OnRedo(object? s, RoutedEventArgs e) => _vm.UndoEdit(true);
    private void OnElementPicker(object? s, RoutedEventArgs e) => OpenElementPicker();
    private void OpenElementPicker()
    {
        _vm.ElementPickerOpen = true;
        MarkChosenElement();
        PickerSearch.Text = "";
        PickerSearch.Focus();
    }
    private void OnPickerClose(object? s, RoutedEventArgs e) => _vm.ElementPickerOpen = false;
    public void MarkChosenForTest() => MarkChosenElement();
    private void OnPickerBackdrop(object? s, PointerPressedEventArgs e) => _vm.ElementPickerOpen = false;
    private void OnPickerCell(object? s, RoutedEventArgs e) { if ((s as Control)?.Tag is string sym) { _vm.BuildElement = sym; MarkChosenElement(); } }
    private void OnPickerType(object? s, KeyEventArgs e)
    {
        if (e.Key == Key.Escape) { _vm.ElementPickerOpen = false; return; }
        if (e.Key == Key.Enter) { _vm.EditTool = 1; _vm.ElementPickerOpen = false; return; }
        if (_vm.PickerJump(PickerSearch.Text ?? "")) MarkChosenElement();
    }
    private void OnPickerReplace(object? s, RoutedEventArgs e) { _vm.ReplacePickedElement(); _vm.ElementPickerOpen = false; }
    private void OnPickerPlace(object? s, RoutedEventArgs e) { _vm.EditTool = 1; _vm.ElementPickerOpen = false; }
    /// <summary>The chosen element's cell amber (a class, so the common/rare colours stay for the rest).</summary>
    private void MarkChosenElement()
    {
        foreach (var b in this.GetVisualDescendants().OfType<Button>().Where(b => b.Classes.Contains("pt")))
            b.Classes.Set("chosen", b.Tag as string == _vm.BuildElement);
    }

    // ---------------------------------------------------------------- selection & stereo (design/boards/SelectionStereo)

    private void OnSelOp(object? s, RoutedEventArgs e) { if ((s as Control)?.Tag is string op) _vm.RunSelect(op); }
    private void OnSelApply(object? s, RoutedEventArgs e) => _vm.RunSelect("replace");
    private void OnSelClear(object? s, RoutedEventArgs e) => _vm.ClearDocSelection();
    private void OnSaveSet(object? s, RoutedEventArgs e) => _vm.SaveSelectionAsSet();
    private void OnUseSet(object? s, RoutedEventArgs e) { if ((s as Control)?.Tag is ViewModels.NamedSet n) _vm.UseNamedSet(n); }
    private async void OnMakeIso(object? s, RoutedEventArgs e) => await _vm.MakeTactic(true);
    private async void OnMakeSyndio(object? s, RoutedEventArgs e) => await _vm.MakeTactic(false);

    // ---------------------------------------------------------------- interactions & checks (design/boards/Interactions)

    private async void OnIxFix(object? s, RoutedEventArgs e) { if ((s as Control)?.Tag is ViewModels.CheckIssue i) await _vm.FixIssue(i); }
    private async void OnIxFixAll(object? s, RoutedEventArgs e) => await _vm.FixAllSafe();
    private async void OnIxExport(object? s, RoutedEventArgs e)
    {
        var file = await StorageProvider.SaveFilePickerAsync(new FilePickerSaveOptions
        {
            Title = "Save the H-bonds and clashes", SuggestedFileName = "interactions.csv",
            FileTypeChoices = [new FilePickerFileType("CSV") { Patterns = ["*.csv"] }],
        });
        if (file?.TryGetLocalPath() is { } path) { File.WriteAllText(path, _vm.InteractionsCsv()); _vm.Status = "Saved " + path; }
    }

    // ---------------------------------------------------------------- appearance (design/boards/Appearance)

    private void OnAppStyle(object? s, RoutedEventArgs e)
    {
        if ((s as Control)?.Tag is string t && int.TryParse(t, out var k)) _vm.AppStyle = k;
    }
    private void OnRemoveLayer(object? s, RoutedEventArgs e) { if ((s as Control)?.Tag is ViewModels.StyleLayer l) _vm.RemoveStyleLayer(l); }
    private void OnAppearanceReset(object? s, RoutedEventArgs e) => _vm.ResetAppearance();
    private void OnAppearanceExport(object? s, RoutedEventArgs e) => _vm.OpenFigure();
    private void OnAppearanceApply(object? s, RoutedEventArgs e) => _vm.AppearanceOpen = false;

    // ---------------------------------------------------------------- mouse

    private void OnPointerPressed(object? sender, PointerPressedEventArgs e)
    {
        _vm.StopFly();   // the hand wins over a camera flight
        _hostField = sender as Control ?? ViewHost;
        var p = e.GetCurrentPoint(_host);
        _press = _last = p.Position;
        _dragging = true;
        _moved = false;
        _pan = p.Properties.IsRightButtonPressed || e.KeyModifiers.HasFlag(KeyModifiers.Alt);
        _addPick = e.KeyModifiers.HasFlag(KeyModifiers.Shift) || e.KeyModifiers.HasFlag(KeyModifiers.Meta) || _vm.MeasureTool;
        _lassoPts = null;
        _moveAtoms = null;
        if (_host == ViewHost && _vm.IsStudio && _vm.Document is { } doc && !_vm.Busy && !_pan)
        {
            if (_vm.EditTool == 4) _lassoPts = new List<Point> { p.Position };
            else if (_vm.EditTool == 5 && _haveLast)
            {
                try
                {
                    EnsurePickBuffer();
                    var hit = doc.Pick((int)(p.Position.X * _scaling), (int)(p.Position.Y * _scaling));
                    var set = _vm.MoveSet(hit);
                    if (set.Length > 0)
                    {
                        // the view plane from the visible atoms' projections (a sample of them)
                        var n = (int)doc.Summary().Atoms;
                        var proj = doc.ProjectAtoms(_lastCam, _lastOpt, n);
                        var fit = new List<(double, double, double, double, double)>();
                        var step = Math.Max(1, n / 600);
                        for (var i = 0; i < n; i += step)
                        {
                            if (proj[3 * i + 2] <= 0) continue;
                            var a = doc.Atom(i);
                            fit.Add((a.X, a.Y, a.Z, proj[3 * i] / _scaling, proj[3 * i + 1] / _scaling));
                        }
                        if (fit.Count >= 4) { _moveAtoms = set; _moveFit = fit; }
                    }
                }
                catch { _moveAtoms = null; }
            }
        }
        e.Pointer.Capture(_host);
        _host.Focus();
    }

    private void OnPointerMoved(object? sender, PointerEventArgs e)
    {
        if (!_dragging && _vm.LensHold && _vm.Document != null && !_vm.Busy)   // L held: the lens follows the atom under the cursor
        {
            var at = e.GetPosition(_host);
            EnsurePickBuffer();
            var hit = _vm.Document.Pick((int)(at.X * _scaling), (int)(at.Y * _scaling));
            if (hit >= 0) { _vm.MoveLens(hit); RequestRender(); }
            return;
        }
        if (!_dragging || _vm.Document == null) return;   // the camera also turns while a run goes (tools wait for it)
        var pos = e.GetPosition(_host);
        if (_lassoPts != null || _moveAtoms != null)   // the lasso or move tool owns the drag: the camera stays
        {
            if (Math.Abs(pos.X - _press.X) + Math.Abs(pos.Y - _press.Y) > 3) _moved = true;
            if (_lassoPts != null)
            {
                if (_lassoPts.Count == 0 || Math.Abs(pos.X - _lassoPts[^1].X) + Math.Abs(pos.Y - _lassoPts[^1].Y) > 2) _lassoPts.Add(pos);
                Labels.SetLasso(_lassoPts);
            }
            else Labels.SetArrow((_press, pos));
            _last = pos;
            return;
        }
        var d = pos - _last;
        _last = pos;
        if (Math.Abs(pos.X - _press.X) + Math.Abs(pos.Y - _press.Y) > 3) _moved = true;
        if (!_moved) return;
        if (_pan)
        {
            // Pan in Å: approximate px→Å from the current fit (view height ≈ 2.2 × half-extent at zoom 1).
            if (!_vm.Busy) { var s = Summary(); _panSpan = Math.Max(1.0, Math.Max(s.CellA, Math.Max(s.CellB, s.CellC))); }   // a run holds the document: the last span
            var span = _panSpan;
            var perPx = span * 1.6 / Math.Max(1, _host.Bounds.Height) / _vm.Camera.Zoom;
            _vm.Camera.PanX += d.X * perPx;
            _vm.Camera.PanY -= d.Y * perPx;
        }
        else
        {
            _vm.Camera.Yaw += d.X * 0.008;
            _vm.Camera.Pitch = Math.Clamp(_vm.Camera.Pitch + d.Y * 0.008, -Math.PI / 2, Math.PI / 2);
        }
        RequestViewRender();
    }

    private Interop.CapsSummary Summary() => _vm.Document?.Summary() ?? default;

    private void OnPointerReleased(object? sender, PointerReleasedEventArgs e)
    {
        e.Pointer.Capture(null);
        if (_dragging && _moved && _lassoPts != null && _vm.Document is { } ldoc)
        {
            // the visible atoms whose centres are inside the lasso polygon
            var poly = _lassoPts;
            var inside = new List<int>();
            if (poly.Count >= 3 && _haveLast)
            {
                EnsurePickBuffer();
                var n = (int)ldoc.Summary().Atoms;
                var proj = ldoc.ProjectAtoms(_lastCam, _lastOpt, n);
                for (var i = 0; i < n; ++i)
                    if (proj[3 * i + 2] > 0 && InPolygon(poly, proj[3 * i] / _scaling, proj[3 * i + 1] / _scaling)) inside.Add(i);
            }
            _vm.LassoSelect(inside, e.KeyModifiers.HasFlag(KeyModifiers.Shift));
            _lassoPts = null;
            Labels.SetLasso(new List<Point>());
            _dragging = false;
            RequestRender();
            return;
        }
        if (_dragging && _moved && _moveAtoms != null && _moveFit != null)
        {
            var pos = e.GetPosition(_host);
            var by = ViewModels.MainViewModel.ScreenToWorld(_moveFit, pos.X - _press.X, pos.Y - _press.Y);
            if (by != null) _vm.TranslateAtoms(_moveAtoms, by);
            _moveAtoms = null;
            Labels.SetArrow(null);
            _dragging = false;
            RequestRender();
            return;
        }
        _lassoPts = null;
        _moveAtoms = null;
        Labels.SetLasso(new List<Point>());
        Labels.SetArrow(null);
        if (_dragging && !_moved && _vm.Document != null && !_vm.Busy)
        {
            var pos = e.GetPosition(_host);
            if (_host == ViewHost) EnsurePickBuffer();
            var hit = _vm.Document.Pick((int)(pos.X * _scaling), (int)(pos.Y * _scaling));
            if (!_vm.PickAllowed(hit)) hit = -1;   // "Measurements only inside" the lens
            if (_host == FieldViewHost) { if (hit >= 0) _vm.Field.SelectAtom(hit); }
            else if (_vm.EditTool != 0 && _vm.IsStudio) _vm.ToolClick(hit);
            else _vm.Pick(hit, _addPick);
            RequestRender();
        }
        _dragging = false;
    }

    private void OnWheel(object? sender, PointerWheelEventArgs e)
    {
        if (_vm.Document == null) return;
        _vm.StopFly();
        _vm.Camera.Zoom = Math.Clamp(_vm.Camera.Zoom * Math.Pow(1.12, e.Delta.Y), 0.1, 40);
        RequestViewRender();
    }

    // ---------------------------------------------------------------- export

    private async Task Export(string ext)
    {
        if (_vm.Document == null || _vm.Busy) return;
        var bgName = MainViewModel.Backgrounds[_vm.ExportBackground].ToLowerInvariant();
        var file = await StorageProvider.SaveFilePickerAsync(new FilePickerSaveOptions
        {
            Title = $"Export figure ({ext.ToUpperInvariant()})",
            SuggestedFileName = $"{Path.GetFileNameWithoutExtension(_vm.Document.Path)}_{bgName}.{ext}",
            DefaultExtension = ext,
        });
        var path = file?.TryGetLocalPath();
        if (path == null) return;
        var (w, h) = _vm.ExportSize((int)ViewHost.Bounds.Width, (int)ViewHost.Bounds.Height);
        var opt = _vm.ExportOptions(w, h);
        var cam = _vm.Camera;
        var doc = _vm.Document;
        try
        {
            _vm.Status = $"Exporting {w}×{h} {bgName} {ext.ToUpperInvariant()}…";
            await Task.Run(() => { if (ext == "svg") doc.ExportSvg(cam, opt, path); else doc.ExportPng(cam, opt, path); });
            _vm.Status = $"Wrote {path} · {w}×{h} · {bgName} background";
        }
        catch (Exception ex)
        {
            _vm.Status = "Export failed: " + ex.Message;
        }
    }

    private void OnChecks(object? s, RoutedEventArgs e) => ViewModel.OpenChecks();
    private async void OnExportPng(object? s, RoutedEventArgs e) => await Export("png");
    private void OnExportDialog(object? s, RoutedEventArgs e) => _vm.OpenExportDialog();
    private async void OnExportSvg(object? s, RoutedEventArgs e) => await Export("svg");
    // ---- Analyze › Visualize
    private void OnOpenVisualize(object? s, RoutedEventArgs e) => _vm.OpenVisualize();
    private void OnAddStep(object? s, RoutedEventArgs e) => _vm.StepLibraryOpen = true;
    private void OnCloseStepLibrary(object? s, RoutedEventArgs e) => _vm.StepLibraryOpen = false;
    private void OnPickStep(object? s, RoutedEventArgs e) { if ((s as Control)?.Tag is StepKind k) _vm.AddStep(k.Type); }
    private void OnPipeStep(object? s, RoutedEventArgs e) { if ((s as Control)?.Tag is PipelineRow r) _vm.PipeSelected = r; }
    private void OnPipeSource(object? s, RoutedEventArgs e) => _vm.PipeSelected = null;
    private void OnGroupFold(object? s, RoutedEventArgs e) { if ((s as Control)?.Tag is PipelineRow r) r.GroupFolded?.Invoke(r); }
    private void OnStepUp(object? s, RoutedEventArgs e) { if (_vm.PipeSelected is { } r) _vm.MoveStep(r, -1); }
    private void OnStepDown(object? s, RoutedEventArgs e) { if (_vm.PipeSelected is { } r) _vm.MoveStep(r, 1); }
    private void OnStepRemove(object? s, RoutedEventArgs e) { if (_vm.PipeSelected is { } r) _vm.RemoveStep(r); }
    private async void OnSeries(object? s, RoutedEventArgs e) => await _vm.ComputeSeries();
    private async void OnStepFile(object? s, RoutedEventArgs e)
    {
        if ((s as Control)?.Tag is not StepField f) return;
        var files = await StorageProvider.OpenFilePickerAsync(new FilePickerOpenOptions
        {
            Title = f.Label, AllowMultiple = false,
            FileTypeFilter = [new FilePickerFileType("Python") { Patterns = ["*.py"] }],
        });
        if (files.Count > 0 && files[0].TryGetLocalPath() is { } p) f.Text = p;
    }
    private void OnStepAction(object? s, RoutedEventArgs e) { if ((s as Control)?.Tag is StepField f) _vm.StepAction(f); }
    private async void OnOverlayChoose(object? s, RoutedEventArgs e)
    {
        var files = await StorageProvider.OpenFilePickerAsync(new FilePickerOpenOptions
        {
            Title = "Python overlay", AllowMultiple = false, FileTypeFilter = [new FilePickerFileType("Python") { Patterns = ["*.py"] }],
        });
        if (files.Count > 0 && files[0].TryGetLocalPath() is { } p) { _vm.OvScript = p; _vm.OvPython = true; }
    }
    private void OnOverlayNew(object? s, RoutedEventArgs e) => _vm.NewOverlayScript();
    private void OnComputeSettings(object? s, RoutedEventArgs e) => _vm.OpenComputeSettings();
    private void OnOverlayEdit(object? s, RoutedEventArgs e)
    {
        if (_vm.OvScript.Length == 0) return;
        try { System.Diagnostics.Process.Start(new System.Diagnostics.ProcessStartInfo(_vm.OvScript) { UseShellExecute = true }); }
        catch (Exception ex) { _vm.Status = "Cannot open the script: " + ex.Message; }
    }
    private void OnInspectorPrev(object? s, RoutedEventArgs e) => _vm.InspectorPageStep(-1);
    private void OnInspectorNext(object? s, RoutedEventArgs e) => _vm.InspectorPageStep(1);

    private async void OnPipelineMenu(object? s, RoutedEventArgs e)
    {
        var menu = new ContextMenu
        {
            ItemsSource = new[]
            {
                new MenuItem { Header = "Save pipeline (YAML)…", Command = new RelayCommand(() => { _vm.OpenSavePipeline(); return Task.CompletedTask; }) },
                new MenuItem { Header = "Save pipeline as JSON…", Command = new RelayCommand(SavePipelineAsync) },
                new MenuItem { Header = "Load pipeline…", Command = new RelayCommand(LoadPipelineAsync) },
                new MenuItem { Header = "Colour by… (gallery)", Command = new RelayCommand(() => { _vm.OpenColourBy(); return Task.CompletedTask; }) },
                new MenuItem { Header = "Compare with another file…", Command = new RelayCommand(() => { _vm.OpenCompare(); return Task.CompletedTask; }) },
                new MenuItem { Header = "Run on many files…", Command = new RelayCommand(() => { _vm.OpenBatch(); return Task.CompletedTask; }) },
                new MenuItem { Header = "Figure bundle…", Command = new RelayCommand(() => { _vm.OpenBundle(); return Task.CompletedTask; }) },
                new MenuItem { Header = "Export the result…", Command = new RelayCommand(() => { _vm.OpenExport(); return Task.CompletedTask; }) },
                new MenuItem { Header = "Clear all steps", Command = new RelayCommand(() => { _vm.ClearPipeline(); return Task.CompletedTask; }) },
            },
        };
        menu.Open(s as Control);
        await Task.CompletedTask;
    }

    private async Task SavePipelineAsync()
    {
        var f = await StorageProvider.SaveFilePickerAsync(new FilePickerSaveOptions
        {
            Title = "Save pipeline", DefaultExtension = "json", SuggestedFileName = "pipeline.json",
            FileTypeChoices = [new FilePickerFileType("CAPS pipeline (JSON)") { Patterns = ["*.json"] }],
        });
        if (f?.TryGetLocalPath() is { } path) _vm.SavePipeline(path);
    }

    private async Task LoadPipelineAsync()
    {
        var files = await StorageProvider.OpenFilePickerAsync(new FilePickerOpenOptions
        {
            Title = "Load pipeline", AllowMultiple = false,
            FileTypeFilter = [new FilePickerFileType("CAPS pipeline") { Patterns = ["*.yaml", "*.yml", "*.json"] }],
        });
        if (files.Count > 0 && files[0].TryGetLocalPath() is { } path) _vm.LoadPipeline(path);
    }

    private void OnRenderPage(object? s, RoutedEventArgs e) => _vm.OpenRender();
    private void OnProvenancePage(object? s, RoutedEventArgs e) => _vm.OpenProvenance();
    private async void OnBenchmarkView(object? s, RoutedEventArgs e) => await _vm.BenchmarkView(_pixW > 0 ? _pixW : 1000, _pixH > 0 ? _pixH : 700);
    private void OnUnitTyping(object? s, RoutedEventArgs e) => _vm.OpenUnitTyping();
    private void OnStepCodeRun(object? s, RoutedEventArgs e) { if ((s as Control)?.Tag is StepField f) f.Commit(); }
    private void OnStepCodeKey(object? s, KeyEventArgs e)
    {
        if (e.Key == Key.Enter && e.KeyModifiers.HasFlag(KeyModifiers.Meta) || e.Key == Key.Enter && e.KeyModifiers.HasFlag(KeyModifiers.Control))
            if ((s as Control)?.Tag is StepField f) { f.Commit(); e.Handled = true; }
    }

    private void OnProjectPage(object? s, RoutedEventArgs e) => _vm.OpenProject();

    /// <summary>Compact layout: rail and toolbar icons only, the inspector and project as drawers, the dock folded.</summary>
    private void ApplyCompact()
    {
        var c = _vm.Compact;
        Rail.Classes.Set("compact", c);
        ToolbarRight.Classes.Set("compact", c);
        Body.ColumnDefinitions[0].Width = new GridLength(c ? 52 : 72);
        Body.ColumnDefinitions[2].Width = new GridLength(c ? 40 : 330);
        SideTabs.IsVisible = c;
        // the inspector: its column, or a drawer over the right of the view
        InspectorPanel.IsVisible = _vm.InspectorShown;
        Grid.SetColumn(InspectorPanel, c ? 1 : 2);
        InspectorPanel.Width = c ? 300 : double.NaN;
        InspectorPanel.HorizontalAlignment = c ? Avalonia.Layout.HorizontalAlignment.Right : Avalonia.Layout.HorizontalAlignment.Stretch;
        InspectorPanel.ZIndex = c ? 6 : 0;
        InspectorPanel.BoxShadow = c ? Avalonia.Media.BoxShadows.Parse("-16 0 32 0 #70000000") : default;
        // the project panel: a drawer over the left of the view
        Grid.SetColumnSpan(ProjectPanel, c ? 2 : 1);
        ProjectPanel.HorizontalAlignment = c ? Avalonia.Layout.HorizontalAlignment.Left : Avalonia.Layout.HorizontalAlignment.Stretch;
        ProjectPanel.ZIndex = c ? 6 : 0;
        ProjectPanel.BoxShadow = c ? Avalonia.Media.BoxShadows.Parse("16 0 32 0 #70000000") : default;
        // the curves dock: its tabs only until opened
        AnalysisDock.Height = !c || _vm.DockOpen ? 230 : 46;
    }

    private void OnInspectorDrawer(object? s, RoutedEventArgs e) => _vm.InspectorDrawer = !_vm.InspectorDrawer;
    private void OnProjectDrawer(object? s, RoutedEventArgs e) => _vm.ProjectDrawer = !_vm.ProjectDrawer;
    private void OnDockToggle(object? s, RoutedEventArgs e) => _vm.DockOpen = !_vm.DockOpen;

    /// <summary>Lights the tour step's region (the 3D view when that region is hidden).</summary>
    private void UpdateTour()
    {
        if (!_vm.TourOpen) return;
        Dispatcher.UIThread.Post(() =>
        {
            var name = _vm.TourCurrent.Region;
            Control? target = this.FindControl<Control>(name);
            if (target == null || !target.IsEffectivelyVisible || target.Bounds.Width < 2) target = ViewHost;
            var p = target.TranslatePoint(new Point(0, 0), Tour);
            if (p is { } q) Tour.SetTarget(new Rect(q, target.Bounds.Size));
        }, DispatcherPriority.Background);
    }
    private void OnViewportsPage(object? s, RoutedEventArgs e) => _vm.OpenViewports();
    private void OnRenderBack(object? s, RoutedEventArgs e) => _vm.SetModule(8);
    private void OnRenderStop(object? s, RoutedEventArgs e) => _vm.StopRender();

    private async void OnRenderImage(object? s, RoutedEventArgs e)
    {
        if (_vm.Document == null) return;
        var f = await StorageProvider.SaveFilePickerAsync(new FilePickerSaveOptions
        {
            Title = "Render image", DefaultExtension = "png",
            SuggestedFileName = $"{Path.GetFileNameWithoutExtension(_vm.Document.Path)}_frame{_vm.Frame}.png",
            FileTypeChoices = [new FilePickerFileType("PNG") { Patterns = ["*.png"] }],
        });
        if (f?.TryGetLocalPath() is not { } path) return;
        try { await _vm.RenderOut(_ => path, false, FigureDrawing.SaveRenderPng); }
        catch (Exception ex) { _vm.Status = "Render failed: " + ex.Message; }
    }

    private async void OnRenderMovie(object? s, RoutedEventArgs e)
    {
        if (_vm.Document == null) return;
        var dirs = await StorageProvider.OpenFolderPickerAsync(new FolderPickerOpenOptions { Title = "Folder for the frames", AllowMultiple = false });
        if (dirs.Count == 0 || dirs[0].TryGetLocalPath() is not { } dir) return;
        var stem = Path.GetFileNameWithoutExtension(_vm.Document.Path);
        try { await _vm.RenderOut(fr => Path.Combine(dir, $"{stem}_{fr:D5}.png"), true, FigureDrawing.SaveRenderPng); }
        catch (Exception ex) { _vm.Status = "Render failed: " + ex.Message; }
    }

    private void OnFigurePage(object? s, RoutedEventArgs e)
    {
        _vm.SetViewAspect(ViewHost.Bounds.Width, ViewHost.Bounds.Height);
        _vm.OpenFigure();
    }
    // Field › coverage: select the atoms a force field has no type for, assign a force field that covers the structure
    private void OnUntypedSelect(object? s, Avalonia.Interactivity.RoutedEventArgs e)
    {
        if ((s as Control)?.Tag is ViewModels.UntypedGroupRow g) ViewModel.Field.SelectGroup(g);
    }
    private async void OnUseForceField(object? s, Avalonia.Interactivity.RoutedEventArgs e)
    {
        if ((s as Control)?.Tag is string id) await ViewModel.Field.UseForceField(id);
    }
    private async void OnCheckCoverage(object? s, Avalonia.Interactivity.RoutedEventArgs e) => await ViewModel.Field.CheckCoverage();
}
