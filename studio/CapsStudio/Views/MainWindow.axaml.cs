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
    private bool _pan, _addPick;

    private readonly DispatcherTimer _playTimer = new() { Interval = TimeSpan.FromMilliseconds(125) };

    public MainViewModel ViewModel => _vm;
    public void SelectAnalysisTab(int k) => AnalysisTabs.SelectedIndex = k;
    public ICommand OpenCommand { get; }
    public ICommand SaveCommand { get; }
    public bool HasSamples => _samples != null;

    public MainWindow()
    {
        OpenCommand = new RelayCommand(OpenDialog);
        InitializeComponent();
        DataContext = _vm;
        _vm.InitProtocol();
        _vm.LoadReactionSet();
        KeyBindings.Add(new KeyBinding { Gesture = new KeyGesture(Key.O, KeyModifiers.Meta), Command = OpenCommand });
        KeyBindings.Add(new KeyBinding { Gesture = new KeyGesture(Key.O, KeyModifiers.Control), Command = OpenCommand });
        SaveCommand = new RelayCommand(() => _vm.HasDocument && _vm.Idle ? SaveAs("data", "LAMMPS data") : Task.CompletedTask);
        KeyBindings.Add(new KeyBinding { Gesture = new KeyGesture(Key.S, KeyModifiers.Meta), Command = SaveCommand });
        KeyBindings.Add(new KeyBinding { Gesture = new KeyGesture(Key.S, KeyModifiers.Control), Command = SaveCommand });
        _vm.RenderRequested += RequestRender;
        _vm.PropertyChanged += (_, e) =>
        {
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
            if (e.PropertyName == nameof(MainViewModel.ChainCurve)) ChainPlot.SetData(_vm.ChainCurve);
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
        _playTimer.Tick += (_, _) => _vm.StepFrame(1);
        AddHandler(KeyDownEvent, OnKey, RoutingStrategies.Tunnel);
        DragDrop.SetAllowDrop(ViewHost, true);
        ViewHost.AddHandler(DragDrop.DropEvent, OnDrop);
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

    private async Task OpenDialog()
    {
        var files = await StorageProvider.OpenFilePickerAsync(new FilePickerOpenOptions
        {
            Title = "Open structure or trajectory",
            AllowMultiple = true,
            FileTypeFilter =
            [
                new FilePickerFileType("Structures and trajectories") { Patterns = ["*.data", "*.lmp", "*.lammpstrj", "*.dump", "*.gro", "*.pdb", "*.ent", "*.xyz", "*.extxyz"] },
                new FilePickerFileType("All files") { Patterns = ["*"] },
            ],
        });
        OpenMany(files.Select(f => f.TryGetLocalPath()).OfType<string>().ToList());
    }

    private void OpenMany(List<string> paths)
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

    private void OnOpenSample(object? sender, RoutedEventArgs e)
    {
        if (_samples != null) TryOpen(Path.Combine(_samples, "ps_melt.lammpstrj"), Path.Combine(_samples, "ps_melt.data"));
    }

    private void OnPlay(object? s, RoutedEventArgs e) => TogglePlay();

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
    private void OnThemeDark(object? s, RoutedEventArgs e) { Tokens.Use(false); RequestRender(); }
    private void OnThemeLight(object? s, RoutedEventArgs e) { Tokens.Use(true); RequestRender(); }
    private void OnCommandPalette(object? s, RoutedEventArgs e) => _vm.Status = "Command palette: coming with the command layer (design board CommandPalette)";

    // ---- Analyze › Properties
    private async void OnAnalyzeRun(object? s, RoutedEventArgs e) { if (_vm.Idle) await _vm.Analyze.Run(); }
    private void OnAnalyzeCancel(object? s, RoutedEventArgs e) => _vm.Analyze.Cancel();

    private async void OnAnalyzeExport(object? s, RoutedEventArgs e)
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
    private async void OnFieldOverride(object? s, RoutedEventArgs e) { if (!_vm.Busy) await _vm.Field.ApplyOverride(); }
    private async void OnFieldResetOverride(object? s, RoutedEventArgs e) { if (!_vm.Busy) await _vm.Field.ResetOverride(); }
    private async void OnFieldAddRule(object? s, RoutedEventArgs e) { if (!_vm.Busy) await _vm.Field.AddRule(); }
    private async void OnFieldRemoveRules(object? s, RoutedEventArgs e) { if (!_vm.Busy) await _vm.Field.RemoveRules(); }

    private async void OnFieldImport(object? s, RoutedEventArgs e)
    {
        if (_vm.Busy) return;
        var files = await StorageProvider.OpenFilePickerAsync(new FilePickerOpenOptions
        {
            Title = "Import force-field parameters",
            AllowMultiple = false,
            FileTypeFilter = [new FilePickerFileType("CAPS force field or moltemplate") { Patterns = ["*.json", "*.lt"] }],
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

    private async void OnSaveTrajectory(object? s, RoutedEventArgs e)
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

    private void TogglePlay()
    {
        if (_vm.Document == null || !_vm.HasFrames) return;
        _vm.IsPlaying = !_vm.IsPlaying;
        if (_vm.IsPlaying) _playTimer.Start(); else _playTimer.Stop();
    }

    // Viewer keys; ignored while typing in a text box or choosing in a list.
    private void OnKey(object? sender, KeyEventArgs e)
    {
        if (_vm.Document == null || _vm.Busy || FocusManager?.GetFocusedElement() is TextBox or ComboBox) return;
        switch (e.Key)
        {
            case Key.Left: _vm.StepFrame(-1); e.Handled = true; break;
            case Key.Right: _vm.StepFrame(1); e.Handled = true; break;
            case Key.Space: TogglePlay(); e.Handled = true; break;
            case Key.R when e.KeyModifiers == KeyModifiers.None: _vm.ResetView(); e.Handled = true; break;
            case Key.Escape: _vm.ClearSelection(); e.Handled = true; break;
        }
    }

    private void OnResetView(object? s, RoutedEventArgs e) => _vm.ResetView();
    private void OnViewTop(object? s, RoutedEventArgs e) => _vm.SetView(0, Math.PI / 2);
    private void OnViewFront(object? s, RoutedEventArgs e) => _vm.SetView(0, 0);
    private void OnViewSide(object? s, RoutedEventArgs e) => _vm.SetView(Math.PI / 2, 0);

    // ---------------------------------------------------------------- rendering

    private void RequestRender()
    {
        _requested++;
        if (!_busy) _ = RenderLoop();
    }

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
                if (_vm.Busy) { _rendered = ticket; break; }   // the core is busy with this document; keep the last image
                var pw = (int)(w * _scaling);
                var ph = (int)(h * _scaling);
                var cam = _vm.Camera;
                // Retina already gives 2 samples per point; supersample only on 1× screens.
                var opt = field ? _vm.FieldViewOptions(pw, ph, _scaling >= 1.5 ? 1 : 2) : _vm.ViewOptions(pw, ph, _scaling >= 1.5 ? 1 : 2);
                var buf = new byte[pw * ph * 4];
                var sw = Stopwatch.StartNew();
                await Task.Run(() => doc.Render(cam, opt, buf));
                sw.Stop();
                var bmp = new WriteableBitmap(new PixelSize(pw, ph), new Vector(96 * _scaling, 96 * _scaling), PixelFormat.Rgba8888, AlphaFormat.Unpremul);
                using (var fb = bmp.Lock())
                {
                    for (var y = 0; y < ph; y++)
                        System.Runtime.InteropServices.Marshal.Copy(buf, y * pw * 4, fb.Address + y * fb.RowBytes, pw * 4);
                }
                var old = image.Source as IDisposable;
                image.Source = bmp;
                old?.Dispose();
                _pixW = pw; _pixH = ph;
                RenderStat.Text = $"{pw}×{ph} px · {sw.ElapsedMilliseconds} ms";
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

    // ---------------------------------------------------------------- mouse

    private void OnPointerPressed(object? sender, PointerPressedEventArgs e)
    {
        _hostField = sender as Control ?? ViewHost;
        var p = e.GetCurrentPoint(_host);
        _press = _last = p.Position;
        _dragging = true;
        _moved = false;
        _pan = p.Properties.IsRightButtonPressed || e.KeyModifiers.HasFlag(KeyModifiers.Alt);
        _addPick = e.KeyModifiers.HasFlag(KeyModifiers.Shift) || e.KeyModifiers.HasFlag(KeyModifiers.Meta) || _vm.MeasureTool;
        e.Pointer.Capture(_host);
        _host.Focus();
    }

    private void OnPointerMoved(object? sender, PointerEventArgs e)
    {
        if (!_dragging || _vm.Document == null || _vm.Busy) return;
        var pos = e.GetPosition(_host);
        var d = pos - _last;
        _last = pos;
        if (Math.Abs(pos.X - _press.X) + Math.Abs(pos.Y - _press.Y) > 3) _moved = true;
        if (!_moved) return;
        if (_pan)
        {
            // Pan in Å: approximate px→Å from the current fit (view height ≈ 2.2 × half-extent at zoom 1).
            var s = Summary();
            var span = Math.Max(1.0, Math.Max(s.CellA, Math.Max(s.CellB, s.CellC)));
            var perPx = span * 1.6 / Math.Max(1, _host.Bounds.Height) / _vm.Camera.Zoom;
            _vm.Camera.PanX += d.X * perPx;
            _vm.Camera.PanY -= d.Y * perPx;
        }
        else
        {
            _vm.Camera.Yaw += d.X * 0.008;
            _vm.Camera.Pitch = Math.Clamp(_vm.Camera.Pitch + d.Y * 0.008, -Math.PI / 2, Math.PI / 2);
        }
        RequestRender();
    }

    private Interop.CapsSummary Summary() => _vm.Document?.Summary() ?? default;

    private void OnPointerReleased(object? sender, PointerReleasedEventArgs e)
    {
        e.Pointer.Capture(null);
        if (_dragging && !_moved && _vm.Document != null && !_vm.Busy)
        {
            var pos = e.GetPosition(_host);
            var hit = _vm.Document.Pick((int)(pos.X * _scaling), (int)(pos.Y * _scaling));
            if (_host == FieldViewHost) { if (hit >= 0) _vm.Field.SelectAtom(hit); }
            else _vm.Pick(hit, _addPick);
            RequestRender();
        }
        _dragging = false;
    }

    private void OnWheel(object? sender, PointerWheelEventArgs e)
    {
        if (_vm.Document == null || _vm.Busy) return;
        _vm.Camera.Zoom = Math.Clamp(_vm.Camera.Zoom * Math.Pow(1.12, e.Delta.Y), 0.1, 40);
        RequestRender();
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

    private async void OnExportPng(object? s, RoutedEventArgs e) => await Export("png");
    private async void OnExportSvg(object? s, RoutedEventArgs e) => await Export("svg");
}
