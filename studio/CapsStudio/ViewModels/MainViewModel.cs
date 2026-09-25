using System.Collections.ObjectModel;
using System.ComponentModel;
using System.Globalization;
using System.Runtime.CompilerServices;
using CapsStudio.Interop;

namespace CapsStudio.ViewModels;

public abstract class ObservableObject : INotifyPropertyChanged
{
    public event PropertyChangedEventHandler? PropertyChanged;

    protected bool Set<T>(ref T field, T value, [CallerMemberName] string? name = null)
    {
        if (EqualityComparer<T>.Default.Equals(field, value)) return false;
        field = value;
        PropertyChanged?.Invoke(this, new PropertyChangedEventArgs(name));
        OnChanged(name);
        return true;
    }

    protected virtual void OnChanged(string? name) { }

    protected void Raise([CallerMemberName] string? name = null) => PropertyChanged?.Invoke(this, new PropertyChangedEventArgs(name));
}

public sealed record Row(string Key, string Value);
public sealed record MoleculeRow(string Molecule, string Atoms, string Mass, string Rg, string Kappa2);
/// <summary>A structure block of a Pack input: name, where it goes, how many, its colour and file.</summary>
/// <summary>A pre-flight check: what, and ok / check / fail.</summary>
public sealed record CheckRow(string Title, string State)
{
    public string Icon => State == "ok" ? "check" : State == "fail" ? "xcircle" : "alert";
    public Avalonia.Media.IBrush Brush => CapsStudio.Tokens.Brush(State == "ok" ? "OkB" : State == "fail" ? "ErrB" : "WarnB");
}
/// <summary>A convergence criterion: what, the rule, the latest value, pass / not yet.</summary>
public sealed record CriterionRow(string Title, string Rule, string Now, string State)
{
    public string Icon => State == "pass" ? "check" : "alert";
    public Avalonia.Media.IBrush Brush => CapsStudio.Tokens.Brush(State == "pass" ? "OkB" : "WarnB");
}
public sealed record PackItem(string Name, string Detail, string Count, string Colour, string File)
{
    public Avalonia.Media.IBrush Brush => Avalonia.Media.Brush.Parse(Colour);
}

public sealed partial class MainViewModel : ObservableObject
{
    public static readonly string[] Styles = ["Ball & stick", "Space filling", "Sticks", "No hydrogens", "Backbone"];
    public static readonly string[] ColourModes = ["Element", "Molecule", "Type", "Distance to molecule centre"];
    public static readonly string[] Backgrounds = ["Dark", "White", "Transparent"];
    public static readonly string[] ViewBackgrounds = ["Dark", "White (paper)"];
    public static readonly string[] RdfPairs = ["C – C", "C – H", "H – H", "all – all"];
    private static readonly (int A, int B)[] RdfElements = [(6, 6), (6, 1), (1, 1), (0, 0)];
    public static readonly string[] SizePresets = ["1920 × 1080 (slide)", "2008 × 1130 (85 mm @ 600 dpi)", "4016 × 2259 (170 mm @ 600 dpi)", "Viewport size"];

    private CapsDocument? _doc;
    private string _title = "No file open";
    private string _status = "Open a LAMMPS data or dump file, a GROMACS .gro, a PDB, XYZ, mol2 or CIF file";
    private int _style;
    private int _colour = 1;
    private bool _outlines = true, _depthCue = true, _showCell = true, _perspective;
    private int _exportBackground = 2;
    private int _sizePreset;
    private int _frame;
    private int _frames = 1;
    private readonly List<int> _selection = new();
    private string _pickedTitle = "Nothing picked";
    private string _measure = "";
    private int _rdfPair;
    private bool _rdfInter = true;
    private bool _showRdf = true;
    private (double R, double G)[] _rdf = [];
    private string _rdfNote = "";

    public CapsCamera Camera = new() { Yaw = 0.55, Pitch = 0.40, Zoom = 1.0 };

    /// <summary>CAPS Field: the force-field assignment of the open structure.</summary>
    public FieldViewModel Field { get; }
    /// <summary>Analyze › Properties: calculations over the frames of the open trajectory.</summary>
    public AnalyzeViewModel Analyze { get; }

    public MainViewModel()
    {
        Field = new FieldViewModel(() => _doc, s => Status = s, () =>
        {
            // types and charges in the document changed: summary, inspector and viewer follow
            RefreshSummary();
            RefreshSelection();
            RefreshLegend();
            FieldInfoText = "";
            RenderRequested?.Invoke();
        });
        Analyze = new AnalyzeViewModel(() => _doc, s => Status = s, running => { _analyzing = running; RaiseBusy(); });
        Field.PropertyChanged += (_, e) => { if (e.PropertyName == nameof(FieldViewModel.RunLine)) Raise(nameof(ForceFieldLine)); };
    }

    /// <summary>Which force field Relax, Dynamics and Equilibrate will use.</summary>
    public string ForceFieldLine => Field.RunLine;

    public event Action? RenderRequested;

    public CapsDocument? Document
    {
        get => _doc;
        private set
        {
            if (Set(ref _doc, value))
            {
                // styles, surfaces and labels belong to a document: a new one starts plain
                AppLayers.Clear(); _appColour = -1; _appSurface = 0; _labelTexts = null;
                NamedSets.Clear(); _selCount = 0; Raise(nameof(SelectedCount)); Raise(nameof(SelectedChip)); Dyads.Clear();
                foreach (var n in new[] { nameof(AppColour), nameof(AppSurface), nameof(AppHasSurface), nameof(AppChip), nameof(ShowAppLegend) }) Raise(n);
                RaiseAppearanceVisibility();
                Raise(nameof(HasDocument));
                Raise(nameof(NoDocument));
                Raise(nameof(ShowEmpty));
                Raise(nameof(CanRelax));
                Raise(nameof(CanRun));
                Raise(nameof(CanEquilibrate));
            }
        }
    }

    public bool HasDocument => _doc != null;
    public bool NoDocument => _doc == null;

    public string Title { get => _title; set => Set(ref _title, value); }
    public string Status { get => _status; set => Set(ref _status, value); }

    public int StyleIndex { get => _style; set { if (Set(ref _style, value)) { Raise(nameof(StyleText)); RenderRequested?.Invoke(); } } }
    // toolbar texts (Main board: "Ball & stick", "Colour: element", "Perspective")
    public string StyleText => Styles[Math.Clamp(_style, 0, Styles.Length - 1)];
    public string ColourText => "Colour: " + (_appColour == 4 ? "partial charge" : ColourModes[Math.Clamp(ColourIndex, 0, ColourModes.Length - 1)].ToLowerInvariant());
    public string ProjectionText => _perspective ? "Perspective" : "Orthographic";
    private bool _measureTool;
    /// <summary>Measure tool: clicks add atoms to the measurement (as ⇧ click does).</summary>
    public bool MeasureTool { get => _measureTool; set => Set(ref _measureTool, value); }
    public int ColourIndex
    {
        get => _colour;
        set { if (Set(ref _colour, value)) { Raise(nameof(ColourText)); Raise(nameof(OvLegendNote)); Raise(nameof(FigLegendAvailable)); RefreshLegend(); RenderRequested?.Invoke(); } }
    }

    private string _legendLo = "", _legendHi = "";
    public bool ShowLegend => _doc != null && _colour == 3 && _module != 19;
    public string LegendLo { get => _legendLo; private set => Set(ref _legendLo, value); }
    public string LegendHi { get => _legendHi; private set => Set(ref _legendHi, value); }

    private void RefreshLegend()
    {
        Raise(nameof(ShowLegend));
        if (_doc == null) return;
        var (lo, hi) = _doc.PropertyRange();
        LegendLo = lo.ToString("0.0", CultureInfo.InvariantCulture) + " Å";
        LegendHi = hi.ToString("0.0", CultureInfo.InvariantCulture) + " Å";
    }

    private bool _playing;
    public bool IsPlaying { get => _playing; set { if (Set(ref _playing, value)) { Raise(nameof(PlayLabel)); Raise(nameof(PlayIcon)); } } }
    public string PlayLabel => _playing ? "Pause" : "Play";

    public void StepFrame(int d)
    {
        if (_doc == null || _frames < 2) return;
        Frame = ((_frame + d) % _frames + _frames) % _frames;
    }

    public void ClearSelection()
    {
        _selection.Clear();
        RefreshSelection();
        RenderRequested?.Invoke();
    }
    public bool Outlines { get => _outlines; set { if (Set(ref _outlines, value)) RenderRequested?.Invoke(); } }
    public bool DepthCue { get => _depthCue; set { if (Set(ref _depthCue, value)) RenderRequested?.Invoke(); } }
    public bool ShowCell { get => _showCell; set { if (Set(ref _showCell, value)) RenderRequested?.Invoke(); } }

    private bool _wrap;
    public bool Wrap
    {
        get => _wrap;
        set
        {
            if (!Set(ref _wrap, value) || _doc == null) return;
            _doc.SetWrap(value);
            RefreshSelection();
            Status = value ? "Atoms folded into the cell · bonds that cross a face are hidden" : "Molecules shown whole";
            RenderRequested?.Invoke();
        }
    }
    public bool Perspective
    {
        get => _perspective;
        set { if (Set(ref _perspective, value)) { Raise(nameof(ProjectionText)); Camera.Perspective = value ? 1 : 0; RenderRequested?.Invoke(); } }
    }

    public int ExportBackground { get => _exportBackground; set => Set(ref _exportBackground, value); }

    private int _viewBackground;
    public int ViewBackground
    {
        get => _viewBackground;
        set { if (Set(ref _viewBackground, value)) { Raise(nameof(ViewIsLight)); RenderRequested?.Invoke(); } }
    }
    public bool ViewIsLight => _viewBackground == 1;
    public int SizePreset { get => _sizePreset; set => Set(ref _sizePreset, value); }

    public int Frames { get => _frames; private set { if (Set(ref _frames, value)) Raise(nameof(HasFrames)); } }
    public bool HasFrames => _frames > 1;
    public int FrameMax => Math.Max(0, _frames - 1);

    public int Frame
    {
        get => _frame;
        set
        {
            if (_doc == null || !Set(ref _frame, Math.Clamp(value, 0, FrameMax))) return;
            _doc.SetFrame(_frame);
            Raise(nameof(FrameLabel));
            if (_module == 32) { TrajUpdateValues(); TrajectoryChanged?.Invoke(); }
            if (_ixOpen && _ixLive) RunInteractions();
            RefreshSummary();
            RefreshSelection();
            RefreshRdf();
            RefreshMolecules();
            RefreshLegend();
            FocusOnFrame();
            if (IsVisualize) RefreshPipeline();
            if (IsViewports) RenderViewports();
            RenderRequested?.Invoke();
        }
    }

    public string FrameLabel => $"frame {_frame} / {FrameMax}";

    public int Picked => _selection.Count > 0 ? _selection[^1] : -1;
    public bool HasPicked => _selection.Count > 0;
    public IReadOnlyList<int> Selection => _selection;
    public string MeasureText { get => _measure; private set { if (Set(ref _measure, value)) Raise(nameof(HasMeasure)); } }
    // Atom inspector and viewport HUD
    private string _pElem = "", _pX = "", _pY = "", _pZ = "", _pType = "", _pCharge = "", _pMol = "", _pWhy = "", _pFf = "";
    public string PickedElement { get => _pElem; private set => Set(ref _pElem, value); }
    public string PickedX { get => _pX; private set => Set(ref _pX, value); }
    public string PickedY { get => _pY; private set => Set(ref _pY, value); }
    public string PickedZ { get => _pZ; private set => Set(ref _pZ, value); }
    public string PickedType { get => _pType; private set => Set(ref _pType, value); }
    public string PickedCharge { get => _pCharge; private set => Set(ref _pCharge, value); }
    public string PickedMolecule { get => _pMol; private set => Set(ref _pMol, value); }
    public string PickedWhy { get => _pWhy; private set => Set(ref _pWhy, value); }
    public string PickedFf { get => _pFf; private set => Set(ref _pFf, value); }
    public bool HasAtom => _doc != null && Picked >= 0;
    public string HudSelection => _selection.Count == 1 ? "1 atom selected" : $"{_selection.Count} atoms selected";
    public bool HasSelection => _selection.Count > 0;
    public bool HasMeasure => _measure.Length > 0;

    public int RdfPairIndex { get => _rdfPair; set { if (Set(ref _rdfPair, value)) RefreshRdf(); } }
    public bool RdfInter { get => _rdfInter; set { if (Set(ref _rdfInter, value)) RefreshRdf(); } }
    public bool ShowRdf { get => _showRdf; set => Set(ref _showRdf, value); }
    public (double R, double G)[] RdfCurve { get => _rdf; private set => Set(ref _rdf, value); }
    public string RdfNote { get => _rdfNote; private set => Set(ref _rdfNote, value); }
    public string PickedTitle { get => _pickedTitle; set => Set(ref _pickedTitle, value); }

    public ObservableCollection<Row> SummaryRows { get; } = new();
    public ObservableCollection<string> Notes { get; } = new();
    public ObservableCollection<Row> PickedRows { get; } = new();
    public ObservableCollection<Row> NeighbourRows { get; } = new();
    public ObservableCollection<MoleculeRow> MoleculeRows { get; } = new();
    private string _moleculeNote = "";
    public string MoleculeNote { get => _moleculeNote; private set => Set(ref _moleculeNote, value); }

    private (double X, double Y)[] _chainCurve = [];
    private string _chainNote = "";
    public (double X, double Y)[] ChainCurve { get => _chainCurve; private set => Set(ref _chainCurve, value); }
    public string ChainNote { get => _chainNote; private set => Set(ref _chainNote, value); }

    private void RefreshChains()
    {
        if (_doc == null) { ChainCurve = []; ChainNote = ""; return; }
        try
        {
            var (n, r, chains, b2) = _doc.InternalDistances();
            ChainCurve = n.Select((v, k) => ((double)v, r[k])).ToArray();
            ChainNote = chains == 0 ? "no chains of four or more heavy atoms"
                : string.Format(CultureInfo.InvariantCulture, "{0} backbones · ⟨b²⟩ {1:F3} Å² · plateau → C∞ for equilibrated long chains · frame {2}", chains, b2, _frame);
        }
        catch (Exception e) { ChainCurve = []; ChainNote = e.Message; }
    }

    private void RefreshMolecules()
    {
        RefreshChains();
        MoleculeRows.Clear();
        if (_doc == null) return;
        var inv = CultureInfo.InvariantCulture;
        var m = _doc.Molecules();
        foreach (var x in m.Take(500))
            MoleculeRows.Add(new(x.Molecule.ToString(inv), x.Atoms.ToString(inv), x.Mass.ToString("F1", inv), x.Rg.ToString("F3", inv), x.Kappa2.ToString("F3", inv)));
        if (m.Length > 0)
            MoleculeNote = string.Format(inv, "{0} molecules · mean Rg {1:F3} Å · mean κ² {2:F3}{3}", m.Length, m.Average(x => x.Rg), m.Average(x => x.Kappa2),
                m.Length > 500 ? " · first 500 listed" : "");
    }

    public void Open(string path, string? topology = null)
    {
        if (OpensProgressively(path)) { _ = OpenProgressive(path, topology); return; }
        Show(CapsDocument.Open(path, topology), System.IO.Path.GetFileName(path));
        if (_doc?.Path == path) Remember(path, topology);
    }

    // ---------------------------------------------------------------- modules
    private int _module = 8;   // 0 Grow, 1 Analyze, 2 Relax, 3 Dynamics, 4 Equilibrate, 5 Pack, 6 React, 7 Field, 8 Studio, 9 Molecule, 10 Settings, 11 Jobs, 12 Bench, 13 Polymer
    public bool IsGrow => _module == 0;
    public bool IsAnalyze => _module == 1;
    public bool IsRelax => _module == 2;
    public bool IsDynamics => _module == 3;
    public bool IsEquilibrate => _module == 4;
    public bool IsPack => _module == 5;
    public bool IsReact => _module == 6;
    public bool IsField => _module == 7;
    /// <summary>Studio: the workspace with the 3D view and the inspector.</summary>
    public bool IsStudio => _module == 8;
    private static readonly string[] Crumbs = ["Grow › Amorphous cell", "Analyze › Properties", "Relax › Minimise", "Dynamics › Run",
        "Equilibrate › Protocol", "Pack › Molecules & regions", "React › Crosslinking", "Field › Typing report", "Studio", "Studio › Molecule", "Settings", "Jobs", "Bench", "Builders › Polymer", "Builders › Surface", "Builders › Nanostructure", "Builders › Polymer › Blend", "Studio › File checks", "Export › Figure", "Studio › Render", "Analyze › Visualize", "Export › Data", "Analyze › Batch", "Analyze › Compare", "Analyze › Visualize › Colour by", "Studio › Viewports", "Export › Figure bundle", "Open file", "Analyze › Visualize › Save pipeline", "Builders › Crystal", "Builders › Biomolecule", "Builders › Solvation", "Studio › Trajectory", "Studio › Torsion scan"];
    /// <summary>Where the user is (top bar).</summary>
    public string Crumb => _module == 8 ? "" : Crumbs[_module];
    /// <summary>Where calculations run (top bar).</summary>
    public string ComputeText { get; } = $"Local · {Environment.ProcessorCount} threads · CPU";
    /// <summary>Kept for scripts and tests: the property calculations are the Analyze module.</summary>
    public bool AnalyzeProperties { get => _module == 1; set { if (value) SetModule(1); else if (_module == 1) SetModule(8); } }
    public bool IsProperties => _module == 1;
    public void SetModule(int m)
    {
        var was = _module;
        if (!Set(ref _module, m, nameof(Module))) return;
        if (was is 20 or 21 && m != 20) SuspendPipeline();
        if (m == 20) ApplyPipeline();
        Raise(nameof(IsGrow));
        Raise(nameof(IsAnalyze));
        Raise(nameof(IsRelax));
        Raise(nameof(IsDynamics));
        Raise(nameof(IsEquilibrate));
        Raise(nameof(IsPack));
        Raise(nameof(IsReact));
        Raise(nameof(IsField));
        Raise(nameof(IsStudio));
        Raise(nameof(IsMolecule));
        Raise(nameof(IsStudioRail));
        Raise(nameof(IsSettings));
        Raise(nameof(IsJobs));
        Raise(nameof(IsBench));
        Raise(nameof(IsPolymer));
        Raise(nameof(IsSurface));
        Raise(nameof(IsNano));
        Raise(nameof(IsBlend));
        Raise(nameof(ShowEmpty));
        Raise(nameof(IsChecks));
        Raise(nameof(IsFigure));
        Raise(nameof(IsRender));
        Raise(nameof(IsVisualize));
        Raise(nameof(HasTimeline));
        Raise(nameof(IsExport));
        Raise(nameof(IsBatch));
        Raise(nameof(IsCompare));
        Raise(nameof(IsColourBy));
        Raise(nameof(IsViewports));
        Raise(nameof(IsBundle));
        Raise(nameof(IsOpenPreview));
        Raise(nameof(IsSavePipeline));
        Raise(nameof(IsCrystal));
        Raise(nameof(IsBio));
        Raise(nameof(IsSolvation));
        Raise(nameof(IsTrajectory));
        Raise(nameof(IsTorsion));
        RaiseAppearanceVisibility();
        Raise(nameof(IsAnalyzeRail));
        Raise(nameof(ShowPipeLegend));
        Raise(nameof(ShowAnalysisPanel));
        Raise(nameof(ShowLegend));
        Raise(nameof(Crumb));
        Raise(nameof(IsProperties));
        RenderRequested?.Invoke();   // the Field page has its own view
    }
    public int Module => _module;

    // ---------------------------------------------------------------- Grow
    public static readonly string[] Tacticities = ["Atactic", "Isotactic", "Syndiotactic"];
    private int _growChains = 10, _growDp = 8, _growTact, _growSeed = 1;
    private double _growDensity = 0.4, _growBox = 0, _growScale = 1.0;
    private bool _growUseBox, _growCurve = true, _growing;
    private string _growLog = "Choose chains and density, then Build. Growth checks every contact against all atoms already in the cell, across the periodic faces.";
    private CancellationTokenSource? _growCancel;
    private bool _grownUnsaved;

    public int GrowChains { get => _growChains; set => Set(ref _growChains, Math.Clamp(value, 1, 2000)); }
    public int GrowDp { get => _growDp; set => Set(ref _growDp, Math.Clamp(value, 2, 2000)); }
    public int GrowTacticity { get => _growTact; set => Set(ref _growTact, value); }
    public int GrowSeed { get => _growSeed; set => Set(ref _growSeed, Math.Max(0, value)); }
    public double GrowDensity { get => _growDensity; set => Set(ref _growDensity, Math.Clamp(value, 0.01, 2.0)); }
    public double GrowBox { get => _growBox; set => Set(ref _growBox, Math.Max(0, value)); }
    public double GrowScale { get => _growScale; set => Set(ref _growScale, Math.Clamp(value, 0.5, 1.2)); }
    public bool GrowUseBox { get => _growUseBox; set { if (Set(ref _growUseBox, value)) Raise(nameof(GrowUseDensity)); } }
    public bool GrowUseDensity => !_growUseBox;
    public bool GrowCurve { get => _growCurve; set => Set(ref _growCurve, value); }
    public bool Growing { get => _growing; private set { if (Set(ref _growing, value)) { Raise(nameof(NotGrowing)); RaiseBusy(); } } }
    public bool NotGrowing => !_growing;
    /// <summary>No build or minimisation is running (the document can be replaced or edited).</summary>
    public bool Idle => !_growing && !_relaxing && !_mdRunning && !_eqRunning && !_packing && !_reacting && !_analyzing;
    /// <summary>The core is working on the open document (Relax or Dynamics): no rendering or edits until it is done.</summary>
    public bool Busy => _relaxing || _mdRunning || _eqRunning || _reacting || _analyzing;
    private bool _analyzing;
    private void RaiseBusy()
    {
        Raise(nameof(Idle));
        Raise(nameof(Busy));
        Raise(nameof(CanRelax));
        Raise(nameof(CanRun));
        Raise(nameof(CanEquilibrate));
        Raise(nameof(CanReact));
    }
    public string GrowLog { get => _growLog; private set => Set(ref _growLog, value); }
    // live statistics (Grow board: chains grown, restarts, elapsed, overall progress)
    private int _growDone, _growRestarts;
    private double _growElapsed;
    public int GrowDone { get => _growDone; private set { if (Set(ref _growDone, value)) { Raise(nameof(GrowProgress)); Raise(nameof(GrowDoneText)); } } }
    public int GrowRestarts { get => _growRestarts; private set => Set(ref _growRestarts, value); }
    public double GrowElapsed { get => _growElapsed; private set { if (Set(ref _growElapsed, value)) Raise(nameof(GrowElapsedText)); } }
    public double GrowProgress => _growChains > 0 ? (double)_growDone / _growChains : 0;
    public string GrowDoneText => $"{_growDone} / {_growChains}";
    public string GrowElapsedText => string.Format(CultureInfo.InvariantCulture, "{0:F1} s", _growElapsed);
    public string GrowTacticityName => Tacticities[Math.Clamp(_growTact, 0, 2)].ToLowerInvariant();
    public string GrowSizeLabel => _growUseBox ? "Box edge (Å)" : "Target density (g/cm³)";
    public string GrowAtomsText => (GrowVaries ? "≈ " : "") + (_growChains * GrowChainSize().Atoms).ToString("N0", CultureInfo.InvariantCulture);
    /// <summary>Random and gradient copolymers draw each chain's sequence, so sizes vary from chain to chain.</summary>
    private bool GrowVaries => _growSpec != null && (_growSpec.Contains("\"random\"") || _growSpec.Contains("\"gradient\""));
    /// <summary>The same build from the command line.</summary>
    public string GrowCommand => string.Format(CultureInfo.InvariantCulture,
        "caps grow -o {7}_{0}x{1}.data --chains {0} --dp {1} {2} --tacticity {3} --seed {4} --scale {5}{6}{8}",
        _growChains, _growDp, _growUseBox ? $"--box {_growBox}" : $"--density {_growDensity}", Tacticities[_growTact].ToLowerInvariant(), _growSeed, _growScale,
        _growCurve ? "" : " --trans", _growSpec == null ? "PS" : "polymer", GrowUnitsArg());
    private string GrowUnitsArg()
    {
        if (_growSpec == null) return "";
        var j = System.Text.Json.Nodes.JsonNode.Parse(_growSpec)!;
        var units = string.Join(",", (j["units"] as System.Text.Json.Nodes.JsonArray ?? []).Select(u => (string?)u!["smiles"]));
        var seq = (string?)j["sequence"] ?? "homopolymer";
        var extra = seq switch
        {
            "random" => " --weights " + string.Join(",", (j["weights"] as System.Text.Json.Nodes.JsonArray ?? []).Select(w => ((double?)w ?? 1).ToString(CultureInfo.InvariantCulture))),
            "block" => " --blocks " + string.Join(",", (j["blocks"] as System.Text.Json.Nodes.JsonArray ?? []).Select(w => (int?)w ?? 1)),
            "pattern" => " --pattern " + (string?)j["pattern"],
            _ => "",
        };
        return $" --units '{units}' --sequence {seq}{extra}";
    }
    public bool GrownUnsaved { get => _grownUnsaved; private set => Set(ref _grownUnsaved, value); }

    // decimal views for NumericUpDown
    public decimal? GrowChainsD { get => _growChains; set { GrowChains = (int)(value ?? 1); Raise(); } }
    public decimal? GrowDpD { get => _growDp; set { GrowDp = (int)(value ?? 2); Raise(); } }
    public decimal? GrowSeedD { get => _growSeed; set { GrowSeed = (int)(value ?? 0); Raise(); } }
    public decimal? GrowDensityD { get => (decimal)_growDensity; set { GrowDensity = (double)(value ?? 0.4m); Raise(); } }
    public decimal? GrowBoxD { get => (decimal)_growBox; set { GrowBox = (double)(value ?? 0m); Raise(); } }
    public decimal? GrowScaleD { get => (decimal)_growScale; set { GrowScale = (double)(value ?? 1m); Raise(); } }
    private bool _growAutoScale = true;
    /// <summary>Repeat-unit chains: lower the contact scale (0.85 … 0.6) when a chain cannot be placed.</summary>
    public bool GrowAutoScale { get => _growAutoScale; set => Set(ref _growAutoScale, value); }

    public string GrowEstimate
    {
        get
        {
            var inv = CultureInfo.InvariantCulture;
            var (perAtoms, perMass) = GrowChainSize();
            var atoms = _growChains * perAtoms;
            var mass = _growChains * perMass;
            var box = _growUseBox ? _growBox : Math.Cbrt(mass / 6.02214076e23 / _growDensity) * 1e8;
            var rho = _growUseBox && _growBox > 0 ? mass / 6.02214076e23 / Math.Pow(_growBox * 1e-8, 3) : _growDensity;
            return string.Format(inv, "{0:N0} atoms · {1:N0} g/mol per chain · box {2:F2} Å · {3:F3} g/cm³", atoms, mass / _growChains, box, rho);
        }
    }

    protected override void OnChanged(string? name)
    {
        if (name != null && name.StartsWith("Grow") && name != nameof(GrowEstimate) && name != nameof(GrowLog) && name != nameof(GrowCommand) && name != nameof(GrowAtomsText)
            && name != nameof(GrowDoneText) && name != nameof(GrowProgress))
        {
            Raise(nameof(GrowEstimate));
            Raise(nameof(GrowCommand));
            Raise(nameof(GrowAtomsText));
            Raise(nameof(GrowDoneText));
            Raise(nameof(GrowTacticityName));
            Raise(nameof(GrowComponentName));
            Raise(nameof(GrowSizeLabel));
        }
    }

    public async Task Grow()
    {
        if (_growing) return;
        Growing = true;
        _growCancel = new CancellationTokenSource();
        var token = _growCancel.Token;
        var o = new CapsGrowOpts
        {
            Chains = _growChains, Dp = _growDp, Tacticity = _growTact, Seed = (ulong)_growSeed,
            Box = _growUseBox ? _growBox : 0, Density = _growUseBox ? 0 : _growDensity, ContactScale = _growScale, Curve = _growCurve ? 1 : 0,
        };
        var spec = _growSpec;
        if (spec != null && _growAutoScale) o.ContactScale = -_growScale;   // caps_grow_chains: start here, lower it when crowded
        if (spec != null)
        {
            var sj = System.Text.Json.Nodes.JsonNode.Parse(spec)!.AsObject();
            sj["dp"] = _growDp;
            spec = sj.ToJsonString();
        }
        var stem = spec == null ? "PS" : string.Concat(_growSpecName.Where(char.IsLetterOrDigit).Take(16));
        var label = $"{stem}_{_growChains}x{_growDp}_{Tacticities[_growTact].ToLowerInvariant()}_seed{_growSeed}";
        GrowLog = "Growing…";
        GrowDone = 0;
        GrowRestarts = 0;
        GrowElapsed = 0;
        Status = $"Growing {_growChains} chains of {(spec == null ? "polystyrene" : _growSpecName)}, DP {_growDp}…";
        var sw = System.Diagnostics.Stopwatch.StartNew();
        try
        {
            var lastUi = 0L;
            (CapsDocument Doc, string Report)? result = null;
            var seed = _growSeed;
            var failures = new List<string>();
            // A seed that jams is not a property of the request: try up to three consecutive seeds, and say so.
            for (var attempt = 0; attempt < 3 && result == null; attempt++, seed++)
            {
                var oa = o;
                oa.Seed = (ulong)seed;
                var s0 = seed;
                try
                {
                    Func<int, int, int, bool> onProgress = (d, t, r) =>
                    {
                        if (sw.ElapsedMilliseconds - lastUi > 150)
                        {
                            lastUi = sw.ElapsedMilliseconds;
                            Avalonia.Threading.Dispatcher.UIThread.Post(() =>
                            {
                                GrowLog = $"Growing (seed {s0})… {d} of {t} chains finished · {r} restarts · {sw.Elapsed.TotalSeconds:F1} s";
                                GrowDone = d;
                                GrowRestarts = r;
                                GrowElapsed = sw.Elapsed.TotalSeconds;
                            });
                        }
                        return !token.IsCancellationRequested;
                    };
                    var lbl = label.Replace($"seed{_growSeed}", $"seed{s0}");
                    result = await Task.Run(() => spec == null ? CapsDocument.Grow(oa, onProgress, lbl) : CapsDocument.GrowChains(spec, oa, onProgress, lbl));
                }
                catch (InvalidOperationException e) when (e.Message != "cancelled")
                {
                    failures.Add($"seed {s0}: {e.Message}");
                }
            }
            sw.Stop();
            if (result == null) throw new InvalidOperationException(string.Join("\n", failures) + "\nLower the density, the chain length or the contact scale.");
            var (doc, report) = result.Value;
            var used = seed - 1;
            var name = label.Replace($"seed{_growSeed}", $"seed{used}");
            Show(doc, name + " (unsaved)");
            GrownUnsaved = true;
            GrowDone = _growChains;
            GrowElapsed = sw.Elapsed.TotalSeconds;
            GrowLog = (failures.Count > 0 ? string.Join("\n", failures) + $"\nused seed {used} instead\n" : "") + report + $"\nbuilt in {sw.Elapsed.TotalSeconds:F2} s";
            Status = $"Grown {name} · save it as LAMMPS data, PDB or XYZ";
        }
        catch (Exception e)
        {
            GrowLog = e.Message == "cancelled" ? "Cancelled." : "Could not grow.\n" + e.Message;
            Status = e.Message == "cancelled" ? "Cancelled" : "Could not grow the cell — see the Grow panel";
        }
        finally
        {
            Growing = false;
        }
    }

    public void CancelGrow() => _growCancel?.Cancel();

    // ---------------------------------------------------------------- Relax
    public static readonly string[] Minimisers = ["Steepest descent", "Polak–Ribière conjugate gradient", "L-BFGS (m = 10)", "FIRE"];
    private int _relaxMethod = 2, _relaxIterations = 5000;
    private double _relaxFtol = 0.5, _relaxDensity = 1.05, _relaxStep = 0.06, _relaxPressure = 1.0, _relaxCutoff = 10.0;
    private bool _relaxPushoff = true, _relaxCompress = true, _relaxBox, _relaxCoulomb = true, _relaxing;
    private string _relaxLog = "Relaxes the structure in the viewer with GAFF (C and H in this version): capped-force push-off, " +
                               "compression to a target density, then minimisation to the force tolerance.";
    private string _fieldInfo = "";
    private CancellationTokenSource? _relaxCancel;
    private readonly List<(double X, double Y)> _relaxEnergy = new(), _relaxForce = new();

    public int RelaxMethod { get => _relaxMethod; set => Set(ref _relaxMethod, value); }
    public bool RelaxPushoff { get => _relaxPushoff; set => Set(ref _relaxPushoff, value); }
    public bool RelaxCompress { get => _relaxCompress; set => Set(ref _relaxCompress, value); }
    public bool RelaxBox { get => _relaxBox; set => Set(ref _relaxBox, value); }
    public bool RelaxCoulomb { get => _relaxCoulomb; set => Set(ref _relaxCoulomb, value); }
    public bool Relaxing { get => _relaxing; private set { if (Set(ref _relaxing, value)) RaiseBusy(); } }
    public bool CanRelax => _doc != null && Idle;
    public string RelaxLog { get => _relaxLog; private set => Set(ref _relaxLog, value); }
    public string FieldInfoText { get => _fieldInfo; private set { if (Set(ref _fieldInfo, value)) Raise(nameof(HasFieldInfo)); } }
    public bool HasFieldInfo => _fieldInfo.Length > 0;
    public (double X, double Y)[] RelaxEnergyCurve => _relaxEnergy.ToArray();
    public (double X, double Y)[] RelaxForceCurve => _relaxForce.ToArray();
    public event Action? RelaxCurvesChanged;

    public decimal? RelaxFtolD { get => (decimal)_relaxFtol; set { _relaxFtol = Math.Clamp((double)(value ?? 0.5m), 0.001, 100); Raise(); } }
    public decimal? RelaxIterationsD { get => _relaxIterations; set { _relaxIterations = Math.Clamp((int)(value ?? 5000), 10, 1000000); Raise(); } }
    public decimal? RelaxDensityD { get => (decimal)_relaxDensity; set { _relaxDensity = Math.Clamp((double)(value ?? 1.05m), 0.05, 3); Raise(); } }
    public decimal? RelaxStepD { get => (decimal)_relaxStep; set { _relaxStep = Math.Clamp((double)(value ?? 0.06m), 0.005, 0.5); Raise(); } }
    public decimal? RelaxPressureD { get => (decimal)_relaxPressure; set { _relaxPressure = (double)(value ?? 1m); Raise(); } }
    public decimal? RelaxCutoffD { get => (decimal)_relaxCutoff; set { _relaxCutoff = Math.Clamp((double)(value ?? 10m), 4, 30); Raise(); } }

    public void CheckField()
    {
        if (_doc == null || _relaxing) return;
        try { FieldInfoText = _doc.FieldInfo(); }
        catch (Exception e) { FieldInfoText = "Cannot type this structure: " + e.Message; }
    }

    public async Task Relax()
    {
        if (_doc == null || !Idle) return;
        var doc = _doc;
        Relaxing = true;
        IsPlaying = false;
        _relaxCancel = new CancellationTokenSource();
        var token = _relaxCancel.Token;
        var o = new CapsRelaxOpts
        {
            Method = _relaxMethod, Ftol = _relaxFtol, MaxIterations = _relaxIterations,
            TargetDensity = _relaxCompress ? _relaxDensity : 0, CompressStep = _relaxStep,
            Pushoff = _relaxPushoff ? 1 : 0, RelaxBox = _relaxBox ? 1 : 0, Pressure = _relaxPressure,
            Cutoff = _relaxCutoff, Coulomb = _relaxCoulomb ? 1 : 0,
        };
        _relaxEnergy.Clear();
        _relaxForce.Clear();
        RelaxCurvesChanged?.Invoke();
        RelaxLog = "Typing and minimising…";
        Status = $"Relaxing {Title} with {Minimisers[_relaxMethod]}…";
        var sw = System.Diagnostics.Stopwatch.StartNew();
        var inv = CultureInfo.InvariantCulture;
        try
        {
            var lastUi = 0L;
            var total = 0;
            int lastStage = -1, lastIt = 0;
            // Every progress point is kept (worker side); the UI copies the series at most ten times a second.
            var series = new List<(double X, double E, double F)>();
            var finished = false;   // progress lines queued before the end must not overwrite the final report
            void Publish(string? line)
            {
                (double, double, double)[] copy;
                lock (series) copy = series.ToArray();
                Avalonia.Threading.Dispatcher.UIThread.Post(() =>
                {
                    _relaxEnergy.Clear();
                    _relaxForce.Clear();
                    foreach (var (x, e, f) in copy)
                    {
                        _relaxEnergy.Add((x, e));
                        _relaxForce.Add((x, f));
                    }
                    if (line != null && !finished) RelaxLog = line;
                    RelaxCurvesChanged?.Invoke();
                });
            }
            var (converged, report) = await Task.Run(() =>
            {
                var r = doc.Relax(o, (st, n, it, e, f, d) =>
                {
                    // cumulative iterations across stages for the plots
                    if (st != lastStage) { lastStage = st; lastIt = 0; }
                    total += Math.Max(0, it - lastIt);
                    lastIt = it;
                    lock (series) series.Add((total, e, Math.Log10(Math.Max(f, 1e-6))));
                    if (sw.ElapsedMilliseconds - lastUi > 100)
                    {
                        lastUi = sw.ElapsedMilliseconds;
                        Publish(string.Format(inv, "stage {0} of {1} · iteration {2} · E {3:F1} kcal/mol · |F|max {4:F3} · {5:F3} g/cm³ · {6:F1} s",
                            st, n, it, e, f, d, sw.Elapsed.TotalSeconds));
                    }
                    return !token.IsCancellationRequested;
                });
                Publish(null);
                return r;
            });
            sw.Stop();
            finished = true;
            RelaxLog = report + string.Format(inv, "\nfinished in {0:F1} s{1}", sw.Elapsed.TotalSeconds, converged ? "" : " · force tolerance not reached");
            AfterRun(doc, " · relaxed");
            Status = converged ? $"Relaxed · {Frames} frames (start, each stage, final) · save it as LAMMPS data with the force field"
                               : "Relax stopped before the force tolerance — see the Relax panel";
        }
        catch (Exception e)
        {
            var cancelled = e.Message.Contains("cancelled");
            RelaxLog = cancelled ? "Cancelled; the structure is unchanged." : "Could not relax.\n" + e.Message;
            Status = cancelled ? "Relax cancelled" : "Could not relax — see the Relax panel";
        }
        finally
        {
            Relaxing = false;
            RelaxCurvesChanged?.Invoke();
        }
    }

    public void CancelRelax() => _relaxCancel?.Cancel();

    // ---------------------------------------------------------------- Dynamics
    public static readonly string[] Ensembles = ["NVE (no thermostat)", "NVT", "NPT (isotropic)"];
    public static readonly string[] Thermostats = ["Bussi velocity rescaling", "Langevin (BAOAB)"];
    public static readonly string[] Barostats = ["Stochastic cell rescaling", "Berendsen (early relaxation only)"];
    private int _mdEnsemble = 1, _mdThermostat, _mdBarostat, _mdSeed = 1;
    private double _mdDt = 1.0, _mdTemp = 300, _mdTauT = 100, _mdPressure = 1, _mdTauP = 1000;
    private long _mdSteps = 20000;
    private int _mdFrameEvery = 1000;
    private bool _mdNewVelocities, _mdRunning;
    private string _mdLog = "Runs molecular dynamics on the structure in the viewer with the Relax force field (GAFF, C and H). " +
                            "Relax first: dynamics from a grown cell with overlaps is unstable.";
    private CancellationTokenSource? _mdCancel;
    private readonly List<CapsThermo> _thermo = new();

    public int MdEnsemble { get => _mdEnsemble; set { if (Set(ref _mdEnsemble, value)) { Raise(nameof(MdHasThermostat)); Raise(nameof(MdHasBarostat)); } } }
    public bool MdHasThermostat => _mdEnsemble >= 1;
    public bool MdHasBarostat => _mdEnsemble == 2;
    public int MdThermostat { get => _mdThermostat; set => Set(ref _mdThermostat, value); }
    public int MdBarostat { get => _mdBarostat; set => Set(ref _mdBarostat, value); }
    public bool MdNewVelocities { get => _mdNewVelocities; set => Set(ref _mdNewVelocities, value); }
    public bool MdRunning { get => _mdRunning; private set { if (Set(ref _mdRunning, value)) RaiseBusy(); } }
    public bool CanRun => _doc != null && Idle;
    public string MdLog { get => _mdLog; private set => Set(ref _mdLog, value); }
    public IReadOnlyList<CapsThermo> Thermo => _thermo;
    public event Action? ThermoChanged;

    public decimal? MdDtD { get => (decimal)_mdDt; set { _mdDt = Math.Clamp((double)(value ?? 1m), 0.1, 5); Raise(); Raise(nameof(MdEstimate)); } }
    public decimal? MdStepsD { get => _mdSteps; set { _mdSteps = Math.Clamp((long)(value ?? 20000), 0, 1_000_000_000); Raise(); Raise(nameof(MdEstimate)); } }
    public decimal? MdTempD { get => (decimal)_mdTemp; set { _mdTemp = Math.Clamp((double)(value ?? 300m), 1, 5000); Raise(); } }
    public decimal? MdTauTD { get => (decimal)_mdTauT; set { _mdTauT = Math.Clamp((double)(value ?? 100m), 1, 1e6); Raise(); } }
    public decimal? MdPressureD { get => (decimal)_mdPressure; set { _mdPressure = (double)(value ?? 1m); Raise(); } }
    public decimal? MdTauPD { get => (decimal)_mdTauP; set { _mdTauP = Math.Clamp((double)(value ?? 1000m), 10, 1e7); Raise(); } }
    public decimal? MdFrameEveryD { get => _mdFrameEvery; set { _mdFrameEvery = Math.Clamp((int)(value ?? 1000), 1, 1_000_000); Raise(); Raise(nameof(MdEstimate)); } }
    public decimal? MdSeedD { get => _mdSeed; set { _mdSeed = Math.Max(0, (int)(value ?? 1)); Raise(); } }
    // ---- pre-flight and export (Dynamics board)
    public ObservableCollection<CheckRow> MdPreflight { get; } = new();
    private string _mdPreflightSummary = "", _mdDeck = "";
    public string MdPreflightSummary { get => _mdPreflightSummary; private set => Set(ref _mdPreflightSummary, value); }
    /// <summary>The LAMMPS input for this run (setup from the core, ensemble lines from the settings).</summary>
    public string MdDeck { get => _mdDeck; private set => Set(ref _mdDeck, value); }

    public void RefreshPreflight()
    {
        MdPreflight.Clear();
        if (_doc == null) { MdPreflightSummary = ""; MdDeck = ""; return; }
        var inv = CultureInfo.InvariantCulture;
        var s = _doc.Summary();
        // force field
        if (Field.Assigned)
            MdPreflight.Add(new CheckRow(Field.Complete ? "All atoms typed, no missing parameters" : "The Field assignment is incomplete: runs are blocked",
                Field.Complete ? "ok" : "fail"));
        else
        {
            string ff;
            try { ff = _doc.FieldInfo(); } catch (Exception e) { ff = "error: " + e.Message; }
            var ok = !ff.StartsWith("error", StringComparison.Ordinal) && !ff.Contains("Cannot", StringComparison.Ordinal);
            var name = ok ? ff.Split('\n')[0] : "";
            MdPreflight.Add(new CheckRow(ok ? $"All atoms typed with {name} (none assigned in Field)" : "The built-in force fields cannot type this structure: assign one in Field", ok ? "ok" : "fail"));
        }
        // charge
        var q = s.HasCharges != 0 ? s.TotalCharge : 0;
        MdPreflight.Add(new CheckRow(string.Format(inv, "Net charge {0:+0.000;−0.000;0.000} e", q), Math.Abs(q) < 1e-3 ? "ok" : "check"));
        // box against the cut-off
        if (s.CellValid != 0)
        {
            var w = Math.Min(s.CellA, Math.Min(s.CellB, s.CellC));
            MdPreflight.Add(new CheckRow(string.Format(inv, "Box {0:F1} Å ≥ 2 r_c ({1:F0} Å) in every direction", w, 2 * _relaxCutoff), w >= 2 * _relaxCutoff ? "ok" : "check"));
        }
        else MdPreflight.Add(new CheckRow("No periodic cell: the run is in vacuum", "check"));
        if (_mdEnsemble == 2 && s.CellValid == 0) MdPreflight.Add(new CheckRow("NPT needs a periodic cell", "fail"));
        // time step
        MdPreflight.Add(new CheckRow(string.Format(inv, "Δt {0:0.##} fs with hydrogens, no bond constraints", _mdDt), _mdDt <= 1.0 ? "ok" : _mdDt <= 2.0 ? "check" : "fail"));
        // velocities
        MdPreflight.Add(new CheckRow(_mdNewVelocities ? string.Format(inv, "New velocities at {0:0} K (seed {1})", _mdTemp, _mdSeed) : "Velocities from the structure, or drawn at the target if it has none", "ok"));
        var fails = MdPreflight.Count(r => r.State == "fail");
        var checks = MdPreflight.Count(r => r.State == "check");
        MdPreflightSummary = $"{MdPreflight.Count - fails - checks} / {MdPreflight.Count} ok";
        // LAMMPS deck
        try
        {
            var setup = _doc.LammpsInput("system.data");
            var steps = _mdSteps;
            var ens = _mdEnsemble switch
            {
                0 => "fix 1 all nve",
                1 => _mdThermostat == 1 ? string.Format(inv, "fix 1 all nve\nfix 2 all langevin {0:0.##} {0:0.##} {1:0.##} {2}", _mdTemp, _mdTauT, _mdSeed + 1)
                                         : string.Format(inv, "fix 1 all nve\nfix 2 all temp/csvr {0:0.##} {0:0.##} {1:0.##} {2}", _mdTemp, _mdTauT, _mdSeed + 1),
                _ => string.Format(inv, "fix 1 all nve\nfix 2 all temp/csvr {0:0.##} {0:0.##} {1:0.##} {2}\nfix 3 all press/berendsen iso {3:0.##} {3:0.##} {4:0.##} modulus 22222",
                    _mdTemp, _mdTauT, _mdSeed + 1, _mdPressure, _mdTauP),   // modulus 1/β for β = 4.5e-5 atm⁻¹, as CAPS's barostat
            };
            MdDeck = "# LAMMPS input written by CAPS Studio: the same force field and settings as this Dynamics run\n" + setup +
                     (_mdNewVelocities ? string.Format(inv, "velocity all create {0:0.##} {1} mom yes rot yes dist gaussian\n", _mdTemp, _mdSeed) : "") +
                     string.Format(inv, "timestep {0:0.###}\n{1}\nthermo {2}\ndump d all custom {3} traj.lammpstrj id mol type xu yu zu\nrun {4}\n",
                         _mdDt, ens, Math.Max(1, _mdFrameEvery / 10), _mdFrameEvery, steps);
        }
        catch (Exception e) { MdDeck = "# cannot write the LAMMPS input: " + e.Message; }
    }

    public string MdEstimate => string.Format(CultureInfo.InvariantCulture, "{0:0.###} ps · {1:N0} frames recorded",
        _mdSteps * _mdDt / 1000, _mdSteps / Math.Max(1, _mdFrameEvery) + 1);

    public async Task RunMd()
    {
        if (_doc == null || !Idle) return;
        var doc = _doc;
        MdRunning = true;
        IsPlaying = false;
        _mdCancel = new CancellationTokenSource();
        var token = _mdCancel.Token;
        var o = new CapsMdOpts
        {
            Dt = _mdDt, Steps = _mdSteps, Temperature = _mdTemp,
            Thermostat = _mdEnsemble == 0 ? 0 : _mdThermostat + 1, TauT = _mdTauT,
            Barostat = _mdEnsemble == 2 ? _mdBarostat + 1 : 0, Pressure = _mdPressure, TauP = _mdTauP,
            NewVelocities = _mdNewVelocities ? 1 : 0, Seed = (ulong)_mdSeed,
            ThermoEvery = (int)Math.Clamp(_mdSteps / 400, 10, 1000), FrameEvery = _mdFrameEvery,
            Cutoff = _relaxCutoff, Coulomb = _relaxCoulomb ? 1 : 0, Tail = 1,
        };
        _thermo.Clear();
        ThermoChanged?.Invoke();
        MdLog = "Starting…";
        Status = $"Running {Ensembles[_mdEnsemble]} dynamics on {Title}…";
        var sw = System.Diagnostics.Stopwatch.StartNew();
        var inv = CultureInfo.InvariantCulture;
        var rows = new List<CapsThermo>();
        var lastUi = 0L;
        var finished = false;   // progress lines queued before the end must not overwrite the final report
        void Publish(string? line)
        {
            CapsThermo[] copy;
            lock (rows) copy = rows.ToArray();
            Avalonia.Threading.Dispatcher.UIThread.Post(() =>
            {
                _thermo.Clear();
                _thermo.AddRange(copy);
                if (line != null && !finished) MdLog = line;
                ThermoChanged?.Invoke();
            });
        }
        try
        {
            var report = await Task.Run(() =>
            {
                var r = doc.Md(o, (row, n) =>
                {
                    lock (rows) rows.Add(row);
                    if (sw.ElapsedMilliseconds - lastUi > 150)
                    {
                        lastUi = sw.ElapsedMilliseconds;
                        var frac = n > 0 ? (double)row.Step / n : 1;
                        var eta = frac > 0.001 ? sw.Elapsed.TotalSeconds * (1 - frac) / frac : 0;
                        Publish(string.Format(inv, "step {0:N0} of {1:N0} · {2:F2} ps · T {3:F1} K · P {4:F0} atm · ρ {5:F4} g/cm³ · E {6:F1} kcal/mol · {7:F0} s left",
                            row.Step, n, row.TimePs, row.Temperature, row.Pressure, row.Density, row.Total, eta));
                    }
                    return !token.IsCancellationRequested;
                });
                Publish(null);
                return r;
            });
            sw.Stop();
            finished = true;
            MdLog = report;
            AfterRun(doc, " · MD");
            Status = $"Dynamics finished · {Frames} frames · save the trajectory or the final structure (with velocities)";
        }
        catch (Exception e)
        {
            finished = true;
            var cancelled = e.Message.Contains("cancelled");
            MdLog = cancelled ? "Cancelled; the structure is unchanged." : "Could not run dynamics.\n" + e.Message;
            Status = cancelled ? "Dynamics cancelled" : "Could not run dynamics — see the Dynamics panel";
        }
        finally
        {
            MdRunning = false;
            ThermoChanged?.Invoke();
        }
    }

    public void CancelMd() => _mdCancel?.Cancel();

    // ---------------------------------------------------------------- Equilibrate
    public static readonly string[] Protocols = ["Larsen et al. 21-step (2011)", "Simulated annealing", "MD push-off", "Custom (edit the text)"];
    private static readonly string[] ProtocolNames = ["larsen21", "annealing", "pushoff", ""];
    private int _eqProtocol, _eqCycles = 3, _eqMaxBlocks = 20;
    private double _eqTFinal = 300, _eqTMax = 600, _eqPFinal = 1, _eqPMax = 49346.2, _eqScale = 1, _eqTLow = 300, _eqTHigh = 600, _eqRamp = 50, _eqHold = 50, _eqBlock = 20;
    private bool _eqUntil, _eqRunning;
    private string _eqText = "", _eqLog = "Runs a published equilibration protocol as a chain of Dynamics stages (thermostat, barostat, cut-off and electrostatics as in the Dynamics and Relax panels).";
    private CancellationTokenSource? _eqCancel;

    public int EqProtocol { get => _eqProtocol; set { if (Set(ref _eqProtocol, value)) { Raise(nameof(EqIsAnnealing)); Raise(nameof(EqIsNamed)); RegenerateProtocol(); } } }
    public bool EqIsAnnealing => _eqProtocol == 1;
    public bool EqIsNamed => _eqProtocol != 3;
    public string EqText { get => _eqText; set { if (Set(ref _eqText, value)) Raise(nameof(EqTotal)); } }
    public string EqLog { get => _eqLog; private set => Set(ref _eqLog, value); }
    public bool EqUntilConverged { get => _eqUntil; set => Set(ref _eqUntil, value); }
    public bool EqRunning { get => _eqRunning; private set { if (Set(ref _eqRunning, value)) RaiseBusy(); } }
    public bool CanEquilibrate => _doc != null && Idle;

    private decimal? D(ref double f, decimal? v, double lo, double hi, [System.Runtime.CompilerServices.CallerMemberName] string? name = null)
    {
        f = Math.Clamp((double)(v ?? (decimal)f), lo, hi);
        Raise(name);
        RegenerateProtocol();
        return v;
    }
    public decimal? EqTFinalD { get => (decimal)_eqTFinal; set => D(ref _eqTFinal, value, 1, 5000); }
    public decimal? EqTMaxD { get => (decimal)_eqTMax; set => D(ref _eqTMax, value, 1, 5000); }
    public decimal? EqPFinalD { get => (decimal)_eqPFinal; set => D(ref _eqPFinal, value, -1e5, 1e6); }
    public decimal? EqPMaxD { get => (decimal)_eqPMax; set => D(ref _eqPMax, value, 1, 1e6); }
    public decimal? EqScaleD { get => (decimal)_eqScale; set => D(ref _eqScale, value, 0.001, 100); }
    public decimal? EqTLowD { get => (decimal)_eqTLow; set => D(ref _eqTLow, value, 1, 5000); }
    public decimal? EqTHighD { get => (decimal)_eqTHigh; set => D(ref _eqTHigh, value, 1, 5000); }
    public decimal? EqRampD { get => (decimal)_eqRamp; set => D(ref _eqRamp, value, 0.1, 1e6); }
    public decimal? EqHoldD { get => (decimal)_eqHold; set => D(ref _eqHold, value, 0.1, 1e6); }
    public decimal? EqCyclesD { get => _eqCycles; set { _eqCycles = Math.Clamp((int)(value ?? 3), 1, 100); Raise(); RegenerateProtocol(); } }
    public decimal? EqBlockD { get => (decimal)_eqBlock; set { _eqBlock = Math.Clamp((double)(value ?? 20m), 0.1, 1e6); Raise(); } }
    public decimal? EqMaxBlocksD { get => _eqMaxBlocks; set { _eqMaxBlocks = Math.Clamp((int)(value ?? 20), 1, 1000); Raise(); } }

    /// <summary>Total simulated time of the protocol text (stage durations in ps, ns or fs).</summary>
    public string EqTotal
    {
        get
        {
            double ps = 0;
            var stages = 0;
            foreach (var line in _eqText.Split('\n'))
            {
                var body = line.Split('#')[0];
                var m = System.Text.RegularExpressions.Regex.Match(body, @"([0-9.]+)\s*(ps|ns|fs)\b", System.Text.RegularExpressions.RegexOptions.IgnoreCase);
                if (!m.Success || !double.TryParse(m.Groups[1].Value, NumberStyles.Float, CultureInfo.InvariantCulture, out var v)) continue;
                stages++;
                ps += m.Groups[2].Value.ToLowerInvariant() switch { "ns" => v * 1000, "fs" => v / 1000, _ => v };
            }
            var scaled = _eqScale != 1 && EqIsNamed ? string.Format(CultureInfo.InvariantCulture, " · durations × {0:0.###} (shortened, not the published schedule)", _eqScale) : "";
            if (_eqScale == 1 || !EqIsNamed) scaled = "";
            return string.Format(CultureInfo.InvariantCulture, "{0} stages · {1:N1} ps{2}", stages, ps, scaled);
        }
    }

    private void RegenerateProtocol()
    {
        if (!EqIsNamed) return;
        var p = new CapsProtocolParams
        {
            TFinal = _eqTFinal, TMax = _eqTMax, PFinal = _eqPFinal, PMax = _eqPMax, TimeScale = _eqScale,
            Cycles = _eqCycles, TLow = _eqTLow, THigh = _eqTHigh, RampPs = _eqRamp, HoldPs = _eqHold,
        };
        try { EqText = CapsDocument.ProtocolText(ProtocolNames[_eqProtocol], p); }
        catch (Exception e) { EqLog = e.Message; }
    }

    public void InitProtocol() => RegenerateProtocol();

    public async Task RunEquilibrate()
    {
        if (_doc == null || !Idle) return;
        var doc = _doc;
        EqRunning = true;
        IsPlaying = false;
        _eqCancel = new CancellationTokenSource();
        var token = _eqCancel.Token;
        var o = new CapsEquilOpts
        {
            Dt = _mdDt, Thermostat = _mdThermostat + 1, Barostat = _mdBarostat + 1, TauT = _mdTauT, TauP = _mdTauP, Seed = (ulong)_mdSeed,
            Cutoff = _relaxCutoff, Coulomb = _relaxCoulomb ? 1 : 0, Tail = 1,
            FramePs = 10, ThermoPs = 0.5, UntilConverged = _eqUntil ? 1 : 0, BlockPs = _eqBlock, MaxBlocks = _eqMaxBlocks,
        };
        _thermo.Clear();
        ThermoChanged?.Invoke();
        EqLog = "Starting…";
        Status = $"Equilibrating {Title}: {Protocols[_eqProtocol]}…";
        var inv = CultureInfo.InvariantCulture;
        var sw = System.Diagnostics.Stopwatch.StartNew();
        var rows = new List<CapsThermo>();
        var lastUi = 0L;
        var finished = false;
        var text = _eqText;
        void Publish(string? line)
        {
            CapsThermo[] copy;
            lock (rows) copy = rows.ToArray();
            Avalonia.Threading.Dispatcher.UIThread.Post(() =>
            {
                _thermo.Clear();
                _thermo.AddRange(copy);
                if (line != null && !finished) EqLog = line;
                ThermoChanged?.Invoke();
            });
        }
        try
        {
            var (converged, report) = await Task.Run(() =>
            {
                var r = doc.Equilibrate(text, o, (st, n, label, row) =>
                {
                    lock (rows) rows.Add(row);
                    if (sw.ElapsedMilliseconds - lastUi > 150)
                    {
                        lastUi = sw.ElapsedMilliseconds;
                        Publish(string.Format(inv, "stage {0} of {1}: {2}\n{3:F2} ps · T {4:F1} K · P {5:F0} atm · ρ {6:F4} g/cm³ · {7:F0} s so far",
                            st, n, label, row.TimePs, row.Temperature, row.Pressure, row.Density, sw.Elapsed.TotalSeconds));
                    }
                    return !token.IsCancellationRequested;
                });
                Publish(null);
                return r;
            });
            finished = true;
            EqLog = report;
            AfterRun(doc, " · equilibrated");
            LoadEqChecks(doc.EquilibrateChecks());
            Status = _eqUntil && !converged ? "Protocol finished; the convergence checks did not pass — see the Equilibrate panel"
                                            : $"Equilibrated · {Frames} frames · save the trajectory or the final structure";
        }
        catch (Exception e)
        {
            finished = true;
            var cancelled = e.Message.Contains("cancelled");
            EqLog = cancelled ? "Cancelled; the structure is unchanged." : "Could not equilibrate.\n" + e.Message;
            Status = cancelled ? "Equilibration cancelled" : "Could not equilibrate — see the Equilibrate panel";
        }
        finally
        {
            EqRunning = false;
            ThermoChanged?.Invoke();
        }
    }

    public void CancelEquilibrate() => _eqCancel?.Cancel();

    // ---- convergence (Convergence board): criteria and block means
    public ObservableCollection<CriterionRow> EqCriteria { get; } = new();
    private string _eqCriteriaText = "", _eqPhase = "set up";
    public string EqCriteriaText { get => _eqCriteriaText; private set => Set(ref _eqCriteriaText, value); }
    /// <summary>Mean Rg per production block (Å).</summary>
    public (double X, double Y)[] EqRgBlocks { get; private set; } = [];
    public event Action? EqChecksChanged;

    private void LoadEqChecks(string json)
    {
        EqCriteria.Clear();
        EqRgBlocks = [];
        var inv = CultureInfo.InvariantCulture;
        if (json.Length > 0)
        {
            using var js = System.Text.Json.JsonDocument.Parse(json);
            var root = js.RootElement;
            foreach (var c in root.GetProperty("checks").EnumerateArray())
            {
                var q = c.GetProperty("quantity").GetString() ?? "";
                var ok = c.GetProperty("ok").GetBoolean();
                var ch = c.GetProperty("change").GetDouble();
                var tol = c.GetProperty("tolerance").GetDouble();
                var energy = q.StartsWith("potential", StringComparison.Ordinal);
                var rule = energy ? string.Format(inv, "change < {0:0.###} kcal/mol per atom between blocks, twice", tol)
                                  : string.Format(inv, "change < {0:0.#} % between blocks, twice", 100 * tol);
                var now = energy ? string.Format(inv, "last {0:0.####}", ch) : string.Format(inv, "last {0:0.##} %", 100 * ch);
                EqCriteria.Add(new CriterionRow(char.ToUpperInvariant(q[0]) + q[1..], rule, now, ok ? "pass" : "not yet"));
                if (q.Contains("Rg", StringComparison.Ordinal))
                    EqRgBlocks = c.GetProperty("blocks").EnumerateArray().Select((v, k) => ((double)(k + 1), v.GetDouble())).ToArray();
            }
            var met = EqCriteria.Count(r => r.State == "pass");
            EqCriteriaText = $"{met} of {EqCriteria.Count} criteria met";
        }
        else EqCriteriaText = "";
        EqChecksChanged?.Invoke();
    }

    // ---------------------------------------------------------------- Pack
    private double _packX = 40, _packY = 40, _packZ = 40, _packTol = 2.0;
    private int _packCount = 100, _packSeed = 1;
    private bool _packPeriodic = true, _packing;
    private string _packText = "", _packBaseDir = Environment.CurrentDirectory;
    private string _packLog = "Write or open a packmol-style input. Every molecule is a rigid body; overlaps below the tolerance are " +
                              "removed by minimisation, and a cell that misses the tolerance is never produced.";
    private CancellationTokenSource? _packCancel;

    public decimal? PackXD { get => (decimal)_packX; set { _packX = Math.Clamp((double)(value ?? 40m), 5, 10000); Raise(); } }
    public decimal? PackYD { get => (decimal)_packY; set { _packY = Math.Clamp((double)(value ?? 40m), 5, 10000); Raise(); } }
    public decimal? PackZD { get => (decimal)_packZ; set { _packZ = Math.Clamp((double)(value ?? 40m), 5, 10000); Raise(); } }
    public decimal? PackTolD { get => (decimal)_packTol; set { _packTol = Math.Clamp((double)(value ?? 2m), 0.5, 10); Raise(); } }
    public decimal? PackCountD { get => _packCount; set { _packCount = Math.Clamp((int)(value ?? 100), 1, 10_000_000); Raise(); } }
    public decimal? PackSeedD { get => _packSeed; set { _packSeed = Math.Max(0, (int)(value ?? 1)); Raise(); } }
    public bool PackPeriodic { get => _packPeriodic; set => Set(ref _packPeriodic, value); }
    public bool Packing { get => _packing; private set { if (Set(ref _packing, value)) RaiseBusy(); } }
    public string PackText { get => _packText; set { if (Set(ref _packText, value)) ParsePackText(); } }

    // the input read back as the board's "Molecules & regions" list, and the objective's settings
    public ObservableCollection<PackItem> PackItems { get; } = new();
    private string _packTolText = "2.0 Å", _packCellText = "no periodic cell";
    public string PackTolText { get => _packTolText; private set => Set(ref _packTolText, value); }
    public string PackCellText { get => _packCellText; private set => Set(ref _packCellText, value); }
    public bool HasPackItems => PackItems.Count > 0;
    private static readonly string[] PackColours = ["#F0A83C", "#6CC4D8", "#DE775D", "#9B7AD5", "#7DC884", "#D6AC5C", "#E9ECEF", "#2271DB"];

    private void ParsePackText()
    {
        PackItems.Clear();
        var inv = CultureInfo.InvariantCulture;
        string? file = null;
        var number = "1";
        var fixedMol = false;
        var constraints = new List<string>();
        foreach (var raw in _packText.Split('\n'))
        {
            var line = raw.Split('#')[0].Trim();
            if (line.Length == 0) continue;
            var w = line.Split((char[]?)null, StringSplitOptions.RemoveEmptyEntries);
            var key = w[0].ToLowerInvariant();
            if (file == null)
            {
                if (key == "tolerance" && w.Length > 1) PackTolText = w[1] + " Å";
                else if (key == "pbc" && w.Length >= 7)
                    PackCellText = string.Format(inv, "{0} × {1} × {2} Å, periodic", w[4], w[5], w[6]);
                else if (key == "structure" && w.Length > 1) { file = line[(line.IndexOf(' ') + 1)..].Trim(); number = "1"; fixedMol = false; constraints.Clear(); }
            }
            else if (key == "end") 
            {
                var name = Path.GetFileNameWithoutExtension(file);
                var colour = PackColours[PackItems.Count % PackColours.Length];
                PackItems.Add(new PackItem(name, constraints.Count > 0 ? string.Join(" · ", constraints) : "anywhere in the cell",
                    fixedMol ? "fixed" : $"× {number}", colour, file));
                file = null;
            }
            else if (key == "number" && w.Length > 1) number = w[1];
            else if (key == "fixed") fixedMol = true;
            else if (key is "inside" or "outside" or "over" or "below" or "above") constraints.Add(line);
        }
        Raise(nameof(HasPackItems));
    }

    // live convergence and the guarantee (Pack board)
    public List<(double X, double Y)> PackCurve { get; } = new();
    public event Action? PackCurveChanged;
    private int _packLoop, _packBad;
    private string _packDmin = "—";
    public int PackLoop { get => _packLoop; private set => Set(ref _packLoop, value); }
    public int PackBad { get => _packBad; private set => Set(ref _packBad, value); }
    public string PackDmin { get => _packDmin; private set => Set(ref _packDmin, value); }
    public string PackBaseDir { get => _packBaseDir; set { if (Set(ref _packBaseDir, value)) Raise(nameof(PackBaseNote)); } }
    public string PackBaseNote => "relative structure paths are read from " + _packBaseDir;
    public string PackLog { get => _packLog; private set => Set(ref _packLog, value); }

    private string BoxText() => string.Format(CultureInfo.InvariantCulture, "0. 0. 0. {0:0.###} {1:0.###} {2:0.###}", _packX, _packY, _packZ);

    /// <summary>A new input with the cell, tolerance and seed from the fields.</summary>
    public void NewPackInput()
    {
        PackText = string.Format(CultureInfo.InvariantCulture, "# CAPS Pack input (packmol syntax)\ntolerance {0:0.###}\nseed {1}\n", _packTol, _packSeed) +
                   (_packPeriodic ? $"pbc {BoxText()}\n" : "") + "\n";
    }

    /// <summary>Append a structure block for a molecule file, placed inside the whole cell.</summary>
    public void AddPackStructure(string path)
    {
        if (_packText.Trim().Length == 0) NewPackInput();
        PackText = _packText.TrimEnd() + $"\n\nstructure {path}\n  number {_packCount}\n  inside box {BoxText()}\nend structure\n";
    }

    public void AddPackExample(string samples)
    {
        NewPackInput();
        AddPackStructure(Path.Combine(samples, "water.pdb"));
    }

    public void LoadPackInput(string path)
    {
        PackText = File.ReadAllText(path);
        PackBaseDir = Path.GetDirectoryName(path) ?? ".";
    }

    public async Task RunPack()
    {
        if (!Idle || _packText.Trim().Length == 0) return;
        Packing = true;
        _packCancel = new CancellationTokenSource();
        var token = _packCancel.Token;
        var text = _packText;
        var baseDir = _packBaseDir;
        PackLog = "Packing…";
        Status = "Packing…";
        PackCurve.Clear();
        PackCurveChanged?.Invoke();
        PackLoop = 0;
        PackBad = 0;
        PackDmin = "—";
        var sw = System.Diagnostics.Stopwatch.StartNew();
        var inv = CultureInfo.InvariantCulture;
        var lastUi = 0L;
        try
        {
            var (doc, report) = await Task.Run(() => CapsDocument.Pack(text, baseDir, (loop, loops, f, bad) =>
            {
                if (sw.ElapsedMilliseconds - lastUi > 150)
                {
                    lastUi = sw.ElapsedMilliseconds;
                    var line = loop == 0 ? string.Format(inv, "placing molecules… {0:F1} s", sw.Elapsed.TotalSeconds)
                        : string.Format(inv, "round {0} of at most {1} · penalty {2:G3} · {3} molecules in violation · {4:F1} s", loop, loops, f, bad, sw.Elapsed.TotalSeconds);
                    Avalonia.Threading.Dispatcher.UIThread.Post(() =>
                    {
                        if (!_packing) return;
                        PackLog = line;
                        if (loop > 0)
                        {
                            PackLoop = loop;
                            PackBad = bad;
                            PackCurve.Add((loop, Math.Log10(Math.Max(f, 1e-12))));
                            PackCurveChanged?.Invoke();
                        }
                    });
                }
                return !token.IsCancellationRequested;
            }, "packed"));
            var s = doc.Summary();
            Show(doc, $"packed_{s.Molecules}_molecules (unsaved)");
            GrownUnsaved = true;
            Packing = false;
            PackLog = report;
            PackBad = 0;
            var m = System.Text.RegularExpressions.Regex.Match(report, @"smallest distance between molecules ([0-9.]+) Å");
            if (m.Success) PackDmin = m.Groups[1].Value + " Å";
            Status = $"Packed {s.Molecules:N0} molecules ({s.Atoms:N0} atoms) · save it, or relax and run dynamics";
        }
        catch (Exception e)
        {
            Packing = false;
            var cancelled = e.Message.Contains("cancelled");
            PackLog = cancelled ? "Cancelled." : "Could not pack.\n" + e.Message;
            Status = cancelled ? "Packing cancelled" : "Could not pack — see the Pack panel";
        }
        finally
        {
            Packing = false;
        }
    }

    public void CancelPack() => _packCancel?.Cancel();

    // ---------------------------------------------------------------- React
    public static readonly string[] ReactionSets = ["C–C crosslink (saturated carbons, H₂ leaves)", "Epoxy–amine (primary + secondary)",
        "Sulfur cure of diene rubber (H–S–S–H donors → C–S–S–C)", "Peroxide cure of diene rubber (allylic C–C)", "Custom (edit the text)"];
    private int _rxSet, _rxCycles = 50, _rxPerCycle = 5, _rxSeed = 1, _rxRelaxIt = 500;
    private double _rxTarget = 1.0, _rxCapture, _rxMdPs = 2, _rxTemp = 500, _rxFa = 2, _rxFb = 4, _rxRatio = 1;
    private bool _rxRelax = true, _reacting;
    private string _rxText = "", _rxLog = "Forms bonds between groups closer than the capture distance, retypes, relaxes and (optionally) runs short " +
                                        "dynamics, cycle after cycle (Polymatic cycle; REACTER-style capture and probability).";
    private CancellationTokenSource? _rxCancel;
    private readonly List<CapsReactCycle> _rxRows = new();
    public event Action? ReactChanged;
    public void RaiseGel() => Raise(nameof(RxGelText));
    public IReadOnlyList<CapsReactCycle> ReactRows => _rxRows;

    public int RxSet { get => _rxSet; set { if (Set(ref _rxSet, value)) { LoadReactionSet(); Raise(nameof(RxShowInsert)); } } }
    // curatives inserted into the cell before a cure (sulfur donors)
    public bool RxShowInsert => _rxSet == 2;
    private string _rxInsertSmiles = "SS";
    private decimal _rxInsertCount = 20;
    public string RxInsertSmiles { get => _rxInsertSmiles; set => Set(ref _rxInsertSmiles, value); }
    public decimal RxInsertCount { get => _rxInsertCount; set => Set(ref _rxInsertCount, Math.Clamp(Math.Round(value), 1, 100000)); }
    public async Task InsertCurative()
    {
        if (_doc == null || !Idle) return;
        var doc = _doc;
        var smiles = _rxInsertSmiles;
        var n = (int)_rxInsertCount;
        try
        {
            Status = $"Inserting {n} × {smiles}…";
            var rep = await Task.Run(() => doc.InsertMolecules(smiles, n, 2.0, (ulong)_rxSeed));
            Field.Reset();
            AfterRun(doc, $" · +{n} {smiles}");
            RxLog = rep;
            Status = $"{n} × {smiles} inserted · react with the sulfur cure";
        }
        catch (Exception e) { RxLog = "Could not insert: " + e.Message; Status = "Could not insert the curative"; }
    }
    public string RxText { get => _rxText; set => Set(ref _rxText, value); }
    public string RxLog { get => _rxLog; private set => Set(ref _rxLog, value); }
    public bool RxRelax { get => _rxRelax; set => Set(ref _rxRelax, value); }
    public bool Reacting { get => _reacting; private set { if (Set(ref _reacting, value)) RaiseBusy(); } }
    public bool CanReact => _doc != null && Idle;
    public decimal? RxCyclesD { get => _rxCycles; set { _rxCycles = Math.Clamp((int)(value ?? 50), 1, 100000); Raise(); } }
    public decimal? RxPerCycleD { get => _rxPerCycle; set { _rxPerCycle = Math.Clamp((int)(value ?? 5), 1, 100000); Raise(); } }
    public decimal? RxSeedD { get => _rxSeed; set { _rxSeed = Math.Max(0, (int)(value ?? 1)); Raise(); } }
    public decimal? RxTargetD { get => (decimal)_rxTarget; set { _rxTarget = Math.Clamp((double)(value ?? 1m), 0.001, 1); Raise(); } }
    public decimal? RxCaptureD { get => (decimal)_rxCapture; set { _rxCapture = Math.Clamp((double)(value ?? 0m), 0, 20); Raise(); } }
    public decimal? RxMdPsD { get => (decimal)_rxMdPs; set { _rxMdPs = Math.Clamp((double)(value ?? 0m), 0, 10000); Raise(); } }
    public decimal? RxTempD { get => (decimal)_rxTemp; set { _rxTemp = Math.Clamp((double)(value ?? 300m), 1, 5000); Raise(); } }
    public decimal? RxFaD { get => (decimal)_rxFa; set { _rxFa = Math.Max(1, (double)(value ?? 2m)); Raise(); Raise(nameof(FloryText)); Raise(nameof(RxAlphaC)); } }
    public decimal? RxFbD { get => (decimal)_rxFb; set { _rxFb = Math.Max(1, (double)(value ?? 4m)); Raise(); Raise(nameof(FloryText)); Raise(nameof(RxAlphaC)); } }
    public decimal? RxRatioD { get => (decimal)_rxRatio; set { _rxRatio = Math.Clamp((double)(value ?? 1m), 0.01, 1); Raise(); Raise(nameof(FloryText)); Raise(nameof(RxAlphaC)); } }
    /// <summary>α_c = 1/√(r (fA−1)(fB−1)), or "—".</summary>
    public string RxAlphaC
    {
        get
        {
            var d = _rxRatio * (_rxFa - 1) * (_rxFb - 1);
            return d > 0 ? (1 / Math.Sqrt(d)).ToString("F3", CultureInfo.InvariantCulture) : "—";
        }
    }
    /// <summary>Simulated gel point: the conversion where the reduced weight-average mass peaks (cluster analysis).</summary>
    public string RxGelText
    {
        get
        {
            if (_rxRows.Count < 3) return "—";
            var best = _rxRows.MaxBy(r => r.ReducedMw);
            var last = _rxRows[^1];
            if (best.Cycle == last.Cycle) return string.Format(CultureInfo.InvariantCulture, "not reached (α {0:F3})", last.Conversion);
            return string.Format(CultureInfo.InvariantCulture, "α {0:F3}", best.Conversion);
        }
    }
    public string FloryText
    {
        get
        {
            var d = _rxRatio * (_rxFa - 1) * (_rxFb - 1);
            return d > 0 ? string.Format(CultureInfo.InvariantCulture, "Flory–Stockmayer α_c = 1/√(r (fA−1)(fB−1)) = {0:F3}", 1 / Math.Sqrt(d))
                         : "Flory–Stockmayer: no gel point (a functionality of 1)";
        }
    }

    public void LoadReactionSet()
    {
        try
        {
            RxText = _rxSet switch
            {
                0 => CapsDocument.ReactionTemplate("cc_crosslink"),
                1 => CapsDocument.ReactionTemplate("epoxy_amine_primary") + "\n" + CapsDocument.ReactionTemplate("epoxy_amine_secondary"),
                2 => CapsDocument.ReactionTemplate("sulfur_allylic"),
                3 => CapsDocument.ReactionTemplate("peroxide_allylic"),
                _ => _rxText,
            };
        }
        catch (Exception e) { RxLog = e.Message; }
    }

    public async Task RunReact()
    {
        if (_doc == null || !Idle || _rxText.Trim().Length == 0) return;
        var doc = _doc;
        Reacting = true;
        IsPlaying = false;
        _rxCancel = new CancellationTokenSource();
        var token = _rxCancel.Token;
        var o = new CapsReactOpts
        {
            Seed = (ulong)_rxSeed, MaxCycles = _rxCycles, MaxPerCycle = _rxPerCycle, TargetConversion = _rxTarget, Capture = _rxCapture,
            Relax = _rxRelax ? 1 : 0, RelaxIterations = _rxRelaxIt, MdPs = _rxRelax ? _rxMdPs : 0, Temperature = _rxTemp, Cutoff = _relaxCutoff, Coulomb = _relaxCoulomb ? 1 : 0,
        };
        _rxRows.Clear();
        ReactChanged?.Invoke();
        RxLog = "Finding reactive pairs…";
        Status = $"Reacting {Title}…";
        var text = _rxText;
        var inv = CultureInfo.InvariantCulture;
        var sw = System.Diagnostics.Stopwatch.StartNew();
        var finished = false;
        try
        {
            var report = await Task.Run(() => doc.React(text, o, row =>
            {
                Avalonia.Threading.Dispatcher.UIThread.Post(() =>
                {
                    _rxRows.Add(row);
                    if (!finished)
                        RxLog = string.Format(inv, "cycle {0}: {1} reactions ({2} in all) · conversion {3:F3} · {4} clusters, largest {5:F1} % · {6:F0} s",
                            row.Cycle, row.Reactions, row.Total, row.Conversion, row.Clusters, 100 * row.LargestFraction, sw.Elapsed.TotalSeconds);
                    ReactChanged?.Invoke();
                });
                return !token.IsCancellationRequested;
            }));
            finished = true;
            RxLog = report + "\n" + FloryText;
            AfterRun(doc, " · reacted");
            Status = "Reaction run finished · save the network (LAMMPS data carries the force field when every atom is typed)";
        }
        catch (Exception e)
        {
            finished = true;
            var cancelled = e.Message.Contains("cancelled");
            RxLog = cancelled ? "Cancelled; the structure is unchanged." : "Could not react.\n" + e.Message;
            Status = cancelled ? "Reaction run cancelled" : "Could not react — see the React panel";
        }
        finally
        {
            Reacting = false;
            ReactChanged?.Invoke();
        }
    }

    public void CancelReact() => _rxCancel?.Cancel();

    public void SaveTrajectory(string path)
    {
        if (_doc == null) return;
        _doc.SaveTrajectory(path);
        Status = $"Saved {Frames} frames to {path}";
    }

    // The document now holds the run's record (relaxation stages or MD frames); show its last frame, keep the camera.
    private void AfterRun(CapsDocument doc, string suffix)
    {
        Field.LoadReport(doc);   // runs keep the assignment (React replaces the topology and ends it)
        var s = doc.Summary();
        Frames = (int)Math.Max(1, s.Frames);
        Raise(nameof(FrameMax));
        _frame = Frames - 1;
        Raise(nameof(Frame));
        Raise(nameof(FrameLabel));
        var t = Title.Replace(" (unsaved)", "");
        Title = (t.EndsWith(suffix) ? t : t + suffix) + " (unsaved)";
        GrownUnsaved = true;
        _selection.Clear();
        RefreshSelection();
        RefreshSummary();
        RefreshRdf();
        RefreshLegend();
        RefreshMolecules();
        Notes.Clear();
        foreach (var n in doc.Notes()) Notes.Add(n);
        LoadFileChecks();
        FieldInfoText = "";
        RenderRequested?.Invoke();
    }

    public void SaveDocument(string path)
    {
        if (_doc == null) return;
        _doc.Save(path);
        GrownUnsaved = false;
        Title = System.IO.Path.GetFileName(path);
        Status = $"Saved {path}";
        Remember(path, null);
    }

    private void Show(CapsDocument doc, string title)
    {
        if (Busy) { doc.Dispose(); Status = "Wait for the run to finish (or cancel it) before opening another structure"; return; }
        if (_wrap) doc.SetWrap(true);
        Document?.Dispose();
        Document = doc;
        ClearFocus();
        if (IsVisualize) Avalonia.Threading.Dispatcher.UIThread.Post(ApplyPipeline);
        Field.Reset();
        Analyze.Load("");
        SyncHeld();
        Title = title;
        FieldInfoText = "";
        GrownUnsaved = false;
        var s = doc.Summary();
        Frames = (int)Math.Max(1, s.Frames);
        Raise(nameof(FrameMax));
        _frame = 0;
        Raise(nameof(Frame));
        Raise(nameof(FrameLabel));
        Camera = new CapsCamera { Yaw = 0.55, Pitch = 0.40, Zoom = 1.0, Perspective = _perspective ? 1 : 0 };
        _selection.Clear();
        RefreshSelection();
        RefreshSummary();
        RefreshRdf();
        RefreshLegend();
        RefreshMolecules();
        IsPlaying = false;
        Notes.Clear();
        foreach (var n in doc.Notes()) Notes.Add(n);
        LoadFileChecks();
        var look = FileChecks.Count(c => c.NeedsLook);
        Status = $"Opened {Title} · {s.Atoms.ToString("N0", CultureInfo.InvariantCulture)} atoms · {s.Format}" + (look > 0 ? $" · {look} file check{(look == 1 ? "" : "s")} need a look" : "");
        RenderRequested?.Invoke();
    }

    private void RefreshSummary()
    {
        if (_doc == null) return;
        var s = _doc.Summary();
        Analyze.OnDocument(Title, s);
        var inv = CultureInfo.InvariantCulture;
        SummaryRows.Clear();
        StatusCounts = string.Format(inv, "{0:N0} atoms · {1:N0} bonds{2}", s.Atoms, s.Bonds, s.Frames > 1 ? $" · {s.Frames:N0} frames" : "");
        StatusCell = "Cell: none (non-periodic)";
        var mols = string.Format(inv, s.Molecules == 1 ? "{0:N0} molecule" : "{0:N0} molecules", s.Molecules);
        HudInfo = s.CellValid != 0 ? string.Format(inv, "{0} · {1:F3} g/cm³", mols, s.Density) : mols;
        SummaryRows.Add(new("Format", s.Format));
        SummaryRows.Add(new("Atoms", s.Atoms.ToString("N0", inv)));
        SummaryRows.Add(new("Bonds", $"{s.Bonds.ToString("N0", inv)} · {(s.BondsFromFile != 0 ? "from file" : "perceived")}"));
        SummaryRows.Add(new("Molecules", s.Molecules.ToString("N0", inv)));
        SummaryRows.Add(new("Frames", s.Frames.ToString(inv)));
        if (s.CellValid != 0)
        {
            SummaryRows.Add(new("Cell", string.Format(inv, "{0:F2} × {1:F2} × {2:F2} Å", s.CellA, s.CellB, s.CellC)));
            StatusCell = string.Format(inv, "Cell: {0:F2} × {1:F2} × {2:F2} Å", s.CellA, s.CellB, s.CellC);
            SummaryRows.Add(new("Volume", string.Format(inv, "{0:N0} Å³", s.Volume)));
            SummaryRows.Add(new("Density", string.Format(inv, "{0:F4} g/cm³", s.Density)));
        }
        SummaryRows.Add(new("Mass", string.Format(inv, "{0:N1} g/mol", s.TotalMass)));
        if (s.HasCharges != 0) SummaryRows.Add(new("Total charge", string.Format(inv, "{0:+0.0e+0;-0.0e+0;0} e", s.TotalCharge)));
    }

    /// <summary>Click selects one atom; with add (Shift/Cmd) atoms are appended, up to four, for measurements.</summary>
    public void Pick(int index, bool add = false)
    {
        if (!add) _selection.Clear();
        if (index >= 0)
        {
            if (add && _selection.Contains(index)) _selection.Remove(index);
            else
            {
                _selection.Add(index);
                if (_selection.Count > 4) _selection.RemoveAt(0);
            }
        }
        RefreshSelection();
    }

    // status bar
    private string _statusCounts = "", _statusCell = "";
    public string StatusCounts { get => _statusCounts; private set => Set(ref _statusCounts, value); }
    public string StatusCell { get => _statusCell; private set => Set(ref _statusCell, value); }
    public string StatusSelection => $"Selection {_selection.Count}";
    private string _hudInfo = "";
    public string HudInfo { get => _hudInfo; private set => Set(ref _hudInfo, value); }

    private void RefreshSelection()
    {
        Raise(nameof(StatusSelection));
        Raise(nameof(HudSelection));
        Raise(nameof(HasSelection));
        Raise(nameof(HasAtom));
        Raise(nameof(Picked));
        Raise(nameof(HasPicked));
        PickedRows.Clear();
        NeighbourRows.Clear();
        var index = Picked;
        if (_doc == null || index < 0)
        {
            PickedTitle = "Nothing picked";
            MeasureText = "";
            return;
        }
        var a = _doc.Atom(index);
        var inv = CultureInfo.InvariantCulture;
        PickedTitle = $"Atom {a.Id} · {a.ElementSymbol}" + (string.IsNullOrEmpty(a.Name) ? "" : $" · {a.Name}");
        // Atom inspector (Main board): element, coordinates, type, charge and why the type
        PickedElement = a.ElementSymbol;
        PickedX = a.X.ToString("F3", inv);
        PickedY = a.Y.ToString("F3", inv);
        PickedZ = a.Z.ToString("F3", inv);
        PickedType = string.IsNullOrEmpty(a.Name) ? a.Type.ToString(inv) : a.Name;
        PickedCharge = string.Format(inv, "q {0:+0.0000;−0.0000;0} e", a.Charge);
        PickedMolecule = a.Mol.ToString(inv);
        if (Field.Assigned)
        {
            PickedWhy = Field.Explain(index);
            PickedFf = Field.ForceFieldName;
        }
        else
        {
            PickedWhy = "Types come from the file. Assign a force field in Field to see the rule behind each type.";
            PickedFf = "file types";
        }
        Raise(nameof(HasAtom));
        PickedRows.Add(new("Molecule", a.Mol.ToString(inv)));
        PickedRows.Add(new("Type", a.Type.ToString(inv)));
        PickedRows.Add(new("Charge", string.Format(inv, "{0:+0.0000;-0.0000;0} e", a.Charge)));
        PickedRows.Add(new("Position", string.Format(inv, "{0:F2}  {1:F2}  {2:F2} Å", a.X, a.Y, a.Z)));
        foreach (var (i, d) in _doc.Neighbours(index, 4))
        {
            var b = _doc.Atom(i);
            NeighbourRows.Add(new($"{b.Id} · {b.ElementSymbol} · mol {b.Mol}", string.Format(inv, "{0:F3} Å", d)));
        }
        if (_selection.Count >= 2)
        {
            var v = _doc.Measure(_selection.ToArray());
            var ids = string.Join("–", _selection.Select(i => _doc.Atom(i).Id));
            MeasureText = _selection.Count switch
            {
                2 => string.Format(inv, "Distance {0}: {1:F3} Å", ids, v),
                3 => string.Format(inv, "Angle {0}: {1:F2}°", ids, v),
                _ => string.Format(inv, "Dihedral {0}: {1:F2}°", ids, v),
            };
            Status = MeasureText;
        }
        else
        {
            MeasureText = "";
            Status = $"Picked atom {a.Id} ({a.ElementSymbol}) in molecule {a.Mol} · Shift-click more atoms to measure";
        }
    }

    private void RefreshRdf()
    {
        if (_doc == null) { RdfCurve = []; return; }
        var s = _doc.Summary();
        if (s.CellValid == 0) { RdfCurve = []; RdfNote = "g(r) needs a periodic cell"; return; }
        var half = Math.Min(s.CellA, Math.Min(s.CellB, s.CellC)) / 2;
        var rmax = Math.Min(12.0, Math.Floor(half));
        var (ea, eb) = RdfElements[_rdfPair];
        RdfCurve = _doc.Rdf(ea, eb, rmax, 0.2, _rdfInter);
        RdfNote = string.Format(CultureInfo.InvariantCulture, "{0} · {1} · r ≤ {2:F0} Å · 0.2 Å bins · frame {3}",
            RdfPairs[_rdfPair], _rdfInter ? "between molecules" : "all pairs", rmax, _frame);
    }

    public CapsRenderOpts ViewOptions(int w, int h, int supersample) => new()
    {
        Width = w, Height = h, Supersample = supersample,
        Background = _viewBackground,
        ColourBy = _colour, Style = _style,
        Outlines = _outlines ? 1 : 0, ShowCell = _showCell ? 1 : 0,
        Highlight0 = _selection.Count > 0 ? _selection[0] : -1,
        Highlight1 = _selection.Count > 1 ? _selection[1] : -1,
        Highlight2 = _selection.Count > 2 ? _selection[2] : -1,
        Highlight3 = _selection.Count > 3 ? _selection[3] : -1,
        Focus = _focusAtom >= 0 ? _focusAtom + 1 : 0,
        AmbientOcclusion = _module == 19 ? (_renderAo ? 1 : 0) : (_viewAo ? 1 : 0),
        DepthCue = _module == 19 ? (_renderDepth ? 1 : 0) : (_depthCue ? 1 : 0),
    };

    /// <summary>The Field page's view: coloured by force-field type, ball and stick, the selected row's atom highlighted.</summary>
    public CapsRenderOpts FieldViewOptions(int w, int h, int supersample) => new()
    {
        Width = w, Height = h, Supersample = supersample,
        Background = _viewBackground,
        ColourBy = 2, Style = _style == 3 ? 0 : _style,
        Outlines = _outlines ? 1 : 0, DepthCue = _depthCue ? 1 : 0, ShowCell = 0,
        Highlight0 = Field.SelectedRow?.Index ?? -1, Highlight1 = -1, Highlight2 = -1, Highlight3 = -1,
    };

    public (int W, int H) ExportSize(int viewW, int viewH) => _sizePreset switch
    {
        0 => (1920, 1080),
        1 => (2008, 1130),
        2 => (4016, 2259),
        _ => (viewW, viewH),
    };

    public CapsRenderOpts ExportOptions(int w, int h)
    {
        var o = ViewOptions(w, h, 2);
        o.Background = _exportBackground;
        o.Highlight0 = o.Highlight1 = o.Highlight2 = o.Highlight3 = -1;
        o.Focus = 0;
        return o;
    }

    public void ResetView()
    {
        Camera.Yaw = 0.55; Camera.Pitch = 0.40; Camera.Zoom = 1; Camera.PanX = 0; Camera.PanY = 0;
        RenderRequested?.Invoke();
    }

    public void SetView(double yaw, double pitch)
    {
        Camera.Yaw = yaw; Camera.Pitch = pitch; Camera.PanX = 0; Camera.PanY = 0;
        RenderRequested?.Invoke();
    }
}
