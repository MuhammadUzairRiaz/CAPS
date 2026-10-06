using System.Collections.ObjectModel;
using System.ComponentModel;
using System.Globalization;
using System.Text.Json.Nodes;
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
public sealed record PackItem(string Name, string Detail, string Count, string Colour, string File, int Row = 0, int FfChoice = 0)
{
    /// <summary>The copies as a number (a fixed molecule: none to edit).</summary>
    public decimal? CountValue => Count.StartsWith("× ", StringComparison.Ordinal) && int.TryParse(Count[2..], out var n) ? n : null;
    public bool IsFixed => Count == "fixed";
    public Avalonia.Media.IBrush Brush => Avalonia.Media.Brush.Parse(Colour);
}

public sealed partial class MainViewModel : ObservableObject
{
    public static readonly string[] Styles = ["Ball & stick", "Space filling", "Sticks", "No hydrogens", "Backbone"];
    public static readonly string[] ColourModes = ["Element", "Molecule", "Type", "Distance to molecule centre"];
    public static readonly string[] Backgrounds = ["Dark", "White", "Transparent"];
    public static readonly string[] ViewBackgrounds = ["Theme", "Dark", "White (paper)"];
    // g(r) pairs from the elements present: the commonest heavy elements with themselves and each other, with hydrogen,
    // H – H, all – all (a water box offers O – O, O – H …, a polystyrene melt C – C, C – H …)
    public System.Collections.ObjectModel.ObservableCollection<string> RdfPairs { get; } = new(["C – C", "C – H", "H – H", "all – all"]);
    private (int A, int B)[] _rdfElements = [(6, 6), (6, 1), (1, 1), (0, 0)];
    private (CapsDocument? Doc, long Atoms) _rdfPairsFor;

    private void EnsureRdfPairs()
    {
        if (_doc == null) return;
        var atoms = _doc.Summary().Atoms;
        if (_rdfPairsFor.Doc == _doc && _rdfPairsFor.Atoms == atoms) return;
        _rdfPairsFor = (_doc, atoms);
        var counts = new Dictionary<int, (int N, string Sym)>();
        var stride = (int)Math.Max(1, atoms / 20000);   // a sample is enough to rank the elements
        for (var i = 0; i < atoms; i += stride)
        {
            var a = _doc.Atom(i);
            if (a.Element <= 0) continue;
            counts[a.Element] = (counts.GetValueOrDefault(a.Element).N + 1, a.ElementSymbol);
        }
        var heavy = counts.Where(kv => kv.Key != 1).OrderByDescending(kv => kv.Value.N).Take(3).ToList();
        var hasH = counts.ContainsKey(1);
        var pairs = new List<(int, int, string)>();
        for (var x = 0; x < heavy.Count; x++)
            for (var y = x; y < heavy.Count; y++)
                pairs.Add((heavy[x].Key, heavy[y].Key, $"{heavy[x].Value.Sym} – {heavy[y].Value.Sym}"));
        if (hasH)
        {
            foreach (var h in heavy.Take(2)) pairs.Add((h.Key, 1, $"{h.Value.Sym} – H"));
            pairs.Add((1, 1, "H – H"));
        }
        pairs.Add((0, 0, "all – all"));
        var keep = _rdfPair < RdfPairs.Count ? RdfPairs[_rdfPair] : "";
        _rdfElements = pairs.Select(p => (p.Item1, p.Item2)).ToArray();
        RdfPairs.Clear();
        foreach (var p in pairs) RdfPairs.Add(p.Item3);
        _rdfPair = Math.Max(0, RdfPairs.IndexOf(keep));
        Raise(nameof(RdfPairIndex));
    }
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
        ProjectItems.CollectionChanged += OnProjectItemsChanged;   // an open project is saved as its structures change
        PropertyChanged += (_, e) =>
        {
            if (e.PropertyName == nameof(ShowStudioTabs) || e.PropertyName == nameof(IsStudio)) { Raise(nameof(ShowBrushPanel)); Raise(nameof(ShowRulesPanel)); }   // the inspector's panels take turns
            else if (e.PropertyName == nameof(Document)) { if (_brushOpen) RefillBrush(); if (_rulesOpen) LoadRulePairs(); ApplyLook(); }   // each structure in the user's look
        };
        Field = new FieldViewModel(() => _doc, s => Status = s, () =>
        {
            // types and charges in the document changed: summary, inspector and viewer follow
            RefreshSummary();
            RefreshSelection();
            RefreshLegend();
            FieldInfoText = "";
            RenderRequested?.Invoke();
            RaiseGrowField();
            RefreshSteps();
        });
        Field.GroupSuggestions = FieldGroupSuggestions;
        Field.ApplyWaterModel = id => RunEdit(new { op = "water_model", model = id }) != null;
        Analyze = new AnalyzeViewModel(() => _doc, s => Status = s, running => { _analyzing = running; RaiseBusy(); });
        Analyze.ViscosityComputed += eta => DfEta = (decimal)eta;   // the Yeh–Hummer correction takes the Green–Kubo η
        Field.PropertyChanged += (_, e) =>
        {
            if (e.PropertyName == nameof(FieldViewModel.RunLine)) Raise(nameof(ForceFieldLine));
            // Pack and Grow follow the Field's force field until they are given their own: their boxes show it
            if (e.PropertyName == nameof(FieldViewModel.FfIndex)) { if (_packFf < 0) Raise(nameof(PackFfIndex)); if (_growFf < 0) Raise(nameof(GrowFfIndex)); }
        };
        Field.Recorder = Record;
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
                QueueSelBar();   // the bar, and the hidden / ghosted counts, of the new document
                PackDocumentChanged();   // Packing's last result belongs to the structure it made
                _rxAimDc = true;         // a new structure: the crosslinking panel sets the target again
                _tagOnly = null;
                Probes.Clear(); ProbeCards.Clear(); Raise(nameof(HasProbes)); Raise(nameof(HasProbeCards));   // probes are a structure's atoms
                if (value == null) RefreshLayers();   // no structure: no layers (a closed one's rows gone)
                RaiseStepDots();                      // the tab dots are the active structure's
                QueueRxSites();          // the chains' reactive sites of this structure
                if (_visionPreview != 0) try { value?.SetVision(_visionPreview); } catch { /* an older core */ }
                foreach (var n in new[] { nameof(AppColour), nameof(AppSurface), nameof(AppHasSurface), nameof(AppChip), nameof(ShowAppLegend) }) Raise(n);
                RaiseAppearanceVisibility();
                RefreshAppColumns();
                if (value == null) { _pipeDone.Clear(); _pipeBuild = ""; PipelineSteps.Clear(); }   // no structure: the strip empties
                Raise(nameof(ShowPipelineStrip));
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

    public string Title
    {
        get => _title;
        set { if (Set(ref _title, value) && _activeItem != null && value.Length > 0) _activeItem.Name = value; }
    }
    public string Status { get => _status; set => Set(ref _status, value); }

    public int StyleIndex { get => _style; set { if (Set(ref _style, value)) { Raise(nameof(StyleText)); Raise(nameof(DisplayStatus)); Raise(nameof(DsStyle)); RenderRequested?.Invoke(); } } }
    // toolbar texts (Main board: "Ball & stick", "Colour: element", "Perspective")
    public string StyleText => Styles[Math.Clamp(_style, 0, Styles.Length - 1)];
    public string ColourText => "Colour: " + (_appColour == 4 ? "partial charge" : ColourModes[Math.Clamp(ColourIndex, 0, ColourModes.Length - 1)].ToLowerInvariant());
    public string ProjectionText => _perspective ? (Math.Abs(ViewFov - 35) < 0.5 ? "Perspective" : $"Perspective {ViewFov:0}°") : "Orthographic";
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
        Raise(nameof(CanAddRestraint));
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
    /// <summary>The perspective view's field of view, degrees (35 by default; narrower flattens, wider exaggerates depth).</summary>
    public double ViewFov
    {
        get => Camera.FovDeg > 0 ? Camera.FovDeg : 35;
        set
        {
            var v = Math.Clamp(value, 10, 120);
            if (Math.Abs(v - ViewFov) < 1e-9) return;
            Camera.FovDeg = v;
            Raise();
            Raise(nameof(ProjectionText));
            RenderRequested?.Invoke();
            if (IsViewports) RenderViewports();
        }
    }

    public bool Perspective
    {
        get => _perspective;
        set { if (Set(ref _perspective, value)) { Raise(nameof(ProjectionText)); Camera.Perspective = value ? 1 : 0; RenderRequested?.Invoke(); } }
    }

    public int ExportBackground { get => _exportBackground; set => Set(ref _exportBackground, value); }

    // the view's background: the mode chosen (0 follow the theme, 1 dark, 2 white) and what it gives now (0 dark, 1 white),
    // which every renderer of the view reads
    private int _viewBackground, _viewBgMode;
    public int ViewBackground
    {
        get => _viewBackground;
        private set { if (Set(ref _viewBackground, value)) { Raise(nameof(ViewIsLight)); RenderRequested?.Invoke(); } }
    }
    public int ViewBackgroundMode
    {
        get => _viewBgMode;
        set
        {
            if (!Set(ref _viewBgMode, Math.Clamp(value, 0, 2))) { UpdateViewBackground(); return; }
            Raise(nameof(ViewBgTheme)); Raise(nameof(ViewBgDark)); Raise(nameof(ViewBgWhite));
            UpdateViewBackground();
        }
    }
    public bool ViewBgTheme { get => _viewBgMode == 0; set { if (value) ViewBackgroundMode = 0; } }
    public bool ViewBgDark { get => _viewBgMode == 1; set { if (value) ViewBackgroundMode = 1; } }
    public bool ViewBgWhite { get => _viewBgMode == 2; set { if (value) ViewBackgroundMode = 2; } }
    /// <summary>Recomputed when the mode or the theme (the system's too) changes.</summary>
    public void UpdateViewBackground()
    {
        var themeLight = Avalonia.Application.Current?.ActualThemeVariant == Avalonia.Styling.ThemeVariant.Light;
        ViewBackground = _viewBgMode switch { 1 => 0, 2 => 1, _ => themeLight ? 1 : 0 };
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
            if (_module == 34 && _syncFrame) SplitFrameB();
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
    /// <summary>The rotational isomeric state reference (polyethylene, Flory) for alkane chains, else empty.</summary>
    private (double X, double Y)[] _risCurve = [];
    public (double X, double Y)[] RisCurve { get => _risCurve; private set => Set(ref _risCurve, value); }
    public string ChainNote { get => _chainNote; private set => Set(ref _chainNote, value); }

    private void RefreshChains()
    {
        if (_doc == null) { ChainCurve = []; ChainNote = ""; return; }
        try
        {
            var (n, r, chains, b2) = _doc.InternalDistances();
            // alkanes (every molecule C_nH_2n+2): the polyethylene RIS reference at the Dynamics temperature, dashed
            var ris = Array.Empty<(double, double)>();
            if (chains > 0 && n.Length > 1 && IsAlkane(_doc))
            {
                var nmax = (int)n.Max();
                var c = new double[nmax];
                if (Native.RisCn(_mdTemp, nmax, c) == nmax) ris = n.Select(v => ((double)v, c[Math.Clamp(v, 1, nmax) - 1])).ToArray();
            }
            RisCurve = ris;
            Raise(nameof(ChainReference)); Raise(nameof(EqTargetNote));
            ChainCurve = n.Select((v, k) => ((double)v, r[k])).ToArray();
            ChainNote = chains == 0 ? "no chains of four or more heavy atoms"
                : string.Format(CultureInfo.InvariantCulture, "{0} backbones · ⟨b²⟩ {1:F3} Å² · plateau → C∞ for equilibrated long chains · frame {2}", chains, b2, _frame)
                  + (ris.Length > 0 ? string.Format(CultureInfo.InvariantCulture, " · dashed: RIS polyethylene at {0:F0} K (Flory), C_{1} = {2:F2}", _mdTemp, n.Max(), ris[^1].Item2) : "");
        }
        catch (Exception e) { ChainCurve = []; ChainNote = e.Message; }
    }

    /// <summary>Every molecule an acyclic alkane (C_nH_2n+2, n ≥ 4): the RIS polyethylene reference applies.</summary>
    private static bool IsAlkane(CapsDocument doc)
    {
        var s = doc.Summary();
        if (s.Atoms > 200000) return false;
        var c = new Dictionary<long, (int C, int H)>();
        for (var i = 0; i < s.Atoms; ++i)
        {
            var a = doc.Atom(i);
            if (a.Element != 6 && a.Element != 1) return false;
            c.TryGetValue(a.Mol, out var t);
            c[a.Mol] = a.Element == 6 ? (t.C + 1, t.H) : (t.C, t.H + 1);
        }
        return c.Count > 0 && c.Values.All(t => t.C >= 4 && t.H == 2 * t.C + 2);
    }

    private void RefreshMolecules()
    {
        RefreshChains();
        MoleculeRows.Clear();
        RefreshLayers();
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
        if (IsReaderFile(path)) { OpenReader(path); return; }   // a log, input or table: read, not opened as a structure
        if (path.EndsWith(".capstable", StringComparison.OrdinalIgnoreCase)) { Status = StudyLoad(path); return; }
        var frames = _pendingFrames;   // a frame selection from the Open page, for this open only
        _pendingFrames = null;
        if (OpensProgressively(path)) { _ = OpenProgressive(path, topology, frames); return; }
        Show(frames is { } f ? CapsDocument.OpenStaged(path, topology, 0, null, f) : CapsDocument.Open(path, topology), System.IO.Path.GetFileName(path));
        if (_doc?.Path == path) { Remember(path, topology); RecordOpen(path, topology, frames); LoadKeptForceField(path); }
    }

    /// <summary>Several trajectory files of one run (dumps written in parts, restarts) opened as one trajectory.</summary>
    public void OpenJoined(IReadOnlyList<string> paths, string? topology = null)
    {
        var doc = CapsDocument.OpenJoined(paths, topology);
        var stem = System.IO.Path.GetFileNameWithoutExtension(paths[0]);
        Show(doc, $"{stem} + {paths.Count - 1} more (joined)");
        Status = $"{paths.Count} files joined in time order · {doc.Summary().Frames} frames";
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
    private static readonly string[] Crumbs = ["Polymer cell › Amorphous cell", "Analyze › Properties", "Minimise", "Dynamics › Run",
        "Equilibrate › Protocol", "Packing › Molecules & regions", "React › Crosslinking", "Force field › Typing report", "Studio", "Studio › Molecule", "Settings", "Jobs", "Bench", "Builders › Polymer", "Builders › Surface", "Builders › Nanostructure", "Builders › Polymer › Blend", "Studio › File checks", "Export › Figure", "Studio › Render", "Analyze › Visualize", "Export › Data", "Analyze › Batch", "Analyze › Compare", "Analyze › Visualize › Colour by", "Studio › Viewports", "Export › Figure bundle", "Open file", "Analyze › Visualize › Save pipeline", "Builders › Crystal", "Builders › Biomolecule", "Builders › Solvation", "Studio › Trajectory", "Studio › Torsion scan", "Studio › Split view", "Studio › Fragment library", "Studio › Macro recorder", "Jobs › Provenance", "Analyze › Mechanics", "Analyze › Scattering", "Analyze › Free volume", "Theory manual", "Project", "Jobs › Sweep", "Builders › Coarse-grained", "React › Template editor", "Settings › Colour vision", "Analyze › Glass transition", "Analyze › Interface", "Analyze › Diffusion", "Studio › Charges", "Studio › Periodic box", "Analyze › Orientation", "Jobs › Recipes", "Export › Figure composer", "Analyze › Chains", "Packing › Density calculator", "Analyze › Surface area", "Studio › Unit cell",
        "Polymer cell › Polydispersity", "Builders › Copolymer", "Analyze › Solvent screen", "Builders › Polymer › Tacticity", "Analyze › Blend phase diagram", "Dynamics › Electrostatics",
        "Studio › Display styles", "Studio › Add hydrogens", "Studio › Model resolution", "Export", "Force field › Type by hand", "Analyze › Adsorption locator", "Analyze › Sorption", "Polymer cell › Mesoscale (DPD)", "Studio › Reader", "Analyze › Study table"];
    /// <summary>Where the user is (top bar).</summary>
    public string Crumb => _module == 8 ? "" : Crumbs[_module];
    /// <summary>Where calculations run (top bar).</summary>
    public string ComputeText
    {
        get
        {
            var threads = _settings.Threads > 0 ? _settings.Threads : Math.Min(16, Environment.ProcessorCount);
            var hosts = _settings.Hosts.Count;
            var remote = Jobs.Count(j => j.IsRemote && j.IsRunning);
            return $"Local · {threads} threads · {(_settings.GpuView ? "GPU view" : "CPU")}" +
                   (hosts > 0 ? $" · {hosts} host{(hosts == 1 ? "" : "s")}" : "") + (remote > 0 ? $" · {remote} remote job{(remote == 1 ? "" : "s")}" : "");
        }
    }
    /// <summary>The top bar's compute button: Settings › Compute &amp; remote.</summary>
    public void OpenComputeSettings() { SettingsTab = 3; SetModule(10); }
    /// <summary>Kept for scripts and tests: the property calculations are the Analyze module.</summary>
    public bool AnalyzeProperties { get => _module == 1; set { if (value) SetModule(1); else if (_module == 1) SetModule(8); } }
    public bool IsProperties => _module == 1;
    /// <summary>Each page's own set-up (its Open… method), run when the page is reached any other way: the rail, the
    /// command palette, a Back button — so no page shows empty because it was not opened through its button.</summary>
    private static readonly Dictionary<int, Action<MainViewModel>> PageOpeners = new()
    {
        [9] = v => v.OpenBuilder(),
        [13] = v => { if (v._polyDoc == null) _ = v.BuildPolyPreview(); },   // the one-chain preview, once
        [14] = v => v.OpenSurface(),
        [15] = v => v.OpenNano(),
        [16] = v => v.OpenBlend(),
        [17] = v => v.OpenChecks(),
        [18] = v => v.OpenFigure(),
        [19] = v => v.OpenRender(),
        [20] = v => v.OpenVisualize(),
        [21] = v => v.OpenExport(),
        [22] = v => v.OpenBatch(),
        [23] = v => v.OpenCompare(),
        [24] = v => v.OpenColourBy(),
        [25] = v => v.OpenViewports(),
        [26] = v => v.OpenBundle(),
        [28] = v => v.OpenSavePipeline(),
        [29] = v => v.OpenCrystal(),
        [30] = v => v.OpenBio(),
        [31] = v => v.OpenSolvation(),
        [32] = v => v.OpenTrajectory(),
        [33] = v => v.OpenTorsion(),
        [34] = v => v.OpenSplit(),
        [35] = v => v.OpenFragments(),
        [36] = v => v.OpenMacro(),
        [37] = v => v.OpenProvenance(),
        [38] = v => v.OpenMechanics(),
        [39] = v => v.OpenScattering(),
        [40] = v => v.OpenFreeVolume(),
        [41] = v => v.OpenManual(),
        [42] = v => v.OpenProject(),
        [43] = v => v.OpenSweep(),
        [44] = v => v.OpenCg(),
        [45] = v => v.OpenTemplateEditor(),
        [46] = v => v.OpenColourVision(),
        [47] = v => v.OpenGlass(),
        [48] = v => v.OpenInterface(),
        [49] = v => v.OpenDiffusion(),
        [50] = v => v.OpenCharges(),
        [51] = v => v.OpenPeriodic(),
        [52] = v => v.OpenOrientation(),
        [53] = v => v.OpenRecipes(),
        [54] = v => v.OpenComposer(),
        [55] = v => v.OpenChainStats(),
        [56] = v => v.OpenDensityCalc(),
        [57] = v => v.OpenSurfaceArea(),
        [58] = v => v.OpenCellEditor(),
        [59] = v => v.OpenPolydispersity(),
        [60] = v => v.OpenCopolymer(),
        [61] = v => v.OpenSolventScreen(),
        [62] = v => v.OpenTacticityStats(),
        [63] = v => v.OpenBlendPhase(),
        [64] = v => v.OpenElectrostatics(),
        [65] = v => v.OpenDisplayStyles(),
        [66] = v => _ = v.OpenAddHydrogens(),
        [67] = v => v.OpenModelResolution(),
        [68] = v => v.OpenExportCenter(),
        [69] = v => v.OpenUnitTyping(),
        [70] = v => v.OpenAdsorption(),
        [71] = v => v.OpenSorption(),
        [72] = v => v.OpenDpd(),
        [73] = v => { if (v.ReaderPath.Length > 0) v.SetModule(73); },
        [74] = v => v.OpenStudy(),
    };

    public void SetModule(int m, [System.Runtime.CompilerServices.CallerMemberName] string caller = "")
    {
        var was = _module;
        if (!Set(ref _module, m, nameof(Module))) return;
        if (!caller.StartsWith("Open", StringComparison.Ordinal) && PageOpeners.TryGetValue(m, out var open))
            Avalonia.Threading.Dispatcher.UIThread.Post(() => { if (_module == m) open(this); });
        if (was is 19 or 20 or 21 && m is not (19 or 20)) SuspendPipeline();   // Render (Visualize › Render) keeps the pipeline
        if (m == 20) ApplyPipeline();
        if (m == 19 || was == 19) RenderRequested?.Invoke();   // the view takes (or gives back) the render background
        if (was == 50) EndChargePreview();   // the view shows the structure's own charges again
        if (was == 40 && _doc != null) { try { _doc.Voids("{\"clear\":true}"); _fvVoids = 0; Raise(nameof(FvHasVoids)); RenderRequested?.Invoke(); } catch { } }
        Raise(nameof(IsGrow));
        Raise(nameof(IsAnalyze));
        Raise(nameof(IsRelax));
        Raise(nameof(IsDynamics));
        Raise(nameof(IsEquilibrate));
        Raise(nameof(IsPack));
        Raise(nameof(IsReact));
        Raise(nameof(IsField));
        Raise(nameof(IsStudio)); Raise(nameof(HasSelBar));
        Raise(nameof(IsMolecule));
        Raise(nameof(IsStudioRail));
        Raise(nameof(IsSettings));
        Raise(nameof(IsJobs));
        // Jobs opens on a job (a running one first), never on an empty detail beside a full list
        if (m == 11 && SelectedJob == null && Jobs.Count > 0) SelectedJob = Jobs.FirstOrDefault(j => j.IsRunning) ?? Jobs[0];
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
        Raise(nameof(IsSplit));
        Raise(nameof(IsFragments));
        Raise(nameof(IsMacro));
        Raise(nameof(IsProvenance));
        Raise(nameof(IsMechanics));
        Raise(nameof(IsManual));
        Raise(nameof(IsProject));
        Raise(nameof(IsSweep));
        Raise(nameof(IsCg));
        Raise(nameof(IsTemplate));
        Raise(nameof(IsColourVision));
        Raise(nameof(IsGlass));
        Raise(nameof(IsInterfacePage));
        Raise(nameof(IsDiffusion));
        Raise(nameof(IsCharges));
        Raise(nameof(IsAdsorption));
        Raise(nameof(IsSorption));
        Raise(nameof(IsDpd));
        Raise(nameof(IsReader));
        Raise(nameof(IsStudy));
        Raise(nameof(IsPeriodic));
        Raise(nameof(IsOrientation));
        Raise(nameof(IsRecipes));
        Raise(nameof(IsComposer));
        Raise(nameof(IsChainStats));
        Raise(nameof(IsDensityCalc));
        Raise(nameof(IsSurfaceArea));
        Raise(nameof(IsCellEditor));
        Raise(nameof(IsPolydispersity)); Raise(nameof(IsCopolymer)); Raise(nameof(IsSolventScreen)); Raise(nameof(IsTacticityStats)); Raise(nameof(IsBlendPhase)); Raise(nameof(IsElectrostatics));
        Raise(nameof(IsDisplayStyles)); Raise(nameof(IsAddHydrogens)); Raise(nameof(IsModelResolution)); Raise(nameof(IsUnitTyping)); Raise(nameof(ShowLensPanel));
        if (was == 57 && m != 57) ClearSurfaceColour();
        if (m != 51) LeavePeriodic();
        Raise(nameof(ShowLodPanel));
        Raise(nameof(ShowHistoryPanel));
        Raise(nameof(ProjectPanelShown));
        Raise(nameof(IsScattering));
        Raise(nameof(IsFreeVolume));
        SyncFocusChips();
        RaiseAppearanceVisibility();
        Raise(nameof(IsAnalyzeRail));
        Raise(nameof(ShowPipeLegend));
        Raise(nameof(ShowAnalysisPanel));
        Raise(nameof(ShowLegend));
        Raise(nameof(Crumb));
        Raise(nameof(IsProperties));
        Raise(nameof(IsExportCenter)); Raise(nameof(IsBuildRail)); Raise(nameof(IsExportRail));
        if (m == 68) RefreshEngines();
        RefreshSteps();
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
    public int GrowSeed { get => _growSeed; set { if (Set(ref _growSeed, Math.Max(1, value)) && _growSeedC != null && _growSeedC.Value != _growSeed) _growSeedC.Value = _growSeed; } }
    public double GrowDensity { get => _growDensity; set => Set(ref _growDensity, Math.Clamp(value, 0.01, 2.0)); }
    public double GrowBox { get => _growBox; set => Set(ref _growBox, Math.Max(0, value)); }
    public double GrowScale { get => _growScale; set => Set(ref _growScale, Math.Clamp(value, 0.5, 1.2)); }
    public bool GrowUseBox { get => _growUseBox; set { if (Set(ref _growUseBox, value)) Raise(nameof(GrowUseDensity)); } }
    // region shape: 0 cubic, 1 slab with vacuum, 2 inside a cylinder along z, 3 around a cylinder (a fibre's place)
    public static readonly string[] GrowShapes = ["Cubic, periodic", "Slab with vacuum", "Cylinder (pore)", "Around a cylinder (fibre)"];
    private int _growShape;
    private double _growSlabH = 30, _growSlabVac = 30, _growCylR = 10, _growCylLen;
    public int GrowShape
    {
        get => _growShape;
        set
        {
            if (!Set(ref _growShape, Math.Clamp(value, 0, 3))) return;
            if (_growShape > 0 && _growUseBox) GrowUseBox = false;   // the region's edges follow from the density
            foreach (var n in new[] { nameof(GrowHasRegion), nameof(GrowRegionALabel), nameof(GrowRegionBLabel), nameof(GrowRegionAD), nameof(GrowRegionBD), nameof(GrowRegionNote), nameof(GrowEstimate) }) Raise(n);
        }
    }
    public bool GrowHasRegion => _growShape > 0;
    // growth method: 0 roomiest trial, 1 Rosenbluth soft spheres, 2 Rosenbluth UFF Lennard-Jones
    private int _growMethod;
    private double _growMethodT = 450;
    public int GrowMethod
    {
        get => _growMethod;
        set { if (Set(ref _growMethod, Math.Clamp(value, 0, 2))) { Raise(nameof(GrowMethodIs0)); Raise(nameof(GrowMethodIs1)); Raise(nameof(GrowMethodIs2)); } }
    }
    public bool GrowMethodIs0 => _growMethod == 0;
    public bool GrowMethodIs1 => _growMethod == 1;
    public bool GrowMethodIs2 => _growMethod == 2;
    public decimal GrowMethodTempD { get => (decimal)_growMethodT; set { _growMethodT = Math.Clamp((double)value, 100, 2000); Raise(); } }
    private static readonly string[] GrowMethodIds = ["trials", "rosenbluth", "rosenbluth_lj"];
    // orientation (design/boards/Grow): isotropic, or an aligning field along x, y or z of strength s (kT)
    public static readonly string[] GrowOrientations = ["Isotropic", "Oriented along x", "Oriented along y", "Oriented along z"];
    private int _growOrient;
    private double _growOrientS = 4;
    public int GrowOrient { get => _growOrient; set { if (Set(ref _growOrient, Math.Clamp(value, 0, 3))) Raise(nameof(GrowOriented)); } }
    public bool GrowOriented => _growOrient > 0;
    public decimal GrowOrientStrengthD { get => (decimal)_growOrientS; set { _growOrientS = Math.Clamp((double)value, 0, 50); Raise(); } }
    private System.Text.Json.Nodes.JsonObject? GrowOrientationJson() =>
        _growOrient == 0 ? null : new() { ["axis"] = "xyz"[_growOrient - 1].ToString(), ["strength"] = _growOrientS };
    public string GrowRegionALabel => _growShape == 1 ? "Film thickness (Å)" : "Cylinder radius (Å)";
    public string GrowRegionBLabel => _growShape == 1 ? "Vacuum, above + below (Å)" : "Length along z (Å, 0: from the density)";
    public decimal GrowRegionAD
    {
        get => (decimal)(_growShape == 1 ? _growSlabH : _growCylR);
        set { var v = Math.Clamp((double)value, 2, 500); if (_growShape == 1) _growSlabH = v; else _growCylR = v; Raise(); Raise(nameof(GrowEstimate)); }
    }
    public decimal GrowRegionBD
    {
        get => (decimal)(_growShape == 1 ? _growSlabVac : _growCylLen);
        set { var v = Math.Clamp((double)value, 0, 1000); if (_growShape == 1) _growSlabVac = v; else _growCylLen = v; Raise(); Raise(nameof(GrowEstimate)); }
    }
    public string GrowRegionNote => _growShape switch
    {
        1 => "x and y from the density; the film sits in the middle of the cell",
        2 => "the cell is 2 Å wider than the cylinder; along z periodic",
        3 => "left empty for a fibre; a cube sized from the density",
        _ => "",
    };
    private System.Text.Json.Nodes.JsonObject? GrowRegionJson() => _growShape switch
    {
        1 => new() { ["shape"] = "slab", ["thickness"] = _growSlabH, ["vacuum"] = _growSlabVac },
        2 => new() { ["shape"] = "cylinder", ["radius"] = _growCylR, ["length"] = _growCylLen },
        3 => new() { ["shape"] = "around_cylinder", ["radius"] = _growCylR, ["length"] = _growCylLen },
        _ => null,
    };
    public bool GrowUseDensity => !_growUseBox;
    public bool GrowCurve { get => _growCurve; set => Set(ref _growCurve, value); }
    public bool Growing { get => _growing; private set { if (Set(ref _growing, value)) { Raise(nameof(NotGrowing)); RaiseBusy(); } } }
    public bool NotGrowing => !_growing;
    /// <summary>No build or minimisation is running (the document can be replaced or edited).</summary>
    public bool Idle => !_growing && !_relaxing && !_mdRunning && !_eqRunning && !_packing && !_reacting && !_analyzing && !_adsRunning && !_sorbRunning && !_dpdRunning;
    /// <summary>The core is working on the open document (Relax or Dynamics): no rendering or edits until it is done.</summary>
    public bool Busy => _relaxing || _mdRunning || _eqRunning || _reacting || _analyzing || _recipeRunning || _adsRunning || _sorbRunning || _dpdRunning;
    private bool _analyzing;
    private void RaiseBusy()
    {
        if (!Busy) QueueProjectSave();   // a run ended: the project is saved once the work is still
        Raise(nameof(Idle));
        Raise(nameof(Busy));
        Raise(nameof(CanRelax));
        Raise(nameof(CanRun));
        Raise(nameof(CanEquilibrate));
        Raise(nameof(CanReact));
        Raise(nameof(HasSelBar));
        if (!Busy && _leavePeriodicPending) LeavePeriodic();
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
    /// <summary>The same build as a CAPS recipe (caps run, Start › From a recipe): build, type, grow, export.</summary>
    public string GrowRecipe()
    {
        var inv = CultureInfo.InvariantCulture;
        var name = (_growSpec == null ? "PS" : "polymer") + $"_{_growChains}x{_growDp}";
        var sb = new System.Text.StringBuilder();
        sb.Append("# CAPS recipe from Grow · run: caps run ").Append(name).Append(".yaml\nrecipe: 1\nname: ").Append(name).Append("\nbuild:\n  polymer:\n");
        var j = _growSpec == null ? null : System.Text.Json.Nodes.JsonNode.Parse(_growSpec);
        var units = j?["units"] is System.Text.Json.Nodes.JsonArray ua ? ua.Select(u => (string?)u!["smiles"] ?? "").ToList() : ["*CC(*)c1ccccc1"];
        sb.Append("    units: [").Append(string.Join(", ", units.Select(u => "\"" + u + "\""))).Append("]\n");
        var seq = (string?)j?["sequence"] ?? "homopolymer";
        sb.Append("    sequence: ").Append(seq).Append('\n');
        if (seq == "block" && j?["blocks"] is System.Text.Json.Nodes.JsonArray ba) sb.Append("    blocks: [").Append(string.Join(", ", ba.Select(b => (int?)b ?? 1))).Append("]\n");
        if (seq is "random" or "shuffled" && j?["weights"] is System.Text.Json.Nodes.JsonArray wa) sb.Append("    weights: [").Append(string.Join(", ", wa.Select(w => ((double?)w ?? 1).ToString(inv)))).Append("]\n");
        if (seq == "pattern") sb.Append("    pattern: ").Append((string?)j?["pattern"] ?? "A").Append('\n');
        if (seq == "terminal" && j?["weights"] is System.Text.Json.Nodes.JsonArray tw)
            sb.Append(inv, $"    weights: [{string.Join(", ", tw.Select(w => ((double?)w ?? 0.5).ToString(inv)))}]\n    r1: {((double?)j["r1"] ?? 1).ToString(inv)}\n    r2: {((double?)j["r2"] ?? 1).ToString(inv)}\n");
        sb.Append(inv, $"    dp: {_growDp}\n    chains: {_growChains}\n    tacticity: {Tacticities[_growTact].ToLowerInvariant()}\n");
        if (_growChainDp != null && _growChainDp.Length == _growChains)   // the lengths drawn here, one per chain
            sb.Append("    chain_dp: [").Append(string.Join(", ", _growChainDp)).Append("]   # ").Append(_growChainDpText).Append('\n');
        if (_growTact == 0 && _growStereo is { } gs)
            sb.Append(gs.Kind == 1 ? FormattableString.Invariant($"    p_mr: {gs.A}\n    p_rm: {gs.B}\n") : FormattableString.Invariant($"    pm: {gs.Pm}\n"));
        if ((string?)j?["architecture"] is { } arch and not "linear")   // star, comb, branched
        {
            sb.Append("    architecture: ").Append(arch).Append('\n');
            if (arch == "star") sb.Append(inv, $"    arms: {(int?)j!["arms"] ?? 4}\n");
            else if (arch == "dendrimer") sb.Append(inv, $"    arms: {(int?)j!["arms"] ?? 4}\n    arm_dp: {(int?)j["arm_dp"] ?? 5}\n    generations: {(int?)j["generations"] ?? 2}\n");
            else sb.Append(inv, $"    arm_dp: {(int?)j!["arm_dp"] ?? 5}\n").Append(arch == "comb" ? $"    spacing: {(int?)j["spacing"] ?? 4}\n"
                                                                                               : $"    branch_probability: {((double?)j["branch_probability"] ?? 0.1).ToString(inv)}\n");
        }
        sb.Append("type: { forcefield: default }\n");
        sb.Append("grow:\n").Append(_growUseBox && _growShape == 0 ? $"  box: {_growBox.ToString(inv)}\n" : $"  density: {_growDensity.ToString(inv)}\n");
        if (GrowRegionJson() is { } region) sb.Append("  region: ").Append(region.ToJsonString().Replace("\"", "").Replace(",", ", ").Replace(":", ": ")).Append('\n');
        sb.Append(inv, $"  seed: {_growSeed}\n  contact_scale: {(_growAutoScale ? "auto" : _growScale.ToString(inv))}\n  curve: {(_growCurve ? "true" : "false")}\n");
        if (_growMethod > 0) sb.Append(inv, $"  method: {GrowMethodIds[_growMethod]}\n  temperature: {_growMethodT}\n");
        if (_growLookahead > 1) sb.Append(inv, $"  lookahead: {_growLookahead}\n");
        if (_growOrient > 0) sb.Append(inv, $"  orientation: {{ axis: {"xyz"[_growOrient - 1]}, strength: {_growOrientS} }}\n");
        sb.Append("relax: { method: lbfgs, fmax: 1.0 }\nexport: [lammps, pdb]\n");
        return sb.ToString();
    }

    /// <summary>The same build in Python (the caps package; notebooks).</summary>
    public string GrowPython()
    {
        var inv = CultureInfo.InvariantCulture;
        var j = _growSpec == null ? null : System.Text.Json.Nodes.JsonNode.Parse(_growSpec);
        var units = j?["units"] is System.Text.Json.Nodes.JsonArray ua ? ua.Select(u => (string?)u!["smiles"] ?? "").ToList() : ["*CC(*)c1ccccc1"];
        var smiles = units.Count == 1 ? $"\"{units[0]}\"" : "[" + string.Join(", ", units.Select(u => $"\"{u}\"")) + "]";
        var seq = (string?)j?["sequence"] ?? "homopolymer";
        var extra = seq == "homopolymer" ? "" : $", sequence=\"{seq}\"";
        if ((string?)j?["architecture"] is { } arch and not "linear")
            extra += arch == "star" ? $", architecture=\"star\", arms={(int?)j!["arms"] ?? 4}"
                   : arch == "dendrimer" ? $", architecture=\"dendrimer\", arms={(int?)j!["arms"] ?? 4}, arm_dp={(int?)j["arm_dp"] ?? 5}, generations={(int?)j["generations"] ?? 2}"
                   : string.Format(inv, ", architecture=\"{0}\", arm_dp={1}{2}", arch, (int?)j!["arm_dp"] ?? 5,
                                   arch == "comb" ? $", spacing={(int?)j["spacing"] ?? 4}" : string.Format(inv, ", branch_probability={0}", (double?)j["branch_probability"] ?? 0.1));
        if (_growMethod > 0) extra += string.Format(inv, ", method=\"{0}\", method_temperature={1}", GrowMethodIds[_growMethod], _growMethodT);
        if (_growOrient > 0) extra += string.Format(inv, ", orientation={{\"axis\": \"{0}\", \"strength\": {1}}}", "xyz"[_growOrient - 1], _growOrientS);
        if (GrowRegionJson() is { } region)
            extra += ", region={" + string.Join(", ", region.Select(kv => $"\"{kv.Key}\": " + (kv.Value is System.Text.Json.Nodes.JsonValue v && v.TryGetValue<string>(out var sv) ? $"\"{sv}\"" : kv.Value!.ToJsonString()))) + "}";
        return "import caps\n\n" + string.Format(inv, "cell = caps.polymer({0}, dp={1}, chains={2}, tacticity=\"{3}\", seed={4}, density={5}{6})\n",
                   smiles, _growDp, _growChains, Tacticities[_growTact].ToLowerInvariant(), _growSeed, _growDensity.ToString(inv), extra) +
               "cell.relax(ftol=1.0)\ncell.save(\"cell.data\")\ncell.view()\n";
    }

    private string GrowUnitsArg()
    {
        if (_growSpec == null) return "";
        var j = System.Text.Json.Nodes.JsonNode.Parse(_growSpec)!;
        var units = string.Join(",", (j["units"] as System.Text.Json.Nodes.JsonArray ?? []).Select(u => (string?)u!["smiles"]));
        var seq = (string?)j["sequence"] ?? "homopolymer";
        var extra = seq switch
        {
            "random" or "shuffled" => " --weights " + string.Join(",", (j["weights"] as System.Text.Json.Nodes.JsonArray ?? []).Select(w => ((double?)w ?? 1).ToString(CultureInfo.InvariantCulture))),
            "block" => " --blocks " + string.Join(",", (j["blocks"] as System.Text.Json.Nodes.JsonArray ?? []).Select(w => (int?)w ?? 1)),
            "pattern" => " --pattern " + (string?)j["pattern"],
            _ => "",
        };
        return $" --units '{units}' --sequence {seq}{extra}";
    }
    public bool GrownUnsaved { get => _grownUnsaved; private set => Set(ref _grownUnsaved, value); }

    // decimal views for NumericUpDown
    public decimal? GrowChainsD { get => _growChains; set { GrowChains = (int)(value ?? 1); Raise(); } }
    public decimal? GrowDpD { get => _growDp; set { if (value == null) return; var before = _growDp; GrowDp = (int)value; Raise(); if (_growDp != before) PolyChanged(); } }   // an emptied box keeps the DP; the sequence strip and composition follow at once   // the sequence strip and composition follow at once
    public decimal? GrowSeedD { get => _growSeed; set { if (value == null) return; GrowSeed = (int)value; Raise(); } }
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
            var room = mass / 6.02214076e23 / _growDensity * 1e24;   // Å³ the chains fill at the target density
            var cell = _growShape switch   // as the grower sizes the region (polymer.cpp)
            {
                1 => string.Format(inv, "cell {0:F2} × {0:F2} × {1:F1} Å (film {2:F0} Å)", Math.Sqrt(room / _growSlabH), _growSlabH + _growSlabVac, _growSlabH),
                2 => _growCylLen > 0 ? string.Format(inv, "cylinder r {0:F1} Å × {1:F1} Å", _growCylR, _growCylLen)
                                     : string.Format(inv, "cylinder r {0:F1} Å, {1:F1} Å long", _growCylR, room / (Math.PI * _growCylR * _growCylR)),
                3 => string.Format(inv, "around a cylinder of r {0:F1} Å", _growCylR),
                _ => string.Format(inv, "box {0:F2} Å", box),
            };
            return string.Format(inv, "{0:N0} atoms · {1:N0} g/mol per chain · {2} · {3:F3} g/cm³", atoms, mass / _growChains, cell, rho);
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
        _growSeed = GrowSeedChoice.Take();
        var growScript = Recording ? GrowPython() : null;   // the macro line, with this run's seed
        var o = new CapsGrowOpts
        {
            Chains = _growChains, Dp = _growDp, Tacticity = _growTact, Seed = (ulong)_growSeed,
            Box = _growUseBox ? _growBox : 0, Density = _growUseBox ? 0 : _growDensity, ContactScale = _growScale, Curve = _growCurve ? 1 : 0,
        };
        // every cell goes through the repeat-unit grower (polystyrene as a styrene unit), so the live view and the
        // per-chain lengths (Polydispersity) and stereo model (Tacticity) apply to all of them
        string? spec = GrowSpecWithStatistics(_growSpec) ?? GrowSpecObject().ToJsonString();
        if (spec != null && _growAutoScale) o.ContactScale = -_growScale;   // caps_grow_chains: start here, lower it when crowded
        if (spec != null)
        {
            var sj = System.Text.Json.Nodes.JsonNode.Parse(spec)!.AsObject();
            sj["dp"] = _growDp;
            sj["trials"] = _growTrials;
            if (_growLookahead > 1) sj["lookahead"] = _growLookahead;
            if (GrowRegionJson() is { } region) { sj["region"] = region; o.Box = 0; o.Density = _growDensity; }
            if (_growMethod > 0) { sj["method"] = GrowMethodIds[_growMethod]; sj["temperature"] = _growMethodT; }
            if (GrowOrientationJson() is { } orient) sj["orientation"] = orient;
            spec = sj.ToJsonString();
        }
        var stem = spec == null ? "PS" : string.Concat(_growSpecName.Where(char.IsLetterOrDigit).Take(16));
        var dpTag = _growChainDp != null && _growChainDp.Length == _growChains ? $"Nn{_growChainDp.Average():0}" : _growDp.ToString(CultureInfo.InvariantCulture);
        var label = $"{stem}_{_growChains}x{dpTag}_{Tacticities[_growTact].ToLowerInvariant()}_seed{_growSeed}";
        GrowLog = "Growing…";
        GrowUnitsText = GrowMarginText = GrowDensityNowText = "—";
        GrowUnitFraction = 0;
        GrowDone = 0;
        GrowRestarts = 0;
        GrowElapsed = 0;
        Status = $"Growing {_growChains} chains of {(spec == null ? "polystyrene" : _growSpecName)}, DP {_growDp}…";
        var sw = System.Diagnostics.Stopwatch.StartNew();
        var grown = false;
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
                                if (!_growing) return;   // arrived after the run ended: the finished numbers stay
                                GrowLog = $"Growing (seed {s0})… {d} of {t} chains finished · {r} restarts · {sw.Elapsed.TotalSeconds:F1} s";
                                GrowDone = d;
                                GrowRestarts = r;
                                GrowElapsed = sw.Elapsed.TotalSeconds;
                            });
                        }
                        return !token.IsCancellationRequested;
                    };
                    var lbl = label.Replace($"seed{_growSeed}", $"seed{s0}");
                    var ticket = ++_growLiveTicket;
                    Action<CapsDocument, string> onLive = (live, stats) => Avalonia.Threading.Dispatcher.UIThread.Post(() =>
                    {
                        if (ticket != _growLiveTicket || !_growing) { live.Dispose(); return; }
                        try { live.SetWrap(_growWrap); } catch { }   // folded into the cell (as Amorphous Cell shows them) or whole, as chosen
                        var old = GrowLiveDoc;
                        GrowLiveDoc = live;
                        old?.Dispose();
                        GrowLiveStats(stats);
                    });
                    result = await Task.Run(() => spec == null ? CapsDocument.Grow(oa, onProgress, lbl) : CapsDocument.GrowChains(spec, oa, onProgress, lbl, _growLiveView ? onLive : null));
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
            if (GrowSmall.Count > 0)   // solvents and gases into the free space the chains left
            {
                GrowLog = "Inserting " + string.Join(", ", GrowSmall.Select(r => $"{r.CountD:0} {r.Name}")) + "…";
                try { report += "\n" + await Task.Run(() => InsertGrowSmall(doc, (ulong)used)); }
                catch (Exception e) { report += "\nsmall molecules not inserted: " + e.Message; }
            }
            var name = label.Replace($"seed{_growSeed}", $"seed{used}");
            Show(doc, name + " (unsaved)");
            if (_activeItem != null) _activeItem.BuildSettings = GrowSettingsJson();   // Edit brings these back
            ShowGrownWrap();   // the grown cell as the live view showed it
            MarkPipeline("Grow", spec == null ? "polystyrene" : _growSpecName);
            if (growScript != null) RecordScript(growScript);
            grown = true;
            AfterGrowStatistics(doc);
            var densityNote = SuggestRelaxDensity(spec == null ? "Polystyrene" : _growSpecName);
            GrownUnsaved = true;
            GrowDone = _growChains;
            // the last live snapshot can land after the run ends and is dropped: the finished cell gives the final numbers
            if (double.IsFinite(_growUnitsTotal)) { GrowUnitsText = $"{_growUnitsTotal:0} / {_growUnitsTotal:0}"; GrowUnitFraction = 1; }
            GrowDensityNowText = doc.Summary().Density.ToString("0.00", CultureInfo.InvariantCulture) + " g/cm³" + (_growShape > 0 ? " (cell)" : "");
            GrowElapsed = sw.Elapsed.TotalSeconds;
            GrowLog = (failures.Count > 0 ? string.Join("\n", failures) + $"\nused seed {used} instead\n" : "") + report + $"\nbuilt in {sw.Elapsed.TotalSeconds:F2} s" + densityNote;
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
            ++_growLiveTicket;
            var old = GrowLiveDoc;
            GrowLiveDoc = null;   // the finished cell is the document now
            old?.Dispose();
        }
        if (grown) await AssignForBuilder(pack: false);
    }

    public void CancelGrow() => _growCancel?.Cancel();

    // ---------------------------------------------------------------- Relax
    public static readonly string[] Minimisers = ["Steepest descent", "Polak–Ribière conjugate gradient", "L-BFGS (m = 10)", "FIRE"];
    private int _relaxMethod = 2, _relaxIterations = 5000;
    private double _relaxFtol = 0.5, _relaxDensity = 1.05, _relaxStep = 0.06, _relaxPressure = 1.0, _relaxCutoff = 10.0;
    private bool _relaxPushoff = true, _relaxCompress = true, _relaxBox, _relaxCoulomb = true, _relaxing;
    private string _relaxLog = "Minimises the chosen structure with its assigned force field: capped-force push-off, " +
                               "compression to a target density, then minimisation to the force tolerance.";
    private string _fieldInfo = "";
    private CancellationTokenSource? _relaxCancel;
    private readonly List<(double X, double Y)> _relaxEnergy = new(), _relaxForce = new();

    public int RelaxMethod { get => _relaxMethod; set => Set(ref _relaxMethod, value); }
    public bool RelaxPushoff { get => _relaxPushoff; set { if (Set(ref _relaxPushoff, value)) Raise(nameof(RelaxPushoffMdEnabled)); } }
    // push-off by MD first, the LJ force cap ramped (Auhl et al. 2003): the board's "Force cap" and "λ ramp over … ps"
    private bool _relaxPushoffMd;
    private decimal _relaxCap = 500, _relaxRampPs = 20, _relaxPushoffT = 300;
    public bool RelaxPushoffMd { get => _relaxPushoffMd; set { if (Set(ref _relaxPushoffMd, value)) Raise(nameof(RelaxPushoffMdEnabled)); } }
    public bool RelaxPushoffMdEnabled => _relaxPushoff && _relaxPushoffMd;
    public decimal RelaxCap { get => _relaxCap; set { if (Set(ref _relaxCap, Math.Clamp(value, 5, 5000))) Raise(nameof(RelaxPushoffText)); } }
    public decimal RelaxRampPs { get => _relaxRampPs; set { if (Set(ref _relaxRampPs, Math.Clamp(value, 1, 1000))) Raise(nameof(RelaxPushoffText)); } }
    public decimal RelaxPushoffT { get => _relaxPushoffT; set => Set(ref _relaxPushoffT, Math.Clamp(value, 10, 2000)); }
    public string RelaxPushoffText => $"Auhl et al., J. Chem. Phys. 119, 12718 (2003): LJ forces capped, the cap raised to {_relaxCap:0} kcal/mol/Å" +
                                      $" (minimisation stages 5 → 20 → 100 → {_relaxCap:0}); with MD first, NVT with the cap raised over {_relaxRampPs:0.#} ps";
    public bool RelaxCompress { get => _relaxCompress; set => Set(ref _relaxCompress, value); }
    public bool RelaxBox { get => _relaxBox; set => Set(ref _relaxBox, value); }
    // how the box relaxes: every axis together (the volume), each on its own, only z (a film's thickness), only x and y
    public static readonly string[] RelaxBoxModes = ["Isotropic (the volume)", "Each axis on its own", "Only z (a film or slab)", "Only x and y (fixed thickness)"];
    private int _relaxBoxMode;
    public int RelaxBoxMode { get => _relaxBoxMode; set => Set(ref _relaxBoxMode, Math.Clamp(value, 0, 3)); }
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
    /// <summary>Relax compresses a grown cell towards a target density: when data/reference/polymers.json has the grown
    /// polymer ("Polyethylene" → "Polyethylene, amorphous"), its measured density (mid-range, 0.01 g/cm³) becomes the
    /// target, so a PE cell is not squeezed to polystyrene's 1.05. Returns a line for the Grow log, or "".</summary>
    private string SuggestRelaxDensity(string polymer)
    {
        try
        {
            if (Paths.References is not { } f || string.IsNullOrWhiteSpace(polymer)) return "";
            var mats = JsonNode.Parse(File.ReadAllText(f))?["materials"] as JsonArray;
            var m = mats?.OfType<JsonObject>().FirstOrDefault(x =>
            {
                var n = (string?)x["name"] ?? "";
                return n.Equals(polymer, StringComparison.OrdinalIgnoreCase) || n.StartsWith(polymer + ",", StringComparison.OrdinalIgnoreCase);
            });
            if (m?["values"]?["density"] is not JsonObject d || d["lo"] is null || d["hi"] is null) return "";
            var mid = Math.Round(((double)d["lo"]! + (double)d["hi"]!) / 2, 2);
            _relaxDensity = mid;
            Raise(nameof(RelaxDensityD));
            return string.Format(CultureInfo.InvariantCulture, "\nRelax target density set to {0:0.00} g/cm³ ({1}, {2})", mid, (string?)m["name"], (string?)d["source"]);
        }
        catch { return ""; }
    }

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

    // stop criteria (design/boards/Relax): the relative energy change and, with box relaxation, the stress left
    public static readonly string[] RelaxEtols = ["1e-6", "1e-8", "1e-10", "1e-12"];
    private int _relaxEtol = 1;
    private double _relaxPressureTol = 100;
    public int RelaxEtolIndex { get => _relaxEtol; set => Set(ref _relaxEtol, Math.Clamp(value, 0, 3)); }
    public decimal RelaxPressureTolD { get => (decimal)_relaxPressureTol; set { _relaxPressureTol = Math.Clamp((double)value, 1, 100000); Raise(); } }

    // atoms held in place besides the held molecule (Relax board: Fixed atoms · Add from selection)
    private int _fixedCount;
    public bool HasFixedAtoms => _fixedCount > 0;
    public string FixedAtomsText => _fixedCount == 0 ? "Fixed atoms: none (besides a held molecule)"
        : $"Fixed atoms: {_fixedCount:N0} held " + (FixedAxes == 7 ? "in place" : "along " + string.Join(", ", new[] { "x", "y", "z" }.Where((_, k) => (FixedAxes >> k & 1) == 1)) + " only");
    // which coordinates of the fixed atoms are held (x 1, y 2, z 4): a substrate that slides in its plane holds z only
    private int FixedAxes => _doc?.FixedAxes() ?? 7;
    private void SetFixedAxis(int bit, bool on)
    {
        if (_doc == null) return;
        var a = on ? FixedAxes | bit : FixedAxes & ~bit;
        if (a == 0) { Status = "At least one axis stays held (free the atoms instead)"; Raise(nameof(FixedX)); Raise(nameof(FixedY)); Raise(nameof(FixedZ)); return; }
        _doc.SetFixedAxes(a);
        Raise(nameof(FixedX)); Raise(nameof(FixedY)); Raise(nameof(FixedZ)); Raise(nameof(FixedAtomsText));
    }
    public bool FixedX { get => (FixedAxes & 1) != 0; set => SetFixedAxis(1, value); }
    public bool FixedY { get => (FixedAxes & 2) != 0; set => SetFixedAxis(2, value); }
    public bool FixedZ { get => (FixedAxes & 4) != 0; set => SetFixedAxis(4, value); }
    public void HoldSelection()
    {
        if (_doc == null) return;
        var sel = new HashSet<int>(_selection);
        try
        {
            if (System.Text.Json.Nodes.JsonNode.Parse(_doc.SelectionJson())?["indices"] is System.Text.Json.Nodes.JsonArray a)
                foreach (var x in a) sel.Add((int)((double?)x ?? -1));
        }
        catch { }
        sel.Remove(-1);
        if (sel.Count == 0) { Status = "Select atoms first (click, ⇧-click, lasso or Select by query), then Hold selection"; return; }
        sel.UnionWith(_doc.FixedAtoms());
        _fixedCount = _doc.SetFixedAtoms(sel);
        Raise(nameof(HasFixedAtoms)); Raise(nameof(FixedAtomsText));
        Status = $"{_fixedCount:N0} atoms held in place in Relax, Dynamics and Equilibrate (a freeze group in GROMACS files)";
    }
    public void FreeFixedAtoms()
    {
        if (_doc == null) return;
        _fixedCount = _doc.SetFixedAtoms(Array.Empty<int>());
        Raise(nameof(HasFixedAtoms)); Raise(nameof(FixedAtomsText));
        Status = "No atoms fixed (a held molecule stays held)";
    }
    private void SyncFixed()
    {
        try { _mdRigid = _doc?.RigidMolecules() ?? ""; } catch { _mdRigid = ""; }
        Raise(nameof(MdRigid));
        MdRigidNote = _mdRigid.Length == 0 ? "" : "rigid in LAMMPS (fix rigid/nvt/small); CAPS's own run moves their atoms freely";
        _fixedCount = _doc?.FixedAtoms().Length ?? 0;
        Raise(nameof(HasFixedAtoms)); Raise(nameof(FixedAtomsText)); Raise(nameof(FixedX)); Raise(nameof(FixedY)); Raise(nameof(FixedZ));
    }

    /// <summary>The Relax settings as the core takes them (also captured when a run is queued).</summary>
    private CapsRelaxOpts RelaxOptions() => new()
    {
        Method = _relaxMethod, Ftol = _relaxFtol, MaxIterations = _relaxIterations,
        TargetDensity = _relaxCompress ? _relaxDensity : 0, CompressStep = _relaxStep,
        Pushoff = _relaxPushoff ? 1 : 0, RelaxBox = _relaxBox ? 1 : 0, Pressure = _relaxPressure,
        BoxAnisotropic = _relaxBoxMode > 0 ? 1 : 0, BoxAxes = _relaxBoxMode switch { 2 => 4, 3 => 3, _ => 7 },
        PushoffRampPs = _relaxPushoff && _relaxPushoffMd ? (double)_relaxRampPs : 0, PushoffCap = (double)_relaxCap, PushoffTemperature = (double)_relaxPushoffT,
        Cutoff = _relaxCutoff, Coulomb = _relaxCoulomb ? 1 : 0,
        Etol = double.Parse(RelaxEtols[_relaxEtol], CultureInfo.InvariantCulture), PressureTol = _relaxPressureTol,
    };

    public async Task Relax() => await Relax(null);

    private async Task Relax(CapsRelaxOpts? preset)
    {
        if (_doc == null || !Idle || BlockedByField("Relax")) return;
        PrepareRunTarget("minimised");
        var doc = _doc!;
        ApplyRestraints();
        Relaxing = true;
        IsPlaying = false;
        _relaxCancel = new CancellationTokenSource();
        var token = _relaxCancel.Token;
        var o = preset ?? RelaxOptions();
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
                    WaitIfPaused(token);
                    return !token.IsCancellationRequested;
                });
                Publish(null);
                return r;
            });
            sw.Stop();
            finished = true;
            Record(string.Format(inv, "doc.relax(ftol={0}, method=\"{1}\", max_iterations={2})", _relaxFtol, _relaxMethod switch { 0 => "sd", 1 => "cg", 3 => "fire", _ => "lbfgs" }, _relaxIterations));
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

    public void CancelRelax() { _relaxCancel?.Cancel(); ResumeRun(); }

    // ---------------------------------------------------------------- Dynamics
    public static readonly string[] Ensembles = ["NVE (no thermostat)", "NVT", "NPT (isotropic)", "NPH (Berendsen, no thermostat)"];
    public static readonly string[] Thermostats = ["Bussi velocity rescaling", "Langevin (BAOAB)", "Nosé–Hoover chain (as LAMMPS fix nvt)"];
    public static readonly string[] Barostats = ["Stochastic cell rescaling", "Berendsen (early relaxation only)", "MTK · Nosé–Hoover (as LAMMPS fix npt)"];
    private int _mdEnsemble = 1, _mdThermostat, _mdBarostat, _mdSeed = 1;
    private double _mdDt = 1.0, _mdTemp = 300, _mdTauT = 100, _mdPressure = 1, _mdTauP = 1000;
    private long _mdSteps = 20000;
    private int _mdFrameEvery = 1000;
    private bool _mdNewVelocities, _mdRunning;
    private string _mdLog = "Runs molecular dynamics on the structure in the viewer with the Relax force field (GAFF, C and H). " +
                            "Relax first: dynamics from a grown cell with overlaps is unstable.";
    private CancellationTokenSource? _mdCancel;
    private readonly List<CapsThermo> _thermo = new();

    public int MdEnsemble
    {
        get => _mdEnsemble;
        set
        {
            if (!Set(ref _mdEnsemble, value)) return;
            if (value == 3) MdBarostat = 1;   // NPH shows the barostat it runs
            Raise(nameof(MdHasThermostat)); Raise(nameof(MdHasBarostat)); Raise(nameof(MdCanChooseBarostat)); Raise(nameof(MdCitation)); Raise(nameof(MdIntegratorText)); Raise(nameof(MdRespaAllowed));
            if (!MdRespaAllowed) MdRespa = 0;
            RefreshPreflight();
        }
    }
    /// <summary>The papers behind the chosen thermostat and barostat.</summary>
    public string MdCitation => string.Join("; ", new[]
    {
        MdHasThermostat ? (_mdThermostat == 1 ? "Leimkuhler & Matthews, Appl. Math. Res. Express 2013, 34 (BAOAB)" : _mdThermostat == 2 ? "Martyna, Klein & Tuckerman, J. Chem. Phys. 97, 2635 (1992)" : "Bussi et al., J. Chem. Phys. 126, 014101 (2007)") : null,
        MdHasBarostat ? (_mdEnsemble == 3 || _mdBarostat == 1 ? "Berendsen et al., J. Chem. Phys. 81, 3684 (1984)" : _mdBarostat == 2 ? "Martyna, Tobias & Klein, J. Chem. Phys. 101, 4177 (1994)" : "Bernetti & Bussi, J. Chem. Phys. 153, 114107 (2020)") : null,
        _mdEnsemble == 0 ? "Swope et al., J. Chem. Phys. 76, 637 (1982)" : null,
        _mdRespa > 0 ? "r-RESPA: Tuckerman, Berne & Martyna, J. Chem. Phys. 97, 1990 (1992)" : null,
        _mdConstraints > 0 ? (_mdConstraintSolver == 1 ? "LINCS: Hess et al., J. Comput. Chem. 18, 1463 (1997); RATTLE: Andersen, J. Comput. Phys. 52, 24 (1983)"
                                                       : "SHAKE/RATTLE: Ryckaert et al., J. Comput. Phys. 23, 327 (1977); Andersen, J. Comput. Phys. 52, 24 (1983)") : null,
    }.Where(x => x != null));
    public bool MdHasThermostat => _mdEnsemble is 1 or 2;
    public bool MdHasBarostat => _mdEnsemble is 2 or 3;
    /// <summary>NPH runs the Berendsen barostat: stochastic cell rescaling needs a thermostat.</summary>
    public bool MdCanChooseBarostat => _mdEnsemble == 2;
    public int MdThermostat
    {
        get => _mdThermostat;
        set
        {
            if (!Set(ref _mdThermostat, value)) return;
            if (value is 1 or 2) MdRespa = 0;
            if (value != 2 && _mdBarostat == 2) MdBarostat = 0;   // MTK is the pressure half of Nosé–Hoover
            Raise(nameof(MdCitation)); Raise(nameof(MdIntegratorText)); Raise(nameof(MdRespaAllowed)); RefreshPreflight();
        }
    }
    // r-RESPA: 0 off, 1 two inner steps, 2 four (bonded forces every Δt/2 or Δt/4); with Bussi or no thermostat
    public static readonly string[] RespaChoices = ["r-RESPA · off", "r-RESPA · bonded ×2", "r-RESPA · bonded ×4"];
    private int _mdRespa;
    public int MdRespa { get => _mdRespa; set { if (Set(ref _mdRespa, Math.Clamp(value, 0, 2))) { Raise(nameof(MdCitation)); Raise(nameof(MdIntegratorText)); Raise(nameof(MdConstraintsAllowed)); RefreshPreflight(); } } }
    public bool MdRespaAllowed => !(MdHasThermostat && _mdThermostat is 1 or 2) && _mdConstraints == 0;
    // bond constraints (SHAKE or LINCS, RATTLE velocities): 0 none, 1 bonds to hydrogen with rigid water, 2 every bond; r-RESPA is the alternative
    public static readonly string[] ConstraintChoices = ["None: every bond flexible", "Bonds to hydrogen, rigid water", "All bonds"];
    public static readonly string[] ConstraintSolvers = ["SHAKE / RATTLE (as LAMMPS)", "LINCS (as GROMACS)"];
    private int _mdConstraintSolver;
    public int MdConstraintSolver { get => _mdConstraintSolver; set { if (Set(ref _mdConstraintSolver, Math.Clamp(value, 0, 1))) { Raise(nameof(MdCitation)); Raise(nameof(MdIntegratorText)); RefreshPreflight(); } } }
    public bool MdHasConstraints => _mdConstraints > 0;
    private int _mdConstraints;
    public int MdConstraints
    {
        get => _mdConstraints;
        set
        {
            if (!Set(ref _mdConstraints, Math.Clamp(value, 0, 2))) return;
            if (value > 0) { MdRespa = 0; if (_mdDt < 2) MdDtD = 2; if (_mdBarostat == 2) MdBarostat = 0; }   // what constraints are for: 2 fs steps; MTK runs without them
            Raise(nameof(MdRespaAllowed)); Raise(nameof(MdCitation)); Raise(nameof(MdIntegratorText)); Raise(nameof(MdConstraintsAllowed)); Raise(nameof(MdHasConstraints));
            RefreshPreflight();
        }
    }
    public bool MdConstraintsAllowed => _mdRespa == 0;
    private int RespaSteps => _mdRespa == 0 ? 1 : _mdRespa == 1 ? 2 : 4;
    public int MdBarostat
    {
        get => _mdBarostat;
        set
        {
            if (!Set(ref _mdBarostat, value)) return;
            if (value == 2) { if (_mdThermostat != 2) MdThermostat = 2; if (_mdConstraints > 0) MdConstraints = 0; }   // MTK: Nosé–Hoover, no constraints
            Raise(nameof(MdCitation)); Raise(nameof(MdIntegratorText)); RefreshPreflight();
        }
    }
    public bool MdNewVelocities { get => _mdNewVelocities; set => Set(ref _mdNewVelocities, value); }
    public bool MdRunning { get => _mdRunning; private set { if (Set(ref _mdRunning, value)) RaiseBusy(); } }
    public bool CanRun => _doc != null && Idle;
    public string MdLog { get => _mdLog; private set => Set(ref _mdLog, value); }
    public IReadOnlyList<CapsThermo> Thermo => _thermo;
    public event Action? ThermoChanged;
    // ---- live view of a running MD or equilibration: the positions so far, about four times a second
    private CapsDocument? _runLiveDoc;
    private string _runLiveText = "";
    private int _runLiveTicket;
    public CapsDocument? RunLiveDoc { get => _runLiveDoc; private set { if (Set(ref _runLiveDoc, value)) Raise(nameof(RunLiveShown)); } }
    public bool RunLiveShown => _runLiveDoc != null;
    public string RunLiveText { get => _runLiveText; private set => Set(ref _runLiveText, value); }
    /// <summary>A receiver for caps_set_live snapshots on the run's thread: shown on the UI thread, the previous one disposed.</summary>
    private Action<CapsDocument, string> LiveReceiver(string what)
    {
        var ticket = ++_runLiveTicket;
        var inv = CultureInfo.InvariantCulture;
        return (live, stats) => Avalonia.Threading.Dispatcher.UIThread.Post(() =>
        {
            if (ticket != _runLiveTicket) { live.Dispose(); return; }
            var old = RunLiveDoc;
            RunLiveDoc = live;
            old?.Dispose();
            try
            {
                var j = System.Text.Json.Nodes.JsonNode.Parse(stats);
                var t = (double?)j?["time_ps"] ?? 0;
                var rho = (double?)j?["density"];
                RunLiveText = string.Format(inv, "{0} · {1:F2} ps · {2:N0} atoms", what, t, (double?)j?["atoms"] ?? 0) +
                              (rho is double r ? string.Format(inv, " · ρ {0:F4} g/cm³", r) : "");
            }
            catch { RunLiveText = what; }
        });
    }
    /// <summary>The live view keeps the last snapshot of the run until the next run starts.</summary>
    private void StartLive() { _runLiveTicket++; var old = RunLiveDoc; RunLiveDoc = null; old?.Dispose(); RunLiveText = ""; }

    public decimal? MdDtD { get => (decimal)_mdDt; set { _mdDt = Math.Clamp((double)(value ?? 1m), 0.1, 5); Raise(); Raise(nameof(MdEstimate)); } }
    public decimal? MdStepsD { get => _mdSteps; set { _mdSteps = Math.Clamp((long)(value ?? 20000), 0, 1_000_000_000); Raise(); Raise(nameof(MdEstimate)); } }
    public decimal? MdTempD { get => (decimal)_mdTemp; set { _mdTemp = Math.Clamp((double)(value ?? 300m), 1, 5000); Raise(); } }
    public decimal? MdTauTD { get => (decimal)_mdTauT; set { _mdTauT = Math.Clamp((double)(value ?? 100m), 1, 1e6); Raise(); } }
    public decimal? MdPressureD { get => (decimal)_mdPressure; set { _mdPressure = (double)(value ?? 1m); Raise(); } }
    public decimal? MdTauPD { get => (decimal)_mdTauP; set { _mdTauP = Math.Clamp((double)(value ?? 1000m), 10, 1e7); Raise(); } }
    public decimal? MdFrameEveryD { get => _mdFrameEvery; set { _mdFrameEvery = Math.Clamp((int)(value ?? 1000), 1, 1_000_000); Raise(); Raise(nameof(MdEstimate)); } }
    public decimal? MdSeedD { get => _mdSeed; set { if (value == null) return; _mdSeed = Math.Max(1, (int)value); if (_mdSeedC != null && _mdSeedC.Value != _mdSeed) _mdSeedC.Value = _mdSeed; Raise(); } }
    // ---- pre-flight and export (Dynamics board)
    public ObservableCollection<CheckRow> MdPreflight { get; } = new();
    private string _mdPreflightSummary = "", _mdDeck = "";
    public string MdPreflightSummary { get => _mdPreflightSummary; private set => Set(ref _mdPreflightSummary, value); }
    /// <summary>fix press/berendsen's coupling as this run's: iso, aniso, one axis or two (the others fixed). LAMMPS has no
    /// Berendsen for the tilts: the full shape writes aniso and says how fix npt tri would do it.</summary>
    private string BerendsenCoupling(IFormatProvider inv)
    {
        string Ax(string a) => string.Format(inv, "{0} {1:0.##} {1:0.##} {2:0.##}", a, _mdPressure, _mdTauP);
        return _mdCoupling switch
        {
            1 => Ax("aniso"),
            2 => Ax("z"),
            3 => Ax("x") + " " + Ax("y") + " couple none",
            4 => Ax("aniso") + string.Format(inv, "\n# full shape: CAPS also relaxes the tilts (Berendsen); in LAMMPS, fix npt ... tri {0:0.##} {0:0.##} {1:0.##} (Nosé–Hoover) does", _mdPressure, _mdTauP),
            _ => Ax("iso"),
        };
    }

    /// <summary>The LAMMPS input for this run (setup from the core, ensemble lines from the settings).</summary>
    public string MdDeck { get => _mdDeck; private set => Set(ref _mdDeck, value); }

    private Avalonia.Threading.DispatcherTimer? _preflightTimer;
    private int _preflightTicket;

    /// <summary>The Dynamics pre-flight and the LAMMPS / GROMACS input for these settings: asked for often (every setting
    /// change), done once the settings rest, off the UI thread (typing and writing the input of a 180 000-atom cell
    /// takes seconds), and never while a run holds the document.</summary>
    public void RefreshPreflight()
    {
        if (_doc == null) { MdPreflight.Clear(); MdPreflightSummary = ""; MdDeck = ""; return; }
        if (Busy) return;
        if (_preflightTimer == null)
        {
            _preflightTimer = new Avalonia.Threading.DispatcherTimer { Interval = TimeSpan.FromMilliseconds(250) };
            _preflightTimer.Tick += (_, _) => { _preflightTimer!.Stop(); _ = RefreshPreflightNow(); };
        }
        _preflightTimer.Stop();
        _preflightTimer.Start();
    }

    /// <summary>The pre-flight at once (tests).</summary>
    public Task PreflightNow() { _preflightTimer?.Stop(); return RefreshPreflightNow(); }

    private async Task RefreshPreflightNow()
    {
        var doc = _doc;
        if (doc == null || Busy) return;
        var ticket = ++_preflightTicket;
        var inv = CultureInfo.InvariantCulture;
        var s = doc.Summary();
        var assigned = Field.Assigned;
        var gromacs = _mdGromacs;
        if (MdPreflight.Count == 0) MdPreflightSummary = "checking…";
        if (MdDeck.Length == 0) MdDeck = gromacs ? "; writing the GROMACS run parameters…" : "# writing the LAMMPS input…";
        var (ff, deck) = await Task.Run(() =>
        {
            var f = "";
            if (!assigned)
                try { f = doc.FieldInfo(); } catch (Exception e) { f = "error: " + e.Message; }
            string d;
            try { d = gromacs ? GromacsDeck(doc) : LammpsDeck(doc); }
            catch (Exception e) { d = gromacs ? "; cannot write the GROMACS files: " + e.Message : "# cannot write the LAMMPS input: " + e.Message; }
            return (f, d);
        });
        if (ticket != _preflightTicket || !ReferenceEquals(doc, _doc)) return;   // superseded
        MdPreflight.Clear();
        // force field
        if (assigned)
            MdPreflight.Add(new CheckRow(Field.Complete ? "All atoms typed, no missing parameters" : "The Field assignment is incomplete: runs are blocked",
                Field.Complete ? "ok" : "fail"));
        else
        {
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
        if (_mdEnsemble >= 2 && s.CellValid == 0) MdPreflight.Add(new CheckRow((_mdEnsemble == 2 ? "NPT" : "NPH") + " needs a periodic cell", "fail"));
        if (_mdEnsemble == 3)
            MdPreflight.Add(new CheckRow("NPH: Berendsen barostat without a thermostat; Berendsen scaling does not conserve the enthalpy exactly, so watch the drift", "check"));
        // time step
        var inner = _mdDt / RespaSteps;   // the step the bonded forces (C–H stretches) see
        if (_mdConstraints > 0)
            MdPreflight.Add(new CheckRow(string.Format(inv, "Δt {0:0.##} fs with {1} held at their lengths ({2})", _mdDt,
                    _mdConstraints == 1 ? "bonds to hydrogen" : "all bonds", _mdConstraintSolver == 1 ? "LINCS" : "SHAKE/RATTLE"),
                _mdDt <= (_mdConstraints == 1 ? 2.0 : 2.5) ? "ok" : _mdDt <= 3.0 ? "check" : "fail"));
        else
            MdPreflight.Add(new CheckRow(RespaSteps > 1
                    ? string.Format(inv, "Δt {0:0.##} fs, bonded forces every {1:0.###} fs (r-RESPA), no bond constraints", _mdDt, inner)
                    : string.Format(inv, "Δt {0:0.##} fs with hydrogens, no bond constraints", _mdDt),
                inner <= 1.0 && _mdDt <= 4.0 ? "ok" : inner <= 2.0 ? "check" : "fail"));
        // velocities
        MdPreflight.Add(new CheckRow(_mdNewVelocities ? string.Format(inv, "New velocities at {0:0} K (seed {1})", _mdTemp, _mdSeed) : "Velocities from the structure, or drawn at the target if it has none", "ok"));
        var fails = MdPreflight.Count(r => r.State == "fail");
        var checks = MdPreflight.Count(r => r.State == "check");
        MdPreflightSummary = $"{MdPreflight.Count - fails - checks} / {MdPreflight.Count} ok";
        MdDeck = deck;
    }

    /// <summary>The LAMMPS input for this run (setup from the core, ensemble lines from the settings).</summary>
    private string LammpsDeck(CapsDocument doc)
    {
        var inv = CultureInfo.InvariantCulture;
        var setup = doc.LammpsInput("system.data");
        var steps = _mdSteps;
        var ens = _mdEnsemble switch
        {
            0 => "fix 1 all nve",
            1 => _mdThermostat == 1 ? string.Format(inv, "fix 1 all nve\nfix 2 all langevin {0:0.##} {0:0.##} {1:0.##} {2}", _mdTemp, _mdTauT, _mdSeed + 1)
               : _mdThermostat == 2 ? string.Format(inv, "fix 1 all nvt temp {0:0.##} {0:0.##} {1:0.##}", _mdTemp, _mdTauT)
                                     : string.Format(inv, "fix 1 all nve\nfix 2 all temp/csvr {0:0.##} {0:0.##} {1:0.##} {2}", _mdTemp, _mdTauT, _mdSeed + 1),
            2 when _mdBarostat == 2 => string.Format(inv, "fix 1 all npt temp {0:0.##} {0:0.##} {1:0.##} iso {2:0.##} {2:0.##} {3:0.##}", _mdTemp, _mdTauT, _mdPressure, _mdTauP),
            3 => string.Format(inv, "fix 1 all nve\nfix 3 all press/berendsen {0} modulus 22222", BerendsenCoupling(inv)),
            _ => (_mdThermostat == 2 ? string.Format(inv, "fix 1 all nvt temp {0:0.##} {0:0.##} {1:0.##}", _mdTemp, _mdTauT)
                  : _mdThermostat == 1 ? string.Format(inv, "fix 1 all nve\nfix 2 all langevin {0:0.##} {0:0.##} {1:0.##} {2}", _mdTemp, _mdTauT, _mdSeed + 1)
                  : string.Format(inv, "fix 1 all nve\nfix 2 all temp/csvr {0:0.##} {0:0.##} {1:0.##} {2}", _mdTemp, _mdTauT, _mdSeed + 1)) +
                 string.Format(inv, "\nfix 3 all press/berendsen {0} modulus 22222", BerendsenCoupling(inv)),   // modulus 1/β for β = 4.5e-5 atm⁻¹, as CAPS's barostat
        };
        if (_mdRigid.Length > 0 && setup.Contains("group           rigid molecule"))
        {   // the bodies first (LAMMPS wants rigid fixes before any box-changing fix), the rest integrated on their own
            var body = _mdEnsemble == 0 ? "rigid/small molecule" : string.Format(inv, "rigid/nvt/small molecule temp {0:0.##} {0:0.##} {1:0.##}", _mdTemp, _mdTauT);
            ens = "group mobile subtract all rigid\nfix 0 rigid " + body + "\n" + ens.Replace("fix 1 all ", "fix 1 mobile ").Replace("fix 2 all ", "fix 2 mobile ");
        }
        return "# LAMMPS input written by CAPS Studio: the same force field and settings as this Dynamics run\n" + setup +
               (_mdNewVelocities ? string.Format(inv, "velocity all create {0:0.##} {1} mom yes rot yes dist gaussian\n", _mdTemp, _mdSeed) : "") +
               (RespaSteps > 1 ? $"run_style respa 2 {RespaSteps} bond 1 angle 1 dihedral 1 improper 1 pair 2 kspace 2\n" : "") +
               (_mdConstraints > 0 ? doc.LammpsShake(_mdConstraints) : "") +
               (_mdFieldOn ? string.Format(inv, "fix efield all efield {0:0.######} {1:0.######} {2:0.######}\n", _mdEx, _mdEy, _mdEz) : "") +
               string.Format(inv, "timestep {0:0.###}\n{1}\nthermo {2}\ndump d all custom {3} traj.lammpstrj id mol type xu yu zu\nrun {4}\n",
                   _mdDt, ens, Math.Max(1, _mdFrameEvery / 10), _mdFrameEvery, steps);
    }

    public string MdEstimate => string.Format(CultureInfo.InvariantCulture, "{0:0.###} ps · {1:N0} frames recorded",
        _mdSteps * _mdDt / 1000, _mdSteps / Math.Max(1, _mdFrameEvery) + 1);

    public async Task RunMd() => await RunMd(null);

    // checkpoints (design/boards/FailedJob): a run that failed or was stopped keeps its last checkpoint
    private bool _mdCanContinue;
    private string _mdCheckpointText = "";
    private long _mdCheckpointStep, _mdCheckpointSteps;
    public bool MdCanContinue { get => _mdCanContinue; private set => Set(ref _mdCanContinue, value); }
    public string MdCheckpointText { get => _mdCheckpointText; private set => Set(ref _mdCheckpointText, value); }
    public string MdContinueLabel => string.Format(CultureInfo.InvariantCulture, "Continue from step {0:N0}", _mdCheckpointStep);

    private void RefreshMdCheckpoint(CapsDocument doc)
    {
        MdCanContinue = false;
        MdCheckpointText = "";
        try
        {
            var c = System.Text.Json.Nodes.JsonNode.Parse(doc.Checkpoint("info"))!;
            if (c["has"]?.GetValue<bool>() != true || (string?)c["kind"] != "md" || (string?)c["ended"] == "finished") return;
            _mdCheckpointStep = (long)c["step"]!.GetValue<double>();
            _mdCheckpointSteps = (long)c["steps"]!.GetValue<double>();
            MdCheckpointText = string.Format(CultureInfo.InvariantCulture, "The run {0} after its checkpoint at step {1:N0} of {2:N0} ({3:0.##} ps): that state and its velocities are kept",
                (string?)c["ended"] == "failed" ? "failed" : "was stopped", _mdCheckpointStep, _mdCheckpointSteps, c["time_ps"]!.GetValue<double>());
            MdCanContinue = _mdCheckpointSteps > _mdCheckpointStep;
            Raise(nameof(MdContinueLabel));
        }
        catch { /* no checkpoint */ }
    }

    /// <summary>Continue from the checkpoint: its state becomes the last frame and the remaining steps run from it with
    /// the same velocities (a run that failed can be continued with a smaller time step).</summary>
    public async Task ContinueMd()
    {
        if (_doc == null || !Idle || !_mdCanContinue) return;
        try { _doc.Checkpoint("restore"); }
        catch (Exception e) { Status = "Could not continue: " + e.Message; MdCanContinue = false; return; }
        await RunMd((_mdCheckpointStep, _mdCheckpointSteps - _mdCheckpointStep));
    }

    // pressure coupling per axis (Berendsen: each axis from its own diagonal pressure) and the checkpoint interval
    public string MdIntegratorText => (MdHasThermostat && _mdThermostat == 1 ? "BAOAB (Langevin)" : "Velocity Verlet") + (RespaSteps > 1 ? " · r-RESPA" : "");
    public static readonly string[] MdCouplings = ["Isotropic", "Each axis on its own (Berendsen)", "Only z (Berendsen)", "Only x and y (Berendsen)", "Full shape: lengths and tilts (Berendsen)"];
    private int _mdCoupling;
    private double _mdCheckpointPs;
    public int MdCoupling { get => _mdCoupling; set => Set(ref _mdCoupling, Math.Clamp(value, 0, 4)); }
    public decimal MdCheckpointPsD { get => (decimal)_mdCheckpointPs; set { _mdCheckpointPs = Math.Max(0, (double)value); Raise(); } }

    // rigid bodies in the LAMMPS deck (filler particles): molecule ranges, kept by the document for every LAMMPS export
    private string _mdRigid = "", _mdRigidNote = "";
    public string MdRigid
    {
        get => _mdRigid;
        set
        {
            if (!Set(ref _mdRigid, (value ?? "").Trim())) return;
            try { var n = _doc?.SetRigidMolecules(_mdRigid) ?? 0; MdRigidNote = n == 0 ? "" : $"{n} molecule{(n == 1 ? "" : "s")} rigid in LAMMPS (fix rigid/nvt/small); CAPS's own run moves their atoms freely"; }
            catch (Exception e) { MdRigidNote = e.Message; }
            RefreshPreflight();
        }
    }
    public string MdRigidNote { get => _mdRigidNote; private set => Set(ref _mdRigidNote, value); }

    // an external electric field (V/Å) on the partial charges: poling, field-driven ion transport, dielectric response
    private bool _mdFieldOn;
    private double _mdEx, _mdEy, _mdEz = 0.1;
    public bool MdFieldOn { get => _mdFieldOn; set { if (Set(ref _mdFieldOn, value)) Raise(nameof(MdFieldText)); } }
    public decimal MdExD { get => (decimal)_mdEx; set { _mdEx = (double)value; Raise(); Raise(nameof(MdFieldText)); } }
    public decimal MdEyD { get => (decimal)_mdEy; set { _mdEy = (double)value; Raise(); Raise(nameof(MdFieldText)); } }
    public decimal MdEzD { get => (decimal)_mdEz; set { _mdEz = (double)value; Raise(); Raise(nameof(MdFieldText)); } }
    /// <summary>|E| in V/Å and V/m, and the force it puts on one elementary charge.</summary>
    public string MdFieldText
    {
        get
        {
            var e = Math.Sqrt(_mdEx * _mdEx + _mdEy * _mdEy + _mdEz * _mdEz);
            return string.Format(CultureInfo.InvariantCulture, "|E| {0:0.###} V/Å = {1:0.##E+0} V/m · {2:0.##} kcal/mol/Å on a charge of 1 e", e, e * 1e10, e * 23.0605);
        }
    }

    /// <summary>The Dynamics settings as the core takes them (also captured when a run is queued).</summary>
    private CapsMdOpts MdOptions((long Offset, long Steps)? resume) => new()
    {
        BoxAnisotropic = _mdCoupling > 0 ? 1 : 0, BoxAxes = _mdCoupling switch { 2 => 4, 3 => 3, _ => 7 }, FullShape = _mdCoupling == 4 ? 1 : 0,
        EfieldX = _mdFieldOn ? _mdEx : 0, EfieldY = _mdFieldOn ? _mdEy : 0, EfieldZ = _mdFieldOn ? _mdEz : 0,
        CheckpointEvery = _mdCheckpointPs > 0 ? Math.Max(1, (long)Math.Round(_mdCheckpointPs * 1000 / Math.Max(0.01, _mdDt))) : 0,
        Dt = _mdDt, Steps = resume?.Steps ?? _mdSteps, Temperature = _mdTemp, StepOffset = resume?.Offset ?? 0,
        Thermostat = _mdEnsemble is 1 or 2 ? _mdThermostat + 1 : 0, TauT = _mdTauT,
        Barostat = _mdEnsemble == 2 ? _mdBarostat + 1 : _mdEnsemble == 3 ? 2 : 0, Pressure = _mdPressure, TauP = _mdTauP,   // NPH: Berendsen
        NewVelocities = resume == null && _mdNewVelocities ? 1 : 0, Seed = (ulong)(_mdSeed = MdSeedChoice.Take()),
        ThermoEvery = (int)Math.Clamp(_mdSteps / 400, 10, 1000), FrameEvery = _mdFrameEvery,
        Cutoff = _relaxCutoff, Coulomb = _relaxCoulomb ? 1 : 0, Tail = TailFlag, Respa = RespaSteps, Constraints = _mdConstraints, ConstraintAlgorithm = _mdConstraintSolver,
    };

    private async Task RunMd((long Offset, long Steps)? resume, CapsMdOpts? preset = null, int? ensemble = null)
    {
        if (_doc == null || !Idle || BlockedByField("Dynamics")) return;
        if (resume == null) PrepareRunTarget("MD");
        MdCanContinue = false;
        MdCheckpointText = "";
        var doc = _doc!;
        MdRunning = true;
        IsPlaying = false;
        _mdCancel = new CancellationTokenSource();
        var token = _mdCancel.Token;
        var o = preset ?? MdOptions(resume);
        var ens = ensemble ?? _mdEnsemble;
        _thermo.Clear();
        ThermoChanged?.Invoke();
        MdLog = "Starting…";
        Status = $"Running {Ensembles[ens]} dynamics on {Title}…";
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
            StartLive();
            var live = LiveReceiver(new[] { "NVE", "NVT", "NPT", "NPH" }[ens] + " dynamics");
            var report = await Task.Run(() =>
            {
                var r = doc.Md(o, (row, n) =>
                {
                    lock (rows) rows.Add(row);
                    if (sw.ElapsedMilliseconds - lastUi > 150)
                    {
                        lastUi = sw.ElapsedMilliseconds;
                        var frac = n > 0 ? (double)(row.Step - o.StepOffset) / n : 1;
                        var eta = frac > 0.001 ? sw.Elapsed.TotalSeconds * (1 - frac) / frac : 0;
                        Publish(string.Format(inv, "step {0:N0} of {1:N0} · {2:F2} ps · T {3:F1} K · P {4:F0} atm · ρ {5:F4} g/cm³ · E {6:F1} kcal/mol · {7:F0} s left",
                            row.Step, n + o.StepOffset, row.TimePs, row.Temperature, row.Pressure, row.Density, row.Total, eta));
                    }
                    WaitIfPaused(token);
                    return !token.IsCancellationRequested;
                }, live);
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
            RefreshMdCheckpoint(doc);
            var kept = MdCanContinue ? "\n" + MdCheckpointText + " — Continue runs the remaining steps from it." : "";
            MdLog = (cancelled ? "Cancelled; the structure is unchanged." : "Could not run dynamics.\n" + e.Message) + kept;
            Status = cancelled ? "Dynamics cancelled" + (MdCanContinue ? $" · {MdContinueLabel} is ready" : "") : "Could not run dynamics — see the Dynamics panel";
        }
        finally
        {
            MdRunning = false;
            ThermoChanged?.Invoke();
        }
    }

    public void CancelMd() { _mdCancel?.Cancel(); ResumeRun(); }

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
    public bool EqRunning { get => _eqRunning; private set { if (Set(ref _eqRunning, value)) { RaiseBusy(); Raise(nameof(EqCanDecide)); } } }
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
    public decimal? EqBlockD { get => (decimal)_eqBlock; set { _eqBlock = Math.Clamp((double)(value ?? 20m), 0.1, 1e6); Raise(); Raise(nameof(EqExtendLabel)); } }
    public decimal? EqMaxBlocksD { get => _eqMaxBlocks; set { _eqMaxBlocks = Math.Clamp((int)(value ?? 20), 1, 1000); Raise(); Raise(nameof(EqExtendLabel)); } }

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

    /// <summary>The protocol's stages on the time axis (from, to, ensemble · stage, NPT tinted): NVT stages hold the volume,
    /// so the density is flat there by construction.</summary>
    public (double From, double To, string Label, bool Shade)[] EqStageBands { get; private set; } = [];

    public async Task RunEquilibrate() => await RunEquilibrate(null, null);

    /// <summary>Extend (design/boards/Convergence): production NPT blocks again, from where the run ended, at the
    /// protocol's final conditions, until the checks pass or the blocks run out.</summary>
    public async Task ExtendEquilibrate()
    {
        var (t, p) = FinalNpt(_eqText);
        var text = string.Format(CultureInfo.InvariantCulture, "npt {0:0.###} ps T {1:0.##} P {2:0.###} atm   # extend: production blocks from the end of the last run", _eqBlock, t, p);
        _eqAccepted = false;
        await RunEquilibrate(text, true);
    }

    /// <summary>The temperature and pressure of a protocol's last NPT stage (the Dynamics temperature and 1 atm without one).</summary>
    internal (double T, double P) FinalNpt(string text)
    {
        var inv = CultureInfo.InvariantCulture;
        double t = _mdTemp, p = 1;
        foreach (var raw in text.Split('\n'))
        {
            var w = raw.Split('#')[0].Split(' ', StringSplitOptions.RemoveEmptyEntries);
            if (w.Length == 0 || !w[0].Equals("npt", StringComparison.OrdinalIgnoreCase)) continue;
            for (var k = 1; k + 1 < w.Length; k++)
            {
                if (w[k] == "T" && double.TryParse(w[k + 1], NumberStyles.Float, inv, out var tv)) t = tv;
                if (w[k] == "P" && double.TryParse(w[k + 1], NumberStyles.Float, inv, out var pv)) p = pv;
            }
        }
        return (t, p);
    }

    private bool _eqAccepted;
    public bool EqCanDecide => EqCriteria.Count > 0 && !_eqRunning && !_eqAccepted && EqCriteria.Any(r => r.State != "pass");
    public string EqExtendLabel => string.Format(CultureInfo.InvariantCulture, "Extend {0:0.#} ps", _eqBlock * _eqMaxBlocks);

    /// <summary>Accept now: the cell is taken as equilibrated with the criteria as they stand; the decision and the
    /// unmet criteria go into the provenance manifest.</summary>
    public void AcceptEquilibration()
    {
        if (_doc == null || EqCriteria.Count == 0) return;
        var unmet = EqCriteria.Where(r => r.State != "pass").Select(r => r.Title).ToArray();
        var met = EqCriteria.Count - unmet.Length;
        _doc.ProvenanceNote(new System.Text.Json.Nodes.JsonObject
        {
            ["engine"] = "equilibrate.accepted",
            ["summary"] = $"accepted as equilibrated with {met} of {EqCriteria.Count} criteria met",
            ["params"] = new System.Text.Json.Nodes.JsonObject
            {
                ["criteria met"] = $"{met} of {EqCriteria.Count}",
                ["not met"] = unmet.Length > 0 ? string.Join(", ", unmet) : "none",
                ["decided by"] = "the user (Accept now)",
            },
        }.ToJsonString());
        _eqAccepted = true;
        EqCriteriaText = $"{met} of {EqCriteria.Count} criteria met · accepted";
        Raise(nameof(EqCanDecide));
        Status = "Accepted as equilibrated; the decision and the unmet criteria are in the provenance";
    }

    /// <summary>The equilibration settings as the core takes them (also captured when a run is queued).</summary>
    private CapsEquilOpts EqOptions(bool? until) => new()
    {
        Dt = _mdDt, Thermostat = _mdThermostat + 1, Barostat = _mdBarostat + 1, TauT = _mdTauT, TauP = _mdTauP, Seed = (ulong)(_mdSeed = MdSeedChoice.Take()),
        Cutoff = _relaxCutoff, Coulomb = _relaxCoulomb ? 1 : 0, Tail = TailFlag,
        FramePs = 10, ThermoPs = 0.5, UntilConverged = (until ?? _eqUntil) ? 1 : 0, BlockPs = _eqBlock, MaxBlocks = _eqMaxBlocks, Constraints = _mdConstraints, ConstraintAlgorithm = _mdConstraintSolver,
    };

    private (string Text, CapsEquilOpts Opts, bool Target) _eqRun;   // the last run, for the macro recorder

    private async Task RunEquilibrate(string? protocol, bool? until, CapsEquilOpts? preset = null, double[]? presetTarget = null)
    {
        if (_doc == null || !Idle || BlockedByField("Equilibrate")) return;
        PrepareRunTarget("equilibrated");
        var doc = _doc!;
        EqRunning = true;
        IsPlaying = false;
        _eqCancel = new CancellationTokenSource();
        var token = _eqCancel.Token;
        var o = preset ?? EqOptions(until);
        // the internal-distance target curve, pinned for the run
        var target = presetTarget ?? EqTargetCurve();
        var pin = target.Length > 0 ? System.Runtime.InteropServices.GCHandle.Alloc(target, System.Runtime.InteropServices.GCHandleType.Pinned) : default;
        if (pin.IsAllocated) { o.InternalTarget = pin.AddrOfPinnedObject(); o.InternalTargetN = target.Length; }
        _thermo.Clear();
        ThermoChanged?.Invoke();
        EqLog = "Starting…";
        Status = $"Equilibrating {Title}: {Protocols[_eqProtocol]}…";
        var inv = CultureInfo.InvariantCulture;
        var sw = System.Diagnostics.Stopwatch.StartNew();
        var rows = new List<CapsThermo>();
        var lastUi = 0L;
        var finished = false;
        var text = protocol ?? _eqText;
        _eqRun = (text, o, target.Length > 0);
        if (protocol == null || preset != null) _eqAccepted = false;
        // each stage's ensemble from the protocol text (its first word: nvt, npt, …); production blocks are NPT
        var ens = text.Split('\n').Select(l => l.Split('#')[0].Trim()).Where(l => l.Length > 0)
                      .Select(l => l.Split(' ', StringSplitOptions.RemoveEmptyEntries)[0].ToUpperInvariant()).ToArray();
        var starts = new List<(int Stage, double T)>();
        EqStageBands = [];
        void Publish(string? line)
        {
            CapsThermo[] copy;
            (int Stage, double T)[] st0;
            lock (rows) { copy = rows.ToArray(); st0 = starts.ToArray(); }
            var end = copy.Length > 0 ? copy[^1].TimePs : 0;
            var bands = st0.Select((b, k) =>
            {
                var e = b.Stage - 1 < ens.Length ? ens[b.Stage - 1] : "NPT";
                return (b.T, k + 1 < st0.Length ? st0[k + 1].T : end, e, e is "NPT" or "NPH");
            }).ToArray();
            Avalonia.Threading.Dispatcher.UIThread.Post(() =>
            {
                _thermo.Clear();
                _thermo.AddRange(copy);
                EqStageBands = bands;
                if (line != null && !finished) EqLog = line;
                ThermoChanged?.Invoke();
            });
        }
        try
        {
            StartLive();
            var live = LiveReceiver("Equilibration");
            var (converged, report) = await Task.Run(() =>
            {
                var r = doc.Equilibrate(text, o, (st, n, label, row) =>
                {
                    lock (rows)
                    {
                        rows.Add(row);
                        if (starts.Count == 0 || starts[^1].Stage != st) starts.Add((st, starts.Count == 0 ? 0 : row.TimePs));
                    }
                    if (sw.ElapsedMilliseconds - lastUi > 150)
                    {
                        lastUi = sw.ElapsedMilliseconds;
                        Publish(string.Format(inv, "stage {0} of {1}: {2}\n{3:F2} ps · T {4:F1} K · P {5:F0} atm · ρ {6:F4} g/cm³ · {7:F0} s so far",
                            st, n, label, row.TimePs, row.Temperature, row.Pressure, row.Density, sw.Elapsed.TotalSeconds));
                    }
                    WaitIfPaused(token);
                    return !token.IsCancellationRequested;
                }, live);
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
            if (pin.IsAllocated) pin.Free();
            EqRunning = false;
            ThermoChanged?.Invoke();
        }
    }

    public void CancelEquilibrate() { _eqCancel?.Cancel(); ResumeRun(); }

    // ---- configurational-bias Monte Carlo (Equilibrate › Chain ends): Siepmann & Frenkel regrowth with the Field assignment
    private int _cbMoves = 2000, _cbTrials = 8, _cbTorsions = 4;
    private double _cbTemp = 300;
    private string _cbNote = "";
    public decimal CbMovesD { get => _cbMoves; set => Set(ref _cbMoves, (int)Math.Clamp(value, 10, 10_000_000)); }
    public decimal CbTrialsD { get => _cbTrials; set => Set(ref _cbTrials, (int)Math.Clamp(value, 2, 64)); }
    public decimal CbTorsionsD { get => _cbTorsions; set => Set(ref _cbTorsions, (int)Math.Clamp(value, 1, 12)); }
    public decimal CbTempD { get => (decimal)_cbTemp; set => Set(ref _cbTemp, Math.Clamp((double)value, 1, 5000)); }
    public string CbNote { get => _cbNote; private set => Set(ref _cbNote, value); }

    public async Task RunCbmc()
    {
        if (_doc == null || !Idle || BlockedByField("CBMC")) return;
        PrepareRunTarget("regrown");
        var doc = _doc!;
        EqRunning = true;
        IsPlaying = false;
        _eqCancel = new CancellationTokenSource();
        var token = _eqCancel.Token;
        var inv = CultureInfo.InvariantCulture;
        var o = new CapsCbmcOpts
        {
            Moves = _cbMoves, Trials = _cbTrials, MaxTorsions = _cbTorsions, Temperature = _cbTemp,
            Cutoff = Math.Min(9.0, _relaxCutoff), Coulomb = _relaxCoulomb ? 1 : 0, Seed = (ulong)(_mdSeed = MdSeedChoice.Take()),
        };
        CbNote = "Starting…";
        Status = $"Regrowing chain ends of {Title} (CBMC)…";
        var sw = System.Diagnostics.Stopwatch.StartNew();
        long lastUi = 0;
        try
        {
            var json = await Task.Run(() => doc.Cbmc(o, (done, n, acc) =>
            {
                if (sw.ElapsedMilliseconds - lastUi > 150)
                {
                    lastUi = sw.ElapsedMilliseconds;
                    var t = string.Format(inv, "{0} of {1} moves · {2} accepted ({3:F0} %) · {4:F0} s", done, n, acc, 100.0 * acc / Math.Max(1, done), sw.Elapsed.TotalSeconds);
                    Avalonia.Threading.Dispatcher.UIThread.Post(() => CbNote = t);
                }
                return !token.IsCancellationRequested;
            }));
            using var j = System.Text.Json.JsonDocument.Parse(json);
            var r = j.RootElement;
            double D(string k) => r.TryGetProperty(k, out var v) ? v.GetDouble() : 0;
            CbNote = string.Format(inv, "{0:F0} of {1:F0} regrowths accepted ({2:F1} %) · {3:F0} chain ends · ⟨R²⟩ {4:F0} → {5:F0} Å² · ΔE {6:F1} kcal/mol · {7:F1} s",
                                   D("accepted"), D("attempted"), 100 * D("acceptance"), D("chains"), D("r2_before"), D("r2_after"), D("energy_change"), D("seconds"));
            AfterRun(doc, " · CBMC");
            Status = "Chain ends regrown · " + CbNote;
        }
        catch (Exception e)
        {
            var cancelled = e.Message.Contains("cancelled");
            CbNote = cancelled ? "Cancelled; the structure is unchanged." : "Could not run CBMC: " + e.Message;
            Status = cancelled ? "CBMC cancelled" : "Could not run CBMC — see the Equilibrate panel";
        }
        finally
        {
            EqRunning = false;
        }
    }

    // ---- the internal-distance target: production blocks also stop only when ⟨R²(n)⟩/(n⟨b²⟩) is within the tolerance of it
    public static readonly string[] EqTargetChoices = ["None: between blocks only", "RIS polyethylene (alkane cells)", "A curve from a file (n, value)"];
    private int _eqTarget;
    private (double X, double Y)[] _eqTargetFile = [];
    private string _eqTargetFileName = "";
    public int EqTarget
    {
        get => _eqTarget;
        set { if (Set(ref _eqTarget, Math.Clamp(value, 0, 2))) { Raise(nameof(EqTargetNote)); Raise(nameof(EqTargetIsFile)); Raise(nameof(ChainReference)); } }
    }
    public bool EqTargetIsFile => _eqTarget == 2;
    /// <summary>The protocol's final temperature (the RIS reference is taken there).</summary>
    private double EqEndTemperature => _eqProtocol == 1 ? _eqTLow : _eqTFinal;
    public string EqTargetNote => _eqTarget switch
    {
        1 => _risCurve.Length > 0 ? string.Format(CultureInfo.InvariantCulture, "Flory's three-state model at {0:F0} K (the protocol's end); every n within the tolerance", EqEndTemperature)
                                  : "Only for cells of alkanes (every molecule CₙH₂ₙ₊₂): this structure has none, so no target is used",
        2 => _eqTargetFile.Length > 0 ? $"{_eqTargetFileName}: {_eqTargetFile.Length} points, n {_eqTargetFile[0].X:0}–{_eqTargetFile[^1].X:0}" : "Load a file: two columns, n and ⟨R²(n)⟩/(n⟨b²⟩) (# comments)",
        _ => "",
    };
    /// <summary>The curve the chain plot compares with: the loaded target, else the RIS reference for alkanes.</summary>
    public (double X, double Y)[] ChainReference => _eqTarget == 2 && _eqTargetFile.Length > 0 ? _eqTargetFile : _risCurve;

    /// <summary>Reads a target curve: lines of n and the value (whitespace or comma separated); others are skipped.</summary>
    public void LoadEqTarget(string path)
    {
        var inv = CultureInfo.InvariantCulture;
        var pts = new SortedDictionary<int, double>();
        foreach (var raw in File.ReadLines(path))
        {
            var line = raw.Split('#')[0].Trim();
            var f = line.Split([' ', '\t', ',', ';'], StringSplitOptions.RemoveEmptyEntries);
            if (f.Length < 2 || !double.TryParse(f[0], NumberStyles.Float, inv, out var n) || !double.TryParse(f[1], NumberStyles.Float, inv, out var v)) continue;
            if (n < 1 || n > 1_000_000 || Math.Abs(n - Math.Round(n)) > 1e-9 || !(v > 0) || !double.IsFinite(v)) continue;
            pts[(int)Math.Round(n)] = v;
        }
        if (pts.Count == 0) { Status = "No (n, value) rows in " + Path.GetFileName(path); return; }
        _eqTargetFile = pts.Select(kv => ((double)kv.Key, kv.Value)).ToArray();
        _eqTargetFileName = Path.GetFileName(path);
        EqTarget = 2;
        Raise(nameof(EqTargetNote)); Raise(nameof(ChainReference));
        Status = $"Target curve: {_eqTargetFile.Length} points from {_eqTargetFileName}";
    }

    /// <summary>The target indexed by n (0 where there is none), or empty.</summary>
    private double[] EqTargetCurve()
    {
        if (_eqTarget == 2 && _eqTargetFile.Length > 0)
        {
            var t = new double[(int)_eqTargetFile[^1].X + 1];
            foreach (var (n, v) in _eqTargetFile) t[(int)n] = v;
            return t;
        }
        if (_eqTarget == 1 && _risCurve.Length > 0)
        {
            var nmax = (int)_risCurve.Max(p => p.X);
            var c = new double[nmax];
            if (Native.RisCn(EqEndTemperature, nmax, c) != nmax) return [];
            var t = new double[nmax + 1];
            Array.Copy(c, 0, t, 1, nmax);
            return t;
        }
        return [];
    }

    // ---- convergence (Convergence board): criteria and block means
    public ObservableCollection<CriterionRow> EqCriteria { get; } = new();
    private string _eqCriteriaText = "", _eqPhase = "set up";
    public string EqCriteriaText { get => _eqCriteriaText; private set => Set(ref _eqCriteriaText, value); }
    /// <summary>Mean Rg per production block (Å).</summary>
    public (double X, double Y)[] EqRgBlocks { get; private set; } = [];
    /// <summary>Block means of the density and the potential energy at each production block's middle (ps), over the time series.</summary>
    public (double X, double Y)[] EqDensityBlocks { get; private set; } = [];
    public (double X, double Y)[] EqEnergyBlocks { get; private set; } = [];
    public event Action? EqChecksChanged;

    private void LoadEqChecks(string json)
    {
        EqCriteria.Clear();
        EqRgBlocks = [];
        EqDensityBlocks = [];
        EqEnergyBlocks = [];
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
                // the blocks are the last ones of the run: their middles, counted back from its end
                var means = c.GetProperty("blocks").EnumerateArray().Select(v => v.GetDouble()).ToArray();
                var bps = root.TryGetProperty("block_ps", out var bp) ? bp.GetDouble() : _eqBlock;
                var tEnd = _thermo.Count > 0 ? _thermo[^1].TimePs : means.Length * bps;
                var mids = means.Select((v, k) => (tEnd - (means.Length - k - 0.5) * bps, v)).ToArray();
                if (q.StartsWith("density", StringComparison.Ordinal)) EqDensityBlocks = mids;
                else if (energy)   // per atom in the checks, the whole cell on the plot
                {
                    var natoms = _doc?.Summary().Atoms ?? 0;
                    EqEnergyBlocks = natoms > 0 ? mids.Select(m => (m.Item1, m.v * natoms)).ToArray() : [];
                }
            }
            var met = EqCriteria.Count(r => r.State == "pass");
            EqCriteriaText = $"{met} of {EqCriteria.Count} criteria met";
        }
        else EqCriteriaText = "";
        Raise(nameof(EqCanDecide));
        EqChecksChanged?.Invoke();
    }

    // ---------------------------------------------------------------- Pack
    private double _packX = 40, _packY = 40, _packZ = 40, _packTol = 2.0;
    private int _packCount = 100, _packSeed = 1;
    private bool _packPeriodic = true, _packing;
    private string _packText = "", _packBaseDir = Environment.CurrentDirectory;
    private string _packLog = "Add molecules (or open an input: CAPS or packmol syntax). Every molecule is a rigid body; overlaps closer than the distance are " +
                              "removed by minimisation, and a cell that misses the tolerance is never produced.";
    private CancellationTokenSource? _packCancel;

    public decimal? PackXD { get => (decimal)_packX; set { _packX = Math.Clamp((double)(value ?? 40m), 5, 10000); Raise(); } }
    public decimal? PackYD { get => (decimal)_packY; set { _packY = Math.Clamp((double)(value ?? 40m), 5, 10000); Raise(); } }
    public decimal? PackZD { get => (decimal)_packZ; set { _packZ = Math.Clamp((double)(value ?? 40m), 5, 10000); Raise(); } }
    public decimal? PackTolD { get => (decimal)_packTol; set { _packTol = Math.Clamp((double)(value ?? 2m), 0.5, 10); Raise(); } }
    public decimal? PackCountD { get => _packCount; set { _packCount = Math.Clamp((int)(value ?? 100), 1, 10_000_000); Raise(); } }
    public decimal? PackSeedD { get => _packSeed; set { if (value == null) return; _packSeed = Math.Max(1, (int)value); if (_packSeedC != null && _packSeedC.Value != _packSeed) _packSeedC.Value = _packSeed; Raise(); } }
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

    /// <summary>The CAPS form of Pack input (cell, distance, molecule … end); packmol's is read too.</summary>
    public static bool IsCapsPack(string text)
    {
        foreach (var raw in text.Split('\n'))
        {
            var w = raw.Split('#')[0].Trim().Split((char[]?)null, StringSplitOptions.RemoveEmptyEntries);
            if (w.Length == 0) continue;
            var k = w[0].ToLowerInvariant();
            if (k is "cell" or "distance" or "molecule" or "compress" or "save" or "loops" or "iterations") return true;
            if (k is "tolerance" or "structure" or "pbc" or "output" or "filetype" or "nloop" or "maxit") return false;
        }
        return true;
    }

    private void ParsePackText()
    {
        PackItems.Clear();
        var inv = CultureInfo.InvariantCulture;
        string? file = null;
        var number = "1";
        var fixedMol = false;
        var ffChoice = 0;
        var constraints = new List<string>();
        PackCellText = "no periodic cell";
        foreach (var raw in _packText.Split('\n'))
        {
            var line = raw.Split('#')[0].Trim();
            if (line.Length == 0) continue;
            var w = line.Split((char[]?)null, StringSplitOptions.RemoveEmptyEntries);
            var key = w[0].ToLowerInvariant();
            if (file == null)
            {
                if (key is "tolerance" or "distance" && w.Length > 1) PackTolText = w[1] + " Å";
                else if (key == "pbc" && w.Length >= 7)
                    PackCellText = string.Format(inv, "{0} × {1} × {2} Å, periodic", w[4], w[5], w[6]);
                else if (key == "cell" && w.Length == 4) PackCellText = string.Format(inv, "{0} × {1} × {2} Å, periodic", w[1], w[2], w[3]);
                else if (key == "cell" && w.Length == 8) PackCellText = string.Format(inv, "{0} × {1} × {2} Å, periodic", w[5], w[6], w[7]);
                else if (key is "structure" or "molecule" && w.Length > 1) { file = line[(line.IndexOf(' ') + 1)..].Trim(); number = "1"; fixedMol = false; ffChoice = 0; constraints.Clear(); }
            }
            else if (key == "end" && !(w.Length > 1 && w[1].Equals("atoms", StringComparison.OrdinalIgnoreCase)))
            {
                var name = Path.GetFileNameWithoutExtension(file);
                var colour = PackColours[PackItems.Count % PackColours.Length];
                PackItems.Add(new PackItem(name, constraints.Count > 0 ? string.Join(" · ", constraints) : "anywhere in the cell",
                    fixedMol ? "fixed" : $"× {number}", colour, file, PackItems.Count, ffChoice));
                file = null;
            }
            else if (key is "number" or "count" && w.Length > 1) number = w[1];
            else if (key == "fixed") fixedMol = true;
            else if (key is "forcefield" or "water" && w.Length > 1) ffChoice = PackFfChoiceOf(key, line[(line.IndexOf(' ') + 1)..].Trim());
            else if (key is "inside" or "outside" or "over" or "below" or "above" or "in" or "not" or "atoms") constraints.Add(line);
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
        PackText = string.Format(CultureInfo.InvariantCulture, "# CAPS Pack\ndistance  {0:0.###}         # Å: no two molecules closer\nseed      {1}\n", _packTol, _packSeed) +
                   (_packPeriodic ? string.Format(CultureInfo.InvariantCulture, "cell      {0:0.###} {1:0.###} {2:0.###}   # Å, periodic\n", _packX, _packY, _packZ) : "") + "\n";
    }

    /// <summary>Append a structure block for a molecule file, placed inside the whole cell.</summary>
    public void AddPackStructure(string path)
    {
        if (_packText.Trim().Length == 0) NewPackInput();
        var inv = CultureInfo.InvariantCulture;
        PackText = _packText.TrimEnd() + (IsCapsPack(_packText)
            ? string.Format(inv, "\n\nmolecule  {0}\n  count   {1}\n  in      box from 0 0 0 to {2:0.###} {3:0.###} {4:0.###}\nend\n", path, _packCount, _packX, _packY, _packZ)
            : $"\n\nstructure {path}\n  number {_packCount}\n  inside box {BoxText()}\nend structure\n");
    }

    public void AddPackExample(string samples)
    {
        NewPackInput();
        AddPackStructure(Path.Combine(samples, "water.pdb"));
    }

    public void LoadPackInput(string path)
    {
        var text = File.ReadAllText(path);
        // a packmol file shown in the CAPS form (it runs the same; Save writes it back in either form)
        if (!IsCapsPack(text)) try { text = CapsDocument.PackConvert(text, true); } catch { }
        PackText = text;
        PackBaseDir = Path.GetDirectoryName(path) ?? ".";
    }

    public async Task RunPack() => await RunPack(null, null);

    /// <summary>Packs the page's input (or a captured one: a queued run keeps the input as it was when queued).</summary>
    private async Task RunPack(string? presetText, string? presetBase)
    {
        if (!Idle || (presetText ?? _packText).Trim().Length == 0) return;
        string text;
        try
        {
            text = presetText ?? PackTextToRun();
            // new each time: this run's seed in the input's seed line (a fixed seed is the one the input says)
            if (presetText == null && PackSeedChoice.Fresh) text = WithSeedLine(text, _packSeed = PackSeedChoice.Take());
        }
        catch (Exception e) { PackLog = "Could not pack.\n" + e.Message; Status = "Could not pack — see the Pack panel"; return; }
        // packing again from the last packed cell: the new cell replaces it (the same host, or the cell itself was open)
        var replaces = presetText == null && _packResult != null && (ReferenceEquals(_doc, _packResult) || (_packRunHost != null && ReferenceEquals(_packRunHost, _packHost)))
            ? _packResult : null;
        var runHost = _packRunHost;
        var runHostTitle = runHost == null ? "" : ReferenceEquals(runHost, _doc) ? Title.Replace(" (unsaved)", "") : _packHostTitle;
        // the macro line, from the input as it is before packing (the packed cell becomes the open structure)
        var packScript = Recording && presetText == null ? PackPython() : null;
        PackDone = false;
        Packing = true;
        var packed = false;
        _packCancel = new CancellationTokenSource();
        var token = _packCancel.Token;
        var baseDir = presetBase ?? _packBaseDir;
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
            _packResult = doc;
            _packHost = runHost;
            _packHostTitle = runHostTitle;
            var baseName = runHost != null && runHostTitle.Length > 0 ? runHostTitle + " + packed" : $"packed_{s.Molecules}_molecules";
            Show(doc, $"{baseName} (unsaved)");
            if (_activeItem != null) _activeItem.BuildSettings = PackSettingsJson(runHost != null ? runHostTitle : "");   // Edit brings these back
            if (replaces != null && !ReferenceEquals(replaces, doc) && ProjectItems.FirstOrDefault(i => ReferenceEquals(i.Doc, replaces)) is { } oldItem)
            {
                ProjectItems.Remove(oldItem);   // the earlier packing: replaced, not kept beside the new one
                if (!replaces.LongRunning) replaces.Dispose();
                RaiseProject();
            }
            GrownUnsaved = true;
            Packing = false;
            PackLog = report;
            PackBad = 0;
            var m = System.Text.RegularExpressions.Regex.Match(report, @"smallest distance between molecules ([0-9.]+) Å");
            if (m.Success) PackDmin = m.Groups[1].Value + " Å";
            Status = $"Packed {s.Molecules:N0} molecules ({s.Atoms:N0} atoms){(replaces != null ? " · the earlier packed cell replaced" : "")} · export it to LAMMPS or GROMACS, or minimise / run dynamics first";
            Raise(nameof(PackStartNote));
            MarkPipeline("Pack");
            if (packScript != null) RecordScript(packScript);
            packed = true;
            PackDone = true;
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
        if (packed) await AssignForBuilder(pack: true);
    }

    public void CancelPack() => _packCancel?.Cancel();

    // ---------------------------------------------------------------- React
    // short names in the list; what each does in RxSetNotes, shown under the choice
    public static readonly string[] ReactionSets = ["C–C crosslink (CH₂–CH₂)", "Epoxy–amine cure", "Sulfur vulcanization", "Peroxide cure (allylic)",
        "Silane coupling", "ENR + carboxylic acid", "ENR + MAH", "Esterification", "C–C crosslink (any C–H)", "Custom"];
    public static readonly string[] RxSetNotes =
    [
        "Two CH₂ sites become CH–CH; each CH₂ reacts once, two H leave (H₂ in radiation crosslinking, the peroxide's alcohols in a peroxide cure — the same network). PE, EPDM, the CH₂ of polydienes.",
        "Epoxide ring opened by a primary, then a secondary amine (each N–H once): the cure of epoxy resins.",
        "Accelerated sulfur cure of a diene rubber: H–Sx–H donors (inserted below, dosed in phr) bond to allylic C–H carbons; two in turn make C–Sx–C bridges.",
        "Peroxide cure of a diene rubber: allylic carbons of two chains joined C–C.",
        "Silane coupling agent (TESPT / TESPD) polysulfide to the allylic C–H of a diene rubber → C–S.",
        "ENR epoxide opened by a carboxylic acid at its tertiary carbon → β-hydroxy ester, nothing leaves; a diacid (maleic acid) bridges two epoxides of two chains — one molecule per link.",
        "MAH is opened by an OH (to a half-ester acid, which then opens an ENR epoxide): it needs OH groups — ring-opened ENR, an alcohol, PBS ends or moisture; a fresh ENR cell has none. For one molecule bridging two epoxides choose ENR + carboxylic acid with maleic acid (MAH once hydrolysed).",
        "COOH + OH → ester + H₂O (polyesters, PBS chain extension).",
        "C–C crosslink of any sp³ C–H carbons (CH₃, CH₂ or CH — polypropylene's tertiary CH); two H leave.",
        "The reaction text as you write it (template syntax, or from the library and the editor).",
    ];
    public string RxSetNote => _rxSet >= 0 && _rxSet < RxSetNotes.Length ? RxSetNotes[_rxSet] : "";
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

    public int RxSet { get => _rxSet; set { if (Set(ref _rxSet, value)) { LoadReactionSet(); Raise(nameof(RxShowInsert)); Raise(nameof(RxSetNote)); RxCrosslinkerDefault(); SyncRxChoice(); } } }
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
    public string RxText { get => _rxText; set { if (Set(ref _rxText, value)) { RefreshRxReactions(); QueueRxSites(); } } }
    public string RxLog { get => _rxLog; private set => Set(ref _rxLog, value); }
    public bool RxRelax { get => _rxRelax; set { if (Set(ref _rxRelax, value)) Raise(nameof(RxMdEnabled)); } }
    public bool Reacting { get => _reacting; private set { if (Set(ref _reacting, value)) RaiseBusy(); } }
    public bool CanReact => _doc != null && Idle;
    public decimal? RxCyclesD { get => _rxCycles; set { _rxCycles = Math.Clamp((int)(value ?? 50), 1, 100000); Raise(); } }
    public decimal? RxPerCycleD { get => _rxPerCycle; set { _rxPerCycle = Math.Clamp((int)(value ?? 5), 1, 100000); Raise(); } }
    public decimal? RxSeedD { get => _rxSeed; set { if (value == null) return; _rxSeed = Math.Max(1, (int)value); if (_rxSeedC != null && _rxSeedC.Value != _rxSeed) _rxSeedC.Value = _rxSeed; Raise(); } }
    public decimal? RxTargetD { get => (decimal)_rxTarget; set { _rxTarget = Math.Clamp((double)(value ?? 1m), 0.001, 1); Raise(); } }
    public decimal? RxCaptureD { get => (decimal)_rxCapture; set { _rxCapture = Math.Clamp((double)(value ?? 0m), 0, 20); Raise(); } }
    public decimal? RxMdPsD { get => (decimal)_rxMdPs; set { _rxMdPs = Math.Clamp((double)(value ?? 0m), 0, 10000); Raise(); } }
    // protocol: Polymatic cycles (react → relax → optional MD) or REACTER-style checks during one continuous NVT run
    private bool _rxDuringMd;
    public bool RxDuringMd
    {
        get => _rxDuringMd;
        set
        {
            if (!Set(ref _rxDuringMd, value)) return;
            if (value && _rxMdPs > 1) RxMdPsD = 0.1m;   // checks every 0.1 ps, as bond/react's Nevery
            Raise(nameof(RxPolymatic)); Raise(nameof(RxMdEnabled)); Raise(nameof(RxMdLabel));
        }
    }
    public bool RxPolymatic { get => !_rxDuringMd; set { if (value) RxDuringMd = false; } }
    public bool RxMdEnabled => _rxRelax || _rxDuringMd;
    public string RxMdLabel => _rxDuringMd ? "Check every (ps)" : "MD per cycle (ps)";
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
                4 => CapsDocument.ReactionTemplate("polysulfide_allylic"),
                5 => CapsDocument.ReactionTemplate("enr_acid_ester"),
                6 => CapsDocument.ReactionTemplate("anhydride_alcohol") + "\n" + CapsDocument.ReactionTemplate("enr_acid_ester"),
                7 => CapsDocument.ReactionTemplate("ester_condensation"),
                8 => CapsDocument.ReactionTemplate("cc_crosslink_any"),
                _ => _rxText,
            };
        }
        catch (Exception e) { RxLog = e.Message; }
    }

    /// <summary>The React settings as the core takes them (also captured when a run is queued).</summary>
    private CapsReactOpts ReactOptions() => new()
    {
        Seed = (ulong)(_rxSeed = RxSeedChoice.Take()), MaxCycles = _rxCycles, MaxPerCycle = _rxPerCycle, TargetConversion = _rxTarget, Capture = _rxCapture,
        Relax = _rxRelax ? 1 : 0, RelaxIterations = _rxRelaxIt, MdPs = _rxRelax || _rxDuringMd ? _rxMdPs : 0, Temperature = _rxTemp, Cutoff = _relaxCutoff, Coulomb = _relaxCoulomb ? 1 : 0,
        DuringMd = _rxDuringMd ? 1 : 0,
        FieldMode = _rxUseField ? 0 : 1, BetweenChains = _rxBetween ? 1 : 0, KeepByproducts = _rxKeepBy ? 1 : 0,
        Selection = RxSeveral && _rxByWeights ? 1 : 0,
        Weights = RxSeveral && _rxByWeights ? string.Join(",", RxReactions.Select(r => r.WeightD.ToString(CultureInfo.InvariantCulture))) : "",
        AutoCapture = _rxAutoCapture ? 1 : 0, CaptureMax = _rxCaptureMax, CaptureStep = 0.5,
        TargetKind = _rxTargetKind, TargetValue = _rxTargetKind == 0 ? 0 : _rxTargetValue,
        SitesPerChain = (int)_rxSitesPer,
    };

    public async Task RunReact() => await RunReact(null, null);

    private async Task RunReact(CapsReactOpts? preset, string? presetText)
    {
        if (_doc == null || !Idle || (presetText ?? _rxText).Trim().Length == 0) return;
        // the crosslinker the asked-for degree still lacks goes in first (a higher degree, more molecules), then the force field
        if (preset == null && RxTopUp && RxTargetKind == 5 && RxCanInsertX) await InsertRxCrosslinker();
        if (preset == null && !await RxAssignField()) return;   // the reaction's force field, complete for the start
        var doc = _doc;
        var rxSnap = preset == null ? PageSnapshot("react") : null;   // the reaction as set on the page, kept with the network
        Reacting = true;
        IsPlaying = false;
        _rxCancel = new CancellationTokenSource();
        var token = _rxCancel.Token;
        var o = preset ?? ReactOptions();
        _rxRows.Clear();
        StartLive();
        RxLiveText = "";
        RxNetworkText = "";
        ReactChanged?.Invoke();
        RxLog = "Finding reactive pairs…";
        Status = $"Reacting {Title}…";
        var text = presetText ?? _rxText;
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
                    OnReactCurves();
                    if (!finished)
                        RxLog = string.Format(inv, "cycle {0}: {1} reactions ({2} in all) · conversion {3:F3} · {4} clusters, largest {5:F1} % · {6:F0} s",
                            row.Cycle, row.Reactions, row.Total, row.Conversion, row.Clusters, 100 * row.LargestFraction, sw.Elapsed.TotalSeconds);
                    ReactChanged?.Invoke();
                });
                return !token.IsCancellationRequested;
            }, ReactLiveReceiver()));
            finished = true;
            var failedAt = System.Text.RegularExpressions.Regex.Match(report, @"failed at cycle (\d+)");
            RxLog = report + "\n" + FloryText;
            AfterRun(doc, " · reacted");
            if (rxSnap != null && _activeItem != null && ReferenceEquals(_activeItem.Doc, doc)) _activeItem.ReactSettings = rxSnap;
            LoadReactSummary(doc);
            Raise(nameof(RxFieldText));
            Raise(nameof(RxBrFfChoices));
            Status = failedAt.Success
                ? $"React failed at cycle {failedAt.Groups[1].Value}; the structure after cycle {int.Parse(failedAt.Groups[1].Value, inv) - 1} is kept — React again to continue from it"
                : "Reaction run finished · save the network (LAMMPS data carries the force field when every atom is typed)";
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
            QueueRxSites();
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
        // the macro recorder: the run as the caps package would do it (relax records itself)
        if (suffix == " · MD") RecordScript(MdPython());
        else if (suffix == " · reacted") RecordScript(ReactPython());
        else if (suffix == " · equilibrated") RecordScript(EqPython());
        RefreshAppColumns();   // a run's record has its own columns (or none)
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
        Raise(nameof(CanAddRestraint));
        RefreshSelection();
        RefreshSummary();
        RefreshRdf();
        RefreshLegend();
        RefreshMolecules();
        Notes.Clear();
        foreach (var n in doc.Notes()) Notes.Add(n);
        LoadFileChecks();
        FieldInfoText = "";
        MarkPipeline(suffix.Contains("relaxed") ? "Relax" : suffix.Contains("equilibrated") ? "Equilibrate" : suffix.Contains("MD") ? "Dynamics" : "Run");
        RenderRequested?.Invoke();
    }

    public void SaveDocument(string path)
    {
        if (_doc == null) return;
        _doc.Save(path);
        Record($"doc.save({PyStr(path)})");
        GrownUnsaved = false;
        Title = System.IO.Path.GetFileName(path);
        Status = $"Saved {path}";
        Remember(path, null);
    }

    private void Show(CapsDocument doc, string title)
    {
        if (Busy) { doc.Dispose(); Status = "Wait for the run to finish (or cancel it) before opening another structure"; return; }
        if (_wrap) doc.SetWrap(true);
        // a new structure joins the project; the one that was active stays there (pick it again from the tabs or the list)
        var item = AddProjectItem(doc, title);
        Document = doc;
        RestraintsFollow(doc);
        ClearFocus();
        if (IsVisualize) Avalonia.Threading.Dispatcher.UIThread.Post(ApplyPipeline);
        Field.Reset();
        try { Field.LoadReport(doc); } catch { Field.Reset(); }   // a structure built with its own model (Kremer–Grest) comes assigned
        _pipeAutoFf = false;
        Analyze.Load("");
        SyncHeld(); SyncFixed();
        Title = title;
        PipelineNewDocument(title);
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
        Raise(nameof(CanAddRestraint));
        RefreshSelection();
        RefreshSummary();
        RefreshRdf();
        RefreshLegend();
        RefreshMolecules();
        IsPlaying = false;
        Notes.Clear();
        foreach (var n in doc.Notes()) Notes.Add(n);
        if (IsProvenance) LoadProvenance();
        AutoLod(doc);
        AutoStyle(doc);
        LoadFileChecks();
        UpdateItemInfo(item);
        var look = FileChecks.Count(c => c.NeedsLook);
        Status = $"Opened {Title} · {s.Atoms.ToString("N0", CultureInfo.InvariantCulture)} atoms · {s.Format}" + (look > 0 ? $" · {look} file check{(look == 1 ? "" : "s")} need a look" : "");
        RenderRequested?.Invoke();
    }

    private void RefreshSummary()
    {
        RequestProperties();
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
        Raise(nameof(CanAddRestraint));
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
        QueueSelBar();
        RequestProperties();
        RefreshMolInspector();
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
        if (CoordinationOf(index) is { Length: > 0 } coord) PickedRows.Add(new("Coordination", coord));
        foreach (var (i, d) in _doc.Neighbours(index, 4))
        {
            var b = _doc.Atom(i);
            NeighbourRows.Add(new($"{b.Id} · {b.ElementSymbol} · mol {b.Mol}", string.Format(inv, "{0:F3} Å", d)));
        }
        if (_selection.Count >= 2)
        {
            var v = _doc.Measure(_selection.ToArray());
            var ids = string.Join("–", _selection.Select(i => _doc.Atom(i).Id));
            Raise(nameof(CanAddRestraint));
            MeasureText = _selection.Count switch
            {
                2 => string.Format(inv, "Distance {0}: {1:F3} Å", ids, v),
                3 => string.Format(inv, "Angle {0}: {1:F2}°", ids, v),
                _ => string.Format(inv, "Dihedral {0}: {1:F2}°", ids, v),
            };
            if (_selection.Count == 4 && OutOfPlane(_selection.ToArray()) is { Length: > 0 } oop) MeasureText += "\n" + oop;
            MeasureTarget = _selection.Count == 2 ? v.ToString("F3", inv) : v.ToString("F2", inv);
            Raise(nameof(CanSetMeasure));
            Status = MeasureText;
        }
        else
        {
            MeasureText = "";
            Raise(nameof(CanAddRestraint));
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
        EnsureRdfPairs();
        var (ea, eb) = _rdfElements[Math.Clamp(_rdfPair, 0, _rdfElements.Length - 1)];
        var sampled = s.Atoms > 20000 ? " · 20 000 sampled centres" : "";
        RdfNote = string.Format(CultureInfo.InvariantCulture, "{0} · {1} · r ≤ {2:F0} Å · 0.2 Å bins · frame {3}{4}",
            RdfPairs[Math.Clamp(_rdfPair, 0, RdfPairs.Count - 1)], _rdfInter ? "between molecules" : "all pairs", rmax, _frame, sampled);
        if (s.Atoms < 200_000) { RdfCurve = _doc.Rdf(ea, eb, rmax, 0.2, _rdfInter); return; }
        // large cells: off the UI thread
        var doc = _doc;
        var inter = _rdfInter;
        Task.Run(() => { try { return doc.Rdf(ea, eb, rmax, 0.2, inter); } catch { return []; } })
            .ContinueWith(t => Avalonia.Threading.Dispatcher.UIThread.Post(() => { if (_doc == doc) RdfCurve = t.Result; }));
    }

    public CapsRenderOpts ViewOptions(int w, int h, int supersample)
    {
        var o = new CapsRenderOpts
        {
        Width = w, Height = h, Supersample = supersample,
        Background = _module == 19 ? _renderBg : _viewBackground,   // Render: the preview on the image's own background
        ColourBy = _colour, Style = _style,
        Outlines = OutlineLevel, ShowCell = _showCell ? 1 : 0,
        Highlight0 = _selection.Count > 0 ? _selection[0] : -1,
        Highlight1 = _selection.Count > 1 ? _selection[1] : -1,
        Highlight2 = _selection.Count > 2 ? _selection[2] : -1,
        Highlight3 = _selection.Count > 3 ? _selection[3] : -1,
        Focus = _focusAtom >= 0 ? _focusAtom + 1 : 0,
        AmbientOcclusion = _module == 19 ? (_renderAo ? 1 : 0) : (_viewAo ? 1 : 0),
        DepthCue = _module == 19 ? (_renderDepth ? 1 : 0) : (_depthCue ? 1 : 0),
        };
        ApplyLod(ref o);
        return o;
    }

    /// <summary>The Field page's view: coloured by force-field type, ball and stick, the selected row's atom highlighted.</summary>
    public CapsRenderOpts FieldViewOptions(int w, int h, int supersample) => new()
    {
        Width = w, Height = h, Supersample = supersample,
        Background = _viewBackground,
        ColourBy = 2, Style = _style == 3 ? 0 : _style,
        Outlines = OutlineLevel, DepthCue = _depthCue ? 1 : 0, ShowCell = 0,
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
        o.LodNear = o.LodFar = 0;   // exports in full detail
        o.Background = _exportBackground;
        o.Highlight0 = o.Highlight1 = o.Highlight2 = o.Highlight3 = -1;
        o.Focus = 0;
        return o;
    }

    public void ResetView() => FlyTo(new CapsCamera { Yaw = 0.55, Pitch = 0.40, Zoom = 1, Perspective = Camera.Perspective });

    public void SetView(double yaw, double pitch) => FlyTo(new CapsCamera { Yaw = yaw, Pitch = pitch, Zoom = Camera.Zoom, Perspective = Camera.Perspective });
}
