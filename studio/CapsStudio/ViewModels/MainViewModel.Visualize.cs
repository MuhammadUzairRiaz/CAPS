using System.Collections.ObjectModel;
using System.ComponentModel;
using System.Globalization;
using System.Text.Json.Nodes;
using Avalonia.Media;
using CapsStudio.Interop;

namespace CapsStudio.ViewModels;

/// <summary>One field of a step's editor: a number, text, expression, switch or choice, stored in the step's parameters.</summary>
public sealed class StepField : INotifyPropertyChanged
{
    public event PropertyChangedEventHandler? PropertyChanged;
    private void Raise(string n) => PropertyChanged?.Invoke(this, new PropertyChangedEventArgs(n));
    public string Key { get; init; } = "";
    public string Label { get; init; } = "";
    public string Kind { get; init; } = "text";         // text | expression | number | bool | choice | vector
    public string Hint { get; init; } = "";
    public string[] Choices { get; init; } = [];
    public Action? Changed;
    private string _text = "";
    private bool _on;
    public string Text { get => _text; set { if (_text == value) return; _text = value; Raise(nameof(Text)); Changed?.Invoke(); } }
    public bool On { get => _on; set { if (_on == value) return; _on = value; Raise(nameof(On)); Changed?.Invoke(); } }
    public int ChoiceIndex
    {
        get => Math.Max(0, Array.IndexOf(Choices, _text));
        set { if (value >= 0 && value < Choices.Length) Text = Choices[value]; Raise(nameof(ChoiceIndex)); }
    }
    public bool IsText => Kind is "text" or "expression" or "number" or "vector" or "file" or "matrix";
    public bool IsFile => Kind == "file";
    public bool IsBool => Kind == "bool";
    public bool IsChoice => Kind == "choice";
    public bool IsNote => Kind == "note";
    public bool IsCode => Kind == "code";                 // several lines, sent with Run (not per keystroke)
    public bool IsAction => Kind == "action";             // a button (Label) that fills the step from the view
    private string _draft = "";
    public string Draft { get => _draft; set { if (_draft == value) return; _draft = value; Raise(nameof(Draft)); Raise(nameof(DraftChanged)); } }
    public bool DraftChanged => _draft != _text;
    public void Commit() { Text = _draft; Raise(nameof(DraftChanged)); }
    private string _output = "";
    /// <summary>A code field's console: what the last run printed.</summary>
    public string Output { get => _output; set { if (_output == value) return; _output = value; Raise(nameof(Output)); Raise(nameof(HasOutput)); } }
    public bool HasOutput => _output.Length > 0;                 // a line of explanation under the fields (Hint)
    public bool ShowLabel => Kind is not ("bool" or "note" or "action");
    public bool IsMono => Kind is "expression" or "number" or "vector" or "file" or "matrix";
}

/// <summary>A step in the pipeline list.</summary>
public sealed class PipelineRow : INotifyPropertyChanged
{
    public event PropertyChangedEventHandler? PropertyChanged;
    private void Raise(string n) => PropertyChanged?.Invoke(this, new PropertyChangedEventArgs(n));
    public string Type { get; init; } = "";
    public string Title { get; init; } = "";
    public string Icon { get; init; } = "sliders";
    public JsonObject Params { get; set; } = new();
    public Action? Toggled;
    private bool _enabled = true, _selected;
    private string _summary = "", _level = "ok";
    public bool Enabled { get => _enabled; set { if (_enabled == value) return; _enabled = value; Raise(nameof(Enabled)); Toggled?.Invoke(); } }
    public bool Selected { get => _selected; set { _selected = value; Raise(nameof(Selected)); } }
    public string Summary { get => _summary; set { _summary = value; Raise(nameof(Summary)); } }
    /// <summary>What the step printed (a Python step's console).</summary>
    public string Output { get; set; } = "";
    public string Level { get => _level; set { _level = value; Raise(nameof(Level)); Raise(nameof(SummaryBrush)); } }
    public IBrush SummaryBrush => Tokens.Brush(_level switch { "error" => "ErrB", "warning" => "WarnB", _ => "DimB" });

    // groups (design/boards/PipelineGroups): a "group" parameter the core ignores; consecutive steps of one group
    // show under one header that turns them all on or off and folds them away
    public string Group => Params["group"] is JsonValue v && v.TryGetValue<string>(out var g) ? g : "";
    public bool HasGroup => Group.Length > 0;
    private bool _header, _visible = true, _groupOn = true, _collapsed;
    private int _groupCount;
    public bool GroupHeader { get => _header; set { _header = value; Raise(nameof(GroupHeader)); } }
    public int GroupCount { get => _groupCount; set { _groupCount = value; Raise(nameof(GroupCount)); Raise(nameof(GroupLabel)); } }
    public string GroupLabel => $"{GroupCount} step{(GroupCount == 1 ? "" : "s")}";
    public bool RowVisible { get => _visible; set { _visible = value; Raise(nameof(RowVisible)); } }
    public bool Collapsed { get => _collapsed; set { _collapsed = value; Raise(nameof(Collapsed)); Raise(nameof(GroupChevron)); } }
    public string GroupChevron => _collapsed ? "chevr" : "chev";
    public Action<PipelineRow, bool>? GroupToggled;
    public Action<PipelineRow>? GroupFolded;
    public bool GroupOn { get => _groupOn; set { if (_groupOn == value) return; _groupOn = value; Raise(nameof(GroupOn)); GroupToggled?.Invoke(this, value); } }
    public void SetGroupOn(bool on) { _groupOn = on; Raise(nameof(GroupOn)); }
    public void RaiseGroup() { Raise(nameof(Group)); Raise(nameof(HasGroup)); }
}

public sealed record StepKind(string Type, string Title, string About, string Group, string Icon);
public sealed record TableRow(string[] Cells, bool Selected);
public sealed record LegendSwatch(string Label, IBrush Brush);

/// <summary>Analyze › Visualize (design/boards/VisPipeline, PipelineSteps, DataInspector): a pipeline of non-destructive
/// steps run on every shown frame, the data inspector (particles, bonds, attributes, tables) and each step's editor.</summary>
public sealed partial class MainViewModel
{
    public bool IsVisualize => _module == 20;

    public ObservableCollection<PipelineRow> PipelineRows { get; } = new();
    public ObservableCollection<StepField> StepFields { get; } = new();
    public ObservableCollection<Row> PipeAttributes { get; } = new();
    public ObservableCollection<TableRow> InspectorRows { get; } = new();
    public ObservableCollection<string> InspectorColumns { get; } = new();
    public ObservableCollection<LegendSwatch> PipeLegendEntries { get; } = new();
    public ObservableCollection<string> PipeTables { get; } = new();

    /// <summary>The Add step library, grouped as on the PipelineSteps board (only steps CAPS runs).</summary>
    public static readonly StepKind[] StepLibrary =
    [
        new("colour_coding", "Colour coding", "any property, categorical or continuous", "Colour & style", "eye"),
        new("assign_colour", "Assign colour", "to the selection", "Colour & style", "eye"),
        new("transparency", "Transparency", "the selection or by property: see through fillers", "Colour & style", "eye"),
        new("particle_radius", "Particle radius", "the selection or by property", "Colour & style", "atom"),
        new("select_expression", "Expression selection", "Type == 2 && Position.Z > 13", "Select", "filter"),
        new("manual_selection", "Manual selection", "atoms picked or lassoed in the view", "Select", "filter"),
        new("expand_selection", "Expand selection", "by bonds or distance", "Select", "filter"),
        new("invert_selection", "Invert selection", "selected ↔ not selected", "Select", "filter"),
        new("clear_selection", "Clear selection", "nothing selected", "Select", "filter"),
        new("slice", "Slice", "slab by normal and width", "Modify", "scissors"),
        new("delete_selected", "Delete selected", "remove the selected particles", "Modify", "scissors"),
        new("replicate", "Replicate", "periodic images", "Modify", "copy"),
        new("wrap", "Wrap into cell", "fold positions into the cell", "Modify", "cube"),
        new("compute_property", "Compute property", "expression per particle", "Modify", "terminal"),
        new("coordination", "Coordination & RDF", "neighbours within a cutoff, g(r)", "Measure", "chart"),
        new("cluster", "Cluster analysis", "by bonds or cutoff, sizes, Rg", "Measure", "layers"),
        new("histogram", "Histogram", "distribution of a property", "Measure", "chart"),
        new("binning", "Spatial binning", "1-D profile along an axis", "Measure", "chart"),
        new("molecule_shape", "Molecule shape", "Rg, κ², asphericity per molecule", "Structure", "atom"),
        new("topology", "Topology distributions", "bond lengths, angles, dihedrals", "Structure", "bond"),
        new("displacements", "Displacements", "vs a reference frame · MSD", "Trajectory", "move"),
        new("msd", "Mean-square displacement", "MSD(τ) of atoms and chain centres, D", "Trajectory", "chart"),
        new("scatter", "Scatter plot", "one property against another · Pearson r", "Measure", "chart"),
        new("smooth", "Smooth trajectory", "positions averaged over frames", "Trajectory", "history"),
        new("unwrap", "Unwrap", "molecules whole across the boundary", "Modify", "cube"),
        new("affine_transform", "Affine transformation", "strain, shear or rotate particles and cell", "Modify", "move"),
        new("freeze_property", "Freeze property", "values at a reference frame, on every frame", "Trajectory", "pin"),
        new("orientation", "Chain orientation", "P₂ per atom, S, director, local crystallinity", "Structure", "grow"),
        new("cna", "Common neighbour analysis", "FCC, HCP, BCC, icosahedral (adaptive CNA)", "Structure", "atom"),
        new("ptm", "Polyhedral template matching", "structure, orientation and strain per particle", "Structure", "atom"),
        new("combine", "Combine datasets", "add the particles of a second file", "Modify", "layers"),
        new("wigner_seitz", "Wigner–Seitz defects", "vacancies and interstitials against a reference lattice", "Measure", "atom"),
        new("centrosymmetry", "Centrosymmetry", "Kelchner's parameter: defects, surfaces", "Structure", "atom"),
        new("create_bonds", "Create bonds", "from distances or a cutoff", "Visual", "link"),
        new("python", "Python step", "your script with an @step function (caps.pipeline API)", "Automate", "terminal"),
        new("primitive_paths", "Primitive paths", "chains pulled tight without crossing · N_e", "Structure", "bond"),
        new("voids", "Voids & pores", "accessible volume for a probe, voids by size", "Structure", "atom"),
        new("voronoi", "Voronoi analysis", "exact cells: volumes, faces, Voronoi index", "Structure", "hex"),
        new("density_field", "Density field", "smoothed mass density, profile, slice", "Structure", "layers"),
        new("vectors", "Vectors", "end-to-end, dipoles, displacements, velocities", "Visual", "move"),
        new("trajectory_lines", "Trajectory lines", "paths of chain centres or particles", "Visual", "history"),
    ];
    public static readonly string[] StepGroups = ["Colour & style", "Select", "Modify", "Structure", "Measure", "Trajectory", "Visual", "Automate"];

    private PipelineRow? _pipeSel;
    private bool _stepLibraryOpen, _pipeLegendVisible = true;
    private int _inspectorTab, _inspectorOffset, _pipeTable;
    private string _inspectorFilter = "", _inspectorNote = "", _pipeError = "", _stepSearch = "";
    private JsonNode? _pipeResult;
    private int _pipeGen;
    private bool _showTableAfterApply;
    private const int InspectorPage = 200;

    public void OpenVisualize()
    {
        if (_doc == null) { SetModule(20); return; }
        if (PipelineRows.Count == 0)
        {
            // the board's starting pipeline: colour chains, lighten hydrogens
            AddStep("colour_coding", select: true);
        }
        SetModule(20);
        ApplyPipeline();
    }

    public PipelineRow? PipeSelected
    {
        get => _pipeSel;
        set
        {
            if (_pipeSel == value) return;
            if (_pipeSel != null) _pipeSel.Selected = false;
            _pipeSel = value;
            if (_pipeSel != null) _pipeSel.Selected = true;
            Raise();
            Raise(nameof(HasPipeSelected));
            Raise(nameof(PipeStepTitle));
            Raise(nameof(PipeStepNumber));
            BuildStepFields();
            AddGroupField();
            ShowTableOf(_pipeSel?.Type);
            RefreshExpressions();
        }
    }

    /// <summary>Selecting a measuring step shows its table in the data inspector.</summary>
    private void ShowTableOf(string? type)
    {
        var name = type switch
        {
            "scatter" => "scatter", "coordination" => "rdf", "cluster" => "clusters", "histogram" => "histogram", "binning" => (_pipeSel != null && double.TryParse(_pipeSel.Params["axis2"]?.ToString(), NumberStyles.Float, CultureInfo.InvariantCulture, out var a2) && a2 >= 0 ? "binning2d" : "binning"),
            "molecule_shape" => "molecules", "cna" => "structures", "ptm" => "structures", "wigner_seitz" => "defects", "wrap" => "outside", "unwrap" => "images", "topology" => "ranges", "voids" => "voids", "voronoi" => (_pipeSel != null && ((string?)_pipeSel.Params["method"] ?? "").StartsWith("exact", StringComparison.Ordinal) ? "voronoi_indices" : "voronoi"), "density_field" => "density_profile",
            "msd" => "msd", "vectors" => "vectors", "displacements" => "displacements", "trajectory_lines" => "paths", "primitive_paths" => "primitive_paths", _ => null,
        };
        if (name == null || _pipeResult?["tables"] is not JsonArray ts) return;
        for (int k = 0; k < ts.Count; ++k)
            if ((string?)ts[k]?["name"] == name) { PipeTable = k; InspectorTab = 3; return; }
    }
    public bool HasPipeSelected => _pipeSel != null;
    public string PipeStepTitle => _pipeSel?.Title ?? "Data source";
    public string PipeStepNumber => _pipeSel == null ? "" : $"step {PipelineRows.Count - PipelineRows.IndexOf(_pipeSel)}";
    public string PipeError { get => _pipeError; private set { if (Set(ref _pipeError, value)) Raise(nameof(HasPipeError)); } }
    public bool HasPipeError => _pipeError.Length > 0;
    public string PipeSourceTitle => _doc == null ? "No file open" : $"{FormatName(_doc.Summary().Format)} · {Title.Replace(" (unsaved)", "")}";
    private static string FormatName(string f) => f switch
    {
        "lammps-dump" => "LAMMPS dump", "lammps-data" => "LAMMPS data", "gro" => "GROMACS .gro", "pdb" => "PDB", "mol2" => "Tripos mol2",
        "cif" => "CIF", "xyz" => "XYZ", "" => "Built in CAPS", _ => f,
    };
    public string PipeSourceDetail => _doc == null ? "" : string.Format(CultureInfo.InvariantCulture, "{0:N0} particles · frame {1} / {2}", _doc.Summary().Atoms, _frame, Math.Max(0, _frames - 1));
    /// <summary>The per-frame cache: this frame from it or computed, and how many frames it holds.</summary>
    public string PipeCacheText => _pipeResult?["cache"] is JsonObject c
        ? $"{((bool?)c["hit"] == true ? "from the cache" : "computed")} · {(int?)c["frames"] ?? 0} frame{((int?)c["frames"] == 1 ? "" : "s")} cached" : "";
    public string PipeStatus { get; private set; } = "";

    public bool StepLibraryOpen { get => _stepLibraryOpen; set { if (Set(ref _stepLibraryOpen, value)) { StepSearch = ""; Raise(nameof(StepLibraryRows)); } } }
    public string StepSearch { get => _stepSearch; set { if (Set(ref _stepSearch, value)) Raise(nameof(StepLibraryRows)); } }
    public IEnumerable<StepKind> StepLibraryRows => StepLibrary.Where(k => _stepSearch.Length == 0 || (k.Title + " " + k.About + " " + k.Group).Contains(_stepSearch, StringComparison.OrdinalIgnoreCase));

    public bool PipeLegendVisible { get => _pipeLegendVisible; set { if (Set(ref _pipeLegendVisible, value)) Raise(nameof(ShowPipeLegend)); } }
    private bool _pipeHasLegend, _pipeLegendContinuous;
    private string _pipeLegendTitle = "", _pipeLegendLo = "", _pipeLegendHi = "";
    public bool ShowPipeLegend => IsVisualize && _pipeHasLegend && _pipeLegendVisible;
    public bool PipeLegendContinuous => _pipeLegendContinuous;
    public bool PipeLegendCategorical => !_pipeLegendContinuous;
    public string PipeLegendTitle => _pipeLegendTitle;
    public string PipeLegendLo => _pipeLegendLo;
    public string PipeLegendHi => _pipeLegendHi;
    public bool PipeLegendDiverging { get; private set; }

    public int InspectorTab { get => _inspectorTab; set { if (Set(ref _inspectorTab, value)) { _inspectorOffset = 0; LoadInspector(); } } }
    public string InspectorFilter { get => _inspectorFilter; set { if (Set(ref _inspectorFilter, value)) { _inspectorOffset = 0; LoadInspector(); } } }
    public string InspectorNote { get => _inspectorNote; private set => Set(ref _inspectorNote, value); }
    public bool InspectorShowsTable => _inspectorTab is 0 or 1;
    public bool IsParticlesTab => _inspectorTab == 0;
    public bool InspectorShowsAttributes => _inspectorTab == 2;
    public bool InspectorShowsTables => _inspectorTab == 3;
    public int PipeTable { get => _pipeTable; set { if (value >= 0 && Set(ref _pipeTable, value)) { _pipeYCol = 1; Raise(nameof(PipeYColumn)); Raise(nameof(PipeTableName)); LoadPipeTable(); } } }
    public ObservableCollection<string> PipeYColumns { get; } = new();
    /// <summary>The table and plotted column by name: the pickers bind to these, so they survive the lists being rebuilt.</summary>
    public string? PipeTableName
    {
        get => _pipeTable >= 0 && _pipeTable < PipeTables.Count ? PipeTables[_pipeTable] : null;
        set { var i = value == null ? -1 : PipeTables.IndexOf(value); if (i >= 0) PipeTable = i; }
    }
    public string? PipeYColumnName
    {
        get => _pipeYCol - 1 >= 0 && _pipeYCol - 1 < PipeYColumns.Count ? PipeYColumns[_pipeYCol - 1] : null;
        set { var i = value == null ? -1 : PipeYColumns.IndexOf(value); if (i >= 0) PipeYColumn = i; }
    }
    private int _pipeYCol = 1;
    /// <summary>Which column the plot shows against the first (index into PipeYColumns + 1).</summary>
    public int PipeYColumn { get => _pipeYCol - 1; set { if (value >= 0 && value + 1 != _pipeYCol) { _pipeYCol = value + 1; Raise(); Raise(nameof(PipeYColumnName)); LoadPipeTable(); } } }
    private JsonObject? _series;
    private bool _seriesRunning;
    public bool SeriesRunning { get => _seriesRunning; private set { if (Set(ref _seriesRunning, value)) Raise(nameof(SeriesIdle)); } }
    public bool SeriesIdle => !_seriesRunning;

    /// <summary>Runs the pipeline on every frame (or every stride-th) and adds the attributes as the Time series table.</summary>
    public async Task ComputeSeries()
    {
        if (_doc == null || SeriesRunning) return;
        var doc = _doc;
        var stride = Math.Max(1, _frames / 500);
        SeriesRunning = true;
        try
        {
            var text = await Task.Run(() => doc.PipelineSeries(stride, f =>
            {
                Avalonia.Threading.Dispatcher.UIThread.Post(() => Status = $"Time series · {f * 100:0} % of the frames");
                return true;
            }));
            _series = JsonNode.Parse(text) as JsonObject;
            RefreshPipeline();
            RefreshTimeline();
            var idx = PipeTables.Count - 1;
            Avalonia.Threading.Dispatcher.UIThread.Post(() =>
            {
                PipeTable = idx;
                InspectorTab = 3;
                if (SparkAttribute != null) PipeYColumnName = SparkAttribute;   // plot the attribute, not the timestep
            });
            Status = $"Time series over {((JsonArray?)_series?["rows"])?.Count ?? 0} frames" + (stride > 1 ? $" (every {stride}th)" : "");
        }
        catch (Exception e) { Status = "Time series failed: " + e.Message; }
        finally { SeriesRunning = false; }
    }
    public bool HasPipeTables => PipeTables.Count > 0;
    public double[] PipeTableX { get; private set; } = [];
    public (double[] X, double[] Y, double[] Z)? PipeTableHeat { get; private set; }
    public string PipeTableZLabel { get; private set; } = "";
    public double[] PipeTableY { get; private set; } = [];
    public string PipeTableXLabel { get; private set; } = "";
    public bool PipeTableScatter { get; private set; }
    public string PipeTableYLabel { get; private set; } = "";
    public event Action? PipeTableChanged;

    // ---------------------------------------------------------------- steps

    public void AddStep(string type, bool select = true)
    {
        var kind = StepLibrary.FirstOrDefault(k => k.Type == type);
        if (kind == null) return;
        var row = new PipelineRow { Type = type, Title = kind.Title, Icon = kind.Icon, Params = DefaultParams(type) };
        if (type == "manual_selection") row.Params["atoms"] = ViewSelectionText();   // what is selected in the view now
        if (_pipeSel?.Group is { Length: > 0 } g) row.Params["group"] = g;   // a new step joins the selected step's group
        WireRow(row);
        // inserted above the selected step, so it runs after it (the list runs bottom to top)
        var at = _pipeSel != null ? PipelineRows.IndexOf(_pipeSel) : 0;
        PipelineRows.Insert(Math.Max(0, at), row);
        StepLibraryOpen = false;
        if (select) PipeSelected = row;
        _showTableAfterApply = select;
        ApplyPipeline();
    }

    public void RemoveStep(PipelineRow row)
    {
        var i = PipelineRows.IndexOf(row);
        if (i < 0) return;
        PipelineRows.Remove(row);
        PipeSelected = PipelineRows.Count == 0 ? null : PipelineRows[Math.Min(i, PipelineRows.Count - 1)];
        ApplyPipeline();
    }

    public void MoveStep(PipelineRow row, int d)
    {
        var i = PipelineRows.IndexOf(row);
        var j = i + d;
        if (i < 0 || j < 0 || j >= PipelineRows.Count) return;
        PipelineRows.Move(i, j);
        Raise(nameof(PipeStepNumber));
        ApplyPipeline();
    }

    private readonly HashSet<string> _foldedGroups = new();

    private void WireRow(PipelineRow row)
    {
        row.Toggled = ApplyPipeline;
        row.GroupToggled = (r, on) =>
        {
            foreach (var x in PipelineRows.Where(x => x.Group == r.Group)) x.Enabled = on;
            ApplyPipeline();
        };
        row.GroupFolded = r =>
        {
            if (!_foldedGroups.Remove(r.Group)) _foldedGroups.Add(r.Group);
            RefreshGroups();
        };
    }

    /// <summary>Headers on the first step of each run of one group, and folded groups hidden behind their header.</summary>
    private void RefreshGroups()
    {
        string? prev = null;
        for (int k = 0; k < PipelineRows.Count; ++k)
        {
            var r = PipelineRows[k];
            r.RaiseGroup();
            var g = r.Group;
            var first = g.Length > 0 && g != prev;
            r.GroupHeader = first;
            if (first)
            {
                var run = PipelineRows.Skip(k).TakeWhile(x => x.Group == g).ToList();
                r.GroupCount = run.Count;
                r.SetGroupOn(run.Any(x => x.Enabled));
                r.Collapsed = _foldedGroups.Contains(g);
            }
            r.RowVisible = g.Length == 0 || !_foldedGroups.Contains(g);
            prev = g;
        }
    }

    public void ClearPipeline()
    {
        PipelineRows.Clear();
        PipeSelected = null;
        ApplyPipeline();
    }

    /// <summary>Writes the steps as JSON (the same file runs with caps pipeline FILE --steps).</summary>
    public void SavePipeline(string path)
    {
        var j = JsonNode.Parse(PipelineJson())!.AsObject();
        j["caps_pipeline"] = 1;
        System.IO.File.WriteAllText(path, j.ToJsonString(new System.Text.Json.JsonSerializerOptions { WriteIndented = true }));
        Status = $"Saved the pipeline ({PipelineRows.Count} steps) to {path}";
    }

    /// <summary>Python steps a loaded pipeline brought in switched off (they run code from the file).</summary>
    public int PipelineHeldPython { get; private set; }

    public void LoadPipeline(string path)
    {
        try
        {
            var text = System.IO.File.ReadAllText(path);
            if (text.TrimStart().StartsWith("caps_pipeline")) text = CapsDocument.PipelineFromYaml(text);   // YAML (SavePipeline)
            var j = JsonNode.Parse(text)!;
            var steps = j is JsonArray a ? a : (JsonArray)j["steps"]!;
            PipelineRows.Clear();
            var held = 0;
            foreach (var st in steps)
            {
                if (st is not JsonObject o) continue;
                var type = (string?)o["type"] ?? "";
                var kind = StepLibrary.FirstOrDefault(k => k.Type == type);
                var prm = (JsonObject)JsonNode.Parse(o.ToJsonString())!;
                var enabled = (bool?)prm["enabled"] ?? true;
                // a Python step runs the file's code (or a .py it names): it arrives switched off until you turn it on
                if (type == "python" && enabled) { enabled = false; ++held; }
                prm.Remove("type");
                prm.Remove("enabled");
                var row = new PipelineRow { Type = type, Title = kind?.Title ?? type, Icon = kind?.Icon ?? "sliders", Params = prm };
                row.Enabled = enabled;
                WireRow(row);
                PipelineRows.Add(row);
            }
            PipelineOutputs.Clear();
            if (j is JsonObject jo && jo["outputs"] is JsonArray outs)
                foreach (var x in outs.OfType<JsonObject>()) { var row = PipelineOutputRow.From(x); WireOutput(row); PipelineOutputs.Add(row); }
            Raise(nameof(HasPipelineOutputs));
            PipeSelected = PipelineRows.FirstOrDefault();
            _showTableAfterApply = true;
            ApplyPipeline();
            Status = $"Loaded {PipelineRows.Count} steps from {path}" +
                     (held == 1 ? " · its Python step is switched off: it runs code from the file — read it, then switch it on"
                      : held > 1 ? $" · its {held} Python steps are switched off: they run code from the file — read them, then switch them on" : "");
            PipelineHeldPython = held;
        }
        catch (Exception e) { Status = "Could not load the pipeline: " + e.Message; }
    }

    private static JsonObject DefaultParams(string type) => type switch
    {
        "select_expression" => new JsonObject { ["expression"] = "Element == \"H\"" },
        "expand_selection" => new JsonObject { ["mode"] = "bonds", ["iterations"] = 1, ["cutoff"] = 3.0 },
        "slice" => new JsonObject { ["normal"] = new JsonArray(0, 0, 1), ["width"] = 12.0, ["invert"] = false, ["select_only"] = false },
        "colour_coding" => new JsonObject { ["property"] = "Molecule", ["mode"] = "auto", ["map"] = "viridis", ["lighten_h"] = true, ["only_selected"] = false },
        "assign_colour" => new JsonObject { ["colour"] = "#E5484D", ["keep_selection"] = false },
        "cluster" => new JsonObject { ["mode"] = "cutoff", ["cutoff"] = 3.3, ["heavy_only"] = true, ["unit"] = "molecules", ["sort_by_size"] = true, ["only_selected"] = false, ["colour"] = true, ["sweep"] = true },
        "coordination" => new JsonObject { ["cutoff"] = 5.0, ["rmax"] = 10.0, ["bins"] = 200, ["element_a"] = 6, ["element_b"] = 6, ["inter_only"] = true, ["only_selected"] = false, ["average_frames"] = false, ["every"] = 1 },
        "create_bonds" => new JsonObject { ["mode"] = "pairs", ["pairs"] = "C-C 1.70, C-H 1.25", ["tolerance"] = 0.45, ["cutoff"] = 1.6, ["keep_file"] = true, ["inter_only"] = false, ["only_selected"] = false },
        "compute_property" => new JsonObject { ["name"] = "Custom", ["expression"] = "Position.Z", ["only_selected"] = false, ["neighbours"] = false, ["neighbour_mode"] = "cutoff", ["cutoff"] = 3.0, ["neighbour_expression"] = "" },
        "wrap" => new JsonObject { ["mode"] = "atoms" },
        "transparency" => new JsonObject { ["mode"] = "selected", ["value"] = 0.7 },
        "particle_radius" => new JsonObject { ["mode"] = "selected", ["value"] = 1.2 },
        "freeze_property" => new JsonObject { ["property"] = "Position.Z", ["frame"] = 0 },
        "cna" => new JsonObject { ["only_selected"] = false },
        "ptm" => new JsonObject { ["rmsd_max"] = 0.1, ["fcc"] = true, ["hcp"] = true, ["bcc"] = true, ["ico"] = true, ["sc"] = false, ["only_selected"] = false },
        "combine" => new JsonObject { ["path"] = "", ["frame"] = 0 },
        "wigner_seitz" => new JsonObject { ["reference"] = "frame", ["frame"] = 0, ["path"] = "", ["output"] = "particles" },
        "centrosymmetry" => new JsonObject { ["neighbours"] = 12 },
        "orientation" => new JsonObject { ["axis"] = "director", ["radius"] = 5.0, ["angle"] = 10.0, ["neighbours"] = 8 },
        "affine_transform" => new JsonObject { ["strain"] = new JsonArray(0.1, 0.0, 0.0), ["target"] = "all" },
        "unwrap" => new JsonObject { ["method"] = "bonds" },
        "replicate" => new JsonObject { ["nx"] = 2, ["ny"] = 2, ["nz"] = 1, ["adjust_cell"] = true },
        "histogram" => new JsonObject { ["property"] = "Charge", ["bins"] = 40, ["stack_by"] = "none", ["only_selected"] = false },
        "molecule_shape" => new JsonObject { ["glyphs"] = true },
        "binning" => new JsonObject { ["property"] = "Mass", ["axis"] = 2, ["bins"] = 50, ["reduction"] = "density", ["axis2"] = -1, ["bins2"] = 50 },
        "topology" => new JsonObject { ["bins"] = 60, ["colour_states"] = true },
        "displacements" => new JsonObject { ["reference"] = "first", ["frame"] = 0 },
        "manual_selection" => new JsonObject { ["atoms"] = "", ["mode"] = "replace" },
        "smooth" => new JsonObject { ["window"] = 5, ["kind"] = "centred", ["unwrap"] = true, ["positions"] = true, ["properties"] = false, ["mark"] = true },
        "vectors" => new JsonObject { ["property"] = "end_to_end", ["scale"] = 1.0, ["radius"] = 0.3 },
        "python" => new JsonObject { ["file"] = "", ["code"] = PythonStepTemplate },
        "msd" => new JsonObject { ["heavy_only"] = true, ["every"] = 1, ["timestep_fs"] = 1.0 },
        "scatter" => new JsonObject { ["x"] = "DistanceToCOM", ["y"] = "Charge", ["only_selected"] = false },
        "voids" => new JsonObject { ["probe"] = 1.4, ["grid"] = 0.5, ["show"] = true },
        "voronoi" => new JsonObject { ["method"] = "exact", ["grid"] = 0.5, ["face_area_min"] = 0.0, ["only_selected"] = false },
        "density_field" => new JsonObject { ["grid"] = 0.8, ["sigma"] = 1.5, ["axis"] = 2, ["position"] = 0.5, ["slice"] = true, ["isosurface"] = false, ["opacity"] = 0.45, ["colour"] = "#6FA8DC" },
        "trajectory_lines" => new JsonObject { ["particles"] = "centres", ["from"] = 0, ["radius"] = 0.12 },
        "primitive_paths" => new JsonObject { ["radius"] = 0.3, ["show_chains"] = false },
        _ => new JsonObject(),
    };

    /// <summary>A Python step to start from: what a step reads and writes.</summary>
    public const string PythonStepTemplate =
        "from caps.pipeline import step\n\n" +
        "@step(name=\"Heavy atoms per molecule\")\n" +
        "def modify(frame, data):\n" +
        "    # data.particles: one list per property (\"Element\", \"Position\", \"Molecule Identifier\", \"Charge\", \"Backbone\", …)\n" +
        "    heavy = [1 if e != \"H\" else 0 for e in data.particles[\"Element\"]]\n" +
        "    per = {}\n" +
        "    for m, h in zip(data.particles[\"Molecule Identifier\"], heavy):\n" +
        "        per[int(m)] = per.get(int(m), 0) + h\n" +
        "    data.particles[\"Heavy\"] = heavy          # a new property, for colour coding and expressions\n" +
        "    data.attributes[\"Heavy atoms\"] = sum(heavy)\n" +
        "    data.tables[\"heavy_per_molecule\"] = per\n" +
        "    print(f\"frame {frame}: {sum(heavy)} heavy atoms in {len(per)} molecules\")   # shows in the console\n";

    /// <summary>The fields of the selected step's editor, from its parameters.</summary>
    private void BuildStepFields()
    {
        StepFields.Clear();
        if (_pipeSel == null) return;
        var p = _pipeSel.Params;
        string S(string k, string d = "") => p[k] is JsonValue v ? (v.TryGetValue<string>(out var s) ? s : v.ToJsonString()) : d;
        bool B(string k, bool d = false) => p[k] is JsonValue v ? v.TryGetValue<bool>(out var b) && b : d;
        var props = PipeProperties();
        void Add(StepField f)
        {
            f.Changed = () => WriteField(f);
            StepFields.Add(f);
        }
        void Text(string key, string label, string kind = "text", string hint = "", string d = "") => Add(new StepField { Key = key, Label = label, Kind = kind, Hint = hint, Text = S(key, d) });
        void Note(string text) => StepFields.Add(new StepField { Kind = "note", Hint = text });
        void Bool(string key, string label, bool d = false) => Add(new StepField { Key = key, Label = label, Kind = "bool", On = B(key, d) });
        void Choice(string key, string label, string[] choices)
        {
            var value = S(key, choices[0]);
            if (!choices.Contains(value)) choices = [.. choices, value];   // a property made by a step, before its first run
            Add(new StepField { Key = key, Label = label, Kind = "choice", Choices = choices, Text = value });
        }
        switch (_pipeSel.Type)
        {
            case "select_expression": Text("expression", "Expression", "expression", "Type == 2 && Position.Z > 13 · Element == \"O\""); Bool("legend", "Show the selection count in the legend"); break;
            case "expand_selection": Choice("mode", "Across", ["bonds", "cutoff"]); Text("iterations", "Steps", "number"); Text("cutoff", "Cutoff (Å)", "number"); break;
            case "slice":
                Add(new StepField { Key = "normal", Label = "Normal", Kind = "vector", Hint = "x y z", Text = p["normal"] is JsonArray a ? string.Join(" ", a.Select(x => x?.ToString())) : "0 0 1" });
                Text("distance", "Distance (Å)", "number", "blank: through the cell centre"); Text("width", "Width (Å)", "number");
                Bool("invert", "Invert (remove the slab)"); Bool("select_only", "Select instead of delete"); break;
            case "colour_coding":
                Choice("property", "Property", props); Choice("mode", "Map", ["auto", "categorical", "continuous"]); Choice("map", "Continuous colours", ["viridis", "diverging"]);
                Text("start", "Start", "number", "blank: minimum"); Text("end", "End", "number", "blank: maximum");
                Bool("lighten_h", "Lighten hydrogens"); Bool("only_selected", "Only selected"); break;
            case "assign_colour": Text("colour", "Colour (#RRGGBB)", "text"); Bool("keep_selection", "Keep selection"); break;
            case "cluster":
                Choice("mode", "Neighbours", ["bonds", "cutoff"]); Text("cutoff", "Cutoff (Å)", "number"); Choice("unit", "Unit", ["atoms", "molecules"]);
                Bool("heavy_only", "Cutoff between heavy atoms"); Bool("sort_by_size", "Sort by size"); Bool("colour", "Colour by cluster"); Bool("sweep", "Cutoff sweep (cutoff mode)"); Bool("only_selected", "Only selected");
                Bool("periodic", "Periodic (links across the cell's faces, minimum image)", true); Bool("unwrap", "Unwrap clusters (made whole through their links, shown whole)"); break;
            case "coordination":
                Text("element_a", "A · element number (0: any)", "number"); Text("element_b", "B · element number (0: any)", "number"); Bool("inter_only", "Only different molecules");
                Text("rmax", "g(r) out to r max (Å)", "number"); Text("bins", "Bins", "number"); Text("cutoff", "Coordination cutoff (Å)", "number");
                Bool("average_frames", "Average g(r) over the frames"); Text("every", "Every n-th frame", "number"); Bool("only_selected", "Only selected"); break;
            case "topology": Text("bins", "Bins", "number");
                Bool("bonds", "Bonds", true); Bool("angles", "Angles", true); Bool("dihedrals", "Dihedrals", true);
                Choice("dihedral_set", "Dihedrals over", ["all", "backbone"]); Text("trans", "Trans when |φ| above (°)", "number", "blank: 120");
                Bool("colour_states", "Colour the backbone by dihedral state (t · g+ · g−)"); break;
            case "displacements": Choice("reference", "Reference", ["first", "previous", "frame"]); Text("frame", "Reference frame", "number"); Bool("subtract_drift", "Subtract system drift"); break;
            case "manual_selection":
                Text("atoms", "Atoms (indices in the frame: 0 5 12-40)", "text", "click, ⇧-click or lasso atoms in the view");
                Add(new StepField { Key = "use_view", Label = "Use the view's selection", Kind = "action" });
                Choice("mode", "Mode", ["replace", "add", "subtract"]);
                Note("The same atoms on every frame. Lasso: the view toolbar's lasso tool (⇧ adds).");
                break;
            case "smooth":
                Text("window", "Window (frames)", "number"); Choice("kind", "Kind", ["centred", "trailing", "gaussian"]);
                Bool("unwrap", "Unwrap before averaging (each atom followed to this frame's image)", true);
                Bool("positions", "Apply to positions", true); Bool("properties", "Apply to per-atom properties (the steps below, run on each frame of the window)");
                Bool("mark", "Mark frames as smoothed", true);
                Note("Averaging removes vibration, but chains that turn within the window average to slightly shrunken shapes; smoothed frames are marked so bond and angle analysis can tell.");
                break;
            case "python":
                Add(new StepField { Key = "file", Label = "Script (.py with an @step function)", Kind = "file", Hint = "blank: the step typed below", Text = S("file") });
                var code = S("code", PythonStepTemplate);
                Add(new StepField { Key = "code", Label = "Or type the step (Run sends it; a file above takes precedence)", Kind = "code", Text = code, Draft = code, Output = _pipeSel.Output.TrimEnd() });
                break;
            case "msd":
                Bool("heavy_only", "Heavy atoms only"); Text("every", "Every n-th atom", "number"); Text("max_lag", "Longest lag (frames)", "number", "blank: half the frames");
                Text("fit_from", "Fit from lag (frames)", "number", "blank: a quarter of the lags"); Text("fit_to", "Fit to lag (frames)", "number", "blank: the longest");
                Text("timestep_fs", "Timestep (fs) for D in cm²/s", "number");
                Bool("remove_drift", "Remove system drift (paths relative to the system's centre of mass)");
                Bool("per_molecule", "Per molecule curves (the first twelve centres)"); break;
            case "scatter": Choice("x", "x", props); Choice("y", "y", props); Bool("only_selected", "Only selected"); break;
            case "voids": Text("probe", "Probe radius (Å)", "number"); Text("grid", "Grid (Å)", "number"); Text("min_volume", "Only voids above (Å³)", "number", "blank or 0: every void");
                Bool("show", "Show void points, coloured by void"); Note("The void surface as a mesh: Free volume page"); break;
            case "voronoi":
                Choice("method", "Method", ["exact", "exact_radical", "grid", "radical"]);
                Text("face_area_min", "Smallest face counted (Å²)", "number", "exact: faces below this stay in the volume, not in the index");
                Text("grid", "Grid (Å)", "number", "grid methods only");
                Bool("only_selected", "Only selected (exact)");
                Bool("periodic_x", "Periodic x (exact)", true); Bool("periodic_y", "Periodic y (exact)", true); Bool("periodic_z", "Periodic z (exact)", true);
                Note("Exact: each atom's cell clipped by the bisecting planes of its neighbours (radical: the power planes weighted by van der Waals radii, Gellatly & Finney 1982) — AtomicVolume, Coordination (faces), Max Face Order, Cell Surface Area and Voronoi Index.3–6; the cells tile the box exactly. Grid: voxels to the nearest atom");
                break;
            case "density_field":
                Text("grid", "Grid (Å)", "number"); Text("sigma", "Smoothing σ (Å)", "number"); Choice("axis", "Slice normal (0 x · 1 y · 2 z)", ["0", "1", "2"]);
                Text("position", "Slice position (0–1 of the cell)", "number"); Bool("slice", "Show the slice", true);
                Bool("isosurface", "Isosurface (the boundary of the denser regions)");
                Text("level", "Level (g/cm³)", "number", "blank: the cell's mean density"); Text("opacity", "Opacity (0–1)", "number"); Text("colour", "Colour (#RRGGBB)", "text");
                Note("Above the level: the volume share, area and specific area are global attributes; a low level traces the free-volume network, the mean one a blend's domains.");
                break;
            case "vectors":
                Choice("property", "Vector", ["end_to_end", "dipole", "displacement", "velocity", "force", "columns"]); Text("scale", "Scale (dipole, displacement, velocity, force)", "number");
                Text("x", "Columns: x", "text", "the dump's per-atom columns (Vector: columns)"); Text("y", "Columns: y", "text"); Text("z", "Columns: z", "text");
                Text("radius", "Arrow radius (Å)", "number"); Bool("arrowheads", "Arrowheads", true);
                Choice("anchor", "The particle at the arrow's", ["tail", "middle", "head"]);
                Choice("colour_by", "Colour by", ["auto", "length", "fixed"]); Bool("flip", "Flip direction (end-to-end)");
                Note("Force: the dump's fx fy fz columns (dump … fx fy fz), kcal/mol/Å × scale; auto colour: per molecule, or the particle's colour coding"); break;
            case "trajectory_lines":
                Choice("particles", "Trace", ["centres", "selected"]); Text("from", "From frame", "number"); Text("to", "To frame", "number", "blank: the last");
                Text("stride", "Every n-th frame", "number", "blank: about 200 steps"); Text("radius", "Line radius (Å)", "number");
                Choice("colour", "Colour", ["molecule", "time"]); Bool("fade", "Fade with age"); Bool("up_to_current", "Show only up to the current frame"); break;
            case "create_bonds":
                Choice("mode", "Mode", ["pairs", "perceive", "cutoff"]); Text("pairs", "Cutoff by pair (Å)", "text", "C-C 1.70, C-H 1.25 (a pair not listed: never bonded)");
                Text("tolerance", "Tolerance over covalent radii (Å)", "number"); Text("cutoff", "One cutoff (Å)", "number");
                Bool("keep_file", "Keep file bonds (compare with them)"); Bool("inter_only", "Only between different molecules"); Bool("replace", "Replace the bonds"); Bool("only_selected", "Only selected"); break;
            case "compute_property":
                Text("name", "Output property"); Text("expression", "Expression", "expression", "e.g. sqrt(Position.X^2 + Position.Y^2)"); Bool("only_selected", "Only selected");
                Bool("neighbours", "Neighbour terms (added: Σ over the neighbours)");
                Choice("neighbour_mode", "Neighbours", ["cutoff", "bonds"]); Text("cutoff", "Cutoff r_c (Å)", "number");
                Text("neighbour_expression", "Neighbour term", "expression", "e.g. Charge / Distance · exp(-Distance^2) · Element == 8");
                Note("In the neighbour term, names are the neighbour's properties; Distance and Delta.X/.Y/.Z are the pair's (from the atom to the neighbour, minimum image).");
                break;
            case "transparency" or "particle_radius":
                var what = _pipeSel.Type == "transparency" ? "Transparency (0 opaque … 1 invisible)" : "Radius (Å)";
                Choice("mode", "Set", ["selected", "property"]);
                if (S("mode", "selected") == "selected") { Text("value", what, "number"); Bool("keep_selection", "Keep the selection"); Note("Select the particles with a step below this one (e.g. Expression selection: Element == \"Si\")"); }
                else
                {
                    Choice("property", "Property", props); Text("start", "From value", "number", "blank: minimum"); Text("end", "To value", "number", "blank: maximum");
                    Text("low", what + " at the start", "number"); Text("high", what + " at the end", "number"); Bool("only_selected", "Only selected");
                }
                break;
            case "freeze_property":
                Choice("property", "Property", props); Text("frame", "Reference frame", "number"); Text("output", "Output property", "text", "blank: <property> frozen");
                Note("The steps below this one run on the reference frame; particles are matched by identifier. Colour by the frozen property to follow where atoms started."); break;
            case "ptm":
                Text("rmsd_max", "RMSD cutoff", "number", "larger admits more distorted neighbourhoods", "0.1");   // the core's defaults when a loaded step leaves them out
                Bool("fcc", "FCC", true); Bool("hcp", "HCP", true); Bool("bcc", "BCC", true); Bool("ico", "Icosahedral", true); Bool("sc", "Simple cubic"); Bool("only_selected", "Only selected");
                Note("Larsen, Schmidt & Schiøtz, Modelling Simul. Mater. Sci. Eng. 24, 055007 (2016). Adds Structure Type (0 other, 1 FCC, 2 HCP, 3 BCC, 4 icosahedral, 5 SC), RMSD, Interatomic Distance, Shear Strain and Orientation.W/X/Y/Z (lattice → template quaternion); the RMSD table shows where to set the cutoff");
                break;
            case "wigner_seitz":
                Choice("reference", "Reference lattice", ["frame", "file"]);
                Text("frame", "Reference frame", "number");
                Text("path", "Reference file", "text", "the perfect lattice (reference: file)");
                Choice("output", "Output", ["particles", "sites"]);
                Note("Each particle belongs to its nearest reference site (minimum image): an empty site is a vacancy, each particle beyond the first on a site an interstitial. Particles: Site Index and Occupancy of their site; sites: the reference sites with their Occupancy (0 = vacancy) — select Occupancy == 0 to see the vacancies");
                break;
            case "combine":
                Text("path", "File", "text", "a LAMMPS data/dump, .gro, .pdb, .xyz, mol2 or CIF file");
                Text("frame", "Frame of that file", "number");
                Note("The second file's particles and bonds are appended as they are in it (same coordinates); they are not selected and take zeros for properties set above");
                break;
            case "cna": Bool("only_selected", "Only selected"); Note("Adds Structure Type (0 other, 1 FCC, 2 HCP, 3 BCC, 4 icosahedral): colour by it to see grains, stacking faults and surfaces"); break;
            case "centrosymmetry": Text("neighbours", "Nearest neighbours (12 FCC, 8 BCC)", "number"); Note("Adds Centrosymmetry (Å²): zero in a perfect lattice, large at surfaces and defects"); break;
            case "orientation":
                Choice("axis", "P₂ against", ["director", "x", "y", "z"]); Text("radius", "Neighbour radius (Å)", "number"); Text("angle", "Aligned within (°)", "number");
                Text("neighbours", "Aligned neighbours for crystalline", "number");
                Note("Adds Orientation (P₂ of the backbone chord through each atom) and Crystalline (1/0) — colour by either"); break;
            case "affine_transform":
                Add(new StepField { Key = "strain", Label = "Strain (εxx εyy εzz)", Kind = "vector", Hint = "0.1 0 0", Text = p["strain"] is JsonArray st ? string.Join(" ", st.Select(x => x?.ToString())) : "" });
                Add(new StepField { Key = "matrix", Label = "Matrix (9 numbers, row by row; blank: identity)", Kind = "matrix", Hint = "1 0.2 0  0 1 0  0 0 1",
                                    Text = p["matrix"] is JsonArray mt ? string.Join(" ", mt.Select(x => x?.ToString())) : "" });
                Choice("target", "Transform", ["all", "particles", "cell"]); Bool("only_selected", "Only selected particles"); break;
            case "wrap": Choice("mode", "Fold", ["atoms", "molecules"]); Note("atoms: each atom into the cell (bonds cross faces) · molecules: each molecule whole, its centre of mass in the cell"); break;
            case "unwrap":
                Choice("method", "Method", ["bonds", "images", "nojump"]);
                Note("bonds: molecules whole along their bonds · images: the file's image flags · nojump: each atom followed through the frames (MSD, diffusion)"); break;
            case "replicate": Text("nx", "Images along a", "number"); Text("ny", "Images along b", "number"); Text("nz", "Images along c", "number"); Bool("adjust_cell", "Enlarge the cell"); Bool("unique_ids", "Unique IDs (the copies numbered after the original)", true); Bool("bonds", "Operate on bonds", true); Bool("vectors", "Operate on arrows, paths and glyphs", true); break;
            case "primitive_paths": Bool("show_chains", "Show the chains too"); Text("radius", "Line radius (Å)", "number"); Text("max_steps", "Minimisation steps at most", "number", "blank: 200 000"); break;
            case "molecule_shape": Choice("group_by", "Group by", ["molecule", "unit"]); Choice("atoms", "Atoms", ["all", "heavy", "backbone"]);
                Bool("glyphs", "Principal-axis glyphs", true); Bool("dim", "Dim chains (the glyphs stand out)"); break;
            case "histogram": Choice("property", "Property", props); Text("bins", "Bins", "number"); Choice("stack_by", "Stack by", ["none", "Type", "Element", "Molecule"]); Bool("only_selected", "Only selected"); break;
            case "binning":
                Choice("property", "Property", props); Choice("axis", "Along (0 x, 1 y, 2 z)", ["0", "1", "2"]); Text("bins", "Bins", "number"); Choice("reduction", "Reduction", ["density", "mean", "sum"]);
                Choice("axis2", "Second axis: a map (−1 none)", ["-1", "0", "1", "2"]); Text("bins2", "Bins along the second axis", "number");
                Note("With a second axis the result is a map, one cell per pair of bins, drawn as a heat map in the data table view (density: each cell's column through the box)");
                break;
        }
    }

    private void AddGroupField()
    {
        if (_pipeSel == null) return;
        var f = new StepField { Key = "group", Label = "Group (steps with the same name fold together)", Kind = "text", Hint = "e.g. Prepare · Structure · Look",
                                Text = _pipeSel.Group };
        f.Changed = () => WriteField(f);
        StepFields.Add(f);
        var b = new StepField { Key = "branch", Label = "Branch (only in that view; the shared steps below run once)", Kind = "text", Hint = "blank: shared · e.g. Figure · Numbers",
                                Text = _pipeSel.Params["branch"] is JsonValue v && v.TryGetValue<string>(out var br) ? br : "" };
        b.Changed = () => WriteField(b);
        StepFields.Add(b);
    }

    // branches (design/boards/PipelineGroups): steps with a "branch" run only in that branch's view; the shared steps
    // below them run once per frame and every branch starts from their cached result
    private string _pipeBranch = "";
    private bool _hadBranches;
    public List<string> PipeBranches => PipelineRows.Select(r => r.Params["branch"] is JsonValue v && v.TryGetValue<string>(out var b) ? b : "")
                                                    .Where(b => b.Length > 0).Distinct().Prepend("Shared steps only").ToList();
    public bool HasPipeBranches => PipeBranches.Count > 1;
    public int PipeBranchIndex
    {
        get => Math.Max(0, PipeBranches.IndexOf(_pipeBranch.Length == 0 ? "Shared steps only" : _pipeBranch));
        set
        {
            var list = PipeBranches;
            var b = value <= 0 || value >= list.Count ? "" : list[value];
            if (b == _pipeBranch) return;
            _pipeBranch = b;
            Raise();
            ApplyPipeline();
        }
    }
    private void RefreshBranches()
    {
        // the shown branch goes when its last step leaves it; a first branch is shown as soon as it exists
        var list = PipeBranches;
        if (_pipeBranch.Length > 0 && !list.Contains(_pipeBranch)) _pipeBranch = "";
        if (_pipeBranch.Length == 0 && list.Count > 1 && !_hadBranches) _pipeBranch = list[1];
        _hadBranches = list.Count > 1;
        Raise(nameof(PipeBranches)); Raise(nameof(HasPipeBranches)); Raise(nameof(PipeBranchIndex));
    }

    internal string[] PipeProperties()
    {
        var list = new List<string> { "Molecule", "Type", "Element", "Charge", "Mass", "Position.X", "Position.Y", "Position.Z", "DistanceToCOM", "Identifier", "Selection",
                                      "Cluster", "Coordination", "Displacement", "MoleculeRg", "MoleculeKappa2" };
        if (_pipeResult?["properties"] is JsonArray a)
            foreach (var x in a) if (x?.GetValue<string>() is { } s && !list.Contains(s)) list.Add(s);
        return list.ToArray();
    }

    private void WriteField(StepField f)
    {
        if (_pipeSel == null) return;
        var p = _pipeSel.Params;
        switch (f.Kind)
        {
            case "bool": p[f.Key] = f.On; break;
            case "number":
                if (string.IsNullOrWhiteSpace(f.Text)) p.Remove(f.Key);
                else if (double.TryParse(f.Text, NumberStyles.Float, CultureInfo.InvariantCulture, out var v)) p[f.Key] = v;
                else return;
                break;
            case "vector":
                var parts = f.Text.Split([' ', ',', ';'], StringSplitOptions.RemoveEmptyEntries).Select(x => double.TryParse(x, NumberStyles.Float, CultureInfo.InvariantCulture, out var d) ? d : double.NaN).ToArray();
                if (parts.Length != 3 || parts.Any(double.IsNaN)) return;
                p[f.Key] = new JsonArray(parts.Select(x => (JsonNode)x).ToArray());
                break;
            case "choice" when f.Key == "axis" && int.TryParse(f.Text, NumberStyles.Integer, CultureInfo.InvariantCulture, out var axis): p[f.Key] = axis; break;
            case "matrix":
                var m = f.Text.Split([' ', ',', ';'], StringSplitOptions.RemoveEmptyEntries).Select(x => double.TryParse(x, NumberStyles.Float, CultureInfo.InvariantCulture, out var d) ? d : double.NaN).ToArray();
                if (m.Length == 0) { p.Remove(f.Key); break; }
                if (m.Length != 9 || m.Any(double.IsNaN)) return;
                p[f.Key] = new JsonArray(m.Select(x => (JsonNode)x).ToArray());
                break;
            case "text" when f.Key == "group" && f.Text.Trim().Length == 0: p.Remove("group"); break;
            case "text" when f.Key == "branch" && f.Text.Trim().Length == 0: p.Remove("branch"); break;
            case "text" when f.Key == "branch": p["branch"] = f.Text.Trim(); break;
            default: p[f.Key] = f.Text; break;
        }
        ApplyPipeline();
        if (f.Key == "mode" && _pipeSel.Type is "transparency" or "particle_radius")   // other fields for the other mode
            Avalonia.Threading.Dispatcher.UIThread.Post(BuildStepFields);
    }

    /// <summary>The view's selection (lassoed and clicked atoms, frame indices) as text with ranges: "0 5 12-40".</summary>
    public string ViewSelectionText()
    {
        var set = new SortedSet<int>(_selection);
        try
        {
            if (_doc != null && JsonNode.Parse(_doc.SelectionJson())?["indices"] is JsonArray a)
                foreach (var x in a) set.Add((int)((double?)x ?? -1));
        }
        catch (Exception) { }
        set.Remove(-1);
        var parts = new List<string>();
        int? start = null, prev = null;
        foreach (var i in set.Append(int.MinValue))
        {
            if (start != null && i == prev + 1) { prev = i; continue; }
            if (start != null) parts.Add(prev == start ? $"{start}" : prev == start + 1 ? $"{start} {prev}" : $"{start}-{prev}");
            start = prev = i;
        }
        return string.Join(" ", parts);
    }

    /// <summary>A step's action button: Manual selection takes the view's selection.</summary>
    public void StepAction(StepField f)
    {
        if (f.Key != "use_view") return;
        var text = ViewSelectionText();
        if (StepFields.FirstOrDefault(x => x.Key == "atoms") is { } atoms) atoms.Text = text;
        Status = text.Length == 0 ? "Nothing selected in the view: click, ⇧-click or lasso atoms first" : $"Manual selection: {text.Split(' ').Length} ranges from the view";
    }

    public string PipelineJson() => new JsonObject
    {
        ["steps"] = new JsonArray(PipelineRows.Select(r =>
        {
            var o = (JsonObject)JsonNode.Parse(r.Params.ToJsonString())!;
            o["type"] = r.Type;
            o["enabled"] = r.Enabled;
            return (JsonNode)o;
        }).ToArray()),
        ["outputs"] = PipelineOutputsJson(),
        ["branch"] = _pipeBranch,
    }.ToJsonString();

    /// <summary>Sends the steps to the core, which runs them on the shown frame; then the list, legend and inspector refresh.</summary>
    public void ApplyPipeline()
    {
        if (_doc == null) return;
        ++_pipeGen;
        RefreshBranches();
        try
        {
            _doc.SetPipeline(PipelineRows.Count == 0 ? "" : PipelineJson());
            PipeError = "";
        }
        catch (Exception e) { PipeError = e.Message; }
        RefreshPipeline();
        RefreshGroups();
        RenderRequested?.Invoke();
    }

    /// <summary>Clears the pipeline in the core (leaving Visualize keeps the steps for when it opens again).</summary>
    private void SuspendPipeline()
    {
        try { _doc?.SetPipeline(""); } catch { }
        _pipeHasLegend = false;
        Raise(nameof(ShowPipeLegend));
    }

    /// <summary>A step made a grid (Density field): it can be written as .cube, .vtk or .npy.</summary>
    public bool PipeHasGrid => _pipeResult?["grid"] != null;
    public void ExportGrid(string path)
    {
        if (_doc == null) return;
        try { _doc.PipelineExportGrid(path); Status = $"Grid written: {System.IO.Path.GetFileName(path)}"; }
        catch (Exception e) { Status = "Grid export: " + e.Message; }
    }

    private void RefreshPipeline()
    {
        PipeAttributes.Clear();
        _pipeResult = null;
        if (_doc == null) return;
        string text;
        try { text = _doc.PipelineResult(); } catch { text = ""; }
        _pipeResult = text.Length > 0 ? JsonNode.Parse(text) : null;
        Raise(nameof(PipeHasGrid)); Raise(nameof(PipeCacheText));
        var inv = CultureInfo.InvariantCulture;
        if (_pipeResult?["steps"] is JsonArray steps)
            for (int k = 0; k < steps.Count && k < PipelineRows.Count; ++k)
            {
                PipelineRows[k].Summary = (string?)steps[k]?["summary"] ?? "";
                PipelineRows[k].Level = (string?)steps[k]?["level"] ?? "ok";
                PipelineRows[k].Output = (string?)steps[k]?["output"] ?? "";
            }
        if (_pipeSel != null)
            foreach (var f in StepFields.Where(f => f.IsCode)) f.Output = _pipeSel.Output.TrimEnd();
        if (_pipeResult?["attributes"] is JsonArray attrs)
            foreach (var a in attrs)
            {
                var name = (string?)a?["name"] ?? "";
                var v = (double?)a?["value"] ?? 0;
                var text2 = name switch
                {
                    "CellVolume" => v.ToString("N0", inv) + " Å³",
                    "Mass" => v.ToString("N1", inv) + " g/mol",
                    "Density" => v.ToString("F3", inv) + " g/cm³",
                    "TotalCharge" => v.ToString("+0.00e+0;-0.00e+0;0", inv) + " e",
                    _ => Math.Abs(v - Math.Round(v)) < 1e-9 && Math.Abs(v) < 1e15 ? v.ToString("0", inv) : v.ToString("G6", inv),
                };
                PipeAttributes.Add(new Row(name, text2));
            }
        // legend
        _pipeHasLegend = _pipeResult?["legend"] is JsonObject;
        PipeLegendEntries.Clear();
        if (_pipeResult?["legend"] is JsonObject l)
        {
            _pipeLegendContinuous = (bool?)l["continuous"] ?? false;
            _pipeLegendTitle = (string?)l["property"] ?? "";
            _pipeLegendLo = ((double?)l["lo"] ?? 0).ToString("G4", inv);
            _pipeLegendHi = ((double?)l["hi"] ?? 0).ToString("G4", inv);
            PipeLegendDiverging = (string?)l["map"] == "diverging";
            if (l["entries"] is JsonArray e)
                foreach (var x in e.Take(12))
                    PipeLegendEntries.Add(new LegendSwatch((_pipeLegendTitle == "Molecule" ? "mol " : "") + ((string?)x?["label"] ?? ""),
                                                           new SolidColorBrush(Color.Parse((string?)x?["colour"] ?? "#888888"))));
        }
        foreach (var n in new[] { nameof(ShowPipeLegend), nameof(PipeLegendContinuous), nameof(PipeLegendCategorical), nameof(PipeLegendTitle),
                                  nameof(PipeLegendLo), nameof(PipeLegendHi), nameof(PipeLegendDiverging), nameof(PipeSourceTitle), nameof(PipeSourceDetail) })
            Raise(n);
        // tables
        // the titles are rebuilt only when they change, so the picker keeps its choice
        var titles = _pipeResult?["tables"] is JsonArray tables ? tables.Select(t => (string?)t?["title"] ?? "").ToList() : new List<string>();
        if (_series?["rows"] is JsonArray sr) titles.Add($"Time series · {sr.Count} frames");
        if (!titles.SequenceEqual(PipeTables))
        {
            var keep = _pipeTable;
            PipeTables.Clear();
            foreach (var t in titles) PipeTables.Add(t);
            _pipeTable = Math.Clamp(keep, 0, Math.Max(0, PipeTables.Count - 1));
            Avalonia.Threading.Dispatcher.UIThread.Post(() => { Raise(nameof(PipeTable)); Raise(nameof(PipeTableName)); });
        }
        Raise(nameof(HasPipeTables));
        LoadPipeTable();
        if (_showTableAfterApply) { _showTableAfterApply = false; ShowTableOf(_pipeSel?.Type); }
        var particles = (double?)_pipeResult?["particles"] ?? _doc.Summary().Atoms;
        PipeStatus = string.Format(inv, "{0:N0} particles · {1} step{2}", particles, PipelineRows.Count, PipelineRows.Count == 1 ? "" : "s");
        Raise(nameof(PipelineAtFrameText));
        TimelineChanged?.Invoke();
        Raise(nameof(PipeStatus));
        LoadInspector();
    }

    private void LoadPipeTable()
    {
        PipeTableX = []; PipeTableY = []; PipeTableHeat = null;
        var coreTables = _pipeResult?["tables"] as JsonArray;
        var ncore = coreTables?.Count ?? 0;
        var t = _pipeTable < ncore ? coreTables![_pipeTable] as JsonObject : _pipeTable == ncore ? _series : null;
        if (t != null && t["rows"] is JsonArray rows && t["columns"] is JsonArray cols && cols.Count >= 2)
        {
            var ys = cols.Skip(1).Select(c => (string?)c ?? "").ToList();
            if (!ys.SequenceEqual(PipeYColumns))
            {
                PipeYColumns.Clear();
                foreach (var y in ys) PipeYColumns.Add(y);
                Avalonia.Threading.Dispatcher.UIThread.Post(() => { Raise(nameof(PipeYColumn)); Raise(nameof(PipeYColumnName)); });
            }
            var yc = Math.Clamp(_pipeYCol, 1, cols.Count - 1);
            PipeTableX = rows.Select(r => (double?)r?[0] ?? 0).ToArray();
            PipeTableY = rows.Select(r => (double?)r?[yc] ?? double.NaN).ToArray();
            PipeTableXLabel = (string?)cols[0] ?? "";
            PipeTableYLabel = (string?)cols[yc] ?? "";
            PipeTableScatter = (bool?)t["points"] ?? false;
            // a map over two axes (Spatial binning, second axis): x, y and the chosen value column as a heat map
            PipeTableHeat = null;
            if ((string?)t["name"] == "binning2d" && cols.Count >= 3)
            {
                var zc = Math.Clamp(Math.Max(_pipeYCol, 2), 2, cols.Count - 1);
                PipeTableHeat = (rows.Select(r => (double?)r?[0] ?? 0).ToArray(), rows.Select(r => (double?)r?[1] ?? 0).ToArray(),
                                 rows.Select(r => (double?)r?[zc] ?? double.NaN).ToArray());
                PipeTableYLabel = (string?)cols[1] ?? "";
                PipeTableZLabel = (string?)cols[zc] ?? "";
            }
            if (_inspectorTab == 3)
            {
                // a text column (the molecules a cluster holds …) goes after the first
                var labels = t["labels"] as JsonArray;
                InspectorColumns.Clear();
                for (var c = 0; c < cols.Count; ++c)
                {
                    InspectorColumns.Add((string?)cols[c] ?? "");
                    if (c == 0 && labels != null) InspectorColumns.Add((string?)t["label_column"] ?? "");
                }
                InspectorRows.Clear();
                var inv = CultureInfo.InvariantCulture;
                var ri = 0;
                foreach (var r in rows.Take(InspectorPage))
                {
                    var cells = ((JsonArray)r!).Select(x => x is null ? "–" : ((double?)x ?? 0).ToString("G5", inv)).ToList();
                    if (labels != null) cells.Insert(1, (string?)labels[ri] ?? "");
                    InspectorRows.Add(new TableRow(cells.ToArray(), false));
                    ++ri;
                }
                InspectorNote = $"{rows.Count} rows · {t["title"]}";
            }
        }
        Raise(nameof(PipeTableXLabel));
        Raise(nameof(PipeTableYLabel));
        PipeTableChanged?.Invoke();
    }

    /// <summary>Every row of the inspector's tab as CSV: the particles matching the filter, the bonds, the attributes, or
    /// the chosen data table (all its rows, not the page shown).</summary>
    public void ExportInspectorCsv(string path)
    {
        if (_doc == null) return;
        var inv = CultureInfo.InvariantCulture;
        static string Q(string v) => v.IndexOfAny([',', '"', '\n', '\r']) >= 0 ? "\"" + v.Replace("\"", "\"\"") + "\"" : v;
        var sb = new System.Text.StringBuilder();
        var n = 0;
        try
        {
            if (_inspectorTab == 2)
            {
                sb.Append("Attribute,Value\n");
                foreach (var a in PipeAttributes) { sb.Append(Q(a.Key)).Append(',').Append(Q(a.Value)).Append('\n'); ++n; }
            }
            else if (_inspectorTab == 3)
            {
                var coreTables = _pipeResult?["tables"] as JsonArray;
                var ncore = coreTables?.Count ?? 0;
                var t = _pipeTable < ncore ? coreTables![_pipeTable] as JsonObject : _pipeTable == ncore ? _series : null;
                if (t?["rows"] is not JsonArray rows || t["columns"] is not JsonArray cols) { Status = "No table to export"; return; }
                var labels = t["labels"] as JsonArray;
                var head = cols.Select(c => (string?)c ?? "").ToList();
                if (labels != null) head.Insert(1, (string?)t["label_column"] ?? "Label");
                sb.Append(string.Join(",", head.Select(Q))).Append('\n');
                var ri = 0;
                foreach (var r in rows)
                {
                    var cells = ((JsonArray)r!).Select(x => x is null ? "" : ((double?)x ?? 0).ToString("R", inv)).ToList();
                    if (labels != null) cells.Insert(1, (string?)labels[ri] ?? "");
                    sb.Append(string.Join(",", cells.Select(Q))).Append('\n');
                    ++ri; ++n;
                }
            }
            else
            {
                var j = JsonNode.Parse(_inspectorTab == 0 ? _doc.PipelineParticles(_inspectorFilter, 0, int.MaxValue / 2) : _doc.PipelineBonds(0, int.MaxValue / 2))!;
                sb.Append(string.Join(",", ((JsonArray)j["columns"]!).Select(c => Q((string?)c ?? "")))).Append('\n');
                foreach (var r in (JsonArray)j["rows"]!)
                {
                    sb.Append(string.Join(",", ((JsonArray)r!["cells"]!).Select(x => Q((string?)x ?? "")))).Append('\n');
                    ++n;
                }
            }
            System.IO.File.WriteAllText(path, sb.ToString());
            Status = $"Wrote {n:N0} rows to {System.IO.Path.GetFileName(path)}";
        }
        catch (Exception e) { Status = "Could not export: " + e.Message; }
    }

    /// <summary>The data inspector's current tab: particles (filtered), bonds, attributes or a data table.</summary>
    public void LoadInspector()
    {
        Raise(nameof(InspectorShowsTable)); Raise(nameof(InspectorShowsAttributes)); Raise(nameof(InspectorShowsTables)); Raise(nameof(IsParticlesTab));
        if (_doc == null || !IsVisualize) return;
        if (_inspectorTab == 3) { LoadPipeTable(); return; }
        if (_inspectorTab == 2) { InspectorNote = $"{PipeAttributes.Count} attributes of frame {_frame}"; return; }
        string text;
        try { text = _inspectorTab == 0 ? _doc.PipelineParticles(_inspectorFilter, _inspectorOffset, InspectorPage) : _doc.PipelineBonds(_inspectorOffset, InspectorPage); }
        catch (Exception e) { InspectorNote = "Filter: " + e.Message; InspectorRows.Clear(); return; }
        var j = JsonNode.Parse(text)!;
        InspectorColumns.Clear();
        foreach (var c in (JsonArray)j["columns"]!) InspectorColumns.Add((string?)c ?? "");
        InspectorRows.Clear();
        var picked = Picked;
        foreach (var r in (JsonArray)j["rows"]!)
        {
            var cells = ((JsonArray)r!["cells"]!).Select(x => (string?)x ?? "").ToArray();
            var origin = (int?)r["origin"] ?? -1;
            InspectorRows.Add(new TableRow(cells, _inspectorTab == 0 && (origin == picked || ((bool?)r["selected"] ?? false))));
        }
        var total = (int)((double?)j["total"] ?? 0);
        var shownTo = Math.Min(total, _inspectorOffset + InspectorPage);
        InspectorNote = total == 0 ? "No rows" + (_inspectorFilter.Length > 0 ? " match the filter" : "")
            : $"Rows {_inspectorOffset + 1}–{shownTo} of {total}" + (_inspectorTab == 0 && _inspectorFilter.Length > 0 ? $" matching the filter ({(int)((double?)j["in_frame"] ?? 0)} in frame)" : "");
        Raise(nameof(InspectorCanPrev)); Raise(nameof(InspectorCanNext));
        _inspectorTotal = total;
    }
    private int _inspectorTotal;
    public bool InspectorCanPrev => _inspectorOffset > 0;
    public bool InspectorCanNext => _inspectorOffset + InspectorPage < _inspectorTotal;
    public void InspectorPageStep(int d)
    {
        _inspectorOffset = Math.Max(0, _inspectorOffset + d * InspectorPage);
        LoadInspector();
    }

    // ---- Expression selection (design/boards/ExpressionSelect): saved expressions with their counts on this frame
    public ObservableCollection<ExpressionRow> SavedExpressions { get; } = new();
    public bool PipeIsExpression => _pipeSel?.Type == "select_expression";
    public bool PipeCanMakeReal => _pipeSel?.Type is "replicate" or "delete_selected" or "slice";
    private string _exprMeaning = "";
    public string ExprMeaning { get => _exprMeaning; private set => Set(ref _exprMeaning, value); }
    public void RefreshExpressions()
    {
        Raise(nameof(PipeIsExpression)); Raise(nameof(PipeCanMakeReal));
        SavedExpressions.Clear();
        if (_doc == null || !PipeIsExpression) return;
        JsonObject? types = null;
        foreach (var e in _settings.PipelineExpressions)
        {
            try
            {
                var j = JsonNode.Parse(_doc.ExpressionCount(e))!;
                types ??= j["types"] as JsonObject;
                SavedExpressions.Add(new(e, (bool?)j["ok"] == true ? ((double?)j["count"] ?? 0).ToString("0", CultureInfo.InvariantCulture) : "error"));
            }
            catch { SavedExpressions.Add(new(e, "error")); }
        }
        // what the types in the current expression are (Type == 2 → ca)
        var cur = _pipeSel?.Params["expression"]?.ToString() ?? "";
        var m = System.Text.RegularExpressions.Regex.Matches(cur, @"Type\s*==\s*(\d+)");
        ExprMeaning = types == null || m.Count == 0 ? "" : string.Join(" · ", m.Select(x => $"Type {x.Groups[1].Value} means {(string?)types[x.Groups[1].Value] ?? "an unlabelled type"}").Distinct());
    }
    public void UseExpression(string e)
    {
        if (_pipeSel == null) return;
        _pipeSel.Params["expression"] = e;
        BuildStepFields();
        ApplyPipeline();
        RefreshExpressions();
    }
    public void SaveExpression()
    {
        var cur = _pipeSel?.Params["expression"]?.ToString() ?? "";
        if (cur.Length == 0 || _settings.PipelineExpressions.Contains(cur)) return;
        _settings.PipelineExpressions.Add(cur);
        Changed("Saved expressions");
        RefreshExpressions();
    }
    /// <summary>The pipeline's particles (the replicas, what is left after deleting) as a new document.</summary>
    public void MakeReal()
    {
        if (_doc == null) return;
        try
        {
            var d = _doc.MaterializePipeline(Title.Replace(" (unsaved)", "") + " (made real, unsaved)");
            var n = d.Summary().Atoms;
            ClearPipeline();
            Show(d, Title.Replace(" (unsaved)", "") + " (made real, unsaved)");
            GrownUnsaved = true;
            Status = $"Made real: {n:N0} atoms with unique identifiers · the original file is unchanged";
        }
        catch (Exception e) { Status = "Make real: " + e.Message; }
    }

    /// <summary>A box dragged on a plot of the data inspector (a scatter's two properties, or one property against the
    /// histogram's bins) selects what lies in it: a "Brush" expression step at the top of the pipeline.</summary>
    public void ApplyBrush(double x0, double x1, double y0, double y1)
    {
        var brush = PipelineRows.FirstOrDefault(r => r.Type == "select_expression" && (bool?)r.Params["brush"] == true);
        if (double.IsNaN(x0)) { if (brush != null) { PipelineRows.Remove(brush); ApplyPipeline(); Status = "Brush cleared"; } return; }
        var xn = PipeTableXLabel;
        var yn = PipeTableYLabel;
        if (!PipeTableScatter || xn.Length == 0 || yn.Length == 0 || xn.Contains(' ') || yn.Contains(' ')) { Status = "Brush: drag on a scatter plot of two properties"; return; }
        var inv = CultureInfo.InvariantCulture;
        var expr = string.Format(inv, "{0} >= {1:G6} && {0} <= {2:G6} && {3} >= {4:G6} && {3} <= {5:G6}", xn, x0, x1, yn, y0, y1);
        if (brush == null)
        {
            brush = new PipelineRow { Type = "select_expression", Title = "Brush", Icon = "filter", Params = new JsonObject { ["expression"] = expr, ["brush"] = true } };
            WireRow(brush);
            PipelineRows.Insert(0, brush);
        }
        else brush.Params["expression"] = expr;
        ApplyPipeline();
        Status = string.Format(inv, "Brushed: {0} from {1:G3} to {2:G3}, {3} from {4:G3} to {5:G3} · the viewport shows them selected", xn, x0, x1, yn, y0, y1);
    }
}

public sealed record ExpressionRow(string Expression, string Count);
