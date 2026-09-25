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
    public bool IsText => Kind is "text" or "expression" or "number" or "vector";
    public bool IsBool => Kind == "bool";
    public bool IsChoice => Kind == "choice";
    public bool IsMono => Kind is "expression" or "number" or "vector";
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
    public string Level { get => _level; set { _level = value; Raise(nameof(Level)); Raise(nameof(SummaryBrush)); } }
    public IBrush SummaryBrush => Tokens.Brush(_level switch { "error" => "ErrB", "warning" => "WarnB", _ => "DimB" });
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
        new("select_expression", "Expression selection", "Type == 2 && Position.Z > 13", "Select", "filter"),
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
        new("create_bonds", "Create bonds", "from distances or a cutoff", "Visual", "link"),
        new("voids", "Voids & pores", "accessible volume for a probe, voids by size", "Structure", "atom"),
        new("voronoi", "Voronoi volumes", "volume per atom (grid or radical)", "Structure", "hex"),
        new("density_field", "Density field", "smoothed mass density, profile, slice", "Structure", "layers"),
        new("vectors", "Vectors", "end-to-end, dipoles, displacements, velocities", "Visual", "move"),
        new("trajectory_lines", "Trajectory lines", "paths of chain centres or particles", "Visual", "history"),
    ];
    public static readonly string[] StepGroups = ["Colour & style", "Select", "Modify", "Structure", "Measure", "Trajectory", "Visual"];

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
            ShowTableOf(_pipeSel?.Type);
        }
    }

    /// <summary>Selecting a measuring step shows its table in the data inspector.</summary>
    private void ShowTableOf(string? type)
    {
        var name = type switch
        {
            "scatter" => "scatter", "coordination" => "rdf", "cluster" => "clusters", "histogram" => "histogram", "binning" => "binning",
            "molecule_shape" => "molecules", "topology" => "bonds", "voids" => "voids", "voronoi" => "voronoi", "density_field" => "density_profile",
            "msd" => "msd", "vectors" => "vectors", _ => null,
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
            var idx = PipeTables.Count - 1;
            Avalonia.Threading.Dispatcher.UIThread.Post(() => { PipeTable = idx; InspectorTab = 3; });
            Status = $"Time series over {((JsonArray?)_series?["rows"])?.Count ?? 0} frames" + (stride > 1 ? $" (every {stride}th)" : "");
        }
        catch (Exception e) { Status = "Time series failed: " + e.Message; }
        finally { SeriesRunning = false; }
    }
    public bool HasPipeTables => PipeTables.Count > 0;
    public double[] PipeTableX { get; private set; } = [];
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
        row.Toggled = ApplyPipeline;
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

    public void LoadPipeline(string path)
    {
        try
        {
            var text = System.IO.File.ReadAllText(path);
            if (text.TrimStart().StartsWith("caps_pipeline")) text = CapsDocument.PipelineFromYaml(text);   // YAML (SavePipeline)
            var j = JsonNode.Parse(text)!;
            var steps = j is JsonArray a ? a : (JsonArray)j["steps"]!;
            PipelineRows.Clear();
            foreach (var st in steps)
            {
                if (st is not JsonObject o) continue;
                var type = (string?)o["type"] ?? "";
                var kind = StepLibrary.FirstOrDefault(k => k.Type == type);
                var prm = (JsonObject)JsonNode.Parse(o.ToJsonString())!;
                var enabled = (bool?)prm["enabled"] ?? true;
                prm.Remove("type");
                prm.Remove("enabled");
                var row = new PipelineRow { Type = type, Title = kind?.Title ?? type, Icon = kind?.Icon ?? "sliders", Params = prm };
                row.Enabled = enabled;
                row.Toggled = ApplyPipeline;
                PipelineRows.Add(row);
            }
            PipeSelected = PipelineRows.FirstOrDefault();
            _showTableAfterApply = true;
            ApplyPipeline();
            Status = $"Loaded {PipelineRows.Count} steps from {path}";
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
        "cluster" => new JsonObject { ["mode"] = "bonds", ["cutoff"] = 3.2, ["sort_by_size"] = true, ["only_selected"] = false, ["colour"] = false },
        "coordination" => new JsonObject { ["cutoff"] = 3.2, ["bins"] = 200, ["element_a"] = 0, ["element_b"] = 0, ["only_selected"] = false },
        "compute_property" => new JsonObject { ["name"] = "Custom", ["expression"] = "Position.Z", ["only_selected"] = false },
        "replicate" => new JsonObject { ["nx"] = 2, ["ny"] = 2, ["nz"] = 1, ["adjust_cell"] = true },
        "histogram" => new JsonObject { ["property"] = "Charge", ["bins"] = 40, ["only_selected"] = false },
        "binning" => new JsonObject { ["property"] = "Mass", ["axis"] = 2, ["bins"] = 50, ["reduction"] = "density" },
        "topology" => new JsonObject { ["bins"] = 60 },
        "displacements" => new JsonObject { ["reference"] = "first", ["frame"] = 0 },
        "smooth" => new JsonObject { ["window"] = 5 },
        "vectors" => new JsonObject { ["property"] = "end_to_end", ["scale"] = 1.0, ["radius"] = 0.3 },
        "msd" => new JsonObject { ["heavy_only"] = true, ["every"] = 1, ["timestep_fs"] = 1.0 },
        "scatter" => new JsonObject { ["x"] = "DistanceToCOM", ["y"] = "Charge", ["only_selected"] = false },
        "voids" => new JsonObject { ["probe"] = 1.4, ["grid"] = 0.5, ["show"] = true },
        "voronoi" => new JsonObject { ["method"] = "grid", ["grid"] = 0.5 },
        "density_field" => new JsonObject { ["grid"] = 0.8, ["sigma"] = 1.5, ["axis"] = 2, ["position"] = 0.5 },
        "trajectory_lines" => new JsonObject { ["particles"] = "centres", ["from"] = 0, ["radius"] = 0.12 },
        "create_bonds" => new JsonObject { ["mode"] = "perceive", ["tolerance"] = 0.45, ["cutoff"] = 1.6, ["replace"] = false, ["only_selected"] = false },
        _ => new JsonObject(),
    };

    /// <summary>The fields of the selected step's editor, from its parameters.</summary>
    private void BuildStepFields()
    {
        StepFields.Clear();
        if (_pipeSel == null) return;
        var p = _pipeSel.Params;
        string S(string k, string d = "") => p[k] is JsonValue v ? (v.TryGetValue<string>(out var s) ? s : v.ToJsonString()) : d;
        bool B(string k) => p[k] is JsonValue v && v.TryGetValue<bool>(out var b) && b;
        var props = PipeProperties();
        void Add(StepField f)
        {
            f.Changed = () => WriteField(f);
            StepFields.Add(f);
        }
        void Text(string key, string label, string kind = "text", string hint = "") => Add(new StepField { Key = key, Label = label, Kind = kind, Hint = hint, Text = S(key) });
        void Bool(string key, string label) => Add(new StepField { Key = key, Label = label, Kind = "bool", On = B(key) });
        void Choice(string key, string label, string[] choices)
        {
            var value = S(key, choices[0]);
            if (!choices.Contains(value)) choices = [.. choices, value];   // a property made by a step, before its first run
            Add(new StepField { Key = key, Label = label, Kind = "choice", Choices = choices, Text = value });
        }
        switch (_pipeSel.Type)
        {
            case "select_expression": Text("expression", "Expression", "expression", "Type == 2 && Position.Z > 13 · Element == \"O\""); break;
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
                Bool("heavy_only", "Cutoff between heavy atoms"); Bool("sort_by_size", "Sort by size"); Bool("colour", "Colour clusters"); Bool("only_selected", "Only selected"); break;
            case "coordination":
                Text("cutoff", "Cutoff (Å)", "number"); Text("bins", "Bins", "number"); Text("element_a", "Central element (0: any)", "number"); Text("element_b", "Neighbour element (0: any)", "number");
                Bool("inter_only", "Only different molecules"); Bool("only_selected", "Only selected"); break;
            case "topology": Text("bins", "Bins", "number"); break;
            case "displacements": Choice("reference", "Reference", ["first", "previous", "frame"]); Text("frame", "Reference frame", "number"); break;
            case "smooth": Text("window", "Window (frames, centred)", "number"); break;
            case "msd":
                Bool("heavy_only", "Heavy atoms only"); Text("every", "Every n-th atom", "number"); Text("max_lag", "Longest lag (frames)", "number", "blank: half the frames");
                Text("timestep_fs", "Timestep (fs) for D in cm²/s", "number"); break;
            case "scatter": Choice("x", "x", props); Choice("y", "y", props); Bool("only_selected", "Only selected"); break;
            case "voids": Text("probe", "Probe radius (Å)", "number"); Text("grid", "Grid (Å)", "number"); Bool("show", "Show void points, coloured by void"); break;
            case "voronoi": Choice("method", "Method", ["grid", "radical"]); Text("grid", "Grid (Å)", "number"); break;
            case "density_field":
                Text("grid", "Grid (Å)", "number"); Text("sigma", "Smoothing σ (Å)", "number"); Choice("axis", "Slice normal (0 x · 1 y · 2 z)", ["0", "1", "2"]);
                Text("position", "Slice position (0–1 of the cell)", "number"); break;
            case "vectors":
                Choice("property", "Vector", ["end_to_end", "dipole", "displacement", "velocity"]); Text("scale", "Scale (dipole, displacement, velocity)", "number");
                Text("radius", "Arrow radius (Å)", "number"); break;
            case "trajectory_lines":
                Choice("particles", "Trace", ["centres", "selected"]); Text("from", "From frame", "number"); Text("to", "To frame", "number", "blank: the last");
                Text("stride", "Every n-th frame", "number", "blank: about 200 steps"); Text("radius", "Line radius (Å)", "number"); break;
            case "create_bonds":
                Choice("mode", "Mode", ["perceive", "cutoff"]); Text("tolerance", "Tolerance over covalent radii (Å)", "number"); Text("cutoff", "Cutoff (Å)", "number");
                Bool("replace", "Replace the bonds"); Bool("only_selected", "Only selected"); break;
            case "compute_property": Text("name", "Output property"); Text("expression", "Expression", "expression", "e.g. sqrt(Position.X^2 + Position.Y^2)"); Bool("only_selected", "Only selected"); break;
            case "replicate": Text("nx", "Images along a", "number"); Text("ny", "Images along b", "number"); Text("nz", "Images along c", "number"); Bool("adjust_cell", "Enlarge the cell"); break;
            case "histogram": Choice("property", "Property", props); Text("bins", "Bins", "number"); Bool("only_selected", "Only selected"); break;
            case "binning": Choice("property", "Property", props); Choice("axis", "Along", ["0", "1", "2"]); Text("bins", "Bins", "number"); Choice("reduction", "Reduction", ["density", "mean", "sum"]); break;
        }
    }

    private string[] PipeProperties()
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
            case "choice" when f.Key == "axis": p[f.Key] = int.Parse(f.Text, CultureInfo.InvariantCulture); break;
            default: p[f.Key] = f.Text; break;
        }
        ApplyPipeline();
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
    }.ToJsonString();

    /// <summary>Sends the steps to the core, which runs them on the shown frame; then the list, legend and inspector refresh.</summary>
    public void ApplyPipeline()
    {
        if (_doc == null) return;
        ++_pipeGen;
        try
        {
            _doc.SetPipeline(PipelineRows.Count == 0 ? "" : PipelineJson());
            PipeError = "";
        }
        catch (Exception e) { PipeError = e.Message; }
        RefreshPipeline();
        RenderRequested?.Invoke();
    }

    /// <summary>Clears the pipeline in the core (leaving Visualize keeps the steps for when it opens again).</summary>
    private void SuspendPipeline()
    {
        try { _doc?.SetPipeline(""); } catch { }
        _pipeHasLegend = false;
        Raise(nameof(ShowPipeLegend));
    }

    private void RefreshPipeline()
    {
        PipeAttributes.Clear();
        _pipeResult = null;
        if (_doc == null) return;
        string text;
        try { text = _doc.PipelineResult(); } catch { text = ""; }
        _pipeResult = text.Length > 0 ? JsonNode.Parse(text) : null;
        var inv = CultureInfo.InvariantCulture;
        if (_pipeResult?["steps"] is JsonArray steps)
            for (int k = 0; k < steps.Count && k < PipelineRows.Count; ++k)
            {
                PipelineRows[k].Summary = (string?)steps[k]?["summary"] ?? "";
                PipelineRows[k].Level = (string?)steps[k]?["level"] ?? "ok";
            }
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
        Raise(nameof(PipeStatus));
        LoadInspector();
    }

    private void LoadPipeTable()
    {
        PipeTableX = []; PipeTableY = [];
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
            if (_inspectorTab == 3)
            {
                InspectorColumns.Clear();
                foreach (var c in cols) InspectorColumns.Add((string?)c ?? "");
                InspectorRows.Clear();
                var inv = CultureInfo.InvariantCulture;
                foreach (var r in rows.Take(InspectorPage))
                    InspectorRows.Add(new TableRow(((JsonArray)r!).Select(x => x is null ? "–" : ((double?)x ?? 0).ToString("G5", inv)).ToArray(), false));
                InspectorNote = $"{rows.Count} rows · {t["title"]}";
            }
        }
        Raise(nameof(PipeTableXLabel));
        Raise(nameof(PipeTableYLabel));
        PipeTableChanged?.Invoke();
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
}
