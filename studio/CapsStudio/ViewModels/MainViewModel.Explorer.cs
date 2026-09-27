using System;
using System.Collections.ObjectModel;
using System.Globalization;
using System.Linq;
using System.Text.Json.Nodes;
using Avalonia.Threading;

namespace CapsStudio.ViewModels;

/// <summary>A row of the Properties explorer: a group header, or a property and its value.</summary>
public sealed record PropRow(string Group, string Name, string Value, bool IsHeader = false)
{
    public bool IsValue => !IsHeader;
}

/// <summary>One output of a job, as the files of a Materials Studio job folder: the result structure, its trajectory,
/// charts, the report, the provenance, a checkpoint.</summary>
public sealed record JobOutput(string Label, string Icon, Action Open);

/// <summary>The Properties explorer and the job-output tree of the project (as Materials Studio's Properties and Project
/// explorers): what the selection is, in numbers; and under each structure, a folder per job with what it produced.</summary>
public partial class MainViewModel
{
    // ---------------------------------------------------------------- Properties explorer
    public ObservableCollection<PropRow> PropertyRows { get; } = new();
    private PropRow[] _allProps = [];
    private string _propFilter = "", _propTitle = "Properties";
    public string PropertyFilter { get => _propFilter; set { if (Set(ref _propFilter, value)) FilterProps(); } }
    public string PropertyTitle { get => _propTitle; private set => Set(ref _propTitle, value); }
    public bool HasProperties => PropertyRows.Count > 0;
    private int _leftTab;   // 0 Properties, 1 Fragments
    public int LeftTab
    {
        get => _leftTab;
        set
        {
            if (!Set(ref _leftTab, Math.Clamp(value, 0, 1))) return;
            Raise(nameof(ShowPropertiesTab)); Raise(nameof(ShowFragmentsTab));
            if (_leftTab == 0 && _propsDirty) RefreshProperties();   // caught up when shown again
        }
    }
    private bool _propsDirty;
    public bool ShowPropertiesTab => _leftTab == 0;
    public bool ShowFragmentsTab => _leftTab == 1;

    private DispatcherTimer? _propTimer;
    /// <summary>Asks for a refresh; many changes in a burst (a pick, the frame, the force field) make one.</summary>
    private void RequestProperties()
    {
        if (_propTimer == null)
        {
            _propTimer = new DispatcherTimer { Interval = TimeSpan.FromMilliseconds(60) };
            _propTimer.Tick += (_, _) => { _propTimer!.Stop(); RefreshProperties(); };
        }
        if (_leftTab != 0) { _propsDirty = true; return; }   // nothing computed while the panel shows Fragments
        _propTimer.Stop();
        _propTimer.Start();
    }

    /// <summary>The picked atom, or the selection, or the structure (with its cell and force field).</summary>
    public void RefreshProperties()
    {
        _propsDirty = false;
        var rows = new System.Collections.Generic.List<PropRow>();
        var inv = CultureInfo.InvariantCulture;
        void H(string g) => rows.Add(new PropRow(g, g, "", true));
        void P(string g, string n, string v) { if (v.Length > 0) rows.Add(new PropRow(g, n, v)); }
        string F(double v, string fmt) => v.ToString(fmt, inv);
        if (_doc == null) { _allProps = []; PropertyTitle = "Properties"; FilterProps(); return; }
        try
        {
            if (Picked >= 0)
            {
                var a = JsonNode.Parse(_doc.AtomProperties(Picked))!;
                if (a["ok"]?.GetValue<bool>() == true)
                {
                    var el = (string?)a["element"] ?? "";
                    PropertyTitle = $"Atom {Picked + 1} · {el}";
                    H("Atom");
                    P("Atom", "Index", (Picked + 1).ToString(inv));
                    P("Atom", "Id in file", F(a["id"]!.GetValue<double>(), "0"));
                    P("Atom", "Element", el);
                    P("Atom", "Name", (string?)a["name"] ?? "");
                    P("Atom", "Force-field type", (string?)a["ff_type"] ?? "not assigned");
                    P("Atom", "Charge", F(a["charge"]!.GetValue<double>(), "+0.0000;−0.0000;0.0000") + " e");
                    P("Atom", "Mass", F(a["mass"]!.GetValue<double>(), "0.000") + " g/mol");
                    P("Atom", "Molecule", F(a["molecule"]!.GetValue<double>(), "0"));
                    if (a["residue"] is { } rs) P("Atom", "Residue", ((string?)a["resname"] ?? "") + " " + F(rs.GetValue<double>(), "0"));
                    H("Position");
                    var x = a["xyz"]!.AsArray().Select(v => v!.GetValue<double>()).ToArray();
                    P("Position", "XYZ (Å)", $"{F(x[0], "0.000")}  {F(x[1], "0.000")}  {F(x[2], "0.000")}");
                    if (a["fractional"] is JsonArray fr)
                    {
                        var f = fr.Select(v => v!.GetValue<double>()).ToArray();
                        P("Position", "Fractional", $"{F(f[0], "0.0000")}  {F(f[1], "0.0000")}  {F(f[2], "0.0000")}");
                    }
                    var nb = a["neighbours"]!.AsArray();
                    H($"Bonded to ({nb.Count})");
                    foreach (var n in nb)
                        P($"Bonded to ({nb.Count})", $"{(string?)n!["element"]} {F(n["index"]!.GetValue<double>() + 1, "0")}",
                          F(n["distance"]!.GetValue<double>(), "0.000") + " Å" + (n["order"]!.GetValue<double>() is var o and > 0 ? o switch { 2 => " · double", 3 => " · triple", 4 => " · aromatic", _ => "" } : ""));
                }
            }
            if (rows.Count == 0)
            {
                var s = JsonNode.Parse(_doc.StructureInfo())!;
                if (s["ok"]?.GetValue<bool>() != true) { _allProps = []; FilterProps(); return; }
                var sel = (int)s["selected"]!.GetValue<double>();
                PropertyTitle = sel > 0 ? $"Structure · {sel:N0} selected" : "Structure";
                H("Structure");
                P("Structure", "Name", Title.Replace(" (unsaved)", ""));
                P("Structure", "Formula", (string?)s["formula"] ?? "");
                P("Structure", "Atoms", F(s["atoms"]!.GetValue<double>(), "N0"));
                P("Structure", "Bonds", F(s["bonds"]!.GetValue<double>(), "N0") + (s["bonds_from_file"]!.GetValue<bool>() ? " · from the file" : " · perceived"));
                P("Structure", "Molecules", F(s["molecules"]!.GetValue<double>(), "N0"));
                P("Structure", "Mass", F(s["mass"]!.GetValue<double>(), "N1") + " g/mol");
                P("Structure", "Net charge", s["has_charges"]!.GetValue<bool>() ? F(s["charge"]!.GetValue<double>(), "+0.000;−0.000;0.000") + " e" : "no charges");
                var frames = (int)s["frames"]!.GetValue<double>();
                if (frames > 1) P("Structure", "Frame", $"{(int)s["frame"]!.GetValue<double>() + 1} of {frames}");
                P("Structure", "Format", (string?)s["format"] ?? "");
                if (s["cell"] is JsonObject c)
                {
                    H("Lattice");
                    P("Lattice", "a  b  c (Å)", $"{F(c["a"]!.GetValue<double>(), "0.000")}  {F(c["b"]!.GetValue<double>(), "0.000")}  {F(c["c"]!.GetValue<double>(), "0.000")}");
                    P("Lattice", "α  β  γ (°)", $"{F(c["alpha"]!.GetValue<double>(), "0.00")}  {F(c["beta"]!.GetValue<double>(), "0.00")}  {F(c["gamma"]!.GetValue<double>(), "0.00")}");
                    P("Lattice", "Volume", F(c["volume"]!.GetValue<double>(), "N1") + " Å³");
                    P("Lattice", "Density", F(c["density"]!.GetValue<double>(), "0.0000") + " g/cm³");
                    var per = c["periodic"]!.AsArray().Select(v => v!.GetValue<bool>()).ToArray();
                    P("Lattice", "Periodic", per.All(p => p) ? "x y z" : string.Join(" ", new[] { "x", "y", "z" }.Where((_, k) => per[k])) is { Length: > 0 } p ? p : "none");
                    P("Lattice", "Positions", s["unwrapped"]!.GetValue<bool>() ? "unwrapped (molecules whole)" : "as read (wrapped)");
                }
                else { H("Lattice"); P("Lattice", "Cell", "none (a molecule in vacuum)"); }
                H("Force field");
                if (s["forcefield"] is JsonObject ff)
                {
                    P("Force field", "Name", (string?)ff["name"] ?? "");
                    P("Force field", "State", ff["complete"]!.GetValue<bool>() ? "complete" : $"incomplete · {F(ff["missing"]!.GetValue<double>(), "0")} missing terms");
                    P("Force field", "Atom types", F(ff["types"]!.GetValue<double>(), "0"));
                }
                else P("Force field", "Name", "not assigned");
                if (s["composition"] is JsonObject comp && comp.Count > 0)
                {
                    H("Composition");
                    var total = s["atoms"]!.GetValue<double>();
                    foreach (var kv in comp)
                        P("Composition", kv.Key, F(kv.Value!.GetValue<double>(), "N0") + $"  ({F(100 * kv.Value!.GetValue<double>() / Math.Max(1, total), "0.0")} %)");
                }
            }
        }
        catch (ObjectDisposedException) { rows.Clear(); }
        catch (Exception e) { rows.Clear(); rows.Add(new PropRow("", "Error", e.Message)); }
        _allProps = rows.ToArray();
        FilterProps();
    }

    private void FilterProps()
    {
        PropertyRows.Clear();
        var q = _propFilter.Trim();
        if (q.Length == 0) { foreach (var r in _allProps) PropertyRows.Add(r); }
        else
        {
            // a header stays when any of its rows match
            foreach (var g in _allProps.GroupBy(r => r.Group))
            {
                var hits = g.Where(r => !r.IsHeader && (r.Name.Contains(q, StringComparison.OrdinalIgnoreCase) || r.Value.Contains(q, StringComparison.OrdinalIgnoreCase))).ToList();
                if (hits.Count == 0) continue;
                if (g.FirstOrDefault(r => r.IsHeader) is { } h) PropertyRows.Add(h);
                foreach (var r in hits) PropertyRows.Add(r);
            }
        }
        Raise(nameof(HasProperties));
    }

    // ---------------------------------------------------------------- job output tree
    /// <summary>A job joins the structure it ran on (Grow and Pack: the cell they made) as a folder of outputs.</summary>
    private void AttachJob(Job job, ProjectItem? item)
    {
        if (item == null || job.Item != null) return;
        job.Item = item;
        item.Jobs.Insert(0, job);
    }

    /// <summary>What the finished job left, as openable items.</summary>
    private void BuildJobOutputs(Job job)
    {
        job.Outputs.Clear();
        var item = job.Item;
        var result = _activeItem;   // a run that keeps the original writes its result to a new structure
        if (job.Status is "done" or "stopped" or "cancelled" && result != null && job.Kind != "Analyze")
        {
            var frames = 1L;
            try { frames = result.Doc.Summary().Frames; } catch { }
            job.Outputs.Add(new JobOutput(result == item ? $"Result · {result.Name}" : $"Result → {result.Name}", "cube", () => { Activate(result); SetModule(8); }));
            if (frames > 1 && job.Kind is "Dynamics" or "Equilibrate" or "Relax" or "React")
                job.Outputs.Add(new JobOutput($"Trajectory · {frames:N0} frames", "play", () => { Activate(result); OpenTrajectory(); }));
            job.Outputs.Add(new JobOutput("Provenance", "history", () => { Activate(result); OpenProvenance(); }));
        }
        if (job.HasCurves)
            job.Outputs.Add(new JobOutput($"Chart · {string.Join(", ", new[] { job.CurveA, job.CurveB }.Where(c => c.Length > 0))}", "chart", () => ShowJob(job)));
        if (job.Kind == "Analyze" && job.IsDone)
            job.Outputs.Add(new JobOutput("Properties and curves", "chart", () => SetModule(1)));
        if (job.Kind == "Dynamics" && MdCanContinue && job.Status is "failed" or "cancelled" or "stopped")
            job.Outputs.Add(new JobOutput(MdContinueLabel.Replace("Continue from", "Checkpoint ·"), "save", () => SetModule(3)));
        job.Outputs.Add(new JobOutput(job.IsFailed ? "Error report" : "Report & log", job.IsFailed ? "alert" : "file", () => ShowJob(job)));
    }

    private void ShowJob(Job job)
    {
        SelectedJob = job;
        SetModule(11);
    }
}
