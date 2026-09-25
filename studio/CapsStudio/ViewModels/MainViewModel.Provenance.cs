using System.Collections.ObjectModel;
using System.Text.Json.Nodes;
using CapsStudio.Interop;

namespace CapsStudio.ViewModels;

/// <summary>A step of a provenance chain as the timeline shows it.</summary>
public sealed record ProvStepRow(int Number, string Engine, string Summary, string Params, string Detail, bool Differs, bool Last)
{
    public bool HasSummary => Summary.Length > 0;
    public bool HasDetail => Detail.Length > 0;
    public bool HasParams => Params.Length > 0;
}

/// <summary>A key and value row (approximations, inputs, differences).</summary>
public sealed record ProvRow(string Key, string Value, string Other = "")
{
    public bool HasOther => Other.Length > 0;
}

/// <summary>Provenance viewer (design/boards/Provenance): everything that produced the open structure, in order — each
/// step's engine, parameters, random-number generator and seed, and the papers it cites — the approximations in force,
/// the inputs with their sha256, the difference from another run's provenance, BibTeX and the manifest itself.</summary>
public sealed partial class MainViewModel
{
    public bool IsProvenance => _module == 37;
    public ObservableCollection<ProvStepRow> ProvSteps { get; } = new();
    public ObservableCollection<ProvRow> ProvApprox { get; } = new();
    public ObservableCollection<ProvRow> ProvInputs { get; } = new();
    public ObservableCollection<ProvRow> ProvDiff { get; } = new();
    public ObservableCollection<string> ProvDiffNotes { get; } = new();

    private string _provJson = "", _provTitle = "", _provChip = "", _provCompare = "", _provDiffSummary = "", _provFooter = "", _provError = "";
    private bool _provJsonOpen, _provHasCompare, _provDeterministic = true;
    public string ProvTitle { get => _provTitle; private set => Set(ref _provTitle, value); }
    public string ProvStepsChip { get => _provChip; private set => Set(ref _provChip, value); }
    public string ProvCompareText { get => _provCompare; private set => Set(ref _provCompare, value); }
    public bool ProvHasCompare { get => _provHasCompare; private set { if (Set(ref _provHasCompare, value)) Raise(nameof(ProvNoCompare)); } }
    public bool ProvNoCompare => !_provHasCompare;
    public string ProvDiffSummary { get => _provDiffSummary; private set => Set(ref _provDiffSummary, value); }
    public string ProvJson { get => _provJson; private set => Set(ref _provJson, value); }
    public bool ProvJsonOpen { get => _provJsonOpen; set => Set(ref _provJsonOpen, value); }
    public string ProvFooter { get => _provFooter; private set => Set(ref _provFooter, value); }
    public string ProvError { get => _provError; private set { if (Set(ref _provError, value)) Raise(nameof(ProvHasError)); } }
    public bool ProvHasError => _provError.Length > 0;
    public bool ProvEmpty => ProvSteps.Count == 0;
    public bool ProvHasSteps => ProvSteps.Count > 0;
    public bool ProvDiffers => _provHasCompare && ProvDiff.Count + ProvDiffNotes.Count > 0;
    public string ProvGenerator => "CAPS " + (typeof(MainViewModel).Assembly.GetName().Version?.ToString(3) ?? "0.1.0") + " · ABI " + Native.AbiVersion();

    private JsonNode? _provA, _provB;
    private string _provBName = "";

    public void OpenProvenance()
    {
        SetModule(37);
        LoadProvenance();
    }

    private void LoadProvenance()
    {
        ProvError = "";
        _provA = null;
        if (_doc != null)
            try { _provA = JsonNode.Parse(_doc.Provenance()); }
            catch (Exception e) { ProvError = e.Message; }
        ProvTitle = "Provenance · " + (_doc != null ? Path.GetFileName(_doc.Path) : "no structure open");
        ProvJson = _provA?.ToJsonString(new System.Text.Json.JsonSerializerOptions { WriteIndented = true, Encoder = System.Text.Encodings.Web.JavaScriptEncoder.UnsafeRelaxedJsonEscaping }) ?? "";
        _provDeterministic = _provA?["deterministic"]?.GetValue<bool>() ?? true;
        ProvFooter = $"caps-manifest/1.0 · {(_provDeterministic ? "reproducible with the same inputs, seeds and thread count" : "not reproducible")}";
        ProvApprox.Clear();
        if (_provA?["approximations"] is JsonObject ap) foreach (var kv in ap) ProvApprox.Add(new ProvRow(kv.Key, kv.Value?.GetValue<string>() ?? ""));
        if (ProvApprox.Count > 0 && !ProvApprox.Any(r => r.Key == "Estimated parameters")) ProvApprox.Add(new ProvRow("Estimated parameters", "none recorded"));
        ProvInputs.Clear();
        if (_provA?["inputs"] is JsonArray ins)
            foreach (var i in ins)
            {
                var h = i?["sha256"]?.GetValue<string>() ?? "";
                ProvInputs.Add(new ProvRow(i?["name"]?.GetValue<string>() ?? "", h.Length == 64 ? h[..4] + "…" + h[^4..] : h));
            }
        if (_provB != null) Compare(); else FillSteps(new HashSet<int>());
    }

    private void FillSteps(HashSet<int> differing)
    {
        ProvSteps.Clear();
        if (_provA?["steps"] is JsonArray steps)
            for (int k = 0; k < steps.Count; k++)
            {
                var st = steps[k]!;
                var ps = st["params"] as JsonObject;
                var edits = st["engine"]?.GetValue<string>() == "edit.builder";
                var parts = ps == null ? [] : ps.Where(kv => !kv.Key.StartsWith('#')).Select(kv => $"{kv.Key} {kv.Value?.GetValue<string>()}").ToList();
                if (edits && ps != null) parts = ps.Where(kv => kv.Key.StartsWith('#')).TakeLast(3).Select(kv => kv.Value?.GetValue<string>() ?? "").ToList();
                var detail = new List<string>();
                if (st["rng"]?.GetValue<string>() is { Length: > 0 } rng) detail.Add(rng);
                if (st["cites"] is JsonArray c && c.Count > 0) detail.Add("cites " + string.Join(" · ", c.Select(x => x?.GetValue<string>())));
                ProvSteps.Add(new ProvStepRow(k + 1, st["engine"]?.GetValue<string>() ?? "", st["summary"]?.GetValue<string>() ?? "", string.Join(" · ", parts),
                                              string.Join(" · ", detail), differing.Contains(k), k == steps.Count - 1));
            }
        var n = ProvSteps.Count;
        ProvStepsChip = $"{n} step{(n == 1 ? "" : "s")}" + (_provHasCompare ? $" · {differing.Count} differ{(differing.Count == 1 ? "s" : "")}" : "");
        Raise(nameof(ProvEmpty)); Raise(nameof(ProvHasSteps)); Raise(nameof(ProvDiffers));
    }

    /// <summary>Compares with the provenance saved beside another file (run B).</summary>
    public string? CompareProvenanceWith(string path)
    {
        var json = CapsDocument.ProvenanceFile(path);
        var b = JsonNode.Parse(json);
        if (b?["ok"]?.GetValue<bool>() != true) { ProvError = b?["error"]?.GetValue<string>() ?? "no provenance"; return ProvError; }
        _provB = b;
        _provBName = Path.GetFileName(path);
        ProvError = "";
        Compare();
        return null;
    }

    public void ClearProvenanceCompare()
    {
        _provB = null;
        ProvHasCompare = false;
        ProvCompareText = "";
        ProvDiff.Clear();
        ProvDiffNotes.Clear();
        FillSteps(new HashSet<int>());
    }

    private void Compare()
    {
        if (_provA == null || _provB == null) return;
        var d = JsonNode.Parse(CapsDocument.ProvenanceCompare(_provA.ToJsonString(), _provB.ToJsonString()))!;
        ProvHasCompare = true;
        ProvCompareText = "run B · " + _provBName;
        ProvDiff.Clear();
        ProvDiffNotes.Clear();
        var differing = new HashSet<int>();
        if (d["rows"] is JsonArray rows)
            foreach (var r in rows)
            {
                var step = r!["step"]!.GetValue<int>();
                differing.Add(step);
                var key = r["key"]!.GetValue<string>();
                ProvDiff.Add(new ProvRow($"{r["engine"]!.GetValue<string>()} {(key == "rng" ? "seed" : key)}", Seedless(r["a"]!.GetValue<string>(), key), Seedless(r["b"]!.GetValue<string>(), key)));
            }
        if (d["notes"] is JsonArray notes) foreach (var n in notes) ProvDiffNotes.Add(n!.GetValue<string>());
        var sameInputs = d["same_inputs"]?.GetValue<bool>() ?? false;
        var onlySeeds = ProvDiff.Count > 0 && ProvDiff.All(r => r.Key.EndsWith(" seed")) && ProvDiffNotes.Count == 0;
        ProvDiffSummary = ProvDiff.Count == 0 && ProvDiffNotes.Count == 0
            ? "Same inputs, versions, parameters and seeds: the two runs should be identical."
            : onlySeeds && sameInputs
                ? "Same inputs, versions and parameters; only the seed changed, so differences between the two structures are statistical."
                : "The runs differ in the parameters above" + (sameInputs ? "." : ", and their inputs differ.");
        FillSteps(differing);
    }

    private static string Seedless(string v, string key) => key == "rng" && v.Contains("seed ") ? v[(v.LastIndexOf("seed ", StringComparison.Ordinal) + 5)..] : v;

    public string ProvenanceBibtex() => _provA == null ? "" : CapsDocument.ProvenanceBibtex(_provA.ToJsonString());
}
