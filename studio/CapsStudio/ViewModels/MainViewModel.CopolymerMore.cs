using System.Globalization;
using System.Text;
using System.Text.Json.Nodes;
using CapsStudio.Interop;

namespace CapsStudio.ViewModels;

/// <summary>Builders › Copolymer, beyond two units in the terminal model (core polystats.hpp): a third monomer with its
/// reactivity ratios, the penultimate model, and the feed's drift with conversion in a batch, each chain made at its
/// own conversion or the conversion running along every chain.</summary>
public sealed partial class MainViewModel
{
    public static string[] CoModels { get; } = ["Terminal (last unit)", "Penultimate (last two units)"];
    public static string[] CoDrifts { get; } = ["No drift (low conversion)", "Drift between chains (batch)", "Drift along each chain (living)"];
    private int _coModel, _coThird, _coDrift;
    private decimal _coR1p = 0.52m, _coR2p = 0.46m, _coF3 = 0.2m, _coR13 = 1, _coR31 = 1, _coR23 = 1, _coR32 = 1, _coConversion = 0.8m;
    private string _coMoreText = "";
    public int CoModel { get => _coModel; set { if (Set(ref _coModel, Math.Clamp(value, 0, 1))) { Raise(nameof(CoPenultimate)); Raise(nameof(CoThirdAllowed)); Raise(nameof(CoHasThird)); Raise(nameof(CoM3Name)); CoRecompute(true); } } }
    public bool CoPenultimate => _coModel == 1;
    public bool CoThirdAllowed => _coModel == 0;
    public decimal CoR1p { get => _coR1p; set { if (Set(ref _coR1p, Math.Clamp(value, 0, 100))) CoRecompute(true); } }
    public decimal CoR2p { get => _coR2p; set { if (Set(ref _coR2p, Math.Clamp(value, 0, 100))) CoRecompute(true); } }
    /// <summary>The third monomer: 0 none, k the monomer k − 1 of the list.</summary>
    public int CoThird { get => _coThird; set { if (Set(ref _coThird, Math.Clamp(value, 0, CoMonomers.Length))) { Raise(nameof(CoHasThird)); Raise(nameof(CoM3Name)); CoRecompute(true); } } }
    public string[] CoThirdNames => CoMonomers.Select(m => m.Name).Prepend("none").ToArray();
    public bool CoHasThird => _coThird > 0 && _coModel == 0;
    public string CoM3Name => CoHasThird ? CoMonomers[_coThird - 1].Name : "";
    public decimal CoF3 { get => _coF3; set { if (Set(ref _coF3, Math.Clamp(value, 0.01m, 0.98m))) CoRecompute(true); } }
    public decimal CoR13 { get => _coR13; set { if (Set(ref _coR13, Math.Clamp(value, 0.001m, 100))) CoRecompute(true); } }
    public decimal CoR31 { get => _coR31; set { if (Set(ref _coR31, Math.Clamp(value, 0.001m, 100))) CoRecompute(true); } }
    public decimal CoR23 { get => _coR23; set { if (Set(ref _coR23, Math.Clamp(value, 0.001m, 100))) CoRecompute(true); } }
    public decimal CoR32 { get => _coR32; set { if (Set(ref _coR32, Math.Clamp(value, 0.001m, 100))) CoRecompute(true); } }
    public int CoDrift { get => _coDrift; set { if (Set(ref _coDrift, Math.Clamp(value, 0, 2))) { Raise(nameof(CoHasDrift)); CoRecompute(true); } } }
    public bool CoHasDrift => _coDrift > 0;
    public decimal CoConversion { get => _coConversion; set { if (Set(ref _coConversion, Math.Clamp(value, 0.01m, 0.99m))) CoRecompute(true); } }
    public string CoMoreText { get => _coMoreText; private set { if (Set(ref _coMoreText, value)) Raise(nameof(CoHasMoreText)); } }
    public bool CoHasMoreText => _coMoreText.Length > 0;
    private bool CoExtended => CoHasThird || _coModel == 1 || _coDrift > 0;
    private string[] CoUnitShorts => CoHasThird ? [CoMonomers[_coM1].Short, CoMonomers[_coM2].Short, CoMonomers[_coThird - 1].Short] : [CoMonomers[_coM1].Short, CoMonomers[_coM2].Short];

    private double[] CoFeed()
    {
        var f1 = (double)_coF1;
        if (!CoHasThird) return [f1, 1 - f1];
        var f3 = (double)_coF3;
        return [f1 * (1 - f3), (1 - f1) * (1 - f3), f3];
    }

    /// <summary>The extra settings into a chain spec (the third unit, the model, the drift).</summary>
    private void CoSpecMore(JsonObject o)
    {
        if (CoHasThird)
        {
            ((JsonArray)o["units"]!).Add(new JsonObject { ["name"] = CoMonomers[_coThird - 1].Short, ["smiles"] = CoMonomers[_coThird - 1].Smiles });
            o["weights"] = new JsonArray(CoFeed().Select(x => (JsonNode)x).ToArray());
            o["reactivity"] = new JsonArray(1.0, (double)_coR1, (double)_coR13, (double)_coR2, 1.0, (double)_coR23, (double)_coR31, (double)_coR32, 1.0);
        }
        if (_coModel == 1)
        {
            o["sequence"] = "penultimate";
            o["penultimate"] = new JsonArray((double)_coR1, (double)_coR1p, (double)_coR2, (double)_coR2p);
        }
        if (_coDrift > 0)
        {
            o["drift"] = _coDrift == 1 ? "chains" : "along";
            o["conversion"] = (double)_coConversion;
        }
    }

    /// <summary>The model of n units against the chain drawn: composition and run lengths per unit, the drift's end.</summary>
    private void CoRecomputeMore()
    {
        var inv = CultureInfo.InvariantCulture;
        var q = new JsonObject { ["feed"] = new JsonArray(CoFeed().Select(x => (JsonNode)x).ToArray()) };
        var spec = CoSpec();
        foreach (var k in new[] { "reactivity", "penultimate" })
            if (spec[k] is JsonNode node) q[k] = node.DeepClone();
        if (!CoHasThird) { q["r1"] = (double)_coR1; q["r2"] = (double)_coR2; }
        if (_coDrift > 0) { q["conversion"] = (double)_coConversion; q["points"] = 21; }
        var r = JsonNode.Parse(CapsDocument.CopolymerModel(q.ToJsonString()))!;
        if ((bool?)r["ok"] != true) { CoError = (string?)r["error"] ?? "cannot compute"; return; }
        // the chain Grow draws for chain 0 at this seed, with the full spec
        var pv = JsonNode.Parse(CapsDocument.ChainPreview(spec.ToJsonString(), (ulong)_coSeed))!;
        CoSequence = (pv["sequence"] as JsonArray ?? []).Select(x => (int)((double?)x ?? 0)).ToArray();
        var names = CoUnitShorts;
        var n = names.Length;
        double[] A(JsonNode? a) => a is JsonArray arr ? arr.Select(x => (double?)x ?? double.NaN).ToArray() : [];
        var comp = A(r["model"]?["composition"]);
        var runs = A(r["model"]?["run_length"]);
        // this chain's own numbers
        var cnt = new double[n];
        var runsN = new double[n];
        for (var i = 0; i < CoSequence.Length; ++i)
        {
            var u = Math.Clamp(CoSequence[i], 0, n - 1);
            cnt[u] += 1;
            if (i == 0 || CoSequence[i] != CoSequence[i - 1]) runsN[u] += 1;
        }
        CoRows.Clear();
        var drift = r["drift"] as JsonArray;
        var cum = drift is { Count: > 0 } ? A(drift[^1]!["cumulative"]) : comp;
        for (var k = 0; k < n; ++k)
            CoRows.Add(new($"Fraction {names[k]}" + (drift != null ? " (whole batch)" : ""), k < cum.Length ? cum[k].ToString("0.000", inv) : "—",
                           CoSequence.Length > 0 ? (cnt[k] / CoSequence.Length).ToString("0.000", inv) : "—"));
        if (drift == null)
            for (var k = 0; k < n; ++k)
                CoRows.Add(new($"Mean run, {names[k]}", k < runs.Length && double.IsFinite(runs[k]) ? runs[k].ToString("0.00", inv) : "—",
                               runsN[k] > 0 ? (cnt[k] / runsN[k]).ToString("0.00", inv) : "—"));
        var sb = new StringBuilder();
        sb.Append(_coModel == 1
            ? "Penultimate model (Merz, Alfrey & Goldfinger, J. Polym. Sci. 1, 75 (1946)): the last two units decide the next; r₁′ and r₂′ apply after the other unit."
            : CoHasThird
                ? "Terminal model of three units (Alfrey & Goldfinger, J. Chem. Phys. 12, 322 (1944))."
                : "Terminal model.");
        if (drift is { Count: > 1 })
        {
            var f0 = A(drift[0]!["feed"]);
            var fe = A(drift[^1]!["feed"]);
            var ie = A(drift[^1]!["instantaneous"]);
            sb.Append($" Drift (Skeist, J. Am. Chem. Soc. 68, 1781 (1946)): by {(double)_coConversion:0%} conversion the feed goes from ");
            sb.Append(string.Join(" / ", f0.Select(x => x.ToString("0.00", inv))) + " to " + string.Join(" / ", fe.Select(x => x.ToString("0.00", inv))));
            sb.Append(", the copolymer formed last " + string.Join(" / ", ie.Select(x => x.ToString("0.00", inv))) + " (" + string.Join(" / ", names) + ").");
            sb.Append(_coDrift == 1 ? " Each chain is made at a conversion drawn evenly up to that." : " Each chain runs through the whole conversion, head to tail.");
        }
        CoMoreText = sb.ToString();
        CoNote = CoMoreText;
        CoStatus = $"{CoSequence.Length} units · " + string.Join(" · ", names.Select((nm, k) => $"{cnt[k]:0} {nm}"));
        CoFooter = $"sequence from the {(_coModel == 1 ? "penultimate" : "terminal")} model{(_coDrift > 0 ? " with drift" : "")}, seed {_coSeed}";
    }

    /// <summary>A grown copolymer cell's statistics for the Grow log: composition and runs, model against the cell.</summary>
    private static string SequenceSummary(CapsDocument doc)
    {
        try
        {
            var r = JsonNode.Parse(doc.SequenceReport());
            if (r?["built"] is not JsonObject b) return "";
            var inv = CultureInfo.InvariantCulture;
            string L(JsonNode? a, string f) => a is JsonArray arr ? string.Join(" ", arr.Select(x => ((double?)x ?? double.NaN).ToString(f, inv))) : "—";
            var sb = new StringBuilder("\n\nsequence of the cell (" + ((string?)r["kind"] ?? "") + ", units A B …):");
            if (r["model"] is JsonObject m) sb.Append("\n  model composition ").Append(L(m["composition"], "0.000")).Append(" · runs ").Append(L(m["run_length"], "0.00"));
            if (r["drift"] is JsonArray d && d.Count > 0) sb.Append("\n  model, whole batch ").Append(L(d[^1]!["cumulative"], "0.000"));
            sb.Append("\n  cell  composition ").Append(L(b["composition"], "0.000")).Append(" · runs ").Append(L(b["run_length"], "0.00"));
            sb.Append($" · {(double?)b["count"] ?? 0:0} units");
            return sb.ToString();
        }
        catch { return ""; }
    }
}
