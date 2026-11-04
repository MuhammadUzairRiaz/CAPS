using System.Collections.ObjectModel;
using System.Globalization;
using System.Text;
using System.Text.Json.Nodes;

namespace CapsStudio.ViewModels;

/// <summary>One parameter the assigned force field lets one fit (caps_fit_terms), ticked to fit it.</summary>
public sealed class FitParamItem
{
    public string Kind { get; init; } = "";
    public string[] Types { get; init; } = [];
    public string Name { get; init; } = "";
    public double Value { get; init; }
    public bool Selected { get; set; }
    public string Label => $"{Kind} {string.Join(" – ", Types)} · {Name} = {Value.ToString("G7", CultureInfo.InvariantCulture)}";
}

/// <summary>Force field › fit to crystal data (core crystalfit.hpp): parameters of the assigned force field refined against the
/// open crystal's observed structure, its lattice energy and elastic constants, by least squares; the fitted values kept in the
/// assignment or written to a force-field file.</summary>
public sealed partial class MainViewModel
{
    private bool _fitRunning, _fitHas, _fitRelaxed = true;
    private string _fitError = "", _fitResult = "", _fitNote = "", _fitTargets = "", _fitUnits = "1";
    private double[] _fitFinal = [];
    private FitParamItem[] _fitUsed = [];
    public ObservableCollection<FitParamItem> FitParams { get; } = [];
    public bool FitRunning { get => _fitRunning; private set { if (Set(ref _fitRunning, value)) RaiseBusy(); } }
    public bool FitHas { get => _fitHas; private set => Set(ref _fitHas, value); }
    public bool FitRelaxed { get => _fitRelaxed; set => Set(ref _fitRelaxed, value); }
    public string FitError { get => _fitError; private set { if (Set(ref _fitError, value)) Raise(nameof(FitHasError)); } }
    public bool FitHasError => _fitError.Length > 0;
    public string FitResult { get => _fitResult; private set => Set(ref _fitResult, value); }
    public string FitNote { get => _fitNote; private set => Set(ref _fitNote, value); }
    /// <summary>One target per line: what, value, optional weight and (energy) unit, as "C11 297" or "energy -941.2 1 eV".</summary>
    public string FitTargets { get => _fitTargets; set => Set(ref _fitTargets, value); }
    public string FitFormulaUnits { get => _fitUnits; set => Set(ref _fitUnits, value); }

    /// <summary>Lists the assigned force field's fittable parameters.</summary>
    public void LoadFitTerms()
    {
        if (_doc == null) return;
        FitError = "";
        try
        {
            var r = JsonNode.Parse(_doc.FitTerms())!;
            FitParams.Clear();
            if (r["ok"]?.GetValue<bool>() != true) { FitError = (string?)r["error"] ?? "no force field"; return; }
            foreach (var t in (JsonArray)r["terms"]!)
            {
                var types = ((JsonArray)t!["types"]!).Select(v => (string?)v ?? "").ToArray();
                var names = ((JsonArray)t["names"]!).Select(v => (string?)v ?? "").ToArray();
                var vals = ((JsonArray)t["values"]!).Select(v => (double?)v ?? 0).ToArray();
                for (var k = 0; k < names.Length; ++k)
                    FitParams.Add(new FitParamItem { Kind = (string?)t["kind"] ?? "", Types = types, Name = names[k], Value = k < vals.Length ? vals[k] : 0 });
            }
            if (FitParams.Count == 0) FitError = "the assigned force field has no pair, bond or charge terms to fit on this structure";
        }
        catch (Exception e) { FitError = e.Message; }
    }

    private JsonArray FitTargetJson()
    {
        var inv = CultureInfo.InvariantCulture;
        var a = new JsonArray();
        foreach (var raw in (_fitTargets ?? "").Split('\n'))
        {
            var w = raw.Split((char[]?)null, StringSplitOptions.RemoveEmptyEntries);
            if (w.Length == 0 || w[0].StartsWith('#')) continue;
            if (w.Length < 2 || !double.TryParse(w[1], NumberStyles.Float, inv, out var v)) throw new FormatException($"target '{raw.Trim()}': what and value, e.g. C11 297");
            var o = new JsonObject { ["kind"] = w[0], ["value"] = v };
            if (w.Length > 2 && double.TryParse(w[2], NumberStyles.Float, inv, out var wt)) o["weight"] = wt;
            var unit = w.Length > 3 ? w[3] : w.Length > 2 && !double.TryParse(w[2], NumberStyles.Float, inv, out _) ? w[2] : "";
            if (unit.Length > 0) o["unit"] = unit;
            a.Add(o);
        }
        return a;
    }

    private static JsonObject FitParamJson(FitParamItem p, double? fixedAt = null)
    {
        var o = new JsonObject { ["kind"] = p.Kind, ["types"] = new JsonArray(p.Types.Select(t => (JsonNode?)JsonValue.Create(t)).ToArray()), ["name"] = p.Name };
        if (fixedAt is double v) { o["lower"] = v; o["upper"] = v; }
        return o;
    }

    private async Task<JsonNode?> FitCall(JsonObject json, string what)
    {
        var doc = _doc!;
        FitRunning = true;
        FitError = "";
        Status = what + "…";
        try
        {
            var text = json.ToJsonString();
            var res = await Task.Run(() => doc.FitCrystal(text, (st, f) =>
            {
                Avalonia.Threading.Dispatcher.UIThread.Post(() => { if (FitRunning) Status = $"{what}: {st} · {f * 100:0} %"; });
                return true;
            }));
            var r = JsonNode.Parse(res)!;
            if (r["ok"]?.GetValue<bool>() != true) { FitError = (string?)r["error"] ?? "fit failed"; Status = what + ": " + FitError; return null; }
            return r;
        }
        catch (Exception e) { FitError = e.Message; return null; }
        finally { FitRunning = false; }
    }

    private JsonObject FitBase()
    {
        var inv = CultureInfo.InvariantCulture;
        if (!double.TryParse(_fitUnits, NumberStyles.Float, inv, out var fu) || fu <= 0) throw new FormatException("formula units: a positive number");
        return new JsonObject { ["relaxed"] = _fitRelaxed, ["formula_units"] = fu, ["cutoff"] = Math.Min(12.0, _relaxCutoff), ["targets"] = FitTargetJson() };
    }

    public async Task RunFit()
    {
        if (_doc == null || !Idle || BlockedByField("Fit")) return;
        JsonObject json;
        FitParamItem[] used;
        try
        {
            used = FitParams.Where(p => p.Selected).ToArray();
            json = FitBase();
            json["parameters"] = new JsonArray(used.Select(p => (JsonNode?)FitParamJson(p)).ToArray());
        }
        catch (Exception e) { FitError = e.Message; return; }
        var r = await FitCall(json, "Fit");
        if (r == null) return;
        var inv = CultureInfo.InvariantCulture;
        var sb = new StringBuilder();
        var ps = (JsonArray)r["parameters"]!;
        if (ps.Count > 0)
        {
            sb.AppendLine($"{"parameter",-36} {"start",14} {"fitted",14} {"± error",11}");
            foreach (var p in ps)
                sb.AppendLine(string.Format(inv, "{0,-36} {1,14:G8} {2,14:G8} {3,11}", (string?)p!["name"], (double?)p["initial"] ?? 0, (double?)p["final"] ?? 0,
                    p["error"] is JsonNode e ? ((double?)e ?? 0).ToString("G3", inv) : "–"));
        }
        var ob = (JsonArray)r["observations"]!;
        if (ob.Count > 0)
        {
            sb.AppendLine($"{"observable",-12} {"target",12} {"start",12} {"fitted",12}");
            foreach (var o in ob)
                sb.AppendLine(string.Format(inv, "{0,-12} {1,12:0.####} {2,12:0.####} {3,12:0.####} {4}", (string?)o!["what"], (double?)o["target"] ?? 0, (double?)o["initial"] ?? 0,
                    (double?)o["final"] ?? 0, (string?)o["unit"]));
        }
        foreach (var st in (JsonArray)r["structures"]!)
        {
            double[] A(JsonNode? n) => n is JsonArray a ? a.Select(v => (double?)v ?? 0).ToArray() : [];
            var f = A(st!["final"]); var obs = A(st["observed"]);
            sb.AppendLine(string.Format(inv, "cell  {0:0.0000} {1:0.0000} {2:0.0000} Å  {3:0.00} {4:0.00} {5:0.00}°  (observed {6:0.0000} {7:0.0000} {8:0.0000}, {9:0.00} {10:0.00} {11:0.00})",
                f[0], f[1], f[2], f[3], f[4], f[5], obs[0], obs[1], obs[2], obs[3], obs[4], obs[5]));
            sb.AppendLine(_fitRelaxed
                ? string.Format(inv, "atoms' RMS shift from observed {0:0.0000} Å · lattice energy {1:0.000} kcal/mol per formula unit", (double?)st["shift_final"] ?? 0, (double?)st["energy_final"] ?? 0)
                : string.Format(inv, "largest force at the observed structure {0:G3} kcal/mol/Å · lattice energy {1:0.000} kcal/mol per formula unit", (double?)st["force_final"] ?? 0, (double?)st["energy_final"] ?? 0));
        }
        sb.Append(string.Format(inv, "S {0:G4} → {1:G4} · {2} iterations ({3}) · {4:0.0} s", (double?)r["ss_initial"] ?? 0, (double?)r["ss_final"] ?? 0,
            (int)((double?)r["iterations"] ?? 0), (string?)r["stopped_by"], (double?)r["seconds"] ?? 0));
        FitResult = sb.ToString();
        FitNote = string.Join(" · ", ((JsonArray)r["notes"]!).Select(v => (string?)v ?? ""));
        _fitUsed = used;
        _fitFinal = ps.Select(p => (double?)p!["final"] ?? 0).ToArray();
        FitHas = true;
        Status = used.Length == 0 ? "Fit: the crystal evaluated with the assigned force field" : $"Fit: {used.Length} parameters, S {((double?)r["ss_final"] ?? 0).ToString("G4", inv)}";
    }

    /// <summary>The last fit's values written into the assignment (apply) or to a force-field file (save): one evaluation with
    /// each parameter held at its fitted value.</summary>
    public async Task KeepFit(string savePath = "")
    {
        if (_doc == null || !Idle || _fitUsed.Length == 0) return;
        JsonObject json;
        try
        {
            json = FitBase();
            json["targets"] = new JsonArray();
            json["iterations"] = 0;
            json["parameters"] = new JsonArray(_fitUsed.Select((p, k) => (JsonNode?)FitParamJson(p, _fitFinal[k])).ToArray());
            if (savePath.Length > 0) json["save"] = savePath;
            else json["apply"] = true;
        }
        catch (Exception e) { FitError = e.Message; return; }
        var r = await FitCall(json, savePath.Length > 0 ? "Save fitted force field" : "Keep fitted values");
        if (r == null) return;
        if (savePath.Length > 0) Status = "Fitted force field written to " + savePath;
        else
        {
            Field.LoadReport(_doc!);
            LoadFitTerms();
            Status = "The fitted values are now the assignment's";
        }
    }
}
