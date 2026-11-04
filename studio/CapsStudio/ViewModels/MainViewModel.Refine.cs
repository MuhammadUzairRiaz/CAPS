using System.Globalization;
using System.Text;
using System.Text.Json.Nodes;

namespace CapsStudio.ViewModels;

/// <summary>Analyze › Scattering › powder: the measured pattern (loaded in 2θ) refined with the open crystal (core refine.hpp):
/// Rietveld or Pawley, the lattice, zero shift and peak widths, the agreement factors and the difference curve.</summary>
public sealed partial class MainViewModel
{
    private bool _refPawley, _refRunning, _refHas;
    private string _refText = "", _refError = "", _expPath = "";
    /// <summary>The loaded measured pattern's file when it is in 2θ ("" otherwise).</summary>
    public string ExperimentPath { get => _expPath; private set { if (Set(ref _expPath, value)) Raise(nameof(CanRefine)); } }
    public bool CanRefine => _expPath.Length > 0;
    public bool RefinePawley { get => _refPawley; set => Set(ref _refPawley, value); }
    public bool RefineRunning { get => _refRunning; private set { if (Set(ref _refRunning, value)) RaiseBusy(); } }
    public bool RefineHas { get => _refHas; private set => Set(ref _refHas, value); }
    public string RefineText { get => _refText; private set => Set(ref _refText, value); }
    public string RefineError { get => _refError; private set { if (Set(ref _refError, value)) Raise(nameof(RefineHasError)); } }
    public bool RefineHasError => _refError.Length > 0;
    private string _idxText = "";
    public string IndexText { get => _idxText; private set { if (Set(ref _idxText, value)) Raise(nameof(IndexHas)); } }
    public bool IndexHas => _idxText.Length > 0;

    /// <summary>The cells that index the measured pattern's lines (successive dichotomy, by M20).</summary>
    public async Task RunIndex()
    {
        if (_expPath.Length == 0 || !Idle) return;
        var lam = new[] { 1.54056, 0.70930, 1.78896, 1.93604, 2.28970, 0.559421 }[Math.Clamp(_powAnode, 0, 5)];
        var json = new JsonObject { ["data"] = _expPath, ["wavelength"] = lam, ["max_length"] = 15 }.ToJsonString();
        RefineRunning = true;
        RefineError = "";
        Status = "Indexing…";
        try
        {
            var r = JsonNode.Parse(await Task.Run(() => CapsStudio.Interop.CapsDocument.Index(json)))!;
            if (r["ok"]?.GetValue<bool>() != true) { RefineError = (string?)r["error"] ?? "indexing failed"; return; }
            var inv = CultureInfo.InvariantCulture;
            var sols = (JsonArray)r["solutions"]!;
            var sb = new StringBuilder($"{((JsonArray)r["lines"]!).Count} lines · ");
            if (sols.Count == 0) sb.Append("no cell indexes them (try a larger tolerance or more unindexed lines from the command line)");
            foreach (var x in sols.Take(5))
            {
                var c = ((JsonArray)x!["cell"]!).Select(v => (double?)v ?? 0).ToArray();
                sb.AppendLine().Append(string.Format(inv, "{0,-12} {1:0.0000} {2:0.0000} {3:0.0000} Å {4:0.00} {5:0.00} {6:0.00}° · M20 {7:0.0}", (string?)x["system"], c[0], c[1], c[2], c[3], c[4], c[5], (double?)x["M20"] ?? 0));
            }
            IndexText = sb.ToString();
            Status = sols.Count > 0 ? $"Indexed: {(string?)sols[0]!["system"]}" : "Indexing: no cell";
        }
        catch (Exception e) { RefineError = e.Message; }
        finally { RefineRunning = false; }
    }

    public ((double X, double Y)[] Obs, (double X, double Y)[] Calc, (double X, double Y)[] Diff) RefineFit { get; private set; } = ([], [], []);

    public async Task RunRefine()
    {
        if (_doc == null || !Idle || _expPath.Length == 0) return;
        var json = new JsonObject
        {
            ["data"] = _expPath, ["mode"] = _refPawley ? "pawley" : "rietveld", ["anode"] = PowderAnodeIds[_powAnode], ["kalpha2"] = _powKa2,
            ["W"] = Math.Max(_powW, 1e-4), ["Y"] = 0.02, ["two_theta_min"] = Math.Min(_powFrom, _powTo - 1), ["two_theta_max"] = _powTo,
            ["refine"] = new JsonArray("U", "W", "Y", "zero"),
        }.ToJsonString();
        var doc = _doc!;
        RefineRunning = true;
        RefineError = "";
        Status = "Refining…";
        try
        {
            var res = await Task.Run(() => doc.Refine(json, (st, f) =>
            {
                Avalonia.Threading.Dispatcher.UIThread.Post(() => { if (RefineRunning) Status = $"Refining: {st}"; });
                return true;
            }));
            var r = JsonNode.Parse(res)!;
            if (r["ok"]?.GetValue<bool>() != true) { RefineError = (string?)r["error"] ?? "the refinement failed"; Status = "Refine: " + RefineError; return; }
            double[] A(JsonNode? n) => n is JsonArray a ? a.Select(v => (double?)v ?? 0).ToArray() : [];
            var t = A(r["two_theta"]); var yo = A(r["y_obs"]); var yc = A(r["y_calc"]);
            var off = yo.Length > 0 ? -0.1 * yo.Max() : 0;
            RefineFit = (t.Zip(yo).ToArray(), t.Zip(yc).ToArray(), t.Select((x, i) => (x, yo[i] - yc[i] + off)).ToArray());
            var inv = CultureInfo.InvariantCulture;
            var sb = new StringBuilder();
            sb.AppendLine(string.Format(inv, "{0} · Rwp {1:0.0000} · Rp {2:0.0000} · Rexp {3:0.0000} · χ² {4:0.###}", _refPawley ? "Pawley" : "Rietveld",
                (double?)r["rwp"] ?? 0, (double?)r["rp"] ?? 0, (double?)r["rexp"] ?? 0, (double?)r["chi2"] ?? 0));
            foreach (var p in (JsonArray)r["phases"]!)
            {
                var l = A(p!["lattice"]); var e = A(p["lattice_esd"]);
                sb.AppendLine(string.Format(inv, "a b c {0:0.00000} {1:0.00000} {2:0.00000} Å (± {3:0.0e0}) · α β γ {4:0.000} {5:0.000} {6:0.000}° · V {7:0.000} Å³", l[0], l[1], l[2], e[0], l[3], l[4], l[5],
                    (double?)p["volume"] ?? 0));
            }
            sb.Append(string.Format(inv, "U {0:G4} W {1:G4} Y {2:G4} · zero {3:0.0000}°", (double?)r["U"] ?? 0, (double?)r["W"] ?? 0, (double?)r["Y"] ?? 0, (double?)r["zero"] ?? 0));
            if (((double?)r["size"] ?? 0) > 0) sb.Append(string.Format(inv, " · crystallites {0:0.#} nm", ((double?)r["size"] ?? 0) / 10));
            RefineText = sb.ToString();
            RefineHas = true;
            Status = string.Format(inv, "Refined: Rwp {0:0.0000}", (double?)r["rwp"] ?? 0);
            ScatterChanged?.Invoke();
        }
        catch (Exception e) { RefineError = e.Message; }
        finally { RefineRunning = false; }
    }
}
