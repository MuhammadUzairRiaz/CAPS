using System.Globalization;
using System.Text.Json.Nodes;

namespace CapsStudio.ViewModels;

/// <summary>Thermal conductivity › lattice vibrations (core phonons.hpp): the crystal's phonons with its force field —
/// dispersion through the zone, the density of states per element, sound speeds along the reciprocal axes and the
/// harmonic heat capacity — the atoms first relaxed in the supercell's periodic copies.</summary>
public sealed partial class MainViewModel
{
    private bool _phRunning, _phHas;
    private string _phError = "", _phSummary = "", _phSound = "", _phNote = "";
    public bool PhRunning { get => _phRunning; private set { if (Set(ref _phRunning, value)) RaiseBusy(); } }
    public bool PhHas { get => _phHas; private set => Set(ref _phHas, value); }
    public string PhError { get => _phError; private set { if (Set(ref _phError, value)) Raise(nameof(PhHasError)); } }
    public bool PhHasError => _phError.Length > 0;
    public string PhSummary { get => _phSummary; private set => Set(ref _phSummary, value); }
    public string PhSound { get => _phSound; private set => Set(ref _phSound, value); }
    public string PhNote { get => _phNote; private set => Set(ref _phNote, value); }
    /// <summary>The bands (x along the path, one array per band), the labelled points, the DOS and its per-element parts.</summary>
    public (double[] X, double[][] Bands, (double X, string Label)[] Ticks) PhDispersion { get; private set; } = ([], [], []);
    public (double[] W, double[] Dos, (string Label, double[] Y)[] Parts) PhDos { get; private set; } = ([], [], []);
    public event Action? PhononsChanged;

    public async Task RunPhonons()
    {
        if (_doc == null || !Idle || BlockedByField("Phonons")) return;
        var inv = CultureInfo.InvariantCulture;
        var json = new JsonObject { ["cutoff"] = Math.Min(12.0, _relaxCutoff), ["temperatures"] = new JsonArray(100, 298.15, 500) }.ToJsonString();
        var doc = _doc!;
        PhRunning = true;
        PhError = "";
        Status = "Phonons…";
        try
        {
            var res = await Task.Run(() => doc.Phonons(json, (st, f) =>
            {
                Avalonia.Threading.Dispatcher.UIThread.Post(() => { if (PhRunning) Status = $"Phonons: {st} · {f * 100:0} %"; });
                return true;
            }));
            var r = JsonNode.Parse(res)!;
            if (r["ok"]?.GetValue<bool>() != true) { PhError = (string?)r["error"] ?? "phonons failed"; Status = "Phonons: " + PhError; return; }
            double[] A(JsonNode? n) => n is JsonArray a ? a.Select(v => (double?)v ?? 0).ToArray() : [];
            var x = A(r["path_x"]);
            var rows = ((JsonArray)r["bands"]!).Select(A).ToArray();
            var nb = rows.Length > 0 ? rows[0].Length : 0;
            var bands = Enumerable.Range(0, nb).Select(b => rows.Select(row => row[b]).ToArray()).ToArray();
            var tx = A(r["tick_x"]);
            var tl = ((JsonArray)r["tick_labels"]!).Select(v => (string?)v ?? "").ToArray();
            PhDispersion = (x, bands, tx.Zip(tl).Select(p => (p.First, p.Second)).ToArray());
            var labels = ((JsonArray)r["pdos_labels"]!).Select(v => (string?)v ?? "").ToArray();
            var parts = ((JsonArray)r["pdos"]!).Select(A).ToArray();
            PhDos = (A(r["dos_w"]), A(r["dos"]), labels.Zip(parts).Select(p => (p.First, p.Second)).ToArray());
            var snd = ((JsonArray)r["sound"]!).Select(A).ToArray();
            string[] ax = ["a*", "b*", "c*"];
            PhSound = string.Join("\n", snd.Select((v, i) => $"along {ax[i]}: " + string.Join(" · ", v.Select(s => (s / 1000).ToString("0.00", inv))) + " km/s"));
            var T = A(r["T"]); var cv = A(r["Cv"]); var S = A(r["S"]); var zpe = A(r["zpe"]);
            PhSummary = string.Format(inv, "lowest {0:0.0} cm⁻¹ · zero-point {1:0.00} kcal/mol per cell", (double?)r["min_wavenumber"] ?? 0, zpe.FirstOrDefault()) + "\n" +
                        string.Join("\n", T.Select((t, i) => string.Format(inv, "{0:0} K: Cv {1:0.00}, S {2:0.00} cal/mol/K per cell", t, cv[i], S[i])));
            PhNote = string.Join(" · ", ((JsonArray)r["notes"]!).Select(v => (string?)v ?? ""));
            PhHas = true;
            Status = $"Phonons: {nb} branches · lowest {((double?)r["min_wavenumber"] ?? 0).ToString("0.0", inv)} cm⁻¹";
            PhononsChanged?.Invoke();
        }
        catch (Exception e) { PhError = e.Message; }
        finally { PhRunning = false; }
    }
}
