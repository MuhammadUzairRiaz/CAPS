using System.Globalization;
using System.Text.Json.Nodes;
using CapsStudio.Interop;

namespace CapsStudio.ViewModels;

/// <summary>Polymer cell › Mesoscale › density fields (core ddft.hpp): the DPD page's molecules and χ as Gaussian-chain
/// density fields on a grid — dynamic mean-field density functional theory, for morphologies larger and cheaper than beads.
/// The molecules' shares are their counts × lengths.</summary>
public sealed partial class MainViewModel
{
    public static readonly string[] FieldGrids = ["2D · 64 × 64", "2D · 128 × 128", "3D · 32 × 32 × 32", "3D · 48 × 48 × 48"];
    private static readonly int[][] FieldGridSizes = [[64, 64, 1], [128, 128, 1], [32, 32, 32], [48, 48, 48]];
    private int _fieldGrid;
    private decimal _fieldSteps = 3000, _fieldNoise = 100;
    private bool _fieldRunning, _fieldHas;
    private string _fieldNote = "", _fieldError = "";
    public int FieldGrid { get => _fieldGrid; set => Set(ref _fieldGrid, Math.Clamp(value, 0, FieldGrids.Length - 1)); }
    public decimal FieldSteps { get => _fieldSteps; set => Set(ref _fieldSteps, Math.Clamp(Math.Round(value), 10, 1_000_000)); }
    public decimal FieldNoise { get => _fieldNoise; set => Set(ref _fieldNoise, Math.Clamp(value, 0, 100000)); }
    public bool FieldRunning { get => _fieldRunning; private set { if (Set(ref _fieldRunning, value)) RaiseBusy(); } }
    public bool FieldHas { get => _fieldHas; private set => Set(ref _fieldHas, value); }
    public string FieldNote { get => _fieldNote; private set => Set(ref _fieldNote, value); }
    public string FieldError { get => _fieldError; private set { if (Set(ref _fieldError, value)) Raise(nameof(FieldHasError)); } }
    public bool FieldHasError => _fieldError.Length > 0;
    /// <summary>The last frame's mid-plane of the first bead type: (x, y, θ) cells, and the order parameter over τ.</summary>
    public (double[] X, double[] Y, double[] Z, string Type) FieldSlice { get; private set; } = ([], [], [], "");
    public (double X, double Y)[] FieldOrder { get; private set; } = [];
    public event Action? FieldsChanged;

    /// <summary>"A5B5" → 10 beads; "AAB" → 3.</summary>
    public static int SequenceLength(string s)
    {
        int n = 0;
        for (int i = 0; i < s.Length;)
        {
            if (!char.IsLetter(s[i])) return 0;
            int j = i + 1, k = 0;
            while (j < s.Length && char.IsDigit(s[j])) k = k * 10 + (s[j++] - '0');
            n += j > i + 1 ? k : 1;
            i = j;
        }
        return n;
    }

    public async Task RunFields()
    {
        if (!Idle) return;
        var inv = CultureInfo.InvariantCulture;
        var mols = new JsonArray();
        foreach (var sp in DpdSpecies)
        {
            var n = SequenceLength(sp.Sequence.Trim());
            if (n == 0) { FieldError = $"'{sp.Sequence}': a sequence of letters with counts, e.g. A5B5"; return; }
            mols.Add(new JsonObject { ["name"] = sp.Name, ["sequence"] = sp.Sequence.Trim(), ["fraction"] = (double)sp.Count * n });
        }
        var chi = new JsonObject();
        foreach (var c in DpdChis) if (c.Pair.Length == 2) chi[c.Pair] = (double)c.Chi;
        var g = FieldGridSizes[_fieldGrid];
        var json = new JsonObject
        {
            ["molecules"] = mols, ["chi"] = chi, ["grid"] = new JsonArray(g[0], g[1], g[2]), ["kappa"] = 20, ["noise"] = (double)_fieldNoise,
            ["dt"] = 0.01, ["steps"] = (int)_fieldSteps, ["frame_every"] = Math.Max(1, (int)_fieldSteps / 30),
        }.ToJsonString();
        FieldRunning = true;
        FieldError = "";
        Status = "Density fields…";
        try
        {
            var res = await Task.Run(() => CapsDocument.DensityFields(json, f =>
            {
                Avalonia.Threading.Dispatcher.UIThread.Post(() => { if (FieldRunning) Status = $"Density fields · {f * 100:0} %"; });
                return true;
            }));
            var r = JsonNode.Parse(res)!;
            if (r["ok"]?.GetValue<bool>() != true) { FieldError = (string?)r["error"] ?? "density fields failed"; Status = "Density fields: " + FieldError; return; }
            var gr = ((JsonArray)r["grid"]!).Select(v => (int)((double?)v ?? 1)).ToArray();
            var slices = (JsonArray)r["slices"]!;
            var last = ((JsonArray)slices[^1]!)[0] as JsonArray;
            var xs = new List<double>(); var ys = new List<double>(); var zs = new List<double>();
            for (int j = 0; j < gr[1]; ++j)
                for (int i = 0; i < gr[0]; ++i) { xs.Add(i + 0.5); ys.Add(j + 0.5); zs.Add((double?)last![j * gr[0] + i] ?? 0); }
            var types = ((JsonArray)r["types"]!).Select(v => (string?)v ?? "").ToArray();
            FieldSlice = (xs.ToArray(), ys.ToArray(), zs.ToArray(), types.FirstOrDefault() ?? "A");
            var times = ((JsonArray)r["times"]!).Select(v => (double?)v ?? 0).ToArray();
            var order = ((JsonArray)r["order"]!).Select(v => (double?)((JsonArray)v!)[0] ?? 0).ToArray();
            FieldOrder = times.Zip(order).Select(p => (p.First, p.Second)).ToArray();
            var spacing = ((JsonArray)r["spacing"]!).Select(v => (double?)v ?? 0).LastOrDefault();
            FieldNote = string.Format(inv, "spacing of the strongest wave of {0}: {1:0.0} grid spacings · ", types.FirstOrDefault(), spacing) +
                        string.Join(" · ", ((JsonArray)r["notes"]!).Select(v => (string?)v ?? "")) + string.Format(inv, " · {0:0.0} s", (double?)r["seconds"] ?? 0);
            FieldHas = true;
            Status = $"Density fields: τ {times.LastOrDefault().ToString("0.#", inv)}";
            FieldsChanged?.Invoke();
        }
        catch (Exception e) { FieldError = e.Message; }
        finally { FieldRunning = false; }
    }
}
