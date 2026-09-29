using System.Collections.ObjectModel;
using System.Globalization;
using System.Text.Json.Nodes;

namespace CapsStudio.ViewModels;

/// <summary>An adsorbate to place: its SMILES and how many molecules.</summary>
public sealed class AdsorbateRow : ObservableObject
{
    private string _smiles = "";
    private decimal _count = 1;
    public string Smiles { get => _smiles; set => Set(ref _smiles, value ?? ""); }
    public decimal Count { get => _count; set => Set(ref _count, Math.Clamp(value, 1, 500)); }
}

public sealed record AdsComponentRow(string Name, string Molecules, string DeDn);
public sealed record AdsConfigRow(string Rank, string Energy, string Cycle);

/// <summary>Analyze › Adsorption locator: where adsorbate molecules sit on a substrate, by Monte Carlo simulated
/// annealing with the assigned force field (core adsorption.hpp; the method of Materials Studio's Adsorption Locator).
/// The open structure is the substrate; the adsorbates are built from SMILES and added after it.</summary>
public sealed partial class MainViewModel
{
    public bool IsAdsorption => _module == 70;
    private bool _adsRunning;
    private int _adsRegion = 1, _adsCycles = 3, _adsSteps = 20000;
    private double _adsTHigh = 1e4, _adsTLow = 100;
    private string _adsNote = "", _adsEnergy = "—", _adsSplit = "", _adsError = "", _adsProgress = "";
    private CancellationTokenSource? _adsCts;
    public ObservableCollection<AdsorbateRow> AdsRows { get; } = new() { new AdsorbateRow { Smiles = "O=C=O", Count = 4 } };
    public ObservableCollection<AdsComponentRow> AdsComponents { get; } = new();
    public ObservableCollection<AdsConfigRow> AdsConfigs { get; } = new();
    public (double X, double Y)[] AdsHistogram { get; private set; } = [];
    public static readonly string[] AdsRegions = ["Anywhere in the cell", "Above the substrate's top face"];
    public int AdsRegion { get => _adsRegion; set => Set(ref _adsRegion, Math.Clamp(value, 0, 1)); }
    public decimal AdsCyclesD { get => _adsCycles; set => Set(ref _adsCycles, (int)Math.Clamp(value, 1, 100)); }
    public decimal AdsStepsD { get => _adsSteps; set => Set(ref _adsSteps, (int)Math.Clamp(value, 100, 10_000_000)); }
    public decimal AdsTHighD { get => (decimal)_adsTHigh; set => Set(ref _adsTHigh, Math.Clamp((double)value, 1, 1e7)); }
    public decimal AdsTLowD { get => (decimal)_adsTLow; set => Set(ref _adsTLow, Math.Clamp((double)value, 1, 1e6)); }
    public bool AdsRunning { get => _adsRunning; private set { if (Set(ref _adsRunning, value)) { RaiseBusy(); Raise(nameof(AdsIdle)); } } }
    public bool AdsIdle => !_adsRunning;
    public string AdsNote { get => _adsNote; private set => Set(ref _adsNote, value); }
    public string AdsEnergy { get => _adsEnergy; private set => Set(ref _adsEnergy, value); }
    public string AdsSplit { get => _adsSplit; private set => Set(ref _adsSplit, value); }
    public string AdsError { get => _adsError; private set { if (Set(ref _adsError, value)) Raise(nameof(AdsHasError)); } }
    public bool AdsHasError => _adsError.Length > 0;
    public string AdsProgress { get => _adsProgress; private set => Set(ref _adsProgress, value); }
    public event Action? AdsorptionChanged;

    public void OpenAdsorption() => SetModule(70);
    public void AddAdsorbate() => AdsRows.Add(new AdsorbateRow { Smiles = "O", Count = 1 });
    public void RemoveAdsorbate(AdsorbateRow r) { if (AdsRows.Count > 1) AdsRows.Remove(r); }
    public void CancelAdsorption() => _adsCts?.Cancel();

    /// <summary>Adds the adsorbates to the open structure and anneals them on it; the result is new frames (undoable).</summary>
    public async Task RunAdsorption()
    {
        if (_doc == null || !Idle || BlockedByField("Adsorption locator")) return;
        var rows = AdsRows.Where(r => r.Smiles.Trim().Length > 0).ToArray();
        if (rows.Length == 0) { AdsError = "Give at least one adsorbate (SMILES)"; return; }
        PrepareRunTarget("adsorbed");
        var doc = _doc!;
        var atomsBefore = doc.Summary().Atoms;
        var inv = CultureInfo.InvariantCulture;
        var json = new JsonObject
        {
            ["adsorbates"] = new JsonArray(rows.Select(r => (JsonNode)new JsonObject { ["smiles"] = r.Smiles.Trim(), ["count"] = (int)r.Count }).ToArray()),
            ["cycles"] = _adsCycles, ["steps"] = _adsSteps, ["t_high"] = _adsTHigh, ["t_low"] = _adsTLow,
            ["region"] = _adsRegion == 1 ? "above" : "cell", ["cutoff"] = Math.Min(12.0, _relaxCutoff), ["coulomb"] = _relaxCoulomb, ["seed"] = AdsSeed.Take(),
        }.ToJsonString();
        AdsRunning = true;
        AdsError = "";
        _adsCts = new CancellationTokenSource();
        var ct = _adsCts.Token;
        Status = "Locating adsorption sites…";
        try
        {
            var res = await Task.Run(() => doc.Adsorption(json, (st, f) =>
            {
                Avalonia.Threading.Dispatcher.UIThread.Post(() => AdsProgress = $"{st} · {f * 100:0} %");
                return !ct.IsCancellationRequested;
            }));
            var r = JsonNode.Parse(res)!;
            if (r["ok"]?.GetValue<bool>() != true)
            {
                AdsError = (string?)r["error"] ?? "the locator failed";
                Status = "Adsorption locator: " + AdsError;
                if (doc.Summary().Atoms != atomsBefore) AfterEdit("Adsorbates added; the locator stopped (undo with ⌘Z)");
                return;
            }
            double D(string k) => r[k]?.GetValue<double>() ?? 0;
            AdsEnergy = D("adsorption_energy").ToString("0.00", inv) + " kcal/mol";
            AdsSplit = $"adsorbate–substrate {D("adsorbate_substrate").ToString("0.00", inv)} · adsorbate–adsorbate {D("adsorbate_adsorbate").ToString("0.00", inv)} kcal/mol";
            AdsComponents.Clear();
            foreach (var c in (JsonArray)r["components"]!)
                AdsComponents.Add(new AdsComponentRow((string?)c!["name"] ?? "", ((int)c["molecules"]!.GetValue<double>()).ToString(inv), c["de_dn"]!.GetValue<double>().ToString("0.00", inv)));
            AdsConfigs.Clear();
            var k = 0;
            foreach (var c in (JsonArray)r["configs"]!)
                AdsConfigs.Add(new AdsConfigRow((++k).ToString(inv), c!["energy"]!.GetValue<double>().ToString("0.00", inv), (int)c["cycle"]!.GetValue<double>() is var cy && cy > 0 ? cy.ToString(inv) : "best"));
            var edges = ((JsonArray)r["hist_edges"]!).Select(x => x!.GetValue<double>()).ToArray();
            var counts = ((JsonArray)r["hist_counts"]!).Select(x => x!.GetValue<double>()).ToArray();
            AdsHistogram = counts.Select((c, i) => ((edges[i] + edges[i + 1]) / 2, c)).ToArray();
            AdsNote = string.Join(" · ", ((JsonArray)r["notes"]!).Select(x => (string?)x ?? "")) +
                      string.Format(inv, " · acceptance {0:0} % · {1:0.0} s · {2}", 100 * D("acceptance"), D("seconds"), (string?)r["forcefield"] ?? "");
            AfterRun(doc, " · adsorbed");
            Status = $"Adsorption energy {AdsEnergy} · the lowest configuration is shown (earlier frames: the other cycles)";
            AdsorptionChanged?.Invoke();
        }
        catch (Exception e) { AdsError = e.Message; }
        finally { AdsRunning = false; AdsProgress = ""; }
    }
}
