using System.Collections.ObjectModel;
using System.Globalization;
using System.Text;
using System.Text.Json.Nodes;
using CapsStudio.Interop;

namespace CapsStudio.ViewModels;

/// <summary>A DPD species: its name, bead sequence (compact "A5B5" or written out) and how many molecules.</summary>
public sealed class DpdSpeciesRow : ObservableObject
{
    private string _name = "", _sequence = "A5B5";
    private decimal _count = 150;
    public string Name { get => _name; set => Set(ref _name, value ?? ""); }
    public string Sequence { get => _sequence; set => Set(ref _sequence, value ?? ""); }
    public decimal Count { get => _count; set => Set(ref _count, Math.Clamp(value, 1, 1_000_000)); }
}

/// <summary>A χ between two bead types.</summary>
public sealed class DpdChiRow : ObservableObject
{
    private string _pair = "AB";
    private decimal _chi = 1.0m;
    public string Pair { get => _pair; set => Set(ref _pair, (value ?? "").ToUpperInvariant()); }
    public decimal Chi { get => _chi; set => Set(ref _chi, Math.Clamp(value, -20, 50)); }
}

/// <summary>Polymer cell › Mesoscale (DPD): coarse-grained beads by dissipative particle dynamics (core dpd.hpp,
/// Groot–Warren) — block copolymers and blends from χ: microdomains, their spacing, the segregation over the run. The
/// run's frames open as a new structure.</summary>
public sealed partial class MainViewModel
{
    public bool IsDpd => _module == 72;
    private bool _dpdRunning;
    private decimal _dpdSteps = 20000, _dpdDensity = 3, _dpdRc = 6.46m;
    private string _dpdError = "", _dpdProgress = "", _dpdNote = "", _dpdOrder = "—", _dpdSpacing = "—", _dpdKt = "—";
    private CancellationTokenSource? _dpdCts;
    public ObservableCollection<DpdSpeciesRow> DpdSpecies { get; } = new() { new DpdSpeciesRow { Name = "A-b-B diblock", Sequence = "A5B5", Count = 150 } };
    public ObservableCollection<DpdChiRow> DpdChis { get; } = new() { new DpdChiRow { Pair = "AB", Chi = 4.3m } };
    public (double X, double Y)[] DpdSq { get; private set; } = [];
    public (double X, double Y)[] DpdOrderSeries { get; private set; } = [];
    public decimal DpdSteps { get => _dpdSteps; set => Set(ref _dpdSteps, Math.Clamp(Math.Round(value), 100, 100_000_000)); }
    public decimal DpdDensity { get => _dpdDensity; set => Set(ref _dpdDensity, Math.Clamp(value, 1, 10)); }
    public decimal DpdRc { get => _dpdRc; set => Set(ref _dpdRc, Math.Clamp(value, 1, 50)); }
    public bool DpdRunning { get => _dpdRunning; private set { if (Set(ref _dpdRunning, value)) RaiseBusy(); } }
    public string DpdError { get => _dpdError; private set { if (Set(ref _dpdError, value)) Raise(nameof(DpdHasError)); } }
    public bool DpdHasError => _dpdError.Length > 0;
    public string DpdProgress { get => _dpdProgress; private set => Set(ref _dpdProgress, value); }
    public string DpdNote { get => _dpdNote; private set => Set(ref _dpdNote, value); }
    public string DpdOrder { get => _dpdOrder; private set => Set(ref _dpdOrder, value); }
    public string DpdSpacing { get => _dpdSpacing; private set => Set(ref _dpdSpacing, value); }
    public string DpdKt { get => _dpdKt; private set => Set(ref _dpdKt, value); }
    public event Action? DpdChanged;

    public void OpenDpd() => SetModule(72);
    public void AddDpdSpecies() => DpdSpecies.Add(new DpdSpeciesRow { Name = "solvent", Sequence = "C", Count = 300 });
    public void RemoveDpdSpecies(DpdSpeciesRow r) { if (DpdSpecies.Count > 1) DpdSpecies.Remove(r); }
    public void AddDpdChi() => DpdChis.Add(new DpdChiRow { Pair = "AC", Chi = 0.5m });
    public void RemoveDpdChi(DpdChiRow r) => DpdChis.Remove(r);
    public void CancelDpd() => _dpdCts?.Cancel();

    /// <summary>"A5B5" → "AAAAABBBBB"; a sequence written out stays as it is.</summary>
    public static string ExpandSequence(string s)
    {
        var sb = new StringBuilder();
        s = s.Replace(" ", "").ToUpperInvariant();
        for (var i = 0; i < s.Length;)
        {
            var c = s[i++];
            if (c < 'A' || c > 'Z') throw new FormatException($"'{c}' is not a bead type (A–Z)");
            var j = i;
            while (j < s.Length && char.IsDigit(s[j])) ++j;
            var n = j > i ? int.Parse(s[i..j], CultureInfo.InvariantCulture) : 1;
            sb.Append(c, n);
            i = j;
        }
        return sb.ToString();
    }

    public async Task RunDpd()
    {
        if (!Idle) return;
        var inv = CultureInfo.InvariantCulture;
        JsonArray species;
        try
        {
            species = new JsonArray(DpdSpecies.Select(r => (JsonNode)new JsonObject
            {
                ["name"] = r.Name, ["sequence"] = ExpandSequence(r.Sequence), ["count"] = (int)r.Count,
            }).ToArray());
        }
        catch (Exception e) { DpdError = e.Message; return; }
        var chi = new JsonObject();
        foreach (var c in DpdChis.Where(c => c.Pair.Length == 2)) chi[c.Pair] = (double)c.Chi;
        var json = new JsonObject
        {
            ["species"] = species, ["chi"] = chi, ["density"] = (double)_dpdDensity, ["steps"] = (long)_dpdSteps,
            ["rc_angstrom"] = (double)_dpdRc, ["seed"] = DpdSeed.Take(),
        }.ToJsonString();
        DpdRunning = true;
        DpdError = "";
        _dpdCts = new CancellationTokenSource();
        var ct = _dpdCts.Token;
        Status = "Running DPD…";
        try
        {
            var title = string.Join(" + ", DpdSpecies.Select(r => r.Name)) + " (DPD)";
            var (doc, rep) = await Task.Run(() => CapsDocument.Dpd(json, (st, f) =>
            {
                Avalonia.Threading.Dispatcher.UIThread.Post(() => DpdProgress = $"{st} · {f * 100:0} %");
                return !ct.IsCancellationRequested;
            }, title));
            var r = JsonNode.Parse(rep)!;
            if (doc == null || r["ok"]?.GetValue<bool>() != true) { DpdError = (string?)r["error"] ?? "the run failed"; Status = "DPD: " + DpdError; return; }
            double D(string k) => r[k]?.GetValue<double>() ?? 0;
            DpdOrder = D("order").ToString("0.00", inv);
            DpdSpacing = D("spacing") > 0 ? $"{D("spacing").ToString("0.0", inv)} r_c · {(D("spacing") * (double)_dpdRc).ToString("0", inv)} Å" : "—";
            DpdKt = $"{D("kT").ToString("0.000", inv)} ± {D("kT_error").ToString("0.000", inv)}";
            var q = ((JsonArray)r["q"]!).Select(x => x!.GetValue<double>()).ToArray();
            var sq = ((JsonArray)r["sq"]!).Select(x => x!.GetValue<double>()).ToArray();
            DpdSq = q.Zip(sq).Select(t => (t.First, t.Second)).ToArray();
            DpdOrderSeries = ((JsonArray)r["order_series"]!).Select(p => (p![0]!.GetValue<double>(), p[1]!.GetValue<double>())).ToArray();
            DpdNote = string.Join(" · ", ((JsonArray)r["notes"]!).Select(x => (string?)x ?? ""));
            DpdRunning = false;   // the run is over: the frames can open (Show waits for no run)
            Show(doc, title);
            GrownUnsaved = true;
            SetModule(72);
            Status = $"DPD: ψ {DpdOrder} · spacing {DpdSpacing} · the frames are open as a new structure";
            DpdChanged?.Invoke();
        }
        catch (Exception e) { DpdError = e.Message; }
        finally { DpdRunning = false; DpdProgress = ""; }
    }
}
