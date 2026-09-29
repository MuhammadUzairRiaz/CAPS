using System.Collections.ObjectModel;
using System.Globalization;
using System.Text.Json.Nodes;

namespace CapsStudio.ViewModels;

public sealed record SorbRow(string Pressure, string Loading, string MolKg, string Cm3, string Heat);

/// <summary>Analyze › Sorption: how much of a gas or small molecule the open structure takes up (core sorption.hpp):
/// Widom insertion (excess chemical potential, Henry constant, solubility coefficient) and grand-canonical Monte Carlo
/// at each pressure (the isotherm, the isosteric heat). The structure is held fixed and is not changed.</summary>
public sealed partial class MainViewModel
{
    public bool IsSorption => _module == 71;
    private bool _sorbRunning;
    private string _sorbate = "O=C=O", _sorbPressures = "10, 50, 100, 500, 1000", _sorbError = "", _sorbProgress = "", _sorbNote = "";
    private double _sorbT = 300;
    private int _sorbInsert = 100000, _sorbSteps = 200000;
    private string _sorbW = "—", _sorbMu = "—", _sorbHenry = "—", _sorbS = "—";
    private CancellationTokenSource? _sorbCts;
    public static readonly string[] SorbPresets = ["O=C=O · CO₂", "C · CH₄", "N#N · N₂", "O=O · O₂", "O · H₂O", "[He] · He", "[H][H] · H₂"];
    public ObservableCollection<SorbRow> SorbRows { get; } = new();
    public (double X, double Y)[] SorbIsotherm { get; private set; } = [];
    public string Sorbate { get => _sorbate; set => Set(ref _sorbate, value ?? ""); }
    public string SorbPreset { get => ""; set { if (value?.Split(" · ")[0] is { Length: > 0 } smi) Sorbate = smi; } }
    public string SorbPressures { get => _sorbPressures; set => Set(ref _sorbPressures, value ?? ""); }
    public decimal SorbTD { get => (decimal)_sorbT; set => Set(ref _sorbT, Math.Clamp((double)value, 1, 5000)); }
    public decimal SorbInsertD { get => _sorbInsert; set => Set(ref _sorbInsert, (int)Math.Clamp(value, 0, 100_000_000)); }
    public decimal SorbStepsD { get => _sorbSteps; set => Set(ref _sorbSteps, (int)Math.Clamp(value, 1000, 1_000_000_000)); }
    public bool SorbRunning { get => _sorbRunning; private set { if (Set(ref _sorbRunning, value)) RaiseBusy(); } }
    public string SorbError { get => _sorbError; private set { if (Set(ref _sorbError, value)) Raise(nameof(SorbHasError)); } }
    public bool SorbHasError => _sorbError.Length > 0;
    public string SorbProgress { get => _sorbProgress; private set => Set(ref _sorbProgress, value); }
    public string SorbNote { get => _sorbNote; private set => Set(ref _sorbNote, value); }
    public string SorbW { get => _sorbW; private set => Set(ref _sorbW, value); }
    public string SorbMu { get => _sorbMu; private set => Set(ref _sorbMu, value); }
    public string SorbHenry { get => _sorbHenry; private set => Set(ref _sorbHenry, value); }
    public string SorbS { get => _sorbS; private set => Set(ref _sorbS, value); }
    public event Action? SorptionChanged;

    public void OpenSorption() => SetModule(71);
    public void CancelSorption() => _sorbCts?.Cancel();

    public async Task RunSorption()
    {
        if (_doc == null || !Idle || BlockedByField("Sorption")) return;
        var inv = CultureInfo.InvariantCulture;
        var pressures = new JsonArray();
        foreach (var t in _sorbPressures.Split([',', ' ', ';'], StringSplitOptions.RemoveEmptyEntries))
            if (double.TryParse(t, NumberStyles.Float, inv, out var p) && p > 0) pressures.Add(p);
        var json = new JsonObject
        {
            ["sorbate"] = _sorbate.Trim(), ["temperature"] = _sorbT, ["insertions"] = _sorbInsert, ["pressures_kpa"] = pressures,
            ["steps"] = _sorbSteps, ["cutoff"] = Math.Min(12.0, _relaxCutoff), ["coulomb"] = _relaxCoulomb, ["seed"] = SorbSeed.Take(),
        }.ToJsonString();
        var doc = _doc!;
        SorbRunning = true;
        SorbError = "";
        _sorbCts = new CancellationTokenSource();
        var ct = _sorbCts.Token;
        Status = $"Sorption of {_sorbate} in {Title}…";
        try
        {
            var res = await Task.Run(() => doc.Sorption(json, (st, f) =>
            {
                Avalonia.Threading.Dispatcher.UIThread.Post(() => SorbProgress = $"{st} · {f * 100:0} %");
                return !ct.IsCancellationRequested;
            }));
            var r = JsonNode.Parse(res)!;
            if (r["ok"]?.GetValue<bool>() != true) { SorbError = (string?)r["error"] ?? "sorption failed"; Status = "Sorption: " + SorbError; return; }
            double D(string k) => r[k]?.GetValue<double>() ?? 0;
            SorbW = $"{D("widom_w").ToString("0.####E+0", inv)} ± {D("widom_error").ToString("0.#E+0", inv)}";
            SorbMu = D("mu_ex").ToString("0.00", inv) + " kcal/mol";
            SorbHenry = D("henry_mol_kg_kpa").ToString("0.###E+0", inv) + " mol/(kg·kPa)";
            SorbS = D("solubility").ToString("0.###", inv) + " cm³(STP)/(cm³·atm)";
            SorbRows.Clear();
            var pts = new List<(double, double)>();
            foreach (var p in (JsonArray)r["isotherm"]!)
            {
                double P(string k) => p![k]!.GetValue<double>();
                SorbRows.Add(new SorbRow(P("pressure_kpa").ToString("0.##", inv), $"{P("loading").ToString("0.##", inv)} ± {P("loading_error").ToString("0.##", inv)}",
                    P("mol_per_kg").ToString("0.###", inv), P("cm3stp_per_cm3").ToString("0.##", inv), P("heat").ToString("0.0", inv)));
                pts.Add((P("pressure_kpa"), P("mol_per_kg")));
            }
            SorbIsotherm = pts.ToArray();
            SorbNote = string.Join(" · ", ((JsonArray)r["notes"]!).Select(x => (string?)x ?? "")) + " · " + ((string?)r["forcefield"] ?? "");
            Status = $"Sorption of {_sorbate}: S {SorbS}" + (SorbRows.Count > 0 ? $" · {SorbRows.Count} isotherm points" : "");
            SorptionChanged?.Invoke();
        }
        catch (Exception e) { SorbError = e.Message; }
        finally { SorbRunning = false; SorbProgress = ""; }
    }
}
