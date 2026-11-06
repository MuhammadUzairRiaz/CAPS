using System.Globalization;
using System.Text;
using System.Text.Json.Nodes;

namespace CapsStudio.ViewModels;

/// <summary>Thermal conductivity › Thermal expansion (core qha.hpp): phonons at strained cells; the minimum of the free
/// energy at each temperature gives the cell, the expansion coefficients, the bulk modulus and C_p.</summary>
public sealed partial class MainViewModel
{
    public static readonly string[] QhaModes = ["As a whole (volume)", "Each axis on its own (a, b, c)"];
    private int _qhaMode;
    private decimal _qhaRange = 0.02m;
    private string _qhaTemps = "0, 100, 200, 300", _qhaResult = "", _qhaNote = "", _qhaError = "";
    private bool _qhaRunning;
    public int QhaMode { get => _qhaMode; set => Set(ref _qhaMode, Math.Clamp(value, 0, 1)); }
    public decimal QhaRange { get => _qhaRange; set => Set(ref _qhaRange, value); }
    public string QhaTemps { get => _qhaTemps; set => Set(ref _qhaTemps, value ?? ""); }
    public string QhaResult { get => _qhaResult; private set => Set(ref _qhaResult, value); }
    public string QhaNote { get => _qhaNote; private set => Set(ref _qhaNote, value); }
    public string QhaError { get => _qhaError; private set { if (Set(ref _qhaError, value)) Raise(nameof(QhaHasError)); } }
    public bool QhaHasError => _qhaError.Length > 0;
    public bool QhaRunning { get => _qhaRunning; private set { if (Set(ref _qhaRunning, value)) RaiseBusy(); } }

    public async Task RunQha()
    {
        if (_doc == null || !Idle || BlockedByField("Thermal expansion")) return;
        var inv = CultureInfo.InvariantCulture;
        var temps = new JsonArray();
        foreach (var t in _qhaTemps.Split([',', ' ', ';'], StringSplitOptions.RemoveEmptyEntries))
        {
            if (!double.TryParse(t, NumberStyles.Float, inv, out var v) || v < 0) { QhaError = $"not a temperature: {t}"; return; }
            temps.Add(v);
        }
        if (temps.Count == 0) { QhaError = "give the temperatures in K"; return; }
        var json = new JsonObject
        {
            ["mode"] = _qhaMode == 1 ? "axes" : "volume", ["range"] = (double)_qhaRange, ["temperatures"] = temps, ["cutoff"] = Math.Min(12.0, _relaxCutoff),
        }.ToJsonString();
        var doc = _doc!;
        QhaRunning = true;
        QhaError = "";
        Status = "Thermal expansion…";
        try
        {
            var res = await Task.Run(() => doc.Qha(json, (st, f) =>
            {
                Avalonia.Threading.Dispatcher.UIThread.Post(() => { if (QhaRunning) Status = $"Thermal expansion: {st} · {f * 100:0} %"; });
                return true;
            }));
            var r = JsonNode.Parse(res)!;
            if (r["ok"]?.GetValue<bool>() != true) { QhaError = (string?)r["error"] ?? "the run failed"; Status = "Thermal expansion: " + QhaError; return; }
            var sb = new StringBuilder();
            sb.AppendLine("  T/K     a/Å     b/Å     c/Å     V/Å³  αa,b,c /10⁻⁶ K⁻¹   B/GPa     Cp");
            foreach (var row in (JsonArray)r["rows"]!)
            {
                double D(string k) => row![k]!.GetValue<double>();
                double[] L(string k) => ((JsonArray)row![k]!).Select(x => x!.GetValue<double>()).ToArray();
                var l = L("length");
                var al = L("alpha");
                sb.AppendLine(string.Format(inv, "{0,5:0} {1,7:0.000} {2,7:0.000} {3,7:0.000} {4,8:0.00} {5,6:0} {6,5:0} {7,5:0} {8,7:0.00} {9,6:0.0}{10}",
                    D("T"), l[0], l[1], l[2], D("volume"), al[0] * 1e6, al[1] * 1e6, al[2] * 1e6, D("B"), D("Cp"),
                    row!["extrapolated"]!.GetValue<bool>() ? " *" : ""));
            }
            QhaResult = sb.ToString().TrimEnd();
            var n = ((JsonArray)r["samples"]!).Count;
            QhaNote = $"{n} strained cells · C_p in cal/mol/K per cell · * extrapolated · " + string.Join(" · ", ((JsonArray)r["notes"]!).Select(v => (string?)v ?? ""));
            Status = $"Thermal expansion: {n} strained cells done";
        }
        catch (Exception e) { QhaError = e.Message; }
        finally { QhaRunning = false; }
    }
}
