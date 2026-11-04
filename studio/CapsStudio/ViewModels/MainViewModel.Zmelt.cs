using System.Globalization;
using System.Text;
using System.Text.Json.Nodes;

namespace CapsStudio.ViewModels;

/// <summary>Thermal page › melting point by the Z method (core zmethod.hpp): NVE runs of the open crystal from its perfect
/// lattice at rising energy; the temperature against the energy draws a Z whose bottom is the melting point.</summary>
public sealed partial class MainViewModel
{
    private bool _zmRunning, _zmHas;
    private string _zmError = "", _zmText = "", _zmNote = "", _zmTemps = "300, 600, 900, 1200";
    private decimal _zmSteps = 10000, _zmBisect = 3, _zmDt = 1;
    private CancellationTokenSource? _zmCancel;
    public bool ZmRunning { get => _zmRunning; private set { if (Set(ref _zmRunning, value)) RaiseBusy(); } }
    public bool ZmHas { get => _zmHas; private set => Set(ref _zmHas, value); }
    public string ZmError { get => _zmError; private set { if (Set(ref _zmError, value)) Raise(nameof(ZmHasError)); } }
    public bool ZmHasError => _zmError.Length > 0;
    public string ZmText { get => _zmText; private set => Set(ref _zmText, value); }
    public string ZmNote { get => _zmNote; private set => Set(ref _zmNote, value); }
    public string ZmTemperatures { get => _zmTemps; set => Set(ref _zmTemps, value); }
    public decimal? ZmSteps { get => _zmSteps; set { if (value != null) Set(ref _zmSteps, Math.Clamp(value.Value, 500, 1000000)); } }
    public decimal? ZmBisect { get => _zmBisect; set { if (value != null) Set(ref _zmBisect, Math.Clamp(value.Value, 0, 8)); } }
    public decimal? ZmDt { get => _zmDt; set { if (value != null) Set(ref _zmDt, Math.Clamp(value.Value, 0.1m, 10)); } }
    /// <summary>The runs by energy: E (kcal/mol per atom), mean T (K), melted.</summary>
    public (double[] E, double[] T, bool[] Melted, double Tm, double Tls) ZmCurve { get; private set; } = ([], [], [], 0, 0);
    public event Action? ZmeltChanged;

    public void CancelZmelt() => _zmCancel?.Cancel();

    public async Task RunZmelt()
    {
        if (_doc == null || !Idle || BlockedByField("Z-method melting")) return;
        var inv = CultureInfo.InvariantCulture;
        var temps = new JsonArray();
        foreach (var t in _zmTemps.Split([',', ' ', ';'], StringSplitOptions.RemoveEmptyEntries))
        {
            if (!double.TryParse(t, NumberStyles.Float, inv, out var v) || v <= 0) { ZmError = $"'{t}' is not a temperature in K"; return; }
            temps.Add(v);
        }
        if (temps.Count == 0) { ZmError = "Give the starting temperatures, in K, from well below to above melting."; return; }
        var json = new JsonObject
        {
            ["temperatures"] = temps, ["steps"] = (long)_zmSteps, ["dt"] = (double)_zmDt, ["bisect"] = (int)_zmBisect, ["cutoff"] = Math.Min(12.0, _relaxCutoff),
        }.ToJsonString();
        var doc = _doc!;
        _zmCancel = new CancellationTokenSource();
        var token = _zmCancel.Token;
        ZmRunning = true;
        ZmError = "";
        Status = "Z method…";
        try
        {
            var res = await Task.Run(() => doc.Zmelt(json, (st, f) =>
            {
                Avalonia.Threading.Dispatcher.UIThread.Post(() => { if (ZmRunning) Status = $"Z method: {st} · {f * 100:0} %"; });
                return !token.IsCancellationRequested;
            }));
            var r = JsonNode.Parse(res)!;
            if (r["ok"]?.GetValue<bool>() != true) { ZmError = (string?)r["error"] ?? "failed"; Status = "Z method: " + ZmError; return; }
            var pts = (JsonArray)r["points"]!;
            var sb = new StringBuilder();
            sb.AppendLine($"{"T0 K",7}{"E kcal/mol",13}{"T K",9}{"P atm",10}{"order",7}  state");
            foreach (var p in pts)
                sb.AppendLine(string.Format(inv, "{0,7:0}{1,13:0.0000}{2,9:0.0}{3,10:0}{4,7:0.00}  {5}", (double?)p!["t0"] ?? 0, (double?)p["energy"] ?? 0, (double?)p["temperature"] ?? 0,
                    (double?)p["pressure"] ?? 0, (double?)p["order"] ?? 0, (bool?)p["melted"] == true ? "melted" : "crystal"));
            var tm = (double?)r["t_melt"] ?? 0;
            var tls = (double?)r["t_superheat"] ?? 0;
            if (tm > 0) sb.Append(string.Format(inv, "melting T_m {0:0} K · superheating limit T_LS {1:0} K (at the crystal's volume)", tm, tls));
            ZmText = sb.ToString().TrimEnd();
            ZmNote = string.Join(" · ", ((JsonArray)r["notes"]!).Select(v => (string?)v ?? "").Where(v => !v.StartsWith("superheating limit")));
            ZmCurve = (pts.Select(p => (double?)p!["energy"] ?? 0).ToArray(), pts.Select(p => (double?)p!["temperature"] ?? 0).ToArray(),
                       pts.Select(p => (bool?)p!["melted"] == true).ToArray(), tm, tls);
            ZmHas = true;
            Status = tm > 0 ? string.Format(inv, "Z method: T_m {0:0} K", tm) : "Z method: no run melted";
            ZmeltChanged?.Invoke();
        }
        catch (Exception ex) { ZmError = ex.Message.Contains("cancel") ? "stopped" : ex.Message; }
        finally { ZmRunning = false; _zmCancel = null; }
    }
}
