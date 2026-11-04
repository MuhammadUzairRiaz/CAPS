using System.Globalization;
using System.Text;
using System.Text.Json.Nodes;

namespace CapsStudio.ViewModels;

/// <summary>Thermal page › dielectric response (core dielectric.hpp): the open crystal's static and high-frequency dielectric
/// tensors, refractive indices, proper piezoelectric constants and reflectivity from its force field's lattice dynamics.</summary>
public sealed partial class MainViewModel
{
    private bool _dlRunning, _dlHas;
    private string _dlError = "", _dlText = "", _dlNote = "";
    public bool DlRunning { get => _dlRunning; private set { if (Set(ref _dlRunning, value)) RaiseBusy(); } }
    public bool DlHas { get => _dlHas; private set => Set(ref _dlHas, value); }
    public string DlError { get => _dlError; private set { if (Set(ref _dlError, value)) Raise(nameof(DlHasError)); } }
    public bool DlHasError => _dlError.Length > 0;
    public string DlText { get => _dlText; private set => Set(ref _dlText, value); }
    public string DlNote { get => _dlNote; private set => Set(ref _dlNote, value); }
    /// <summary>Reflectivity along the three principal axes of ε⁰.</summary>
    public (double[] W, double[][] R) DlReflectivity { get; private set; } = ([], []);
    public event Action? DielectricChanged;

    public async Task RunDielectric()
    {
        if (_doc == null || !Idle || BlockedByField("Dielectric response")) return;
        var json = new JsonObject { ["cutoff"] = Math.Min(12.0, _relaxCutoff), ["damping"] = 10, ["points"] = 400 }.ToJsonString();
        var doc = _doc!;
        DlRunning = true;
        DlError = "";
        Status = "Dielectric response…";
        try
        {
            var res = await Task.Run(() => doc.Dielectric(json, (st, f) =>
            {
                Avalonia.Threading.Dispatcher.UIThread.Post(() => { if (DlRunning) Status = $"Dielectric response: {st} · {f * 100:0} %"; });
                return true;
            }));
            var r = JsonNode.Parse(res)!;
            if (r["ok"]?.GetValue<bool>() != true) { DlError = (string?)r["error"] ?? "failed"; Status = "Dielectric response: " + DlError; return; }
            double[] A(JsonNode? n) => n is JsonArray a ? a.Select(v => (double?)v ?? 0).ToArray() : [];
            var inv = CultureInfo.InvariantCulture;
            var e0 = A(r["eps0_principal"]); var ei = A(r["epsinf_principal"]); var n = A(r["refractive"]);
            var sb = new StringBuilder();
            sb.AppendLine(string.Format(inv, "ε⁰ {0:0.###} {1:0.###} {2:0.###} · ε^∞ {3:0.###} {4:0.###} {5:0.###}", e0[0], e0[1], e0[2], ei[0], ei[1], ei[2]));
            if ((bool?)r["shells"] == true) sb.AppendLine(string.Format(inv, "refractive indices {0:0.###} {1:0.###} {2:0.###}", n[0], n[1], n[2]));
            var e = ((JsonArray)r["e"]!).Select(A).ToArray();
            var emax = e.SelectMany(x => x).Select(Math.Abs).DefaultIfEmpty(0).Max();
            if (emax > 1e-5)
            {
                sb.AppendLine("piezoelectric e (C/m²), xx yy zz yz xz xy:");
                foreach (var row in e) sb.AppendLine("  " + string.Join("  ", row.Select(v => v.ToString("0.0000", inv).PadLeft(8))));
                var d = ((JsonArray)r["d"]!).Select(A).ToArray();
                sb.AppendLine("d (pC/N):");
                foreach (var row in d) sb.AppendLine("  " + string.Join("  ", row.Select(v => v.ToString("0.000", inv).PadLeft(8))));
            }
            else sb.AppendLine("no piezoelectricity (a centrosymmetric crystal)");
            var modes = A(r["modes"]);
            if (modes.Length > 0) sb.Append(string.Format(inv, "zone-centre modes {0:0} … {1:0} cm⁻¹", modes.Min(), modes.Max()));
            DlText = sb.ToString();
            DlNote = string.Join(" · ", ((JsonArray)r["notes"]!).Select(v => (string?)v ?? ""));
            DlReflectivity = (A(r["w"]), ((JsonArray)r["reflectivity"]!).Select(A).ToArray());
            DlHas = true;
            Status = string.Format(inv, "Dielectric response: ε⁰ {0:0.##}", e0.Average());
            DielectricChanged?.Invoke();
        }
        catch (Exception ex) { DlError = ex.Message; }
        finally { DlRunning = false; }
    }
}
