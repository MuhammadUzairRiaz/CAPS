using System.Globalization;
using System.Text;
using System.Text.Json.Nodes;

namespace CapsStudio.ViewModels;

/// <summary>Analyze › Scattering › atomic displacement parameters (core adp.hpp): U_ij, U_iso and B per atom or per site
/// of the unit cell from the open trajectory; the view coloured by U_iso; a CIF of the mean positions with U_ij.</summary>
public sealed partial class MainViewModel
{
    private decimal _adpNa = 1, _adpNb = 1, _adpNc = 1;
    private string _adpText = "", _adpError = "";
    private bool _adpRunning;
    private double[] _adpUiso = [];
    public decimal AdpNa { get => _adpNa; set => Set(ref _adpNa, Math.Clamp(value, 1, 20)); }
    public decimal AdpNb { get => _adpNb; set => Set(ref _adpNb, Math.Clamp(value, 1, 20)); }
    public decimal AdpNc { get => _adpNc; set => Set(ref _adpNc, Math.Clamp(value, 1, 20)); }
    public string AdpText { get => _adpText; private set { if (Set(ref _adpText, value)) Raise(nameof(AdpHas)); } }
    public bool AdpHas => _adpText.Length > 0;
    public string AdpError { get => _adpError; private set { if (Set(ref _adpError, value)) Raise(nameof(AdpHasError)); } }
    public bool AdpHasError => _adpError.Length > 0;
    public bool AdpRunning { get => _adpRunning; private set { if (Set(ref _adpRunning, value)) RaiseBusy(); } }

    private string AdpJson(string cif = "") => new JsonObject
    {
        ["supercell"] = new JsonArray((int)_adpNa, (int)_adpNb, (int)_adpNc), ["probability"] = 0.5, ["cif"] = cif,
    }.ToJsonString();

    public async Task RunAdp(string cif = "")
    {
        if (_doc == null || !Idle) return;
        var doc = _doc;
        AdpRunning = true;
        AdpError = "";
        try
        {
            var json = AdpJson(cif);
            var r = JsonNode.Parse(await Task.Run(() => doc.Adp(json)))!;
            if ((bool?)r["ok"] != true) { AdpError = (string?)r["error"] ?? "failed"; return; }
            var inv = CultureInfo.InvariantCulture;
            var atoms = (JsonArray)r["atoms"]!;
            _adpUiso = atoms.Select(a => (double?)a!["u_iso"] ?? 0).ToArray();
            var sites = (JsonArray)r["sites"]!;
            var rows = sites.Count > 0 ? sites : atoms;
            var sb = new StringBuilder();
            sb.AppendLine(string.Format(inv, "{0:0} frames · {1} {2} · mean U_iso {3:0.0000} Å², B {4:0.00} Å²", (double?)r["frames"] ?? 0, rows.Count, sites.Count > 0 ? "sites" : "atoms",
                                        _adpUiso.Average(), 8 * Math.PI * Math.PI * _adpUiso.Average()));
            if (sites.Count > 0)
            {
                sb.AppendLine($"{"site",-6}{"U_iso Å²",10}{"B Å²",8}{"λmax/λmin",11}");
                for (var k = 0; k < sites.Count && k < 24; ++k)
                    sb.AppendLine(string.Format(inv, "{0,-6}{1,10:0.0000}{2,8:0.00}{3,11:0.00}", k + 1, (double?)sites[k]!["u_iso"] ?? 0, (double?)sites[k]!["b"] ?? 0, (double?)sites[k]!["anisotropy"] ?? double.NaN));
            }
            foreach (var n in (JsonArray)r["notes"]!) sb.AppendLine((string?)n ?? "");
            if (r["cif"] is JsonNode c) sb.Append("wrote " + Path.GetFileName((string?)c ?? ""));
            AdpText = sb.ToString().TrimEnd();
            Status = "Displacement parameters: mean B " + (8 * Math.PI * Math.PI * _adpUiso.Average()).ToString("0.00", inv) + " Å²";
        }
        catch (Exception e) { AdpError = e.Message; }
        finally { AdpRunning = false; }
    }

    /// <summary>The view coloured by each atom's U_iso (viridis: mobile atoms bright).</summary>
    public async Task ColourByMobility()
    {
        if (_doc == null) return;
        if (_adpUiso.Length != _doc.Summary().Atoms) await RunAdp();
        if (_adpUiso.Length == 0 || _doc == null) return;
        _doc.SetAtomValues(_adpUiso, 0);
        RenderRequested?.Invoke();
        Status = "Coloured by U_iso (viridis, mobile atoms bright); clear it from Colour by";
    }
}
