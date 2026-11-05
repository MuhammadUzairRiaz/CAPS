using System;
using System.Globalization;
using System.Linq;
using System.Text.Json.Nodes;
using CapsStudio.Interop;

namespace CapsStudio.ViewModels;

// Unit cell page › Lattice planes (core miller.hpp through caps_miller): the planes (hkl) cut to the cell and drawn in the
// view, their spacing, and the camera turned to look straight down their normal.
public partial class MainViewModel
{
    private decimal _millerH = 1, _millerK = 1, _millerL = 1, _millerShift;
    public decimal MillerH { get => _millerH; set => Set(ref _millerH, value); }
    public decimal MillerK { get => _millerK; set => Set(ref _millerK, value); }
    public decimal MillerL { get => _millerL; set => Set(ref _millerL, value); }
    /// <summary>Where the planes sit, in units of d (0: through the origin).</summary>
    public decimal MillerShift { get => _millerShift; set => Set(ref _millerShift, value); }
    private string _millerText = "", _millerError = "";
    public string MillerText { get => _millerText; private set => Set(ref _millerText, value); }
    public string MillerError { get => _millerError; private set { if (Set(ref _millerError, value)) Raise(nameof(MillerHasError)); } }
    public bool MillerHasError => _millerError.Length > 0;
    private bool _millerShown;
    public bool MillerShown { get => _millerShown; private set => Set(ref _millerShown, value); }
    private double _millerYaw, _millerPitch;

    /// <summary>Draws the planes (hkl) in the view and reports their spacing; view: also turn the camera down the normal.</summary>
    public void ShowMiller(bool view = false)
    {
        if (_doc == null) return;
        var inv = CultureInfo.InvariantCulture;
        var json = new JsonObject { ["h"] = (int)_millerH, ["k"] = (int)_millerK, ["l"] = (int)_millerL, ["shift"] = (double)_millerShift, ["show"] = true }.ToJsonString();
        try
        {
            var r = JsonNode.Parse(_doc.Miller(json))!;
            if (r["ok"]?.GetValue<bool>() != true) { MillerError = (string?)r["error"] ?? "no planes"; MillerShown = false; RenderRequested?.Invoke(); return; }
            MillerError = "";
            MillerShown = true;
            _millerYaw = r["yaw"]!.GetValue<double>();
            _millerPitch = r["pitch"]!.GetValue<double>();
            var nrm = ((JsonArray)r["normal"]!).Select(x => x!.GetValue<double>().ToString("0.000", inv));
            MillerText = $"{(string?)r["label"]}  d = {r["d"]!.GetValue<double>().ToString("0.0000", inv)} Å · {(int)r["planes"]!.GetValue<double>()} plane(s) in the cell · normal ({string.Join(", ", nrm)})";
            Status = "Lattice planes " + MillerText;
            if (view) SetView(_millerYaw, _millerPitch);
            RenderRequested?.Invoke();
        }
        catch (Exception e) { MillerError = e.Message; }
    }
    public void ViewOntoMiller() => ShowMiller(true);
    public void HideMiller()
    {
        if (_doc == null) return;
        try { _doc.Miller("{\"clear\": true}"); } catch { }
        MillerShown = false;
        MillerText = "";
        RenderRequested?.Invoke();
    }
}
