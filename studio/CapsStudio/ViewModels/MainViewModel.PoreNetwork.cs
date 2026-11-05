using System.Globalization;
using System.Text.Json.Nodes;

namespace CapsStudio.ViewModels;

/// <summary>Analyze › Free volume › Pore network (core pores.hpp): the probe's channels and pockets, the largest included
/// and pore-limiting diameters, accessible and non-accessible volume and surface, from the current frame.</summary>
public sealed partial class MainViewModel
{
    private string _pnDiameters = "—", _pnVolume = "—", _pnSurface = "—", _pnClusters = "", _pnNote = "";
    private bool _pnWorking;
    private decimal _pnSamples = 2000;
    public string PnDiameters { get => _pnDiameters; private set => Set(ref _pnDiameters, value); }
    public string PnVolume { get => _pnVolume; private set => Set(ref _pnVolume, value); }
    public string PnSurface { get => _pnSurface; private set => Set(ref _pnSurface, value); }
    public string PnClusters { get => _pnClusters; private set => Set(ref _pnClusters, value); }
    public string PnNote { get => _pnNote; private set => Set(ref _pnNote, value); }
    public bool PnWorking { get => _pnWorking; private set { if (Set(ref _pnWorking, value)) Raise(nameof(PnIdle)); } }
    public bool PnIdle => !_pnWorking;
    public decimal PnSamples { get => _pnSamples; set => Set(ref _pnSamples, Math.Clamp(Math.Round(value), 50, 100000)); }

    /// <summary>The pore network of the open frame with the page's probe, grid and radii.</summary>
    public async Task RunPoreNetwork()
    {
        if (_doc == null || _pnWorking) return;
        var doc = _doc;
        var inv = CultureInfo.InvariantCulture;
        var opts = new JsonObject
        {
            ["probe"] = (double)_fvProbe, ["grid"] = Math.Min(0.25, (double)_fvGrid), ["samples"] = (int)_pnSamples, ["radii"] = FvRadiiIds[_fvRadii],
        }.ToJsonString();
        PnWorking = true;
        PnNote = "Finding channels and pockets…";
        try
        {
            var json = await Task.Run(() => doc.Pores(opts));
            var j = JsonNode.Parse(json)!;
            if (j["ok"]?.GetValue<bool>() != true) { PnNote = j["error"]?.GetValue<string>() ?? "the pore network failed"; return; }
            double D(string k) => j[k]?.GetValue<double>() ?? 0;
            var pa = ((JsonArray)j["pld_axis"]!).Select(x => x!.GetValue<double>()).ToArray();
            PnDiameters = string.Format(inv, "LCD {0:0.00} Å · PLD {1:0.00} Å (a {2:0.00}, b {3:0.00}, c {4:0.00}) · Dif {5:0.00} Å", D("lcd"), D("pld"), pa[0], pa[1], pa[2], D("dif"));
            PnVolume = string.Format(inv, "channels {0:0.0} % ({1:0.000} cm³/g) · pockets {2:0.0} % ({3:0.000} cm³/g) · geometric {4:0.0} %",
                100 * D("av_fraction"), D("av"), 100 * D("nav_fraction"), D("nav"), 100 * D("geometric_fraction"));
            PnSurface = string.Format(inv, "accessible {0:0.0} Å² · {1:0} m²/cm³ · {2:0} m²/g; non-accessible {3:0.0} Å² · {4:0} m²/g",
                D("asa"), D("asa_volume"), D("asa_mass"), D("nasa"), D("nasa_mass"));
            var cl = (JsonArray)j["clusters"]!;
            var channels = cl.Count(c => c!["channel"]!.GetValue<bool>());
            PnClusters = $"{channels} channel(s) along {(int)D("channel_dims")} axis(es), {cl.Count - channels} pocket(s)";
            PnNote = string.Join(" · ", ((JsonArray)j["notes"]!).Select(x => (string?)x ?? ""));
        }
        catch (Exception e) { PnNote = e.Message; }
        finally { PnWorking = false; }
    }
}
