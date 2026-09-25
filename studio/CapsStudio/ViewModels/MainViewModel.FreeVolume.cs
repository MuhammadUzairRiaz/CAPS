using System.Globalization;
using System.Text.Json.Nodes;

namespace CapsStudio.ViewModels;

/// <summary>Analyze › Free volume (design/boards/FreeVolume): probe insertion on a grid — the fractional free volume
/// (Bondi), the share a probe reaches, the largest empty spheres drawn in the cell, and the pore size distribution.</summary>
public sealed partial class MainViewModel
{
    public bool IsFreeVolume => _module == 40;
    public ResultCell FvFfv { get; } = new("Fractional free volume", "Bondi: 1 − 1.3 V_w/V");
    public ResultCell FvAccessible { get; } = new("Accessible to probe", "probe 1.4 Å");
    private decimal _fvProbe = 1.4m, _fvGrid = 0.5m;
    private int _fvCount = 40;
    public decimal FvProbe { get => _fvProbe; set { if (Set(ref _fvProbe, Math.Clamp(value, 0m, 5m))) Analyze.ProbeD = value; } }
    public decimal FvGrid { get => _fvGrid; set { if (Set(ref _fvGrid, Math.Clamp(value, 0.2m, 2m))) Analyze.GridD = value; } }
    public decimal FvCount { get => _fvCount; set => Set(ref _fvCount, (int)Math.Clamp(value, 1, 500)); }
    private string _fvChip = "", _fvStatus = "", _fvFooter = "";
    public string FvChip { get => _fvChip; private set => Set(ref _fvChip, value); }
    public string FvStatus { get => _fvStatus; private set => Set(ref _fvStatus, value); }
    public string FvFooter { get => _fvFooter; private set => Set(ref _fvFooter, value); }
    public bool FvHasVoids => _fvVoids > 0;
    private int _fvVoids;
    public (double X, double Y)[] FvPsd { get; private set; } = [];
    public bool FvHasPsd => FvPsd.Length > 1;
    public event Action? FvChanged;

    public void OpenFreeVolume()
    {
        SetModule(40);
        FvChanged?.Invoke();
    }

    /// <summary>The voids (drawn in the view) and the Analyze free volume and pore size distribution.</summary>
    public async Task RunFreeVolume()
    {
        if (_doc == null || Analyze.Working) return;
        var doc = _doc;
        var inv = CultureInfo.InvariantCulture;
        var opts = string.Format(inv, "{{\"grid\":{0},\"probe\":{1},\"count\":{2},\"show\":true}}", _fvGrid, _fvProbe, _fvCount);
        FvStatus = "Finding voids…";
        string json;
        try { json = await Task.Run(() => doc.Voids(opts)); }
        catch (Exception e) { FvStatus = e.Message; return; }
        var j = JsonNode.Parse(json)!;
        if (j["ok"]?.GetValue<bool>() != true) { FvStatus = j["error"]?.GetValue<string>() ?? "voids failed"; return; }
        _fvVoids = ((JsonArray)j["spheres"]!).Count;
        var acc = j["accessible_probe"]!.GetValue<double>();
        FvAccessible.Value = (acc * 100).ToString("0.0", inv) + " %";
        FvAccessible.Caption = $"probe {_fvProbe.ToString("0.##", inv)} Å · point probe {(j["accessible_point"]!.GetValue<double>() * 100).ToString("0.0", inv)} %";
        var g = (JsonArray)j["grid"]!;
        FvChip = $"{_fvVoids} largest voids · biggest r = {j["largest"]!.GetValue<double>().ToString("0.0", inv)} Å";
        FvFooter = $"Grid {g[0]} × {g[1]} × {g[2]} points · {_fvGrid.ToString("0.0#", inv)} Å · current frame";
        Raise(nameof(FvHasVoids));
        RenderRequested?.Invoke();
        FvChanged?.Invoke();
        // the Analyze engine: Bondi FFV and the pore size distribution over the chosen frames
        FvStatus = "Pore size distribution…";
        var chips = Analyze.Groups.SelectMany(c => c.Chips).ToList();
        var was = chips.ToDictionary(c => c, c => c.IsOn);
        foreach (var c in chips) c.IsOn = c.Id is "ffv" or "psd";
        Analyze.ProbeD = _fvProbe;
        Analyze.GridD = _fvGrid;
        try { await Analyze.Run(); }
        finally { foreach (var (c, on) in was) c.IsOn = on; }
        var ffv = Analyze.Results.FirstOrDefault(r => r.Id == "ffv");
        var bondi = ffv?.Extra.FirstOrDefault(e => e.Key.StartsWith("Bondi FFV", StringComparison.Ordinal));
        FvFfv.Value = bondi is { Key: not null } b ? b.Value.ToString("0.000", inv) : "—";
        var psdCard = Analyze.Results.FirstOrDefault(r => r.Id == "psd");
        var psd = psdCard == null ? null : Analyze.Curves.FirstOrDefault(c => c.Property == psdCard.Name);
        FvPsd = psd == null ? [] : psd.X.Zip(psd.Y).Where(p => double.IsFinite(p.Second)).ToArray();
        FvStatus = Analyze.Log;
        Raise(nameof(FvHasPsd));
        FvChanged?.Invoke();
    }

    public void ExportVoids(string path)
    {
        if (_doc == null) return;
        try { _doc.VoidsPdb(path); Status = $"Wrote {Path.GetFileName(path)} · {_fvVoids} voids (radius in the B-factor)"; }
        catch (Exception e) { Status = e.Message; }
    }
}
