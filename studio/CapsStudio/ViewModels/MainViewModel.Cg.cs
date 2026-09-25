using System.Globalization;
using System.Text.Json.Nodes;
using CapsStudio.Interop;

namespace CapsStudio.ViewModels;

/// <summary>Coarse-grained melts (design/boards/CoarseGrained): Kremer–Grest bead-spring chains built as random walks at a
/// reduced density, drawn, and written for LAMMPS with the push-off and the FENE + WCA run.</summary>
public sealed partial class MainViewModel
{
    public bool IsCg => _module == 44;
    private decimal _cgChains = 50, _cgBeads = 100, _cgDensity = 0.85m, _cgKTheta = 0, _cgSeed = 1, _cgPush = 20000, _cgRun = 100000;
    public decimal CgChains { get => _cgChains; set { if (Set(ref _cgChains, Math.Clamp(Math.Round(value), 1, 20000))) CgPreview(); } }
    public decimal CgBeads { get => _cgBeads; set { if (Set(ref _cgBeads, Math.Clamp(Math.Round(value), 2, 100000))) CgPreview(); } }
    public decimal CgDensity { get => _cgDensity; set { if (Set(ref _cgDensity, Math.Clamp(value, 0.06m, 1.2m))) CgPreview(); } }
    public decimal CgKTheta { get => _cgKTheta; set { if (Set(ref _cgKTheta, Math.Clamp(value, 0, 50))) CgPreview(); } }
    public decimal CgSeed { get => _cgSeed; set { if (Set(ref _cgSeed, Math.Max(1, Math.Round(value)))) CgPreview(); } }
    public decimal CgPushSteps { get => _cgPush; set => Set(ref _cgPush, Math.Clamp(Math.Round(value), 100, 1e8m)); }
    public decimal CgRunSteps { get => _cgRun; set => Set(ref _cgRun, Math.Clamp(Math.Round(value), 100, 1e10m)); }
    public string CgBoxText => ((double)(_cgChains * _cgBeads / _cgDensity) is var v && v > 0 ? Math.Cbrt(v) : 0).ToString("0.00", CultureInfo.InvariantCulture);
    private string _cgHud = "", _cgInfo = "", _cgError = "";
    public string CgHud { get => _cgHud; private set => Set(ref _cgHud, value); }
    public string CgInfo { get => _cgInfo; private set => Set(ref _cgInfo, value); }
    public string CgError { get => _cgError; private set => Set(ref _cgError, value); }
    public string CgFooter => $"{_cgChains * _cgBeads:N0} beads · L = {CgBoxText} σ (from ρσ³ = {_cgDensity.ToString("0.00", CultureInfo.InvariantCulture)})";
    private CapsDocument? _cgDoc;
    public CapsDocument? CgDoc { get => _cgDoc; private set { var old = _cgDoc; if (Set(ref _cgDoc, value)) { CgViewChanged?.Invoke(); old?.Dispose(); } } }
    public event Action? CgViewChanged;
    private int _cgTicket;

    public string CgOptions() => new JsonObject { ["chains"] = (int)_cgChains, ["beads"] = (int)_cgBeads, ["density"] = (double)_cgDensity, ["k_theta"] = (double)_cgKTheta, ["seed"] = (int)_cgSeed }.ToJsonString();

    public void OpenCg()
    {
        SetModule(44);
        CgPreview();
    }

    /// <summary>The melt as built (instant: random walks), for the preview.</summary>
    public void CgPreview()
    {
        Raise(nameof(CgBoxText));
        Raise(nameof(CgFooter));
        if (!IsCg) return;
        var ticket = ++_cgTicket;
        var opts = CgOptions();
        Task.Run(() =>
        {
            try { var (d, r) = CapsDocument.KgBuild(opts, "Kremer–Grest melt"); return (Doc: (CapsDocument?)d, Report: r, Error: (string?)null); }
            catch (Exception e) { return (Doc: (CapsDocument?)null, Report: "", Error: (string?)e.Message); }
        }).ContinueWith(t => Avalonia.Threading.Dispatcher.UIThread.Post(() =>
        {
            var (doc, rep, err) = t.Result;
            if (ticket != _cgTicket) { doc?.Dispose(); return; }
            CgError = err ?? "";
            if (doc == null) return;
            CgDoc = doc;
            var inv = CultureInfo.InvariantCulture;
            var j = JsonNode.Parse(rep)!;
            CgHud = $"{_cgChains:0} × {_cgBeads:0} beads · bond 0.97 σ";
            CgInfo = $"⟨R²⟩/(N−1) = {j["r2_per_bond"]!.GetValue<double>().ToString("0.00", inv)} σ² per bond · closest pair {j["closest"]!.GetValue<double>().ToString("0.00", inv)} σ — random walks overlap until the push-off opens them";
        }));
    }

    /// <summary>The melt as the Studio document.</summary>
    public void BuildCg()
    {
        try
        {
            var (doc, _) = CapsDocument.KgBuild(CgOptions(), "Kremer–Grest melt");
            Show(doc, $"KG_melt_{_cgChains:0}x{_cgBeads:0}");
            GrownUnsaved = true;
            Status = "Kremer–Grest melt built · export it for LAMMPS from Builders › Coarse-grained (σ = 1 Å in CAPS files)";
            SetModule(8);
        }
        catch (Exception e) { CgError = e.Message; }
    }

    public string ExportCg(string stem)
    {
        if (_cgDoc == null) return "Nothing built";
        _cgDoc.KgLammps(CgOptions(), stem, (double)_cgPush, (double)_cgRun);
        return $"Wrote {System.IO.Path.GetFileName(stem)}.data and {System.IO.Path.GetFileName(stem)}.in · run: lmp -in {System.IO.Path.GetFileName(stem)}.in";
    }
}
