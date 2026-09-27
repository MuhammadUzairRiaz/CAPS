using System.Globalization;
using System.Text.Json.Nodes;
using CapsStudio.Interop;

namespace CapsStudio.ViewModels;

/// <summary>Nanostructure builder (design/boards/NanoBuilder, Nanoparticle): graphene sheets, carbon nanotubes and
/// nanoparticles cut from crystals — rubber fillers — alone or embedded in a grown polymer matrix.</summary>
public sealed partial class MainViewModel
{
    public bool IsNano => _module == 15;
    // silane coupling agents (silica fillers in rubber), grafted as an undoable edit
    public static readonly string[] SilaneChoices = ["TESPT (Si69, tetrasulfide)", "TESPD (Si75, disulfide)", "MPTES (mercaptopropyl)", "APTES (aminopropyl)", "VTES (vinyl)", "OCTEO (octyl, hydrophobic)"];
    private static readonly string[] SilaneIds = ["TESPT", "TESPD", "MPTES", "APTES", "VTES", "OCTEO"];
    private int _silanePick;
    private double _silaneFraction = 0.25, _silaneSpacing = 5;
    public int SilanePick { get => _silanePick; set => Set(ref _silanePick, Math.Clamp(value, 0, SilaneIds.Length - 1)); }
    public decimal? SilaneFractionD { get => (decimal)_silaneFraction; set { _silaneFraction = Math.Clamp((double)(value ?? 0.25m), 0.01, 1); Raise(); } }
    public decimal? SilaneSpacingD { get => (decimal)_silaneSpacing; set { _silaneSpacing = Math.Clamp((double)(value ?? 5m), 2, 30); Raise(); } }
    public void GraftSilane()
    {
        if (_doc == null) return;
        if (RunEdit(new { op = "graft", silane = SilaneIds[_silanePick], fraction = _silaneFraction, min_spacing = _silaneSpacing, seed = 1 }) is { } r)
            Status = (r["what"]?.GetValue<string>() ?? "Grafted") + " · relax before dynamics";
    }

    public static readonly string[] ParticleShapes = ["Sphere", "Cuboctahedron", "Octahedron", "Cube", "Fibre", "Truncated octahedron", "Icosahedron"];

    public void OpenNano()
    {
        LoadSurface();   // the crystal library and the film polymers
        SetModule(15);
        NanoPreview();
    }

    // 0 sheet, 1 nanotube, 2 particle, 3 pore
    private int _nanoKind = 1;
    public int NanoKind
    {
        get => _nanoKind;
        set { if (Set(ref _nanoKind, value)) { if (value == 3) NanoMatrix = false; RaisePore(); NanoPreview(); } }
    }
    public bool NanoIsSheet => _nanoKind == 0;
    public bool NanoIsTube => _nanoKind == 1;
    public bool NanoIsParticle => _nanoKind == 2;

    // tube
    private decimal _tubeN = 10, _tubeM = 10, _tubeLength = 25;
    private bool _nanoPeriodic = true;
    public decimal TubeN { get => _tubeN; set { if (Set(ref _tubeN, Math.Clamp(Math.Round(value), 1, 60))) { if (_tubeM > _tubeN) { _tubeM = _tubeN; Raise(nameof(TubeM)); } RaiseNano(); NanoPreview(); } } }
    public decimal TubeM { get => _tubeM; set { if (Set(ref _tubeM, Math.Clamp(Math.Round(value), 0, _tubeN))) { RaiseNano(); NanoPreview(); } } }
    public decimal TubeLength { get => _tubeLength; set { if (Set(ref _tubeLength, Math.Clamp(value, 3, 500))) NanoPreview(); } }
    // C–C bond length (sheets and tubes) and concentric walls (armchair or zigzag tubes)
    private decimal _nanoCc = 1.42m, _tubeWalls = 1;
    public decimal NanoCc { get => _nanoCc; set { if (Set(ref _nanoCc, Math.Clamp(value, 1.30m, 1.60m))) { RaiseNano(); NanoPreview(); } } }
    public bool TubeMultiWalled
    {
        get => _tubeWalls > 1;
        set { TubeWalls = value ? Math.Max(2, _tubeWalls) : 1; Raise(); }
    }
    public decimal TubeWalls { get => _tubeWalls; set { if (Set(ref _tubeWalls, Math.Clamp(Math.Round(value), 1, 6))) { Raise(nameof(TubeMultiWalled)); RaiseNano(); NanoPreview(); } } }
    /// <summary>"(5,5)@(10,10)@(15,15) · outer d = 20.35 Å", the walls as the core steps them (about 3.4 Å apart).</summary>
    private string TubeWallsTitle()
    {
        var a = Math.Sqrt(3) * (double)_nanoCc;
        var arm = TubeKind == 0;
        var dn = Math.Max(1, (int)Math.Round(2 * 3.4 * Math.PI / (arm ? a * Math.Sqrt(3) : a)));
        var n = (int)_tubeN;
        var names = Enumerable.Range(0, (int)_tubeWalls).Select(k => arm ? $"({n + k * dn},{n + k * dn})" : $"({n + k * dn},0)");
        var nOut = n + ((int)_tubeWalls - 1) * dn;
        var d = arm ? a * Math.Sqrt(3) * nOut / Math.PI : a * nOut / Math.PI;
        return string.Join("@", names) + string.Format(CultureInfo.InvariantCulture, " · outer d = {0:F2} Å", d);
    }
    /// <summary>Concentric walls need armchair or zigzag tubes (their periods along the axis match).</summary>
    public bool TubeWallsAllowed => TubeKind != 2;
    public bool NanoPeriodic { get => _nanoPeriodic; set { if (Set(ref _nanoPeriodic, value)) { Raise(nameof(NanoCap)); RaiseNano(); NanoPreview(); } } }
    public bool NanoCap { get => !_nanoPeriodic; set => NanoPeriodic = !value; }
    /// <summary>0 armchair (m = n), 1 zigzag (m = 0), 2 chiral.</summary>
    public int TubeKind
    {
        get => _tubeM == _tubeN ? 0 : _tubeM == 0 ? 1 : 2;
        set
        {
            if (value == 0) TubeM = _tubeN;
            else if (value == 1) TubeM = 0;
            else if (_tubeM == _tubeN || _tubeM == 0) TubeM = Math.Max(1, Math.Round(_tubeN / 2));
        }
    }
    public bool TubeArmchair { get => TubeKind == 0; set { if (value) TubeKind = 0; } }
    public bool TubeZigzag { get => TubeKind == 1; set { if (value) TubeKind = 1; } }
    public bool TubeChiral { get => TubeKind == 2; set { if (value) TubeKind = 2; } }
    private double[] TubeGeometry()
    {
        double n = (double)_tubeN, m = (double)_tubeM, a = Math.Sqrt(3) * (double)_nanoCc;
        var ch = a * Math.Sqrt(n * n + n * m + m * m);
        var dr = Gcd((int)(2 * m + n), (int)(2 * n + m));
        return [ch / Math.PI, Math.Atan2(Math.Sqrt(3) * m, 2 * n + m) * 180 / Math.PI, Math.Sqrt(3) * ch / dr];
    }
    private static int Gcd(int a, int b) { while (b != 0) (a, b) = (b, a % b); return Math.Abs(a); }
    public string TubeAngleText => TubeGeometry()[1].ToString("F2", CultureInfo.InvariantCulture) + "°";
    public string TubeDiameterText => TubeGeometry()[0].ToString("F2", CultureInfo.InvariantCulture) + " Å";
    public string TubeTText => TubeGeometry()[2].ToString("F2", CultureInfo.InvariantCulture) + " Å";

    // sheet
    private decimal _sheetLx = 20, _sheetLy = 20, _sheetLayers = 1;
    public decimal SheetLx { get => _sheetLx; set { if (Set(ref _sheetLx, Math.Clamp(value, 5, 300))) NanoPreview(); } }
    public decimal SheetLy { get => _sheetLy; set { if (Set(ref _sheetLy, Math.Clamp(value, 5, 300))) NanoPreview(); } }
    public decimal SheetLayers { get => _sheetLayers; set { if (Set(ref _sheetLayers, Math.Clamp(Math.Round(value), 1, 20))) NanoPreview(); } }

    // particle
    private int _particleCrystal, _particleShape;
    private decimal _particleRadius = 12;
    private bool _particleOnAtom = true, _particlePassivate = true;
    public int ParticleCrystal { get => _particleCrystal; set { if (Set(ref _particleCrystal, value)) NanoPreview(); } }
    public int ParticleShape { get => _particleShape; set { if (Set(ref _particleShape, value)) { Raise(nameof(ParticleIsFibre)); RaiseNano(); NanoPreview(); } } }
    public bool ParticleIsFibre => _particleShape == 4;
    private decimal _fibreLength = 22;
    public decimal FibreLength { get => _fibreLength; set { if (Set(ref _fibreLength, Math.Clamp(value, 3, 500))) NanoPreview(); } }
    public decimal ParticleRadius { get => _particleRadius; set { if (Set(ref _particleRadius, Math.Clamp(value, 3, 60))) NanoPreview(); } }
    public int ParticleCentre { get => _particleOnAtom ? 0 : 1; set { _particleOnAtom = value == 0; Raise(); NanoPreview(); } }
    public bool ParticlePassivate { get => _particlePassivate; set { if (Set(ref _particlePassivate, value)) NanoPreview(); } }
    public static readonly string[] ParticleCentres = ["On an atom", "On the cell centre"];

    // matrix
    private bool _nanoMatrix;
    private decimal _matrixChains = 10, _matrixDp = 20, _matrixDensity = 0.9m;
    private FilmPolymer? _matrixPolymer;
    public bool NanoMatrix { get => _nanoMatrix; set { if (Set(ref _nanoMatrix, value)) RaiseNano(); } }
    public decimal MatrixChains { get => _matrixChains; set => Set(ref _matrixChains, Math.Clamp(Math.Round(value), 1, 2000)); }
    public decimal MatrixDp { get => _matrixDp; set => Set(ref _matrixDp, Math.Clamp(Math.Round(value), 2, 2000)); }
    public decimal MatrixDensity { get => _matrixDensity; set => Set(ref _matrixDensity, Math.Clamp(value, 0.1m, 2.0m)); }
    public FilmPolymer? MatrixPolymer { get => _matrixPolymer ?? _surfFilmItem; set => Set(ref _matrixPolymer, value); }

    public string NanoTitle => _nanoKind switch
    {
        3 => PoreTitle,
        0 => $"Graphene · {_sheetLayers} layer{(_sheetLayers > 1 ? "s" : "")}",
        1 => TubeWallsAllowed && _tubeWalls > 1 ? TubeWallsTitle()
                                                : $"({_tubeN},{_tubeM}) {(TubeKind == 0 ? "armchair" : TubeKind == 1 ? "zigzag" : "chiral")} · d = {TubeGeometry()[0]:F2} Å",
        _ => $"{(_particleCrystal < Crystals.Count ? Crystals[_particleCrystal].Name : "crystal")} {ParticleShapes[_particleShape].ToLowerInvariant()} · r = {_particleRadius:0.#} Å",
    };
    public string NanoAxisText => _nanoKind == 3 ? (_poreType == 0 ? (_poreVacuum ? "vacuum above the walls" : "periodic in x, y, z") : "periodic in x, y, z") : _nanoKind == 1 ? (_nanoPeriodic ? "periodic along z" : "capped ends") : _nanoKind == 0 ? (_nanoPeriodic ? "periodic in the plane" : "flake")
        : _particleShape == 4 ? "fibre · periodic along z" : "cut from the crystal";
    public string NanoBuildText => _nanoMatrix ? "Build composite" : _nanoKind switch { 0 => "Build sheet", 1 => "Build nanotube", 3 => "Build pore", _ => "Build particle" };
    public string NanoBuildIcon => _nanoKind switch { 2 => "atom", 0 => "hex", 3 => "ring", _ => "layers" };

    private void RaiseNano()
    {
        foreach (var n in new[] { nameof(NanoIsSheet), nameof(NanoIsTube), nameof(NanoIsParticle), nameof(TubeKind), nameof(TubeArmchair), nameof(TubeZigzag), nameof(TubeChiral),
                                  nameof(TubeAngleText), nameof(TubeDiameterText), nameof(TubeTText), nameof(NanoTitle), nameof(NanoAxisText), nameof(NanoBuildText), nameof(NanoBuildIcon), nameof(TubeWallsAllowed) })
            Raise(n);
    }

    private string NanoOptions()
    {
        if (_nanoKind == 3) return PoreOptions();
        var o = new JsonObject { ["kind"] = _nanoKind switch { 0 => "sheet", 1 => "tube", _ => "particle" }, ["periodic"] = _nanoPeriodic ? 1 : 0 };
        switch (_nanoKind)
        {
            case 0: o["lx"] = (double)_sheetLx; o["ly"] = (double)_sheetLy; o["layers"] = (int)_sheetLayers; o["cc"] = (double)_nanoCc; break;
            case 1:
                o["n"] = (int)_tubeN; o["m"] = (int)_tubeM; o["length"] = (double)_tubeLength; o["cc"] = (double)_nanoCc;
                o["walls"] = TubeWallsAllowed ? (int)_tubeWalls : 1;
                break;
            default:
                o["crystal"] = _particleCrystal < Crystals.Count ? Crystals[_particleCrystal].File : "";
                o["shape"] = ParticleShapes[_particleShape].ToLowerInvariant();
                o["radius"] = (double)_particleRadius;
                o["on_atom"] = _particleOnAtom ? 1 : 0;
                o["passivate"] = _particlePassivate ? 1 : 0;
                o["length"] = (double)_fibreLength;
                break;
        }
        return o.ToJsonString();
    }

    private CapsDocument? _nanoDoc;
    public CapsDocument? NanoDoc { get => _nanoDoc; private set { if (Set(ref _nanoDoc, value)) NanoViewChanged?.Invoke(); } }
    public event Action? NanoViewChanged;
    private int _nanoTicket;
    private string _nanoSummary = "", _nanoLog = "", _nanoError = "";
    public string NanoSummary { get => _nanoSummary; private set => Set(ref _nanoSummary, value); }
    public string NanoLog { get => _nanoLog; private set => Set(ref _nanoLog, value); }
    public string NanoError { get => _nanoError; private set { if (Set(ref _nanoError, value)) Raise(nameof(NanoHasError)); } }
    public bool NanoHasError => _nanoError.Length > 0;
    private bool _nanoBuilding;
    public bool NanoBuilding { get => _nanoBuilding; private set { if (Set(ref _nanoBuilding, value)) Raise(nameof(NanoIdle)); } }
    public bool NanoIdle => !_nanoBuilding;

    /// <summary>The filler alone, for the preview (milliseconds).</summary>
    public void NanoPreview()
    {
        if (!_surfLoaded) return;
        RaiseNano();
        var ticket = ++_nanoTicket;
        var opts = NanoOptions();
        var title = NanoTitle;
        var pore = _nanoKind == 3;
        Task.Run(() =>
        {
            try
            {
                var (d, r) = pore ? CapsDocument.PoreBuild(opts, title) : CapsDocument.NanoBuild(opts, title);
                return (Doc: (CapsDocument?)d, Report: r, Error: (string?)null);
            }
            catch (Exception e) { return (Doc: (CapsDocument?)null, Report: "", Error: (string?)e.Message); }
        }).ContinueWith(t => Avalonia.Threading.Dispatcher.UIThread.Post(() =>
        {
            var (doc, rep, err) = t.Result;
            if (ticket != _nanoTicket) { doc?.Dispose(); return; }
            if (err != null || doc == null) { NanoError = err ?? "cannot build"; return; }
            NanoError = "";
            var old = _nanoDoc;
            NanoDoc = doc;
            old?.Dispose();
            NanoLog = pore ? PoreReportText(rep) : rep;
            var s = doc.Summary();
            NanoSummary = string.Format(CultureInfo.InvariantCulture, "{0:N0} atoms · {1:N0} bonds", s.Atoms, s.Bonds);
        }));
    }

    /// <summary>Builds the filler (or the composite) and opens it as the Studio document.</summary>
    public async Task BuildNano()
    {
        if (NanoBuilding) return;
        NanoBuilding = true;
        var opts = NanoOptions();
        var title = NanoTitle;
        try
        {
            if (_nanoKind == 3)
            {
                var (doc, rep) = await Task.Run(() => CapsDocument.PoreBuild(opts, title));
                Show(doc, title);
                NanoLog = PoreReportText(rep);
                HoldPick = 1;
                HoldOn = true;   // the walls stay where they are
                RelaxCompress = false;
                Status = "Pore built · the walls (molecule 1) are held in Relax and Dynamics";
            }
            else if (!_nanoMatrix)
            {
                var (doc, rep) = await Task.Run(() => CapsDocument.NanoBuild(opts, title));
                Show(doc, title);
                NanoLog = rep;
                Status = "Built · " + (rep.Split('\n').FirstOrDefault() ?? "");
            }
            else
            {
                var poly = MatrixPolymer ?? FilmPolymers.FirstOrDefault() ?? throw new InvalidOperationException("no polymer for the matrix");
                var o = JsonNode.Parse(opts)!.AsObject();
                o["matrix"] = new JsonObject { ["chains"] = (int)_matrixChains, ["density"] = (double)_matrixDensity };
                var spec = JsonNode.Parse(poly.Spec)!.AsObject();
                spec["dp"] = (int)_matrixDp;
                var g = new CapsGrowOpts { Chains = (int)_matrixChains, Dp = (int)_matrixDp, Tacticity = 0, Seed = 1, Density = 0, ContactScale = 1.0, Curve = 1 };
                var (optsText, specText) = (o.ToJsonString(), spec.ToJsonString());
                var (doc, rep) = await Task.Run(() => CapsDocument.NanoEmbed(optsText, specText, g, (d, t, r) =>
                {
                    Avalonia.Threading.Dispatcher.UIThread.Post(() => Status = $"Growing the matrix · {d} of {t} chains · {r} restarts");
                    return true;
                }, title + " composite"));
                Show(doc, $"{title} in {poly.Name.Split(" (")[0]}");
                GrownUnsaved = true;
                RelaxCompress = false;   // compression would scale the filler with the matrix
                NanoLog = rep;
                Status = "Composite built · the filler (molecule 1) is held in place in Relax and Dynamics";
            }
            SetModule(8);
        }
        catch (Exception e) { NanoError = e.Message; Status = "Could not build: " + e.Message; }
        finally { NanoBuilding = false; }
    }
}
