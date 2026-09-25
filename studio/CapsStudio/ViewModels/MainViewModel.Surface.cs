using System.Collections.ObjectModel;
using System.Globalization;
using System.Text.Json;
using System.Text.Json.Nodes;
using CapsStudio.Interop;

namespace CapsStudio.ViewModels;

/// <summary>A bulk crystal of data/crystals/catalogue.json (or a CIF the user opened).</summary>
public sealed record CrystalEntry(string Id, string Name, string File, string Use, int H, int K, int L)
{
    public override string ToString() => Name;
}

/// <summary>A polymer for the film: a library rubber, or the Polymer builder's chain.</summary>
public sealed record FilmPolymer(string Name, string Spec)
{
    public override string ToString() => Name;
}

/// <summary>One layer of the interface stack (top → bottom), as the Surface board shows it.</summary>
public sealed record StackLayer(string Number, string Title, string Height, string Colour);

/// <summary>Surface &amp; interface builder (design/boards/SurfaceBuilder): cleave a crystal along (hkl), choose the
/// termination, passivate, and grow a polymer film on the surface (fibre–rubber interfaces).</summary>
public sealed partial class MainViewModel
{
    public bool IsSurface => _module == 14;

    public ObservableCollection<CrystalEntry> Crystals { get; } = new();
    private List<string> _surfTerms = new();
    /// <summary>The terminations of the plane (replaced as a whole so the list's selection follows).</summary>
    public List<string> SurfTerminations { get => _surfTerms; private set { _surfTerms = value; Raise(); } }
    public ObservableCollection<StackLayer> SurfStack { get; } = new();
    public ObservableCollection<FilmPolymer> FilmPolymers { get; } = new();
    private bool _surfLoaded;

    public void LoadSurface()
    {
        if (_surfLoaded) return;
        _surfLoaded = true;
        try
        {
            var dir = Paths.Crystals;
            if (dir != null)
            {
                using var js = JsonDocument.Parse(File.ReadAllText(Path.Combine(dir, "catalogue.json")));
                foreach (var e in js.RootElement.GetProperty("crystals").EnumerateArray())
                {
                    var hkl = e.TryGetProperty("hkl", out var h) ? h.EnumerateArray().Select(x => x.GetInt32()).ToArray() : [0, 0, 1];
                    Crystals.Add(new CrystalEntry(e.GetProperty("id").GetString()!, e.GetProperty("name").GetString()!, Path.Combine(dir, e.GetProperty("file").GetString()!),
                        e.TryGetProperty("use", out var u) ? u.GetString() ?? "" : "", hkl[0], hkl[1], hkl[2]));
                }
            }
        }
        catch (Exception e) { SurfError = "Cannot read the crystal library: " + e.Message; }
        LoadPolymerLibrary();
        foreach (var p in PolymerLibrary.Where(p => p.Rubber && !p.Copolymer))
            FilmPolymers.Add(new FilmPolymer(p.Name, new JsonObject
            {
                ["units"] = new JsonArray(new JsonObject { ["name"] = p.Name, ["smiles"] = p.Smiles }),
                ["sequence"] = "homopolymer",
            }.ToJsonString()));
        var nr = FilmPolymers.ToList().FindIndex(f => f.Name.Contains("natural rubber", StringComparison.OrdinalIgnoreCase));
        if (FilmPolymers.Count > 0) SurfFilmItem = FilmPolymers[nr >= 0 ? nr : 0];
        _surfCrystal = Crystals.Count > 0 ? 0 : -1;
        if (_surfCrystal >= 0) UseCrystalIndices(Crystals[0]);
        Raise(nameof(SurfCrystal));
        SurfRefresh();
    }

    public void OpenSurface()
    {
        LoadSurface();
        // the Polymer builder's chain, when it has a valid one, is offered first
        if (!PolyHasError && PolyUnits.Count > 0 && (FilmPolymers.Count == 0 || FilmPolymers[0].Name != "Polymer builder chain"))
        {
            var keep = _surfFilmItem;
            FilmPolymers.Insert(0, new FilmPolymer("Polymer builder chain", PolySpecJson()));
            _surfFilmItem = null;
            SurfFilmItem = keep ?? FilmPolymers[0];
        }
        SetModule(14);
        SurfPreview();
    }

    // ---------------------------------------------------------------- cleave
    private int _surfCrystal = -1;
    public int SurfCrystal
    {
        get => _surfCrystal;
        set { if (Set(ref _surfCrystal, value) && value >= 0 && value < Crystals.Count) { UseCrystalIndices(Crystals[value]); SurfRefresh(); SurfPreview(); } }
    }
    private void UseCrystalIndices(CrystalEntry c)
    {
        _surfH = c.H; _surfK = c.K; _surfL = c.L;
        Raise(nameof(SurfH)); Raise(nameof(SurfK)); Raise(nameof(SurfL));
    }
    public string SurfCif => _surfCrystal >= 0 && _surfCrystal < Crystals.Count ? Crystals[_surfCrystal].File : "";
    public string SurfCrystalName => _surfCrystal >= 0 && _surfCrystal < Crystals.Count ? Crystals[_surfCrystal].Name : "crystal";
    public string SurfCrystalUse => _surfCrystal >= 0 && _surfCrystal < Crystals.Count ? Crystals[_surfCrystal].Use : "";

    /// <summary>Adds a CIF file to the list and selects it.</summary>
    public void UseCif(string path)
    {
        LoadSurface();
        var name = Path.GetFileNameWithoutExtension(path);
        Crystals.Add(new CrystalEntry("file:" + path, name + " (file)", path, "opened from a file", 0, 0, 1));
        SurfCrystal = Crystals.Count - 1;
    }

    private decimal _surfH, _surfK, _surfL = 1, _surfLayers = 3, _surfVacuum = 15, _surfNa = 1, _surfNb = 1;
    public decimal SurfH { get => _surfH; set { if (Set(ref _surfH, Math.Clamp(Math.Round(value), -9, 9))) { SurfRefresh(); SurfPreview(); } } }
    public decimal SurfK { get => _surfK; set { if (Set(ref _surfK, Math.Clamp(Math.Round(value), -9, 9))) { SurfRefresh(); SurfPreview(); } } }
    public decimal SurfL { get => _surfL; set { if (Set(ref _surfL, Math.Clamp(Math.Round(value), -9, 9))) { SurfRefresh(); SurfPreview(); } } }
    public decimal SurfLayers { get => _surfLayers; set { if (Set(ref _surfLayers, Math.Clamp(Math.Round(value), 1, 30))) { RaiseStack(); SurfPreview(); } } }
    public decimal SurfVacuum { get => _surfVacuum; set { if (Set(ref _surfVacuum, Math.Clamp(value, 0, 200))) { RaiseStack(); SurfPreview(); } } }
    public decimal SurfNa { get => _surfNa; set { if (Set(ref _surfNa, Math.Clamp(Math.Round(value), 1, 20))) { _surfAutoCell = false; SurfPreview(); } } }
    public decimal SurfNb { get => _surfNb; set { if (Set(ref _surfNb, Math.Clamp(Math.Round(value), 1, 20))) { _surfAutoCell = false; SurfPreview(); } } }
    /// <summary>Until the user sets the supercell, it is chosen so each surface edge is at least 20 Å (room for chains).</summary>
    private bool _surfAutoCell = true;
    private int _surfTermination;
    public int SurfTermination { get => _surfTermination; set { if (value < 0) return; if (Set(ref _surfTermination, value)) { RaiseStack(); SurfPreview(); } } }
    private bool _surfPassivate = true, _surfOrthogonal = true;
    public bool SurfPassivate { get => _surfPassivate; set { if (Set(ref _surfPassivate, value)) SurfPreview(); } }
    public bool SurfOrthogonal { get => _surfOrthogonal; set { if (Set(ref _surfOrthogonal, value)) SurfPreview(); } }
    private double _surfD;
    private string _surfBulk = "", _surfError = "", _surfLog = "", _surfMatch = "";
    public string SurfBulkText { get => _surfBulk; private set => Set(ref _surfBulk, value); }
    public string SurfError { get => _surfError; private set { if (Set(ref _surfError, value)) Raise(nameof(SurfHasError)); } }
    public bool SurfHasError => _surfError.Length > 0;
    public string SurfLog { get => _surfLog; private set => Set(ref _surfLog, value); }
    /// <summary>Lattice matching row of the slab: the rectangular cell and its shear.</summary>
    public string SurfMatchCell { get => _surfMatch; private set => Set(ref _surfMatch, value); }
    private string _surfShear = "0.00 %";
    public string SurfShear { get => _surfShear; private set => Set(ref _surfShear, value); }
    public string SurfTitle => $"{SurfCrystalShort} ({(int)_surfH}{(int)_surfK}{(int)_surfL})";
    private string SurfCrystalShort => SurfCrystalName.Split(' ')[0];
    public string SurfPassivateText => SurfCrystalName.Contains("quartz", StringComparison.OrdinalIgnoreCase) ? "Hydroxylate dangling bonds (Si–OH)" : "Passivate dangling bonds (O–H, M–OH, H)";

    /// <summary>Terminations of the chosen plane (fast: no slab is built).</summary>
    private void SurfRefresh()
    {
        if (SurfCif.Length == 0) return;
        try
        {
            var r = JsonNode.Parse(CapsDocument.SurfaceTerminations(SurfCif, (int)_surfH, (int)_surfK, (int)_surfL))!;
            if ((bool?)r["ok"] != true) { SurfTerminations = new(); SurfError = (string?)r["error"] ?? "cannot cleave"; return; }
            SurfError = "";
            var terms = r["terminations"] as JsonArray ?? [];
            int k = 0;
            var list = new List<string>();
            foreach (var t in terms)
            {
                var top = PolyUnit.Sub((string?)t!["top"] ?? "");
                var bottom = PolyUnit.Sub((string?)t["bottom"] ?? "");
                list.Add(string.Format(CultureInfo.InvariantCulture, "{0}-terminated · {1} of {2} · {3:F1} bonds/nm² cut{4}", top, ++k, terms.Count,
                    (double?)t["bonds_per_nm2"] ?? 0, bottom != top ? " · bottom " + bottom : ""));
            }
            _surfTermination = -1;
            SurfTerminations = list;
            _surfAutoCell = true;
            _surfD = (double?)r["d"] ?? 0;
            var cell = (r["cell"] as JsonArray ?? []).Select(x => (double?)x ?? 0).ToArray();
            SurfBulkText = string.Format(CultureInfo.InvariantCulture, "{0} · {1} atoms per cell · {2:F3} g/cm³ · d({3}{4}{5}) = {6:F3} Å", PolyUnit.Sub((string?)r["formula"] ?? ""),
                (int?)r["atoms"] ?? 0, (double?)r["density"] ?? 0, (int)_surfH, (int)_surfK, (int)_surfL, _surfD);
            _surfTermination = 0;
            Raise(nameof(SurfTermination));
            // the list control clears its selection when its items are replaced: select the first again afterwards
            Avalonia.Threading.Dispatcher.UIThread.Post(() => { if (_surfTermination <= 0) { _surfTermination = 0; Raise(nameof(SurfTermination)); RaiseStack(); } });
        }
        catch (Exception e) { SurfError = e.Message; }
        Raise(nameof(SurfTitle)); Raise(nameof(SurfCrystalUse)); Raise(nameof(SurfPassivateText));
        RaiseStack();
    }

    private string SurfOptions(bool forInterface) => new JsonObject
    {
        ["h"] = (int)_surfH, ["k"] = (int)_surfK, ["l"] = (int)_surfL, ["layers"] = (int)_surfLayers, ["termination"] = Math.Max(0, _surfTermination),
        ["vacuum"] = forInterface ? 10 : (double)_surfVacuum, ["orthogonal"] = forInterface || _surfOrthogonal ? 1 : 0, ["max_strain"] = 0.02,
        ["na"] = (int)_surfNa, ["nb"] = (int)_surfNb, ["passivate"] = _surfPassivate ? 1 : 0,
    }.ToJsonString();

    // ---------------------------------------------------------------- the slab preview
    private CapsDocument? _surfDoc;
    public CapsDocument? SurfDoc { get => _surfDoc; private set { if (Set(ref _surfDoc, value)) SurfViewChanged?.Invoke(); } }
    public event Action? SurfViewChanged;
    private int _surfTicket;
    private string _surfSummary = "";
    public string SurfSummary { get => _surfSummary; private set => Set(ref _surfSummary, value); }

    /// <summary>Builds the slab alone for the preview (milliseconds); the film is grown by Build interface.</summary>
    public void SurfPreview()
    {
        if (!_surfLoaded || SurfCif.Length == 0 || SurfTerminations.Count == 0) return;
        var ticket = ++_surfTicket;
        var cif = SurfCif;
        var opts = SurfOptions(false);
        var name = SurfTitle;
        Task.Run(() =>
        {
            try
            {
                var (d, r) = CapsDocument.SurfaceBuild(cif, opts, name + " slab");
                return (Doc: (CapsDocument?)d, Report: r, Error: (string?)null);
            }
            catch (Exception e) { return (Doc: (CapsDocument?)null, Report: "", Error: (string?)e.Message); }
        }).ContinueWith(t => Avalonia.Threading.Dispatcher.UIThread.Post(() =>
        {
            var (doc, rep, err) = t.Result;
            if (ticket != _surfTicket) { doc?.Dispose(); return; }
            if (err != null || doc == null) { SurfError = err ?? "cannot build the slab"; return; }
            SurfError = "";
            var old = _surfDoc;
            SurfDoc = doc;
            old?.Dispose();
            SurfLog = rep;
            var s = doc.Summary();
            if (_surfAutoCell && s.CellA > 0 && s.CellB > 0)
            {
                var na = (decimal)Math.Max(1, Math.Ceiling(20.0 / (s.CellA / (double)_surfNa)));
                var nb = (decimal)Math.Max(1, Math.Ceiling(20.0 / (s.CellB / (double)_surfNb)));
                if (na != _surfNa || nb != _surfNb)
                {
                    _surfNa = na; _surfNb = nb;
                    Raise(nameof(SurfNa)); Raise(nameof(SurfNb));
                    doc.Dispose();
                    SurfPreview();
                    return;
                }
            }
            SurfSummary = string.Format(CultureInfo.InvariantCulture, "{0:N0} atoms · cell {1:F2} × {2:F2} × {3:F2} Å", s.Atoms, s.CellA, s.CellB, s.CellC);
            SurfMatchCell = $"{(int)_surfNa} × {(int)_surfNb}";
            var shear = rep.Split('\n').FirstOrDefault(l => l.Contains("sheared by"));
            SurfShear = shear != null ? shear.Replace("surface cell sheared by ", "").Replace(" to make it rectangular", "") : "0.00 %";
            RaiseStack();
        }));
    }

    // ---------------------------------------------------------------- film and stack
    private int _surfFilmPolymer;
    public int SurfFilmPolymer
    {
        get => _surfFilmPolymer;
        set { if (value >= 0 && value < FilmPolymers.Count) SurfFilmItem = FilmPolymers[value]; }
    }
    private FilmPolymer? _surfFilmItem;
    /// <summary>The film polymer (bound by item: the list grows when the Polymer builder's chain is offered).</summary>
    public FilmPolymer? SurfFilmItem
    {
        get => _surfFilmItem;
        set
        {
            if (value == null || !Set(ref _surfFilmItem, value)) return;
            _surfFilmPolymer = FilmPolymers.IndexOf(value);
            Raise(nameof(SurfFilmPolymer));
            RaiseStack();
        }
    }
    private bool _surfFilm = true;
    public bool SurfFilm { get => _surfFilm; set { if (Set(ref _surfFilm, value)) { RaiseStack(); Raise(nameof(SurfBuildText)); } } }
    private decimal _filmThickness = 25, _filmDensity = 0.9m, _filmVacuum, _filmGap = 1, _filmDp = 12;
    public decimal FilmThickness { get => _filmThickness; set { if (Set(ref _filmThickness, Math.Clamp(value, 5, 200))) RaiseStack(); } }
    public decimal FilmDensity { get => _filmDensity; set { if (Set(ref _filmDensity, Math.Clamp(value, 0.1m, 2.0m))) RaiseStack(); } }
    public decimal FilmVacuum { get => _filmVacuum; set { if (Set(ref _filmVacuum, Math.Clamp(value, 0, 200))) RaiseStack(); } }
    public decimal FilmGap { get => _filmGap; set { if (Set(ref _filmGap, Math.Clamp(value, 0, 10))) RaiseStack(); } }
    public decimal FilmDp { get => _filmDp; set { if (Set(ref _filmDp, Math.Clamp(Math.Round(value), 2, 500))) RaiseStack(); } }
    public string SurfBuildText => _surfFilm ? "Build interface" : "Build slab";
    public string SurfStackChip => _surfFilm ? $"Interface stack · 2 layers" : "Slab · 1 layer";
    public string SurfStackKind => _surfFilm ? (_filmVacuum > 0 ? "slab + film + vacuum" : "slab + film · periodic") : "slab + vacuum";
    public string FilmName => _surfFilmItem?.Name ?? "polymer";

    private void RaiseStack()
    {
        SurfStack.Clear();
        var term = _surfTermination >= 0 && _surfTermination < SurfTerminations.Count ? SurfTerminations[_surfTermination].Split(" · ")[0] : "";
        var slabH = _surfD * (double)_surfLayers;
        if (_surfFilm)
        {
            if (_filmVacuum > 0) SurfStack.Add(new StackLayer("3", "Vacuum", $"{_filmVacuum:0.0} Å", "#5B8DEF"));
            else SurfStack.Add(new StackLayer("3", "Periodic image of the surface", "—", "#5B8DEF"));
            SurfStack.Add(new StackLayer("2", $"{ShortName(FilmName)} film · amorphous (CAPS Grow)", $"{_filmThickness:0.0} Å", "#B9BEC4"));
        }
        else SurfStack.Add(new StackLayer("2", "Vacuum", $"{_surfVacuum:0.0} Å", "#5B8DEF"));
        SurfStack.Add(new StackLayer("1", $"{SurfTitle} · {term}", string.Format(CultureInfo.InvariantCulture, "{0:F1} Å", slabH), "#D6A45E"));
        Raise(nameof(SurfStackChip)); Raise(nameof(SurfStackKind)); Raise(nameof(FilmName)); Raise(nameof(FilmMatchText)); Raise(nameof(SurfTitle));
    }
    private static string ShortName(string n) => n.Split(" (")[0];
    public string FilmMatchText => _surfFilm ? "fits cell" : "—";

    // ---------------------------------------------------------------- build
    private bool _surfBuilding;
    public bool SurfBuilding { get => _surfBuilding; private set { if (Set(ref _surfBuilding, value)) Raise(nameof(SurfIdle)); } }
    public bool SurfIdle => !_surfBuilding;

    /// <summary>Builds the slab, or grows the film on it, and opens the result as the Studio document.</summary>
    public async Task BuildSurface()
    {
        if (SurfBuilding || SurfCif.Length == 0 || SurfTerminations.Count == 0) return;
        SurfBuilding = true;
        var cif = SurfCif;
        var title = SurfTitle;
        try
        {
            if (_surfAutoCell)
            {
                // the supercell the preview would choose (each edge ≥ 20 Å), in case the preview has not finished
                _surfNa = _surfNb = 1;
                var one = SurfOptions(false);
                var (probe, _) = await Task.Run(() => CapsDocument.SurfaceBuild(cif, one, "probe"));
                var ps = probe.Summary();
                probe.Dispose();
                _surfNa = (decimal)Math.Max(1, Math.Ceiling(20.0 / ps.CellA));
                _surfNb = (decimal)Math.Max(1, Math.Ceiling(20.0 / ps.CellB));
                Raise(nameof(SurfNa)); Raise(nameof(SurfNb));
            }
            if (!_surfFilm)
            {
                var opts = SurfOptions(false);
                var (doc, rep) = await Task.Run(() => CapsDocument.SurfaceBuild(cif, opts, title + " slab"));
                Show(doc, $"{title} slab · {(int)_surfLayers} layers");
                SurfLog = rep;
                Status = "Slab built · " + (rep.Split('\n').FirstOrDefault() ?? "");
            }
            else
            {
                var spec = JsonNode.Parse((_surfFilmItem ?? FilmPolymers[0]).Spec)!.AsObject();
                spec["dp"] = (int)_filmDp;
                var options = new JsonObject
                {
                    ["crystal"] = cif,
                    ["slab"] = JsonNode.Parse(SurfOptions(true)),
                    ["film"] = new JsonObject { ["thickness"] = (double)_filmThickness, ["density"] = (double)_filmDensity, ["gap"] = (double)_filmGap, ["vacuum"] = (double)_filmVacuum },
                }.ToJsonString();
                var name = FilmName;
                var specText = spec.ToJsonString();
                var o = new CapsGrowOpts { Chains = 0, Dp = (int)_filmDp, Tacticity = 0, Seed = 1, Density = 0, ContactScale = 1.0, Curve = 1 };
                var (doc, rep) = await Task.Run(() => CapsDocument.InterfaceBuild(options, specText, o, (d, t, r) =>
                {
                    Avalonia.Threading.Dispatcher.UIThread.Post(() => Status = $"Growing the film · {d} of {t} chains · {r} restarts");
                    return true;
                }, title + " interface"));
                Show(doc, $"{title} + {ShortName(name)} film");
                GrownUnsaved = true;
                SurfLog = rep;
                Status = "Interface built · the surface (molecule 1) is held in place in Relax";
            }
            SetModule(8);
        }
        catch (Exception e) { SurfError = e.Message; Status = "Could not build: " + e.Message; }
        finally { SurfBuilding = false; }
    }
}
