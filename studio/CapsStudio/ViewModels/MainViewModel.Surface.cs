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

/// <summary>One layer of the interface stack (top → bottom), as the Surface board shows it; Extra ≥ 0: an added layer
/// (removable).</summary>
public sealed record StackLayer(string Number, string Title, string Height, string Colour, int Extra = -1, LayerOps? Ops = null)
{
    public bool Removable => Extra >= 0;
}

/// <summary>What is done to an added layer before it goes on the stack (B9): turned upside down, shifted in the plane.</summary>
public sealed class LayerOps : ObservableObject
{
    private bool _flip;
    private string _shift = "";
    public bool Flip { get => _flip; set => Set(ref _flip, value); }
    /// <summary>"dx dy" in Å.</summary>
    public string Shift { get => _shift; set => Set(ref _shift, value ?? ""); }
    public (double X, double Y) ShiftXY()
    {
        var p = _shift.Split([' ', ',', ';'], StringSplitOptions.RemoveEmptyEntries);
        double V(int k) => k < p.Length && double.TryParse(p[k], NumberStyles.Float, CultureInfo.InvariantCulture, out var v) ? v : 0;
        return (V(0), V(1));
    }
}

/// <summary>A row of the lattice-matching table: a layer, its repeats, and the strain it takes to fit the slab's cell.</summary>
public sealed record MatchRow(string Layer, string Supercell, string StrainA, string StrainB, bool TooStrained);

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
    // molecular crystals (polymer crystals, cellulose, organic solids): the surfaces made of whole molecules
    private bool _surfWhole;
    public bool SurfWhole { get => _surfWhole; set { if (Set(ref _surfWhole, value)) SurfPreview(); } }
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
        ["vacuum"] = forInterface ? 10 : (double)_surfVacuum, ["orthogonal"] = forInterface || _surfOrthogonal ? 1 : 0, ["max_strain"] = (double)_surfMaxStrain / 100,
        ["na"] = (int)_surfNa, ["nb"] = (int)_surfNb, ["passivate"] = _surfPassivate ? 1 : 0, ["whole_molecules"] = _surfWhole ? 1 : 0,
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
            var s = doc.Summary();
            if (_surfAutoCell && s.CellA > 0 && s.CellB > 0)
            {
                var na = (decimal)Math.Max(1, Math.Ceiling(20.0 / (s.CellA / (double)_surfNa)));
                var nb = (decimal)Math.Max(1, Math.Ceiling(20.0 / (s.CellB / (double)_surfNb)));
                if (na != _surfNa || nb != _surfNb)
                {   // a larger surface cell first: this slab is not shown (the preview keeps the last one until then)
                    _surfNa = na; _surfNb = nb;
                    Raise(nameof(SurfNa)); Raise(nameof(SurfNb));
                    doc.Dispose();
                    SurfPreview();
                    return;
                }
            }
            var old = _surfDoc;
            SurfDoc = doc;
            old?.Dispose();
            SurfLog = rep;
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
    // a grafted brush in place of the free film (core polymer.hpp build_brush)
    public static readonly string[] FilmKinds = ["Free film grown against the surface", "Brush: chains grafted to surface sites"];
    private int _filmKind;
    private decimal _brushSigma = 0.3m, _brushSpacing = 4m;
    private string _brushSite = "O";
    public int FilmKind { get => _filmKind; set { if (Set(ref _filmKind, Math.Clamp(value, 0, 1))) Raise(nameof(FilmIsBrush)); } }
    public bool FilmIsBrush => _filmKind == 1;
    public decimal BrushSigma { get => _brushSigma; set => Set(ref _brushSigma, Math.Clamp(value, 0.01m, 10m)); }
    public decimal BrushSpacing { get => _brushSpacing; set => Set(ref _brushSpacing, Math.Clamp(value, 1m, 30m)); }
    public string BrushSite { get => _brushSite; set => Set(ref _brushSite, (value ?? "O").Trim()); }
    public decimal FilmThickness { get => _filmThickness; set { if (Set(ref _filmThickness, Math.Clamp(value, 5, 200))) RaiseStack(); } }
    public decimal FilmDensity { get => _filmDensity; set { if (Set(ref _filmDensity, Math.Clamp(value, 0.1m, 2.0m))) RaiseStack(); } }
    public decimal FilmVacuum { get => _filmVacuum; set { if (Set(ref _filmVacuum, Math.Clamp(value, 0, 200))) RaiseStack(); } }
    public decimal FilmGap { get => _filmGap; set { if (Set(ref _filmGap, Math.Clamp(value, 0, 10))) RaiseStack(); } }
    public decimal FilmDp { get => _filmDp; set { if (Set(ref _filmDp, Math.Clamp(Math.Round(value), 2, 500))) RaiseStack(); } }
    public string SurfBuildText => _surfFilm || _surfExtra.Count > 0 ? "Build interface" : "Build slab";
    public string SurfStackChip => (_surfFilm ? 2 : 1) + _surfExtra.Count is var n && n > 1 ? $"Interface stack · {n} layers" : "Slab · 1 layer";
    public string SurfStackKind => _surfFilm ? (_filmVacuum > 0 ? "slab + film + vacuum" : "slab + film · periodic") : "slab + vacuum";
    public string FilmName => _surfFilmItem?.Name ?? "polymer";

    // lattice matching: the slab sets the lateral cell; a grown film is made in it; an added layer is repeated to come
    // closest and stretched to fit (the same choice stack_layers makes, first layer as it is)
    private decimal _surfMaxStrain = 2;
    public decimal SurfMaxStrain { get => _surfMaxStrain; set { if (Set(ref _surfMaxStrain, Math.Clamp(value, 0, 20))) { RaiseStack(); SurfPreview(); } } }
    public ObservableCollection<MatchRow> SurfMatchRows { get; } = new();
    private string _surfMatchNote = "";
    public string SurfMatchNote { get => _surfMatchNote; private set => Set(ref _surfMatchNote, value); }

    private void RefreshMatchRows()
    {
        SurfMatchRows.Clear();
        var inv = CultureInfo.InvariantCulture;
        var max = (double)_surfMaxStrain / 100;
        var shearText = _surfShear;
        var shear = double.TryParse(shearText.Replace("%", "").Trim(), NumberStyles.Float, inv, out var sh) ? Math.Abs(sh) / 100 : 0;
        SurfMatchRows.Add(new MatchRow(SurfTitle + " · reference", SurfMatchCell.Length > 0 ? SurfMatchCell : "—", shearText, "0.00 %", shear > max + 1e-12));
        if (_surfFilm) SurfMatchRows.Add(new MatchRow($"{ShortName(FilmName)} film", "grown in it", "0.00 %", "0.00 %", false));
        double A = 0, B = 0;
        if (_surfDoc != null) { var s = _surfDoc.Summary(); A = s.CellA; B = s.CellB; }
        var worst = 0.0;
        foreach (var e in _surfExtra)
        {
            if (A <= 0 || B <= 0) { SurfMatchRows.Add(new MatchRow(e.Name, "—", "—", "—", false)); continue; }
            var c = e.Doc.Summary();
            static int Best(double L, double T) { var best = 1; for (var k = 1; k <= 6; ++k) if (Math.Abs(T / (k * L) - 1) < Math.Abs(T / (best * L) - 1)) best = k; return best; }
            var na = Best(c.CellA, A); var nb = Best(c.CellB, B);
            double sa = A / (na * c.CellA) - 1, sb = B / (nb * c.CellB) - 1;
            worst = Math.Max(worst, Math.Max(Math.Abs(sa), Math.Abs(sb)));
            SurfMatchRows.Add(new MatchRow(e.Name, $"{na} × {nb}", (100 * sa).ToString("+0.00;-0.00", inv) + " %", (100 * sb).ToString("+0.00;-0.00", inv) + " %",
                                           Math.Abs(sa) > max || Math.Abs(sb) > max));
        }
        SurfMatchNote = worst > max ? $"A layer is strained by more than {_surfMaxStrain:0.#} %: an amorphous polymer relaxes it away; for a crystal pick another slab size (the supercell above) or another plane"
                      : _surfExtra.Count > 0 ? "Strains before stacking; the slab may be repeated too when that fits the layers better (the report after Build lists the repeats)" : "";
    }

    private void RaiseStack()
    {
        RefreshMatchRows();
        SurfStack.Clear();
        var term = _surfTermination >= 0 && _surfTermination < SurfTerminations.Count ? SurfTerminations[_surfTermination].Split(" · ")[0] : "";
        var slabH = _surfD * (double)_surfLayers;
        var n = 1 + (_surfFilm ? 1 : 0) + _surfExtra.Count;   // layers under the top one
        var vacuum = _surfFilm ? _filmVacuum : _surfVacuum;
        if (vacuum > 0) SurfStack.Add(new StackLayer($"{n + 1}", "Vacuum", $"{vacuum:0.0} Å", "#5B8DEF"));
        else SurfStack.Add(new StackLayer($"{n + 1}", "Periodic image of the surface", "—", "#5B8DEF"));
        for (var k = _surfExtra.Count - 1; k >= 0; --k)
        {
            var e = _surfExtra[k];
            SurfStack.Add(new StackLayer($"{(_surfFilm ? 3 : 2) + k}", e.Name, e.Height, e.Colour, k, e.Ops));
        }
        if (_surfFilm) SurfStack.Add(new StackLayer("2", $"{ShortName(FilmName)} film · amorphous (CAPS Grow)", $"{_filmThickness:0.0} Å", "#B9BEC4"));
        SurfStack.Add(new StackLayer("1", $"{SurfTitle} · {term}", string.Format(CultureInfo.InvariantCulture, "{0:F1} Å", slabH), "#D6A45E"));
        Raise(nameof(SurfStackChip)); Raise(nameof(SurfStackKind)); Raise(nameof(FilmName)); Raise(nameof(FilmMatchText)); Raise(nameof(SurfTitle));
        Raise(nameof(SurfHasExtra));
    }

    // ---------------------------------------------------------------- added layers (design/boards/SurfaceBuilder "Add layer")
    private sealed record ExtraLayer(string Name, CapsDocument Doc, string Height, string Colour)
    {
        public LayerOps Ops { get; } = new();
    }
    private bool _surfShareStrain;
    /// <summary>The lateral cell the mean of the layers' (match "average"), every layer strained a little; else the
    /// first layer's (the slab's) and the others stretched to it.</summary>
    public bool SurfShareStrain { get => _surfShareStrain; set => Set(ref _surfShareStrain, value); }
    private readonly List<ExtraLayer> _surfExtra = new();
    public bool SurfHasExtra => _surfExtra.Count > 0;
    private decimal _surfStackGap = 2.5m;
    /// <summary>Å between one layer's top atoms and the next one's lowest (added layers).</summary>
    public decimal SurfStackGap { get => _surfStackGap; set => Set(ref _surfStackGap, Math.Clamp(value, 0, 20)); }

    private void AddExtra(string name, CapsDocument doc, string colour)
    {
        var s = doc.Summary();
        if (!(s.CellA > 0 && s.CellB > 0 && s.CellC > 0)) { doc.Dispose(); SurfError = $"{name}: a layer needs a periodic cell"; return; }
        SurfError = "";
        _surfExtra.Add(new ExtraLayer(name, doc, string.Format(CultureInfo.InvariantCulture, "{0:F1} Å cell", s.CellC), colour));
        RaiseStack();
        Raise(nameof(SurfBuildText));
    }

    /// <summary>From CAPS Grow cell: the open structure (a grown polymer cell) as the layer above the surface, in place
    /// of growing a film there.</summary>
    public void AddLayerFromOpen()
    {
        if (Document == null) { SurfError = "Open or grow a polymer cell first (Grow), then add it here"; return; }
        var name = string.IsNullOrWhiteSpace(Title) ? "grown cell" : Title;
        AddExtra($"{ShortName(name)} · from the open structure", Document.FrameCopy(name), "#B9BEC4");
        if (_surfFilm) { SurfFilm = false; Status = "The open cell goes on the surface in place of a grown film"; }
    }

    /// <summary>Add layer › this slab again (a film between two surfaces).</summary>
    public void AddSlabLayer()
    {
        if (SurfCif.Length == 0 || SurfTerminations.Count == 0) return;
        try
        {
            var (doc, _) = CapsDocument.SurfaceBuild(SurfCif, SurfOptions(false), SurfTitle + " slab");
            AddExtra($"{SurfTitle} slab", doc, "#D6A45E");
        }
        catch (Exception e) { SurfError = e.Message; }
    }

    /// <summary>Add layer › a structure file with a rectangular periodic cell.</summary>
    public void AddLayerFile(string path)
    {
        try { AddExtra(Path.GetFileName(path), CapsDocument.Open(path), "#9C8FD6"); }
        catch (Exception e) { SurfError = $"{Path.GetFileName(path)}: {e.Message}"; }
    }

    public void RemoveSurfLayer(int extra)
    {
        if (extra < 0 || extra >= _surfExtra.Count) return;
        _surfExtra[extra].Doc.Dispose();
        _surfExtra.RemoveAt(extra);
        RaiseStack();
        Raise(nameof(SurfBuildText));
    }

    /// <summary>The built base (slab or slab + film) with the added layers stacked on it; null when there are none.</summary>
    private (CapsDocument Doc, string Log) StackOnto(CapsDocument baseDoc, string title)
    {
        var docs = new List<CapsDocument> { baseDoc };
        docs.AddRange(_surfExtra.Select(e => e.Doc));
        var names = new JsonArray(new[] { (JsonNode)title }.Concat(_surfExtra.Select(e => (JsonNode)e.Name)).ToArray());
        var flips = new JsonArray(new[] { (JsonNode)false }.Concat(_surfExtra.Select(e => (JsonNode)e.Ops.Flip)).ToArray());
        var shifts = new JsonArray(new[] { (JsonNode)new JsonArray(0.0, 0.0) }.Concat(_surfExtra.Select(e => { var (x, y) = e.Ops.ShiftXY(); return (JsonNode)new JsonArray(x, y); })).ToArray());
        var opts = new JsonObject { ["names"] = names, ["gap"] = (double)_surfStackGap, ["vacuum"] = (double)(_surfFilm ? _filmVacuum : _surfVacuum),
                                    ["match"] = _surfShareStrain ? "average" : "both", ["flips"] = flips, ["shifts"] = shifts }.ToJsonString();
        var (doc, rep) = CapsDocument.Stack(docs, opts, title + " stack");
        var r = JsonNode.Parse(rep);
        if (doc == null) throw new InvalidOperationException((string?)r?["error"] ?? "cannot stack the layers");
        return (doc, string.Join("\n", (r?["notes"] as JsonArray ?? []).Select(x => (string?)x ?? "")));
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
                SurfLog = rep;
                if (_surfExtra.Count > 0)
                {
                    var (stacked, log) = await Task.Run(() => { try { return StackOnto(doc, title); } finally { doc.Dispose(); } });
                    SurfLog = log;
                    Show(stacked, $"{title} + {string.Join(" + ", _surfExtra.Select(e => ShortName(e.Name.Split(" · ")[0])))}");
                    KeepPageSettings("surface");   // Edit brings the page back as it was for this structure
                    GrownUnsaved = true;
                    RelaxCompress = false;
                    Status = "Layers stacked · " + (SurfLog.Split('\n').FirstOrDefault() ?? "");
                }
                else
                {
                    Show(doc, $"{title} slab · {(int)_surfLayers} layers");
                    KeepPageSettings("surface");   // Edit brings the page back as it was for this structure
                    Status = "Slab built · " + (rep.Split('\n').FirstOrDefault() ?? "");
                }
            }
            else
            {
                var spec = JsonNode.Parse((_surfFilmItem ?? FilmPolymers[0]).Spec)!.AsObject();
                spec["dp"] = (int)_filmDp;
                var options = new JsonObject
                {
                    ["crystal"] = cif,
                    ["slab"] = JsonNode.Parse(SurfOptions(true)),
                    ["film"] = new JsonObject { ["thickness"] = (double)_filmThickness, ["density"] = (double)_filmDensity, ["gap"] = (double)_filmGap,
                                                ["vacuum"] = _filmKind == 1 ? Math.Max(10.0, (double)_filmVacuum) : (double)_filmVacuum },
                };
                if (_filmKind == 1) options["brush"] = new JsonObject { ["density"] = (double)_brushSigma, ["site"] = _brushSite, ["min_spacing"] = (double)_brushSpacing };
                var optionsText = options.ToJsonString();
                var name = FilmName;
                var specText = spec.ToJsonString();
                var o = new CapsGrowOpts { Chains = 0, Dp = (int)_filmDp, Tacticity = 0, Seed = (ulong)FilmSeed.Take(), Density = 0, ContactScale = 1.0, Curve = 1 };
                var (doc, rep) = await Task.Run(() => CapsDocument.InterfaceBuild(optionsText, specText, o, (d, t, r) =>
                {
                    Avalonia.Threading.Dispatcher.UIThread.Post(() => Status = $"Growing the film · {d} of {t} chains · {r} restarts");
                    return true;
                }, title + " interface"));
                SurfLog = rep;
                if (_surfExtra.Count > 0)
                {
                    var built = doc;
                    (doc, var log) = await Task.Run(() => { try { return StackOnto(built, $"{title} + {ShortName(name)} film"); } finally { built.Dispose(); } });
                    SurfLog = rep + "\n" + log;
                }
                Show(doc, $"{title} + {ShortName(name)} {(_filmKind == 1 ? "brush" : "film")}" + (_surfExtra.Count > 0 ? $" + {_surfExtra.Count} layer{(_surfExtra.Count == 1 ? "" : "s")}" : ""));
                KeepPageSettings("surface");   // Edit brings the page back as it was for this structure
                GrownUnsaved = true;
                RelaxCompress = false;   // compression would scale the crystal with the film
                Status = "Interface built · the surface (molecule 1) is held in place in Relax";
            }
            SetModule(8);
        }
        catch (Exception e) { SurfError = e.Message; Status = "Could not build: " + e.Message; }
        finally { SurfBuilding = false; }
    }

    // ---------------------------------------------------------------- held molecule (Relax, Dynamics, Equilibrate)
    private long _holdMol;
    private decimal _holdPick = 1;
    /// <summary>Hold molecule HoldPick in place (an interface's surface is molecule 1 and starts held).</summary>
    public bool HoldOn
    {
        get => _holdMol > 0;
        set
        {
            _holdMol = value ? (long)_holdPick : 0;
            try { Document?.SetHeldMolecule(_holdMol); } catch { }
            Raise(); Raise(nameof(HoldText));
            Status = value ? $"Molecule {_holdPick} is held in place in Relax, Dynamics and Equilibrate" : "No atoms held";
        }
    }
    public decimal HoldPick
    {
        get => _holdPick;
        set
        {
            if (!Set(ref _holdPick, Math.Max(1, Math.Round(value)))) return;
            if (_holdMol > 0) HoldOn = true;
        }
    }
    public string HoldText => _holdMol > 0 ? $"molecule {_holdMol} held" : "all atoms move";
    private void SyncHeld()
    {
        _holdMol = Document?.HeldMolecule() ?? 0;
        if (_holdMol > 0) { _holdPick = _holdMol; Raise(nameof(HoldPick)); }
        Raise(nameof(HoldOn)); Raise(nameof(HoldText));
    }
}
