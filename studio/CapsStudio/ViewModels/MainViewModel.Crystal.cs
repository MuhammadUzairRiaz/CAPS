using System.Collections.ObjectModel;
using System.Globalization;
using System.Text.Json.Nodes;
using CapsStudio.Interop;

namespace CapsStudio.ViewModels;

/// <summary>A space-group setting as the crystal builder lists it.</summary>
public sealed record SpaceGroupItem(string Key, int Number, string Hm, string System)
{
    /// <summary>Hermann–Mauguin without spaces, screw axes subscripted, the 1s of monoclinic full symbols dropped.</summary>
    public string Pretty
    {
        get
        {
            var sym = Hm.Split(':')[0];
            var tokens = sym.Split(' ', StringSplitOptions.RemoveEmptyEntries).ToList();
            if (Number is >= 3 and <= 15 && tokens.Count == 4) tokens = tokens.Where((t, i) => i == 0 || t != "1").ToList();
            const string sub = "₀₁₂₃₄₅₆₇₈₉";
            return string.Concat(tokens.Select(t => t.Length >= 2 && char.IsDigit(t[0]) && char.IsDigit(t[1]) ? t[0] + sub[t[1] - '0'].ToString() + t[2..] : t));
        }
    }
    /// <summary>"No. 227:2" where the setting matters (origin choice, axes, a non-standard monoclinic cell).</summary>
    public string NumberText
    {
        get
        {
            var c = Key.IndexOf(':');
            if (c < 0) return $"No. {Number}";
            var code = Key[(c + 1)..];
            if (Number is >= 3 and <= 15) return code is "b" or "b1" ? $"No. {Number}" : $"No. {Key}";
            return code is "1" or "2" or "h" or "r" ? $"No. {Key}" : $"No. {Number}";
        }
    }
    public string Tip => Setting.Length > 0 ? $"{Hm} · {Setting} · Hall setting {Key}" : $"{Hm} · setting {Key}";
    /// <summary>The setting where the symbol alone does not say it: origin choice, axes, the monoclinic unique axis.</summary>
    public string Setting
    {
        get
        {
            var c = Key.IndexOf(':');
            if (c < 0) return "";
            var s = Key[(c + 1)..];
            if (Number is >= 3 and <= 15) return $"unique {s.TrimStart('-')[0]}";   // P2₁/c with b or c unique read the same
            return s switch { "1" => "origin 1", "2" => "origin 2", "h" => "hexagonal axes", "r" => "rhombohedral axes", _ => "" };   // orthorhombic settings: the symbol says it
        }
    }
    public string Title => Pretty;
    public override string ToString() => Title;
}

/// <summary>A site of the asymmetric unit: label, element, fractional x y z, occupancy.</summary>
public sealed class CrystalSiteRow : ObservableObject
{
    private readonly Action _changed;
    private string _label, _element, _x, _y, _z, _occ;
    public CrystalSiteRow(Action changed, string label, string element, double x, double y, double z, double occ = 1)
    {
        _changed = changed;
        _label = label; _element = element;
        _x = F(x); _y = F(y); _z = F(z); _occ = occ.ToString("0.00", CultureInfo.InvariantCulture);
    }
    private static string F(double v) => v.ToString("0.0000", CultureInfo.InvariantCulture);
    public string Label { get => _label; set { if (Set(ref _label, value)) _changed(); } }
    public string Element { get => _element; set { if (Set(ref _element, value)) _changed(); } }
    public string X { get => _x; set { if (Set(ref _x, value)) _changed(); } }
    public string Y { get => _y; set { if (Set(ref _y, value)) _changed(); } }
    public string Z { get => _z; set { if (Set(ref _z, value)) _changed(); } }
    public string Occ { get => _occ; set { if (Set(ref _occ, value)) _changed(); } }
    private string _mult = "";
    /// <summary>Atoms this site gives per cell (its multiplicity), after a build.</summary>
    public string Multiplicity { get => _mult; set => Set(ref _mult, value); }
    public static double Parse(string s)
    {
        s = s.Trim();
        var slash = s.IndexOf('/');
        if (slash > 0 && double.TryParse(s[..slash], NumberStyles.Float, CultureInfo.InvariantCulture, out var p) &&
            double.TryParse(s[(slash + 1)..], NumberStyles.Float, CultureInfo.InvariantCulture, out var q) && q != 0) return p / q;
        return double.TryParse(s, NumberStyles.Float, CultureInfo.InvariantCulture, out var v) ? v : double.NaN;
    }
}

/// <summary>Crystal builder (design/boards/CrystalBuilder): a space group (the 530 settings of the 230 groups, from
/// their Hall symbols), a lattice and an asymmetric unit expanded into a cell; symmetry found, sites snapped onto
/// special positions, primitive cells and supercells; CIF import.</summary>
public sealed partial class MainViewModel
{
    public bool IsCrystal => _module == 29;

    private List<SpaceGroupItem>? _allGroups;
    public ObservableCollection<SpaceGroupItem> CrystalGroups { get; } = new();
    public ObservableCollection<CrystalSiteRow> CrystalSites { get; } = new();
    public static readonly string[] CrystalTolerances = ["0.001 Å", "0.01 Å", "0.05 Å", "0.1 Å"];
    private static readonly double[] CrystalTolValues = [0.001, 0.01, 0.05, 0.1];
    // groups offered before anything is typed: polymer crystals and the fillers and fibres of rubber composites
    private static readonly string[] CrystalFavourites = ["Pnam", "Pnma", "P 1 21/c 1", "F m -3 m", "F d -3 m:2", "P 32 2 1", "P 42/m n m", "P 63 m c", "P 63/m m c", "C 1 2/m 1", "P -1", "I m -3 m"];

    public void OpenCrystal()
    {
        LoadSurface();   // the crystal library
        if (_allGroups == null)
        {
            _allGroups = new();
            try
            {
                foreach (var g in JsonNode.Parse(CapsDocument.SpaceGroups())!.AsArray())
                    _allGroups.Add(new SpaceGroupItem(g!["key"]!.GetValue<string>(), (int)g["number"]!.GetValue<double>(), g["hm"]!.GetValue<string>(), g["system"]!.GetValue<string>()));
            }
            catch (Exception e) { CrystalError = "Cannot list the space groups: " + e.Message; }
            if (CrystalSites.Count == 0) LoadPolyethylene();
        }
        SetModule(29);
        CrystalFilter();
        CrystalPreview();
    }

    /// <summary>Orthorhombic polyethylene (Bunn 1939): carbon from the paper, the CH₂ hydrogens placed at C–H 1.09 Å, H–C–H 107°.</summary>
    private void LoadPolyethylene()
    {
        _quiet = true;
        CrystalGroup = FindGroup("Pnam");
        _cA = 7.40m; _cB = 4.93m; _cC = 2.534m; _cAlpha = _cBeta = _cGamma = 90;
        CrystalSites.Clear();
        CrystalSites.Add(new CrystalSiteRow(CrystalPreview, "C1", "C", 0.0380, 0.0650, 0.25));
        CrystalSites.Add(new CrystalSiteRow(CrystalPreview, "H1", "H", 0.1848, 0.0466, 0.25));
        CrystalSites.Add(new CrystalSiteRow(CrystalPreview, "H2", "H", 0.0068, 0.2811, 0.25));
        _cSuper = "2 × 3 × 2";
        CrystalTitle = "PE_crystal";
        CrystalCite = "Bunn, Trans. Faraday Soc. 35, 482 (1939) — orthorhombic polyethylene; H placed at C–H 1.09 Å, H–C–H 107°";
        _quiet = false;
        RaiseLattice();
        Raise(nameof(CrystalSupercell));
    }

    private SpaceGroupItem? FindGroup(string q)
    {
        if (_allGroups == null) return null;
        var squash = new string(q.Where(c => !char.IsWhiteSpace(c)).ToArray()).ToLowerInvariant();
        return _allGroups.FirstOrDefault(g => g.Key == q) ?? _allGroups.FirstOrDefault(g => new string(g.Hm.Where(c => !char.IsWhiteSpace(c)).ToArray()).ToLowerInvariant() == squash)
               ?? _allGroups.FirstOrDefault(g => new string(g.Hm.Split(':')[0].Where(c => !char.IsWhiteSpace(c)).ToArray()).ToLowerInvariant() == squash);
    }

    // ---------------------------------------------------------------- space group

    private string _sgQuery = "";
    public string CrystalQuery { get => _sgQuery; set { if (Set(ref _sgQuery, value)) CrystalFilter(); } }

    private void CrystalFilter()
    {
        if (_allGroups == null) return;
        var q = _sgQuery.Trim();
        IEnumerable<SpaceGroupItem> hits;
        if (q.Length == 0)
            hits = CrystalFavourites.Select(FindGroup).OfType<SpaceGroupItem>();
        else
        {
            static string Norm(string s) => new string(s.Where(c => !char.IsWhiteSpace(c) && c != '_').ToArray()).ToLowerInvariant()
                .Replace("₁", "1").Replace("₂", "2").Replace("₃", "3").Replace("₄", "4").Replace("₅", "5");
            var nq = Norm(q);
            hits = int.TryParse(q, out var num)
                ? _allGroups.Where(g => g.Number == num)
                : _allGroups.Where(g => Norm(g.Hm).StartsWith(nq) || Norm(g.Pretty).StartsWith(nq) || g.System.StartsWith(q, StringComparison.OrdinalIgnoreCase))
                    .OrderBy(g => Norm(g.Hm) == nq || Norm(g.Pretty) == nq ? 0 : 1);
        }
        var list = hits.Take(12).ToList();
        if (_cGroup != null && !list.Contains(_cGroup) && q.Length == 0) list.Insert(0, _cGroup);
        CrystalGroups.Clear();
        foreach (var g in list) CrystalGroups.Add(g);
        Raise(nameof(CrystalGroupPick));
        CrystalQueryNote = list.Count == 0 ? $"No space group matches “{q}”" : q.Length == 0 ? "" : $"{list.Count}{(list.Count == 12 ? "+" : "")} settings";
    }
    private string _queryNote = "";
    public string CrystalQueryNote { get => _queryNote; private set => Set(ref _queryNote, value); }

    private SpaceGroupItem? _cGroup;
    public SpaceGroupItem? CrystalGroup
    {
        get => _cGroup;
        set
        {
            if (value == null || !Set(ref _cGroup, value)) return;
            ApplyLatticeRules();
            RaiseLattice();
            Raise(nameof(CrystalGroupPick));
            Raise(nameof(CrystalCanPrimitive));
            if (!CrystalCanPrimitive) CrystalPrimitive = false;
            CrystalPreview();
        }
    }
    /// <summary>The list's selection: the chosen group when it is listed (a rebuilt list does not clear the choice).</summary>
    public SpaceGroupItem? CrystalGroupPick { get => _cGroup != null && CrystalGroups.Contains(_cGroup) ? _cGroup : null; set { if (value != null) CrystalGroup = value; } }
    public string CrystalSystemText => _cGroup == null ? "—" : CultureInfo.InvariantCulture.TextInfo.ToTitleCase(_cGroup.System);

    // ---------------------------------------------------------------- lattice

    private decimal _cA = 5, _cB = 5, _cC = 5, _cAlpha = 90, _cBeta = 90, _cGamma = 90;
    private bool _quiet;
    public decimal CrystalA { get => _cA; set { if (Set(ref _cA, Math.Max(0.5m, value))) { ApplyLatticeRules(); RaiseLattice(); CrystalPreview(); } } }
    public decimal CrystalB { get => _cB; set { if (Set(ref _cB, Math.Max(0.5m, value))) { ApplyLatticeRules(); RaiseLattice(); CrystalPreview(); } } }
    public decimal CrystalC { get => _cC; set { if (Set(ref _cC, Math.Max(0.5m, value))) { ApplyLatticeRules(); RaiseLattice(); CrystalPreview(); } } }
    public decimal CrystalAlpha { get => _cAlpha; set { if (Set(ref _cAlpha, Math.Clamp(value, 10, 170))) { ApplyLatticeRules(); RaiseLattice(); CrystalPreview(); } } }
    public decimal CrystalBeta { get => _cBeta; set { if (Set(ref _cBeta, Math.Clamp(value, 10, 170))) { ApplyLatticeRules(); RaiseLattice(); CrystalPreview(); } } }
    public decimal CrystalGamma { get => _cGamma; set { if (Set(ref _cGamma, Math.Clamp(value, 10, 170))) { ApplyLatticeRules(); RaiseLattice(); CrystalPreview(); } } }

    private char UniqueAxis => _cGroup != null && _cGroup.Number is >= 3 and <= 15 && _cGroup.Key.Contains(':') ? _cGroup.Key[(_cGroup.Key.IndexOf(':') + 1)..].TrimStart('-').FirstOrDefault('b') : 'b';
    private bool Rhombohedral => _cGroup != null && _cGroup.Key.EndsWith(":r");
    private string Sys => _cGroup?.System ?? "triclinic";
    // which lengths and angles the crystal system leaves free
    public bool CrystalBFree => Sys is "triclinic" or "monoclinic" or "orthorhombic";
    public bool CrystalCFree => Sys != "cubic" && !Rhombohedral;
    public bool CrystalAlphaFree => Sys == "triclinic" || Rhombohedral || (Sys == "monoclinic" && UniqueAxis == 'a');
    public bool CrystalBetaFree => Sys == "triclinic" || (Sys == "monoclinic" && UniqueAxis == 'b');
    public bool CrystalGammaFree => Sys == "triclinic" || (Sys == "monoclinic" && UniqueAxis == 'c');
    public string CrystalLatticeRule => Sys switch
    {
        "cubic" => "Cubic: a = b = c, α = β = γ = 90°",
        "tetragonal" => "Tetragonal: a = b, α = β = γ = 90°",
        "trigonal" or "hexagonal" when Rhombohedral => "Rhombohedral axes: a = b = c, α = β = γ",
        "trigonal" or "hexagonal" => "Hexagonal axes: a = b, α = β = 90°, γ = 120°",
        "orthorhombic" => "Orthorhombic: α = β = γ = 90°",
        "monoclinic" => $"Monoclinic, unique axis {UniqueAxis}: the other two angles 90°",
        _ => "Triclinic: every length and angle free",
    };

    /// <summary>Lengths and angles the crystal system fixes follow the free ones.</summary>
    private void ApplyLatticeRules()
    {
        if (_cGroup == null) return;
        if (!CrystalBFree) _cB = _cA;
        if (!CrystalCFree) _cC = _cA;
        if (Rhombohedral) { _cBeta = _cGamma = _cAlpha; return; }
        if (Sys is "trigonal" or "hexagonal") { _cAlpha = _cBeta = 90; _cGamma = 120; return; }
        if (!CrystalAlphaFree) _cAlpha = 90;
        if (!CrystalBetaFree) _cBeta = 90;
        if (!CrystalGammaFree) _cGamma = 90;
    }

    private void RaiseLattice()
    {
        foreach (var n in new[] { nameof(CrystalA), nameof(CrystalB), nameof(CrystalC), nameof(CrystalAlpha), nameof(CrystalBeta), nameof(CrystalGamma), nameof(CrystalBFree),
                                  nameof(CrystalCFree), nameof(CrystalAlphaFree), nameof(CrystalBetaFree), nameof(CrystalGammaFree), nameof(CrystalLatticeRule), nameof(CrystalSystemText),
                                  nameof(CrystalHudGroup) })
            Raise(n);
    }

    // ---------------------------------------------------------------- asymmetric unit, tools

    public void AddCrystalSite()
    {
        var el = CrystalSites.LastOrDefault()?.Element ?? "C";
        var k = CrystalSites.Count(s => s.Element == el) + 1;
        CrystalSites.Add(new CrystalSiteRow(CrystalPreview, $"{el}{k}", el, 0, 0, 0));
        CrystalPreview();
    }

    public void RemoveCrystalSite(CrystalSiteRow row)
    {
        CrystalSites.Remove(row);
        CrystalPreview();
    }

    private string _cSuper = "1 × 1 × 1";
    public string CrystalSupercell { get => _cSuper; set { if (Set(ref _cSuper, value)) { Raise(nameof(CrystalHudSuper)); CrystalPreview(); } } }
    private int[] SupercellCounts()
    {
        var parts = _cSuper.Split(['×', 'x', 'X', '*', ',', ' '], StringSplitOptions.RemoveEmptyEntries);
        var n = parts.Select(p => int.TryParse(p, out var v) ? Math.Clamp(v, 1, 50) : 1).ToList();
        while (n.Count < 3) n.Add(n.Count == 0 ? 1 : n[^1]);
        return n.Take(3).ToArray();
    }

    private int _cTol = 1;
    public int CrystalTolerance { get => _cTol; set { if (Set(ref _cTol, Math.Clamp(value, 0, 3))) CrystalPreview(); } }
    private bool _cPrim;
    public bool CrystalPrimitive { get => _cPrim; set { if (Set(ref _cPrim, value)) { Raise(nameof(CrystalHudGroup)); CrystalPreview(); } } }
    /// <summary>A centred lattice has a smaller primitive cell (rhombohedral axes are already primitive).</summary>
    public bool CrystalCanPrimitive => _cGroup != null && _cGroup.Hm.Length > 0 && _cGroup.Hm[0] != 'P' && !Rhombohedral;

    private string _cTitle = "crystal", _cCite = "";
    public string CrystalTitle { get => _cTitle; set => Set(ref _cTitle, value); }
    public string CrystalCite { get => _cCite; private set { if (Set(ref _cCite, value)) Raise(nameof(CrystalHasCite)); } }
    public bool CrystalHasCite => _cCite.Length > 0;

    /// <summary>The spec for the core (caps_crystal_*); null with CrystalError set when a field does not parse.</summary>
    private string? CrystalSpec(bool withSupercell = true)
    {
        if (_cGroup == null) { CrystalError = "Choose a space group"; return null; }
        var sites = new JsonArray();
        foreach (var s in CrystalSites)
        {
            double x = CrystalSiteRow.Parse(s.X), y = CrystalSiteRow.Parse(s.Y), z = CrystalSiteRow.Parse(s.Z), o = CrystalSiteRow.Parse(s.Occ);
            if (double.IsNaN(x) || double.IsNaN(y) || double.IsNaN(z)) { CrystalError = $"Site {s.Label}: x, y and z are fractions (0.25 or 1/4)"; return null; }
            if (string.IsNullOrWhiteSpace(s.Element)) { CrystalError = $"Site {s.Label}: no element"; return null; }
            sites.Add(new JsonObject { ["label"] = s.Label, ["element"] = s.Element.Trim(), ["x"] = x, ["y"] = y, ["z"] = z, ["occupancy"] = double.IsNaN(o) ? 1 : o });
        }
        var sc = withSupercell ? SupercellCounts() : [1, 1, 1];
        return new JsonObject
        {
            ["space_group"] = _cGroup.Key, ["a"] = (double)_cA, ["b"] = (double)_cB, ["c"] = (double)_cC,
            ["alpha"] = (double)_cAlpha, ["beta"] = (double)_cBeta, ["gamma"] = (double)_cGamma,
            ["sites"] = sites, ["supercell"] = new JsonArray(sc[0], sc[1], sc[2]), ["tolerance"] = CrystalTolValues[_cTol],
            ["primitive"] = _cPrim && CrystalCanPrimitive, ["title"] = _cTitle,
        }.ToJsonString();
    }

    /// <summary>Loads a spec the core returned (find symmetry, symmetrize, CIF import) into the fields.</summary>
    private void LoadCrystalSpec(JsonNode spec)
    {
        _quiet = true;
        var g = FindGroup(spec["space_group"]!.GetValue<string>());
        if (g != null) { _cGroup = g; Raise(nameof(CrystalGroup)); Raise(nameof(CrystalCanPrimitive)); }
        decimal D(string k) => Math.Round((decimal)spec[k]!.GetValue<double>(), 4);
        _cA = D("a"); _cB = D("b"); _cC = D("c"); _cAlpha = D("alpha"); _cBeta = D("beta"); _cGamma = D("gamma");
        CrystalSites.Clear();
        foreach (var s in spec["sites"]!.AsArray())
            CrystalSites.Add(new CrystalSiteRow(CrystalPreview, s!["label"]!.GetValue<string>(), s["element"]!.GetValue<string>(), s["x"]!.GetValue<double>(), s["y"]!.GetValue<double>(),
                                                s["z"]!.GetValue<double>(), s["occupancy"]?.GetValue<double>() ?? 1));
        _quiet = false;
        RaiseLattice();
        CrystalFilter();
    }

    /// <summary>Snaps every site onto its special position (the mean of its images that land within 0.3 Å of it).</summary>
    public void CrystalApplySymmetry()
    {
        var spec = CrystalSpec();
        if (spec == null) return;
        var r = JsonNode.Parse(CapsDocument.CrystalSymmetrize(spec, 0.3))!;
        if (r["ok"]?.GetValue<bool>() != true) { CrystalError = r["error"]?.GetValue<string>() ?? "cannot apply the symmetry"; return; }
        var moved = (int)r["moved"]!.GetValue<double>();
        var keep = _cSuper;
        LoadCrystalSpec(r["spec"]!);
        _cSuper = keep;
        CrystalLog = moved == 0 ? "Every site already sits on its position: nothing moved" : $"{moved} site{(moved == 1 ? "" : "s")} moved onto {(moved == 1 ? "its" : "their")} special position{(moved == 1 ? "" : "s")}";
        CrystalPreview();
    }

    /// <summary>The highest-symmetry setting that maps the current crystal onto itself; its asymmetric unit replaces the sites.</summary>
    public async Task CrystalFindSymmetry(string? cif = null)
    {
        string? spec = null;
        if (cif == null) { spec = CrystalSpec(false); if (spec == null) return; }
        var tol = Math.Max(0.01, CrystalTolValues[_cTol]);
        CrystalBusy = true;
        try
        {
            var text = await Task.Run(() => CapsDocument.CrystalFindSymmetry(spec, cif, cif != null ? Math.Max(tol, 0.05) : tol));
            var r = JsonNode.Parse(text)!;
            if (r["ok"]?.GetValue<bool>() != true) { CrystalError = r["error"]?.GetValue<string>() ?? "cannot find the symmetry"; return; }
            var was = _cGroup?.Pretty;
            LoadCrystalSpec(r["spec"]!);
            var found = $"{_cGroup?.Pretty} (No. {(int)r["number"]!.GetValue<double>()}, {r["system"]!.GetValue<string>()}) · {(int)r["operations"]!.GetValue<double>()} operations";
            if (cif != null)
            {
                CrystalTitle = Path.GetFileNameWithoutExtension(cif);
                CrystalCite = $"Imported from {Path.GetFileName(cif)}: {(int)r["atoms"]!.GetValue<double>()} atoms per cell reduced to {CrystalSites.Count} site{(CrystalSites.Count == 1 ? "" : "s")}";
                CrystalLog = "Found " + found;
            }
            else CrystalLog = was == _cGroup?.Pretty ? $"Symmetry confirmed: {found}" : $"Found {found} (was {was})";
            CrystalPreview();
        }
        catch (Exception e) { CrystalError = e.Message; }
        finally { CrystalBusy = false; }
    }

    public Task ImportCrystalCif(string path) => CrystalFindSymmetry(path);

    /// <summary>The crystal library (data/crystals) for the Import menu.</summary>
    public IEnumerable<CrystalEntry> CrystalLibrary => Crystals;

    // ---------------------------------------------------------------- preview, build

    private CapsDocument? _crystalDoc;
    public CapsDocument? CrystalDoc { get => _crystalDoc; private set { if (Set(ref _crystalDoc, value)) CrystalViewChanged?.Invoke(); } }
    public event Action? CrystalViewChanged;
    private int _crystalTicket;
    private string _crystalLog = "", _crystalError = "", _crystalAtoms = "", _crystalVolume = "", _crystalExpand = "";
    public string CrystalLog { get => _crystalLog; private set => Set(ref _crystalLog, value); }
    public string CrystalError { get => _crystalError; private set { if (Set(ref _crystalError, value)) Raise(nameof(CrystalHasError)); } }
    public bool CrystalHasError => _crystalError.Length > 0;
    public string CrystalAtomsText { get => _crystalAtoms; private set => Set(ref _crystalAtoms, value); }
    public string CrystalVolumeText { get => _crystalVolume; private set => Set(ref _crystalVolume, value); }
    /// <summary>"1 site → 4 atoms / cell".</summary>
    public string CrystalExpandText { get => _crystalExpand; private set => Set(ref _crystalExpand, value); }
    private bool _crystalBusy;
    public bool CrystalBusy { get => _crystalBusy; private set { if (Set(ref _crystalBusy, value)) Raise(nameof(CrystalIdle)); } }
    public bool CrystalIdle => !_crystalBusy;
    public string CrystalHudGroup => _cGroup == null ? "No space group" : $"{CrystalSystemText} · {_cGroup.Pretty}{(_cPrim && CrystalCanPrimitive ? " · primitive" : "")}";
    public string CrystalHudSuper { get { var n = SupercellCounts(); return $"supercell {n[0]} × {n[1]} × {n[2]}"; } }

    /// <summary>Builds the crystal in the background for the preview; the newest request wins.</summary>
    public void CrystalPreview()
    {
        if (_quiet || _module != 29) return;
        Raise(nameof(CrystalHudGroup));
        Raise(nameof(CrystalHudSuper));
        var spec = CrystalSpec();
        if (spec == null) return;
        var ticket = ++_crystalTicket;
        var title = _cTitle;
        Task.Run(() =>
        {
            try
            {
                var info = JsonNode.Parse(CapsDocument.CrystalInfo(spec))!;
                if (info["ok"]?.GetValue<bool>() != true) return (Doc: (CapsDocument?)null, Info: info, Report: "", Error: info["error"]?.GetValue<string>() ?? "cannot build");
                var (d, r) = CapsDocument.CrystalBuild(spec, title);
                return (Doc: (CapsDocument?)d, Info: info, Report: r, Error: (string?)null);
            }
            catch (Exception e) { return (Doc: (CapsDocument?)null, Info: (JsonNode)new JsonObject(), Report: "", Error: (string?)e.Message); }
        }).ContinueWith(t => Avalonia.Threading.Dispatcher.UIThread.Post(() =>
        {
            var (doc, info, rep, err) = t.Result;
            if (ticket != _crystalTicket) { doc?.Dispose(); return; }
            if (err != null || doc == null) { CrystalError = err ?? "cannot build"; return; }
            CrystalError = "";
            var old = _crystalDoc;
            CrystalDoc = doc;
            old?.Dispose();
            var mult = info["multiplicity"]!.AsArray().Select(m => (int)m!.GetValue<double>()).ToList();
            for (var k = 0; k < CrystalSites.Count && k < mult.Count; k++) CrystalSites[k].Multiplicity = mult[k].ToString(CultureInfo.InvariantCulture);
            var perCell = (int)info["atoms_per_cell"]!.GetValue<double>();
            var n = CrystalSites.Count;
            CrystalExpandText = $"{n} site{(n == 1 ? "" : "s")} → {perCell} atoms / cell";
            var s = doc.Summary();
            CrystalAtomsText = string.Format(CultureInfo.InvariantCulture, "{0:N0} atoms · {1}", s.Atoms, info["formula"]!.GetValue<string>());
            CrystalVolumeText = string.Format(CultureInfo.InvariantCulture, "V cell {0:0.00} Å³ · {1:0.000} g/cm³", info["volume"]!.GetValue<double>(), info["density"]!.GetValue<double>());
            var notes = info["notes"]!.AsArray().Select(x => x!.GetValue<string>()).ToList();
            CrystalLog = notes.Count > 0 ? string.Join("\n", notes) : string.Format(CultureInfo.InvariantCulture, "{0} · {1} operations · Hall {2}",
                info["hm"]!.GetValue<string>(), (int)info["operations"]!.GetValue<double>(), info["hall"]!.GetValue<string>());
        }));
    }

    /// <summary>Builds the crystal (with its supercell) and opens it as the Studio document.</summary>
    public async Task BuildCrystal()
    {
        var spec = CrystalSpec();
        if (spec == null || CrystalBusy) return;
        CrystalBusy = true;
        var title = _cTitle;
        try
        {
            var (doc, rep) = await Task.Run(() => CapsDocument.CrystalBuild(spec, title));
            Show(doc, title);
            KeepPageSettings("crystal");   // Edit brings the page back as it was for this structure
            GrownUnsaved = true;
            Status = "Built · " + (rep.Split('\n').FirstOrDefault() ?? "");
            SetModule(8);
        }
        catch (Exception e) { CrystalError = e.Message; Status = "Could not build: " + e.Message; }
        finally { CrystalBusy = false; }
    }

    /// <summary>The open structure's atoms as the asymmetric unit (B4: a crystal of the current molecule): Cartesian
    /// positions about their centre put into this cell's fractional coordinates at ¼ ¼ ¼, away from the common special
    /// positions; the space group places its copies. The cell and group are the user's (from the published structure).</summary>
    public void SitesFromOpenMolecule()
    {
        if (_doc == null) { CrystalError = "Open or build the molecule first (the structure in the Studio becomes the asymmetric unit)"; return; }
        var n = (int)_doc.Summary().Atoms;
        if (n == 0 || n > 2000) { CrystalError = n == 0 ? "The open structure has no atoms" : $"{n} atoms: the asymmetric unit is one molecule (or a few)"; return; }
        var at = Enumerable.Range(0, n).Select(i => _doc.Atom(i)).ToArray();
        double cx = at.Average(a => a.X), cy = at.Average(a => a.Y), cz = at.Average(a => a.Z);
        // the cell as cell_parameters builds it: a along x, b in the xy plane
        double A = (double)_cA, B = (double)_cB, Cc = (double)_cC, d2r = Math.PI / 180;
        double ca = Math.Cos((double)_cAlpha * d2r), cb = Math.Cos((double)_cBeta * d2r), cg = Math.Cos((double)_cGamma * d2r), sg = Math.Sin((double)_cGamma * d2r);
        double cyc = (ca - cb * cg) / sg, czc = Math.Sqrt(Math.Max(1e-12, 1 - cb * cb - cyc * cyc));
        double[,] m = { { A, B * cg, Cc * cb }, { 0, B * sg, Cc * cyc }, { 0, 0, Cc * czc } };   // columns a, b, c
        // upper triangular: back substitution for the fractional coordinates
        (double, double, double) Frac(double x, double y, double z)
        {
            var fz = z / m[2, 2];
            var fy = (y - m[1, 2] * fz) / m[1, 1];
            var fx = (x - m[0, 1] * fy - m[0, 2] * fz) / m[0, 0];
            return (fx, fy, fz);
        }
        _quiet = true;
        CrystalSites.Clear();
        var count = new Dictionary<string, int>();
        foreach (var a in at)
        {
            var el = string.IsNullOrWhiteSpace(a.ElementSymbol) ? "C" : a.ElementSymbol.Trim();
            count[el] = count.GetValueOrDefault(el) + 1;
            var (fx, fy, fz) = Frac(a.X - cx, a.Y - cy, a.Z - cz);
            CrystalSites.Add(new CrystalSiteRow(CrystalPreview, el + count[el], el, fx + 0.25, fy + 0.25, fz + 0.25));
        }
        CrystalTitle = Title.Replace(" (unsaved)", "").Replace(' ', '_') + "_crystal";
        CrystalCite = $"asymmetric unit: the open molecule ({n} atoms), centred at ¼ ¼ ¼ of this cell — give the cell and space group of the published structure";
        _quiet = false;
        CrystalPreview();
    }
}
