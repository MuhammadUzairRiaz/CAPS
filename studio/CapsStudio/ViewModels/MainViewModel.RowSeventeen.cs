using System;
using System.Collections.ObjectModel;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Text.Json.Nodes;
using System.Threading.Tasks;
using CapsStudio.Interop;

namespace CapsStudio.ViewModels;

public sealed class SpeciesRow : ObservableObject
{
    public required string Name { get; init; }
    public required string Smiles { get; init; }
    public double Mass { get; init; }
    private decimal _rho;
    private string _nExact = "—", _n = "—", _achieved = "—";
    public decimal Rho { get => _rho; set { if (Set(ref _rho, value)) Changed?.Invoke(); } }
    public string MassText => Mass.ToString("0.000", CultureInfo.InvariantCulture);
    public string NExact { get => _nExact; set => Set(ref _nExact, value); }
    public string N { get => _n; set => Set(ref _n, value); }
    public string Achieved { get => _achieved; set => Set(ref _achieved, value); }
    public int Count;
    public Action? Changed;
}
public sealed record SasaGroupRow(string Part, string Area, string Share);
public sealed record SasaConvRow(string Points, string Area, string Delta);
public sealed record CellPresetRow(string Name, string Volume, string Z, string RhoCalc, string RhoLit, double A, double B, double C, double Al, double Be, double Ga, int ZCount, double FormulaMass);

/// <summary>Row 17 of the design (ChainStats, DensityCalc, SurfaceArea, CellEditor): chain size and stiffness, molecule
/// counts from a box and a density, solvent-accessible areas, and the unit cell edited in place.</summary>
public partial class MainViewModel
{
    // ---------------------------------------------------------------- Analyze › Chains (design/boards/ChainStats)
    public bool IsChainStats => _module == 55;
    public ResultCell CsBonds { get; } = new("Bonds per chain, n", "backbone bonds");
    public ResultCell CsRee { get; } = new("⟨R²ₑₑ⟩½", "first to last backbone atom");
    public ResultCell CsRg { get; } = new("⟨R²g⟩½", "mass-weighted, whole molecules");
    public ResultCell CsRatio { get; } = new("⟨R²ₑₑ⟩ / ⟨R²g⟩", "6 for an ideal Gaussian chain");
    public ResultCell CsCn { get; } = new("C_n at n = N", "whole backbones");
    public ResultCell CsCinf { get; } = new("C∞", "extrapolated: C_n = C∞(1 − a/n)");
    public (double X, double Y)[] CsCnCurve { get; private set; } = [];
    public (double X, double Y)[] CsReeHist { get; private set; } = [];
    private string _csStatus = "";
    public string CsStatus { get => _csStatus; private set => Set(ref _csStatus, value); }
    public event Action? CsChanged;

    public void OpenChainStats() { SetModule(55); CsChanged?.Invoke(); }

    public async Task RunChainStats()
    {
        if (_doc == null || Analyze.Working) return;
        var inv = CultureInfo.InvariantCulture;
        CsStatus = "Chains…";
        await RunChips("ree", "rg", "cn");
        var ree = Analyze.Results.FirstOrDefault(r => r.Id == "ree");
        var rg = Analyze.Results.FirstOrDefault(r => r.Id == "rg");
        var cn = Analyze.Results.FirstOrDefault(r => r.Id == "cn");
        string V(ResultCard? c, string f) => c != null && double.IsFinite(c.Value) ? c.Value.ToString(f, inv) + (double.IsFinite(c.Error) && c.Error > 0 ? " ± " + c.Error.ToString(f, inv) : "") : "—";
        double X(ResultCard? c, string k) => c == null ? double.NaN : c.Extra.Where(e => e.Key.StartsWith(k, StringComparison.Ordinal)).Select(e => e.Value).DefaultIfEmpty(double.NaN).First();
        CsBonds.Value = double.IsFinite(X(cn, "backbone bonds")) ? X(cn, "backbone bonds").ToString("0", inv) : "—";
        CsRee.Value = V(ree, "0.0") + (ree != null ? " Å" : "");
        CsRg.Value = V(rg, "0.0") + (rg != null ? " Å" : "");
        CsRatio.Value = ree != null && rg != null && rg.Value > 0 ? (ree.Value * ree.Value / (rg.Value * rg.Value)).ToString("0.00", inv) : "—";
        CsCn.Value = double.IsFinite(X(cn, "C_N")) ? X(cn, "C_N").ToString("0.00", inv) : "—";
        CsCinf.Value = V(cn, "0.0");
        CsCnCurve = cn == null ? [] : Analyze.Curves.Where(c => c.Property == cn.Name).Select(c => c.X.Zip(c.Y).ToArray()).FirstOrDefault() ?? [];
        CsReeHist = ree == null ? [] : Analyze.Curves.Where(c => c.Property == ree.Name).Select(c => c.X.Zip(c.Y).ToArray()).FirstOrDefault() ?? [];
        CsStatus = Analyze.Log + (cn?.Notes.Length > 0 ? " · " + cn.Notes[0] : "");
        CsChanged?.Invoke();
    }

    // ---------------------------------------------------------------- Pack › Density calculator (design/boards/DensityCalc)
    public bool IsDensityCalc => _module == 56;
    public const double Avogadro = 6.02214076e23;
    private decimal _dcA = 40, _dcB = 40, _dcC = 40, _dcTarget = 0.95m, _dcN = 1000;
    private int _dcSolvent = 1, _dcBoxSpecies;
    public decimal DcA { get => _dcA; set { if (Set(ref _dcA, Math.Clamp(value, 1, 10000))) Recalculate(); } }
    public decimal DcB { get => _dcB; set { if (Set(ref _dcB, Math.Clamp(value, 1, 10000))) Recalculate(); } }
    public decimal DcC { get => _dcC; set { if (Set(ref _dcC, Math.Clamp(value, 1, 10000))) Recalculate(); } }
    public decimal DcTarget { get => _dcTarget; set { if (Set(ref _dcTarget, Math.Clamp(value, 0.01m, 30m))) Recalculate(); } }
    public int DcSolvent { get => _dcSolvent; set { if (value >= 0 && Set(ref _dcSolvent, value)) Recalculate(); } }   // −1: a list being replaced
    public decimal DcN { get => _dcN; set { if (Set(ref _dcN, Math.Clamp(value, 1, 1e9m))) Recalculate(); } }
    public int DcBoxSpecies { get => _dcBoxSpecies; set { if (value >= 0 && Set(ref _dcBoxSpecies, value)) Recalculate(); } }
    public ObservableCollection<SpeciesRow> DcSpecies { get; } = new();
    private string _dcVolume = "", _dcPolymer = "", _dcSolventN = "—", _dcFraction = "—", _dcEdge = "—", _dcNewSmiles = "", _dcError = "";
    public string DcVolume { get => _dcVolume; private set => Set(ref _dcVolume, value); }
    public string DcPolymer { get => _dcPolymer; private set => Set(ref _dcPolymer, value); }
    public string DcSolventN { get => _dcSolventN; private set => Set(ref _dcSolventN, value); }
    public string DcFraction { get => _dcFraction; private set => Set(ref _dcFraction, value); }
    public string DcEdge { get => _dcEdge; private set => Set(ref _dcEdge, value); }
    public string DcNewSmiles { get => _dcNewSmiles; set => Set(ref _dcNewSmiles, value ?? ""); }
    public string DcError { get => _dcError; private set => Set(ref _dcError, value); }
    private string[] _dcNames = [];
    /// <summary>The species' names; a new list only when they change, with the two choices put back (a replaced list clears a
    /// ComboBox's selection).</summary>
    public string[] DcSpeciesNames => _dcNames;
    private void RaiseSpeciesNames()
    {
        var names = DcSpecies.Select(s => s.Name).ToArray();
        if (names.SequenceEqual(_dcNames)) return;
        _dcNames = names;
        Raise(nameof(DcSpeciesNames));
        Raise(nameof(DcSolvent));
        Raise(nameof(DcBoxSpecies));
    }

    public void OpenDensityCalc()
    {
        if (DcSpecies.Count == 0)
        {
            // target densities are literature values at 25 °C, used only as inputs here (CRC Handbook)
            AddSpecies("Water · H₂O", "O", 0.99705m);
            AddSpecies("Toluene · C₇H₈", "Cc1ccccc1", 0.862m);
            AddSpecies("Ethanol · C₂H₆O", "CCO", 0.785m);
        }
        SetModule(56);
        Recalculate();
    }

    private void AddSpecies(string name, string smiles, decimal rho)
    {
        double mass;
        try { mass = JsonNode.Parse(CapsDocument.SmilesInfo(smiles))!["mass"]!.GetValue<double>(); }
        catch (Exception e) { DcError = $"{smiles}: {e.Message}"; return; }
        var row = new SpeciesRow { Name = name, Smiles = smiles, Mass = mass, Rho = rho };
        row.Changed = Recalculate;
        DcSpecies.Add(row);
        RaiseSpeciesNames();
    }

    public void AddSpeciesFromSmiles()
    {
        var smi = _dcNewSmiles.Trim();
        if (smi.Length == 0) return;
        DcError = "";
        AddSpecies(smi, smi, 1.0m);
        DcNewSmiles = "";
        Recalculate();
    }

    public void Recalculate()
    {
        var inv = CultureInfo.InvariantCulture;
        var v = (double)(_dcA * _dcB * _dcC);          // Å³
        var cm3 = v * 1e-24;
        DcVolume = $"{v.ToString("N0", inv)} Å³ = {Sci(cm3)} cm³";
        foreach (var s in DcSpecies)
        {
            var exact = (double)s.Rho * cm3 * Avogadro / s.Mass;
            s.Count = (int)Math.Round(exact);
            s.NExact = exact.ToString("0.0", inv);
            s.N = s.Count.ToString("N0", inv);
            s.Achieved = (s.Count * s.Mass / Avogadro / cm3).ToString("0.0000", inv);
        }
        // polymer (the Grow setup) + solvent to an overall density
        var (_, chainMass) = GrowChainSize();
        var polymerMass = _growChains * chainMass;                                      // g/mol for the whole set of chains
        DcPolymer = $"{_growChains} × {_growSpecName} DP {_growDp} · M = {chainMass.ToString("N1", inv)} g/mol per chain";
        if (_dcSolvent < DcSpecies.Count && chainMass > 0)
        {
            var sol = DcSpecies[_dcSolvent];
            var total = (double)_dcTarget * cm3 * Avogadro;                               // g/mol in the box
            var add = Math.Max(0, (total - polymerMass) / sol.Mass);
            DcSolventN = $"{Math.Round(add).ToString("N0", inv)} {sol.Name.Split(' ')[0].ToLowerInvariant()}";
            var solMass = Math.Round(add) * sol.Mass;
            DcFraction = polymerMass + solMass > 0 ? (polymerMass / (polymerMass + solMass) * 100).ToString("0.0", inv) + " %" : "—";
        }
        // a cubic box for N molecules at the species' density
        if (_dcBoxSpecies < DcSpecies.Count)
        {
            var s = DcSpecies[_dcBoxSpecies];
            var vol = (double)_dcN * s.Mass / Avogadro / (double)s.Rho * 1e24;          // Å³
            DcEdge = Math.Cbrt(vol).ToString("0.000", inv) + " Å";
        }
        RaiseSpeciesNames();
    }

    /// <summary>The box and every species (built from its SMILES) into a new Pack input.</summary>
    public void SendDensityToPack()
    {
        var dir = Path.Combine(Path.GetTempPath(), "caps-density-molecules");
        Directory.CreateDirectory(dir);
        PackXD = _dcA; PackYD = _dcB; PackZD = _dcC;
        PackPeriodic = true;
        NewPackInput();
        foreach (var s in DcSpecies.Where(s => s.Count > 0))
        {
            var path = Path.Combine(dir, string.Concat(s.Smiles.Select(c => char.IsLetterOrDigit(c) ? c : '_')) + ".pdb");
            try
            {
                var (doc, _) = CapsDocument.BuildSmiles(s.Smiles, "uff", 1, 1, s.Name);
                using (doc) doc.Save(path);
            }
            catch (Exception e) { DcError = $"{s.Name}: {e.Message}"; return; }
            PackCountD = s.Count;
            AddPackStructure(path);
        }
        SetModule(5);
        Status = $"Pack input with {DcSpecies.Count(s => s.Count > 0)} species in a {_dcA} × {_dcB} × {_dcC} Å periodic box";
    }

    // ---------------------------------------------------------------- Analyze › Surface area (design/boards/SurfaceArea)
    public bool IsSurfaceArea => _module == 57;
    private decimal _saProbe = 1.4m, _saPoints = 200;
    public decimal SaProbe { get => _saProbe; set => Set(ref _saProbe, Math.Clamp(value, 0m, 5m)); }
    public decimal SaPoints { get => _saPoints; set => Set(ref _saPoints, Math.Clamp(value, 8, 2000)); }
    public ObservableCollection<SasaGroupRow> SaGroups { get; } = new();
    public ObservableCollection<SasaConvRow> SaConvergence { get; } = new();
    public (double X, double Y)[] SaCurve { get; private set; } = [];
    private string _saTotal = "—", _saStatus = "";
    public string SaTotal { get => _saTotal; private set => Set(ref _saTotal, value); }
    public string SaStatus { get => _saStatus; private set => Set(ref _saStatus, value); }
    private double[] _saArea = [];
    public event Action? SaChanged;

    public void OpenSurfaceArea() { SetModule(57); SaChanged?.Invoke(); }

    public async Task RunSurfaceArea()
    {
        if (_doc == null) return;
        var doc = _doc;
        var inv = CultureInfo.InvariantCulture;
        SaStatus = "Shrake–Rupley…";
        var json = new JsonObject { ["probe"] = (double)_saProbe, ["points"] = (int)_saPoints, ["convergence"] = true, ["colour"] = true }.ToJsonString();
        var sw = System.Diagnostics.Stopwatch.StartNew();
        var r = JsonNode.Parse(await Task.Run(() => doc.Sasa(json)))!;
        if (r["ok"]?.GetValue<bool>() != true) { SaStatus = (string?)r["error"] ?? "failed"; return; }
        var total = r["total"]!.GetValue<double>();
        SaTotal = total.ToString("0.0", inv) + " Å²";
        _saArea = ((JsonArray)r["area"]!).Select(x => x!.GetValue<double>()).ToArray();
        SaGroups.Clear();
        SaGroups.Add(new SasaGroupRow("Whole structure", total.ToString("0.0", inv), "100"));
        foreach (var g in (JsonArray)r["groups"]!)
            SaGroups.Add(new SasaGroupRow((string)g!["name"]!, g["area"]!.GetValue<double>().ToString("0.0", inv), (g["share"]!.GetValue<double>() * 100).ToString("0.0", inv)));
        SaConvergence.Clear();
        var conv = (JsonArray?)r["convergence"] ?? new JsonArray();
        foreach (var c in conv)
            SaConvergence.Add(new SasaConvRow(((int)c!["points"]!.GetValue<double>()).ToString(inv), c["total"]!.GetValue<double>().ToString("0.0", inv),
                                              (c["delta"]!.GetValue<double>() * 100).ToString("+0.00;−0.00;0.00", inv) + " %"));
        SaCurve = conv.Select(c => (c!["points"]!.GetValue<double>(), c["total"]!.GetValue<double>())).ToArray();
        SaStatus = $"probe {_saProbe.ToString("0.00", inv)} Å · {(int)_saPoints} points per atom · {sw.Elapsed.TotalSeconds.ToString("0.0", inv)} s · view coloured by exposure";
        RenderRequested?.Invoke();
        SaChanged?.Invoke();
    }

    public void ClearSurfaceColour() { try { _doc?.SetAtomValues(null); } catch { } RenderRequested?.Invoke(); }

    public void ExportSasaCsv(string path)
    {
        if (_doc == null || _saArea.Length == 0) return;
        var inv = CultureInfo.InvariantCulture;
        using var w = new StreamWriter(path);
        w.WriteLine("atom,element,sasa_A2");
        for (var i = 0; i < _saArea.Length; ++i) w.WriteLine($"{i + 1},{_doc.Atom(i).ElementSymbol},{_saArea[i].ToString("0.###", inv)}");
        Status = $"Wrote {Path.GetFileName(path)} · {_saArea.Length} atoms";
    }

    // ---------------------------------------------------------------- Studio › Unit cell (design/boards/CellEditor)
    public bool IsCellEditor => _module == 58;
    private decimal _ceA = 10, _ceB = 10, _ceC = 10, _ceAl = 90, _ceBe = 90, _ceGa = 90, _ceZ = 1, _ceSupA = 2, _ceSupB = 2, _ceSupC = 2;
    private double _ceFormula;
    private bool _ceScale = true;
    public decimal CeA { get => _ceA; set { if (Set(ref _ceA, Math.Clamp(value, 0.5m, 5000))) RaiseCell(); } }
    public decimal CeB { get => _ceB; set { if (Set(ref _ceB, Math.Clamp(value, 0.5m, 5000))) RaiseCell(); } }
    public decimal CeC { get => _ceC; set { if (Set(ref _ceC, Math.Clamp(value, 0.5m, 5000))) RaiseCell(); } }
    public decimal CeAlpha { get => _ceAl; set { if (Set(ref _ceAl, Math.Clamp(value, 20, 160))) RaiseCell(); } }
    public decimal CeBeta { get => _ceBe; set { if (Set(ref _ceBe, Math.Clamp(value, 20, 160))) RaiseCell(); } }
    public decimal CeGamma { get => _ceGa; set { if (Set(ref _ceGa, Math.Clamp(value, 20, 160))) RaiseCell(); } }
    public decimal CeZ { get => _ceZ; set { if (Set(ref _ceZ, Math.Clamp(value, 1, 1000))) RaiseCell(); } }
    public bool CeScale { get => _ceScale; set => Set(ref _ceScale, value); }
    public decimal CeSupA { get => _ceSupA; set => Set(ref _ceSupA, Math.Clamp(value, 1, 10)); }
    public decimal CeSupB { get => _ceSupB; set => Set(ref _ceSupB, Math.Clamp(value, 1, 10)); }
    public decimal CeSupC { get => _ceSupC; set => Set(ref _ceSupC, Math.Clamp(value, 1, 10)); }
    private string _ceFormulaText = "", _ceError = "";
    public string CeFormulaText { get => _ceFormulaText; private set => Set(ref _ceFormulaText, value); }
    public string CeError { get => _ceError; private set => Set(ref _ceError, value); }
    public ObservableCollection<CellPresetRow> CellPresets { get; } = new();

    /// <summary>Volume from the six parameters (triclinic).</summary>
    private double CellVolume()
    {
        double r(decimal d) => (double)d * Math.PI / 180;
        var (ca, cb, cg) = (Math.Cos(r(_ceAl)), Math.Cos(r(_ceBe)), Math.Cos(r(_ceGa)));
        var f = 1 - ca * ca - cb * cb - cg * cg + 2 * ca * cb * cg;
        return (double)(_ceA * _ceB * _ceC) * Math.Sqrt(Math.Max(0, f));
    }
    public string CeVolume => CellVolume().ToString("0.00", CultureInfo.InvariantCulture) + " Å³";
    public string CeDensity => _ceFormula > 0 && CellVolume() > 0 ? ((double)_ceZ * _ceFormula / Avogadro / (CellVolume() * 1e-24)).ToString("0.0000", CultureInfo.InvariantCulture) + " g/cm³" : "—";
    public string CeVectors
    {
        get
        {
            double r(decimal d) => (double)d * Math.PI / 180;
            double a = (double)_ceA, b = (double)_ceB, c = (double)_ceC;
            var (al, be, ga) = (r(_ceAl), r(_ceBe), r(_ceGa));
            var bx = b * Math.Cos(ga); var by = b * Math.Sin(ga);
            var cx = c * Math.Cos(be); var cy = c * (Math.Cos(al) - Math.Cos(be) * Math.Cos(ga)) / Math.Sin(ga);
            var cz = Math.Sqrt(Math.Max(0, c * c - cx * cx - cy * cy));
            string F(double v) => v.ToString("0.000", CultureInfo.InvariantCulture).PadLeft(8);
            return $"a {F(a)} {F(0)} {F(0)}\nb {F(bx)} {F(by)} {F(0)}\nc {F(cx)} {F(cy)} {F(cz)}";
        }
    }
    private void RaiseCell() { foreach (var n in new[] { nameof(CeVolume), nameof(CeDensity), nameof(CeVectors) }) Raise(n); }

    public void OpenCellEditor()
    {
        if (CellPresets.Count == 0)
        {
            // crystallographic cells from the literature (lattice only; densities computed from them)
            void Preset(string name, double a, double b, double c, double al, double be, double ga, int z, double fm, string lit)
            {
                double rr(double d) => d * Math.PI / 180;
                var f = 1 - Math.Pow(Math.Cos(rr(al)), 2) - Math.Pow(Math.Cos(rr(be)), 2) - Math.Pow(Math.Cos(rr(ga)), 2) + 2 * Math.Cos(rr(al)) * Math.Cos(rr(be)) * Math.Cos(rr(ga));
                var v = a * b * c * Math.Sqrt(f);
                CellPresets.Add(new CellPresetRow(name, v.ToString("0.00", CultureInfo.InvariantCulture), z.ToString(CultureInfo.InvariantCulture),
                    (z * fm / Avogadro / (v * 1e-24)).ToString("0.000", CultureInfo.InvariantCulture), lit, a, b, c, al, be, ga, z, fm));
            }
            Preset("PE (Bunn 1939)", 7.40, 4.93, 2.534, 90, 90, 90, 2, 28.054, "≈ 1.00");
            Preset("α-quartz", 4.9134, 4.9134, 5.4052, 90, 90, 120, 3, 60.084, "2.65");
            Preset("PET (Daubeny 1954)", 4.56, 5.94, 10.75, 98.5, 118, 112, 1, 192.17, "1.455");
        }
        if (_doc?.Summary() is { CellValid: not 0 } s)
        {
            // the open cell: its parameters from the summary lengths; the angles from the document's cell vectors when known
            _ceA = (decimal)Math.Round(s.CellA, 4); _ceB = (decimal)Math.Round(s.CellB, 4); _ceC = (decimal)Math.Round(s.CellC, 4);
            try
            {
                if (JsonNode.Parse(_doc.Periodic("{}"))?["angles"] is JsonArray an && an.Count == 3)
                    (_ceAl, _ceBe, _ceGa) = ((decimal)Math.Round(an[0]!.GetValue<double>(), 3), (decimal)Math.Round(an[1]!.GetValue<double>(), 3), (decimal)Math.Round(an[2]!.GetValue<double>(), 3));
            }
            catch { /* orthogonal */ }
            foreach (var n in new[] { nameof(CeAlpha), nameof(CeBeta), nameof(CeGamma) }) Raise(n);
            _ceFormula = s.TotalMass; _ceZ = 1;
            CeFormulaText = $"the open structure · {s.TotalMass.ToString("N1", CultureInfo.InvariantCulture)} g/mol per cell";
            foreach (var n in new[] { nameof(CeA), nameof(CeB), nameof(CeC), nameof(CeZ) }) Raise(n);
        }
        SetModule(58);
        RaiseCell();
    }

    public void UseCellPreset(CellPresetRow p)
    {
        (_ceA, _ceB, _ceC, _ceAl, _ceBe, _ceGa, _ceZ) = ((decimal)p.A, (decimal)p.B, (decimal)p.C, (decimal)p.Al, (decimal)p.Be, (decimal)p.Ga, p.ZCount);
        _ceFormula = p.FormulaMass;
        CeFormulaText = $"{p.Name} · Z = {p.ZCount} × {p.FormulaMass.ToString("0.###", CultureInfo.InvariantCulture)} g/mol";
        foreach (var n in new[] { nameof(CeA), nameof(CeB), nameof(CeC), nameof(CeAlpha), nameof(CeBeta), nameof(CeGamma), nameof(CeZ) }) Raise(n);
        RaiseCell();
    }

    public void ApplyCell()
    {
        if (_doc == null) return;
        var j = new JsonObject { ["a"] = (double)_ceA, ["b"] = (double)_ceB, ["c"] = (double)_ceC, ["alpha"] = (double)_ceAl, ["beta"] = (double)_ceBe, ["gamma"] = (double)_ceGa, ["scale"] = _ceScale };
        try { _doc.SetCell(j.ToJsonString()); CeError = ""; }
        catch (Exception e) { CeError = e.Message; return; }
        AfterEdit($"Cell applied · V = {CeVolume}" + (_ceScale ? " · atoms scaled" : "") + " (undo with ⌘Z)");
    }

    public void MakeSupercell()
    {
        if (_doc == null) return;
        try { _doc.Supercell((int)_ceSupA, (int)_ceSupB, (int)_ceSupC); CeError = ""; }
        catch (Exception e) { CeError = e.Message; return; }
        Frames = 1; Raise(nameof(FrameMax)); _frame = 0; Raise(nameof(Frame));
        AfterEdit($"Supercell {_ceSupA} × {_ceSupB} × {_ceSupC} · {_doc.Summary().Atoms:N0} atoms (undo with ⌘Z)");
        if (_doc.Summary() is { CellValid: not 0 } s) { _ceA = (decimal)Math.Round(s.CellA, 4); _ceB = (decimal)Math.Round(s.CellB, 4); _ceC = (decimal)Math.Round(s.CellC, 4); foreach (var n in new[] { nameof(CeA), nameof(CeB), nameof(CeC) }) Raise(n); RaiseCell(); }
    }
}
