using System.Collections.ObjectModel;
using System.Globalization;
using System.Text;
using System.Text.Json;
using Avalonia.Media;
using CapsStudio.Interop;

namespace CapsStudio.ViewModels;

/// <summary>One calculation that can be switched on (a property id of caps analyze), or one that is not built yet.</summary>
public sealed class CalcChip : ObservableObject
{
    public string Id { get; init; } = "";
    public string Label { get; init; } = "";
    public string Tip { get; init; } = "";
    public bool Available { get; init; } = true;
    private bool _on;
    public bool IsOn { get => _on; set => Set(ref _on, value && Available); }
    private bool _active;
    /// <summary>Shown on the focused page it belongs to (Mechanics, Scattering, Free volume).</summary>
    public bool Active { get => _active; set => Set(ref _active, value); }
}

public sealed record CalcGroup(string Name, CalcChip[] Chips);

/// <summary>An experimental or literature range for one property of one material (data/reference/polymers.json).</summary>
public sealed record RefValue(double Lo, double Hi, string Unit, string Condition, string Source);

public sealed record RefMaterial(string Id, string Name, Dictionary<string, RefValue> Values)
{
    public override string ToString() => Name;
}

/// <summary>A property result as a card: value ± error, method, notes, further numbers and the comparison with experiment.</summary>
public sealed class ResultCard
{
    public string Id { get; init; } = "";
    public string Name { get; init; } = "";
    public double Value { get; init; } = double.NaN;
    public double Error { get; init; } = double.NaN;
    public string Unit { get; init; } = "";
    public string Method { get; init; } = "";
    public string[] Notes { get; init; } = [];
    public (string Key, double Value)[] Extra { get; init; } = [];
    public RefValue? Ref { get; set; }
    // normal modes: the card plays one (1 = the lowest) in a copy of the structure
    public bool IsModes => Id == "modes" && HasValue;
    public bool IsConformers => Id == "conformers" && HasValue;
    public decimal PickMode { get; set; } = 1;
    public decimal ModeCount => (decimal)Math.Max(1, Extra.FirstOrDefault(e => e.Key == "modes").Value);

    private static readonly CultureInfo Inv = CultureInfo.InvariantCulture;
    public bool HasValue => !double.IsNaN(Value);
    private bool HasError => HasValue && !double.IsNaN(Error) && Error > 0;
    /// <summary>Decimals that show the error to two significant figures (the value is rounded to match).</summary>
    private int ErrorDecimals => Math.Clamp(1 - (int)Math.Floor(Math.Log10(Error)), 0, 8);
    // shown in the display units chosen in Settings › Units (the stored value keeps CAPS's units)
    private double F => DisplayUnits.For(Unit).Factor;
    private int ShownDecimals => Math.Clamp(1 - (int)Math.Floor(Math.Log10(Error * F)), 0, 8);
    public string ValueText => !HasValue ? "—" : HasError ? (Value * F).ToString("F" + ShownDecimals, Inv) : Num(Value * F);
    public string ErrorText => HasError ? "± " + (Error * F).ToString("F" + ShownDecimals, Inv) : "";
    public string UnitText => DisplayUnits.For(Unit).Unit.Replace("^½", "½");
    public string NotesText => string.Join("\n", Notes);
    public bool HasNotes => Notes.Length > 0;
    public string ExtraText => string.Join("\n", Extra.Select(e => $"{e.Key}: {Num(e.Value)}"));
    public bool HasExtra => Extra.Length > 0;
    public bool HasRef => Ref != null;
    public string RefText => Ref == null ? "" :
        $"{(Ref.Source.StartsWith("Flory", StringComparison.Ordinal) ? "lit." : "exp.")} {(Ref.Lo == Ref.Hi ? Num(Ref.Lo) : Num(Ref.Lo) + "–" + Num(Ref.Hi))} {Ref.Unit.Replace("^½", "½")}".TrimEnd() +
        (Ref.Condition.Length > 0 ? " · " + Ref.Condition : "");
    public string RefSource => Ref?.Source ?? "";
    /// <summary>Within the range (allowing the error), outside it, or not compared.</summary>
    public string Verdict
    {
        get
        {
            if (Ref == null || !HasValue) return "";
            var e = double.IsNaN(Error) ? 0 : Error;
            var uncertain = e > 0.3 * Math.Abs(Value);
            if (Value + e >= Ref.Lo && Value - e <= Ref.Hi)
                return uncertain ? string.Format(Inv, "overlaps the range, but the error is {0:F0} % of the value", 100 * e / Math.Abs(Value)) : "within the range";
            var d = Value < Ref.Lo ? (Ref.Lo - Value) / Math.Max(Math.Abs(Ref.Lo), 1e-12) : (Value - Ref.Hi) / Math.Max(Math.Abs(Ref.Hi), 1e-12);
            return string.Format(Inv, "{0:F1} % {1} the range", 100 * d, Value < Ref.Lo ? "below" : "above");
        }
    }
    public IBrush VerdictBrush => Verdict.StartsWith("within", StringComparison.Ordinal) ? Good : Warn;   // "overlaps … error" stays amber
    internal static IBrush Good => Tokens.Brush("OkB");
    internal static IBrush Warn => Tokens.Brush("WarnB");

    public static string Num(double v)
    {
        if (double.IsNaN(v)) return "—";
        var a = Math.Abs(v);
        if (a == 0) return "0";
        if (a >= 1e5 || a < 1e-3) return v.ToString("0.###e0", Inv);
        var digits = Math.Max(0, 3 - (int)Math.Floor(Math.Log10(a)));   // four significant figures
        return Math.Round(v, Math.Min(digits, 6)).ToString("0.######", Inv);
    }
}

/// <summary>A curve behind a result.</summary>
public sealed record SeriesItem(string Property, string Label, string XLabel, string YLabel, double[] X, double[] Y, bool LogLog, double? RefY)
{
    /// <summary>Draw the data as points (with Overlay as the line through them: a fit or a smoothed curve).</summary>
    public bool Markers { get; init; }
    public double[]? OverlayX { get; init; }
    public double[]? OverlayY { get; init; }
    public override string ToString() => $"{Property} · {Label}";
}

/// <summary>Analyze › Properties: choose calculations, frames and groups, run them on the open trajectory, and read the
/// results as cards compared with experiment, with the curves behind them. Export CSV and a LaTeX table.</summary>
public sealed class AnalyzeViewModel : ObservableObject
{
    private readonly Func<CapsDocument?> _doc;
    private readonly Action<string> _status;
    private readonly Action<bool> _running;
    private static readonly CultureInfo Inv = CultureInfo.InvariantCulture;

    public AnalyzeViewModel(Func<CapsDocument?> doc, Action<string> status, Action<bool> running)
    {
        _doc = doc;
        _status = status;
        _running = running;
        Groups =
        [
            new("Structure", [Chip("density", "Density", on: true), Chip("rdf", "RDF", on: true), Chip("sq", "S(q)"), Chip("xray", "X-ray"), Chip("electron", "Electron"), Chip("neutron", "Neutron")]),
            new("Chains", [Chip("rg", "Rg", on: true), Chip("ree", "Ree"), Chip("cn", "Cn, C∞"), Chip("persistence", "Persistence"), Chip("orientation", "Orientation"),
                Chip("entanglements", "Entanglements"), Chip("conformation", "Torsions"), Chip("p2r", "P₂(r)")]),
            new("Thermo", [Chip("ced", "CED"), Chip("delta", "δ"), Chip("dielectric", "ε dielectric"), Chip("fluct", "Cp·κT·α"), Chip("modes", "Normal modes"), ConfChip, TgChip]),
            new("Mechanics", [StrainChip, FluctChip, TensileChip, CreepChip]),
            new("Dynamics", [Chip("msd", "MSD"), Chip("diffusion", "D"), Chip("relaxation", "Relaxation"), Chip("vacf", "VACF · VDOS"), Chip("vanhove", "van Hove · α₂"), ViscChip, NemdChip]),
            new("Free volume", [Chip("ffv", "Probe insertion"), Chip("psd", "Pore size")]),
            new("Interface", [Chip("zprofile", "z profile"), Chip("adhesion", "Adhesion"), Chip("interaction", "Filler–matrix"), PullShearChip, PullNormalChip, FrictionChip]),
            new("Rubber network", [Chip("crosslinks", "Crosslink density"), Chip("hbonds", "H-bonds")]),
        ];
        LoadReferences();
        PullShearChip.PropertyChanged += (_, _) => Raise(nameof(PullOn));
        FluctChip.PropertyChanged += (_, _) => Raise(nameof(NvtRunOn));
        ViscChip.PropertyChanged += (_, _) => Raise(nameof(NvtRunOn));
        PullNormalChip.PropertyChanged += (_, _) => Raise(nameof(PullOn));
        foreach (var c in new[] { FluctChip, ViscChip, TensileChip, PullShearChip, PullNormalChip }) c.PropertyChanged += (_, _) => Raise(nameof(ShearOwnT));
    }

    // protocols (their settings show when switched on)
    public CalcChip TgChip { get; } = Chip("tg", "Tg", tip: "Glass transition from a stepwise NPT cooling run of the current frame (a copy: the document is not changed)");
    public CalcChip StrainChip { get; } = Chip("cij_strain", "Cij strain", tip: "Static elastic constants: minimise, strain ±ε in each direction, re-minimise (Theodorou & Suter)");
    /// <summary>A Green–Kubo viscosity came out (mPa·s): the Diffusion page's finite-size correction takes it.</summary>
    public event Action<double>? ViscosityComputed;
    public CalcChip ViscChip { get; } = Chip("viscosity", "Viscosity", tip: "Shear viscosity by Green–Kubo: an NVT run of the current frame (Nosé–Hoover), the pressure tensor's autocorrelation integrated (Daivis & Evans); melts need long runs");
    public bool NvtRunOn => FluctChip.IsOn || ViscChip.IsOn;
    public CalcChip ConfChip { get; } = Chip("conformers", "Conformers", tip: "Conformer search of the frame shown or the selected molecule, in vacuum: minimised starts clustered by heavy-atom RMSD, energies and Boltzmann populations; open them as frames from the card");
    private int _confTrials = 50, _confMethod;
    private double _confWindow = 10, _confRmsd = 0.5;
    private bool _confSel = true;
    public decimal ConfTrialsD { get => _confTrials; set => Set(ref _confTrials, (int)Math.Clamp(value, 1, 5000), nameof(ConfTrialsD)); }
    public int ConfMethod { get => _confMethod; set => Set(ref _confMethod, Math.Clamp(value, 0, 1)); }
    public static readonly string[] ConfMethods = ["Random staggered torsions", "Anneal: 1000 K snapshots quenched"];
    public decimal ConfWindowD { get => (decimal)_confWindow; set => Set(ref _confWindow, Math.Max(0.1, (double)value), nameof(ConfWindowD)); }
    public decimal ConfRmsdD { get => (decimal)_confRmsd; set => Set(ref _confRmsd, Math.Max(0.05, (double)value), nameof(ConfRmsdD)); }
    public bool ConfSelection { get => _confSel; set => Set(ref _confSel, value); }
    /// <summary>The search as the frames call takes it (the same settings and seed as the protocol: the same conformers).</summary>
    public string ConfJson(double temperature) => new System.Text.Json.Nodes.JsonObject
    {
        ["trials"] = _confTrials, ["method"] = _confMethod == 1 ? "anneal" : "torsions", ["selection"] = _confSel, ["window"] = _confWindow,
        ["rmsd"] = _confRmsd, ["temperature"] = temperature, ["seed"] = (long)MechSeed,
    }.ToJsonString();
    public CalcChip FrictionChip { get; } = Chip("friction", "Sliding friction", tip: "One wall molecule slid over the film at a fixed gap on a copy of the current frame: the friction force, the shear stress at the wall and μ = −F_x/F_z");
    private int _frMoving = 2, _frFixed = 1;
    private double _frSpeed = 10, _frPs = 50, _frEq = 10, _frT = 300;   // m/s, ps, ps, K
    public decimal FrMovingD { get => _frMoving; set => Set(ref _frMoving, (int)Math.Max(1, value), nameof(FrMovingD)); }
    public decimal FrFixedD { get => _frFixed; set => Set(ref _frFixed, (int)Math.Max(0, value), nameof(FrFixedD)); }
    public decimal FrSpeedD { get => (decimal)_frSpeed; set { Set(ref _frSpeed, (double)value, nameof(FrSpeedD)); Raise(nameof(FrText)); } }
    public decimal FrPsD { get => (decimal)_frPs; set { Set(ref _frPs, Math.Max(0.1, (double)value), nameof(FrPsD)); Raise(nameof(FrText)); } }
    public decimal FrEqD { get => (decimal)_frEq; set => Set(ref _frEq, Math.Max(0, (double)value), nameof(FrEqD)); }
    public decimal FrTD { get => (decimal)_frT; set => Set(ref _frT, Math.Max(1, (double)value), nameof(FrTD)); }
    public string FrText => string.Format(CultureInfo.InvariantCulture, "{0:0.###} Å/ps · the wall moves {1:0.#} Å in {2:0.#} ps", _frSpeed / 100, _frSpeed / 100 * _frPs, _frPs);
    public CalcChip CreepChip { get; } = Chip("creep", "Creep", tip: "Constant true stress along one axis on a copy of the current frame: strain and creep compliance against time, the creep rate and the lateral contraction");
    private double _crStress = 50, _crT = 300, _crPs = 200, _crEq = 20;
    private int _crAxis = 2;
    public decimal CreepStressD { get => (decimal)_crStress; set { Set(ref _crStress, (double)value, nameof(CreepStressD)); Raise(nameof(CreepText)); } }
    public decimal CreepTD { get => (decimal)_crT; set => Set(ref _crT, Math.Max(1, (double)value), nameof(CreepTD)); }
    public decimal CreepPsD { get => (decimal)_crPs; set { Set(ref _crPs, Math.Max(1, (double)value), nameof(CreepPsD)); Raise(nameof(CreepText)); } }
    public decimal CreepEqD { get => (decimal)_crEq; set => Set(ref _crEq, Math.Max(0, (double)value), nameof(CreepEqD)); }
    public int CreepAxis { get => _crAxis; set => Set(ref _crAxis, Math.Clamp(value, 0, 2)); }
    public static readonly string[] CreepAxes = ["x", "y", "z"];
    public string CreepText => string.Format(CultureInfo.InvariantCulture, "{0} {1:0.##} MPa for {2:0.#} ps ({3:N0} steps at 1 fs)",
        _crStress >= 0 ? "tension" : "compression", Math.Abs(_crStress), _crPs, _crPs * 1000);
    public CalcChip NemdChip { get; } = Chip("nemd", "η(γ̇) shear", tip: "Shear viscosity at imposed shear rates (SLLOD, non-equilibrium MD): η(γ̇), shear thinning and the first normal-stress difference, each rate on a copy of the current frame");
    // NEMD: rates spaced evenly in log γ̇, each its own run from one NVT-equilibrated start
    private double _shLo = 0.01, _shHi = 0.1, _shPs = 50, _shEq = 10, _shT = 300;
    private int _shPoints = 3;
    public decimal ShearLoD { get => (decimal)_shLo; set { Set(ref _shLo, Math.Max(1e-6, (double)value), nameof(ShearLoD)); Raise(nameof(ShearRatesText)); } }
    public decimal ShearHiD { get => (decimal)_shHi; set { Set(ref _shHi, Math.Max(1e-6, (double)value), nameof(ShearHiD)); Raise(nameof(ShearRatesText)); } }
    public decimal ShearPointsD { get => _shPoints; set { Set(ref _shPoints, (int)Math.Clamp(value, 1, 12), nameof(ShearPointsD)); Raise(nameof(ShearRatesText)); } }
    public decimal ShearPsD { get => (decimal)_shPs; set { Set(ref _shPs, Math.Max(0.1, (double)value), nameof(ShearPsD)); Raise(nameof(ShearRatesText)); } }
    public decimal ShearEqD { get => (decimal)_shEq; set => Set(ref _shEq, Math.Max(0, (double)value), nameof(ShearEqD)); }
    public decimal ShearTD { get => (decimal)_shT; set => Set(ref _shT, Math.Max(1, (double)value), nameof(ShearTD)); }
    /// <summary>One temperature goes to every protocol of a run: another protocol's, when one is on with the shear.</summary>
    public bool ShearOwnT => !(PullOn || NvtRunOn || TensileChip.IsOn);
    /// <summary>The rates the sweep runs, in 1/ps and 1/s, and the strain each reaches.</summary>
    public string ShearRatesText
    {
        get
        {
            var inv = CultureInfo.InvariantCulture;
            var n = _shHi > _shLo ? Math.Max(2, _shPoints) : 1;
            var rates = Enumerable.Range(0, n).Select(k => n == 1 ? _shLo : _shLo * Math.Pow(_shHi / _shLo, k / (double)(n - 1))).ToArray();
            return string.Join(" · ", rates.Select(r => r.ToString("0.####", inv))) + " /ps  (" + (rates[0] * 1e12).ToString("0.#E+0", inv) + "–" +
                   (rates[^1] * 1e12).ToString("0.#E+0", inv) + " s⁻¹) · strain " + (rates[0] * _shPs).ToString("0.##", inv) + "–" + (rates[^1] * _shPs).ToString("0.##", inv);
        }
    }
    public CalcChip FluctChip { get; } = Chip("cij_run", "Cij fluct.", tip: "Elastic constants from stress fluctuations: an NVT run of the current frame, the stress sampled at every step (Lutsko; Clavier et al.)");
    public CalcChip PullShearChip { get; } = Chip("pull_shear", "Pull · shear", tip: "Steered MD: the film dragged along x over the held surface (molecule 1); interfacial shear strength and work");
    public CalcChip PullNormalChip { get; } = Chip("pull_normal", "Pull · normal", tip: "Steered MD: the film pulled off the held surface along +z (needs vacuum above the film); peak normal stress and work of separation");
    public CalcChip TensileChip { get; } = Chip("tensile", "Stress–strain", tip: "Uniaxial deformation MD of the current frame: modulus, Poisson ratio, yield");

    private static CalcChip Chip(string id, string label, bool on = false, string tip = "") => new() { Id = id, Label = label, IsOn = on, Tip = tip.Length > 0 ? tip : Tips.GetValueOrDefault(id, "") };
    private static CalcChip Soon(string label, string tip) => new() { Id = "", Label = label, Tip = tip, Available = false };

    private static readonly Dictionary<string, string> Tips = new()
    {
        ["density"] = "Mass over volume, averaged over the frames",
        ["rdf"] = "g(r) of the chosen pair, averaged over the frames; first peak and coordination number",
        ["sq"] = "Total structure factor: direct reciprocal-lattice sum at low q, g(r) transform above",
        ["xray"] = "X-ray I(q) with Cromer–Mann form factors (Faber–Ziman)",
        ["electron"] = "Electron-diffraction I(q) with Peng et al. 1996 elastic scattering factors (Faber–Ziman)",
        ["neutron"] = "Neutron S(q) with coherent scattering lengths (Faber–Ziman)",
        ["rg"] = "Radius of gyration of the chains, √⟨Rg²⟩",
        ["ree"] = "Backbone end-to-end distance, √⟨R²⟩",
        ["cn"] = "Characteristic ratio C_n and C∞ extrapolated in 1/n",
        ["persistence"] = "Persistence length: Flory projection and bond-correlation decay",
        ["ced"] = "Cohesive energy density (E isolated − E bulk)/V with the force field (Field, else GAFF for C and H, UFF otherwise), split into its van der Waals and electrostatic parts",
        ["delta"] = "Hildebrand solubility parameter δ = √CED",
        ["cij_fluct"] = "Elastic constants from stress fluctuations of the saved frames of an NVT run",
        ["msd"] = "Mean-square displacement of atoms and molecule centres, all time origins, drift removed",
        ["diffusion"] = "Diffusion coefficient from the linear part of the molecule-centre MSD (Einstein)",
        ["relaxation"] = "End-to-end and segmental (P2) autocorrelations with KWW fits",
        ["ffv"] = "Free volume by probe insertion on a grid: accessible fraction and Bondi FFV",
        ["psd"] = "Pore size distribution: largest atom-free sphere containing each free point",
        ["orientation"] = "Nematic order S of backbone chords, director, Herman's f along z, local crystallinity, and P₂ against height (orientation near a surface)",
        ["entanglements"] = "Primitive-path analysis (Everaers 2004): chain ends fixed, chains pulled tight without crossing; entanglement length N_e (modified S-coil), M_e, tube step and plateau modulus G_N⁰ at the temperature set for Tg / Cij",
        ["conformation"] = "Backbone torsion distribution: trans (|φ| > 120°) and gauche± fractions, t/g ratio, with the backbone angle and bond distributions",
        ["vacf"] = "Velocity autocorrelation (all time origins), Green–Kubo D, and the vibrational density of states (needs velocities in the frames)",
        ["modes"] = "Harmonic vibrations of the last chosen frame (minimise it tightly first; at most 600 atoms, or a group): wavenumbers, the IR spectrum of the fixed charges, the vibrational density of states, ZPE, S and Cv; play a mode from its card",
        ["fluct"] = "Heat capacity, isothermal compressibility, bulk modulus and thermal expansion from the fluctuations of an equilibrium run (NPT: Cp, κT, B, α; NVT: Cv). Needs the force field and the run's temperature; hundreds of well-spaced frames",
        ["dielectric"] = "Static dielectric constant from the fluctuations of the cell's total dipole (needs charges, the run's temperature and a long equilibrium trajectory), with the dipole autocorrelation",
        ["vanhove"] = "Self part of the van Hove function G_s(r, t) at log-spaced times and the non-Gaussian parameter α₂(t): heterogeneous, glassy dynamics",
        ["p2r"] = "Orientational correlation of backbone chords against their distance: local chain alignment, order near a surface or in a stretched sample",
        ["hbonds"] = "Hydrogen bonds by the Luzar–Chandler criterion (D···A ≤ 3.5 Å, ∠H–D···A ≤ 30°): count per frame and lifetime from the intermittent correlation",
        ["crosslinks"] = "Sulfur bridges (mono-, di-, polysulfidic), pendant groups, crosslink density ν and strand mass Mc = ρ/2ν",
        ["zprofile"] = "Mass density along z for the surface (molecule 1) and the film: first-layer peak and the film's own density",
        ["adhesion"] = "Work of adhesion −(E all − E surface − E film)/area between the surface (molecule 1) and the film, with the force field",
        ["interaction"] = "Filler–matrix interaction energy E all − E filler − E matrix (kcal/mol; the filler is molecule 1: a functionalised tube, sheet or particle in its matrix), with the van der Waals and Coulomb parts",
    };

    public CalcGroup[] Groups { get; }

    // ---------------------------------------------------------------- source
    private int _first, _last = -1, _stride = 1, _pair;
    private double _framePs, _timestepFs = 1, _probe, _grid = 0.4, _fitFrom = 0.2, _fitTo = 0.5;
    private bool _inter;
    public decimal FirstD { get => _first; set => Set(ref _first, (int)Math.Max(0, value), nameof(FirstD)); }
    /// <summary>Last frame, −1 for the last one.</summary>
    public decimal LastD { get => _last; set => Set(ref _last, (int)Math.Max(-1, value), nameof(LastD)); }
    public decimal StrideD { get => _stride; set => Set(ref _stride, (int)Math.Max(1, value), nameof(StrideD)); }
    /// <summary>Time between frames in ps; 0 takes it from the timesteps in the file × the MD time step.</summary>
    public decimal FramePsD { get => (decimal)_framePs; set => Set(ref _framePs, (double)Math.Max(0, value), nameof(FramePsD)); }
    public decimal TimestepFsD { get => (decimal)_timestepFs; set => Set(ref _timestepFs, (double)Math.Max(0.01m, value), nameof(TimestepFsD)); }
    public int PairIndex { get => _pair; set => Set(ref _pair, value); }
    public bool InterOnly { get => _inter; set => Set(ref _inter, value); }
    private double _qmax = 25, _dq = 0.02, _qDirect = 4;
    private int _deuterate;
    public decimal QmaxD { get => (decimal)_qmax; set => Set(ref _qmax, (double)Math.Clamp(value, 1m, 40m), nameof(QmaxD)); }
    public decimal DqD { get => (decimal)_dq; set => Set(ref _dq, (double)Math.Clamp(value, 0.002m, 0.2m), nameof(DqD)); }
    public decimal QDirectD { get => (decimal)_qDirect; set => Set(ref _qDirect, (double)Math.Clamp(value, 0m, 10m), nameof(QDirectD)); }
    /// <summary>Neutron contrast: 0 none, 1 every H → D, 2 aliphatic H (d-backbone), 3 aromatic H (d-ring), 4 H on O/N.</summary>
    public int Deuterate { get => _deuterate; set => Set(ref _deuterate, Math.Clamp(value, 0, 4)); }
    private string _radii = "bondi";
    /// <summary>Free volume and pore radii: bondi, uff or forcefield (the Field assignment).</summary>
    public string Radii { get => _radii; set => Set(ref _radii, value); }
    public decimal ProbeD { get => (decimal)_probe; set => Set(ref _probe, (double)Math.Max(0, value), nameof(ProbeD)); }
    public decimal GridD { get => (decimal)_grid; set => Set(ref _grid, (double)Math.Clamp(value, 0.1m, 2m), nameof(GridD)); }
    public decimal FitFromD { get => (decimal)_fitFrom; set => Set(ref _fitFrom, (double)Math.Clamp(value, 0m, 0.95m), nameof(FitFromD)); }
    public decimal FitToD { get => (decimal)_fitTo; set => Set(ref _fitTo, (double)Math.Clamp(value, 0.05m, 1m), nameof(FitToD)); }
    // protocol settings
    private double _fluctPs = 100, _eqPs = 20;
    /// <summary>Unsampled NPT run before the tensile pull and the first cooling hold (ps).</summary>
    public decimal EqPsD { get => (decimal)_eqPs; set => Set(ref _eqPs, (double)Math.Max(0, value), nameof(EqPsD)); }
    public decimal FluctPsD { get => (decimal)_fluctPs; set => Set(ref _fluctPs, (double)Math.Max(1, value), nameof(FluctPsD)); }
    private double _tgFrom = 500, _tgTo = 200, _tgStep = 20, _tgPs = 100, _cijStrain = 1e-4, _fluctT = 300, _tensRate = 1e-3, _tensMax = 0.1, _tensT = 300;
    private int _cijConfigs = 1, _tensAxis;
    private bool _tensFixed;
    public decimal TgFromD { get => (decimal)_tgFrom; set => Set(ref _tgFrom, (double)value, nameof(TgFromD)); }
    public decimal TgToD { get => (decimal)_tgTo; set => Set(ref _tgTo, (double)value, nameof(TgToD)); }
    public decimal TgStepD { get => (decimal)_tgStep; set => Set(ref _tgStep, (double)Math.Max(1, value), nameof(TgStepD)); }
    public decimal TgPsD { get => (decimal)_tgPs; set => Set(ref _tgPs, (double)Math.Max(1, value), nameof(TgPsD)); }
    // the cooling scan's ensemble and how Tg is read
    private double _tgPressure = 1, _tgTauT = 100, _tgTauP = 1000, _tgAverage = 50, _tgGlassy, _tgRubbery;
    private int _tgBarostat, _tgProperty, _tgFit;
    public static readonly string[] TgBarostats = ["Stochastic cell rescaling (Bernetti & Bussi 2020)", "Berendsen", "MTK (Nosé–Hoover chains)"];
    public static readonly string[] TgProperties = ["Specific volume", "Potential energy per atom"];
    public static readonly string[] TgFits = ["Two lines, break found by the fit", "Two lines through a glassy and a rubbery range"];
    public decimal TgPressureD { get => (decimal)_tgPressure; set => Set(ref _tgPressure, (double)Math.Clamp(value, 0.001m, 1e6m), nameof(TgPressureD)); }
    public int TgBarostat { get => _tgBarostat; set => Set(ref _tgBarostat, Math.Clamp(value, 0, 2)); }
    public decimal TgTauTD { get => (decimal)_tgTauT; set => Set(ref _tgTauT, (double)Math.Clamp(value, 1m, 1e6m), nameof(TgTauTD)); }
    public decimal TgTauPD { get => (decimal)_tgTauP; set => Set(ref _tgTauP, (double)Math.Clamp(value, 10m, 1e7m), nameof(TgTauPD)); }
    /// <summary>Per cent of each hold averaged (the rest, at its start, discarded).</summary>
    public decimal TgAverageD { get => (decimal)_tgAverage; set => Set(ref _tgAverage, (double)Math.Clamp(value, 5m, 100m), nameof(TgAverageD)); }
    public int TgProperty { get => _tgProperty; set => Set(ref _tgProperty, Math.Clamp(value, 0, 1)); }
    public int TgFit { get => _tgFit; set { if (Set(ref _tgFit, Math.Clamp(value, 0, 1))) Raise(nameof(TgFitRanges)); } }
    public bool TgFitRanges => _tgFit == 1;
    /// <summary>The glassy range ends / the rubbery range starts (K; 0: the lowest / highest third of the scan).</summary>
    public decimal TgGlassyMaxD { get => (decimal)_tgGlassy; set => Set(ref _tgGlassy, (double)Math.Max(0, value), nameof(TgGlassyMaxD)); }
    public decimal TgRubberyMinD { get => (decimal)_tgRubbery; set => Set(ref _tgRubbery, (double)Math.Max(0, value), nameof(TgRubberyMinD)); }
    public decimal CijStrainD { get => (decimal)_cijStrain; set => Set(ref _cijStrain, (double)Math.Clamp(value, 1e-6m, 0.01m), nameof(CijStrainD)); }
    public decimal CijConfigsD { get => _cijConfigs; set => Set(ref _cijConfigs, (int)Math.Max(1, value), nameof(CijConfigsD)); }
    public decimal FluctTD { get => (decimal)_fluctT; set => Set(ref _fluctT, (double)Math.Max(1, value), nameof(FluctTD)); }
    public decimal TensRateD { get => (decimal)_tensRate; set => Set(ref _tensRate, (double)Math.Max(1e-6m, value), nameof(TensRateD)); }
    public decimal TensMaxD { get => (decimal)_tensMax; set => Set(ref _tensMax, (double)Math.Clamp(value, 0.005m, 2m), nameof(TensMaxD)); }
    public decimal TensTD { get => (decimal)_tensT; set => Set(ref _tensT, (double)Math.Max(1, value), nameof(TensTD)); }
    public int TensAxis { get => _tensAxis; set => Set(ref _tensAxis, value); }
    public bool TensFixed { get => _tensFixed; set => Set(ref _tensFixed, value); }
    public static readonly string[] Axes = ["x", "y", "z", "x, y, z averaged"];
    public string TensRateText => string.Format(Inv, "{0:0.##e0} s⁻¹", _tensRate * 1e12);
    protected override void OnChanged(string? name) { if (name == nameof(TensRateD)) Raise(nameof(TensRateText)); }

    // pull test (interfaces): distance and rate travel in the tensile fields when the tensile run is off
    private double _pullDist = 10, _pullRate = 2, _pullT = 300, _pullEq = 5;
    public decimal PullDistD { get => (decimal)_pullDist; set => Set(ref _pullDist, (double)Math.Clamp(value, 0.5m, 200m), nameof(PullDistD)); }
    public decimal PullRateD { get => (decimal)_pullRate; set => Set(ref _pullRate, (double)Math.Clamp(value, 0.01m, 100m), nameof(PullRateD)); }
    public decimal PullTD { get => (decimal)_pullT; set => Set(ref _pullT, (double)Math.Max(1, value), nameof(PullTD)); }
    public decimal PullEqD { get => (decimal)_pullEq; set => Set(ref _pullEq, (double)Math.Max(0, value), nameof(PullEqD)); }
    public bool PullOn => PullShearChip.IsOn || PullNormalChip.IsOn;
    private int _pullAxis;
    /// <summary>Shear axis: x, y, or z for pull-out along a fibre.</summary>
    public int PullAxis { get => _pullAxis; set => Set(ref _pullAxis, value); }
    public static readonly string[] PullAxes = ["x (slab)", "y (slab)", "z (fibre pull-out)"];

    /// <summary>Seed of the runs (cooling scans, pulls): replicas differ only in it.</summary>
    public ulong MechSeed { get; set; } = 1;

    public CapsMechOpts MechOptions()
    {
        var pull = PullOn && !TensileChip.IsOn;
        return new CapsMechOpts
        {
            Configurations = _cijConfigs, Strain = _cijStrain,
            Temperature = pull ? _pullT : NvtRunOn && !TensileChip.IsOn ? _fluctT : TensileChip.IsOn ? _tensT : NemdChip.IsOn ? _shT : _fluctT,
            Axis = pull ? _pullAxis : _tensAxis, Rate = pull ? _pullRate : _tensRate, MaxStrain = pull ? _pullDist : _tensMax, LateralFixed = _tensFixed ? 1 : 0,
            TStart = _tgFrom, TEnd = _tgTo, TStep = _tgStep, PsPerStep = _tgPs, RunPs = _fluctPs,
            EquilibratePs = pull ? (_pullEq > 0 ? _pullEq : -1) : _eqPs > 0 ? _eqPs : -1,
            Seed = MechSeed,
            Pressure = _tgPressure, Barostat = _tgBarostat, TauT = _tgTauT, TauP = _tgTauP, AverageFrom = Math.Max(1e-6, 1 - _tgAverage / 100.0),
            TgProperty = _tgProperty, TgFit = _tgFit, GlassyMax = _tgGlassy, RubberyMin = _tgRubbery,
            FrMoving = _frMoving, FrFixed = _frFixed, FrVelocity = _frSpeed / 100, FrPs = _frPs, FrEqPs = _frEq > 0 ? _frEq : -1, FrT = _frT,
            CreepStress = _crStress, CreepT = _crT, CreepPs = _crPs, CreepEqPs = _crEq > 0 ? _crEq : -1, CreepAxis = _crAxis,
            ConfTrials = _confTrials, ConfMethod = _confMethod, ConfSelection = _confSel ? 1 : 0, ConfWindow = _confWindow, ConfRmsd = _confRmsd,
            ShearLo = _shLo, ShearHi = _shHi, ShearPoints = _shPoints, ShearPs = _shPs, ShearEqPs = _shEq > 0 ? _shEq : -1,
        };
    }

    /// <summary>The cooling scan's settings as a recipe's tg keys (core/src/recipe.cpp), the same scan as Run gives.</summary>
    public string TgRecipeKeys(ulong seed)
    {
        var inv = CultureInfo.InvariantCulture;
        string G(double v) => v.ToString("R", inv);
        var keys = new List<string>
        {
            "t_start: " + G(_tgFrom), "t_end: " + G(_tgTo), "t_step: " + G(_tgStep), "ps_per_step: " + G(_tgPs), "equilibrate_ps: " + G(_eqPs > 0 ? _eqPs : -1),
            "pressure: " + G(_tgPressure), "barostat: " + (_tgBarostat == 1 ? "berendsen" : _tgBarostat == 2 ? "mtk" : "crescale"),
            "tau_t: " + G(_tgTauT), "tau_p: " + G(_tgTauP), "average_from: " + G(Math.Max(1e-6, 1 - _tgAverage / 100.0)),
            "property: " + (_tgProperty == 1 ? "energy" : "volume"), "fit: " + (_tgFit == 1 ? "ranges" : "hinge"),
        };
        if (_tgFit == 1) { keys.Add("glassy_max: " + G(_tgGlassy)); keys.Add("rubbery_min: " + G(_tgRubbery)); }
        keys.Add("seed: " + seed.ToString(inv));
        return "{ " + string.Join(", ", keys) + " }";
    }

    public static readonly string[] Pairs = ["all – all", "C – C", "C – H", "H – H", "C – O", "C – N", "O – H"];
    private static readonly (int A, int B)[] PairElements = [(0, 0), (6, 6), (6, 1), (1, 1), (6, 8), (6, 7), (8, 1)];

    private string _sourceText = "No trajectory open";
    public string SourceText { get => _sourceText; private set => Set(ref _sourceText, value); }
    private string _framesText = "";
    public string FramesText { get => _framesText; private set => Set(ref _framesText, value); }

    /// <summary>Called when the document or frame count changes.</summary>
    public void OnDocument(string title, CapsSummary? s)
    {
        if (s is not CapsSummary x) { SourceText = "No trajectory open"; FramesText = ""; return; }
        SourceText = title;
        FramesText = x.Frames > 1 ? string.Format(Inv, "{0:N0} frames · {1:N0} atoms", x.Frames, x.Atoms) : string.Format(Inv, "one structure · {0:N0} atoms · dynamics need frames", x.Atoms);
    }

    // ---------------------------------------------------------------- references
    public ObservableCollection<RefMaterial> References { get; } = new();
    private int _refIndex;
    public int RefIndex { get => _refIndex; set { if (Set(ref _refIndex, value)) ApplyReferences(); } }
    private RefMaterial? SelectedRef => _refIndex > 0 && _refIndex < References.Count ? References[_refIndex] : null;

    private void LoadReferences()
    {
        References.Add(new RefMaterial("", "No comparison", new()));
        var f = Paths.References;
        if (f != null)
        {
            try
            {
                using var js = JsonDocument.Parse(File.ReadAllText(f));
                foreach (var m in js.RootElement.GetProperty("materials").EnumerateArray())
                {
                    var vals = new Dictionary<string, RefValue>();
                    foreach (var v in m.GetProperty("values").EnumerateObject())
                        vals[v.Name] = new RefValue(v.Value.GetProperty("lo").GetDouble(), v.Value.GetProperty("hi").GetDouble(), Str(v.Value, "unit"),
                            Str(v.Value, "condition"), Str(v.Value, "source"));
                    References.Add(new RefMaterial(Str(m, "id"), Str(m, "name"), vals));
                }
            }
            catch (Exception e) { _status("Cannot read " + f + ": " + e.Message); }
        }
    }

    private static string Str(JsonElement e, string k) => e.TryGetProperty(k, out var v) && v.ValueKind == JsonValueKind.String ? v.GetString()! : "";

    private void ApplyReferences()
    {
        var m = SelectedRef;
        var cards = Results.ToList();
        Results.Clear();
        foreach (var c in cards)
        {
            c.Ref = m != null && m.Values.TryGetValue(c.Id, out var r) ? r : null;
            Results.Add(c);
        }
    }

    // ---------------------------------------------------------------- run
    private bool _working;
    public bool Working { get => _working; private set { if (Set(ref _working, value)) { Raise(nameof(NotWorking)); _running(value); } } }
    public bool NotWorking => !_working;
    private string _log = "Choose calculations and run them on the open trajectory.";
    public string Log { get => _log; private set => Set(ref _log, value); }
    private double _progress;
    public double Progress { get => _progress; private set => Set(ref _progress, value); }
    private CancellationTokenSource? _cancel;

    public ObservableCollection<ResultCard> Results { get; } = new();
    public ObservableCollection<SeriesItem> Curves { get; } = new();
    private int _curveIndex = -1;
    public int CurveIndex { get => _curveIndex; set { if (Set(ref _curveIndex, value)) Raise(nameof(Curve)); } }
    public SeriesItem? Curve => _curveIndex >= 0 && _curveIndex < Curves.Count ? Curves[_curveIndex] : null;
    public bool HasResults => Results.Count > 0;
    private string _runInfo = "";
    public string RunInfo { get => _runInfo; private set => Set(ref _runInfo, value); }
    private string _json = "";

    public string[] SelectedIds => Groups.SelectMany(g => g.Chips).Where(c => c.IsOn && c.Available).Select(c => c.Id).ToArray();

    // Groups (design/boards/Analyze "Groups"): the atoms the properties see
    public static readonly string[] GroupChoices = ["All atoms", "The selection", "Without the held molecule", "Molecules…"];
    private int _group;
    private string _groupMolecules = "1";
    public int GroupIndex { get => _group; set { if (Set(ref _group, Math.Clamp(value, 0, 3))) Raise(nameof(GroupIsMolecules)); } }
    public bool GroupIsMolecules => _group == 3;
    public string GroupMolecules { get => _groupMolecules; set => Set(ref _groupMolecules, value ?? ""); }
    private string GroupSpec => _group switch { 1 => "selection", 2 => "exclude-held", 3 => "molecules:" + _groupMolecules.Trim(), _ => "" };

    // interfaces (z profile, adhesion, filler interaction) and Herman's f: the axis they run along and the surface's molecules
    public static readonly string[] AxisChoices = ["x (a)", "y (b)", "z (c)"];
    private int _axis = 2;
    private double _zbin = 0.5;
    private string _surfaceMolecules = "1";
    public int AxisIndex { get => _axis; set => Set(ref _axis, Math.Clamp(value, 0, 2)); }
    public string AxisName => "xyz"[_axis].ToString();
    public decimal ZBinD { get => (decimal)_zbin; set => Set(ref _zbin, (double)Math.Clamp(value, 0.05m, 5m), nameof(ZBinD)); }
    /// <summary>The surface or filler: molecule ids such as "1" or "1-3,7" (empty: molecule 1).</summary>
    public string SurfaceMolecules { get => _surfaceMolecules; set => Set(ref _surfaceMolecules, value ?? ""); }

    public CapsAnalyzeOpts Options()
    {
        var (a, b) = PairElements[Math.Clamp(_pair, 0, PairElements.Length - 1)];
        return new CapsAnalyzeOpts
        {
            First = _first, Last = _last < 0 ? -1 : _last, Stride = _stride, FramePs = _framePs, TimestepFs = _timestepFs, Blocks = 5,
            ElemA = a, ElemB = b, InterOnly = _inter ? 1 : 0, FitFrom = _fitFrom, FitTo = _fitTo, Probe = _probe, Grid = _grid,
            Qmax = _qmax, Dq = _dq, QDirect = _qDirect, Deuterate = _deuterate, Group = GroupSpec, Radii = _radii,
            ZBin = _zbin, Axis = _axis + 1, Surface = _surfaceMolecules.Trim(),
        };
    }

    /// <summary>The result cards again, so they show the display units chosen in Settings.</summary>
    public void RefreshDisplayUnits()
    {
        var cards = Results.ToList();
        Results.Clear();
        foreach (var c in cards) Results.Add(c);
    }

    public async Task Run()
    {
        var doc = _doc();
        if (doc == null || Working) return;
        var ids = SelectedIds;
        if (ids.Length == 0) { Log = "Switch on at least one calculation."; return; }
        Working = true;
        Progress = 0;
        _cancel = new CancellationTokenSource();
        var token = _cancel.Token;
        var o = Options();
        var sw = System.Diagnostics.Stopwatch.StartNew();
        Log = "Starting…";
        _status("Analyzing " + string.Join(", ", ids) + "…");
        try
        {
            if (FluctChip.IsOn && TensileChip.IsOn && Math.Abs(_fluctT - _tensT) > 1e-9)
                Log = "Note: the fluctuation constants and the tensile run share one temperature; the tensile temperature is used";
            var mo = MechOptions();
            var json = await Task.Run(() => doc.Analyze(string.Join(",", ids), o, mo, (what, f) =>
            {
                Avalonia.Threading.Dispatcher.UIThread.Post(() =>
                {
                    if (!_working) return;
                    Progress = f;
                    Log = string.Format(Inv, "{0} · {1:F0} % · {2:F0} s", what, 100 * f, sw.Elapsed.TotalSeconds);
                });
                return !token.IsCancellationRequested;
            }));
            Load(json);
            Log = string.Format(Inv, "{0} properties in {1:F1} s", Results.Count, sw.Elapsed.TotalSeconds);
            _status("Analysis finished · " + RunInfo);
        }
        catch (Exception e)
        {
            var cancelled = e.Message.Contains("cancelled");
            Log = cancelled ? "Cancelled." : "Could not analyse: " + e.Message;
            _status(cancelled ? "Analysis cancelled" : "Analysis failed — see the Analyze page");
        }
        finally { Working = false; Progress = 0; }
    }

    public void Cancel() => _cancel?.Cancel();

    /// <summary>Reads the JSON report of caps_analyze into cards and curves.</summary>
    public void Load(string json)
    {
        _json = json;
        Results.Clear();
        Curves.Clear();
        if (json.Length == 0) { Raise(nameof(HasResults)); return; }
        using var js = JsonDocument.Parse(json);
        var root = js.RootElement;
        RunInfo = string.Format(Inv, "{0} of {1} frames · {2:N0} atoms", root.GetProperty("frames").GetInt32(), root.GetProperty("of").GetInt32(), root.GetProperty("atoms").GetInt32());
        var m = SelectedRef;
        foreach (var p in root.GetProperty("properties").EnumerateArray())
        {
            var id = Str(p, "id");
            var extra = new List<(string, double)>();
            if (p.TryGetProperty("extra", out var ex))
                foreach (var e in ex.EnumerateObject())
                    extra.Add((e.Name, e.Value.ValueKind == JsonValueKind.Number ? e.Value.GetDouble() : double.NaN));
            var notes = p.TryGetProperty("notes", out var ns) ? ns.EnumerateArray().Select(n => n.GetString() ?? "").ToArray() : [];
            var card = new ResultCard
            {
                Id = id, Name = Str(p, "name"), Unit = Str(p, "unit"), Method = Str(p, "method"), Notes = notes, Extra = extra.ToArray(),
                Value = Dbl(p, "value"), Error = Dbl(p, "error"),
                Ref = m != null && m.Values.TryGetValue(id, out var r) ? r : null,
            };
            Results.Add(card);
            if (id == "viscosity" && double.IsFinite(card.Value) && card.Value > 0) ViscosityComputed?.Invoke(card.Value);
            if (p.TryGetProperty("series", out var ss))
            {
                var list = new List<SeriesItem>();
                foreach (var s in ss.EnumerateArray())
                {
                    var x = s.GetProperty("x").EnumerateArray().Select(v => v.ValueKind == JsonValueKind.Number ? v.GetDouble() : double.NaN).ToArray();
                    var y = s.GetProperty("y").EnumerateArray().Select(v => v.ValueKind == JsonValueKind.Number ? v.GetDouble() : double.NaN).ToArray();
                    double? refY = id is "rdf" or "sq" or "xray" or "electron" or "neutron" ? 1.0 : id is "relaxation" ? 0.0 : null;
                    list.Add(new SeriesItem(card.Name, Str(s, "label"), Str(s, "x_label"), Str(s, "y_label"), x, y, id == "msd", refY));
                }
                if (id == "tg" && list.Count >= 2)          // specific volume as points, the two-line fit through them
                    list = [list[0] with { Markers = true, OverlayX = list[1].X, OverlayY = list[1].Y, Label = "specific volume and two-line fit" }, .. list.Skip(2)];
                else if (id == "tensile_modulus" && list.Count >= 2)   // raw stress as points, the smoothed curve through them
                    list = [list[1] with { Markers = true, OverlayX = list[0].X, OverlayY = list[0].Y, Label = "stress–strain (points) and smoothed" }, .. list.Skip(2)];
                foreach (var it in list) Curves.Add(it);
            }
        }
        Raise(nameof(HasResults));
        CurveIndex = -1;
        CurveIndex = Curves.Count == 0 ? -1 : Math.Max(0, Curves.ToList().FindIndex(c => c.X.Length >= 2));
    }

    private static double Dbl(JsonElement e, string k) => e.TryGetProperty(k, out var v) && v.ValueKind == JsonValueKind.Number ? v.GetDouble() : double.NaN;

    // ---------------------------------------------------------------- export
    /// <summary>Writes results.csv (one row per property, with the comparison), one CSV per curve, results.tex (a booktabs
    /// table) and results.json into a folder. Returns the number of files.</summary>
    public int Export(string dir)
    {
        Directory.CreateDirectory(dir);
        var n = 0;
        var sb = new StringBuilder("id,property,value,error,unit,reference_lo,reference_hi,reference_source,method\n");
        foreach (var c in Results)
            sb.Append(Csv(c.Id)).Append(',').Append(Csv(c.Name)).Append(',').Append(c.HasValue ? c.Value.ToString("R", Inv) : "").Append(',')
              .Append(double.IsNaN(c.Error) ? "" : c.Error.ToString("R", Inv)).Append(',').Append(Csv(c.Unit)).Append(',')
              .Append(c.Ref?.Lo.ToString("R", Inv) ?? "").Append(',').Append(c.Ref?.Hi.ToString("R", Inv) ?? "").Append(',')
              .Append(Csv(c.Ref?.Source ?? "")).Append(',').Append(Csv(c.Method)).Append('\n');
        File.WriteAllText(Path.Combine(dir, "results.csv"), sb.ToString());
        n++;
        var used = new HashSet<string>();
        foreach (var s in Curves)
        {
            var name = Slug(s.Property + "_" + s.Label);
            while (!used.Add(name)) name += "_";
            var w = new StringBuilder();
            w.Append(Csv(s.XLabel)).Append(',').Append(Csv(s.YLabel)).Append('\n');
            for (var k = 0; k < s.X.Length; k++) w.Append(s.X[k].ToString("R", Inv)).Append(',').Append(s.Y[k].ToString("R", Inv)).Append('\n');
            File.WriteAllText(Path.Combine(dir, name + ".csv"), w.ToString());
            n++;
        }
        File.WriteAllText(Path.Combine(dir, "results.tex"), Latex());
        n++;
        if (_json.Length > 0) { File.WriteAllText(Path.Combine(dir, "results.json"), _json); n++; }
        return n;
    }

    public string Latex()
    {
        var sb = new StringBuilder();
        sb.Append("% CAPS Analyze · ").Append(SourceText).Append(" · ").Append(RunInfo).Append("\n");
        sb.Append("\\begin{tabular}{llll}\n\\toprule\nProperty & Simulation & Experiment & Unit \\\\\n\\midrule\n");
        foreach (var c in Results.Where(c => c.HasValue))
        {
            var sim = c.ValueText + (c.ErrorText.Length == 0 ? "" : " $\\pm$ " + c.ErrorText[2..]);
            var exp = c.Ref == null ? "--" : (c.Ref.Lo == c.Ref.Hi ? ResultCard.Num(c.Ref.Lo) : ResultCard.Num(c.Ref.Lo) + "--" + ResultCard.Num(c.Ref.Hi));
            sb.Append(Tex(c.Name)).Append(" & ").Append(sim).Append(" & ").Append(exp).Append(" & ").Append(Tex(c.Unit)).Append(" \\\\\n");
        }
        sb.Append("\\bottomrule\n\\end{tabular}\n");
        var sources = Results.Where(c => c.Ref != null).Select(c => c.Ref!.Source).Distinct().ToArray();
        if (sources.Length > 0) sb.Append("% experiment: ").Append(string.Join("; ", sources)).Append('\n');
        return sb.ToString();
    }

    private static string Csv(string s) => s.IndexOfAny([',', '"', '\n']) >= 0 ? "\"" + s.Replace("\"", "\"\"") + "\"" : s;
    private static string Slug(string s) => new(s.Select(ch => char.IsLetterOrDigit(ch) ? char.ToLowerInvariant(ch) : '_').ToArray());
    private static string Tex(string s) => s.Replace("\\", "\\textbackslash{}").Replace("&", "\\&").Replace("%", "\\%").Replace("_", "\\_")
        .Replace("#", "\\#").Replace("^", "\\^{}").Replace("³", "$^3$").Replace("²", "$^2$").Replace("⁻¹", "$^{-1}$").Replace("⁻⁵", "$^{-5}$")
        .Replace("½", "$^{1/2}$").Replace("√", "$\\sqrt{}$").Replace("⟨", "$\\langle$").Replace("⟩", "$\\rangle$").Replace("Å", "\\AA{}")
        .Replace("δ", "$\\delta$").Replace("∞", "$_\\infty$").Replace("τ", "$\\tau$");

    /// <summary>The frame a point of the shown curve belongs to (D9: click a chart point → that frame): curves over
    /// "frame" map directly; curves over time map through the analysis's frame times (frame_ps, or the timesteps × the
    /// time step) from the first analysed frame, as the core lays them out. −1 when the curve is not over the trajectory.</summary>
    public int FrameOfCurveX(double x)
    {
        if (_curveIndex < 0 || _curveIndex >= Curves.Count) return -1;
        var c = Curves[_curveIndex];
        var lbl = c.XLabel.ToLowerInvariant();
        var doc = _doc();
        if (doc == null) return -1;
        var steps = doc.FrameTimesteps();
        if (steps.Length == 0) return -1;
        if (lbl.StartsWith("frame")) return Math.Clamp((int)Math.Round(x), 0, steps.Length - 1);
        if (!lbl.Contains("ps")) return -1;
        double T(int k) => _framePs > 0 ? k * _framePs : steps[k] * _timestepFs * 1e-3;
        var first = Math.Clamp(_first, 0, steps.Length - 1);
        var t = x + T(first);
        int best = -1;
        var bd = double.MaxValue;
        for (var k = first; k < steps.Length; ++k)
        {
            var d = Math.Abs(T(k) - t);
            if (d < bd) { bd = d; best = k; }
        }
        return best;
    }
}
