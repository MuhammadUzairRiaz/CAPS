using System.Globalization;
using System.Text;
using CapsStudio.Interop;

namespace CapsStudio.ViewModels;

/// <summary>Copy as command: a page's run on the open structure as the caps command line runs it — the same option names
/// and defaults, the same engine (the CLI calls the core, the Studio the C API over the same core functions), so a run here
/// and the copied command give the same result. cli/gui_sync.json maps every CLI option to its Studio control and the
/// test cli_gui_sync keeps the two in step.</summary>
public sealed partial class MainViewModel
{
    /// <summary>The open structure's file when it is saved, else null.</summary>
    public string? DocumentPath => _doc?.Path is { Length: > 0 } p && File.Exists(p) && !Title.Contains("(unsaved)") ? p : null;

    /// <summary>The structure as a command-line argument (its file, with --topology when it has one).</summary>
    public string CliStructure
    {
        get
        {
            var p = DocumentPath;
            if (p == null) return "structure.data";
            var topo = TopologyFor(p);
            return FieldViewModel.Quote(p) + (topo is { Length: > 0 } ? " --topology " + FieldViewModel.Quote(topo) : "");
        }
    }

    private static string C(double x) => x.ToString("0.######", CultureInfo.InvariantCulture);

    /// <summary>The start of a copied command: what wrote it and a note when the structure is not saved yet.</summary>
    private StringBuilder CliStart(string what)
    {
        var sb = new StringBuilder($"# {what} · written by CAPS Studio · the same run with the caps command line\n");
        if (DocumentPath == null) sb.Append("# the structure is not saved yet: File › Save, then put its path for structure.data\n");
        return sb;
    }

    /// <summary>The force field part: --ff / --overlay / --charges, or a note when the assignment has no command-line form.</summary>
    private string CliField(StringBuilder sb)
    {
        var f = Field.CliArgs();
        if (f == null)
        {
            sb.Append("# the force field is assigned by group (or a water model, or the file's own): Python's doc.field.save(\"structure.ff.json\") writes it whole; pass --ff structure.ff.json\n");
            return " --ff structure.ff.json";
        }
        return f;
    }

    public string RelaxCommand()
    {
        var sb = CliStart("Relax");
        var method = _relaxMethod switch { 0 => "sd", 1 => "cg", 3 => "fire", _ => "lbfgs" };
        var ff = CliField(sb);
        sb.Append($"caps relax {CliStructure} -o relaxed.data --method {method} --ftol {C(_relaxFtol)} --iterations {_relaxIterations}");
        if (_relaxCompress) sb.Append($" --density {C(_relaxDensity)} --step {C(_relaxStep)}");
        if (_relaxBox) sb.Append($" --box-relax --pressure {C(_relaxPressure)}");
        if (!_relaxPushoff) sb.Append(" --no-pushoff");
        sb.Append($" --cutoff {C(_relaxCutoff)}");
        if (!_relaxCoulomb) sb.Append(" --no-coulomb");
        sb.Append(ff).Append('\n');
        return sb.ToString();
    }

    public string MdCommand()
    {
        var sb = CliStart("Dynamics");
        var ff = CliField(sb);
        string[] th = ["bussi", "langevin", "nose-hoover"], ba = ["crescale", "berendsen", "mtk"], cons = ["none", "h-bonds", "all-bonds"], solver = ["shake", "lincs"];
        var thermostat = _mdEnsemble is 1 or 2 ? th[Math.Clamp(_mdThermostat, 0, 2)] : "none";
        var barostat = _mdEnsemble == 2 ? ba[Math.Clamp(_mdBarostat, 0, 2)] : _mdEnsemble == 3 ? "berendsen" : "none";
        sb.Append($"caps md {CliStructure} -o final.data --steps {_mdSteps} --dt {C(_mdDt)} --temp {C(_mdTemp)} --thermostat {thermostat} --tau-t {C(_mdTauT)}" +
                  $" --barostat {barostat} --pressure {C(_mdPressure)} --tau-p {C(_mdTauP)} --seed {_mdSeed} --dump trajectory.lammpstrj --every {_mdFrameEvery}" +
                  $" --cutoff {C(_relaxCutoff)}");
        if (_mdConstraints > 0) sb.Append($" --constraints {cons[Math.Clamp(_mdConstraints, 0, 2)]} --constraint-solver {solver[Math.Clamp(_mdConstraintSolver, 0, 1)]}");
        if (_mdHmr) sb.Append(" --hmr 3.024");
        if (_mdNewVelocities) sb.Append(" --new-velocities");
        if (!_relaxCoulomb) sb.Append(" --no-coulomb");
        if (_mdDouble) sb.Append(" --precision double");
        if (RespaSteps > 1) sb.Insert(0, "# r-RESPA has no command-line option: the command runs every force at each step\n");
        sb.Append(ff).Append('\n');
        return sb.ToString();
    }

    public string EqCommand()
    {
        var sb = CliStart("Equilibrate");
        var ff = CliField(sb);
        var name = ProtocolNames[Math.Clamp(_eqProtocol, 0, ProtocolNames.Length - 1)];
        if (name.Length == 0)
        {   // a protocol written by hand: its text beside the command
            sb.Append("cat > protocol.txt <<'CAPS_PROTOCOL'\n").Append(EqText.TrimEnd()).Append("\nCAPS_PROTOCOL\n");
            name = "protocol.txt";
        }
        string[] th = ["bussi", "langevin", "nose-hoover"], ba = ["crescale", "berendsen", "mtk"], cons = ["none", "h-bonds", "all-bonds"];
        sb.Append($"caps equilibrate {CliStructure} -o equilibrated.data --protocol {name}");
        if (name == "larsen21") sb.Append($" --tfinal {C(_eqTFinal)} --tmax {C(_eqTMax)} --pfinal {C(_eqPFinal)} --pmax {C(_eqPMax)}");
        if (name == "annealing") sb.Append($" --cycles {_eqCycles} --tlow {C(_eqTLow)} --thigh {C(_eqTHigh)} --ramp {C(_eqRamp)} --hold {C(_eqHold)} --pfinal {C(_eqPFinal)}");
        if (_eqScale != 1 && name != "protocol.txt") sb.Append($" --scale {C(_eqScale)}");
        if (_eqUntil) sb.Append($" --until-converged --block {C(_eqBlock)} --max-blocks {_eqMaxBlocks}");
        sb.Append($" --dt {C(_mdDt)} --thermostat {th[Math.Clamp(_mdThermostat, 0, 2)]} --tau-t {C(_mdTauT)} --barostat {ba[Math.Clamp(_mdBarostat, 0, 2)]} --tau-p {C(_mdTauP)} --seed {_mdSeed}" +
                  $" --dump equilibration.lammpstrj --every-ps 10 --log thermo.csv --cutoff {C(_relaxCutoff)}");
        if (_mdConstraints > 0) sb.Append($" --constraints {cons[Math.Clamp(_mdConstraints, 0, 2)]}");
        if (_mdHmr) sb.Append(" --hmr 3.024");
        if (_mdDouble) sb.Append(" --precision double");
        sb.Append(ff).Append('\n');
        return sb.ToString();
    }

    /// <summary>The React page as caps react: the templates (built-in names, else the text beside the command), every
    /// option of the page under the CLI's name.</summary>
    public string ReactCommand()
    {
        var sb = CliStart("React");
        // built-in templates by name; any other text goes to a file beside the command
        var names = BuiltinTemplateNames(_rxText);
        string templates;
        if (names != null) templates = string.Concat(names.Select(n => " --template " + n));
        else
        {
            sb.Append("cat > reactions.txt <<'CAPS_TEMPLATES'\n").Append(_rxText.TrimEnd()).Append("\nCAPS_TEMPLATES\n");
            templates = " --template reactions.txt";
        }
        var ff = _rxUseField ? CliField(sb) : "";
        sb.Append($"caps react {CliStructure} -o network.data{templates} --cycles {_rxCycles} --per-cycle {_rxPerCycle} --seed {_rxSeed} --temp {C(_rxTemp)}");
        if (_rxTargetKind == 0) sb.Append($" --target {C(_rxTarget)}");
        else sb.Append($" --crosslinks {new[] { "", "links", "per-chain", "density", "mc", "degree" }[Math.Clamp(_rxTargetKind, 1, 5)]}:{C(_rxTargetValue)}");
        if (_rxCapture > 0) sb.Append($" --capture {C(_rxCapture)}");
        if (!_rxRelax) sb.Append(" --no-relax");
        else sb.Append($" --relax-iterations {_rxRelaxIt}");
        if ((_rxRelax || _rxDuringMd) && _rxMdPs > 0) sb.Append($" --md-ps {C(_rxMdPs)}");
        if (_rxDuringMd) sb.Append(" --during-md");
        if (_rxBetween) sb.Append(" --between-chains");
        if (_rxKeepBy) sb.Append(" --keep-byproducts");
        if (_rxRemoveUnreacted) sb.Append(" --remove-unreacted");
        if (_rxAutoCapture) sb.Append($" --auto-capture --capture-max {C(_rxCaptureMax)}");
        if ((int)_rxSitesPer > 0) sb.Append($" --sites-per-chain {(int)_rxSitesPer}");
        if (RxSeveral && _rxByWeights) sb.Append(" --weights " + string.Join(",", RxReactions.Select(r => r.WeightD.ToString(CultureInfo.InvariantCulture))));
        sb.Append($" --cutoff {C(_relaxCutoff)}");
        if (!_relaxCoulomb) sb.Append(" --no-coulomb");
        sb.Append(ff).Append(" --summary network.json\n");
        return sb.ToString();
    }

    /// <summary>The built-in template names whose texts make up the text exactly (in order), else null.</summary>
    private static List<string>? BuiltinTemplateNames(string text)
    {
        static string Norm(string t) => string.Join("\n", t.Replace("\r", "").Split('\n').Select(l => l.TrimEnd()).Where(l => l.Length > 0));
        var want = Norm(text);
        if (want.Length == 0) return null;
        List<(string Name, string Text)> all;
        try
        {
            all = CapsDocument.ReactionTemplate("").Split('\n', StringSplitOptions.RemoveEmptyEntries | StringSplitOptions.TrimEntries)
                .Select(n => (n, Norm(CapsDocument.ReactionTemplate(n)))).Where(x => x.Item2.Length > 0).ToList();
        }
        catch { return null; }
        var names = new List<string>();
        var rest = want;
        while (rest.Length > 0)
        {
            var hit = all.FirstOrDefault(x => rest == x.Text || rest.StartsWith(x.Text + "\n", StringComparison.Ordinal));
            if (hit.Name == null) return null;
            names.Add(hit.Name);
            rest = rest.Length == hit.Text.Length ? "" : rest[(hit.Text.Length + 1)..];
        }
        return names;
    }

    /// <summary>Tensile and Tg (Analyze › Mechanics, Glass transition) as caps tensile / caps tg.</summary>
    public string TensileCommand()
    {
        var a = Analyze;
        var sb = CliStart("Tensile");
        var ff = CliField(sb);
        sb.Append($"caps tensile {CliStructure} -o tensile.data --axis {(a.TensAxis >= 3 ? "xyz" : AxisName(a.TensAxis))} --rate {N(a.TensRateD)} --strain {N(a.TensMaxD)} " +
                  $"--temp {N(a.TensTD)}{(a.TensFixed ? " --fixed-lateral" : "")}{(a.MechDouble ? " --precision double" : "")} --csv curves{ff}\n");
        return sb.ToString();
    }

    public string TgCommand()
    {
        var a = Analyze;
        var sb = CliStart("Glass transition");
        var ff = CliField(sb);
        sb.Append($"caps tg {CliStructure} -o cooled.data --from {N(a.TgFromD)} --to {N(a.TgToD)} --step {N(a.TgStepD)} --ps {N(a.TgPsD)}" +
                  $"{(a.MechDouble ? " --precision double" : "")} --csv tg{ff}\n");
        return sb.ToString();
    }

    // ---- React: unreacted sulfur donors leave the cell afterwards (caps react --remove-unreacted; the core's)
    private bool _rxRemoveUnreacted;
    public bool RxRemoveUnreacted { get => _rxRemoveUnreacted; set => Set(ref _rxRemoveUnreacted, value); }

    // ---- MD and equilibration in double precision (caps md / equilibrate --precision double; mixed is the default)
    private bool _mdDouble;
    public bool MdDoublePrecision { get => _mdDouble; set { if (Set(ref _mdDouble, value)) RefreshPreflight(); } }
}
