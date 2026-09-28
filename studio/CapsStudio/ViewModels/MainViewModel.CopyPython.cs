using System.Globalization;
using System.Text;

namespace CapsStudio.ViewModels;

/// <summary>Copy as Python (design/boards Relax, Dynamics, React): the page's run on the open structure as a script for
/// the caps package — the same settings, runnable with python3 on this or another machine.</summary>
public sealed partial class MainViewModel
{

    /// <summary>import caps, the open structure (its file; a note when it has not been saved) and its force field.</summary>
    private StringBuilder PythonStart(string what)
    {
        var sb = new StringBuilder($"# {what} · written by CAPS Studio · run with python3 (the caps package)\nimport caps\n\n");
        var path = _doc?.Path is { Length: > 0 } p && File.Exists(p) && !Title.Contains("(unsaved)") ? p : null;
        if (path == null) sb.Append("# the structure is not saved yet: File › Save, then put its path here\n");
        sb.Append($"doc = caps.open({PyStr(path ?? "structure.data")}{(path != null && TopologyFor(path) is { Length: > 0 } topo ? ", topology=" + PyStr(topo) : "")})\n");
        if (Field.Assigned && Field.FfIndex >= 0 && Field.FfIndex < Field.Library.Count)
            sb.Append($"doc.field.assign({PyStr(Field.Library[Field.FfIndex].Id)})\n");
        return sb;
    }

    private static string Py(bool b) => b ? "True" : "False";
    private static string Num(double x) => x.ToString("0.######", Inv);

    public string RelaxPython()
    {
        var sb = PythonStart("Relax");
        var method = _relaxMethod switch { 0 => "sd", 1 => "cg", 3 => "fire", _ => "lbfgs" };
        sb.Append($"rc = doc.relax(ftol={Num(_relaxFtol)}, method=\"{method}\", max_iterations={_relaxIterations}, pushoff={Py(_relaxPushoff)}");
        if (_relaxCompress) sb.Append($", density={Num(_relaxDensity)}");
        if (_relaxBox) sb.Append(_relaxBoxMode switch { 2 => ", box_axes=\"z\"", 3 => ", box_axes=\"xy\"", 1 => ", box_axes=\"xyz\"", _ => ", box=True" }).Append($", pressure={Num(_relaxPressure)}");
        if (_relaxPushoff && _relaxPushoffMd) sb.Append($", pushoff_md_ps={Num((double)_relaxRampPs)}, pushoff_cap={Num((double)_relaxCap)}, pushoff_temperature={Num((double)_relaxPushoffT)}");
        sb.Append($", cutoff={Num(_relaxCutoff)}, coulomb={Py(_relaxCoulomb)})\n");
        sb.Append("print(\"converged\" if rc == 0 else \"stopped before the tolerance\")\nprint(doc.report)\ndoc.save(\"relaxed.data\")\n");
        return sb.ToString();
    }

    public string MdPython()
    {
        var sb = PythonStart("Dynamics");
        string[] th = ["bussi", "langevin", "nose-hoover"], ba = ["crescale", "berendsen", "mtk"], cons = ["none", "h-bonds", "all-bonds"], solver = ["shake", "lincs"];
        var thermostat = _mdEnsemble is 1 or 2 ? th[Math.Clamp(_mdThermostat, 0, 2)] : "none";
        var barostat = _mdEnsemble == 2 ? ba[Math.Clamp(_mdBarostat, 0, 2)] : _mdEnsemble == 3 ? "berendsen" : "none";
        sb.Append($"doc.md(steps={_mdSteps}, dt={Num(_mdDt)}, temperature={Num(_mdTemp)}, thermostat=\"{thermostat}\", barostat=\"{barostat}\", pressure={Num(_mdPressure)}, " +
                  $"seed={_mdSeed}, frame_every={_mdFrameEvery}, cutoff={Num(_relaxCutoff)}");
        if (RespaSteps > 1) sb.Append($", respa={RespaSteps}");
        if (_mdConstraints > 0) sb.Append($", constraints=\"{cons[Math.Clamp(_mdConstraints, 0, 2)]}\", constraint_solver=\"{solver[Math.Clamp(_mdConstraintSolver, 0, 1)]}\"");
        sb.Append(")\nprint(doc.report)\ndoc.save(\"final.data\")\ndoc.save_trajectory(\"trajectory.lammpstrj\")\n");
        return sb.ToString();
    }

    /// <summary>The Pack page's input as a script: the packmol text written beside the script, packed by caps.pack.</summary>
    public string PackPython()
    {
        string text;
        try { text = PackTextToRun(); } catch (Exception e) { return "# " + e.Message + "\n"; }
        var sb = new StringBuilder("# Pack · written by CAPS Studio · run with python3 (the caps package)\nimport caps\n\n");
        sb.Append("packmol_input = \"\"\"\n").Append(text.Replace("\"\"\"", "\\\"\\\"\\\"").TrimEnd()).Append("\n\"\"\"\n");
        sb.Append("with open(\"pack.inp\", \"w\") as f:\n    f.write(packmol_input)\n");
        sb.Append($"doc = caps.pack(inp=\"pack.inp\", base_dir={PyStr(PackBaseDir)})   # the structure files are found here\n");
        sb.Append("print(doc.atoms, \"atoms\")\ndoc.save(\"packed.data\")\n");
        return sb.ToString();
    }

    public string ReactPython()
    {
        var sb = PythonStart("React");
        sb.Append("templates = \"\"\"\n").Append(_rxText.Replace("\"\"\"", "\\\"\\\"\\\"").TrimEnd()).Append("\n\"\"\"\n");
        sb.Append($"report = doc.react(templates, cycles={_rxCycles}, per_cycle={_rxPerCycle}, target={Num(_rxTarget)}, capture={Num(_rxCapture)}, relax={Py(_rxRelax)}, " +
                  $"relax_iterations={_rxRelaxIt}, md_ps={Num(_rxMdPs)}, temperature={Num(_rxTemp)}, cutoff={Num(_relaxCutoff)}, seed={_rxSeed})\n");
        sb.Append("print(report)\ndoc.save(\"network.data\")\n");
        return sb.ToString();
    }
}
