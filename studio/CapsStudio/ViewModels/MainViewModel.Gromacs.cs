using System.Globalization;

namespace CapsStudio.ViewModels;

/// <summary>Dynamics › Export to other engines › GROMACS: the topology, coordinates and run parameters for the same force
/// field and settings. The .mdp non-bonded section and the notes come from the core (checked against GROMACS, energies
/// and forces: bench/ff/check_gromacs.py); the run section is written here from the Dynamics settings.</summary>
public sealed partial class MainViewModel
{
    // the engine the deck is written for: LAMMPS, GROMACS, or both into one folder (the preview shows the LAMMPS input)
    private bool _mdGromacs, _mdBoth;
    public bool MdGromacs
    {
        get => _mdGromacs && !_mdBoth;
        set { if (value) { _mdBoth = false; } if (Set(ref _mdGromacs, value) || value) RaiseEngine(); }
    }
    public bool MdLammps { get => !_mdGromacs && !_mdBoth; set { if (value) { _mdBoth = false; _mdGromacs = false; RaiseEngine(); } } }
    public bool MdBoth { get => _mdBoth; set { if (Set(ref _mdBoth, value) && value) { _mdGromacs = false; RaiseEngine(); } } }
    private void RaiseEngine()
    {
        Raise(nameof(MdGromacs)); Raise(nameof(MdLammps)); Raise(nameof(MdBoth)); Raise(nameof(MdDeckFiles)); Raise(nameof(MdParityTip));
        RefreshPreflight();
    }
    public string MdDeckFiles => _mdBoth ? "system.data and system.in for LAMMPS, system.top, system.gro and system.mdp for GROMACS (one folder)"
        : _mdGromacs ? "system.top, system.gro and system.mdp (gmx grompp, then gmx mdrun)" : "system.data and system.in (lmp -in system.in)";
    public string MdParityTip => _mdGromacs
        ? "Energies and forces of the topology match CAPS in GROMACS 2026 (bench/ff/check_gromacs.py): bonded terms, Lennard-Jones and PME Coulomb; the tail correction differs by definition"
        : "Energies and forces of the data file match CAPS in LAMMPS (bench/ff/check_data_lammps.py)";

    // ---------------------------------------------------------------- Compare energies (design/boards/Dynamics)
    // The data file and input CAPS writes, run for zero steps in the LAMMPS on this machine, every energy term set beside
    // CAPS's own for the same frame: the parity the chip claims, checked for this structure and force field.
    private string _parity = "", _parityTable = "";
    private int _parityLevel;   // 0 not checked, 1 matches, 2 differs, 3 could not run
    private bool _parityRunning;
    public string ParityText => _parityLevel switch { 1 => "energies match LAMMPS", 2 => "energies differ", 3 => "not checked", _ => "compare energies" };
    public string ParityDetail { get => _parity; private set => Set(ref _parity, value); }
    public string ParityTable { get => _parityTable; private set { if (Set(ref _parityTable, value)) Raise(nameof(HasParityTable)); } }
    public bool HasParityTable => _parityTable.Length > 0;
    public bool ParityOk => _parityLevel == 1;
    public bool ParityBad => _parityLevel is 2 or 3;
    public bool ParityIdle => !_parityRunning;

    /// <summary>The LAMMPS executable: LAMMPS_EXE, else lmp / lmp_serial / lmp_mpi on the PATH or in ~/.local/bin.</summary>
    private static string? FindLammps()
    {
        if (Environment.GetEnvironmentVariable("LAMMPS_EXE") is { Length: > 0 } e && File.Exists(e)) return e;
        var dirs = (Environment.GetEnvironmentVariable("PATH") ?? "").Split(Path.PathSeparator)
                   .Append(Path.Combine(Environment.GetFolderPath(Environment.SpecialFolder.UserProfile), ".local", "bin"));
        foreach (var d in dirs)
            foreach (var n in new[] { "lmp", "lmp_serial", "lmp_mpi" })
                foreach (var ext in OperatingSystem.IsWindows() ? new[] { ".exe", "" } : [""])
                    if (File.Exists(Path.Combine(d, n + ext))) return Path.Combine(d, n + ext);
        return null;
    }

    public async Task CompareEnergies()
    {
        if (_doc == null || _parityRunning) return;
        var lmp = FindLammps();
        void Done(int level, string detail, string table)
        {
            _parityLevel = level;
            ParityDetail = detail;
            ParityTable = table;
            foreach (var n in new[] { nameof(ParityText), nameof(ParityOk), nameof(ParityBad) }) Raise(n);
        }
        if (lmp == null) { Done(3, "LAMMPS was not found (lmp on the PATH, or set LAMMPS_EXE): nothing to compare with", ""); return; }
        _parityRunning = true;
        Raise(nameof(ParityIdle));
        var doc = _doc;
        var dir = Path.Combine(Path.GetTempPath(), $"caps-parity-{Environment.ProcessId}-{DateTime.Now:HHmmssfff}");
        try
        {
            var (ok, caps, lmpTerms, log) = await Task.Run(() =>
            {
                Directory.CreateDirectory(dir);
                doc.Save(Path.Combine(dir, "system.data"));
                var deck = doc.LammpsInput("system.data");
                // a many-body potential (Tersoff, AIREBO …): LAMMPS's van der Waals includes it, CAPS does not evaluate it —
                // its energy is computed apart and taken out of the comparison
                var mbStyle = deck.Split('\n').Select(l => l.Split(' ', StringSplitOptions.RemoveEmptyEntries))
                                  .FirstOrDefault(w => w.Length > 4 && w[0] == "pair_coeff" && w[1] == "*" && w[2] == "*" && !w[3].StartsWith("coul", StringComparison.Ordinal))?[3];
                File.WriteAllText(Path.Combine(dir, "system.in"), deck + (mbStyle != null ? $"\ncompute mb all pair {mbStyle}" : "") +
                    "\nthermo_style custom step pe ebond eangle edihed eimp evdwl ecoul elong press" + (mbStyle != null ? " c_mb" : "") + "\nthermo_modify format float %.10f\nrun 0\n");
                // metal units (eV): back to kcal/mol by LAMMPS's own factor
                var metalScale = deck.Contains("\nunits           metal", StringComparison.Ordinal) ? 23.060549 : 1.0;
                var caps = System.Text.Json.Nodes.JsonNode.Parse(doc.EnergyTerms())!;
                var psi = new System.Diagnostics.ProcessStartInfo(lmp) { WorkingDirectory = dir, RedirectStandardOutput = true, RedirectStandardError = true, UseShellExecute = false };
                foreach (var a in new[] { "-in", "system.in", "-log", "log.lammps", "-screen", "none" }) psi.ArgumentList.Add(a);
                using var p = System.Diagnostics.Process.Start(psi)!;
                var err = p.StandardError.ReadToEndAsync();
                p.StandardOutput.ReadToEnd();
                if (!p.WaitForExit(600000)) { try { p.Kill(true); } catch { } return (false, caps, (double[]?)null, "LAMMPS took longer than 10 minutes"); }
                var logText = File.Exists(Path.Combine(dir, "log.lammps")) ? File.ReadAllText(Path.Combine(dir, "log.lammps")) : "";
                // the thermo line after "Step PotEng E_bond E_angle E_dihed E_impro E_vdwl E_coul E_long Press"
                var lines = logText.Split('\n');
                var h = Array.FindIndex(lines, l => l.TrimStart().StartsWith("Step") && l.Contains("PotEng") && l.Contains("E_vdwl"));
                if (p.ExitCode != 0 || h < 0 || h + 1 >= lines.Length)
                {
                    var why = lines.LastOrDefault(l => l.StartsWith("ERROR")) ?? err.Result.Split('\n').FirstOrDefault(l => l.Length > 0) ?? $"exit {p.ExitCode}";
                    return (false, caps, (double[]?)null, why.Trim());
                }
                var v = lines[h + 1].Split(' ', StringSplitOptions.RemoveEmptyEntries).Select(x => double.TryParse(x, NumberStyles.Float, CultureInfo.InvariantCulture, out var d) ? d : double.NaN).ToArray();
                if (v.Length >= 10)
                {
                    for (var k = 1; k <= 8; ++k) v[k] *= metalScale;
                    if (mbStyle != null && v.Length >= 11)
                    {
                        var mb = v[10] * metalScale;
                        v[6] -= mb;   // van der Waals and the total without the many-body energy CAPS does not compute
                        v[1] -= mb;
                    }
                }
                return (true, caps, v, (mbStyle != null ? $"{mbStyle} taken out of LAMMPS's energy (CAPS does not evaluate it)" : "") +
                                       (metalScale != 1 ? (mbStyle != null ? "; " : "") + "LAMMPS ran in metal units, converted to kcal/mol" : ""));
            });
            if (!ok || lmpTerms == null || lmpTerms.Length < 9) { Done(3, "LAMMPS did not run the deck: " + log, ""); return; }
            double C(string k) => (double?)caps[k] ?? double.NaN;
            var pme = (string?)caps["electrostatics"] == "pme";
            var rows = new (string Name, double Caps, double Lmp, bool Compare)[]
            {
                ("bonds", C("bond"), lmpTerms[2], true), ("angles", C("angle"), lmpTerms[3], true), ("dihedrals", C("dihedral"), lmpTerms[4], true),
                ("impropers", C("improper"), lmpTerms[5], true), ("van der Waals", C("vdw"), lmpTerms[6], true),
                ("Coulomb", C("coulomb"), lmpTerms[7] + lmpTerms[8], !pme), ("total", C("total"), lmpTerms[1], !pme),
            };
            var inv = CultureInfo.InvariantCulture;
            var worst = 0.0;
            var sb = new System.Text.StringBuilder("term              CAPS            LAMMPS          Δ\n");
            foreach (var r in rows)
            {
                var d = r.Caps - r.Lmp;
                var rel = Math.Abs(d) / Math.Max(1.0, Math.Abs(r.Lmp));
                if (r.Compare) worst = Math.Max(worst, rel);
                sb.Append(string.Format(inv, "{0,-16}{1,16:F6}{2,16:F6}{3,14:E2}{4}\n", r.Name, r.Caps, r.Lmp, d, r.Compare ? "" : "  (PME here, PPPM there)"));
            }
            if (log.Length > 0) sb.Append(log + "\n");
            var match = worst < 1e-4;
            Done(match ? 1 : 2, match ? $"Every term within {worst:0.0e0} (relative) of LAMMPS {Path.GetFileName(lmp)}, kcal/mol, for this frame and force field"
                                      : $"The largest difference is {worst:0.0e0} relative: see the table (a style LAMMPS computes differently, or a tail or cutoff setting)", sb.ToString());
        }
        catch (Exception e) { Done(3, "Could not compare: " + e.Message, ""); }
        finally
        {
            _parityRunning = false;
            Raise(nameof(ParityIdle));
            try { Directory.Delete(dir, true); } catch { }
        }
    }

    /// <summary>The GROMACS .mdp for this run: the core's non-bonded settings, then the integrator and coupling.</summary>
    private string GromacsDeck(Interop.CapsDocument doc)
    {
        var inv = CultureInfo.InvariantCulture;
        var nb = doc.Gromacs(null);
        var respa = RespaSteps > 1 && !(MdHasThermostat && _mdThermostat == 1);
        var dtPs = _mdDt / 1000 / (respa ? RespaSteps : 1);
        var steps = (long)_mdSteps * (respa ? RespaSteps : 1);
        var every = (long)_mdFrameEvery * (respa ? RespaSteps : 1);
        var sb = new System.Text.StringBuilder();
        sb.Append("; GROMACS run parameters written by CAPS Studio: the same force field and settings as this Dynamics run\n");
        sb.Append(MdHasThermostat && _mdThermostat == 1 ? "integrator               = sd          ; Langevin\n" : "integrator               = md\n");
        sb.Append(string.Format(inv, "dt                       = {0:0.######}\nnsteps                   = {1}\n", dtPs, steps));
        if (respa)
            sb.Append(string.Format(inv, "mts                      = yes         ; r-RESPA: bonded forces every step, non-bonded every {0}\nmts-levels               = 2\nmts-level2-forces        = longrange-nonbonded nonbonded pair\nmts-level2-factor        = {0}\n", RespaSteps));
        sb.Append(string.Format(inv, "nstxout-compressed       = {0}\nnstenergy                = {1}\nnstlog                   = {1}\n", every, Math.Max(1, every / 10)));
        if (_mdConstraints > 0)
        {   // the core's settings say every bond is flexible: this run holds some
            nb = string.Join("\n", nb.Split('\n').Where(l => !l.StartsWith("constraints ", StringComparison.Ordinal)));
            sb.Append(_mdConstraints == 1 ? "constraints              = h-bonds     ; bonds to hydrogen, as this run\n"
                                          : "constraints              = all-bonds   ; as this run\n");
            sb.Append(_mdConstraintSolver == 1 ? "constraint-algorithm     = lincs\nlincs-order              = 4\n"
                                               : "constraint-algorithm     = shake       ; as this run (GROMACS allows it without domain decomposition)\n");
        }
        sb.Append(nb);
        if (_mdFieldOn)   // E0 (V/nm) omega t0 sigma: a static field when omega and sigma are 0; 1 V/Å = 10 V/nm
            foreach (var (ax, e) in new[] { ("x", _mdEx), ("y", _mdEy), ("z", _mdEz) })
                if (e != 0) sb.Append(string.Format(inv, "electric-field-{0}         = {1:0.######} 0 0 0   ; V/nm, static, as this run\n", ax, e * 10));
        // coupling: τ in ps, pressure in bar
        if (MdHasThermostat)
        {
            sb.Append(_mdThermostat == 1 ? "; tcoupl is implied by the sd integrator\n"
                    : _mdThermostat == 2 ? "tcoupl                   = Nose-Hoover\nnh-chain-length          = 3           ; as this run\n"
                    : "tcoupl                   = V-rescale   ; Bussi\n");
            sb.Append(string.Format(inv, "tc-grps                  = System\ntau-t                    = {0:0.####}\nref-t                    = {1:0.##}\n", _mdTauT / 1000, _mdTemp));
        }
        else sb.Append("tcoupl                   = no\n");
        if (MdHasBarostat)
        {
            var berendsen = _mdEnsemble == 3 || _mdBarostat == 1;
            var mttk = !berendsen && _mdBarostat == 2;
            sb.Append(string.Format(inv, "pcoupl                   = {0}\npcoupltype               = isotropic\ntau-p                    = {1:0.####}\nref-p                    = {2:0.#####}\ncompressibility          = 4.5e-5\n",
                berendsen ? "Berendsen" : mttk ? "MTTK" : "C-rescale", _mdTauP / 1000, _mdPressure * 1.01325));
            if (mttk) sb.Append("nstcalcenergy            = 1           ; MTTK needs the energies every step\n");
        }
        else sb.Append("pcoupl                   = no\n");
        sb.Append(string.Format(inv, "gen-vel                  = yes         ; {0}\ngen-temp                 = {1:0.##}\ngen-seed                 = {2}\n",
            _mdNewVelocities ? "new velocities, as this run" : "the .gro carries no velocities, so they are drawn at the target", _mdTemp, _mdSeed));
        return sb.ToString();
    }

    /// <summary>Writes the GROMACS files into a folder: system.top and system.gro from the core, system.mdp with the run.</summary>
    public void SaveGromacs(string dir)
    {
        var stem = System.IO.Path.Combine(dir, "system");
        _doc!.Gromacs(stem);
        // the core's .mdp ends with the freeze groups of the held atoms (freezegrps, freezedim): kept under this run's deck
        var core = System.IO.File.Exists(stem + ".mdp") ? System.IO.File.ReadAllText(stem + ".mdp") : "";
        var at = core.IndexOf("; atoms held in place", StringComparison.Ordinal);
        System.IO.File.WriteAllText(stem + ".mdp", GromacsDeck(_doc!) + (at >= 0 ? "\n" + core[at..] : ""));
    }
}
