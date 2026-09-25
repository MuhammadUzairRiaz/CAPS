using System.Globalization;

namespace CapsStudio.ViewModels;

/// <summary>Dynamics › Export to other engines › GROMACS: the topology, coordinates and run parameters for the same force
/// field and settings. The .mdp non-bonded section and the notes come from the core (checked against GROMACS, energies
/// and forces: bench/ff/check_gromacs.py); the run section is written here from the Dynamics settings.</summary>
public sealed partial class MainViewModel
{
    private bool _mdGromacs;
    public bool MdGromacs
    {
        get => _mdGromacs;
        set { if (Set(ref _mdGromacs, value)) { Raise(nameof(MdLammps)); Raise(nameof(MdDeckFiles)); Raise(nameof(MdParityTip)); RefreshPreflight(); } }
    }
    public bool MdLammps { get => !_mdGromacs; set => MdGromacs = !value; }
    public string MdDeckFiles => _mdGromacs ? "system.top, system.gro and system.mdp (gmx grompp, then gmx mdrun)" : "system.data and system.in (lmp -in system.in)";
    public string MdParityTip => _mdGromacs
        ? "Energies and forces of the topology match CAPS in GROMACS 2026 (bench/ff/check_gromacs.py): bonded terms, Lennard-Jones and PME Coulomb; the tail correction differs by definition"
        : "Energies and forces of the data file match CAPS in LAMMPS (bench/ff/check_data_lammps.py)";

    /// <summary>The GROMACS .mdp for this run: the core's non-bonded settings, then the integrator and coupling.</summary>
    private string GromacsDeck()
    {
        var inv = CultureInfo.InvariantCulture;
        var nb = _doc!.Gromacs(null);
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
        sb.Append(nb);
        // coupling: τ in ps, pressure in bar
        if (MdHasThermostat)
        {
            sb.Append(_mdThermostat == 1 ? "; tcoupl is implied by the sd integrator\n" : "tcoupl                   = V-rescale   ; Bussi\n");
            sb.Append(string.Format(inv, "tc-grps                  = System\ntau-t                    = {0:0.####}\nref-t                    = {1:0.##}\n", _mdTauT / 1000, _mdTemp));
        }
        else sb.Append("tcoupl                   = no\n");
        if (MdHasBarostat)
        {
            var berendsen = _mdEnsemble == 3 || _mdBarostat == 1;
            sb.Append(string.Format(inv, "pcoupl                   = {0}\npcoupltype               = isotropic\ntau-p                    = {1:0.####}\nref-p                    = {2:0.#####}\ncompressibility          = 4.5e-5\n",
                berendsen ? "Berendsen" : "C-rescale", _mdTauP / 1000, _mdPressure * 1.01325));
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
        System.IO.File.WriteAllText(stem + ".mdp", GromacsDeck());
    }
}
