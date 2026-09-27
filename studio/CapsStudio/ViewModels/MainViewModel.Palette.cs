using System.Collections.ObjectModel;

namespace CapsStudio.ViewModels;

/// <summary>A command the palette can run: a title, a stable id (as scripts would name it), an icon, a shortcut.</summary>
public sealed class PaletteCommand
{
    public required string Title { get; init; }
    public required string Id { get; init; }
    public string Icon { get; init; } = "chevr";
    public string Shortcut { get; init; } = "";
    public string Section { get; init; } = "Commands";
    public string Keywords { get; init; } = "";
    public Func<bool> Enabled { get; init; } = () => true;
    public required Action Run { get; init; }
}

/// <summary>A row of the palette: a section header or a command.</summary>
public sealed record PaletteRow(string Title, string Id, string Icon, string Shortcut, bool IsHeader, PaletteCommand? Command)
{
    public bool HasShortcut => Shortcut.Length > 0;
    public bool IsCommand => !IsHeader;
}

/// <summary>The command palette (design/boards/CommandPalette): ⌘K, type to filter, ↑ ↓ to move, ↵ to run.</summary>
public sealed partial class MainViewModel
{
    private readonly List<PaletteCommand> _commands = new();
    private bool _modelCommands;
    public ObservableCollection<PaletteRow> PaletteRows { get; } = new();
    private bool _paletteOpen;
    public bool PaletteOpen { get => _paletteOpen; set { if (Set(ref _paletteOpen, value) && value) { _paletteQuery = ""; Raise(nameof(PaletteQuery)); FilterPalette(); } } }
    private string _paletteQuery = "";
    public string PaletteQuery { get => _paletteQuery; set { if (Set(ref _paletteQuery, value)) FilterPalette(); } }
    private int _paletteIndex = -1;
    public int PaletteIndex { get => _paletteIndex; set { if (Set(ref _paletteIndex, value)) Raise(nameof(PaletteHint)); } }
    /// <summary>The footer: the id of the command under the cursor.</summary>
    public string PaletteHint => _paletteIndex >= 0 && _paletteIndex < PaletteRows.Count && PaletteRows[_paletteIndex].Command is { } c ? "command: " + c.Id : "";

    public void AddCommand(PaletteCommand c) => _commands.Add(c);

    /// <summary>Runs the command with this id (the menu bar and the palette share them); false when there is none or it
    /// is not available now.</summary>
    public bool RunCommand(string id)
    {
        var c = _commands.FirstOrDefault(x => x.Id == id);
        if (c == null || !c.Enabled()) return false;
        c.Run();
        return true;
    }

    /// <summary>A page by its module number, as the rail, the palette and the menu bar reach it (each page's own set-up
    /// runs through SetModule's openers).</summary>
    public void GoModule(int m)
    {
        switch (m)
        {
            case 9: OpenBuilder(); break;
            case 66: _ = OpenAddHydrogens(); break;
            case 90: OpenExportDialog(); break;
            case 91: OpenExportDialog(1); break;
            case 92: StartTour(); break;
            default: SetModule(m); break;
        }
    }

    /// <summary>Commands that only need the view model (the window adds those that need dialogs).</summary>
    private void AddModelCommands()
    {
        var modules = new (string Name, int M, string Icon, string Words)[]
        {
            ("Studio", 8, "hex", "view structure inspect"), ("Molecule builder", 9, "hex", "smiles sketch draw 3d conformers"),
            ("Grow", 0, "grow", "amorphous cell chains polystyrene"), ("Pack", 5, "pack", "packmol solvate mixture"),
            ("Relax", 2, "relax", "minimise minimize energy"), ("Dynamics", 3, "dyn", "md nvt npt lammps"),
            ("Equilibrate", 4, "equil", "protocol 21-step annealing convergence"), ("React", 6, "react", "crosslink cure gel"),
            ("Analyze", 1, "chart", "properties density rdf tg modulus"), ("Field", 7, "tag", "force field typing gaff opls"),
            ("Jobs", 11, "jobs", "runs progress log provenance history"), ("Settings", 10, "gear", "preferences theme palette threads"),
            ("Bench", 12, "bench", "validation benchmark tables paper"),
            ("Polymer builder", 13, "grow", "repeat unit copolymer smiles library rubber"), ("Surface builder", 14, "layers", "slab cleave cif crystal interface film fibre silica graphite"),
            ("Crystal builder", 29, "cube", "space group lattice unit cell asymmetric unit cif symmetry primitive supercell polyethylene quartz"),
            ("Macro recorder", 36, "terminal", "macro record script python replay automate parameter batch"),
            ("Fragment library", 35, "hex", "fragment library ring functional group monomer amino acid solvent ion additive accelerator sulfur tmtd cbs silane attach"),
            ("Chain statistics", 55, "grow", "chain statistics end to end rg radius of gyration characteristic ratio cn c infinity ree"),
            ("Density calculator", 56, "pack", "density calculator molecules count box volume avogadro solvent mixing pack"),
            ("Surface area", 57, "chart", "surface area sasa solvent accessible shrake rupley exposure probe"),
            ("Unit cell", 58, "cube", "unit cell editor lattice parameters a b c alpha beta gamma supercell volume density triclinic"),
            ("Polydispersity", 59, "grow", "polydispersity chain length distribution schulz zimm flory poisson dispersity mn mw pdi"),
            ("Copolymer builder", 60, "hex", "copolymer terminal model mayo lewis reactivity ratios r1 r2 feed composition azeotrope run length sbr nbr"),
            ("Solvent screen", 61, "flask", "solvent screen hildebrand solubility parameter chi flory huggins good solvent non-solvent"),
            ("Tacticity statistics", 62, "hex", "tacticity triads pentads meso racemo dyads bernoulli markov nmr pm fit"),
            ("Blend phase diagram", 63, "chart", "blend phase diagram flory huggins binodal spinodal critical point ucst lcst chi"),
            ("Electrostatics", 64, "dyn", "electrostatics ewald pme beta cut-off tolerance mesh fft grid spacing"),
            ("Display styles", 65, "eye", "display style all atoms no hydrogens backbone tubes polar h automatic by size"),
            ("Add hydrogens", 66, "plus", "add hydrogens protonate valence heavy atoms pdb missing h"),
            ("Type by hand", 69, "tag", "type by hand manual atom types head body tail repeat unit copolymer junction oplsaa assign override force field example"),
            ("Model resolution", 67, "layers", "model resolution united atom coarse grained beads convert backmap trappe"),
            ("Figure composer", 54, "layers", "figure composer panels layout multi panel journal column svg tiff pdf paper"),
            ("Recipes", 53, "file", "recipe yaml pipeline steps reusable run caps run protocol schedule larsen hash template"),
            ("Orientation", 52, "chart", "orientation order parameter nematic s director herman crystallinity chord draw"),
            ("Periodic box", 51, "cube", "periodic box wrap unwrap whole images boundary minimum image pbc centre"),
            ("Partial charges", 50, "tag", "charges partial gasteiger qeq resp am1-bcc chg colour by charge net charge"),
            ("Glass transition", 47, "chart", "glass transition tg cooling scan specific volume two lines fit replicas soldera"),
            ("Interface analysis", 48, "layers", "interface density profile surface normal gap adsorbed layer adhesion film slab"),
            ("Diffusion", 49, "chart", "diffusion msd einstein coefficient yeh hummer finite size viscosity"),
            ("Colour-vision check", 46, "eye", "colour vision colour blind color blindness protanopia deuteranopia tritanopia accessibility palette contrast"),
            ("Reaction template editor", 45, "react", "reaction template editor atom map pre post bond formed broken epoxy amine sulfur cure"),
            ("Coarse-grained melt", 44, "grow", "coarse grained kremer grest bead spring fene wca lammps melt cg"),
            ("Parameter sweep", 43, "jobs", "sweep parameter grid tacticity chain length dp seed replicas batch runs"),
            ("Project home", 42, "folder", "project home folder methods section paper results share zip documents"),
            ("Theory manual", 41, "file", "theory manual method equation reference thermostat barostat pme dsf lbfgs uff qeq gaff citation documentation help"),
            ("Studio tour", 92, "cursor", "tour help first run tutorial introduction guide walkthrough"),
            ("Provenance", 37, "history", "provenance manifest steps seeds citations bibtex reproduce compare run approximations"),
            ("Export image", 90, "download", "export image png 16-bit 4k poster journal figure dpi clipboard provenance"),
            ("Export movie", 91, "play", "export movie animation apng mp4 turntable frames video"),
            ("Split view", 34, "layers", "split side by side two documents compare sync camera light paper theme"),
            ("Torsion scan", 33, "chart", "torsion dihedral scan rotation barrier conformer gauche trans energy profile rotamer"),
            ("Trajectory player", 32, "play", "trajectory frames play movie timeline log lammps thermo temperature density rg end-to-end ree smooth"),
            ("Solvation builder", 31, "flask", "solvate water tip3p tip4p spc ions salt nacl concentration neutralise toluene swelling solvent box"),
            ("Biomolecule builder", 30, "hex", "peptide protein sequence fasta helix sheet amino acid residue biomolecule"),
            ("Nanostructure builder", 15, "atom", "nanotube cnt graphene sheet nanoparticle filler composite carbon black silica"),
            ("Blend builder", 16, "grow", "blend mixture nr br sbr tyre compound two polymers"),
            ("Adsorption locator", 70, "layers", "adsorption adsorbate substrate surface binding site simulated annealing monte carlo filler gas coupling agent sizing"),
            ("File checks", 17, "check", "problems warnings validation report file errors"),
            ("Figure bundle", 26, "download", "bundle zip reproduce provenance sha256 hash data csv figure paper supplementary"),
            ("Four views", 25, "cube", "viewports top front left side perspective ortho orthographic quad layout"),
            ("Colour by", 24, "eye", "colour color gallery element chain type charge position property map"),
            ("Compare", 23, "layers", "compare two cells seeds runs a b difference side by side"),
            ("Batch", 22, "layers", "batch many files seeds pipeline table results csv compare replicas"),
            ("Export data", 21, "download", "save write lammps data dump gromacs gro pdb xyz mol2 coefficients preview"),
            ("Visualize", 20, "eye", "pipeline steps modifiers ovito colour coding selection expression slice cluster rdf coordination data inspector particles bonds attributes tables"),
            ("Render", 19, "eye", "render image movie frames ambient occlusion overlays label legend scale bar tripod"),
            ("Export figure", 18, "download", "figure png svg image picture paper journal slide poster transparent background dpi scale bar"),
        };
        foreach (var (name, m, icon, words) in modules)
            AddCommand(new PaletteCommand
            {
                Title = $"Go to {name}", Id = $"module.open {name.Split(' ')[0].ToLowerInvariant()}", Icon = icon, Section = "Modules", Keywords = words,
                Run = () => { if (m == 9) OpenBuilder(); else if (m == 14) OpenSurface(); else if (m == 15) OpenNano(); else if (m == 16) OpenBlend(); else if (m == 17) OpenChecks(); else if (m == 18) OpenFigure(); else if (m == 19) OpenRender(); else if (m == 20) OpenVisualize(); else if (m == 21) OpenExport(); else if (m == 22) OpenBatch(); else if (m == 23) OpenCompare(); else if (m == 24) OpenColourBy(); else if (m == 25) OpenViewports(); else if (m == 26) OpenBundle(); else if (m == 29) OpenCrystal(); else if (m == 30) OpenBio(); else if (m == 31) OpenSolvation(); else if (m == 32) OpenTrajectory(); else if (m == 33) OpenTorsion(); else if (m == 34) OpenSplit(); else if (m == 35) OpenFragments(); else if (m == 36) OpenMacro(); else if (m == 37) OpenProvenance(); else if (m == 92) StartTour(); else if (m == 41) OpenManual(); else if (m == 42) OpenProject(); else if (m == 43) OpenSweep(); else if (m == 44) OpenCg(); else if (m == 45) OpenTemplateEditor(); else if (m == 46) OpenColourVision(); else if (m == 47) OpenGlass(); else if (m == 48) OpenInterface(); else if (m == 49) OpenDiffusion(); else if (m == 50) OpenCharges(); else if (m == 51) OpenPeriodic(); else if (m == 52) OpenOrientation(); else if (m == 53) OpenRecipes(); else if (m == 54) OpenComposer(); else if (m == 55) OpenChainStats(); else if (m == 56) OpenDensityCalc(); else if (m == 57) OpenSurfaceArea(); else if (m == 58) OpenCellEditor(); else if (m == 59) OpenPolydispersity(); else if (m == 60) OpenCopolymer(); else if (m == 61) OpenSolventScreen(); else if (m == 62) OpenTacticityStats(); else if (m == 63) OpenBlendPhase(); else if (m == 64) OpenElectrostatics(); else if (m == 65) OpenDisplayStyles(); else if (m == 66) _ = OpenAddHydrogens(); else if (m == 67) OpenModelResolution(); else if (m == 90) OpenExportDialog(); else if (m == 91) OpenExportDialog(1); else SetModule(m); },
            });
        AddCommand(new PaletteCommand { Title = "Start page", Id = "start.open", Icon = "cube", Section = "Modules", Keywords = "home recent new",
            Enabled = () => _doc == null, Run = () => SetModule(8) });
        AddCommand(new PaletteCommand { Title = "Close document", Id = "document.close", Icon = "close", Shortcut = "⌘W", Section = "File",
            Enabled = () => _doc != null, Run = CloseDocument });
        AddCommand(new PaletteCommand { Title = "Reset view", Id = "view.reset", Icon = "rotate", Shortcut = "R", Section = "View",
            Enabled = () => _doc != null, Run = ResetView });
        AddCommand(new PaletteCommand { Title = "Check for updates", Id = "app.update", Icon = "download", Section = "Settings",
            Keywords = "update new version release download changelog upgrade", Run = () => _ = CheckForUpdates() });
        AddCommand(new PaletteCommand { Title = "Select by query", Id = "select.query", Icon = "search", Shortcut = "⌘F", Section = "Selection",
            Keywords = "select query smarts element chain within ring stereo and or not find", Enabled = () => _doc != null, Run = () => { SetModule(8); QueryOpen = true; } });
        AddCommand(new PaletteCommand { Title = "Frame selection", Id = "view.frame", Icon = "cube", Shortcut = "F", Section = "View",
            Keywords = "focus fit zoom to selection centre camera fly", Enabled = () => _doc != null, Run = FrameSelection });
        for (var ax = 0; ax < 3; ++ax)
        {
            var a = ax;
            var name = "xyz"[a].ToString();
            AddCommand(new PaletteCommand { Title = $"Rotate selection 90° about {name}", Id = "edit.rotate " + name, Icon = "rotate", Section = "Edit",
                Keywords = "rotate turn selection molecule axis 90 degrees", Enabled = () => _doc != null, Run = () => RotateSelection(a, 90) });
            AddCommand(new PaletteCommand { Title = $"Mirror selection through the plane normal to {name}", Id = "edit.mirror " + name, Icon = "mirror", Section = "Edit",
                Keywords = "mirror reflect enantiomer image flip selection plane", Enabled = () => _doc != null, Run = () => MirrorSelection(a) });
        }
        foreach (var rs in new[] { "R", "S" })
        {
            var c = rs;
            AddCommand(new PaletteCommand { Title = $"Make the picked stereocentre {c}", Id = "edit.make " + c, Icon = "atom", Section = "Edit",
                Keywords = "stereo chirality configuration cip set r s centre", Enabled = () => _doc != null && _selection.Count == 1, Run = () => MakePicked(c) });
        }
        AddCommand(new PaletteCommand { Title = "Full-screen view", Id = "view.fullscreen", Icon = "expand", Section = "View",
            Keywords = "full screen presentation screenshot clean view only visualization window big", Enabled = () => _doc != null || RunLiveDoc != null,
            Run = () => FullViewRequested?.Invoke(RunLiveDoc != null) });
        AddCommand(new PaletteCommand { Title = "Auto-clean while editing", Id = "edit.autoclean", Icon = "wand", Shortcut = "A", Section = "Edit",
            Keywords = "auto clean uff relax geometry while building toggle", Run = ToggleAutoClean });
        AddCommand(new PaletteCommand { Title = "Reduce motion", Id = "settings.motion", Icon = "gear", Section = "Settings",
            Keywords = "motion animation reduce accessibility camera fly cut vestibular", Run = () => SetReduceMotion = Motion.Reduced ? "off" : "on" });
        for (var k = 0; k < Styles.Length; k++)
        {
            var s = k;
            AddCommand(new PaletteCommand { Title = "Display: " + Styles[k], Id = "view.style " + Styles[k].ToLowerInvariant().Replace(' ', '_').Replace("&", "and"), Icon = "eye",
                Section = "View", Keywords = "style representation", Enabled = () => _doc != null, Run = () => StyleIndex = s });
        }
        for (var k = 0; k < ColourModes.Length; k++)
        {
            var c = k;
            AddCommand(new PaletteCommand { Title = "Colour by " + ColourModes[k].ToLowerInvariant(), Id = "view.colour " + ColourModes[k].Split(' ')[0].ToLowerInvariant(),
                Icon = "layers", Section = "View", Keywords = "color", Enabled = () => _doc != null, Run = () => ColourIndex = c });
        }
        AddCommand(new PaletteCommand { Title = "Toggle perspective / orthographic", Id = "view.projection", Icon = "cube", Section = "View",
            Enabled = () => _doc != null, Run = () => Perspective = !Perspective });
        AddCommand(new PaletteCommand { Title = "Theme: Graphite (dark)", Id = "theme dark", Icon = "gear", Section = "View", Keywords = "appearance", Run = () => Tokens.Use(false) });
        AddCommand(new PaletteCommand { Title = "Theme: Paper (light)", Id = "theme light", Icon = "gear", Section = "View", Keywords = "appearance", Run = () => Tokens.Use(true) });

        AddCommand(new PaletteCommand { Title = "Grow an amorphous cell (current settings)", Id = "grow.run", Icon = "grow", Section = "Run",
            Enabled = () => Idle, Run = () => { SetModule(0); _ = Grow(); } });
        AddCommand(new PaletteCommand { Title = "Relax the structure", Id = "relax.run", Icon = "relax", Section = "Run", Keywords = "minimise",
            Enabled = () => CanRelax, Run = () => { SetModule(2); _ = Relax(); } });
        AddCommand(new PaletteCommand { Title = "Run dynamics", Id = "md.run", Icon = "dyn", Section = "Run", Keywords = "md",
            Enabled = () => CanRun, Run = () => { SetModule(3); _ = RunMd(); } });
        AddCommand(new PaletteCommand { Title = "Equilibrate with the protocol", Id = "equilibrate.run", Icon = "equil", Section = "Run",
            Enabled = () => CanEquilibrate, Run = () => { SetModule(4); _ = RunEquilibrate(); } });
        AddCommand(new PaletteCommand { Title = "Crosslink (React)", Id = "react.run", Icon = "react", Section = "Run",
            Enabled = () => CanReact, Run = () => { SetModule(6); _ = RunReact(); } });
        AddCommand(new PaletteCommand { Title = "Run the validation suite", Id = "bench.run --all", Icon = "bench", Section = "Run", Keywords = "benchmark tables",
            Enabled = () => BenchIdle, Run = () => { SetModule(12); LoadBench(); _ = RunBench(true); } });
        AddCommand(new PaletteCommand { Title = "Compute the selected properties", Id = "analyze.run", Icon = "chart", Section = "Run", Keywords = "analyze analyse",
            Enabled = () => _doc != null && Idle, Run = () => { SetModule(1); _ = Analyze.Run(); } });
    }

    private void FilterPalette()
    {
        if (!_modelCommands) { _modelCommands = true; AddModelCommands(); }
        PaletteRows.Clear();
        var q = _paletteQuery.Trim();
        var sections = new List<(string Section, List<(PaletteCommand C, int Score)> Items)>();
        void Put(string section, PaletteCommand c, int score)
        {
            var s = sections.FirstOrDefault(x => x.Section == section);
            if (s.Items == null) { s = (section, new List<(PaletteCommand, int)>()); sections.Add(s); }
            s.Items.Add((c, score));
        }
        // what the text itself is: a SMILES to build, a file to open
        if (q.Length > 0)
        {
            var t = q.Trim('"', '\'');
            if (LooksLikePath(t) && File.Exists(ExpandHome(t)))
                Put("From what you typed", new PaletteCommand { Title = "Open " + Path.GetFileName(t), Id = "document.open", Icon = "folder", Run = () => OpenRequested?.Invoke(ExpandHome(t)) }, 1000);
            else if (IsSmiles(t) && FindModule(t) == null)
                Put("From what you typed", new PaletteCommand { Title = "Build 3D from SMILES " + t, Id = "builder.molecule.open", Icon = "hex", Run = () => OpenBuilder(t) }, 1000);
        }
        foreach (var c in _commands)
        {
            if (!c.Enabled()) continue;
            var score = q.Length == 0 ? 1 : Score(q, c.Title + " " + c.Keywords, c.Id);
            if (score > 0) Put(c.Section, c, score);
        }
        // documents: the recent list
        foreach (var r in Recent)
        {
            var score = q.Length == 0 ? 1 : Score(q, r.Name + " " + r.Folder, "document.open");
            var item = r;
            if (score > 0)
                Put("Documents", new PaletteCommand { Title = r.Name, Id = "document.open " + r.Folder, Icon = "file",
                    Run = () => OpenRequested?.Invoke(item.Path + (item.Topology != null ? "\n" + item.Topology : "")) }, score);
        }
        // best matches first: sections ordered by their best score, commands by score
        var order = q.Length == 0 ? sections : sections.OrderByDescending(s => s.Items.Max(x => x.Score)).ToList();
        var shown = 0;
        foreach (var (section, items) in order)
        {
            if (shown >= 40) break;
            PaletteRows.Add(new PaletteRow(section.ToUpperInvariant(), "", "", "", true, null));
            foreach (var (c, _) in (q.Length == 0 ? items.AsEnumerable() : items.OrderByDescending(x => x.Score)).Take(q.Length == 0 ? 8 : 12))
            {
                PaletteRows.Add(new PaletteRow(c.Title, c.Id, c.Icon, ShortcutShown(c), false, c));
                shown++;
            }
        }
        _paletteIndex = -1;
        PaletteIndex = PaletteRows.Select((r, i) => (r, i)).FirstOrDefault(x => x.r.IsCommand, (null!, -1)).Item2;
    }

    /// <summary>Match: every query word must appear in the title, keywords or id (a word start scores higher); a
    /// word not found as written may match the title's letters in order within a short span ("eqlb" → Equilibrate).</summary>
    private static int Score(string query, string text, string id)
    {
        var hay = (text + " " + id).ToLowerInvariant();
        var title = text.ToLowerInvariant();
        var total = 0;
        foreach (var w in query.ToLowerInvariant().Split(' ', StringSplitOptions.RemoveEmptyEntries))
        {
            var at = hay.IndexOf(w, StringComparison.Ordinal);
            if (at >= 0)
            {
                var start = at == 0 || !char.IsLetterOrDigit(hay[at - 1]);
                total += (start ? 30 : 12) + w.Length + (at < title.Length ? 10 : 0);
                continue;
            }
            if (w.Length < 3) return 0;
            var best = int.MaxValue;
            for (var s0 = 0; s0 < title.Length; s0++)
            {
                if (title[s0] != w[0] || (s0 > 0 && char.IsLetterOrDigit(title[s0 - 1]))) continue;
                var k = 1;
                var e = s0 + 1;
                for (; e < title.Length && k < w.Length; e++) if (title[e] == w[k]) k++;
                if (k == w.Length) best = Math.Min(best, e - s0);
            }
            if (best > 2 * w.Length + 2) return 0;
            total += 6;
        }
        return total;
    }

    public void PaletteMove(int delta)
    {
        if (PaletteRows.Count == 0) return;
        var i = _paletteIndex;
        for (var n = 0; n < PaletteRows.Count; n++)
        {
            i = (i + delta + PaletteRows.Count) % PaletteRows.Count;
            if (PaletteRows[i].IsCommand) { PaletteIndex = i; return; }
        }
    }

    public void PaletteRun(PaletteRow? row = null)
    {
        row ??= _paletteIndex >= 0 && _paletteIndex < PaletteRows.Count ? PaletteRows[_paletteIndex] : null;
        if (row?.Command is not { } c) return;
        PaletteOpen = false;
        try { c.Run(); }
        catch (Exception ex) { Status = $"{c.Title}: {ex.Message}"; }
    }

    /// <summary>The window opens files (a path, or "path\ntopology").</summary>
    public event Action<string>? OpenRequested;
}
