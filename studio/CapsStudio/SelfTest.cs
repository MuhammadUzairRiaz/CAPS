using CapsStudio.Interop;
using CapsStudio.ViewModels;

namespace CapsStudio;

/// <summary>Headless check of the Studio's path into the native core: open, summary, pick, render, export.</summary>
internal static class SelfTest
{
    public static int Run(string[] args)
    {
        var fails = 0;
        void Check(bool ok, string what) { Console.WriteLine($"{(ok ? "ok  " : "FAIL")} {what}"); if (!ok) fails++; }

        Check(Native.AbiVersion() == 20, "native ABI version 20");
        var dir = args.Length > 0 ? args[0] : "samples";
        var outDir = args.Length > 1 ? args[1] : Path.GetTempPath();
        AppSettings.Override = Path.Combine(outDir, "caps-selftest-settings.json");
        if (File.Exists(AppSettings.Override)) File.Delete(AppSettings.Override);
        RecentFiles.Override = Path.Combine(outDir, "caps-selftest-recent");
        if (Directory.Exists(RecentFiles.Override)) Directory.Delete(RecentFiles.Override, true);
        if (File.Exists(MainViewModel.JobsFile)) File.Delete(MainViewModel.JobsFile);
        var vm = new MainViewModel();
        vm.HookJobs();

        // Start: the quick-start box recognises SMILES, files and modules
        Check(MainViewModel.IsSmiles("C=Cc1ccccc1") && MainViewModel.IsSmiles("CC(=O)O[C@@H]1CCCC1") && MainViewModel.IsSmiles("[Na+].[Cl-]")
              && !MainViewModel.IsSmiles("Pack") && !MainViewModel.IsSmiles("C(C") && !MainViewModel.IsSmiles("hello world"), "SMILES grammar");
        vm.QuickText = "C=Cc1ccccc1";
        var k1 = vm.QuickKind;
        vm.QuickText = Path.Combine(dir, "ps_melt.data");
        var k2 = vm.QuickKind;
        vm.QuickText = "equil";
        var k3 = vm.QuickKind;
        vm.QuickGo();
        Check(k1 == 2 && k2 == 1 && k3 == 3 && vm.IsEquilibrate, $"quick start: SMILES {k1}, file {k2}, module {k3} → Equilibrate {vm.IsEquilibrate}");
        vm.SetModule(8);
        vm.QuickText = "";

        vm.Open(Path.Combine(dir, "ps_melt.lammpstrj"), Path.Combine(dir, "ps_melt.data"));
        for (var i = 0; i < 100 && RecentFiles.Load().Count == 0; i++) Thread.Sleep(50);
        var recent = RecentFiles.Load();
        Check(recent.Count == 1 && recent[0].Name == "ps_melt.lammpstrj" && recent[0].Topology != null, $"recent: {string.Join(", ", recent.Select(r => r.Name))}");
        var s = vm.Document!.Summary();
        Check(s.Atoms == 1300 && s.Bonds == 1370 && s.Molecules == 10, $"summary: {s.Atoms} atoms, {s.Bonds} bonds, {s.Molecules} molecules");
        Check(s.Frames == 3 && vm.HasFrames, $"frames: {s.Frames}");
        Check(Math.Abs(s.Density - 0.386) < 0.001, $"density {s.Density:F4} g/cm³");
        Check(s.Format == "lammps-dump", $"format '{s.Format}'");

        // Command palette: filters commands, offers a typed SMILES, runs the chosen one
        vm.PaletteOpen = true;
        vm.PaletteQuery = "colour mol";
        var ids = vm.PaletteRows.Where(r => r.IsCommand).Select(r => r.Id).ToList();
        vm.PaletteRun();
        Check(ids.FirstOrDefault() == "view.colour molecule" && vm.ColourIndex == 1 && !vm.PaletteOpen, "palette: " + string.Join(", ", ids.Take(3)));
        vm.ColourIndex = 0;
        vm.PaletteOpen = true;
        vm.PaletteQuery = "c1ccccc1O";
        Check(vm.PaletteRows.Any(r => r.Id == "builder.molecule.open"), "palette: a typed SMILES builds in 3D");
        vm.PaletteOpen = false;

        var a = vm.Document.Atom(40);
        Check(a.Id == 41 && a.ElementSymbol == "H", $"atom 41: {a.ElementSymbol}, mol {a.Mol}, q {a.Charge:+0.0000}");
        var nb = vm.Document.Neighbours(40, 3);
        Check(nb.Length == 3 && nb[0].Distance > 1.0 && nb[0].Distance < 1.1, $"nearest neighbour {nb[0].Distance:F3} Å (C–H bond)");

        var opt = vm.ViewOptions(640, 400, 1);
        var buf = new byte[640 * 400 * 4];
        vm.Document.Render(vm.Camera, opt, buf);
        var covered = 0;
        for (var k = 3; k < buf.Length; k += 4) covered += buf[k] == 255 ? 1 : 0;
        Check(covered == 640 * 400, "view render is opaque (dark background)");
        var hits = 0;
        for (var y = 150; y < 250; y += 5) for (var x = 250; x < 390; x += 5) hits += vm.Document.Pick(x, y) >= 0 ? 1 : 0;
        Check(hits > 0, $"picking finds atoms near the centre ({hits} hits)");
        vm.Pick(vm.Document.Pick(320, 200) is var p && p >= 0 ? p : 40);
        Check(vm.PickedRows.Count == 4 && vm.NeighbourRows.Count == 4, $"inspector: {vm.PickedTitle}");

        // CAPS Field: assign GAFF2 from the library, override one atom, clear
        var gaff = vm.Field.Library.ToList().FindIndex(x => x.Id == "gaff-amber25-dlfield");
        Check(gaff >= 0, $"force-field library: {vm.Field.LibraryNote}");
        if (gaff >= 0)
        {
            // the core call on its own first, so a failure names its cause
            try
            {
                var file = vm.Field.Library[gaff].File;
                var done = vm.Document!.FieldAssign(file, null, 1);
                Check(vm.Document.FieldReport().Length > 1000, $"field assign in the core: complete {done}, report {vm.Document.FieldReport().Length} chars, {file}");
                vm.Document.FieldClear();
            }
            catch (Exception e) { Check(false, "field assign in the core: " + e.Message); }
            vm.Field.FfIndex = gaff;
            vm.Field.ChargeMode = 1;
            vm.Field.Assign().GetAwaiter().GetResult();
            Check(vm.Field.Assigned && vm.Field.Complete, $"Field: {vm.Field.ForceFieldName} · {vm.Field.TypedText} · {vm.Field.MissingText} {vm.Field.Log}");
            Check(vm.Field.Swatches.Select(x => x.Name).OrderBy(x => x).SequenceEqual(["c3", "ca", "ha", "hc"]), "polystyrene types c3 ca ha hc: " + string.Join(" ", vm.Field.Swatches.Select(x => x.Label)));
            vm.Field.SelectAtom(0);
            Check(vm.Field.WhyText.StartsWith("Rule "), "why: " + vm.Field.WhyText);
            var typedAs = vm.Document!.Atom(0).Type;
            vm.Field.OverrideType = "c3";
            vm.Field.ApplyOverride().GetAwaiter().GetResult();
            Check(vm.Field.SelectedRow?.Overridden == true, $"override: atom 1 set to {vm.Field.SelectedRow?.Type}");
            vm.Field.Clear().GetAwaiter().GetResult();
            Check(!vm.Field.Assigned && vm.ForceFieldLine.StartsWith("Force field: built-in"), "clear: " + vm.ForceFieldLine);
            Check(vm.Document.Atom(0).Name.Length > 0, $"clear restores the file's types (atom 1 {vm.Document.Atom(0).Name}, type {typedAs} while assigned)");
        }
        // UFF from the same library: every atom typed from its bonds (polystyrene: C_3, C_R, H_)
        var uffFf = vm.Field.Library.ToList().FindIndex(x => x.Id == "uff");
        if (uffFf >= 0)
        {
            vm.Field.FfIndex = uffFf;
            vm.Field.ChargeMode = 0;
            vm.Field.Assign().GetAwaiter().GetResult();
            Check(vm.Field.Assigned && vm.Field.Complete && vm.Field.Swatches.Select(x => x.Name).OrderBy(x => x).SequenceEqual(["C_3", "C_R", "H_"]),
                  $"Field UFF: {vm.Field.TypedText} · {string.Join(" ", vm.Field.Swatches.Select(x => x.Label))} {vm.Field.Log}");
            vm.Field.Clear().GetAwaiter().GetResult();
        }
        else Check(false, "force-field library has no UFF entry");

        // Analyze › Properties: density, g(r), Rg, MSD and free volume over the three frames; export
        {
            foreach (var c in vm.Analyze.Groups.SelectMany(g => g.Chips)) c.IsOn = c.Id is "density" or "rdf" or "rg" or "msd" or "ffv";
            vm.Analyze.FramePsD = 1;
            vm.Analyze.RefIndex = vm.Analyze.References.ToList().FindIndex(r => r.Id == "ps-atactic");
            vm.Analyze.Run().GetAwaiter().GetResult();
            var cards = vm.Analyze.Results.ToDictionary(c => c.Id);
            Check(cards.Count == 5, $"analyze: {cards.Count} cards · {vm.Analyze.Log}");
            Check(cards.TryGetValue("density", out var d) && Math.Abs(d.Value - s.Density) < 1e-3, $"analyze density {d?.ValueText} g/cm³ (file {s.Density:F4})");
            Check(d?.HasRef == true && d.Verdict.Contains("below"), $"compared with experiment: {d?.RefText} · {d?.Verdict}");
            Check(cards.TryGetValue("rg", out var rg) && rg.Value > 2 && rg.Value < 20, $"analyze Rg {rg?.ValueText} Å");
            Check(vm.Analyze.Curves.Count >= 4 && vm.Analyze.Curve != null, $"analyze curves: {vm.Analyze.Curves.Count}");
            var adir = Path.Combine(outDir, "analyze_selftest");
            var files = vm.Analyze.Export(adir);
            Check(files >= 6 && File.ReadAllText(Path.Combine(adir, "results.tex")).Contains("\\toprule"), $"analyze export: {files} files");
            Check(!vm.Busy && vm.Idle, "analyze leaves the document idle");

            // protocols on a copy of the current frame: a short tensile run and a five-temperature cooling run
            foreach (var c in vm.Analyze.Groups.SelectMany(g => g.Chips)) c.IsOn = c.Id is "tensile" or "tg";
            vm.Analyze.TensRateD = 0.01m;
            vm.Analyze.TensMaxD = 0.01m;
            vm.Analyze.TensTD = 300;
            vm.Analyze.TgFromD = 500;
            vm.Analyze.TgToD = 300;
            vm.Analyze.TgStepD = 50;
            vm.Analyze.TgPsD = 1;
            vm.Analyze.EqPsD = 1;
            var atoms0 = vm.Document!.Summary().Atoms;
            vm.Analyze.Run().GetAwaiter().GetResult();
            var pc = vm.Analyze.Results.ToDictionary(c => c.Id);
            Check(pc.ContainsKey("tensile_modulus") && pc.ContainsKey("yield") && pc.ContainsKey("tg"), $"protocols: {string.Join(", ", pc.Keys)} · {vm.Analyze.Log}");
            Check(pc.TryGetValue("tg", out var tg) && tg.HasValue, $"tg card: {tg?.ValueText} K ({tg?.NotesText.Split('\n')[0]})");
            Check(vm.Analyze.Curves.Any(c => c.Markers && c.OverlayX != null), "tg and stress–strain curves carry points with a line through them");
            Check(vm.Document.Summary().Atoms == atoms0 && vm.Document.Summary().Frames == s.Frames, "protocols leave the document unchanged");
        }

        vm.Frame = 2;
        var f2 = vm.Document.Atom(10);
        vm.Frame = 0;
        var f0 = vm.Document.Atom(10);
        Check(Math.Abs(f2.X - f0.X - 1.0) < 1e-3, $"frame 2 moves atom 11 by {f2.X - f0.X:F3} Å in x (sample pattern)");

        foreach (var (bg, name) in new[] { (0, "dark"), (1, "white"), (2, "transparent") })
        {
            vm.ExportBackground = bg;
            var o = vm.ExportOptions(960, 540);
            var png = Path.Combine(outDir, $"studio_selftest_{name}.png");
            vm.Document.ExportPng(vm.Camera, o, png);
            var bytes = File.ReadAllBytes(png);
            Check(bytes.Length > 1000 && bytes[1] == (byte)'P', $"export {name} PNG ({bytes.Length / 1024} KB)");
            var svg = Path.ChangeExtension(png, ".svg");
            vm.Document.ExportSvg(vm.Camera, o, svg);
            var text = File.ReadAllText(svg);
            Check(text.Contains("<rect width=\"100%\"") == (bg != 2), $"export {name} SVG background shape {(bg != 2 ? "present" : "absent")}");
        }
        vm.Document.Dispose();

        // Grow: a small cell through the same interop the Grow panel uses, saved and reopened.
        var (grown, report) = CapsDocument.Grow(new CapsGrowOpts { Chains = 4, Dp = 6, Density = 0.35, ContactScale = 1.0, Seed = 5, Curve = 1 }, null, "grow_selftest");
        var gs = grown.Summary();
        Check(gs.Atoms == 4 * (16 * 6 + 2) && gs.Molecules == 4, $"grow: {gs.Atoms} atoms in {gs.Molecules} chains");
        Check(Math.Abs(gs.Density - 0.35) < 1e-6, $"grow: density {gs.Density:F4} g/cm³");
        Check(report.Contains("worst contact margin"), "grow: report returned");
        var saved = Path.Combine(outDir, "grow_selftest.data");
        grown.Save(saved);
        grown.Dispose();
        using (var reopened = CapsDocument.Open(saved))
        {
            var rs = reopened.Summary();
            Check(rs.Atoms == gs.Atoms && rs.Bonds == gs.Bonds && rs.BondsFromFile == 1, $"grow: saved LAMMPS data reopens with {rs.Bonds} bonds");
        }
        var cancelled = false;
        try { CapsDocument.Grow(new CapsGrowOpts { Chains = 50, Dp = 30, Density = 0.4, ContactScale = 1.0, Seed = 1, Curve = 1 }, (d, t, r) => false, "x"); }
        catch (InvalidOperationException e) { cancelled = e.Message == "cancelled"; }
        Check(cancelled, "grow: cancel stops the build");

        // Relax: type, push off, compress and minimise through the Relax panel's interop; save with the force field.
        var (cell, _) = CapsDocument.Grow(new CapsGrowOpts { Chains = 4, Dp = 6, Density = 0.4, ContactScale = 0.85, Seed = 2, Curve = 1 }, null, "relax_selftest");
        var info = cell.FieldInfo();
        Check(info.Contains("GAFF") && info.Contains("c3") && info.Contains("energy"), "relax: force-field summary");
        var calls = 0;
        var (conv, rrep) = cell.Relax(new CapsRelaxOpts { Method = 2, Ftol = 1.0, MaxIterations = 5000, TargetDensity = 0.9, CompressStep = 0.1, Pushoff = 1, Cutoff = 10, Coulomb = 1 },
            (st, n, it, e, f, d) => { calls++; return true; });
        var cs = cell.Summary();
        Check(conv && Math.Abs(cs.Density - 0.9) < 1e-6, $"relax: compressed to {cs.Density:F4} g/cm³, converged {conv}");
        Check(cs.Frames > 3 && calls > 0, $"relax: {cs.Frames} frames recorded, {calls} progress calls");
        Check(rrep.Contains("L-BFGS") && rrep.Contains("push-off"), "relax: report lists the stages");
        var ffdata = Path.Combine(outDir, "relax_selftest.data");
        cell.Save(ffdata);
        var ftext = File.ReadAllText(ffdata);
        Check(ftext.Contains("Angle Coeffs") && ftext.Contains("Dihedrals") && ftext.Contains("lj/cut/coul/dsf"), "relax: LAMMPS data carries the force field");
        using (var back = CapsDocument.Open(ffdata))
        {
            var bs = back.Summary();
            Check(bs.Atoms == cs.Atoms && bs.Bonds == cs.Bonds && Math.Abs(bs.Density - 0.9) < 1e-3, $"relax: saved file reopens ({bs.Atoms} atoms, {bs.Density:F4} g/cm³)");
        }
        // Dynamics on the relaxed cell: NVT, then a continuation that reuses the velocities, then a cancelled run.
        var rows = new List<CapsThermo>();
        var mdOpts = new CapsMdOpts { Dt = 1, Steps = 400, Temperature = 300, Thermostat = 1, TauT = 100, Seed = 3, ThermoEvery = 20, FrameEvery = 100, Cutoff = 10, Coulomb = 1, Tail = 1 };
        var mdrep = cell.Md(mdOpts, (r, n) => { rows.Add(r); return true; });
        var ms = cell.Summary();
        Check(ms.Frames == 5 && rows.Count == 21 && rows[^1].Step == 400, $"md: {ms.Frames} frames, {rows.Count} thermo rows");
        Check(rows.Skip(5).Average(r => r.Temperature) is var tm && tm > 200 && tm < 400, $"md: mean temperature {rows.Skip(5).Average(r => r.Temperature):F0} K");
        var traj = Path.Combine(outDir, "md_selftest.lammpstrj");
        var mddata = Path.Combine(outDir, "md_selftest.data");
        cell.SaveTrajectory(traj);
        cell.Save(mddata);
        Check(File.ReadAllText(mddata).Contains("Velocities"), "md: saved data carries velocities");
        using (var tr = CapsDocument.Open(traj, mddata))
            Check(tr.Summary().Frames == 5, $"md: trajectory reopens with {tr.Summary().Frames} frames");
        var cont = cell.Md(mdOpts with { Steps = 100 }, null);
        Check(cont.Contains("velocities taken"), "md: a second run continues with the same velocities");
        // Equilibrate: a named protocol's text, a short custom protocol with convergence blocks, chain statistics.
        var l21 = CapsDocument.ProtocolText("larsen21", new CapsProtocolParams { TFinal = 300, TMax = 600, PFinal = 1, PMax = 49346.2, TimeScale = 1 });
        Check(l21.Split('\n', StringSplitOptions.RemoveEmptyEntries).Length == 21 && l21.Contains("P 49346.2 atm"), "equilibrate: Larsen 21-step text");
        var eqRows = 0;
        var eqOpts = new CapsEquilOpts { Dt = 1, Thermostat = 1, Barostat = 1, TauT = 100, TauP = 1000, Seed = 2, Cutoff = 10, Coulomb = 1, Tail = 1,
                                         FramePs = 0.5, ThermoPs = 0.1, UntilConverged = 1, BlockPs = 0.3, MaxBlocks = 3, TolDensity = 0.5, TolEnergy = 5, TolRg = 0.5 };
        var (eqConv, eqRep) = cell.Equilibrate("nvt 0.5 ps T 450 # warm\nnpt 0.5 ps T 300 P 1 atm # settle", eqOpts, (st, n, label, r) => { eqRows++; return true; });
        Check(eqConv && eqRows > 5 && eqRep.Contains("warm") && eqRep.Contains("check density"), $"equilibrate: 2 stages + blocks, {eqRows} rows, converged {eqConv}");
        var chainsInfo = cell.InternalDistances();
        Check(chainsInfo.Chains == 4 && chainsInfo.N.Length > 5 && Math.Abs(chainsInfo.Ratio[0] - 1) < 1e-9, $"chains: {chainsInfo.Chains} backbones, {chainsInfo.N.Length} separations");
        var badText = false;
        try { cell.Equilibrate("npt 5 ps T 300", eqOpts, null); }
        catch (InvalidOperationException e) { badText = e.Message.Contains("pressure"); }
        Check(badText, "equilibrate: a protocol line without a pressure is reported");
        // React: C–C crosslinks through the React panel's interop (two cycles, relaxed, no dynamics).
        var atomsBefore = cell.Summary().Atoms;
        var rxRows = 0;
        var rxRep = cell.React(CapsDocument.ReactionTemplate("cc_crosslink"),
            new CapsReactOpts { Seed = 1, MaxCycles = 2, MaxPerCycle = 2, TargetConversion = 1, Relax = 1, RelaxIterations = 300, Cutoff = 10, Coulomb = 1 },
            r => { rxRows++; return true; });
        var rxAfter = cell.Summary();
        Check(rxRows == 2 && rxAfter.Atoms < atomsBefore && (atomsBefore - rxAfter.Atoms) % 2 == 0 && rxRep.Contains("reactions"),
              $"react: {(atomsBefore - rxAfter.Atoms) / 2} crosslinks, {rxAfter.Molecules} molecules left");
        var names = CapsDocument.ReactionTemplate("");
        Check(names.Contains("epoxy_amine_primary") && names.Contains("cc_crosslink"), "react: built-in templates listed");
        var before = cell.Summary().Frames;
        var mdCancelled = false;
        try { cell.Md(mdOpts with { Steps = 100000 }, (r, n) => r.Step < 100); }
        catch (InvalidOperationException e) { mdCancelled = e.Message.Contains("cancelled"); }
        Check(mdCancelled && cell.Summary().Frames == before, "md: cancel leaves the document unchanged");

        var relaxCancelled = false;
        try { cell.Relax(new CapsRelaxOpts { Method = 3, Ftol = 0.001, MaxIterations = 100000, Pushoff = 1, Cutoff = 10, Coulomb = 1 }, (st, n, it, e, f, d) => false); }
        catch (InvalidOperationException e) { relaxCancelled = e.Message.Contains("cancelled"); }
        Check(relaxCancelled && cell.Summary().Frames == before, "relax: cancel leaves the document unchanged");
        cell.Dispose();

        // Pack: packmol-format text through the Pack panel's interop; an impossible request must fail and return nothing.
        var packText = "tolerance 2.0\nseed 3\npbc 0 0 0 22 22 22\nstructure water.pdb\n  number 300\n  inside box 0 0 0 22 22 22\nend structure\n";
        var (packed, packRep) = CapsDocument.Pack(packText, dir, null, "pack_selftest");
        var ps = packed.Summary();
        Check(ps.Atoms == 900 && ps.Molecules == 300 && packRep.Contains("0 pairs closer"), $"pack: {ps.Molecules} waters, report: {packRep.Split('\n')[^1]}");
        packed.Dispose();
        var packFailed = false;
        try { CapsDocument.Pack(packText.Replace("number 300", "number 900"), dir, null, "x"); }
        catch (InvalidOperationException e) { packFailed = e.Message.Contains("could not pack"); }
        Check(packFailed, "pack: an overfull cell fails with a reason and no structure");

        // Molecule builder: SMILES → 3D with GAFF2, then into the Studio
        vm.OpenBuilder("N[C@@H](C)C(=O)O");
        Check(vm.MolOk && vm.MolFormula == "C₃H₇NO₂", $"molecule: {vm.MolFormula} · {vm.MolMass} · {vm.MolStereo} stereocentre");
        vm.BuildMolecule().GetAwaiter().GetResult();
        Check(vm.MolDoc != null && vm.MolDoc.Summary().Atoms == 13 && vm.MolConformers.Count >= 1 && vm.MolConformers[0].Minimised,
              $"molecule: {vm.MolConformers.Count} conformers · {vm.MolMethodUsed} {vm.MolNotes}");
        var sk = CapsDocument.SmilesDepict("c1ccccc1O");
        Check(sk.Contains("\"atoms\"") && CapsDocument.SmilesWrite(sk).Length > 0, "molecule: sketch graph round trip " + CapsDocument.SmilesWrite(sk));
        vm.MolSmiles = "C1CC";
        Check(!vm.MolOk && vm.MolHasError, "molecule: bad SMILES reported: " + vm.MolError);
        // UFF: a siloxane (silicone rubber) cleaned up with every element typed
        var uffIx = vm.CleanChoices.ToList().FindIndex(c => c.File == "uff");
        vm.MolSmiles = "C[Si](C)(C)O[Si](C)(C)O[Si](C)(C)C";
        vm.MolClean = uffIx;
        vm.BuildMolecule().GetAwaiter().GetResult();
        Check(uffIx >= 0 && vm.MolConformers.Count >= 1 && vm.MolConformers[0].Minimised && vm.MolMethodUsed.Contains("UFF"),
              $"molecule: UFF clean-up of a siloxane · {vm.MolMethodUsed} {vm.MolNotes}");
        vm.MolClean = 0;
        vm.MolSmiles = "N[C@@H](C)C(=O)O";
        vm.BuildMolecule().GetAwaiter().GetResult();
        vm.OpenMoleculeInStudio();
        Check(vm.IsStudio && vm.Document?.Summary().Atoms == 13, "molecule: opened in the Studio as the document");

        // Settings: saved to the file, applied to the renderer
        vm.SetPalette = 1;
        vm.SetThreads = 2;
        var reread = AppSettings.Load();
        var okabe = new byte[64 * 64 * 4];
        Check(reread.Palette == 1 && reread.Threads == 2 && File.Exists(AppSettings.FilePath), $"settings: saved to {AppSettings.FilePath}");
        vm.ResetSettings();
        Check(AppSettings.Load().Palette == 0 && vm.SetThreads == 0, "settings: reset to the defaults");

        // Polymer builder: the library, a rubber preset, a grown cell of it
        vm.LoadPolymerLibrary();
        var enr = vm.PolymerLibrary.FirstOrDefault(e => e.Id == "C103");
        Check(vm.PolymerLibrary.Count(e => !e.Copolymer) > 100 && enr != null, $"polymer library: {vm.LibraryCount}");
        if (enr != null)
        {
            vm.UseLibrary(enr, null);
            Check(!vm.PolyHasError && vm.PolyUnits.Count == 2 && vm.PolyStripUnits.Length == vm.GrowDpD, $"ENR-50 chain: {vm.PolyPreview} {vm.PolyError}");
            vm.GrowChainsD = 4;
            vm.GrowDensityD = 0.3m;
            vm.SendPolymerToGrow();
            vm.Grow().GetAwaiter().GetResult();
            Check(vm.Document != null && vm.Document.Summary().Molecules == 4 && vm.GrowComponentName.StartsWith("ENR"), $"grown ENR-50 cell: {vm.Status}");
            // sulfur cure of the rubber: H–S–S–H donors inserted, allylic C–S bonds formed
            vm.RxSet = 2;
            vm.RxInsertCount = 6;
            vm.InsertCurative().GetAwaiter().GetResult();
            Check(vm.Document!.Summary().Molecules == 4 + 6, $"sulfur donors inserted: {vm.Document.Summary().Molecules} molecules · {vm.RxLog.Split('\n')[0]}");
            vm.RxRelax = false;
            vm.RxCyclesD = 6;
            vm.RunReact().GetAwaiter().GetResult();
            var cured = System.Text.RegularExpressions.Regex.Match(vm.RxLog, @"^(\d+) reactions");
            Check(cured.Success && int.Parse(cured.Groups[1].Value) > 0, $"sulfur cure: {vm.RxLog.Split('\n')[0]}");
            vm.RxRelax = true;
            vm.RxSet = 0;
            vm.UsePolystyreneInGrow();
        }

        // Surface builder: quartz (001) terminations, a hydroxylated slab, and a thin rubber film grown on it
        vm.OpenSurface();
        Check(vm.Crystals.Count >= 10 && vm.SurfTerminations.Count == 3 && vm.SurfTerminations[0].StartsWith("O-terminated"),
              $"surface: {vm.Crystals.Count} crystals · {string.Join(" | ", vm.SurfTerminations)} · {vm.SurfError}");
        vm.SurfLayers = 1;
        vm.FilmThickness = 12;
        vm.FilmDensity = 0.6m;
        vm.FilmDp = 6;
        vm.BuildSurface().GetAwaiter().GetResult();
        var iface = vm.Document?.Summary();
        Check(iface is { } isum && isum.Molecules > 1 && vm.Title.Contains("film"), $"interface: {vm.Title} · {iface?.Atoms} atoms · {vm.SurfError} {vm.Status}");
        Check(vm.HoldOn && !vm.RelaxCompress, $"interface: the surface is held ({vm.HoldText}), no compression");
        vm.RelaxFtolD = 5;
        vm.Relax().GetAwaiter().GetResult();
        Check(vm.RelaxLog.Contains("UFF"), "interface relaxed with UFF, the surface held: " + vm.RelaxLog.Split('\n')[0]);

        // Nanostructure builder: a (5,5) tube in a natural-rubber matrix
        vm.OpenNano();
        vm.NanoKind = 1;
        vm.TubeN = 5;
        vm.TubeM = 5;
        vm.TubeLength = 12;
        vm.NanoMatrix = true;
        vm.MatrixChains = 4;
        vm.MatrixDp = 8;
        vm.MatrixDensity = 0.6m;
        vm.BuildNano().GetAwaiter().GetResult();
        var comp = vm.Document?.Summary();
        Check(comp is { } csum && csum.Molecules == 5 && vm.HoldOn, $"nanotube composite: {vm.Title} · {comp?.Atoms} atoms · {vm.NanoError} {vm.Status}");

        // Blend builder: NR / BR 70 : 30
        vm.OpenBlend();
        vm.BlendChains = 4;
        foreach (var r in vm.BlendRows) r.Dp = 8;
        vm.BuildBlend().GetAwaiter().GetResult();
        var bsum = vm.Document?.Summary();
        Check(vm.BlendRows.Count == 2 && bsum is { } blendSum && blendSum.Molecules >= 5 && vm.Title.Contains("blend"), $"blend: {vm.Title} · {bsum?.Molecules} chains · {vm.BlendError}");

        // Crystal builder: polyethylene (Pnam) to start, rutile's space group found from its CIF, a supercell built
        vm.OpenCrystal();
        Check(vm.CrystalGroup?.Number == 62 && vm.CrystalSites.Count == 3 && vm.CrystalGroups.Count >= 8 && vm.CrystalBFree,
              $"crystal: starts from {vm.CrystalGroup?.Title} with {vm.CrystalSites.Count} sites · {vm.CrystalGroups.Count} groups listed · {vm.CrystalError}");
        vm.CrystalQuery = "Fm-3m";
        Check(vm.CrystalGroups.Any(g => g.Number == 225), "crystal: search Fm-3m → " + string.Join(", ", vm.CrystalGroups.Select(g => g.Title)));
        vm.CrystalQuery = "14";
        Check(vm.CrystalGroups.Count >= 9 && vm.CrystalGroups.All(g => g.Number == 14), $"crystal: search 14 → {vm.CrystalGroups.Count} settings");
        vm.CrystalQuery = "";
        var rutileCif = vm.Crystals.FirstOrDefault(c => c.Id == "rutile")?.File ?? "";
        vm.ImportCrystalCif(rutileCif).GetAwaiter().GetResult();
        Check(vm.CrystalGroup?.Number == 136 && vm.CrystalSites.Count == 2 && !vm.CrystalBFree && vm.CrystalCFree,
              $"crystal: rutile.cif → {vm.CrystalGroup?.Title} · {vm.CrystalSites.Count} sites · {vm.CrystalLog} {vm.CrystalError}");
        vm.CrystalSupercell = "2 × 2 × 3";
        vm.BuildCrystal().GetAwaiter().GetResult();
        var xsum = vm.Document?.Summary();
        Check(xsum is { } xs && xs.Atoms == 72 && xs.Bonds == 144, $"crystal: rutile 2 × 2 × 3 → {xsum?.Atoms} atoms · {xsum?.Bonds} bonds · {vm.CrystalError} {vm.Status}");

        // Biomolecule builder: the board's peptide, a β-strand applied to a selection, built with the UFF clean-up
        vm.OpenBio();
        Check(vm.BioCells.Count == 30 && vm.BioSegments.Count == 2 && vm.BioSegments[0] == "α-helix 1–22", $"bio: {vm.BioCountText} · {vm.BioSegmentsText}");
        vm.SelectResidue(vm.BioCells[24], false);
        vm.SelectResidue(vm.BioCells[27], true);
        vm.BioType = 1;
        Check(vm.BioSegments.Count == 4 && vm.BioSegments[2] == "β-strand 25–28", "bio: strand on 25–28 → " + vm.BioSegmentsText);
        vm.BioNTerm = 2;
        vm.BioCTerm = 2;
        vm.BuildPeptide().GetAwaiter().GetResult();
        var psum = vm.Document?.Summary();
        Check(psum is { } pepSum && pepSum.Atoms > 300 && vm.Status.Contains("30 residues") && !vm.BioHasError, $"bio: built {psum?.Atoms} atoms · {vm.Status} {vm.BioError}");

        // Solvation builder: the peptide just built in TIP3P water with 0.15 M NaCl, 6 Å padding
        vm.OpenSolvation();
        vm.SolvWaterModel = "TIP3P";
        vm.SolvPadding = 6;
        Check(vm.SolvUseSolute && vm.SolvShape == 2 && vm.SolvCounts.Contains("Na⁺") && !vm.SolvHasError, $"solvation: plan {vm.SolvCounts} · {vm.SolvBoxText} {vm.SolvError}");
        vm.Solvate().GetAwaiter().GetResult();
        var wsum = vm.Document?.Summary();
        Check(wsum is { } solvSum && solvSum.Atoms > 2000 && vm.Title.Contains("water") && vm.SolvStages.All(st => st.IsDone), $"solvation: {vm.Title} · {wsum?.Atoms} atoms · {vm.SolvError} {vm.Status}");

        // Appearance: a style layer, colour by charge, a surface and R/S labels on the solvated peptide
        vm.SetModule(8);
        vm.AppearanceOpen = true;
        vm.AppTarget = 1;
        vm.AppExpression = "Molecule == 1";
        vm.AppStyle = 3;
        vm.AppTarget = 0;
        Check(vm.ShowAppearance && !vm.ShowStudioTabs && vm.AppLayers.Count == 1 && vm.AppChip.Contains("Molecule == 1"), $"appearance: {vm.AppChip}");
        if (vm.Document is { } adoc)
        {
            adoc.SetAppearance("{\"layers\":[{\"expression\":\"Molecule == 1\",\"style\":\"space_filling\"}],\"colour\":\"charge\",\"surface\":{\"kind\":\"excluded\",\"expression\":\"Molecule == 1\"}}");
            var look = System.Text.Json.Nodes.JsonNode.Parse(adoc.AppearanceInfo())!;
            var area = look["surface"]?["area"]?.GetValue<double>() ?? 0;
            Check(area > 500 && (look["styles"]?["space_filling"]?.GetValue<double>() ?? 0) > 100, $"appearance: excluded surface {area:0} Å² · {look["styles"]?.ToJsonString()}");
            var rs = System.Text.Json.Nodes.JsonNode.Parse(adoc.AtomLabels("rs"))!.AsArray().Count(x => x?.GetValue<string>() == "S");
            Check(rs >= 20, $"appearance: {rs} S centres labelled (L residues)");
            adoc.SetAppearance("{\"active\":false}");
        }
        vm.ResetAppearance();
        vm.AppearanceOpen = false;

        // Interactions & checks: a packed water box has H-bonds and a few O···O clashes; pushing them apart clears them
        {
            var (wbox, _) = CapsDocument.Solvate(null, "{\"shape\":0,\"edge\":15,\"ion_mode\":0,\"water_model\":\"TIP3P\"}", null, "water");
            using (wbox)
            {
                var ix = System.Text.Json.Nodes.JsonNode.Parse(wbox.Interactions("{}"))!;
                var clashes = ix["clashes"]?.GetValue<double>() ?? -1;
                var fixedJ = System.Text.Json.Nodes.JsonNode.Parse(wbox.Edit("{\"op\":\"clean\"}"))!;
                var after = System.Text.Json.Nodes.JsonNode.Parse(wbox.Interactions("{}"))!;
                Check(ix["hbonds"]?.GetValue<double>() > 10 && fixedJ["ok"]?.GetValue<bool>() == true && after["clashes"]?.GetValue<double>() < clashes,
                      $"interactions: {ix["hbonds"]} H-bonds · {clashes} clashes → {after["clashes"]} after clean-up · {fixedJ["error"]}");
            }
        }

        // Compute & remote: a host added, edited and removed (saved in the test settings file)
        {
            var hostsBefore = vm.Hosts.Count;
            vm.AddHost();
            vm.HostHostname = "hpc-login2.example.edu";
            vm.HostUser = "someone";
            var detail = vm.SelectedHost?.Detail ?? "";
            var added = vm.Hosts.Count;
            vm.RemoveHost();
            Check(added == hostsBefore + 1 && detail.Contains("someone@hpc-login2.example.edu") && vm.Hosts.Count == hostsBefore && vm.JobTemplate.Contains("#SBATCH"),
                  $"remote: {hostsBefore} → {added} → {vm.Hosts.Count} hosts · {detail}");
        }

        // Jobs: the runs above were recorded with their log and provenance
        Check(vm.Jobs.Any(j => j.Kind == "Analyze" && j.IsDone && j.Log.Count > 1 && j.Provenance.Any(f => f.Key == "sha256")) && File.Exists(MainViewModel.JobsFile),
              $"jobs: {vm.Jobs.Count} recorded ({string.Join(", ", vm.Jobs.Select(j => j.Id + " " + j.Status))})");

        // Progressive open (VisLoading): frame 0 first, the other frames read into the same document, then the file checks
        {
            var text = File.ReadAllText(Path.Combine(dir, "ps_melt.lammpstrj"));
            var chunks = text.Split("ITEM: TIMESTEP\n", StringSplitOptions.RemoveEmptyEntries);
            var longDump = Path.Combine(Path.GetTempPath(), "caps-selftest-long.lammpstrj");
            using (var w = new StreamWriter(longDump))
                for (int k = 0; k < 40; ++k) w.Write("ITEM: TIMESTEP\n" + (k * 1000) + "\n" + chunks[k % chunks.Length].Split('\n', 2)[1]);
            var progressiveBytes = MainViewModel.ProgressiveBytes;
            MainViewModel.ProgressiveBytes = 0;
            vm.SetModule(8);
            vm.OpenProgressive(longDump, Path.Combine(dir, "ps_melt.data")).GetAwaiter().GetResult();
            Avalonia.Threading.Dispatcher.UIThread.RunJobs();
            Check(vm.Frames == 40 && !vm.IsLoading && vm.LoadStages.Count == 5 && vm.LoadStages.All(st => st.State == "done") && vm.FileChecks.Count > 3,
                  $"progressive open: {vm.Frames} frames · " + string.Join(", ", vm.LoadStages.Select(st => st.Title + " " + st.State)));
            // cancelled at once: frame 0 stays open, reading stops after a frame
            var t = vm.OpenProgressive(longDump, Path.Combine(dir, "ps_melt.data"));
            vm.CancelLoad();
            t.GetAwaiter().GetResult();
            Avalonia.Threading.Dispatcher.UIThread.RunJobs();
            Check(vm.HasDocument && vm.Frames < 40 && vm.LoadStages[3].State == "stopped", $"progressive open cancelled: {vm.Frames} frames kept · {vm.LoadStages[3].Detail}");
            MainViewModel.ProgressiveBytes = progressiveBytes;
            File.Delete(longDump);
        }

        // Nothing open (VisEmpty): Analyze shows where to get a structure
        vm.CloseDocument();
        vm.SetModule(1);
        Check(vm.ShowEmpty && vm.KeyOpen.EndsWith("O") && vm.RecentFew.Count <= 4, $"empty state on Analyze · {vm.RecentCount} recent");
        vm.SetModule(8);
        Check(!vm.ShowEmpty, "Studio shows Start, not the empty state");
        vm.Open(Path.Combine(dir, "ps_melt.lammpstrj"), Path.Combine(dir, "ps_melt.data"));

        // Trajectory player: per-frame series of the melt, a chain's ends, smoothing
        if (vm.Document is { } tdoc)
        {
            var ts = System.Text.Json.Nodes.JsonNode.Parse(tdoc.TrajectorySeries("{\"dt_fs\":2}"))!;
            var cols = ts["columns"]?.AsArray().Select(x => x!.GetValue<string>()).ToList() ?? [];
            var nrows = ts["rows"]?.AsArray().Count ?? 0;
            var ends = ts["ends"]?.AsArray().Select(x => (int)x!.GetValue<double>()).ToArray() ?? [];
            Check(ts["ok"]?.GetValue<bool>() == true && nrows == vm.FrameMax + 1 && cols.Contains("Density (g/cm³)") && cols.Contains("Ree (Å)") && ends.Length == 2 && ends[0] != ends[1],
                  $"trajectory: {nrows} frames · {string.Join(", ", cols)} · ends {string.Join("–", ends)} {ts["error"]}");
            tdoc.SetSmoothing(3);
            tdoc.SetSmoothing(1);
            vm.OpenTrajectory();
            Check(vm.IsTrajectory && vm.TrajFrameText.StartsWith("frame 0"), "trajectory player opens: " + vm.TrajFrameText);
            vm.SetModule(8);
        }

        // Torsion scan: n-butane from SMILES (UFF), relaxed: trans lowest, gauche± above it
        using (var but = CapsDocument.BuildSmiles("CCCC", "uff", 1, 1, "butane").Doc)
        {
            var r = System.Text.Json.Nodes.JsonNode.Parse(but.TorsionScan("{\"atoms\":[0,1,2,3],\"step\":15,\"relax\":true}", null))!;
            var conf = r["conformers"]?.AsArray().Select(c => c!["state"]!.GetValue<string>()).ToList() ?? [];
            Check(r["ok"]?.GetValue<bool>() == true && conf.FirstOrDefault() == "trans" && conf.Contains("gauche+") && conf.Contains("gauche−"),
                  $"torsion: {string.Join(", ", conf)} · barrier {r["barrier"]} · {r["forcefield"]} {r["error"]}");
            but.TorsionShow(0);
            but.TorsionShow(-1);
        }

        // Editing (builder tools, Element picker): place, bond, delete, hydrogens, undo and redo on the single-frame melt
        vm.Open(Path.Combine(dir, "ps_melt.data"));
        {
            var n0 = vm.Document?.Summary().Atoms ?? 0;
            vm.BuildElement = "N";
            vm.EditTool = 1;
            vm.ToolClick(0);
            var n1 = vm.Document?.Summary().Atoms ?? 0;
            vm.EditTool = 3;
            vm.ToolClick((int)n1 - 1);
            var n2 = vm.Document?.Summary().Atoms ?? 0;
            vm.EditTool = 0;
            vm.UndoEdit(false);
            var n3 = vm.Document?.Summary().Atoms ?? 0;
            vm.UndoEdit(true);
            Check(n1 == n0 + 1 && n2 == n0 && n3 == n0 + 1 && vm.EditHistory.Count >= 1, $"edit: place N {n0}→{n1}, delete →{n2}, undo →{n3} · {string.Join(" / ", vm.EditHistory)} {vm.EditError}");
            vm.UndoEdit(false);
            vm.UndoEdit(false);
            var sel = System.Text.Json.Nodes.JsonNode.Parse(vm.Document!.Select("{\"mode\":\"smarts\",\"pattern\":\"c1ccccc1\"}"))!;
            var tac = System.Text.Json.Nodes.JsonNode.Parse(vm.Document!.Tacticity())!;
            Check(sel["count"]?.GetValue<double>() > 100 && tac["centres"]?.GetValue<double>() > 10, $"select: {sel["count"]} ring atoms · tacticity {tac["label"]} m {tac["m"]} r {tac["r"]}");
            vm.Document!.Select("{\"mode\":\"none\"}");
            // Selection & stereo: the panel, a SMARTS selection saved as a set, the melt made isotactic
            vm.SelectionOpen = true;
            vm.SelectMode = 0;
            vm.SelectPattern = "c1ccccc1";
            vm.RunSelect("replace");
            vm.SaveSelectionAsSet();
            vm.MakeTactic(true).GetAwaiter().GetResult();
            Check(vm.ShowSelectionPanel && vm.SelectedCount == 480 && vm.NamedSets.Count == 1 && vm.TacticityLabel == "isotactic" && vm.Dyads.All(d => d.IsMeso),
                  $"selection & stereo: {vm.SelectedChip} · {vm.NamedSets.Count} set · {vm.TacticityLabel} · {vm.DyadCounts} {vm.SelectError}");
            vm.UndoEdit(false);
            vm.SelectionOpen = false;
            // Fragment library: a methyl attached to a carbon (replacing one of its H), a water placed beside the melt
            vm.OpenFragments();
            var nf0 = vm.Document!.Summary().Atoms;
            var methyl = vm.FragmentTiles.Concat(vm.QuickFragments).First(f => f.Name == "Methyl");
            var carbon = Enumerable.Range(0, (int)nf0).First(i => vm.Document!.Atom(i).ElementSymbol == "C");
            vm.Pick(carbon);
            vm.UseFragment(methyl, false).GetAwaiter().GetResult();
            var nf1 = vm.Document!.Summary().Atoms;
            var water = vm.QuickFragments.First(f => f.Name == "Water");
            vm.UseFragment(water, false).GetAwaiter().GetResult();
            var nf2 = vm.Document!.Summary().Atoms;
            Check(vm.QuickFragments.Count == 9 && nf1 == nf0 + 3 && nf2 == nf1 + 3 && vm.FragmentCategories.Count >= 10,
                  $"fragments: {vm.FragmentCategories.Count} categories · methyl {nf0}→{nf1} · water →{nf2} · {vm.FragmentError}");
            vm.UndoEdit(false);
            vm.UndoEdit(false);
            vm.SetModule(8);
        }
        vm.Open(Path.Combine(dir, "ps_melt.lammpstrj"), Path.Combine(dir, "ps_melt.data"));

        // Split view: the melt beside its GROMACS copy, compared row by row
        vm.OpenSplit();
        vm.SetSplitB(Path.Combine(dir, "ps_melt.gro")).GetAwaiter().GetResult();
        Check(vm.IsSplit && vm.SplitHasB && vm.SplitRows.Count == 7 && vm.SplitRows[0].A == vm.SplitRows[0].B, $"split: {string.Join(" · ", vm.SplitRows.Select(r => $"{r.Property} {r.A}/{r.B}"))} {vm.SplitError}");
        vm.CloseSplitB();
        vm.SetModule(8);

        // Macro recorder: opening the melt and placing a water record as Python; a literal becomes a parameter; the script
        // replays with python3 through the caps package and this Studio's library
        {
            MainViewModel.MacroFolderOverride = Path.Combine(outDir, "caps-selftest-macros");
            if (Directory.Exists(MainViewModel.MacroFolderOverride)) Directory.Delete(MainViewModel.MacroFolderOverride, true);
            vm.OpenMacro();
            vm.NewMacro();
            vm.Recording = true;
            var data = Path.GetFullPath(Path.Combine(dir, "ps_melt.data"));
            vm.Open(data);
            var nm0 = vm.Document!.Summary().Atoms;
            vm.UseFragment(vm.QuickFragments.First(f => f.Name == "Water"), false).GetAwaiter().GetResult();
            vm.Record("print(\"atoms\", doc.atoms)");
            vm.Recording = false;
            var promo = vm.PromoteToParameter("\"Water\"");
            vm.RunMacro().GetAwaiter().GetResult();
            for (var i = 0; i < 20; i++) { Avalonia.Threading.Dispatcher.UIThread.RunJobs(); Thread.Sleep(25); }
            Check(vm.IsMacro && vm.RecordedCommands.Count == 3 && promo == null && vm.MacroParameters.Count == 1 && vm.MacroOutput.Contains($"atoms {nm0 + 3}") && vm.MacroOutput.Contains("done"),
                  $"macro: {vm.RecordedCommands.Count} recorded · params {string.Join(",", vm.MacroParameters.Select(p => p.Name + "=" + p.Default))} · {vm.MacroOutput.Replace('\n', ' ').Trim()}");
            vm.SetModule(8);
        }
        vm.Open(Path.Combine(dir, "ps_melt.lammpstrj"), Path.Combine(dir, "ps_melt.data"));

        // Import: the melt as extended XYZ (no bonds): bonds perceived with orders, molecules and the cell; None drops the
        // bonds; the whole file imported as a new document
        {
            var xyz = Path.Combine(dir, "ps_melt.xyz");
            vm.PreviewOpen(xyz);
            vm.WaitImport();
            var perceived = vm.ImportSummary;
            var fragOk = vm.ImportFragmentDoc?.Summary().Atoms > 10;
            vm.ImportBondMode = 2;
            vm.WaitImport();
            var none = vm.ImportSummary;
            vm.ImportBondMode = 0;
            vm.WaitImport();
            vm.ImportChecks = false;
            vm.ConfirmImport().GetAwaiter().GetResult();
            var si = vm.Document!.Summary();
            Check(perceived == "1,370 bonds · 10 molecules · 480 aromatic" && none.StartsWith("0 bonds") && fragOk && !vm.ImportOpen && si.Bonds == 1370 && si.Molecules == 10 && si.CellValid != 0,
                  $"import: {perceived} · none: {none} · document {si.Atoms} atoms {si.Bonds} bonds {si.Molecules} molecules {vm.ImportError}");
        }
        vm.Open(Path.Combine(dir, "ps_melt.lammpstrj"), Path.Combine(dir, "ps_melt.data"));

        // Export dialog: a 16-bit PNG with the provenance manifest and a measurement overlay; an animated PNG of the three
        // frames; a turntable of one frame as a PNG sequence
        {
            vm.OpenExportDialog();
            vm.ExportWidth = 320; vm.ExportHeight = 180;
            vm.ExportFormat = 1;
            var png = Path.Combine(outDir, "caps-selftest-export16.png");
            var wrote = vm.ExportDialogImage(png, (w, h) => { var l = new byte[w * h * 4]; for (int k = 3; k < 40 * 4; k += 4) { l[k - 3] = 255; l[k] = 255; } return l; }).GetAwaiter().GetResult();
            var meta = System.Text.Json.Nodes.JsonNode.Parse(CapsDocument.PngText(png))!;
            var man = System.Text.Json.Nodes.JsonNode.Parse(meta["caps:provenance"]!.GetValue<string>())!;
            var bytes = File.ReadAllBytes(png);
            Check(wrote != null && bytes[24] == 16 && man["schema"]?.GetValue<string>() == "caps-image/1.0" && man["source_sha256"]?.GetValue<string>().Length == 64 && man["frames"]?.GetValue<double>() == 3,
                  $"export image: {vm.ExportResult} · manifest {man["source"]} {man["source_sha256"]?.GetValue<string>()[..12]}…");
            vm.ExportTab = 1;
            var apng = Path.Combine(outDir, "caps-selftest-movie.png");
            vm.ExportDialogMovie(apng).GetAwaiter().GetResult();
            var ab = File.ReadAllBytes(apng);
            var actl = System.Text.Encoding.ASCII.GetString(ab).IndexOf("acTL", StringComparison.Ordinal);
            var frames = actl > 0 ? ab[actl + 4] << 24 | ab[actl + 5] << 16 | ab[actl + 6] << 8 | ab[actl + 7] : 0;
            vm.MovieSource = 1; vm.MovieTurnFrames = 12; vm.MovieFormat = 1;
            var seq = Path.Combine(outDir, "caps-selftest-turntable");
            if (Directory.Exists(seq)) Directory.Delete(seq, true);
            vm.Frame = 1;
            vm.ExportDialogMovie(seq).GetAwaiter().GetResult();
            var seqFiles = Directory.Exists(seq) ? Directory.GetFiles(seq, "frame_*.png").Length : 0;
            Check(frames == 3 && seqFiles == 12 && vm.Frame == 1, $"export movie: animated PNG {frames} frames · turntable {seqFiles} PNGs · {vm.ExportResult}");
            vm.ExportDialogOpen = false;
            vm.Frame = 0;
        }

        // Provenance: an edit recorded after the read, provSaved beside the file, read back on open, compared with itself
        {
            var src = Path.Combine(outDir, "caps-selftest-prov-src.data");
            File.Copy(Path.Combine(dir, "ps_melt.data"), src, true);
            if (File.Exists(src + ".provenance.json")) File.Delete(src + ".provenance.json");
            vm.Open(src);
            vm.Document!.Edit("{\"op\":\"place\",\"smiles\":\"O\",\"name\":\"Water\",\"resname\":\"WAT\"}");
            var provSaved = Path.Combine(outDir, "caps-selftest-prov.data");
            vm.Document.Save(provSaved);
            var side = File.Exists(provSaved + ".provenance.json");
            vm.Open(provSaved);
            vm.OpenProvenance();
            var engines = string.Join(",", vm.ProvSteps.Select(r => r.Engine));
            vm.CompareProvenanceWith(provSaved);
            var same = vm.ProvDiff.Count == 0 && vm.ProvDiffNotes.Count == 0;
            var bib = vm.ProvenanceBibtex();
            Check(side && engines == "io.read,edit.builder" && vm.ProvInputs.Count == 1 && same && vm.IsProvenance && bib.Length == 0,
                  $"provenance: {engines} · inputs {vm.ProvInputs.Count} · compare with itself: {vm.ProvDiffSummary}");
            vm.ClearProvenanceCompare();
            vm.SetModule(8);
        }
        vm.Open(Path.Combine(dir, "ps_melt.lammpstrj"), Path.Combine(dir, "ps_melt.data"));

        // Analyze focus pages: scattering with deuteration contrast, free volume with the voids drawn and cleared on leaving
        {
            vm.OpenScattering();
            vm.IsotopePattern = 2;
            vm.RunScattering().GetAwaiter().GetResult();
            var neutron = vm.Analyze.Results.FirstOrDefault(r => r.Id == "neutron");
            var deut = neutron?.Extra.FirstOrDefault(e => e.Key == "deuterated hydrogens").Value ?? 0;
            Check(vm.IsScattering && vm.ScatterXrayCurve.Length > 10 && vm.ScatterNeutronCurve.Length > 10 && deut > 0 && vm.ScatterLengths.Any(r => r.Key == "²H (D)"),
                  $"scattering: x-ray {vm.ScatterXrayCurve.Length} pts · neutron {vm.ScatterNeutronCurve.Length} pts · {deut} H deuterated · {vm.ScatterText[..Math.Min(60, vm.ScatterText.Length)]}");
            vm.IsotopePattern = 0;
            var exp = Path.Combine(outDir, "caps-selftest-self-overlay.csv");
            File.WriteAllLines(exp, new[] { "q,I" }.Concat(vm.ScatterXrayCurve.Select(p => $"{p.X.ToString(System.Globalization.CultureInfo.InvariantCulture)},{p.Y.ToString(System.Globalization.CultureInfo.InvariantCulture)}")));
            var why = vm.LoadExperiment(exp);
            Check(why == null && vm.ScatterExperiment.Length == vm.ScatterXrayCurve.Length && vm.HasExperiment, $"scattering overlay: {vm.ExperimentName} {why}");
            vm.ClearExperiment();
            vm.OpenFreeVolume();
            vm.RunFreeVolume().GetAwaiter().GetResult();
            var hadVoids = vm.FvHasVoids;
            vm.OpenMechanics();
            Check(hadVoids && vm.FvFfv.HasValue && vm.FvAccessible.HasValue && vm.FvPsd.Length > 3 && !vm.FvHasVoids && vm.IsMechanics,
                  $"free volume: FFV {vm.FvFfv.Value} · accessible {vm.FvAccessible.Value} · {vm.FvChip} · PSD {vm.FvPsd.Length} bins");
            vm.SetModule(8);
        }

        // Keyboard walk (VisAccess): atoms, bonds and molecules, announced
        {
            vm.FocusOn(40);
            var said = vm.Announcement;
            var start = vm.FocusAtom;
            vm.FocusBond();
            var bondedOk = vm.Document!.Bonded(start).Contains(vm.FocusAtom);
            vm.FocusStep(1);
            var m0 = vm.Document.Atom(start).Mol;
            var sameMol = vm.Document.Atom(vm.FocusAtom).Mol == m0;
            vm.FocusMolecule(1);
            var nextMol = vm.Document.Atom(vm.FocusAtom).Mol;
            vm.FocusSelect();
            var selOk = vm.Picked == vm.FocusAtom;
            vm.FocusStep(1);
            vm.FocusMeasure();
            Check(said.StartsWith("Atom 41,") && said.Contains("Nearest: atom") && bondedOk && sameMol && nextMol != m0 && selOk && vm.HasMeasure && vm.Announcement.StartsWith("Distance")
                  && vm.FocusRows.Count == 4 && vm.ViewOptions(64, 64, 1).Focus == vm.FocusAtom + 1,
                  $"keyboard walk: \"{said}\" · bond {bondedOk} · molecule {m0}→{nextMol} · \"{vm.Announcement}\"");
            vm.ReaderVerbosity = 0;
            Check(vm.Announcement.EndsWith($"molecule {vm.Document.Atom(vm.FocusAtom).Mol}."), "brief: " + vm.Announcement);
            vm.ReaderVerbosity = 1;
            vm.ClearFocus();
            vm.ClearSelection();
        }

        // Export › Figure (FigureBackground): journal size, true scale bar, overlay in the SVG, dpi in the PNG
        {
            vm.OpenFigure();
            vm.FigPreset = 0;
            vm.FigAspect = 0;
            var px = vm.FigPixels;
            var ov = vm.FigureOverlayFor(2, 2008, 1130, vm.Document!.ViewScale(vm.Camera, vm.FigureOptions(2, 2008, 1130, 1)));
            var svgPath = Path.Combine(outDir, "caps-selftest-figure.svg");
            vm.FigFormat = 1;
            vm.FigBackground = 2;
            var what = vm.ExportFigure(svgPath, (_, _, _, _, _, _) => { }, CapsStudio.Views.FigureDrawing.AddToSvg).GetAwaiter().GetResult();
            var svgText = File.ReadAllText(svgPath);
            var pngPath = Path.Combine(outDir, "caps-selftest-figure.png");
            vm.Document.ExportPng(vm.Camera, vm.FigureOptions(2, 400, 225, 1), pngPath);
            CapsStudio.Views.FigureDrawing.WriteDpi(pngPath, 600);
            var png = File.ReadAllBytes(pngPath);
            var phys = System.Text.Encoding.ASCII.GetString(png, 37, 4) == "pHYs" && ((png[41] << 24) | (png[42] << 16) | (png[43] << 8) | png[44]) == 23622;
            Check(px == (2008, 1130) && ov.BarPx > 0 && Math.Abs(ov.BarPx / ov.BarAngstrom - vm.Document.ViewScale(vm.Camera, vm.FigureOptions(2, 2008, 1130, 1))) < 1e-9
                  && svgText.Contains("caps-figure-overlay") && svgText.Contains(" Å</text>") && !svgText.Contains("<rect width=\"100%\"") && phys,
                  $"figure: {what} · pHYs {phys}");
            vm.FigFormat = 0;
            vm.SetModule(8);
        }

        // Render (RenderOverlays): tokens, output scale bar, guide inside the view, ambient occlusion in the options
        {
            vm.OpenRender();
            var label = vm.ResolveTokens("[Title] · frame [SourceFrame] · [Particles] atoms · ρ = [Density]", 0);
            var guide = vm.RenderGuide(1000, 700);
            var ropt = vm.RenderOptionsFor(1920, 1080, 1);
            var spec = vm.RenderSpecFor(1920, 1080, vm.Document!.ViewScale(vm.Camera, ropt), true, 0);
            Check(vm.IsRender && label.Contains("1,300 atoms") && label.Contains("0.386") && guide is { } gd && gd.W <= 1000 * 0.861 && gd.H <= 700 * 0.861
                  && (Math.Abs(gd.W - 860) < 1 || Math.Abs(gd.H - 602) < 1) && ropt.AmbientOcclusion == 1 && vm.ViewOptions(100, 100, 1).AmbientOcclusion == 1
                  && spec.BarPx > 0 && spec.Lines.Length == 2,
                  $"render: \"{label}\" · guide {guide?.W:F0} × {guide?.H:F0} in 1000 × 700 · bar {spec.BarLabel}");
            vm.SetModule(8);
            Check(vm.ViewOptions(100, 100, 1).AmbientOcclusion == 0, "the Studio view has no ambient occlusion");
        }

        // Analyze › Visualize (VisPipeline, DataInspector): steps in the core, the inspector reads their result
        {
            vm.OpenVisualize();
            vm.ClearPipeline();
            vm.AddStep("select_expression");   // Element == "H"
            vm.AddStep("delete_selected");     // above it, so it runs after
            var sel = vm.PipelineRows.Last().Summary;
            var particles = vm.PipeAttributes.FirstOrDefault(a => a.Key == "Particles")?.Value;
            var note = vm.InspectorNote;
            vm.InspectorFilter = "Molecule == 2";
            var filtered = vm.InspectorNote;
            var pixels = new byte[64 * 64 * 4];
            vm.Document!.Render(vm.Camera, vm.ViewOptions(64, 64, 1), pixels);
            vm.SetModule(8);
            var cleared = vm.Document.PipelineResult() == "";
            Check(vm.PipelineRows.Count == 2 && sel == "660 of 1300 selected" && particles == "640" && note.StartsWith("Rows 1–200 of 640") && filtered.StartsWith("Rows 1–64 of 64")
                  && cleared && vm.PipelineJson().Contains("\"delete_selected\""),
                  $"visualize: {sel} · {particles} particles · {note} · {filtered} · cleared on leaving {cleared}");
            vm.InspectorFilter = "";
        }

        // Export › Data (ExportData): every format previews and writes; the .gro reads back
        {
            var okFormats = new List<string>();
            foreach (var f in MainViewModel.ExportFormats)
            {
                var j = System.Text.Json.Nodes.JsonNode.Parse(vm.Document!.ExportPreview(f.Id, "{\"coeffs\":true}", 20))!;
                var lines = ((System.Text.Json.Nodes.JsonArray)j["lines"]!).Count;
                if (lines > 2 && (double?)j["bytes"] > 100) okFormats.Add(f.Id);
            }
            vm.OpenExport();
            vm.ExportFormatIndex = 2;
            var gro = Path.Combine(outDir, "caps-selftest-export.gro");
            vm.ExportNow(gro).GetAwaiter().GetResult();
            var back = CapsDocument.Open(gro);
            var sameAtoms = back.Summary().Atoms == vm.Document!.Summary().Atoms;
            back.Dispose();
            Check(okFormats.Count == MainViewModel.ExportFormats.Length && sameAtoms, $"export: {string.Join(", ", okFormats)} · .gro reads back {sameAtoms}");
            vm.ExportFormatIndex = 0;
            vm.SetModule(8);
        }

        // Analyze › Batch (BatchRun): the pipeline on several inputs, a failure kept as a row, results.csv written
        {
            vm.ClearPipeline();
            vm.AddStep("molecule_shape");
            vm.OpenBatch();
            var bad = Path.Combine(outDir, "caps-selftest-bad.data");
            File.WriteAllText(bad, "garbage\n");
            vm.BatchInputs.Clear();
            vm.AddBatchFiles([Path.Combine(dir, "ps_melt.data"), Path.Combine(dir, "water.pdb"), bad]);
            var outBatch = Path.Combine(outDir, "caps-selftest-batch");
            vm.BatchWorkers = 0;   // one worker: with no UI thread nothing serialises the row updates
            vm.RunBatch(outBatch).GetAwaiter().GetResult();
            var states = string.Join(",", vm.BatchInputs.Select(b => b.State));
            var csv = File.Exists(Path.Combine(outBatch, "results.csv")) ? File.ReadAllLines(Path.Combine(outBatch, "results.csv")) : [];
            Check(states == "done,done,failed" && csv.Length == 4 && csv[0].Contains("MoleculeShape.mean_rg") && vm.BatchInputs[0].Attributes["Particles"] == 1300,
                  $"batch: {states} · {csv.Length - 1} rows · {vm.BatchErrors}");
            File.Delete(bad);
            vm.ClearPipeline();
            vm.SetModule(8);
        }

        // Analyze › Compare (CompareCells): the same pipeline on two inputs, B − A per attribute
        {
            vm.ClearPipeline();
            vm.AddStep("molecule_shape");
            vm.OpenCompare();
            vm.SetCompareInput(true, Path.Combine(dir, "ps_melt.data")).GetAwaiter().GetResult();
            vm.SetCompareInput(false, Path.Combine(dir, "ps_melt.gro")).GetAwaiter().GetResult();
            var cmpRows = vm.CompareRows.ToDictionary(r => r.Metric);
            Check(vm.CompareHasB && cmpRows.TryGetValue("Particles", out var pr) && pr.Diff == "+0" && cmpRows.ContainsKey("MoleculeShape.mean_rg") && vm.CompareTables.Contains("Molecule shape"),
                  $"compare: {vm.CompareRows.Count} metrics · tables {string.Join(", ", vm.CompareTables)}");
            vm.ClearPipeline();
            vm.SetModule(8);
        }

        // Colour by (ColourBy): choosing a colouring sets the Visualize colour coding
        {
            vm.ClearPipeline();
            var tile = vm.ColourTiles.First(t => t.Property == "Charge");
            vm.ApplyColourTile(tile);
            var cc = vm.PipelineRows.FirstOrDefault(r => r.Type == "colour_coding");
            Check(vm.IsVisualize && cc != null && (string?)cc.Params["property"] == "Charge" && (string?)cc.Params["map"] == "diverging" && cc.Summary.Contains("diverging"),
                  $"colour by: {cc?.Summary}");
            vm.ClearPipeline();
            vm.SetModule(8);
        }

        // Figure bundle (FigureBundle): figure, data, pipeline, input and provenance; caps reproduce rebuilds it
        {
            vm.ClearPipeline();
            vm.AddStep("molecule_shape");
            vm.OpenBundle();
            vm.BundleInput = true;
            var zip = Path.Combine(outDir, "caps-selftest.caps-bundle.zip");
            var n = vm.ExportBundle(zip).GetAwaiter().GetResult();
            var entries = System.IO.Compression.ZipFile.OpenRead(zip).Entries.Select(e => e.FullName).ToList();
            Check(n >= 7 && entries.Contains("figure.png") && entries.Contains("data/molecules.csv") && entries.Contains("provenance.json") && entries.Any(x => x.StartsWith("input/")),
                  $"bundle: {n} files · {string.Join(", ", entries)}");
            vm.ClearPipeline();
            vm.SetModule(8);
        }

        // Open preview (OpenLammps, OpenGromacs): columns, types and frames before reading
        {
            var j = System.Text.Json.Nodes.JsonNode.Parse(CapsDocument.InspectFile(Path.Combine(dir, "ps_melt.lammpstrj"), Path.Combine(dir, "ps_melt.data")))!;
            var g = System.Text.Json.Nodes.JsonNode.Parse(CapsDocument.InspectFile(Path.Combine(dir, "ps_melt.gro"), null))!;
            Check((string?)j["format"] == "lammps-dump" && (double?)j["frames"] == 3 && ((System.Text.Json.Nodes.JsonArray)j["columns"]!).Count == 7
                  && ((System.Text.Json.Nodes.JsonArray)j["types"]!).Count == 4 && (string?)g["format"] == "gro" && (double?)g["atoms"] == 1300,
                  $"open preview: {j["format"]} · {j["frames"]} frames · {g["format"]} {g["atoms"]} atoms");
        }

        // Save pipeline (SavePipeline): YAML out, and back into Visualize
        {
            vm.ClearPipeline();
            vm.AddStep("coordination");
            vm.AddStep("colour_coding");
            vm.OpenSavePipeline();
            var yamlPath = Path.Combine(outDir, "caps-selftest.caps-pipeline.yaml");
            vm.SavePipelineYaml(yamlPath);
            var jsonBefore = vm.PipelineJson();
            vm.ClearPipeline();
            vm.LoadPipeline(yamlPath);
            Check(File.ReadAllText(yamlPath).StartsWith("caps_pipeline: 1") && vm.PipelineRows.Count == 2 && vm.PipelineJson() == jsonBefore && vm.PipelineYamlHash.Length == 64,
                  $"save pipeline: {vm.PipelineRows.Count} steps back from YAML · sha256 {vm.PipelineYamlHash[..12]}");
            vm.ClearPipeline();
            vm.SetModule(8);
        }

        // Pipeline groups (PipelineGroups): a group's switch turns all its steps off; folding hides them
        {
            vm.ClearPipeline();
            vm.AddStep("coordination");
            vm.PipeSelected!.Params["group"] = "Structure";
            vm.AddStep("cluster");   // joins the selected step's group
            vm.AddStep("colour_coding");
            vm.PipeSelected!.Params["group"] = "Look";
            vm.ApplyPipeline();
            var header = vm.PipelineRows.First(r => r.GroupHeader && r.Group == "Structure");
            header.GroupOn = false;
            var off = vm.PipelineRows.Where(r => r.Group == "Structure").All(r => !r.Enabled);
            header.GroupFolded?.Invoke(header);
            var folded = vm.PipelineRows.Where(r => r.Group == "Structure").All(r => !r.RowVisible);
            Check(header.GroupCount == 2 && off && folded && vm.PipelineJson().Contains("\"group\":\"Structure\""),
                  $"groups: Structure {header.GroupCount} steps · off {off} · folded {folded}");
            header.GroupFolded?.Invoke(header);
            vm.ClearPipeline();
            vm.SetModule(8);
        }

        // Close goes back to Start
        vm.SetModule(1);
        vm.CloseDocument();
        Check(vm.NoDocument && vm.IsStudio && vm.Title == "", "close: back to Start with no document");

        Console.WriteLine(fails == 0 ? "all checks passed" : $"{fails} check(s) failed");
        return fails == 0 ? 0 : 1;
    }
}
