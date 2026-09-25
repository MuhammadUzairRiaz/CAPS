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

        Check(Native.AbiVersion() == 18, "native ABI version 18");
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

        // Close goes back to Start
        vm.SetModule(1);
        vm.CloseDocument();
        Check(vm.NoDocument && vm.IsStudio && vm.Title == "", "close: back to Start with no document");

        Console.WriteLine(fails == 0 ? "all checks passed" : $"{fails} check(s) failed");
        return fails == 0 ? 0 : 1;
    }
}
