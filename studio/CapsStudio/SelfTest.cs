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

        Check(Native.AbiVersion() == 29, "native ABI version 29");
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
        {
            // the GPU view's scene and camera: the same atoms and bonds, and every sphere projected where the CPU image has it
            var gcam = new CapsCamera { Yaw = 0.9, Pitch = 0.3, Zoom = 1.3, PanX = 1.5, PanY = -2, Perspective = 1 };
            var gopt = vm.ViewOptions(800, 500, 1);
            var sc = vm.Document!.RenderScene(gopt);
            var fit = vm.Document!.ViewFit(gcam, gopt);
            vm.Document!.Render(gcam, gopt, new byte[800 * 500 * 4]);
            var proj = vm.Document!.ProjectAtoms(gcam, gopt, 1300);
            double worst = 0;
            for (var k = 0; k < sc.SphereId.Length; k++)
            {
                double dx = sc.Spheres[4 * k] - fit.Cx, dy = sc.Spheres[4 * k + 1] - fit.Cy, dz = sc.Spheres[4 * k + 2] - fit.Cz;
                var rx = dx * fit.CosYaw + dz * fit.SinYaw; var rz = -dx * fit.SinYaw + dz * fit.CosYaw;
                var ry = dy * fit.CosPitch - rz * fit.SinPitch; var rz2 = dy * fit.SinPitch + rz * fit.CosPitch;
                rx += fit.PanX; ry += fit.PanY;
                var kk = fit.Perspective != 0 ? fit.Dist / Math.Max(1e-3, fit.Dist - rz2) : 1;
                var i = sc.SphereId[k];
                worst = Math.Max(worst, Math.Max(Math.Abs(fit.W / 2 + rx * fit.Scale * kk - proj[3 * i]), Math.Abs(fit.H / 2 - ry * fit.Scale * kk - proj[3 * i + 1])));
            }
            Check(sc.SphereId.Length == 1300 && sc.CapsuleRgb.Length == 2 * 1370 && sc.LineRgb.Length == 12 && !sc.CpuOnly && worst < 1e-3,
                  $"GPU scene: {sc.SphereId.Length} spheres, {sc.CapsuleRgb.Length} half-bonds, {sc.LineRgb.Length} cell edges; projected as the CPU image to {worst:G2} px");
            // the view's own fit from the scene (used while a run holds the document) equals the core's
            double fitWorst = 0;
            foreach (var c in new[] { gcam, new CapsCamera { Yaw = -2.1, Pitch = 1.2, Zoom = 0.7 }, new CapsCamera { Yaw = 0.3, Pitch = -0.4, Zoom = 2.5, PanX = -3, Perspective = 1 } })
            {
                var fa = vm.Document!.ViewFit(c, gopt);
                var fb = sc.Fit(c, 800, 500);
                fitWorst = Math.Max(fitWorst, Math.Max(Math.Abs(fa.Scale - fb.Scale) / fa.Scale, Math.Max(Math.Abs(fa.Dist - fb.Dist) / fa.Dist,
                    Math.Max(Math.Abs(fa.Cx - fb.Cx) + Math.Abs(fa.Cy - fb.Cy) + Math.Abs(fa.Cz - fb.Cz), Math.Abs(fa.ZMin - fb.ZMin) + Math.Abs(fa.ZMax - fb.ZMax)))));
            }
            Check(fitWorst < 1e-5, $"GPU view's own camera fit equals the core's (worst {fitWorst:G2})");
        }
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
        var gaff = vm.Field.Library.ToList().FindIndex(x => x.Id == "gaff-amber25");
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
            vm.Field.ChargeMode = 2;   // Gasteiger–Marsili
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
        // Automatic charges: the library's OPLS-AA keeps charges on templates, not types; the default falls back to Gasteiger
        var opls = vm.Field.Library.ToList().FindIndex(x => x.Id == "opls2005");
        if (opls >= 0)
        {
            vm.Field.FfIndex = opls;
            vm.Field.ChargeMode = 0;
            vm.Field.Assign().GetAwaiter().GetResult();
            var rep = vm.Document!.FieldReport();
            Check(vm.Field.Assigned && vm.Field.Complete && rep.Contains("Gasteiger–Marsili charges were used instead", StringComparison.Ordinal),
                  $"Field OPLS-AA 2005 with automatic charges: {vm.Field.TypedText} · {vm.Field.MissingText} {vm.Field.Log}");
            vm.Field.Clear().GetAwaiter().GetResult();
        }
        // Coverage: OPLS-AA 2024 by type number types polystyrene completely and neutrally; the library check lists the others
        var opls24 = vm.Field.Library.ToList().FindIndex(x => x.Id == "oplsaa2024-moltemplate");
        if (opls24 >= 0)
        {
            vm.Field.FfIndex = opls24;
            vm.Field.ChargeMode = 0;
            vm.Field.Assign().GetAwaiter().GetResult();
            var complete24 = vm.Field.Assigned && vm.Field.Complete;
            vm.Field.CheckCoverage().GetAwaiter().GetResult();
            Check(complete24 && vm.Field.Alternatives.Any(a => a.Id == "pcff") && !vm.Field.HasUntypedGroups && !vm.Field.HasBalanceNote,
                  $"Field OPLS-AA 2024 by number: complete {complete24} · {vm.Field.TypedText} · {vm.Field.MissingText} · alternatives {string.Join(", ", vm.Field.Alternatives.Select(a => a.Id))} · {vm.Field.CoverageNote}");
            vm.Field.Clear().GetAwaiter().GetResult();
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
        {
            // compare states: frame 3 is frame 1 moved 1 Å along x (the sample pattern above), so the fit leaves nothing and no fit leaves 1 Å everywhere
            var fit = System.Text.Json.Nodes.JsonNode.Parse(vm.Document.CompareStates("{\"reference\":{\"kind\":\"frame\",\"index\":0},\"moving\":{\"kind\":\"frame\",\"index\":2},\"periodic\":\"no\"}"))!;
            var raw = System.Text.Json.Nodes.JsonNode.Parse(vm.Document.CompareStates("{\"reference\":{\"kind\":\"frame\",\"index\":0},\"moving\":{\"kind\":\"frame\",\"index\":2},\"fit\":\"none\",\"periodic\":\"no\"}"))!;
            var rf = fit["rmsd"]!["all"]!.GetValue<double>();
            var rr = raw["rmsd"]!["all"]!.GetValue<double>();
            vm.StatesOpen = true;
            var panel = vm.CompareRmsdAll;
            vm.StatesOpen = false;
            Check(rf < 1e-6 && Math.Abs(rr - (f2.X - f0.X)) < 1e-3 && panel.EndsWith("Å"), $"compare states: frame 3 on frame 1 RMSD {rf:0.######} Å after the fit, {rr:0.###} Å without · panel {panel}");
        }

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
        // the live view: snapshots of the running positions (the first at once), each a document of the cell's atoms
        var liveN = 0; var liveAtoms = 0L; var liveStats = "";
        cell.Md(mdOpts with { Steps = 50 }, null, (snap, stats) => { liveN++; liveAtoms = snap.Summary().Atoms; liveStats = stats; snap.Dispose(); });
        Check(liveN >= 1 && liveAtoms == cell.Summary().Atoms && liveStats.Contains("density"), $"live MD view: {liveN} snapshot(s) of {liveAtoms} atoms · {liveStats}");
        {
            // while a run holds the structure the window's own calls read a shadow of the shown frame and never wait
            var started = new ManualResetEventSlim();
            var stop = false;
            var run = Task.Run(() => cell.Md(mdOpts with { Steps = 10_000_000 }, (r, n) => { started.Set(); return !Volatile.Read(ref stop); }));
            started.Wait(20000);
            var sw = System.Diagnostics.Stopwatch.StartNew();
            var during = cell.Summary();
            var readMs = sw.Elapsed.TotalMilliseconds;
            var longRun = cell.LongRunning;
            Volatile.Write(ref stop, true);
            try { run.Wait(); } catch { }
            Check(longRun && readMs < 200 && during.Atoms == cell.Summary().Atoms && during.Frames == 1,
                  $"a run holds the structure: the window reads its shadow in {readMs:F1} ms ({during.Atoms} atoms)");
        }
        {
            // bond constraints: 2 fs steps with bonds to hydrogen held; the report says so, the LAMMPS line holds the same bonds
            var shk = cell.Md(mdOpts with { Dt = 2, Steps = 100, Constraints = 1 }, null);
            var line = cell.LammpsShake(1);
            Check(shk.Contains("SHAKE/RATTLE") && line.Contains("shake") && line.Contains(" m 1."), $"bond constraints: {shk.Split('\n').FirstOrDefault(l => l.Contains("SHAKE"))} · {line.Trim()}");
        }
        {
            // checkpoints: an NVE run stopped at step 350 continues from its checkpoint at 300 to the same end as one
            // uninterrupted run (same state and velocities; NVE has no random numbers)
            var nve = mdOpts with { Thermostat = 0, Steps = 1000, CheckpointEvery = 100, FrameEvery = 1000, NewVelocities = 0 };
            using var a1 = cell.Copy("checkpoint a");
            using var b1 = cell.Copy("checkpoint b");
            a1.Md(nve, null);
            var endA = a1.Summary().Frames;
            var posA = Enumerable.Range(0, 50).Select(i => a1.Atom(i)).ToArray();
            try { b1.Md(nve, (r, n) => r.Step < 350); } catch (InvalidOperationException) { }
            var ckpt = System.Text.Json.Nodes.JsonNode.Parse(b1.Checkpoint("info"))!;
            b1.Checkpoint("restore");
            b1.Md(nve with { Steps = 700, StepOffset = 300 }, null);
            var worst = Enumerable.Range(0, 50).Max(i => { var p = b1.Atom(i); return Math.Max(Math.Abs(p.X - posA[i].X), Math.Max(Math.Abs(p.Y - posA[i].Y), Math.Abs(p.Z - posA[i].Z))); });
            Check((string?)ckpt["ended"] == "stopped" && ckpt["step"]!.GetValue<double>() == 300 && worst < 1e-6,
                  $"checkpoint: stopped at 350, continued from step {ckpt["step"]} to 1000; positions agree with an uninterrupted run to {worst:0.0e+0} Å");
        }
        var cont = cell.Md(mdOpts with { Steps = 100 }, null);
        Check(cont.Contains("velocities taken"), "md: a second run continues with the same velocities");
        var respaRows = new List<CapsThermo>();
        var respa = cell.Md(mdOpts with { Dt = 2, Steps = 100, Thermostat = 0, Respa = 4, ThermoEvery = 10 }, (r, n) => { respaRows.Add(r); return true; });
        var band = respaRows.Max(r => r.Total) - respaRows.Min(r => r.Total);
        Check(respa.Contains("r-RESPA") && band < 0.03 * respaRows.Average(r => r.Kinetic), $"md: r-RESPA 2 fs × 4 inner steps, NVE energy band {band:F2} kcal/mol");
        // Equilibrate: a named protocol's text, a short custom protocol with convergence blocks, chain statistics.
        var l21 = CapsDocument.ProtocolText("larsen21", new CapsProtocolParams { TFinal = 300, TMax = 600, PFinal = 1, PMax = 49346.2, TimeScale = 1 });
        Check(l21.Split('\n', StringSplitOptions.RemoveEmptyEntries).Length == 21 && l21.Contains("P 49346.2 atm"), "equilibrate: Larsen 21-step text");
        var eqRows = 0;
        var eqOpts = new CapsEquilOpts { Dt = 1, Thermostat = 1, Barostat = 1, TauT = 100, TauP = 1000, Seed = 2, Cutoff = 10, Coulomb = 1, Tail = 1,
                                         FramePs = 0.5, ThermoPs = 0.1, UntilConverged = 1, BlockPs = 0.3, MaxBlocks = 3, TolDensity = 0.5, TolEnergy = 5, TolRg = 0.5 };
        var (eqConv, eqRep) = cell.Equilibrate("nvt 0.5 ps T 450 # warm\nnpt 0.5 ps T 300 P 1 atm # settle", eqOpts, (st, n, label, r) => { eqRows++; return true; });
        Check(eqConv && eqRows > 5 && eqRep.Contains("warm") && eqRep.Contains("check density"), $"equilibrate: 2 stages + blocks, {eqRows} rows, converged {eqConv}");
        {
            // Convergence board: Extend runs at the protocol's final NPT conditions; Accept now is recorded in the provenance
            var (tf, pf) = vm.FinalNpt("nvt 5 ps T 600\nnpt 800 ps T 310 P 2.5 atm # final\nnvt 1 ps T 900");
            cell.ProvenanceNote("{\"engine\":\"equilibrate.accepted\",\"summary\":\"accepted with 2 of 3 criteria met\",\"params\":{\"not met\":\"Rg\"}}");
            var prov = cell.Provenance();
            Check(tf == 310 && pf == 2.5 && prov.Contains("equilibrate.accepted") && prov.Contains("not met"), $"convergence: extend at T {tf} K, P {pf} atm; the acceptance is in the provenance");
        }
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

        // Architecture: a 4-arm natural-rubber star (one molecule; its atoms as the preview estimates), then linear again
        if (vm.PolymerLibrary.FirstOrDefault(e => e.Id == "P056") is { } nr)
        {
            while (vm.PolyUnits.Count > 1) vm.RemovePolyUnit(vm.PolyUnits.Last());
            vm.UseLibrary(nr, null);
            vm.GrowDpD = 10;
            vm.PolyArch = 2;
            vm.PolyArms = 4;
            var est = System.Text.RegularExpressions.Regex.Match(vm.PolyPreview, @"per molecule: 3 more arms · ([\d,]+) atoms");
            vm.BuildPolyPreview().GetAwaiter().GetResult();
            var sm = vm.PolyDoc?.Summary();
            Check(est.Success && sm != null && sm.Value.Molecules == 1 && sm.Value.Atoms == int.Parse(est.Groups[1].Value.Replace(",", "")),
                  $"star NR: {vm.PolyPreview.Replace('\n', ' ')} · built {sm?.Atoms} atoms in {sm?.Molecules} molecule {vm.PolyError}");
            vm.PolyArch = 0;
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
        // Analyze › Interface: the film's profile above the held quartz, this frame and over the relaxation's frames
        vm.OpenInterface();
        vm.RunInterface().GetAwaiter().GetResult();
        {
            string Row(string q) => vm.IfRows.FirstOrDefault(r => r.Quantity == q)?.ThisFrame ?? "—";
            Check(vm.IsInterfacePage && vm.IfRows.Count == 7 && Row("Surface top") != "—" && Row("Gap to the surface") != "—" && Row("Film plateau ρ") != "—" &&
                  vm.IfSeries.Any(x => x.Label == "film") && vm.IfSeries.Any(x => x.Label == "surface Si") && vm.IfGapBand != null,
                  $"interface page: {string.Join(" · ", vm.IfRows.Select(r => $"{r.Quantity} {r.ThisFrame}/{r.Trajectory}"))} · {vm.IfStatus}");
        }
        vm.SetModule(8);

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
        vm.WaitAppearance();   // the page's own applies run in the background: let them land before setting one directly
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
        vm.CloseAllStructures();
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

        // A united-atom force field (TraPPE-UA): assigning it folds the hydrogens on carbon into their carbons (undoable)
        {
            vm.UsePolystyreneInGrow();
            vm.GrowChainsD = 2; vm.GrowDpD = 6; vm.GrowDensityD = 0.3m;
            vm.GrowAssignField = false;
            vm.Grow().GetAwaiter().GetResult();
            var uaBefore = vm.Document!.Summary().Atoms;
            vm.Field.FfIndex = vm.Field.Library.ToList().FindIndex(e => e.Id == "trappe-ua");
            vm.Field.ChargeMode = 0;
            vm.Field.Assign().GetAwaiter().GetResult();
            var uaAfter = vm.Document!.Summary().Atoms;
            Check(vm.Field.Assigned && uaAfter < uaBefore && vm.Field.Notes.Any(n => n.Contains("united-atom")),
                  $"united-atom force field: {uaBefore} atoms → {uaAfter} sites · {vm.Field.Notes.FirstOrDefault(n => n.Contains("united-atom"))}");
            vm.UndoEdit(false);
            Check(vm.Document!.Summary().Atoms == uaBefore, $"united-atom conversion undone: {vm.Document!.Summary().Atoms} atoms");
            vm.GrowAssignField = true;
        }

        // Materials Studio .car/.mdf: a typed, charged structure keeps its types under a force field that has them
        // (ClayFF's SPC water o*, h*), exports complete, saves as .car and reopens with the same atoms and bonds
        {
            var car = Path.GetFullPath(Path.Combine(dir, "..", "tests", "data", "car", "water_pbc.car"));
            if (File.Exists(car))
            {
                var ffBack = vm.Field.FfIndex;
                var chBack = vm.Field.ChargeMode;
                vm.Open(car);
                var carAtoms = vm.Document!.Summary().Atoms;
                vm.Field.FfIndex = vm.Field.Library.ToList().FindIndex(e => e.Id == "inorganic-clay");
                vm.Field.ChargeMode = 3;   // keep the file's charges
                vm.Field.Assign().GetAwaiter().GetResult();
                Check(vm.Field.Assigned && vm.Field.Complete, $"ClayFF on a .car with its types: {(vm.Field.Complete ? "complete" : "incomplete")} · {vm.Field.FooterTyper}");
                var carOut = Path.Combine(outDir, "caps-selftest-water.car");
                vm.SaveDocument(carOut);
                vm.Open(carOut);
                var back = vm.Document!.Summary();
                Check(carAtoms == 6 && back.Atoms == 6 && back.Bonds == 4 && File.Exists(Path.ChangeExtension(carOut, ".mdf")),
                      $"Materials Studio .car/.mdf: {carAtoms} atoms typed by the file (o*, h*) under ClayFF, saved and reopened: {back.Atoms} atoms, {back.Bonds} bonds");
                vm.Field.FfIndex = ffBack;
                vm.Field.ChargeMode = chBack;
            }
        }

        // Coarse-grained: the builder's MARTINI DPPC template (12 beads, bond lengths from the force field)
        {
            var cgBack = vm.Module;
            vm.OpenBuilder();
            vm.CgFf = 0;
            var hasDppc = vm.CgTemplates.Contains("DPPC");
            vm.CgTemplate = "DPPC";
            vm.BuildBeadsMolecule().GetAwaiter().GetResult();
            Check(hasDppc && vm.MolDoc?.Summary().Atoms == 12, $"coarse-grained builder: {vm.CgTemplates.Count} MARTINI templates, DPPC {vm.MolDoc?.Summary().Atoms} beads · {vm.MolNotes}");
            // Martini 3: thymine (a virtual site at the centre of two beads) built from its template, opened and assigned
            vm.CgFf = Array.FindIndex(MainViewModel.CgForceFields, f => f.File == "martini3.json");
            var m3n = vm.CgTemplates.Count;
            vm.CgTemplate = "THYM";
            vm.BuildBeadsMolecule().GetAwaiter().GetResult();
            var thym = vm.MolDoc?.Summary().Atoms ?? 0;
            vm.OpenMoleculeInStudio();
            var ffBack = vm.Field.FfIndex;
            vm.Field.FfIndex = vm.Field.Library.ToList().FindIndex(e => e.Id == "martini3");
            vm.Field.Assign().GetAwaiter().GetResult();
            Check(m3n > 200 && thym == 5 && vm.Field.Assigned && vm.Field.Complete,
                  $"Martini 3 builder: {m3n} molecules, THYM {thym} beads, field {(vm.Field.Complete ? "complete" : "incomplete")} ({vm.Field.FfIndex}) {vm.Field.Log} · {string.Join(" · ", vm.Field.Notes.Take(3))}");
            vm.Field.FfIndex = ffBack;   // later checks assign their own force fields
            vm.SetModule(cgBack);
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
            // History & snapshots: a snapshot, two edits, an undo, an edit that branches, the snapshot back, a jump
            vm.HistoryOpen = true;
            var h0 = vm.Document!.Summary().Atoms;
            var b0 = vm.HistoryBranches.Count;
            vm.TakeSnapshot("as opened");
            vm.BuildElement = "O"; vm.EditTool = 1;
            vm.ToolClick(-1);
            vm.ToolClick(-1);
            b0 = vm.HistoryBranches.Count;   // the earlier tests' undone steps became a branch on the first edit
            vm.UndoEdit(false);
            var undone = vm.HistoryRows.Count(r => r.IsUndone);
            vm.ToolClick(-1);   // an edit after the undo: the undone step becomes a branch
            vm.EditTool = 0;
            var branches = vm.HistoryBranches.Count;
            var current = vm.HistoryRows.LastOrDefault()?.IsCurrent == true;
            vm.RestoreSnapshot(vm.Snapshots[0]);
            var back = vm.Document!.Summary().Atoms == h0;
            vm.JumpToStep(vm.HistoryRows.First(r => r.What.StartsWith("Restore")).Step - 1);
            var jumped = vm.Document!.Summary().Atoms == h0 + 2;
            Check(vm.ShowHistoryPanel && undone >= 1 && branches == b0 + 1 && current && back && jumped && vm.Snapshots.Count == 1,
                  $"history: {undone} undone · {b0}→{branches} branches · current {current} · snaps {vm.Snapshots.Count} thumb {vm.Snapshots.FirstOrDefault()?.Thumb != null} · {vm.EditError} · snapshot back {back} · jump {jumped} · {vm.HistorySummary}");
            vm.HistoryOpen = false;
        }
        vm.Open(Path.Combine(dir, "ps_melt.lammpstrj"), Path.Combine(dir, "ps_melt.data"));

        // Analyze › Diffusion: the MSD of the three-frame test trajectory; the Yeh–Hummer correction for water in a 3 nm box
        {
            var yh = MainViewModel.YehHummer(298, 0.89, 3.0);
            vm.OpenDiffusion();
            vm.RunDiffusion().GetAwaiter().GetResult();
            Check(vm.IsDiffusion && Math.Abs(yh - 2.32e-10) < 0.01e-10 && vm.DfMsd.Length >= 2 && vm.DfL > 0 && vm.DfCorrectionText.EndsWith("m²/s"),
                  $"diffusion: Yeh–Hummer {yh:0.000e0} m²/s · MSD {vm.DfMsd.Length} lags · L {vm.DfL} nm · {vm.DfStatus} {vm.DfWarning}");
            vm.SetModule(8);
        }

        // Select by query (SmartSelect): counted while typing, applied on Enter, errors in words, saved queries counted
        {
            vm.QueryOpen = true;
            vm.QueryText = "smarts \"c1ccccc1\"";
            var preview = vm.QueryInfo;
            vm.ApplyQuery();
            var selected = vm.SelectedCount;
            vm.QueryText = "within 5 sel";
            var bad = vm.QueryBad && vm.QueryInfo.Contains("of");
            vm.QueryText = "stereo * and chain 1";
            var stereo = vm.QueryInfo;
            var savedCount = vm.SavedQueries.FirstOrDefault(q => q.Name == "Aromatic rings")?.Count;
            Check(preview == "480 atoms · 80 rings" && selected == 480 && bad && savedCount == "480" && stereo.EndsWith("atoms"),
                  $"query: {preview} → {selected} selected · bad {bad} · stereo in chain 1: {stereo} · saved {savedCount}");
            vm.QueryOpen = false;
            vm.ClearDocSelection();
        }

        // Charges: Gasteiger on the melt, by group, applied (undoable), the Field's charges refused without an assignment
        {
            vm.OpenCharges();
            var groups = vm.ChargeGroups.Select(g => g.Group).ToList();
            var ok = vm.IsCharges && !vm.ChargeHasError && vm.ChargeNet.Contains("0.000000") && groups.Contains("C aromatic") && groups.Contains("H on sp³ C") && vm.ChargeHistogram.Length == 21;
            var q0 = vm.Document!.Atom(0).Charge;
            vm.ChargeMethod = 1;   // QEq differs from the Gasteiger charges the melt carries
            vm.ApplyCharges();
            var q1 = vm.Document!.Atom(0).Charge;
            vm.UndoEdit(false);
            var q2 = vm.Document!.Atom(0).Charge;
            vm.ChargeMethod = 2;
            var refused = vm.ChargeHasError && vm.ChargeError.Contains("Field");
            vm.ChargeMethod = 0;
            Check(ok && Math.Abs(q1 - q0) > 1e-6 && Math.Abs(q2 - q0) < 1e-12 && refused,
                  $"charges: {vm.ChargeNet} · largest {vm.ChargeMax} · {string.Join(", ", groups)} · applied {q0:0.000}→{q1:0.000}, undone {q2:0.000} · no field: {refused}");
            vm.AppColour = 0;
            vm.SetModule(8);
        }

        // Periodic box: the melt's crossing chains, pieces when wrapped, a crossing bond measured three ways, images, centring
        {
            vm.OpenPeriodic();
            var bond = vm.PbRows.FirstOrDefault(r => r.What.Contains("crosses"));
            var ok = vm.IsPeriodic && vm.PbChains == "10" && vm.PbCross == "7" && int.Parse(vm.PbPieces) > 10 && bond != null && bond.MinImage == bond.Whole &&
                     double.Parse(bond.Wrapped, System.Globalization.CultureInfo.InvariantCulture) > 10;
            vm.PbShow = 2;
            var popt = new CapsRenderOpts { Width = 120, Height = 90, Supersample = 1, Background = 0, Style = 0, ColourBy = 1, Outlines = 1, DepthCue = 1, ShowCell = 1,
                                           Highlight0 = -1, Highlight1 = -1, Highlight2 = -1, Highlight3 = -1 };
            var withImages = new byte[120 * 90 * 4];
            vm.Document!.Render(vm.Camera, popt, withImages);
            vm.PbShow = 1;
            var plain = new byte[120 * 90 * 4];
            vm.Document!.Render(vm.Camera, popt, plain);
            vm.Pick(0);
            var x0 = vm.Document!.Atom(0).X;
            vm.CentreOnSelection();
            var moved = Math.Abs(vm.Document!.Atom(0).X - x0) > 1e-6;
            vm.UndoEdit(false);
            vm.SetModule(8);
            Check(ok && !withImages.SequenceEqual(plain) && moved, $"periodic: {vm.PbChains} molecules · {vm.PbCross} cross · {vm.PbPieces} pieces · {bond?.What} {bond?.Wrapped}/{bond?.MinImage}/{bond?.Whole} · centred {moved}");
        }

        // Orientation: the estimator check (Bunn PE, S = 1 from 336 chords) and the melt's order parameter
        {
            var (checkS, chords) = MainViewModel.OrientationCheck();
            vm.OpenOrientation();
            vm.RunOrientation().GetAwaiter().GetResult();
            Check(Math.Abs(checkS - 1) < 1e-9 && chords == 336 && vm.IsOrientation && vm.OrS.HasValue && vm.OrCryst.HasValue && vm.OrHasZ,
                  $"orientation: check S {checkS} from {chords} chords · melt S {vm.OrS.Value} · f {vm.OrHerman.Value} · crystallinity {vm.OrCryst.Value}");
            vm.SetModule(8);
        }

        // Figure composer: panels from the project, checked at print size, SVG with vector plots; TIFF and PDF writers
        {
            var chips = vm.Analyze.Groups.SelectMany(g => g.Chips).ToList();
            foreach (var c in chips) c.IsOn = c.Id is "rdf" or "density";
            vm.Analyze.Run().GetAwaiter().GetResult();
            vm.OpenComposer();
            var sources = vm.ComposerSources.Count;
            var svgPath = Path.Combine(outDir, "caps-selftest-figure.svg");
            vm.ComposerPanel = 1;
            vm.ComposerLine = 0.3m;
            var thin = vm.ComposerChecks.Any(c => !c.Ok && c.Text.Contains("0.5 pt"));
            vm.ComposerLine = 1m;
            var curve = vm.Composition.Panels.Count(p => p.Kind == "curve" && p.HasData);
            vm.Status = vm.ExportComposerSvg(svgPath);
            var svg = File.ReadAllText(svgPath);
            var rgba = Enumerable.Repeat((byte)200, 40 * 20 * 4).ToArray();
            var tif = Path.Combine(outDir, "caps-selftest-figure.tiff");
            var pdf = Path.Combine(outDir, "caps-selftest-figure.pdf");
            Views.FigureFiles.WriteTiff(tif, rgba, 40, 20, 600);
            Views.FigureFiles.WritePdf(pdf, rgba, 40, 20, 72, 36);
            var tb = File.ReadAllBytes(tif);
            var pdfText = System.Text.Encoding.Latin1.GetString(File.ReadAllBytes(pdf));
            Check(vm.IsComposer && sources >= 3 && thin && curve >= 1 && svg.StartsWith("<svg") && svg.Contains("<polyline") && svg.Contains("width=\"7in\"") &&
                  tb[0] == 'I' && tb[2] == 42 && tb.Length == 8 + 2 + 12 * 12 + 4 + 6 + 16 + 40 * 20 * 3 && pdfText.StartsWith("%PDF-1.4") && pdfText.Contains("/MediaBox [0 0 72 36]") && pdfText.TrimEnd().EndsWith("%%EOF"),
                  $"composer: {sources} sources · {curve} curve panels · thin flagged {thin} · svg {svg.Length} chars · tiff {tb.Length} bytes");
            vm.SetModule(8);
        }

        // Row 17: chain statistics, the density calculator, SASA, the cell editor, display units, the molecule inspector
        {
            vm.OpenChainStats();
            vm.RunChainStats().GetAwaiter().GetResult();
            var chains = vm.CsRee.HasValue && vm.CsRg.HasValue && vm.CsRatio.HasValue && vm.CsCnCurve.Length > 3 && vm.CsReeHist.Length > 3;
            vm.OpenDensityCalc();
            var water = vm.DcSpecies.FirstOrDefault(sp => sp.Smiles == "O");
            var dens = water != null && water.N == "2,133" && vm.DcVolume.StartsWith("64,000 Å³") && vm.DcEdge.EndsWith("Å");
            vm.OpenSurfaceArea();
            vm.SaPoints = 60;
            vm.RunSurfaceArea().GetAwaiter().GetResult();
            var sasa = vm.SaGroups.Count >= 2 && vm.SaGroups[0].Part == "Whole structure" && vm.SaConvergence.Count == 5 && vm.SaTotal.EndsWith("Å²");
            vm.SetModule(8);
            vm.OpenCellEditor();
            var n0 = vm.Document!.Summary().Atoms;
            var vol0 = vm.CeVolume;
            vm.CeSupA = 2; vm.CeSupB = 1; vm.CeSupC = 1;
            vm.MakeSupercell();
            var doubled = vm.Document!.Summary().Atoms == 2 * n0 && Math.Abs(vm.Document!.Summary().CellA - 66) < 1e-6;
            vm.UndoEdit(false);
            vm.CeA = 34;
            vm.ApplyCell();
            var scaled = Math.Abs(vm.Document!.Summary().CellA - 34) < 1e-6;
            vm.UndoEdit(false);
            var cellOk = vol0 == "35937.00 Å³" && doubled && scaled && Math.Abs(vm.Document!.Summary().CellA - 33) < 1e-6;
            vm.SetModule(8);
            // units: a density result shown in kg/m³ under the SI-derived system, and back
            var card = vm.Analyze.Results.FirstOrDefault(r => r.Unit == "g/cm³");
            vm.SetUnitSystem = 1;
            var si = card == null || card.UnitText == "kg/m³";
            vm.SetUnitSystem = 0;
            var units = si && (card == null || card.UnitText == "g/cm³") && vm.UnitTable.Count == 8;
            vm.Pick(0);
            var mol = vm.MolInfoRows.FirstOrDefault(r => r.Key == "Formula")?.Value;
            Check(chains && dens && sasa && cellOk && units && mol == "C64H66" && vm.MolInfoRows.Any(r => r.Key == "SMILES" && r.Value.Length > 20 && r.Value.Contains("c")),
                  $"row 17: [{chains} {dens} {sasa} {cellOk} {units}] smiles {vm.MolInfoRows.FirstOrDefault(r => r.Key == "SMILES")?.Value} · vol {vm.DcVolume} edge {vm.DcEdge} · chains {vm.CsRee.Value}/{vm.CsRg.Value} ratio {vm.CsRatio.Value} · water N {water?.N} · SASA {vm.SaTotal} · cell {vol0} doubled {doubled} scaled {scaled} · units {units} · molecule {mol}");
        }

        // Row 18: polydispersity, copolymer, solvent screen, tacticity, blend phase diagram, electrostatics
        {
            vm.OpenPolydispersity();
            var pdRows = vm.PdRows.Count == 5 && vm.PdLengths.Length == vm.GrowChains && vm.PdK == "k = 1/(Đ − 1) = 10";
            vm.UsePdLengths();
            var pdUsed = vm.GrowPolydisperse && vm.GrowDispersityText.StartsWith("Đ 1.", StringComparison.Ordinal);
            vm.ClearPdLengths();
            // per-chain lengths reach the core: 3 styrene chains of 5, 7 and 9 units → 21 × 16 + 3 × 2 atoms
            var (pdDoc, pdReport) = CapsDocument.GrowChains("{\"units\":[{\"name\":\"styrene\",\"smiles\":\"[*]CC([*])C1=CC=CC=C1\"}],\"dp\":5,\"chain_dp\":[5,7,9]}",
                new CapsGrowOpts { Chains = 3, Dp = 0, Seed = 3, Density = 0.1, ContactScale = 0.8, Curve = 1 }, null, "pd");
            var pdAtoms = pdDoc.Summary().Atoms == 21 * 16 + 6 && pdDoc.AtomResidues().Max() == 9 && pdReport.Contains("5–9 units");
            pdDoc.Dispose();
            vm.OpenCopolymer();
            vm.BuildCoChain().GetAwaiter().GetResult();
            var coOk = vm.CoRows.Count == 4 && vm.CoRows[0].A == "0.510" && vm.CoAzeotrope is > 0.528 and < 0.53 && vm.CoSequence.Length == 80
                       && vm.CoDoc != null && vm.CoDoc.AtomResidues().Max() == 80;
            vm.CoM1 = 2; vm.CoM2 = 3;   // butadiene / acrylonitrile: NBR ratios filled in
            var nbr = vm.CoR1 == 0.30m && vm.CoR2 == 0.02m && vm.CoPresetNote.Contains("NBR");
            vm.CoM1 = 0; vm.CoM2 = 1;
            vm.OpenSolventScreen();
            var ssOk = vm.SsRows.Count == 8 && vm.SsRows.Count(r => r.Mismatch) == 1 && vm.SsRows.First(r => r.Mismatch).Name == "Acetone"
                       && vm.SsPolymerNames.Length == 3 && vm.SsExample.EndsWith("0.347", StringComparison.Ordinal);
            vm.OpenTacticityStats();
            vm.TsBuild().GetAwaiter().GetResult();
            var tsOk = vm.TsTriads.Count == 3 && vm.TsTriads[0].A == "0.250" && vm.TsDyads.Length == 198 && vm.TsModelPentads.Length == 10;
            vm.TsMeasured = "0.2401, 0.2058, 0.0441, 0.0882, 0.2058, 0.0882, 0.0378, 0.0081, 0.0378, 0.0441";   // Bernoulli P_m = 0.7
            vm.TsFitMeasured();
            var fitOk = vm.TsFit.Contains("P_m = 0.700", StringComparison.Ordinal) && vm.TsFit.Contains("consistent", StringComparison.Ordinal) && vm.TsPm == 0.7m;
            vm.TsBuild().GetAwaiter().GetResult();
            vm.OpenBlendPhase();
            var bpOk = vm.BpResults[0].Value == "0.01457" && vm.BpResults[2].Value == "433.9 K" && vm.BpResults[4].Value == "0.0385 · 0.9932" && vm.BpBinodal.Length > 100;
            vm.OpenElectrostatics();
            vm.EsMethod = 0; vm.EsCutoff = 12; vm.EsTolerance = 2; vm.EsSpacing = 1.2m;
            var esOk = vm.EsBeta == "0.2603 Å⁻¹" && vm.EsTable.Count == 3 && vm.EsTable.Count(r => r.Current) == 1 && vm.EsMesh.Count == 3;
            vm.SetModule(8);
            Check(pdRows && pdUsed && pdAtoms && coOk && nbr && ssOk && tsOk && fitOk && bpOk && esOk,
                  $"row 18: [{pdRows} {pdUsed} {pdAtoms} {coOk} {nbr} {ssOk} {tsOk} {fitOk} {bpOk} {esOk}] k {vm.PdK} · F1 {vm.CoRows.FirstOrDefault()?.A} · mismatches {vm.SsRows.Count(r => r.Mismatch)} · dyads {vm.TsDyads.Length} · fit {vm.TsFit} · χc {vm.BpResults.FirstOrDefault()?.Value} · β {vm.EsBeta} mesh {vm.EsMeshText}");
        }

        // χ by pair contacts: polystyrene against the solvent file (every row gets a χ and a verdict), and a blend fit
        {
            vm.OpenSolventScreen();
            vm.SsPolymer = 0;
            vm.SsComputeContacts().GetAwaiter().GetResult();
            var withC = vm.SsRows.Count(r => r.ChiC != "—");
            var tol = vm.SsRows.FirstOrDefault(r => r.Name == "Toluene");
            var water = vm.SsRows.FirstOrDefault(r => r.Name == "Water");
            var ssC = withC == 8 && tol != null && tol.AgreesC && water != null && water.AgreesC && vm.SsFailText.Contains("Pair contacts", StringComparison.Ordinal);
            vm.OpenBlendPhase();
            var bpNames = vm.BpUnitNames;
            vm.BpUnitA = Array.FindIndex(bpNames, n => n.Contains("styrene", StringComparison.OrdinalIgnoreCase));
            vm.BpUnitB = Array.FindIndex(bpNames, n => n.Contains("isoprene", StringComparison.OrdinalIgnoreCase) || n.Contains("natural rubber", StringComparison.OrdinalIgnoreCase));
            vm.BpFitContacts().GetAwaiter().GetResult();
            // strongly segregated at 300 K: two phases, each within powers of ten of pure (the binodal in logits)
            var coex = vm.BpResults[4].Value;
            var bpC = vm.BpUnitA >= 0 && vm.BpUnitB >= 0 && vm.BpSource == 1 && vm.BpB > 0 && vm.BpFitNote.Contains("χ(T) =", StringComparison.Ordinal)
                      && coex.Contains("E-", StringComparison.Ordinal) && vm.BpInputsChip.StartsWith("fitted", StringComparison.Ordinal);
            vm.BpA = -0.02m; vm.BpB = 15;
            var back = vm.BpSource == 0;
            vm.SetModule(8);
            Check(ssC && bpC && back, $"χ by pair contacts: {withC} solvents, toluene {tol?.ChiC} {tol?.PredictedC}, water {water?.ChiC}; blend {(vm.BpUnitA >= 0 ? bpNames[vm.BpUnitA] : "?")}/{(vm.BpUnitB >= 0 ? bpNames[vm.BpUnitB] : "?")} · {vm.BpFitNote.Split('·')[0]} · coexisting {coex}");
        }

        // Row 19: display styles, lens, add hydrogens, model resolution, live grow
        {
            vm.OpenDisplayStyles();
            var dsOk = vm.DsBackboneCaption.StartsWith("10 tubes through 160", StringComparison.Ordinal) && vm.DsAuto.Count == 4 && vm.DsAuto[0].Now;
            vm.DsStyle = 2;
            var bbStatus = vm.DisplayStatus == "Display: Backbone";
            vm.DsStyle = 0;
            vm.SetModule(8);
            vm.LensOpen = true;
            var lensOk = vm.LensOn && int.TryParse(vm.LensInsideCount.Replace("\u202F", ""), out var inLens) && inLens > 10 && inLens < 1300 && vm.DisplayStatus.Contains("lens", StringComparison.Ordinal);
            vm.LensOn = false;
            vm.LensOpen = false;
            // heavy atoms of a PS 4-mer: 34 H by the valence rules (aromatic C from the geometry)
            var frag = CapsDocument.Open(Path.Combine(dir, "ps_frag.pdb"));
            var plan = System.Text.Json.Nodes.JsonNode.Parse(frag.HydrogenPlan())!;
            var ahOk = (double?)plan["add"] == 34 && (double?)plan["aromatic_bonds"] == 24;
            var fragCopy = frag.Copy("copy");
            fragCopy.Edit("{\"op\":\"add_h\",\"atoms\":\"\"}");
            ahOk = ahOk && fragCopy.Summary().Atoms == 66 && frag.Summary().Atoms == 32;
            fragCopy.Dispose(); frag.Dispose();
            // a hydrogenated LAMMPS data file (no bond orders in it): the rings are aromatic from the geometry, nothing is missing
            using (var melt = CapsDocument.Open(Path.Combine(dir, "ps_melt.data")))
            {
                var mp = System.Text.Json.Nodes.JsonNode.Parse(melt.HydrogenPlan())!;
                var meltOk = (double?)mp["add"] == 0 && (double?)mp["aromatic_bonds"] == 480;
                Check(meltOk, $"add hydrogens on a complete LAMMPS data file: {mp["add"]} to add, {mp["aromatic_bonds"]} aromatic bonds from the geometry");
            }
            var res = System.Text.Json.Nodes.JsonNode.Parse(vm.Document!.ResolutionSummary("{\"per_bead\":5}"))!;
            var mOk = Math.Abs((double)res["all_atom"]!["mass"]! - (double)res["coarse_grained"]!["mass"]!) < 1e-6 && (double)res["united_atom"]!["sites"]! == 640 && (double)res["united_atom"]!["hydrogens"]! == 0;
            // a live grow: snapshots arrive while growing and the finished cell replaces them
            var snaps = 0;
            var (liveDoc, _) = CapsDocument.GrowChains("{\"units\":[{\"name\":\"styrene\",\"smiles\":\"[*]CC([*])C1=CC=CC=C1\"}],\"dp\":20}",
                new CapsGrowOpts { Chains = 6, Dp = 0, Seed = 4, Density = 0.3, ContactScale = 1.0, Curve = 1 }, null, "live", (d, stats) => { snaps++; d.Dispose(); });
            var liveOk = snaps >= 1 && liveDoc.Summary().Atoms == 6 * (20 * 16 + 2);
            liveDoc.Dispose();
            vm.SetModule(8);
            Check(dsOk && bbStatus && lensOk && ahOk && mOk && liveOk,
                  $"row 19: [{dsOk} {bbStatus} {lensOk} {ahOk} {mOk} {liveOk}] backbone '{vm.DsBackboneCaption}' · lens {vm.LensInsideCount} · H plan {plan["add"]} · UA sites {res["united_atom"]!["sites"]} · live snapshots {snaps}");
        }

        // Row 21 steps: expression counts, vector expressions, bonds against the file, replicas made real
        {
            var ec = System.Text.Json.Nodes.JsonNode.Parse(vm.Document!.ExpressionCount("Type == 2 && Position.Z > 13"))!;
            var exprOk = (double?)ec["count"] == 278 && (string?)ec["types"]?["2"] == "ca";
            vm.Document!.SetPipeline("{\"steps\":[{\"type\":\"compute_property\",\"name\":\"D\",\"expression\":\"norm(Position − MoleculeCOM(MoleculeIdentifier)) - DistanceToCOM\"},{\"type\":\"create_bonds\",\"mode\":\"pairs\",\"pairs\":\"C-C 1.70, C-H 1.25\"}]}");
            var pr = System.Text.Json.Nodes.JsonNode.Parse(vm.Document!.PipelineResult())!;
            double A(string k) => (pr["attributes"] as System.Text.Json.Nodes.JsonArray ?? []).Where(a => (string?)a?["name"] == k).Select(a => (double?)a?["value"] ?? double.NaN).DefaultIfEmpty(double.NaN).First();
            var bondsOk = A("CreateBonds.in_both") == 1370 && A("CreateBonds.cutoff_only") == 0 && A("CreateBonds.topology_only") == 0;
            var vecOk = ((string?)pr["steps"]?[0]?["summary"] ?? "").Contains("e-1", StringComparison.Ordinal) || ((string?)pr["steps"]?[0]?["summary"] ?? "").Contains(" 0 … 0", StringComparison.Ordinal);
            vm.Document!.SetPipeline("{\"steps\":[{\"type\":\"replicate\",\"nx\":2,\"ny\":2,\"nz\":2}]}");
            vm.Document!.PipelineResult();
            var real = vm.Document!.MaterializePipeline("real");
            var realOk = real.Summary().Atoms == 10400 && real.Summary().Molecules == 80;
            real.Dispose();
            vm.Document!.SetPipeline(null);
            Check(exprOk && bondsOk && vecOk && realOk, $"row 21 steps: [{exprOk} {bondsOk} {vecOk} {realOk}] {pr["steps"]?[0]?["summary"]} · {pr["steps"]?[1]?["summary"]}");
        }

        // Entanglements: primitive paths of the Kremer–Grest sample (LAMMPS gives N_e 68.1 by the modified S-coil)
        {
            using var kg = CapsDocument.Open(Path.Combine(dir, "kg_melt.data"));
            var ej = System.Text.Json.Nodes.JsonNode.Parse(kg.Analyze("entanglements", new CapsAnalyzeOpts { Last = -1, Stride = 1, Blocks = 5, Grid = 0.4, Qmax = 25, Dq = 0.02, FitFrom = 0.2, FitTo = 0.5, TimestepFs = 1 }, null))!;
            var ne = (double?)ej["properties"]?[0]?["value"] ?? double.NaN;
            Check(Math.Abs(ne - 68.1) < 0.05 * 68.1, $"entanglements: N_e {ne:F1} bonds on the Kremer–Grest sample (LAMMPS 68.1)");
        }

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
            var bondNotice = vm.Notices.FirstOrDefault(n => n.Key == "import.bonds");
            Check(bondNotice != null && bondNotice.Severity == "check" && bondNotice.Body.Contains("1,370 bonds") && bondNotice.Primary == "Review bonds",
                  $"system state: {bondNotice?.Title} · {bondNotice?.Body}");
            vm.Dismiss("import.bonds");
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

        // Pores: a graphite slit with methane packed between the walls, built into the Studio with the walls held
        {
            vm.OpenNano();
            vm.NanoKind = 3;
            vm.PoreType = 0;
            vm.PoreFluidIndex = 0;
            vm.PoreCount = 10;
            vm.BuildNano().GetAwaiter().GetResult();
            var sp = vm.Document!.Summary();
            var prov = System.Text.Json.Nodes.JsonNode.Parse(vm.Document.Provenance())!;
            var last = ((System.Text.Json.Nodes.JsonArray)prov["steps"]!).Last()!["engine"]!.GetValue<string>();
            Check(vm.IsStudio && sp.Atoms == 440 + 10 * 5 && vm.HoldOn && vm.NanoLog.Contains("packed inside the pore") && last == "nano.pore",
                  $"pore: {sp.Atoms} atoms · held {vm.HoldOn} · {last} · {vm.PoreChip}");
            vm.HoldOn = false;
            vm.NanoKind = 1;
        }
        vm.Open(Path.Combine(dir, "ps_melt.lammpstrj"), Path.Combine(dir, "ps_melt.data"));

        // First-run tour: starts once with the first structure (not in self-tests), steps, remembers it was done
        {
            vm.StartTour();
            var first = vm.TourCurrent.Region;
            vm.TourNext(); vm.TourNext(); vm.TourBack();
            var second = vm.TourCurrent.Region;
            for (int k = 0; k < 10 && vm.TourOpen; k++) vm.TourNext();
            vm.MaybeStartTour();
            Check(first == "Rail" && second == "Toolbar" && !vm.TourOpen && AppSettings.Load().TourDone, $"tour: {first} → {second} · done {AppSettings.Load().TourDone}");
        }

        // Compact layout: drawers are exclusive, the project shows only as a drawer
        {
            vm.SetModule(8);
            vm.Compact = true;
            var inspector = vm.InspectorShown && !vm.ProjectPanelShown;
            vm.ProjectDrawer = true;
            var swapped = vm.ProjectPanelShown && !vm.InspectorShown;
            vm.Compact = false;
            Check(inspector && swapped && vm.ProjectPanelShown && vm.InspectorShown, $"compact: inspector drawer {inspector} · project drawer {swapped}");
            vm.ProjectDrawer = false;
            vm.InspectorDrawer = true;
        }

        // Theory manual: every page's references resolve in CAPS's BibTeX table
        {
            vm.OpenManual("csvr");
            var unresolved = new List<string>();
            foreach (var item in vm.ManualNav.Where(n => n.Page != null))
            {
                vm.ShowManualPage(item.Page);
                unresolved.AddRange(vm.ManualRefs.Where(r => !r.Contains(" (")));
            }
            vm.OpenManual("csvr");
            var bib = vm.ManualBibtex();
            Check(vm.IsManual && vm.ManualCount >= 24 && unresolved.Count == 0 && bib.Contains("@article{bussi2007,") && vm.ManualSymbols.Count >= 4,
                  $"manual: {vm.ManualCount} pages · unresolved {string.Join(",", unresolved)} · {vm.ManualRefs.FirstOrDefault()?[..40]}");
            vm.SetModule(8);
        }

        // Project home: two seed-only replicas saved with provenance, the methods written for them
        {
            var proj = Path.Combine(outDir, "caps-selftest-project");
            if (Directory.Exists(proj)) Directory.Delete(proj, true);
            Directory.CreateDirectory(proj);
            foreach (var (name, seed) in new[] { ("cell_a.data", 11UL), ("cell_b.data", 12UL) })
            {
                using var d = CapsDocument.Open(Path.Combine(dir, "ps_melt.data"));
                d.Edit($"{{\"op\":\"place\",\"smiles\":\"O\",\"name\":\"Water {seed}\"}}");
                d.Save(Path.Combine(proj, name));
            }
            // make b differ from a only by a seed: rewrite its sidecar's first place name back and add rng
            var sa = System.Text.Json.Nodes.JsonNode.Parse(File.ReadAllText(Path.Combine(proj, "cell_a.data.provenance.json")))!;
            var sb2 = System.Text.Json.Nodes.JsonNode.Parse(sa.ToJsonString())!;
            sa["steps"]![1]!["rng"] = "mt19937-64 · seed 11";
            sb2["steps"]![1]!["rng"] = "mt19937-64 · seed 12";
            File.WriteAllText(Path.Combine(proj, "cell_a.data.provenance.json"), sa.ToJsonString());
            File.WriteAllText(Path.Combine(proj, "cell_b.data.provenance.json"), sb2.ToJsonString());
            vm.OpenProject(proj);
            for (int k = 0; k < 40; k++) { Thread.Sleep(20); }
            var zip = Path.Combine(outDir, "caps-selftest-project.zip");
            var shared = vm.ShareProject(zip);
            Check(vm.IsProject && vm.ProjectDocs.Count == 2 && vm.ProjectDocs.All(d => d.HasProvenance) && vm.ProjectMethods.Contains("2 independent replicas") && File.Exists(zip),
                  $"project: {vm.ProjectDocs.Count} structures · {vm.ProjectMethodsFor} · {shared}");
            vm.SetModule(8);
        }

        // Parameter sweep: 2 tacticities × 1 DP × 2 seeds of a short polystyrene, each run saved with provenance
        {
            vm.SweepFolder = Path.Combine(outDir, "caps-selftest-sweep");
            if (Directory.Exists(vm.SweepFolder)) Directory.Delete(vm.SweepFolder, true);
            vm.OpenSweep();
            vm.SweepSyn = false;
            vm.SweepDps = "3";
            vm.SweepSeeds = "1, 2";
            vm.SweepChains = 2;
            vm.RunSweep().GetAwaiter().GetResult();
            var files = Directory.GetFiles(vm.SweepFolder, "*.data").Length;
            var sides = Directory.GetFiles(vm.SweepFolder, "*.provenance.json").Length;
            var results = File.Exists(Path.Combine(vm.SweepFolder, "results.json")) && File.ReadAllText(Path.Combine(vm.SweepFolder, "results.json")).Contains("caps-sweep/1");
            Check(files == 4 && sides == 4 && results && vm.SweepResults.Count == 2 && vm.SweepResults.All(r => r.Seeds.StartsWith("2 / 2")) && (vm.SweepPolymer?.Name.StartsWith("Polystyrene") ?? false),
                  $"sweep: {files} cells · {sides} manifests · {string.Join(" | ", vm.SweepResults.Select(r => $"{r.Condition} Rg {r.Rg}"))} · {vm.SweepError}");
            vm.SetModule(8);
        }

        // Coarse-grained melt: 5 × 20 Kremer–Grest beads, the LAMMPS deck written
        {
            vm.CgChains = 5; vm.CgBeads = 20;
            var (kg, rep) = CapsDocument.KgBuild(vm.CgOptions(), "KG");
            var stem = Path.Combine(outDir, "caps-selftest-kg");
            kg.KgLammps(vm.CgOptions(), stem, 1000, 1000);
            var ok = File.Exists(stem + ".data") && File.ReadAllText(stem + ".in").Contains("bond_style fene") && kg.Summary().Atoms == 100 && kg.Summary().Bonds == 95;
            kg.Dispose();
            Check(ok, $"coarse-grained: 100 beads, 95 bonds, deck written · {rep}");
        }

        // Reaction template editor: the epoxy–amine templates pass their checks; a broken edit is flagged
        {
            vm.OpenTemplateEditor("epoxy_amine_primary");
            var pass = vm.TemplateChecksPass && vm.TemplateFormed.Count == 2 && vm.TemplateBroken.Count == 2;
            vm.TestTemplate();
            vm.TemplateText = vm.TemplateText + "\nbreak 4 1\n";
            var flagged = !vm.TemplateChecksPass && vm.TemplateHasError;
            Check(vm.IsTemplate && pass && flagged && vm.TemplateTestText.Contains("reactive sites"), $"template: pass {pass} · edit flagged {flagged} · {vm.TemplateTestText}");
            vm.SetModule(8);
        }

        // Level of detail: tiers counted by the renderer, the memory estimate, exports stay in full detail
        {
            vm.LodOn = true;
            vm.LodNear = 10; vm.LodFar = 20;
            var lodBuf = new byte[400 * 300 * 4];
            vm.Document!.Render(vm.Camera, vm.ViewOptions(400, 300, 1), lodBuf);
            var (near, mid, far, _) = vm.Document.RenderStats();
            var exportOpt = vm.ExportOptions(400, 300);
            vm.RefreshLod();
            Check(near > 0 && mid + far > 0 && exportOpt.LodNear == 0 && vm.LodMemory.Contains("B per atom"), $"level of detail: near {near} · mid {mid} · points {far} · {vm.LodMemory}");
            vm.LodOn = false;
            vm.LodNear = 40; vm.LodFar = 80;
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

        // Colour vision: three palettes × four visions, pairs sorted, the view preview changes pixels and never exports
        {
            vm.OpenColourVision();
            var sorted = vm.VisionPairs.Select(p => double.Parse(p.De, System.Globalization.CultureInfo.InvariantCulture)).ToList();
            var ok = vm.IsColourVision && vm.VisionPalettes.Count == 3 && vm.VisionPalettes.All(p => p.Rows.Length == 4) &&
                     vm.VisionPalettes[0].Rows[0].Swatches.Length == 8 && sorted.SequenceEqual(sorted.OrderBy(x => x));
            var vopt = new CapsRenderOpts { Width = 160, Height = 120, Supersample = 1, Background = 0, Style = 0, ColourBy = 1, Outlines = 1, DepthCue = 1, ShowCell = 1,
                                           Highlight0 = -1, Highlight1 = -1, Highlight2 = -1, Highlight3 = -1 };
            var va = new byte[160 * 120 * 4];
            var vb = new byte[160 * 120 * 4];
            vm.Document!.Render(vm.Camera, vopt, va);
            vm.VisionPreview = 1;
            vm.Document.Render(vm.Camera, vopt, vb);
            var differs = !va.SequenceEqual(vb);
            vm.VisionPreview = 0;
            Check(ok && differs, $"colour vision: {vm.VisionPalettes.Count} palettes · {vm.VisionSummary} · preview changes the view {differs}");
            vm.SetModule(8);
        }

        // Motion: the standard easing, reduced motion cuts, a flight lands on its target, F frames the selection
        {
            var ease = Math.Abs(Motion.Standard(0)) < 1e-9 && Math.Abs(Motion.Standard(1) - 1) < 1e-9 && Motion.Standard(0.5) > 0.75;
            vm.SetReduceMotion = "on";
            vm.Camera = new CapsCamera { Yaw = 1.0, Pitch = 0.1, Zoom = 3 };
            vm.ResetView();
            var cut = Math.Abs(vm.Camera.Yaw - 0.55) < 1e-12 && Math.Abs(vm.Camera.Zoom - 1) < 1e-12;
            vm.SetReduceMotion = "off";
            vm.Camera = new CapsCamera { Yaw = 1.0, Pitch = 0.1, Zoom = 3 };
            vm.ResetView();
            vm.StopFly();
            vm.FlyStep(0.5);
            var mid = vm.Camera.Yaw < 1.0 && vm.Camera.Yaw > 0.55 && vm.Camera.Zoom < 3 && vm.Camera.Zoom > 1;
            vm.FlyStep(1.0);
            var landed = Math.Abs(vm.Camera.Yaw - 0.55) < 1e-12 && Math.Abs(vm.Camera.Pitch - 0.40) < 1e-12 && Math.Abs(vm.Camera.Zoom - 1) < 1e-12;
            var one = vm.Document!.Focus(vm.Camera, Enumerable.Range(0, 130).ToArray());   // molecule 1
            vm.SetReduceMotion = "system";
            Check(ease && cut && mid && landed && one.Zoom > 1.2 && (Math.Abs(one.PanX) + Math.Abs(one.PanY)) > 0.5,
                  $"motion: easing {ease} · reduced cut {cut} · mid-flight {mid} · landed {landed} · focus zoom {one.Zoom:0.00}");
        }

        // Accessibility map: icon buttons are named from their tooltips' first clause
        {
            var n1 = Views.AccessibleNames.From("Export image or movie (⌘E)");
            var n2 = Views.AccessibleNames.From("Interactions & checks (design/boards/Interactions): H-bonds, contacts, clashes and fixes");
            var n3 = Views.AccessibleNames.From("Reset the view (R); drag in the view rotates");
            Check(n1 == "Export image or movie" && n2 == "Interactions & checks" && n3 == "Reset the view", $"accessible names: {n1} | {n2} | {n3}");
        }

        // Update available: a local feed with a newer release opens the dialog, results-altering notes first
        {
            var feed = Path.Combine(outDir, "caps-selftest-release.json");
            File.WriteAllText(feed, new System.Text.Json.Nodes.JsonObject
            {
                ["tag_name"] = "v9.9.0", ["html_url"] = "https://github.com/MuhammadUzairRiaz/CAPS/releases/tag/v9.9.0",
                ["body"] = "## New\n- Trajectory export to XTC\n## Changes that can alter results\n- PME tolerance 1e-5 → 1e-6\n## Fixed\n- A fix",
            }.ToJsonString());
            Environment.SetEnvironmentVariable("CAPS_UPDATE_FEED", feed);
            vm.CheckForUpdates().GetAwaiter().GetResult();
            var first = vm.UpdateSections.FirstOrDefault();
            var ok = vm.UpdateOpen && vm.UpdateVersion == "9.9.0" && vm.UpdateSections.Count == 3 && first!.AltersResults && first.Items[0].StartsWith("PME");
            vm.CloseUpdate();
            File.WriteAllText(feed, "{\"tag_name\":\"v0.0.1\",\"body\":\"\"}");
            vm.CheckForUpdates().GetAwaiter().GetResult();
            Environment.SetEnvironmentVariable("CAPS_UPDATE_FEED", null);
            Check(ok && !vm.UpdateOpen && vm.UpdateText.Contains("latest"), $"update: {vm.UpdateTitle} · {vm.UpdateSections.Count} sections · then {vm.UpdateText}");
        }

        // Start › From a recipe: a small polyethylene recipe runs, exports beside itself and opens as the document
        {
            var rdir = Path.Combine(outDir, "caps-selftest-recipe");
            if (Directory.Exists(rdir)) Directory.Delete(rdir, true);
            Directory.CreateDirectory(rdir);
            var recipe = Path.Combine(rdir, "pe.yaml");
            File.WriteAllText(recipe, "recipe: 1\nname: pe\nbuild:\n  polymer: { smiles: \"*CC*\", dp: 4, chains: 2 }\ntype: { forcefield: default }\n" +
                                      "grow: { density: 0.5, seed: 1 }\nrelax: { fmax: 5 }\nexport: [lammps]\n");
            vm.RunRecipeFile(recipe).GetAwaiter().GetResult();
            var ok = File.Exists(Path.Combine(rdir, "pe.data")) && File.Exists(Path.Combine(rdir, "pe.data.provenance.json")) && vm.Document?.Summary().Atoms == 52;
            Check(ok && vm.RecipeLog.Contains("[4/5] relax") && vm.Status.Contains("done"), $"recipe: {vm.Status} · {vm.RecipeLog.Replace('\n', '|')}");
            // Grow › Save recipe runs as saved (a small PS cell)
            vm.GrowChains = 3; vm.GrowDp = 4;
            File.WriteAllText(recipe, vm.GrowRecipe());
            vm.RunRecipeFile(recipe).GetAwaiter().GetResult();
            var grownAtoms = vm.Document?.Summary().Atoms ?? 0;
            Check(grownAtoms == 3 * (4 * 16 + 2) && File.Exists(Path.Combine(rdir, "PS_3x4.data")) && vm.GrowPython().Contains("caps.polymer(\"*CC(*)c1ccccc1\", dp=4, chains=3"),
                  $"grow recipe: {grownAtoms} atoms · {vm.Status}");
            // Analyze › Glass transition on that small cell: two replicas of a short cooling scan, pooled
            {
                var (gf0, gt0, gs0, gp0, ge0) = (vm.Analyze.TgFromD, vm.Analyze.TgToD, vm.Analyze.TgStepD, vm.Analyze.TgPsD, vm.Analyze.EqPsD);
                vm.Analyze.TgFromD = 450; vm.Analyze.TgToD = 250; vm.Analyze.TgStepD = 50; vm.Analyze.TgPsD = 2; vm.Analyze.EqPsD = 1;
                vm.GtReplicas = 2;
                vm.OpenGlass();
                var schedule = vm.GtSchedule.Length;
                vm.RunGlass().GetAwaiter().GetResult();
                Check(vm.IsGlass && schedule == 2 + 2 * 5 && vm.GtTemperaturesText == "5" && vm.GtPoints.Length == 5 && vm.GtErrors.Length == 5 &&
                      vm.GtRateText == "2.5 × 10¹³ K/s" && vm.GlassRecipe("/tmp/x.data").Contains("properties: [tg]"),
                      $"glass transition: {vm.GtPoints.Length} temperatures · Tg {vm.GtTg.Value} ({vm.GtTg.Caption}) · {vm.GtRateText} · {vm.GtStatus}");
                (vm.Analyze.TgFromD, vm.Analyze.TgToD, vm.Analyze.TgStepD, vm.Analyze.TgPsD, vm.Analyze.EqPsD) = (gf0, gt0, gs0, gp0, ge0);
                vm.GtReplicas = 1;
                vm.SetModule(8);
            }
            // Jobs › Recipes: templates checked, a broken edit flagged, saved and duplicated with their hash, run with the hash in provenance
            {
                var rfolder = Path.Combine(outDir, "caps-selftest-recipes");
                if (Directory.Exists(rfolder)) Directory.Delete(rfolder, true);
                Environment.SetEnvironmentVariable("CAPS_RECIPES", rfolder);
                vm.OpenRecipes();
                var first = vm.SelectedRecipe;
                var templates = vm.Recipes.Count;
                var valid = vm.RecipeOk && vm.RecipeStages.Count == 8 && vm.RecipeSchedule.Count == 21 && vm.RecipeT.Length == 42;
                var good = vm.RecipeText;
                vm.RecipeText = good.Replace("export: [lammps, gromacs]", "export: [lammps, tiff]");
                var flagged = !vm.RecipeOk && vm.RecipeStages.Any(st => !st.Ok) && vm.RecipeValid.Contains("tiff");
                vm.RecipeText = good;
                vm.SaveRecipe();
                vm.DuplicateRecipe();
                var files = Directory.GetFiles(rfolder, "*.yaml").Length;
                vm.SelectedRecipe = vm.Recipes.First(r => r.Name == "Quick molecule");
                vm.RunSelectedRecipe().GetAwaiter().GetResult();
                var prov = vm.Document?.Provenance() ?? "";
                Check(templates == 4 && first != null && valid && flagged && files == 2 && vm.Document?.Summary().Atoms == 9 && prov.Contains("recipe.run") && prov.Contains(vm.RecipeSha),
                      $"recipes: {templates} · valid {valid} · flagged {flagged} · {files} files · {vm.Status}");
                Environment.SetEnvironmentVariable("CAPS_RECIPES", null);
                vm.SetModule(8);
            }
            // a broken recipe says why, with its exit code
            File.WriteAllText(recipe, "build: {molecule: CCO}\nbogus: 1\n");
            vm.RunRecipeFile(recipe).GetAwaiter().GetResult();
            Check(vm.Status.Contains("exit 2") && vm.Status.Contains("bogus"), $"recipe error: {vm.Status}");
        }

        // Relax › distance restraints: two picked carbons of different chains pulled to 4 Å
        {
            vm.Open(Path.Combine(dir, "ps_melt.data"));
            var n = vm.Document!.Summary().Atoms;
            int ra = -1, rb = -1;
            double far = 0;
            for (int i = 0; i < 40; ++i)
                for (int j = 200; j < n; j += 7)
                {
                    var ai = vm.Document.Atom(i); var aj = vm.Document.Atom(j);
                    if (ai.ElementSymbol != "C" || aj.ElementSymbol != "C" || ai.Mol == aj.Mol) continue;
                    var d = vm.Document.Measure([i, j]);
                    if (d > far && d < 8) (far, ra, rb) = (d, i, j);
                }
            vm.Pick(ra);
            vm.Pick(rb, true);
            var canAdd = vm.CanAddRestraint;
            vm.AddMeasuredRestraint();
            vm.RelaxRestraints[0].R0D = 4;
            vm.RelaxRestraints[0].KD = 50;
            vm.RelaxFtolD = 2;
            vm.Relax().GetAwaiter().GetResult();
            var after = vm.Document!.Measure([ra, rb]);
            Check(canAdd && vm.RelaxLog.Contains("target 4.000 Å") && Math.Abs(after - 4) < 0.3,
                  $"relax restraint: {far:F2} → {after:F2} Å (target 4) · {vm.RelaxLog.Split('\n').FirstOrDefault(l => l.StartsWith("restraint"))}");
            vm.RemoveRestraint(vm.RelaxRestraints[0]);
        }

        // Dynamics › Export › GROMACS: the deck switches engine, the files are written and (when gmx is installed) grompp accepts them
        {
            vm.Open(Path.Combine(dir, "ps_melt.data"));
            vm.SetModule(3);
            vm.MdEnsemble = 2;
            vm.MdRespa = 1;
            vm.MdGromacs = true;
            vm.PreflightNow().GetAwaiter().GetResult();
            var deck = vm.MdDeck;
            var gdir = Path.Combine(outDir, "gmx-selftest");
            Directory.CreateDirectory(gdir);
            vm.SaveGromacs(gdir);
            var files = new[] { "system.top", "system.gro", "system.mdp" }.All(f => File.Exists(Path.Combine(gdir, f)));
            var gmx = new[] { "/opt/homebrew/bin/gmx", "/usr/local/bin/gmx", "/usr/bin/gmx" }.FirstOrDefault(File.Exists);
            var grompp = "gmx not installed: not run";
            if (gmx != null)
            {
                var psi = new System.Diagnostics.ProcessStartInfo(gmx, "-quiet grompp -f system.mdp -c system.gro -p system.top -o system.tpr -maxwarn 10")
                    { WorkingDirectory = gdir, RedirectStandardOutput = true, RedirectStandardError = true };
                using var pr = System.Diagnostics.Process.Start(psi)!;
                var err = pr.StandardError.ReadToEnd();
                pr.WaitForExit();
                grompp = pr.ExitCode == 0 ? "grompp ok" : "grompp failed: " + err.Split('\n').LastOrDefault(l => l.Contains("ERROR") || l.Contains("rror")) ;
            }
            Check(deck.Contains("coulombtype") && deck.Contains("C-rescale") && deck.Contains("mts-level2-factor        = 2") && files && !grompp.StartsWith("grompp failed"),
                  $"GROMACS export: deck {(deck.Contains("coulombtype") ? "ok" : deck.Split('\n')[0])}, files {files}, {grompp}");
            vm.MdGromacs = false;
            vm.MdRespa = 0;
            vm.MdEnsemble = 1;
            vm.PreflightNow().GetAwaiter().GetResult();
            Check(vm.MdDeck.Contains("read_data"), "LAMMPS deck back");
        }

        // View tools: the view-plane fit behind Move, lasso selection, a move with undo, a pinned distance
        {
            // an orthographic view: screen = 12 (r·w, u·w) + (400, 300); a drag of (30, −20) px is a move in the r–u plane
            double[] r = [0.6, 0.8, 0], u = [-0.48, 0.36, 0.8];
            var pts = new List<(double, double, double, double, double)>();
            var rng = new Random(3);
            for (var k = 0; k < 20; ++k)
            {
                double x = rng.NextDouble() * 20, y = rng.NextDouble() * 20, z = rng.NextDouble() * 20;
                pts.Add((x, y, z, 400 + 12 * (r[0] * x + r[1] * y + r[2] * z), 300 + 12 * (u[0] * x + u[1] * y + u[2] * z)));
            }
            var dw = MainViewModel.ScreenToWorld(pts, 30, -20)!;
            double sx = 12 * (r[0] * dw[0] + r[1] * dw[1] + r[2] * dw[2]), sy = 12 * (u[0] * dw[0] + u[1] * dw[1] + u[2] * dw[2]);
            double[] nrm = [r[1] * u[2] - r[2] * u[1], r[2] * u[0] - r[0] * u[2], r[0] * u[1] - r[1] * u[0]];
            var off = nrm[0] * dw[0] + nrm[1] * dw[1] + nrm[2] * dw[2];
            var fitOk = Math.Abs(sx - 30) < 1e-6 && Math.Abs(sy + 20) < 1e-6 && Math.Abs(off) < 1e-6;

            vm.Open(Path.Combine(dir, "ps_melt.data"));
            vm.LassoSelect(Enumerable.Range(0, 10).ToArray(), false);
            var lassoOk = vm.SelectedCount == 10;
            vm.ClearDocSelection();
            var x0 = vm.Document!.Atom(0).X;
            vm.TranslateAtoms([0, 1, 2], [1.5, 0, 0]);
            var moved = vm.Document!.Atom(0).X - x0;
            vm.UndoEdit(false);
            var back = vm.Document!.Atom(0).X - x0;
            vm.Pick(0);
            vm.Pick(5, true);
            vm.PinMeasurement();
            var pinOk = vm.Monitors.Count == 1 && vm.Monitors[0].Value.EndsWith(" Å") && vm.Monitors[0].Kind == "d";
            vm.UnpinMonitor(vm.Monitors[0]);
            Check(fitOk && lassoOk && Math.Abs(moved - 1.5) < 1e-9 && Math.Abs(back) < 1e-9 && pinOk && vm.Monitors.Count == 0,
                  $"view tools: [{fitOk} {lassoOk} {moved:F2} {back:F2} {pinOk}] move {dw[0]:F3} {dw[1]:F3} {dw[2]:F3} Å");
        }

        // React › REACTER-style: C–C crosslinks checked during one NVT run on the polystyrene melt
        {
            vm.Open(Path.Combine(dir, "ps_melt.data"));
            vm.RxSet = 0;
            vm.RxDuringMd = true;
            vm.RxMdPsD = 0.05m;
            vm.RxCyclesD = 2;
            vm.RxPerCycleD = 2;
            vm.RxTempD = 400;
            vm.RunReact().GetAwaiter().GetResult();
            Check(vm.RxLog.Contains("REACTER-style") && vm.Document!.Summary().Atoms < 1300, $"react during MD: {vm.RxLog.Split('\n')[0]} · {vm.Document!.Summary().Atoms} atoms");
            vm.RxDuringMd = false;
        }

        // Grow › Add solvent / gas: toluene into the free space of a small polystyrene cell
        {
            vm.UsePolystyreneInGrow();
            vm.GrowChainsD = 3;
            vm.GrowDpD = 6;
            vm.GrowDensityD = 0.3m;
            vm.GrowSmallPick = vm.GrowSmallChoices.FindIndex(c => c.Name == "Toluene");
            vm.AddGrowSmall();
            vm.GrowSmall[0].CountD = 5;
            vm.Grow().GetAwaiter().GetResult();
            var mols = vm.Document?.Summary().Molecules ?? 0;
            Check(mols == 3 + 5 && vm.GrowLog.Contains("5 × Toluene"), $"grow with solvent: {mols} molecules · {vm.GrowLog.Split('\n').FirstOrDefault(l => l.Contains("Toluene"))}");
            vm.RemoveGrowSmall(vm.GrowSmall[0]);
        }

        // The pipeline in order (design/boards/PipelineGrow, ExportCenter): the force field chosen in Grow is assigned when
        // growing finishes, the strip says where the project is, and the Export center writes LAMMPS and GROMACS files
        {
            vm.UsePolystyreneInGrow();
            vm.GrowChainsD = 3;
            vm.GrowDpD = 6;
            vm.GrowDensityD = 0.3m;
            vm.Field.FfIndex = vm.Field.Library.ToList().FindIndex(e => e.Id == "gaff-amber25");
            vm.Field.ChargeMode = 0;
            vm.GrowAssignField = true;
            vm.Grow().GetAwaiter().GetResult();
            var steps = string.Join(" · ", vm.PipelineSteps.Select(s => $"{s.Name} {s.State}"));
            // the strip states what the active structure is and what was done to it, in no imposed order
            Check(vm.Field.Assigned && vm.Field.Complete && vm.PipelineSteps.Count == 6 && vm.PipelineSteps[1].Detail == "Polymer cell"
                  && vm.PipelineSteps[2].State == "done" && vm.PipelineSteps[2].Detail.StartsWith("GAFF") && vm.PipelineNextLabel == ""
                  && vm.ActiveItem?.Origin == "Polymer cell",
                  $"structure status after Grow: force field assigned · {steps}");
            // the project keeps every structure: another one opened, then the grownCell cell picked again with its force field
            {
                var grownCell = vm.ActiveItem!;
                var structuresBefore = vm.ProjectItems.Count;
                vm.Open(Path.Combine(dir, "ps_melt.data"));
                var opened = vm.ActiveItem!;
                var twoKept = vm.ProjectItems.Count == structuresBefore + 1 && opened != grownCell && !vm.Field.Assigned;
                vm.Activate(grownCell);
                Check(twoKept && vm.ActiveItem == grownCell && vm.Field.Assigned && vm.Field.ForceFieldName.StartsWith("GAFF") && vm.Document == grownCell.Doc,
                      $"project: {vm.ProjectItems.Count} structures · back on {grownCell.Name} with {vm.Field.ForceFieldName}");
                // Minimise keeping the original: the result is a new structure of the project, force field kept
                vm.ResultAsNew = true;
                vm.RelaxIterationsD = 30;
                vm.Relax().GetAwaiter().GetResult();
                vm.ResultAsNew = false;
                var result = vm.ActiveItem!;
                Check(result != grownCell && result.Name.Contains("relaxed") && !result.Name.Contains("minimised ·") && vm.Field.Assigned && vm.ProjectItems.Count == structuresBefore + 2
                      && result.History.Contains("minimised") && !grownCell.History.Contains("minimised"),
                      $"minimise into a new structure: {result.Name} · {result.History} · original: '{grownCell.History}'");
                vm.CloseDocument();   // the result leaves the project; another structure becomes active
                vm.Activate(opened);
                vm.CloseDocument();
                Check(vm.ProjectItems.Count == structuresBefore && vm.HasDocument && vm.ActiveItem == grownCell, $"close one structure: {vm.ProjectItems.Count} left, working on {vm.Title}");
            }
            vm.OpenExportCenter();
            vm.RefreshEnginesNow();
            Check(vm.IsExportCenter && vm.IsExportRail && !vm.EngineHasError && vm.EngineLammpsFiles.Count == 2 && vm.EngineGromacsFiles.Count >= 4
                  && vm.EnginePreviewLines.Count > 10 && vm.EngineChecks.Count >= 4, $"export center: {vm.EngineSummary} {vm.EngineError}");
            var pkg = Path.Combine(outDir, "caps-selftest-engines");
            if (Directory.Exists(pkg)) Directory.Delete(pkg, true);
            vm.EngineFolder = pkg;
            vm.WriteEngines().GetAwaiter().GetResult();
            var deck = File.Exists(Path.Combine(pkg, "system.in")) ? File.ReadAllText(Path.Combine(pkg, "system.in")) : "";
            Check(deck.Contains("pair_coeff") && deck.Contains("all npt") && File.Exists(Path.Combine(pkg, "system.itp")) && File.Exists(Path.Combine(pkg, "system_em.mdp")),
                  $"export center: wrote {Directory.GetFiles(pkg).Length} files · {vm.Status}");
            // incomplete force field: the Export center refuses and says where to go
            vm.Field.Clear().GetAwaiter().GetResult();
            vm.RefreshEnginesNow();
            Check(vm.EngineHasError && vm.EngineError.Contains("force field") && vm.PipelineSteps[2].State == "todo", $"export center without a force field: {vm.EngineError}");
            vm.SetModule(8);
        }

        // Add hydrogens by pH: a peptide built at pH 7, stripped of its hydrogens, gets them back at pH 7
        {
            var (pep, _) = CapsDocument.PeptideBuild("{\"sequence\":\"KDEHRYCA\",\"ph\":7,\"cleanup\":false}", "pH peptide");
            var h0 = pep.Summary().Atoms;
            pep.Select("{\"mode\":\"element\",\"pattern\":\"H\"}");
            pep.Edit("{\"op\":\"delete\",\"atoms\":\"selection\"}");
            var heavy = pep.Summary().Atoms;
            pep.SetPh(7);
            var plan = System.Text.Json.Nodes.JsonNode.Parse(pep.HydrogenPlan())!;
            pep.Edit("{\"op\":\"add_h\",\"atoms\":\"\"}");
            var h1 = pep.Summary().Atoms;
            Check(h1 == h0, $"add hydrogens at pH 7: {heavy} heavy → {h1} atoms (built {h0}) · plan net formal charge {(double?)plan["net_charge"]}");
            pep.Dispose();
        }

        // Model resolution › backmap: the melt's beads saved and read back, the atoms carried onto them and relaxed
        {
            vm.Open(Path.Combine(dir, "ps_melt.data"));
            vm.MrPerBead = 5;
            var (cg, _) = vm.Document!.ResolutionConvert("{\"to\":\"coarse-grained\",\"per_bead\":5}", "beads");
            var beadsFile = Path.Combine(outDir, "ps_beads.xyz");
            cg.Save(beadsFile);
            var nBeads = cg.Summary().Atoms;
            cg.Dispose();
            vm.BackmapFrom(beadsFile).GetAwaiter().GetResult();
            Check(vm.Document!.Summary().Atoms == 1300 && vm.Status.StartsWith("Backmapped") && vm.Title.Contains("backmapped"),
                  $"backmap: {nBeads} beads → {vm.Document!.Summary().Atoms} atoms · {vm.Status}");
        }

        // Equilibrate › chain statistics: an alkane cell gets the RIS polyethylene reference (dashed), polystyrene does not
        {
            var (pe, _) = CapsDocument.GrowChains("{\"units\":[{\"name\":\"E\",\"smiles\":\"*CC*\"}],\"dp\":20}",
                new CapsGrowOpts { Chains = 3, Dp = 0, Seed = 2, Density = 0.3, ContactScale = -1.0, Curve = 1 }, null, "pe");
            var peFile = Path.Combine(outDir, "pe_ris.data");
            pe.Save(peFile);
            pe.Dispose();
            vm.Open(peFile);
            var risOk = vm.RisCurve.Length > 5 && vm.ChainNote.Contains("RIS polyethylene");
            vm.Open(Path.Combine(dir, "ps_melt.data"));
            Check(risOk && vm.RisCurve.Length == 0, $"RIS reference: {vm.RisCurve.Length} on PS · alkane note ok {risOk}");
        }

        // Biomolecule › nucleic acid: a DNA strand of 6 nt (D-sugars kept, 3′ phosphate), then RNA
        {
            vm.OpenBio();
            vm.BioNucleic = true;
            vm.NaSequence = "ACGTTA";
            vm.BuildNucleic().GetAwaiter().GetResult();
            var dna = vm.Document!.Summary();
            int Oxygens() { var n = 0; for (var i = 0; i < vm.Document!.Summary().Atoms; ++i) n += vm.Document!.Atom(i).Element == 8 ? 1 : 0; return n; }
            var oDna = Oxygens();
            vm.NaRna = true;   // T becomes U
            var rnaSeq = vm.NaSequence;
            vm.BuildNucleic().GetAwaiter().GetResult();
            var rna = vm.Document!.Summary();
            var oRna = Oxygens();   // one 2′-OH per nucleotide
            Check(dna.Molecules == 1 && dna.Atoms > 150 && rnaSeq == "ACGUUA" && oRna == oDna + 6 && vm.Title.StartsWith("RNA"),
                  $"nucleic acid: DNA {dna.Atoms} atoms ({oDna} O), RNA {rnaSeq} {rna.Atoms} atoms ({oRna} O) · {vm.Status}");
            vm.BioNucleic = false;
            vm.NaRna = false;
        }

        // Close goes back to Start
        vm.SetModule(1);
        vm.CloseAllStructures();
        Check(vm.NoDocument && vm.IsStudio && vm.Title == "" && vm.ProjectItems.Count == 0, "close all structures: back to Start with no document");

        Console.WriteLine(fails == 0 ? "all checks passed" : $"{fails} check(s) failed");
        return fails == 0 ? 0 : 1;
    }
}
