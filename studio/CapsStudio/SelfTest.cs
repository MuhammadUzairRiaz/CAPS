using System.Collections.Generic;
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

        Check(Native.AbiVersion() == 59, "native ABI version 59");
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
        vm.QuickText = "toluene";
        var namedMol = vm.QuickKind == 5 && vm.QuickBadge.StartsWith("Toluene");
        vm.QuickText = "polystyrene";
        var namedPoly = vm.QuickKind == 5 && vm.QuickBadge.Contains("polymer library");
        vm.QuickText = "CCO";
        Check(namedMol && namedPoly && vm.QuickKind == 2, $"quick start names: toluene {namedMol} · polystyrene {namedPoly} · CCO as SMILES {vm.QuickKind == 2}");
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
        Check(vm.PickedRows.Count >= 4 && vm.NeighbourRows.Count == 4, $"inspector: {vm.PickedTitle} · {string.Join(", ", vm.PickedRows.Select(r => r.Key))}");

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
            // the mixing rule in place of the force field's own (Dynamics › Mixing), and back
            {
                var own = vm.Field.MixingRule;
                vm.Field.MixingIndex = 2;
                for (var k = 0; k < 400 && vm.Field.Working; ++k) { Avalonia.Threading.Dispatcher.UIThread.RunJobs(); Thread.Sleep(10); }
                var geo = vm.Field.MixingRule;
                var noted = vm.Field.Log.Contains("geometric") || vm.Document!.FieldReport().Contains("mixing rule geometric in place of");
                vm.Field.MixingIndex = 0;
                for (var k = 0; k < 400 && vm.Field.Working; ++k) { Avalonia.Threading.Dispatcher.UIThread.RunJobs(); Thread.Sleep(10); }
                Check(own.StartsWith("Lorentz") && geo.StartsWith("geometric") && noted && vm.Field.MixingRule == own && vm.Field.MixingIndex == 0,
                      $"mixing override: {own} → {geo} → {vm.Field.MixingRule}");
            }
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
        // Force fields by group: the melt's first five chains and the rest, each GAFF2; then OPLS 2005 for the rest, refused
        // for its 1-4 scaling until the first group's is taken; a crystal under Tersoff's silicon (a literature potential)
        if (gaff >= 0)
        {
            vm.Field.GroupMode = true;
            vm.Field.Groups.Clear();
            var ga = vm.Field.AddGroup("A", "1-5");
            ga.FfIndex = gaff; ga.Charges = 2;
            var gb = vm.Field.AddGroup("B", "rest");
            gb.FfIndex = gaff; gb.Charges = 2;
            vm.Field.AssignGroups().GetAwaiter().GetResult();
            var halves = vm.Field.IsGrouped && vm.Field.Complete;
            var oplsIx = vm.Field.Library.ToList().FindIndex(x => x.Id == "opls2005");
            var refused = false;
            var mixedName = "";
            if (oplsIx >= 0)
            {
                gb.FfIndex = oplsIx; gb.Charges = 0;
                vm.Field.Scaling14 = 2;   // refuse different 1-4 scalings (the default keeps each group's own)
                vm.Field.AssignGroups().GetAwaiter().GetResult();
                refused = vm.Field.Log.Contains("1-4", StringComparison.Ordinal);
                vm.Field.Scaling14 = 1;
                vm.Field.AssignGroups().GetAwaiter().GetResult();
                mixedName = vm.Field.Complete ? vm.Field.ForceFieldName : vm.Field.Log;
            }
            vm.Field.PairsText = "c3 bad";
            vm.Field.AssignGroups().GetAwaiter().GetResult();
            var badPairs = vm.Field.Log.StartsWith("Cross pair line", StringComparison.Ordinal);
            vm.Field.PairsText = "";
            Check(halves && refused && mixedName.Contains(" + ") && badPairs, $"Field by group: halves {halves} · OPLS refused {refused} · {mixedName} · bad pairs refused {badPairs}");
            // a silicon crystal (molecule 1, held) under Tersoff below methane typed by GAFF2
            var mbDir = Path.Combine(Path.GetTempPath(), "caps_selftest_manybody");
            Directory.CreateDirectory(mbDir);
            File.WriteAllText(Path.Combine(mbDir, "Si.tersoff"), "# UNITS: metal CITATION: Tersoff, Phys Rev B, 37, 6991 (1988)\n" +
                              "Si Si Si 3.0 1.0 1.3258 4.8381 2.0417 0.0000 22.956 0.33675 1.3258 95.373 3.0 0.2 3.2394 3264.7\n");
            var pdb = new System.Text.StringBuilder("CRYST1   30.000   30.000   30.000  90.00  90.00  90.00 P 1           1\n");
            (string Rn, int Res, double X, double Y, double Z, string El)[] siRows =
                [("SI", 1, 10, 10, 10, "Si"), ("SI", 1, 12.35, 10, 10, "Si"), ("MET", 2, 11, 10, 14, "C"), ("MET", 2, 11, 10, 15.09, "H"),
                 ("MET", 2, 12.03, 10, 13.64, "H"), ("MET", 2, 10.49, 10.89, 13.64, "H"), ("MET", 2, 10.49, 9.11, 13.64, "H")];
            for (var k = 0; k < siRows.Length; ++k)
                pdb.Append(string.Format(System.Globalization.CultureInfo.InvariantCulture, "HETATM{0,5} {1,-4} {2,3} A{3,4}    {4,8:F3}{5,8:F3}{6,8:F3}  1.00  0.00          {7,2}\n",
                                         k + 1, siRows[k].El, siRows[k].Rn, siRows[k].Res, siRows[k].X, siRows[k].Y, siRows[k].Z, siRows[k].El));
            pdb.Append("END\n");
            var pdbPath = Path.Combine(mbDir, "si_methane.pdb");
            File.WriteAllText(pdbPath, pdb.ToString());
            vm.Open(pdbPath);
            vm.Document!.SetHeldMolecule(1);
            vm.Field.SuggestGroups();
            // the held silicon: the library's Tersoff silicon suggested for it
            var suggested = vm.Field.Groups.Count == 2 && vm.Field.Groups[0].Molecules == "1" && vm.Field.Groups[1].Molecules == "rest" &&
                            vm.Field.Groups[0].IsPotential && vm.Field.Groups[0].Pick > 0 && vm.Field.PotentialLibrary[vm.Field.Groups[0].Pick].Id == "si-tersoff1988";
            if (vm.Field.Groups.Count == 2)
            {
                vm.Field.Groups[1].FfIndex = gaff;
                vm.Field.Groups[1].Charges = 2;
                vm.Field.AssignGroups().GetAwaiter().GetResult();
            }
            var notes = string.Join(" ", vm.Field.Notes);
            Check(suggested && vm.Field.Complete && vm.Field.IsGrouped && notes.Contains("tersoff", StringComparison.Ordinal) && vm.Field.Swatches.Any(x => x.Name == "Si"),
                  $"Field by group with Tersoff silicon: suggested {suggested} (library {vm.Field.PotentialLibrary.Count}, pick {vm.Field.Groups.FirstOrDefault()?.Pick}, kind {vm.Field.Groups.FirstOrDefault()?.Kind}) · {vm.Field.ForceFieldName} · {vm.Field.Log}");
            // LAMMPS agrees term by term once the Tersoff energy (which CAPS does not compute) is taken out
            vm.CompareEnergies().GetAwaiter().GetResult();
            Check(vm.ParityOk || vm.ParityDetail.Contains("not found"), $"compare energies beside Tersoff: {vm.ParityText}\n{vm.ParityTable}");
            // a carbon filler: the library's AIREBO suggested; LAMMPS reads it in metal units only, so the deck is in metal units
            var cc = new System.Text.StringBuilder("CRYST1   30.000   30.000   30.000  90.00  90.00  90.00 P 1           1\n");
            (string Rn, int Res, double X, double Y, double Z, string El)[] ccRows =
                [("CC", 1, 10, 10, 10, "C"), ("CC", 1, 11.42, 10, 10, "C"), ("MET", 2, 11, 10, 14, "C"), ("MET", 2, 11, 10, 15.09, "H"),
                 ("MET", 2, 12.03, 10, 13.64, "H"), ("MET", 2, 10.49, 10.89, 13.64, "H"), ("MET", 2, 10.49, 9.11, 13.64, "H")];
            for (var k = 0; k < ccRows.Length; ++k)
                cc.Append(string.Format(System.Globalization.CultureInfo.InvariantCulture, "HETATM{0,5} {1,-4} {2,3} A{3,4}    {4,8:F3}{5,8:F3}{6,8:F3}  1.00  0.00          {7,2}\n",
                                        k + 1, ccRows[k].El, ccRows[k].Rn, ccRows[k].Res, ccRows[k].X, ccRows[k].Y, ccRows[k].Z, ccRows[k].El));
            cc.Append("END\n");
            var ccPath = Path.Combine(mbDir, "cc_methane.pdb");
            File.WriteAllText(ccPath, cc.ToString());
            // coordination: methane's carbon read as tetrahedral, made square planar exactly, back with undo; out-of-plane angle
            vm.Open(pdbPath);
            vm.Pick(2);
            var coordRow = vm.PickedRows.FirstOrDefault(r => r.Key == "Coordination")?.Value ?? "";
            vm.Pick(3, true); vm.Pick(4, true); vm.Pick(5, true);
            var oop = vm.MeasureText.Contains("Out of plane", StringComparison.Ordinal);
            vm.Pick(2);
            vm.CoordinationIndex = 3;
            vm.ApplyCoordination();
            var planar = vm.PickedRows.FirstOrDefault(r => r.Key == "Coordination")?.Value ?? "";
            Check(coordRow.Contains("4 neighbours · tetrahedral", StringComparison.Ordinal) && oop && planar.Contains("square planar (RMS 0.0°)", StringComparison.Ordinal),
                  $"coordination: {coordRow} → {planar} · out-of-plane shown {oop}");
            vm.Field.Clear().GetAwaiter().GetResult();
            vm.Open(ccPath);
            vm.Document!.SetHeldMolecule(1);
            vm.Field.SuggestGroups();
            var airebo = vm.Field.Groups.Count == 2 && vm.Field.Groups[0].IsPotential && vm.Field.PotentialLibrary[vm.Field.Groups[0].Pick].Style == "airebo";
            if (vm.Field.Groups.Count == 2) { vm.Field.Groups[1].FfIndex = gaff; vm.Field.Groups[1].Charges = 2; }
            vm.Field.AssignGroups().GetAwaiter().GetResult();
            var metalDeck = vm.Field.Complete && vm.Document!.LammpsInput("system.data").Contains("units           metal", StringComparison.Ordinal);
            vm.CompareEnergies().GetAwaiter().GetResult();
            Check(airebo && metalDeck && (vm.ParityOk || vm.ParityDetail.Contains("not found")),
                  $"carbon filler under AIREBO: suggested {airebo} · metal units {metalDeck} · {vm.ParityText}\n{vm.ParityTable}");
            // MEAM from the library: the SiC set's entries in its order, carbon mapped
            var meamIx = vm.Field.PotentialLibrary.FindIndex(p => p.Id == "sic-meam");
            var meamDeck = "";
            if (meamIx > 0 && vm.Field.Groups.Count == 2)
            {
                vm.Field.Groups[0].Pick = meamIx;
                vm.Field.AssignGroups().GetAwaiter().GetResult();
                meamDeck = vm.Field.Complete ? vm.Document!.LammpsInput("system.data") : vm.Field.Log;
            }
            Check(meamDeck.Contains("* * meam library.meam Si C SiC.meam C NULL NULL", StringComparison.Ordinal) && vm.Field.Groups[0].IsMeam && vm.Field.Groups[0].EntriesText == "Si=Si C=C",
                  "MEAM from the library: " + meamDeck.Split('\n').FirstOrDefault(l => l.Contains("meam", StringComparison.Ordinal) && l.StartsWith("pair_coeff", StringComparison.Ordinal)));
            vm.Field.Clear().GetAwaiter().GetResult();
            vm.Field.GroupMode = false;
            vm.Field.Groups.Clear();
            vm.Open(Path.Combine(dir, "ps_melt.lammpstrj"), Path.Combine(dir, "ps_melt.data"));
        }
        // Type by hand: polystyrene's head, body and tail typed by OPLS-AA 2024's rules, one body CH changed by hand
        // (with every equivalent atom), applied to the whole melt by environment; an SBR copolymer example holds every junction
        {
            var o24 = vm.Field.Library.ToList().FindIndex(x => x.Id == "oplsaa2024-moltemplate");
            vm.Field.FfIndex = o24;
            vm.UsePolystyreneInGrow();
            while (vm.PolyUnits.Count > 1) vm.RemovePolyUnit(vm.PolyUnits[^1]);
            if (vm.PolyUnits.Count == 0) vm.AddPolyUnit("Styrene", "*CC(*)c1ccccc1");
            else { vm.PolyUnits[0].Name = "Styrene"; vm.PolyUnits[0].Smiles = "*CC(*)c1ccccc1"; }
            vm.OpenUnitTyping();
            var exampleAtoms = vm.UtAtoms.Count;
            var roles = vm.UtAtoms.Select(a => a.RoleKind).Distinct().Count();
            var bodyCh = vm.UtAtoms.FirstOrDefault(a => a.RoleKind == 1 && a.Auto.StartsWith("515_", StringComparison.Ordinal));
            var applied = "";
            if (bodyCh != null)
            {
                vm.UtSelected = bodyCh;
                vm.UtTypeSelected = vm.UtTypes.FirstOrDefault(t => t.Name.StartsWith("137_", StringComparison.Ordinal));
                vm.AssignUnitType();
                vm.ApplyUnitTyping();
                applied = vm.UtStatus;
            }
            var rep = vm.Document!.FieldReport();
            var count137 = System.Text.RegularExpressions.Regex.Matches(rep, "\"type\": ?\"137_").Count;
            Check(exampleAtoms == 50 && roles == 3 && bodyCh != null && vm.UtUntyped == 0 && count137 == 60 && applied.Contains("atoms of the structure typed from the example"),
                  $"type by hand: {exampleAtoms} example atoms, {roles} roles, 137 on {count137} melt atoms · {applied} {vm.UtStatus}");
            vm.Field.Clear().GetAwaiter().GetResult();
            vm.SetModule(8);
            Check(MainViewModel.DeBruijnPairs("AB") == "AABBA" && MainViewModel.DeBruijnPairs("ABC").Length == 10, $"junction sequence: {MainViewModel.DeBruijnPairs("AB")} · {MainViewModel.DeBruijnPairs("ABC")}");
        }

        // Automatic charges: OPLS 2005's own        // Automatic charges: OPLS 2005's own — charge keys and its bond charge increments (no Gasteiger stand-in)
        var opls = vm.Field.Library.ToList().FindIndex(x => x.Id == "opls2005");
        if (opls >= 0)
        {
            vm.Field.FfIndex = opls;
            vm.Field.ChargeMode = 0;
            vm.Field.Assign().GetAwaiter().GetResult();
            var rep = vm.Document!.FieldReport();
            Check(vm.Field.Assigned && vm.Field.Complete && rep.Contains("charges from bond increments", StringComparison.Ordinal) &&
                  !rep.Contains("Gasteiger", StringComparison.Ordinal) && rep.Contains("\"ck\"", StringComparison.Ordinal),
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
            Check(complete24 && vm.Field.Alternatives.Any(a => a.Id == "pcff-frc") && !vm.Field.HasUntypedGroups && !vm.Field.HasBalanceNote,
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
            // both in the split view: frame 1 on the right, frame 3 on the left
            var st = vm.Document.StateDocument("{\"kind\":\"frame\",\"index\":0}", "frame 1");
            var right10 = st.Atom(10);
            st.Dispose();
            Check(Math.Abs(right10.X - f0.X) < 1e-9, $"state document: frame 1's atom 11 at x {right10.X:F3}");
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
                                         FramePs = 0.5, ThermoPs = 0.1, UntilConverged = 1, BlockPs = 0.3, MaxBlocks = 3, TolDensity = 0.5, TolEnergy = 5, TolRg = 0.5, TolInternal = 0.5 };
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
        {
            // an internal-distance target far from the cell's curve: the blocks never converge, and the check names the curve
            var far = Enumerable.Range(0, 30).Select(k => k == 0 ? 0.0 : 50.0).ToArray();
            var h = System.Runtime.InteropServices.GCHandle.Alloc(far, System.Runtime.InteropServices.GCHandleType.Pinned);
            try
            {
                var (farConv, farRep) = cell.Equilibrate("nvt 0.3 ps T 300 # hold", eqOpts with { MaxBlocks = 2, InternalTarget = h.AddrOfPinnedObject(), InternalTargetN = far.Length }, null);
                Check(!farConv && farRep.Contains("of the target curve"), $"equilibrate: a far internal-distance target keeps it unconverged ({farConv})");
            }
            finally { h.Free(); }
        }
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
        // PBSA from the library; the composition calculator: BS:BA 80:20 mol at DP 25 → 20 and 5 per chain; 80:20 by weight → mole
        if (vm.PolymerLibrary.FirstOrDefault(e => e.Id == "C112") is { } pbsa && vm.PolymerLibrary.Any(e => e.Id == "C113" && e.Copolymer))
        {
            if (File.Exists(MainViewModel.UserPolymerFile)) File.Delete(MainViewModel.UserPolymerFile);
            vm.UseLibrary(pbsa, null);
            var dpBefore = vm.GrowDpD;
            vm.GrowDpD = 12;
            var strip12 = vm.PolyStripUnits.Length;
            vm.GrowDpD = 25;
            // seeds: a fixed seed draws the same sequence every time; new each run draws others, shows the seed, Keep repeats it
            {
                vm.PolySeed.Fresh = false; vm.PolySeed.Value = 7;
                var a1 = string.Concat(vm.PolyStripUnits); vm.PolyRedraw(); var a2 = string.Concat(vm.PolyStripUnits);
                vm.PolySeed.Value = 8; var b1 = string.Concat(vm.PolyStripUnits);
                vm.PolySeed.Fresh = true;
                var draws = new HashSet<string>();
                for (var k = 0; k < 6; ++k) { vm.PolyRedraw(); draws.Add(string.Concat(vm.PolyStripUnits)); }
                var shown = string.Concat(vm.PolyStripUnits);
                var used = vm.PolySeed.Used;
                vm.PolySeed.Keep();
                var kept = string.Concat(vm.PolyStripUnits);
                Check(a1 == a2 && a1 != b1 && draws.Count >= 4 && used > 0 && !vm.PolySeed.Fresh && vm.PolySeed.Value == used && kept == shown,
                      $"polymer seed: fixed repeats ({a1 == a2}), another seed differs ({a1 != b1}), new each run gave {draws.Count}/6 sequences, kept {used} reproduces ({kept == shown})");
                vm.PolySeed.Value = 1;
                // Pack: the seed line follows the picker; new each run writes a fresh seed per pack
                vm.NewPackInput();
                vm.PackSeedChoice.Value = 42;
                Check(vm.PackText.Contains("seed 42\n"), $"pack seed line: {vm.PackText.Split('\n').FirstOrDefault(l => l.StartsWith("seed"))}");
                vm.PackSeedChoice.Value = 1;
            }
            vm.GrowDpD = null;   // the box emptied while typing: the DP stays
            // PBSA ends: the head * is on O, the tail * on the acid carbonyl C — hydrogen / hydroxyl gives HO–…–COOH
            vm.HeadCap = 0; vm.TailCap = Array.IndexOf(MainViewModel.EndGroups, "hydroxyl");
            var endsOk = vm.PolyEndsText.Contains("Head: hydrogen on O → –OH") && vm.PolyEndsText.Contains("Tail: hydroxyl on C=O → acid –COOH") && !vm.PolyEndsText.Contains("⚠");
            vm.TailCap = Array.IndexOf(MainViewModel.EndGroups, "carboxyl");
            vm.HeadCap = Array.IndexOf(MainViewModel.EndGroups, "hydroxyl");
            var warnOk = vm.PolyEndsText.Split('\n').All(l => l.Contains("⚠"));
            vm.HeadCap = 0; vm.TailCap = 0;
            Check(endsOk && warnOk, $"end groups explained: {vm.PolyEndsText.Replace('\n', '|')}");
            Check(strip12 == 12 && vm.PolyStripUnits.Length == 25 && vm.PolyPreview.Contains("A ") , $"sequence follows the DP at once: {strip12} then {vm.PolyStripUnits.Length} units · {vm.PolyPreview}");
            vm.CompBasis = 0;
            vm.PolyUnits[0].Target = 80; vm.PolyUnits[1].Target = 20;
            Check(vm.CompShown && vm.CompRows.Count == 2 && vm.CompRows[0].PerChain.StartsWith("20 · ") && vm.CompRows[1].PerChain.StartsWith("5 · "),
                $"composition 80:20 mol, DP 25: {string.Join(" | ", vm.CompRows.Select(r => r.PerChain))}");
            // by weight: x_BS = (80/172.18) / (80/172.18 + 20/200.23) = 0.8231 (BS C8H12O4, BA C10H16O4 repeat units)
            vm.CompBasis = 1;
            vm.ApplyComposition(true);
            var mBS = vm.PolyUnits[0].Mass; var mBA = vm.PolyUnits[1].Mass;
            var xBS = 80 / mBS / (80 / mBS + 20 / mBA);
            Check(Math.Abs(mBS - 172.18) < 0.05 && Math.Abs(mBA - 200.23) < 0.05 && Math.Abs((double)vm.PolyUnits[0].Weight - xBS) < 1e-4 && vm.PolySequence == 6
                  && vm.PolySpecJson().Contains("\"shuffled\"") && vm.CompRows[0].Wt == "80.0 %",
                $"composition 80:20 by weight → mole {vm.PolyUnits[0].Weight} (expect {xBS:0.0000}), masses {mBS:0.00}/{mBA:0.00}, {vm.CompRows[0].Wt}");
            // saved to your polymers, found in the library, then removed
            vm.PolyName = "My PBSA 80-20 wt";
            var savedPoly = vm.SavePolymerToLibrary();
            var mine = vm.PolymerLibrary.FirstOrDefault(e => e.User && e.Name == "My PBSA 80-20 wt");
            Check(mine != null && mine.Copolymer && File.Exists(MainViewModel.UserPolymerFile) && !MainViewModel.UserPolymerFile.StartsWith(AppSettings.Folder), $"saved to your polymers: {savedPoly}");
            if (mine != null)
            {
                vm.UseLibrary(mine, null);
                Check(vm.PolyUnits.Count == 2 && vm.PolySequence == 6 && Math.Abs((double)vm.PolyUnits[0].Weight - xBS) < 1e-4, $"your PBSA reopened: {vm.PolyPreview}");
                var removed = vm.RemoveUserPolymer(mine);
                Check(!vm.PolymerLibrary.Any(e => e.User), $"your polymer removed: {removed}");
            }
            vm.GrowDpD = dpBefore;
            vm.PolyName = "";
        }
        else Check(false, "PBSA / PBAT in the library");
        if (enr != null)
        {
            vm.UseLibrary(enr, null);
            Check(!vm.PolyHasError && vm.PolyUnits.Count == 2 && vm.PolyStripUnits.Length == vm.GrowDpD, $"ENR-50 chain: {vm.PolyPreview} {vm.PolyError}");
            // end groups: a tert-butyl head (an initiator fragment) and a hydroxyl tail travel with the chain spec
            vm.HeadCap = 3; vm.TailCap = 6;
            var capped = vm.PolySpecJson();
            Check(capped.Contains("\"head_cap\":\"tert-butyl\"") && capped.Contains("\"tail_cap\":\"hydroxyl\"") && !vm.PolyHasError, $"end groups: {capped} {vm.PolyError}");
            vm.HeadCap = 0; vm.TailCap = 0;
            vm.PolyLinkage = 1;
            Check(vm.PolySpecJson().Contains("\"linkage\":\"head-to-head\"") && vm.PolyPreview.Contains(" reversed") && !vm.PolyHasError, $"head-to-head linkage: {vm.PolyPreview} {vm.PolyError}");
            vm.PolyLinkage = 0;
            // your own ENR (the library's units, alternating): guessed a rubber, saved tagged, found at once under Rubbers + search
            vm.PolySequence = 1;
            vm.PolyName = "ENR-50 test alternate";
            vm.LibraryRubberOnly = true;
            vm.LibraryQuery = "ENR";
            var guessed = vm.PolySaveRubber;
            var savedEnr = vm.SavePolymerToLibrary();
            var mineEnr = vm.LibraryShown.FirstOrDefault(e => e.User && e.Name == "ENR-50 test alternate");
            Check(guessed && mineEnr is { Rubber: true } && File.ReadAllText(MainViewModel.UserPolymerFile).Contains("\"rubber\""), $"your ENR saved as a rubber and shown under Rubbers + ENR: {savedEnr}");
            // polystyrene is no rubber: saved untagged it is shown by turning the Rubbers filter off
            if (vm.PolymerLibrary.FirstOrDefault(e => !e.Copolymer && e.Smiles.Contains("c1ccccc1")) is { } psLib)
            {
                vm.PolyUnits.Clear(); vm.AddPolyUnit(psLib.Name, psLib.Smiles);
                vm.PolyName = "PS test";
                var psRubber = vm.PolySaveRubber;
                vm.SavePolymerToLibrary();
                Check(!psRubber && !vm.LibraryRubberOnly && vm.LibraryShown.Any(e => e.User && e.Name == "PS test" && !e.Rubber), "your PS saved, not a rubber, and shown");
                if (vm.PolymerLibrary.FirstOrDefault(e => e.User && e.Name == "PS test") is { } m2) vm.RemoveUserPolymer(m2);
            }
            if (mineEnr != null) vm.RemoveUserPolymer(mineEnr);
            vm.LibraryRubberOnly = false; vm.LibraryQuery = ""; vm.PolyName = "";
            vm.UseLibrary(enr, null);
            vm.GrowChainsD = 4;
            vm.GrowDensityD = 0.3m;
            vm.SendPolymerToGrow();
            vm.Grow().GetAwaiter().GetResult();
            Check(vm.Document != null && vm.Document.Summary().Molecules == 4 && vm.GrowComponentName.StartsWith("ENR"), $"grown ENR-50 cell: {vm.Status}");
            // the grown cell shows as the live view did (wrapped by default), and the switch changes it both ways
            {
                var wrappedAfterGrow = vm.Wrap == vm.GrowWrap;
                vm.GrowUnwrapped = true;
                var unwrapped = !vm.Wrap && !vm.GrowWrap;
                vm.GrowWrap = true;
                Check(wrappedAfterGrow && unwrapped && vm.Wrap, $"grow view wrap: after growing {wrappedAfterGrow}, unwrapped {unwrapped}, back {vm.Wrap}");
                vm.Wrap = false;
            }
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
            // ENR crosslinked through a diacid (succinic acid, the PBS acid end's model): each COOH opens an epoxide at its
            // tertiary carbon, so a diacid reaching two chains is a link between them; links only between chains, to two links
            vm.RxSet = 5;
            vm.RxInsertSmiles = "OC(=O)CCC(=O)O";
            vm.RxInsertCount = 8;
            vm.InsertCurative().GetAwaiter().GetResult();
            vm.RxBetweenChains = true;
            vm.RxAutoCapture = true;
            vm.RxTargetKind = 1;
            vm.RxTargetValueD = 2;
            vm.RxCyclesD = 30;
            vm.RunReact().GetAwaiter().GetResult();
            var linksLine = vm.RxNetworkText.Split('\n')[0];
            Check(linksLine.StartsWith("Links between chains: ") && !linksLine.StartsWith("Links between chains: 0") && vm.RxNetworkText.Contains("Force field during the run"),
                  $"ENR + diacid: {linksLine} · {vm.RxLog.Split('\n')[0]}");
            // several reactions with weights: the reaction list follows the text
            vm.RxSet = 6;
            var several = vm.RxReactions.Count == 2 && vm.RxSeveral;
            vm.RxByWeights = true;
            vm.RxReactions[0].WeightD = 3;
            var weighted = vm.RxReactions.Select(r => r.Name).SequenceEqual(["anhydride_alcohol", "enr_acid_ester"]);
            Check(several && weighted, $"ENR + MAH: reactions {string.Join(", ", vm.RxReactions.Select(r => r.Name + " ×" + r.WeightD))}");
            vm.RxByWeights = false;
            // the reaction library: the cure templates and the worked schemes by category; a three-molecule scheme (ENR–MA–ENR)
            // loads as its steps
            vm.RxLibOpen = true;
            var cats = vm.RxLibCategories.ToList();
            var enrCat = cats.FindIndex(c => c.StartsWith("ENR"));
            if (enrCat >= 0) vm.RxLibCategory = enrCat;
            vm.RxLibSelected = vm.RxLibItems.FirstOrDefault(x => x.Id == "enr_mah_enr");
            vm.UseRxLib(false);
            var steps = vm.RxText.Contains("reaction epoxy_acid");
            vm.RxLibSelected = vm.RxLibItems.FirstOrDefault(x => x.Id == "enr_hydrolysis");
            vm.UseRxLib(true);
            Check(cats.Count >= 10 && steps && vm.RxReactions.Count == 2, $"reaction library: {cats.Count} categories · ENR–MA–ENR as steps {steps} · reactions {string.Join(", ", vm.RxReactions.Select(r => r.Name))}");
            vm.RxLibOpen = false;
            vm.RxTargetKind = 0;
            vm.RxInsertSmiles = "SS";
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
            // a generation-2 dendrimer of the same unit: 3 core arms, 18 branches of 2 units
            vm.GrowDpD = 3;
            vm.PolyArch = 4;
            vm.PolyArms = 3;
            vm.PolyArmDp = 2;
            vm.PolyGenerations = 2;
            var dest = System.Text.RegularExpressions.Regex.Match(vm.PolyPreview, @"per molecule: 18 branches on 3 core arms · ([\d,]+) atoms");
            vm.BuildPolyPreview().GetAwaiter().GetResult();
            var dsm = vm.PolyDoc?.Summary();
            Check(dest.Success && dsm != null && dsm.Value.Molecules == 1 && dsm.Value.Atoms == int.Parse(dest.Groups[1].Value.Replace(",", "")),
                  $"dendrimer NR: {vm.PolyPreview.Replace('\n', ' ')} · built {dsm?.Atoms} atoms in {dsm?.Molecules} molecule {vm.PolyError}");
            vm.PolyArch = 0;
            vm.PolyArmDp = 5;
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
        {   // Type by hand shows the film's polymer (from the structure's history), not the Polymer builder's
            var own = vm.PolymerOfDocument();
            var film = vm.SurfFilmItem?.Name ?? "";
            Check(own != null && own.Contains("units") && film.Length > 0 && own.Contains(System.Text.Json.Nodes.JsonNode.Parse(vm.SurfFilmItem!.Spec)!["units"]![0]!["smiles"]!.GetValue<string>().Replace("\\", "\\\\")),
                  $"the open structure's polymer: {film} · {own}");
        }
        {   // typing by hand inside a by-group assignment: the film's group typed from the example in OPLS-AA 2024, the
            // surface keeps its own (UFF here), cross terms remade; the LAMMPS input lists each group's types with masses
            var ffKeep = vm.Field.FfIndex;
            var o24g = vm.Field.Library.ToList().FindIndex(x => x.Id == "oplsaa2024-moltemplate");
            var uffg = vm.Field.Library.ToList().FindIndex(x => x.Id == "uff");
            vm.Field.Groups.Clear();
            var gSurf = vm.Field.AddGroup("surface", "1");
            gSurf.FfIndex = uffg;
            var gm = vm.Field.AddGroup("film", "rest");
            gm.FfIndex = o24g;
            vm.Field.FfIndex = o24g;
            vm.OpenUnitTyping();
            var exUntyped = vm.UtUntyped;
            if (exUntyped == 0) vm.ApplyUnitTyping();
            for (var k = 0; k < 400 && vm.Field.Working; ++k) { Avalonia.Threading.Dispatcher.UIThread.RunJobs(); Thread.Sleep(25); }
            Avalonia.Threading.Dispatcher.UIThread.RunJobs();
            var deck = vm.Field.Assigned ? vm.Document!.LammpsInput("system.data") : "";
            Check(exUntyped == 0 && vm.Field.Assigned && vm.UtStatus.Contains("film") && vm.UtStatus.Contains("typed from the example") &&
                  deck.Contains("group           surface") && deck.Contains("group           film") && deck.Contains("# atom types by group") && deck.Contains("mass            1 "),
                  $"hand typing by group: example untyped {exUntyped} · {vm.UtStatus} · {vm.Field.ForceFieldName}");
            vm.Field.Clear().GetAwaiter().GetResult();
            vm.Field.Groups.Clear();
            vm.Field.FfIndex = ffKeep;
            vm.SetModule(8);
        }
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
            // the surface chosen by hand: the held molecule, then ids with no atoms (the page says why nothing splits)
            vm.IfUseHeld();
            var heldId = vm.Document!.HeldMolecule();
            var heldOk = vm.Analyze.SurfaceMolecules == heldId.ToString() && vm.IfSetupText.StartsWith($"surface: molecule {heldId}");
            vm.Analyze.SurfaceMolecules = "999";
            vm.Analyze.AxisIndex = 0;
            vm.RunInterface().GetAwaiter().GetResult();
            var noSurface = vm.IfNote.Contains("no atom in the surface") && Row("Surface top") == "—" && vm.IfSetupText == "surface: molecule 999 · along x";
            vm.Analyze.SurfaceMolecules = "1";
            vm.Analyze.AxisIndex = 2;
            Check(heldOk && noSurface, $"interface setup: held molecule {heldId} → '{vm.Analyze.SurfaceMolecules}'; ids with no atoms: {vm.IfNote}");
        }
        vm.SetModule(8);

        // Layer stack: a small grown polystyrene cell from CAPS Grow on the quartz slab, then a second slab on top
        {
            vm.UsePolystyreneInGrow();
            vm.GrowChainsD = 2; vm.GrowDpD = 5; vm.GrowDensityD = 0.4m;
            vm.Grow().GetAwaiter().GetResult();
            var cellAtoms = vm.Document?.Summary().Atoms ?? 0;
            vm.OpenSurface();
            vm.AddLayerFromOpen();
            vm.AddSlabLayer();
            var chip = vm.SurfStackChip;
            var matchRows = string.Join(" | ", vm.SurfMatchRows.Select(r => $"{r.Layer} {r.Supercell} {r.StrainA} {r.StrainB}"));
            Check(vm.SurfMatchRows.Count == 3 && vm.SurfMatchRows[0].StrainB == "0.00 %", "lattice matching rows: " + matchRows);
            vm.BuildSurface().GetAwaiter().GetResult();
            var st = vm.Document?.Summary();
            Check(st is { } ss && ss.Atoms > cellAtoms && ss.Molecules >= 2 + 2 && vm.Document!.Provenance().Contains("build.stack") && chip.Contains("3 layers") && !vm.SurfFilm,
                  $"layer stack: {chip} · {st?.Atoms} atoms, {st?.Molecules} molecules · {vm.SurfLog.Split('\n').FirstOrDefault()} {vm.SurfError}");
            while (vm.SurfHasExtra) vm.RemoveSurfLayer(0);
            vm.SurfFilm = true;
            vm.SetModule(8);
        }

        // Jobs: a Dynamics run paused (it holds its step), a second run queued behind it with its settings, resumed,
        // stopped, and the queued run starting by itself when the first ends
        {
            bool Until(Func<bool> c, int ms)
            {
                var sw = System.Diagnostics.Stopwatch.StartNew();
                while (!c() && sw.ElapsedMilliseconds < ms) { Avalonia.Threading.Dispatcher.UIThread.RunJobs(); Thread.Sleep(20); }
                return c();
            }
            vm.UsePolystyreneInGrow();
            vm.GrowChainsD = 1; vm.GrowDpD = 3; vm.GrowDensityD = 0.3m;
            vm.GrowAssignField = true;
            vm.Grow().GetAwaiter().GetResult();
            vm.SetModule(3);
            var steps = vm.MdStepsD;
            vm.MdStepsD = 5_000_000;
            var run = vm.RunMd();
            var started = Until(() => vm.MdRunning && vm.MdLog.StartsWith("step"), 30000);
            vm.PauseRun();
            Until(() => false, 500);
            var held = vm.MdLog;
            Until(() => false, 800);
            var still = vm.MdLog == held && vm.RunPaused && vm.Jobs.FirstOrDefault(j => j.IsRunning)?.StatusText == "paused";
            vm.MdStepsD = 3_000_000;
            vm.QueueMd();
            var queued = vm.QueuedCount == 1 && vm.Jobs[0].IsQueued && vm.JobsSummary.Contains("1 queued");
            vm.ResumeRun();
            var moved = Until(() => vm.MdLog != held, 10000);
            vm.CancelMd();
            Until(() => run.IsCompleted, 30000);
            var second = Until(() => vm.MdRunning && vm.QueuedCount == 0 && vm.MdLog.Contains("of 3,000,000"), 30000);
            vm.CancelMd();
            Until(() => !vm.MdRunning, 30000);
            Check(started && still && queued && moved && second && !vm.RunPaused,
                  $"jobs: started {started} · paused and held {still} · queued {queued} · resumed {moved} · queued run started {second} · {vm.JobsSummary}");
            vm.MdStepsD = steps;

            // a remote job, end to end against stand-ins on this machine: ssh runs the command in a shell, scp copies, the
            // scheduler is none (the job runs in the background) and caps is the command-line tool of this build
            var cli = Path.GetFullPath(Path.Combine(dir, "..", "build", "cli", "caps"));
            if (!OperatingSystem.IsWindows() && File.Exists(cli))
            {
                var shim = Path.Combine(outDir, "remote-shim");
                var work = Path.Combine(outDir, "remote-work");
                Directory.CreateDirectory(shim);
                if (Directory.Exists(work)) Directory.Delete(work, true);
                File.WriteAllText(Path.Combine(shim, "ssh"), "#!/bin/bash\nwhile [[ \"$1\" == -* ]]; do case \"$1\" in -o|-p) shift 2;; *) shift;; esac; done\nshift\nexec bash -c \"$*\"\n");
                File.WriteAllText(Path.Combine(shim, "scp"), "#!/bin/bash\nargs=()\nwhile [ $# -gt 0 ]; do case \"$1\" in -P) shift 2;; -*) shift;; *) args+=(\"$1\"); shift;; esac; done\n" +
                                                            "dst=\"${args[${#args[@]}-1]}\"; dst=\"${dst#*:}\"\nfor a in \"${args[@]:0:${#args[@]}-1}\"; do src=\"${a#*:}\"; cp -R $src \"$dst\" || exit 1; done\n");
                File.WriteAllText(Path.Combine(shim, "caps"), $"#!/bin/bash\nexec \"{cli}\" \"$@\"\n");
                foreach (var f in new[] { "ssh", "scp", "caps" }) File.SetUnixFileMode(Path.Combine(shim, f), (UnixFileMode)0b111_101_101);
                var path = Environment.GetEnvironmentVariable("PATH") ?? "";
                Environment.SetEnvironmentVariable("PATH", shim + ":" + path);
                vm.AddHost();
                vm.HostName = "stand-in";
                vm.HostHostname = "localhost";
                vm.HostScheduler = "none";
                vm.HostWorkDir = work;
                vm.RunWhereIndex = vm.RunWhereChoices.Count - 1;
                vm.MdStepsD = 200;
                vm.SubmitRemote("Dynamics").GetAwaiter().GetResult();
                var rj = vm.Jobs.FirstOrDefault(j => j.IsRemote);
                var sent = rj is { IsRunning: true } && rj.Where.Contains("stand-in");
                var back = Until(() => { if (rj is { IsRunning: true }) vm.CheckRemote(rj).GetAwaiter().GetResult(); return rj is { IsRunning: false }; }, 120000);
                var result = rj?.Remote is { } rr ? Path.Combine(rr.Local, "out", rr.Stem + ".data") : "";
                Check(sent && back && rj!.IsDone && File.Exists(result) && File.Exists(result + ".provenance.json") && rj.CanOpenRemote,
                      $"remote job: sent {sent} · {rj?.Status} · {rj?.Where} · {string.Join(" | ", rj?.Log.Select(l => l.Text) ?? [])}");
                // the same for Relax (the open structure) and Grow (a new cell from its recipe)
                var others = new List<string>();
                foreach (var kind in new[] { "Relax", "Grow" })
                {
                    vm.GrowChainsD = 1; vm.GrowDpD = 3; vm.GrowDensityD = 0.3m;
                    vm.SubmitRemote(kind).GetAwaiter().GetResult();
                    var job = vm.Jobs.FirstOrDefault(j => j.IsRemote && j.Kind == kind);
                    Until(() => { if (job is { IsRunning: true }) vm.CheckRemote(job).GetAwaiter().GetResult(); return job is { IsRunning: false }; }, 120000);
                    others.Add($"{kind} {job?.Status}");
                }
                Check(others.All(o => o.EndsWith("done")), "remote relax and grow: " + string.Join(" · ", others));
                vm.RunWhereIndex = 0;
                vm.RemoveHost();
                Environment.SetEnvironmentVariable("PATH", path);
                vm.MdStepsD = steps;
            }

            // Copy as Python: the Relax and Dynamics scripts run with python3 on the saved structure; the React one compiles
            {
                var pySaved = Path.Combine(outDir, "caps-selftest-py.data");
                vm.Document!.Save(pySaved);
                vm.MdStepsD = 50;
                var scripts = new[] { ("relax", vm.RelaxPython()), ("md", vm.MdPython()), ("react", vm.ReactPython()), ("pack", vm.PackPython()) };
                var results = new List<string>();
                foreach (var (name, text) in scripts)
                {
                    var file = Path.Combine(outDir, $"caps-selftest-{name}.py");
                    var body = System.Text.RegularExpressions.Regex.Replace(text, @"caps\.open\([^\n]*\)", $"caps.open(\"{pySaved}\")");
                    File.WriteAllText(file, body);
                    var psi = new System.Diagnostics.ProcessStartInfo("python3") { RedirectStandardOutput = true, RedirectStandardError = true, UseShellExecute = false, WorkingDirectory = outDir };
                    if (name is "react" or "pack") { psi.ArgumentList.Add("-m"); psi.ArgumentList.Add("py_compile"); }
                    psi.ArgumentList.Add(file);
                    if (Paths.Python is { } pkg) psi.Environment["PYTHONPATH"] = pkg;
                    psi.Environment["CAPS_LIB"] = MainViewModel.NativeLibraryPath;
                    using var proc = System.Diagnostics.Process.Start(psi)!;
                    var err = proc.StandardError.ReadToEndAsync();
                    proc.StandardOutput.ReadToEnd();
                    proc.WaitForExit(120000);
                    results.Add($"{name} {proc.ExitCode}{(proc.ExitCode != 0 ? " " + err.Result.Trim().Split('\n').LastOrDefault() : "")}");
                }
                Check(results.All(r => r.Split(' ')[1] == "0") && scripts[0].Item2.Contains("doc.relax(") && scripts[1].Item2.Contains("doc.md(steps=50"),
                      "copy as Python: " + string.Join(" · ", results));
                vm.MdStepsD = steps;
            }

            // Compare energies: the LAMMPS files run for zero steps in the LAMMPS on this machine, term by term against CAPS
            {
                vm.CompareEnergies().GetAwaiter().GetResult();
                Check(vm.ParityOk || vm.ParityDetail.Contains("not found"), $"compare energies: {vm.ParityText} · {vm.ParityDetail} · mixing {vm.Field.MixingRule}\n{vm.ParityTable}");
            }

            // the queue beyond Dynamics: a Relax on this structure and a Grow (as its recipe) behind a running MD; stopping
            // the MD starts them one after the other
            {
                vm.MdStepsD = 5_000_000;
                var md = vm.RunMd();
                Until(() => vm.MdRunning && vm.MdLog.StartsWith("step"), 30000);
                vm.QueueRelax();
                vm.GrowChainsD = 1; vm.GrowDpD = 3; vm.GrowDensityD = 0.3m;
                vm.QueueGrow();
                var both = vm.QueuedCount == 2;
                vm.CancelMd();
                Until(() => md.IsCompleted, 30000);
                var relaxJob = vm.Jobs.FirstOrDefault(j => j.Kind == "Relax" && j.Title.StartsWith("Relax"));
                var growJob = vm.Jobs.FirstOrDefault(j => j.Kind == "Recipe");
                var ran = Until(() => relaxJob is { IsDone: true } && growJob is { IsQueued: false, IsRunning: false } && vm.Idle, 180000);
                Check(both && ran && growJob!.IsDone && vm.Title.Contains("cell"),
                      $"queue: relax {relaxJob?.Status} · grow {growJob?.Status} · {vm.Title} · {growJob?.Log.LastOrDefault()?.Text}");
                vm.MdStepsD = steps;
            }
            vm.SetModule(8);
        }

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

        // Nanostructure › Functional groups: a bare (8,8) boron nitride tube, then carboxyls on 4 % of its borons (undoable)
        {
            vm.OpenNano();
            vm.NanoKind = 1;
            vm.NanoMaterial = 1;
            vm.TubeN = 8;
            vm.TubeM = 8;
            vm.TubeLength = 15;
            vm.NanoMatrix = false;
            vm.BuildNano().GetAwaiter().GetResult();
            var bare = vm.Document!.Summary().Atoms;
            vm.FnGroup = "carboxyl";
            vm.FnPattern = 0;
            vm.FnFraction = 0.04m;
            vm.FnElements = "B";
            vm.Functionalize().GetAwaiter().GetResult();
            var grafted = vm.Document!.Summary().Atoms - bare;
            Check(grafted > 0 && grafted % 4 == 0 && vm.Status.Contains("carboxyl groups"), $"functional groups on a BN tube: +{grafted} atoms · {vm.Status}");
            vm.UndoEdit(false);
            Check(vm.Document!.Summary().Atoms == bare, "functional groups undone");
            vm.UndoEdit(true);   // redo: the groups back for the composite
            // the functionalised tube in a natural-rubber matrix: the matrix grows around the structure shown
            var withGroups = vm.Document!.Summary().Atoms;
            vm.NanoMatrix = true;
            vm.NanoAroundShown = true;
            vm.MatrixChains = 3;
            vm.MatrixDp = 6;
            vm.MatrixDensity = 0.5m;
            vm.BuildNano().GetAwaiter().GetResult();
            var composite = vm.Document!.Summary();
            Check(composite.Atoms > withGroups && composite.Molecules == 4 && vm.HoldOn && vm.Document!.Provenance().Contains("nano.embed"),
                  $"composite around the functionalised tube: {composite.Atoms} atoms, {composite.Molecules} molecules · {vm.NanoError} {vm.Status}");
            vm.NanoAroundShown = false;
            vm.NanoMatrix = false;
            vm.NanoMaterial = 0;
        }
        // the builder as a flow (the user's report): graphene sheet settings never built → Graft builds the sheet, then grafts;
        // an element filter the sheet lacks says which elements it has; then Polymer matrix › Build composite wraps it
        {
            vm.OpenNano();
            vm.NanoKind = 0;
            vm.NanoMaterial = 0;
            vm.SheetLayers = 2;
            vm.SheetLx = 18;
            vm.NanoMatrix = false;
            vm.NanoAroundShown = false;
            vm.FnGroup = "*O";
            vm.FnPattern = 0;
            vm.FnFraction = 0.05m;
            vm.FnElements = "B";
            vm.Functionalize().GetAwaiter().GetResult();
            var built = vm.Document?.Summary().Atoms ?? 0;
            var saidText = vm.NanoError;
            var said = saidText.Contains("this structure has C");
            vm.FnElements = "";
            vm.Functionalize().GetAwaiter().GetResult();
            var sheetGrafted = (vm.Document?.Summary().Atoms ?? 0) - built;
            var around = vm.NanoAroundShown;
            vm.NanoMatrix = true;
            vm.MatrixChains = 3; vm.MatrixDp = 6; vm.MatrixDensity = 0.5m;
            vm.BuildNano().GetAwaiter().GetResult();
            var comp2 = vm.Document?.Summary();
            Check(built > 0 && said && sheetGrafted > 0 && around && comp2 is { } c2 && c2.Atoms > built + sheetGrafted && vm.Document!.Provenance().Contains("nano.embed"),
                  $"nano flow: sheet built {built} · B filter said '{saidText}' · +{sheetGrafted} grafted · composite {comp2?.Atoms} · {vm.Status}");
            vm.NanoAroundShown = false;
            vm.NanoMatrix = false;
            vm.FnGroup = "carboxyl";
        }
        // the project tree: Edit opens the builder that made the structure; Delete asks once more, then removes it
        {
            var item = vm.ActiveItem!;
            vm.EditProjectItem(item);
            var toBuilder = vm.IsNano && item.Origin == "Nanostructure builder";
            var itemsBefore = vm.ProjectItems.Count;
            vm.DeleteProjectItem(item);
            var armed = item.DeleteArmed && vm.ProjectItems.Count == itemsBefore;
            vm.DeleteProjectItem(item);
            Check(toBuilder && armed && vm.ProjectItems.Count == itemsBefore - 1 && !vm.ProjectItems.Contains(item),
                  $"project tree: edit → {(vm.IsNano ? "Nanostructure builder" : "module " + vm.Module)} ({item.Origin}) · delete asked {armed} · {itemsBefore} → {vm.ProjectItems.Count}");
            vm.SetModule(8);
        }


        // Mesoscale (DPD): a small A5B5 diblock melt at χN = 43; the frames open as a new structure
        {
            vm.OpenDpd();
            vm.DpdSpecies.Clear();
            vm.DpdSpecies.Add(new DpdSpeciesRow { Name = "diblock", Sequence = "A5B5", Count = 80 });
            vm.DpdSteps = 2000;
            vm.RunDpd().GetAwaiter().GetResult();
            Check(!vm.DpdHasError && vm.IsDpd && vm.DpdOrder != "—" && vm.Document!.Summary().Atoms == 800 && vm.DpdSq.Length > 0
                  && MainViewModel.ExpandSequence("A2B3") == "AABBB" && vm.Document!.Provenance().Contains("dpd.run"),
                  $"DPD: ψ {vm.DpdOrder} · spacing {vm.DpdSpacing} · kT {vm.DpdKt} · {vm.DpdError}");
            vm.SetModule(8);
        }

        // Blend builder: NR / BR 70 : 30
        vm.OpenBlend();
        vm.BlendChains = 4;
        foreach (var r in vm.BlendRows) r.Dp = 8;
        vm.BuildBlend().GetAwaiter().GetResult();
        var bsum = vm.Document?.Summary();
        Check(vm.BlendRows.Count == 2 && bsum is { } blendSum && blendSum.Molecules >= 5 && vm.Title.Contains("blend"), $"blend: {vm.Title} · {bsum?.Molecules} chains · {vm.BlendError}");
        // the same blend grown by configurational bias (Rosenbluth with UFF Lennard-Jones), as Grow offers
        vm.OpenBlend();
        vm.BlendChains = 4;
        foreach (var r in vm.BlendRows) r.Dp = 8;
        vm.BlendGrowMethod = 2;
        vm.BuildBlend().GetAwaiter().GetResult();
        var cbSum = vm.Document?.Summary();
        Check(cbSum is { } cbs && cbs.Molecules >= 4 && vm.BlendError.Length == 0 && vm.BlendGrowBiased,
              $"blend grown by Rosenbluth + UFF LJ: {cbSum?.Molecules} chains · {vm.BlendError}");
        vm.BlendGrowMethod = 0;
        // by chain count: exactly the counts given; by volume: the densities turn volume shares into weight shares
        vm.OpenBlend();
        vm.BlendMode = 2;
        vm.BlendRows[0].Count = 3;
        vm.BlendRows[1].Count = 2;
        vm.BuildBlend().GetAwaiter().GetResult();
        var counted = vm.Document?.Summary().Molecules;
        // a force field per component: both GAFF2 here; after the build Field assigns them by group, one group per component
        {
            vm.OpenBlend();
            var gaffEntry = vm.BlendForceFields.FirstOrDefault(f => f.Id == "gaff-amber25");
            foreach (var r in vm.BlendRows) r.ForceField = gaffEntry;
            vm.BuildBlend().GetAwaiter().GetResult();
            for (int k = 0; k < 200 && vm.Field.Working; ++k) { Thread.Sleep(25); }
            var blendGroups = string.Join(" | ", vm.Field.Groups.Select(g => $"{g.Name}: {g.Molecules}"));
            Check(gaffEntry != null && vm.Field.IsGrouped && vm.Field.Complete && vm.Field.Groups.Count == 2 && vm.Field.Groups[1].Molecules == "4-5" && vm.IsField,
                  $"blend with a force field per component: {blendGroups} · ff {gaffEntry?.Id} · log {vm.BlendLog.Replace("\n", " / ")} · {vm.Field.ForceFieldName} · {vm.Field.Log}");
            foreach (var r in vm.BlendRows) r.ForceField = vm.BlendForceFields[0];
            vm.Field.Clear().GetAwaiter().GetResult();
            vm.Field.Groups.Clear();
            vm.Field.GroupMode = false;
        }
        vm.OpenBlend();
        vm.BlendMode = 1;
        vm.BlendRows[0].Weight = 50;
        vm.BlendRows[1].Weight = 50;
        vm.BlendRows[0].Density = 0.92m;
        vm.BlendRows[1].Density = 0.92m;
        var sameDensity = vm.BlendRows[0].ChainsText;
        vm.BlendRows[1].Density = 1.84m;   // twice as dense: twice the mass for the same volume
        Check(counted == 5 && sameDensity.Contains("vol %") && vm.BlendRows[0].ChainsText.Contains("33.") && vm.BlendRows[1].ChainsText.Contains("66."),
              $"blend by chain count: {counted} chains · by volume: '{sameDensity}' then '{vm.BlendRows[0].ChainsText}' / '{vm.BlendRows[1].ChainsText}'");
        vm.BlendMode = 0;
        // copolymer components: NR with ENR-50 from the library (its preset: isoprene and epoxidised units, 50 : 50), and the
        // chain made in the polymer builder; each chain carries its own units. Cancel goes back to the polymer builder.
        {
            vm.LoadPolymerLibrary();
            vm.SetModule(13);
            vm.OpenBlend();
            var enrEntry = vm.BlendLibrary.FirstOrDefault(e => e.Copolymer && e.Name.StartsWith("ENR-50", StringComparison.Ordinal));
            vm.BlendRows[1].Polymer = enrEntry;
            vm.BlendChains = 3;
            foreach (var r in vm.BlendRows) r.Dp = 10;
            vm.BuildBlend().GetAwaiter().GetResult();
            var enrDoc = vm.Document;
            var oxygen = 0;
            if (enrDoc != null) for (var i = 0; i < (int)enrDoc.Summary().Atoms; ++i) if (enrDoc.Atom(i).ElementSymbol == "O") ++oxygen;
            var enrOk = enrEntry != null && vm.BlendError.Length == 0 && enrDoc != null && oxygen > 0 && vm.Title.Contains("ENR-50", StringComparison.Ordinal);
            // the builder's chain (an alternating A-B copolymer made there) as a third component
            vm.LoadPolymerLibrary();
            vm.SetModule(13);
            vm.OpenBlend();
            var rowsBefore = vm.BlendRows.Count;
            vm.AddBuilderChainToBlend();
            var builderOk = vm.BlendRows.Count == rowsBefore + 1 && vm.BlendRows[^1].Polymer?.Id.StartsWith("builder:", StringComparison.Ordinal) == true && vm.BlendRows[^1].ChainsText.Length > 0;
            vm.RemoveBlendRow(vm.BlendRows[^1]);
            vm.CancelBlend();
            var cancelOk = vm.IsPolymer;
            vm.BlendRows[1].Polymer = vm.BlendLibrary.FirstOrDefault(p => p.Name.StartsWith("Cis-1,4-Polybutadiene", StringComparison.OrdinalIgnoreCase));
            vm.SetModule(8);
            Check(enrOk && builderOk && cancelOk, $"blend with copolymers: [{enrOk} {builderOk} {cancelOk}] {vm.Title} · {oxygen} O · {vm.BlendError}");
        }

        // Packing around a structure, then again with another count: it starts from the structure packed around (not the packed
        // cell) and the new cell replaces the last; another structure opened clears the last result
        {
            vm.Open(Path.Combine(dir, "ps_melt.data"));
            var host = vm.Document!;
            var hostAtoms = (int)host.Summary().Atoms;
            vm.SetModule(5);
            vm.PackStart = 1;
            vm.NewPackInput();
            vm.AddPackMolecule("O=C1OC(=O)C=C1", "maleic anhydride").GetAwaiter().GetResult();
            var row = vm.PackItems.Count - 1;
            vm.SetPackRowCount(row, 4);
            vm.PackAssignField = false;
            vm.RunPack().GetAwaiter().GetResult();
            var first = vm.Document!;
            var firstOk = vm.PackDone && first.Summary().Atoms == hostAtoms + 4 * 9 && vm.PackRepacking && vm.PackStartNote.StartsWith("Packing again", StringComparison.Ordinal);
            var itemsBefore = vm.ProjectItems.Count;
            vm.SetPackRowCount(row, 7);
            vm.RunPack().GetAwaiter().GetResult();
            var second = vm.Document!;
            var secondAtoms = (int)second.Summary().Atoms;
            var againOk = second.Summary().Atoms == hostAtoms + 7 * 9 && vm.ProjectItems.Count == itemsBefore && first.IsDisposed
                          && vm.ProjectItems.Any(i => ReferenceEquals(i.Doc, host));
            // the packed cell (made here, never saved) with a force field, deleted: kept, in Recent, opens with its force field
            vm.Field.FfIndex = vm.Field.Library.ToList().FindIndex(e => e.Id == "gaff-amber25");
            vm.Field.Assign().GetAwaiter().GetResult();
            var packedItem = vm.ProjectItems.First(i => ReferenceEquals(i.Doc, second));
            vm.Activate(vm.ProjectItems.First(i => ReferenceEquals(i.Doc, host)));
            vm.DeleteProjectItem(packedItem);
            vm.DeleteProjectItem(packedItem);
            var kept = RecentFiles.Load().FirstOrDefault(r => r.Path.StartsWith(MainViewModel.KeptFolder, StringComparison.Ordinal));
            var keptOk = kept != null && File.Exists(Path.ChangeExtension(kept.Path, ".ff.json"));
            if (kept != null) vm.Open(kept.Path);
            keptOk &= kept != null && vm.Document?.Summary().Atoms == hostAtoms + 7 * 9 && vm.Field.Assigned;
            // another structure: the page starts fresh
            vm.Open(Path.Combine(dir, "ps_melt.data"));
            var freshOk = !vm.PackDone && vm.PackDmin == "—" && !vm.PackRepacking;
            Check(firstOk && againOk && freshOk && keptOk, $"pack again replaces the last packing; a deleted packed cell kept in Recent: [{firstOk} {againOk} {freshOk} {keptOk}] {kept?.Path} {first.IsDisposed} · {secondAtoms} atoms · {vm.ProjectItems.Count} items · {vm.PackLog.Split('\n')[0]}");
            vm.SetModule(8);
        }

        // packing components each with their own types (as the blend's), the force field the Pack's: groups by type in the
        // LAMMPS input, cross pairs by OPLS's own geometric rule, every number at six decimals
        {
            vm.Open(Path.Combine(dir, "ps_melt.data"));
            vm.SetModule(5);
            vm.PackStart = 1;
            vm.NewPackInput();
            vm.AddPackMolecule("O=C1OC(=O)C=C1", "MAH").GetAwaiter().GetResult();
            vm.SetPackRowCount(vm.PackItems.Count - 1, 3);
            vm.PackFfIndex = vm.Field.Library.ToList().FindIndex(e => e.Id == "opls2005");
            vm.PackAssignField = true;
            vm.RunPack().GetAwaiter().GetResult();
            for (int k = 0; k < 400 && vm.Field.Working; ++k) Thread.Sleep(25);
            var exDir = Path.Combine(outDir, "pack-groups");
            Directory.CreateDirectory(exDir);
            vm.ExportFieldLammps(Path.Combine(exDir, "packed.data")).GetAwaiter().GetResult();
            var lin = File.Exists(Path.Combine(exDir, "packed.in")) ? File.ReadAllText(Path.Combine(exDir, "packed.in")) : "";
            var groupLines = lin.Split('\n').Where(l => l.StartsWith("group ", StringComparison.Ordinal)).ToList();
            var pairs = lin.Split('\n').Where(l => l.StartsWith("pair_coeff", StringComparison.Ordinal)).ToList();
            var sixOk = pairs.Count > 0 && pairs.All(l => System.Text.RegularExpressions.Regex.IsMatch(l, @"^pair_coeff\s+\d+ \d+ \d+\.\d{6} \d+\.\d{6}\s"));
            var grpOk = vm.Field.IsGrouped && groupLines.Count == 2 && groupLines.All(l => l.Contains(" type ", StringComparison.Ordinal)) && groupLines[1].Contains("MAH", StringComparison.Ordinal);
            var mixOk = vm.Field.MixingRule.Contains("geometric", StringComparison.Ordinal) && !vm.Field.MixingRule.Contains("arithmetic", StringComparison.Ordinal);
            Check(sixOk && grpOk && mixOk, $"packed components by group: [{sixOk} {grpOk} {mixOk}] {string.Join(" | ", groupLines)} · {vm.Field.MixingRule} · {pairs.FirstOrDefault()}");
            vm.SetModule(8);
        }

        // React › reactive sites & crosslinking: each chain's sites for the reaction, the links a degree of crosslinking asks
        // for, the crosslinker and its phr (one molecule per link), and a cap on the sites per chain
        {
            // a grown cell: its chains carry repeat-unit numbers
            vm.UsePolystyreneInGrow();
            vm.GrowChainsD = 10;
            vm.GrowDpD = 8;
            vm.GrowDensityD = 0.3m;
            vm.Grow().GetAwaiter().GetResult();
            for (int k = 0; k < 400 && vm.Field.Working; ++k) Thread.Sleep(25);
            vm.SetModule(6);
            vm.RxSet = 5;   // ENR + acid: polystyrene has no epoxide — no sites
            vm.RefreshRxSites();
            var noneOk = vm.RxSitesText.Contains("0 reactive sites each", StringComparison.Ordinal) || vm.RxSitesMax == 0;
            vm.RxSet = 0;   // C–C (CH2–CH2): polystyrene's backbone CH2
            vm.RefreshRxSites();
            var sitesOk = vm.RxHasSites && vm.RxSitesMax > 0 && vm.RxSitesText.StartsWith("10 chains", StringComparison.Ordinal);
            // the degree the panel shows is the run's target without touching it (the 20 % default)
            sitesOk &= vm.RxTargetKind == 5 && Math.Abs((double)vm.RxTargetValueD - 20) < 1e-9;
            vm.RxDcD = 25;
            vm.RxDcD = 20;   // asked for here: the run aims at it
            var calcOk = vm.RxCalcText.Contains("DC 20 %", StringComparison.Ordinal) && !vm.RxHasCrosslinker && vm.RxTargetKind == 5;
            vm.RxSitesPerChainD = 1;
            var capOk = vm.RxCalcText.Contains("allow at most 5 links", StringComparison.Ordinal) || vm.RxCalcText.Contains("up to 5 links", StringComparison.Ordinal);
            var unitsOk = vm.RxCalcText.Contains("of 80 repeat units", StringComparison.Ordinal) && vm.RxCalcText.Contains("→ 8 in this cell", StringComparison.Ordinal);
            // a crosslinker: maleic acid, 116.07 g/mol; phr = 100 N M / m_rubber
            vm.RxCrosslinkerSmiles = "OC(=O)/C=C\\C(=O)O";
            for (int k = 0; k < 200 && !vm.RxCalcText.Contains("maleic", StringComparison.OrdinalIgnoreCase) && !vm.RxCalcText.Contains("116.0", StringComparison.Ordinal); ++k)
            { Avalonia.Threading.Dispatcher.UIThread.RunJobs(); Thread.Sleep(25); }
            var phrOk = vm.RxCalcText.Contains("116.0", StringComparison.Ordinal) && vm.RxPhrXD > 0;
            // phr back to the degree: the same numbers
            var phr = vm.RxPhrXD;
            vm.RxPhrXD = phr;
            var backOk = Math.Abs((double)vm.RxDcD - 20) < 1.5;
            vm.RxSitesPerChainD = 0;
            vm.RxCrosslinkerSmiles = "";
            vm.RxTargetKind = 0;
            Check(noneOk && sitesOk && calcOk && capOk && unitsOk && phrOk && backOk, $"react sites & crosslinking: [{noneOk} {sitesOk} {calcOk} {capOk} {unitsOk} {phrOk} {backOk}] {vm.RxSitesText} · {vm.RxCalcText.Replace("\n", " / ")}");
            vm.SetModule(8);
        }

        // the Force field tab's dot is the active structure's: on with a complete assignment, off once the last structure closes;
        // layers empty with no structure; a homopolymer chain named poly(UNIT)
        {
            while (vm.Document != null) vm.CloseDocument();
            vm.Open(Path.Combine(dir, "ps_melt.data"));
            vm.Field.FfIndex = vm.Field.Library.ToList().FindIndex(e => e.Id == "gaff-amber25");
            vm.Field.Assign().GetAwaiter().GetResult();
            var on = vm.StepFieldDone && vm.StepFieldTip.Contains("every atom typed", StringComparison.Ordinal);
            while (vm.Document != null) vm.CloseDocument();
            var off = !vm.StepFieldDone && !vm.StepFieldWarn && !vm.HasLayers && vm.LayerRows.Count == 0;
            Check(on && off, $"force-field dot and layers follow the structure: on {on} · off after closing {off} · {vm.StepFieldTip}");
        }

        // AMBER prmtop: opened with its restart beside it, its own force field assigned; another force field, then back
        {
            var amber = Path.GetFullPath(Path.Combine(dir, "..", "tests", "data", "amber", "phenol.prmtop"));
            if (File.Exists(amber))
            {
                vm.Open(amber);
                var own = vm.Field.Assigned && vm.Field.Complete && vm.Field.ForceFieldName.Contains("phenol.prmtop", StringComparison.Ordinal) && !vm.Field.CanUseFileForceField;
                var gaffIx = vm.Field.Library.ToList().FindIndex(x => x.Id == "gaff-amber25");
                var amberOffered = false;
                if (gaffIx >= 0)
                {
                    vm.Field.FfIndex = gaffIx;
                    vm.Field.Assign().GetAwaiter().GetResult();
                    amberOffered = vm.Field.CanUseFileForceField;
                    vm.Field.UseFileForceField().GetAwaiter().GetResult();
                }
                Check(own && amberOffered && vm.Field.ForceFieldName.Contains("phenol.prmtop", StringComparison.Ordinal) && !vm.Field.CanUseFileForceField,
                      $"AMBER prmtop opens with its own force field, offered back after GAFF: own {own}, offered {amberOffered} · {vm.Field.ForceFieldName} · {vm.Field.Log}");
            }
            else Check(false, "AMBER test topology missing: " + amber);
        }

        // Coarse-grained Kremer–Grest melt: reduced units, then mapped to real units by σ, T and the bead mass
        {
            vm.CgChains = 10;
            vm.CgBeads = 20;
            vm.OpenCg();
            for (var i = 0; i < 200 && vm.CgDoc == null; i++) { Avalonia.Threading.Dispatcher.UIThread.RunJobs(); Thread.Sleep(25); }
            var cgStem = Path.Combine(outDir, "caps-selftest-kg");
            vm.CgUnits = 0;
            vm.ExportCg(cgStem);
            var lj = File.Exists(cgStem + ".in") ? File.ReadAllText(cgStem + ".in") : "";
            vm.CgUnits = 1;
            vm.CgSigma = 5m; vm.CgTemp = 300m; vm.CgMass = 50m;
            vm.ExportCg(cgStem + "-real");
            var real = File.Exists(cgStem + "-real.in") ? File.ReadAllText(cgStem + "-real.in") : "";
            Check(lj.Contains("units lj") && real.Contains("units real") && real.Contains("bond_coeff 1 " + (30 * 0.0019872067 * 300 / 25).ToString("0.########", System.Globalization.CultureInfo.InvariantCulture)) &&
                  vm.CgMapText.StartsWith("ε = k_B T = 0.5962", StringComparison.Ordinal),
                  $"CG units: reduced and mapped · {vm.CgMapText[..Math.Min(90, vm.CgMapText.Length)]}");
            vm.CgUnits = 0;
            // the melt as the Studio structure comes with its own force field: nothing to assign before exporting
            vm.BuildCg();
            var kgAssigned = vm.Field.Assigned && vm.Field.Complete && vm.Field.ForceFieldName.StartsWith("Kremer–Grest", StringComparison.Ordinal);
            Check(kgAssigned, $"Kremer–Grest melt assigned on its own: {vm.Field.ForceFieldName} · {vm.Field.Log}");
            // reduced units: backmapping is refused (no lengths in Å); mapped to styrene units: 16 atoms a bead, GAFF2 typed
            vm.CgBackmapUnit = vm.CgBackmapPolymers.FindIndex(p => p.Name == "Polystyrene");
            vm.BackmapCg().GetAwaiter().GetResult();
            var reducedRefused = vm.CgError.Contains("reduced units");
            vm.CgUnits = 1;
            vm.CgChains = 4; vm.CgBeads = 12;
            vm.CgSigma = 5.5m; vm.CgTemp = 450m; vm.CgMass = 104.15m;
            vm.BuildCg();
            vm.BackmapCg().GetAwaiter().GetResult();
            for (int k = 0; k < 400 && vm.Field.Working; ++k) Thread.Sleep(25);
            var bm = vm.Document?.Summary();
            Check(reducedRefused && bm is { } bms && bms.Atoms == 4 * 12 * 16 + 8 && vm.Field.Complete && vm.CgBackmapLog.Contains("relaxed"),
                  $"KG backmap to polystyrene: refused in reduced units {reducedRefused} · {bm?.Atoms} atoms · {vm.Field.ForceFieldName} · {vm.CgBackmapLog.Replace("\n", " / ")} {vm.CgError}");
            vm.CgUnits = 0;
            vm.CgChains = 10; vm.CgBeads = 20;
            // MARTINI: PEO chains of SN0 beads with the library's MARTINI 2 polymers, packed and compressed to 1.1 g/cm³
            vm.OpenCg();
            vm.CgModel = 1;
            vm.MtExample = 1;
            vm.MtRepeats = 10; vm.MtChains = 12; vm.MtDensity = 1.1m;
            vm.BuildMartini().GetAwaiter().GetResult();
            var mtDensity = vm.Document?.Summary().Density ?? 0;
            Check(vm.Field.Complete && vm.Field.ForceFieldName.Contains("MARTINI", StringComparison.Ordinal) && Math.Abs(mtDensity - 1.1) < 0.01,
                  $"MARTINI melt: {vm.Field.ForceFieldName} · {mtDensity:0.000} g/cm³ · {vm.Status} {vm.CgError}");
            // From a polymer: polystyrene, backbone + side group, from an all-atom reference melt
            vm.CgModel = 2;
            vm.MpPolymer = vm.CgBackmapPolymers.FindIndex(p => p.Name == "Polystyrene");
            vm.MpScheme = 1;
            vm.MpChains = 6; vm.MpDp = 12; vm.MpDensity = 1.04m; vm.MpTemp = 450;
            vm.BuildMappedCg().GetAwaiter().GetResult();
            var mp = vm.Document?.Summary();
            Check(mp is { } mps && mps.Atoms == 6 * 12 * 2 && vm.Field.Complete && vm.MpLog.Contains("bond STY_B–STY_S") && vm.Field.ForceFieldName.StartsWith("Structure-based CG"),
                  $"CG from a polymer: {mp?.Atoms} beads · complete {vm.Field.Complete} · log {vm.MpLog.Split((char)10).LastOrDefault()} · {vm.Field.ForceFieldName} · {vm.CgError}");
            vm.CgModel = 0;
            vm.SetModule(8);
        }

        // Settings › Accessibility and Python & scripting
        {
            vm.SetHighContrast = true;
            var strong = vm.ViewOptions(64, 64, 1).Outlines;
            var hcBuf = new byte[160 * 120 * 4];
            vm.Document?.Render(vm.Camera, vm.ViewOptions(160, 120, 1), hcBuf);
            vm.SetHighContrast = false;
            var light = vm.ViewOptions(64, 64, 1).Outlines;
            vm.CheckPython().GetAwaiter().GetResult();
            var py = vm.PythonCheck;
            vm.Settings.Shortcuts[vm.PaletteRowsFor("theory manual").FirstOrDefault(r => !r.IsHeader && r.Id.Length > 0)?.Id ?? "x"] = OperatingSystem.IsMacOS() ? "Meta+K" : "Ctrl+K";
            vm.FillShortcutConflicts();
            var conflicts = vm.ShortcutConflicts.Count;
            if (vm.ShortcutConflicts.FirstOrDefault() is { } sc) vm.ResolveShortcutConflict(sc);
            Check(strong == 2 && light == 1 && py.StartsWith("✓", StringComparison.Ordinal) && py.Contains("caps package: ABI", StringComparison.Ordinal) && conflicts == 1 && vm.ShortcutConflicts.Count == 0,
                  $"settings: outlines {strong}/{light} · python {py.Split('\n')[0]} · conflicts {conflicts} → {vm.ShortcutConflicts.Count}");
        }

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
        // the crystal saved for VASP (POSCAR, species grouped: 24 Ti then 48 O) and as a CIF, both read back
        {
            var vp = Path.Combine(outDir, "caps-selftest-rutile.vasp");
            var cp = Path.Combine(outDir, "caps-selftest-rutile.cif");
            vm.SaveDocument(vp);
            vm.SaveDocument(cp);
            var pl = File.Exists(vp) ? File.ReadAllLines(vp) : [];
            var cif = File.Exists(cp) ? File.ReadAllText(cp) : "";
            Check(pl.Length == 8 + 72 && pl[5].Split(' ', StringSplitOptions.RemoveEmptyEntries).SequenceEqual(new[] { "Ti", "O" })
                  && pl[6].Split(' ', StringSplitOptions.RemoveEmptyEntries).SequenceEqual(new[] { "24", "48" }) && pl[7] == "Direct"
                  && cif.Contains("_chemical_formula_sum 'O48 Ti24'") && cif.Contains("_space_group_IT_number 1")
                  && MainViewModel.ExportFormats.Any(f => f.Id == "poscar"),
                  $"crystal saved as POSCAR ({(pl.Length > 6 ? pl[5] + " / " + pl[6] : "none")}) and CIF · {vm.Status}");
        }
        // the rotate tool's drag: on an orthographic view (screen = 10 x + 100, −10 y + 200), a drag from right of the
        // centre to above it turns world +x into +y (+90° about z, whatever sign the axis comes with); X locks the world x axis
        {
            var rng = new Random(4);
            var pts = Enumerable.Range(0, 20).Select(_ => { double x = rng.NextDouble() * 10, y = rng.NextDouble() * 10, z = rng.NextDouble() * 10; return (x, y, z, 10 * x + 100, -10 * y + 200); }).ToList();
            double[] c = [5, 5, 5];   // on screen (150, 150)
            var rot = MainViewModel.RotationFromDrag(pts, c, 200, 150, 150, 100, '\0', false)!.Value;
            // Rodrigues: +x turned by the angle about the axis
            double th = rot.Degrees * Math.PI / 180, ux = rot.Axis[0], uy = rot.Axis[1], uz = rot.Axis[2];
            double[] v = [1, 0, 0];
            var dot = ux * v[0] + uy * v[1] + uz * v[2];
            double[] cr = [uy * v[2] - uz * v[1], uz * v[0] - ux * v[2], ux * v[1] - uy * v[0]];
            double[] r = Enumerable.Range(0, 3).Select(k => v[k] * Math.Cos(th) + cr[k] * Math.Sin(th) + new[] { ux, uy, uz }[k] * dot * (1 - Math.Cos(th))).ToArray();
            var locked = MainViewModel.RotationFromDrag(pts, c, 200, 150, 260, 150, 'x', true)!.Value;
            Check(Math.Abs(r[0]) < 1e-6 && Math.Abs(r[1] - 1) < 1e-6 && Math.Abs(Math.Abs(uz) - 1) < 1e-9 && locked.Axis[0] == 1 && locked.Degrees == 30,
                  $"rotate tool: +x → ({r[0]:0.###}, {r[1]:0.###}, {r[2]:0.###}) by {rot.Degrees:0.#}° about ({ux:0.#}, {uy:0.#}, {uz:0.#}) · locked x {locked.Degrees}°");
        }
        // cell tools on the open structure: the 2 × 2 × 3 rutile folds back to its 6-atom cell, the conventional cell is
        // P 42/m n m (136), then a vacuum slab and a [001] wire; a matrix that is not a lattice map is refused
        {
            vm.FindPrimitiveCell();
            var prim = vm.Document!.Summary().Atoms;
            vm.ConventionalCell();
            var convText = vm.CellToolText;
            vm.MakeVacuumSlab();
            var slab = vm.CellToolText;
            vm.RedefineMatrix = "1/2 0 0\n0 1 0\n0 0 1";
            vm.RedefineLattice();
            var refused = vm.CellToolText;
            vm.WireUvw = "0 0 1"; vm.WireRadius = 6; vm.WireRepeats = 3;
            vm.MakeNanowire();
            Check(prim == 6 && convText.Contains("No. 136") && slab.StartsWith("Vacuum slab") && refused.Contains("does not map") && vm.CellToolText.StartsWith("Nanowire"),
                  $"cell tools: primitive {prim} atoms · {convText} · {slab} · {refused} · {vm.CellToolText}");
        }
        // the Modify toolbar on the picked atoms: an element, a bond order, a geometry; each one undoable edit
        {
            var d = vm.Document!;
            var bonds = System.Text.Json.Nodes.JsonNode.Parse(d.BondLabels("index"))!["pairs"]!.AsArray().Select(x => (int)x!.GetValue<double>()).ToArray();
            var (i0, j0) = (bonds[0], bonds[1]);
            vm.Pick(i0, false);
            vm.ModifyElementPicked("Si");
            var el = System.Text.Json.Nodes.JsonNode.Parse(d.AtomLabels("element"))!.AsArray()[i0]!.GetValue<string>();
            vm.Pick(i0, false); vm.Pick(j0, true);
            vm.BondOrderPicked(2);
            var order = System.Text.Json.Nodes.JsonNode.Parse(d.BondLabels("order"))!["labels"]!.AsArray()[0]!.GetValue<string>();
            vm.Pick(i0, false);
            vm.GeometryPicked("octahedral");
            var geomOk = vm.EditError.Length == 0 || vm.EditError.Contains("ring");   // a crystal's atoms are all in rings: refused, and said why
            Check(el == "Si" && order == "2" && geomOk, $"modify toolbar: element {el} · bond order {order} · geometry {(geomOk ? "ok" : vm.EditError)}");
        }
        // two dump files of one run, chosen together: joined in time order, the boundary frame kept once
        {
            var lines = File.ReadAllLines(Path.Combine(dir, "ps_melt.lammpstrj"));
            var starts = lines.Select((l, k) => (l, k)).Where(x => x.l.StartsWith("ITEM: TIMESTEP")).Select(x => x.k).ToList();
            var pa = Path.Combine(outDir, "caps-selftest-part1.lammpstrj");
            var pb = Path.Combine(outDir, "caps-selftest-part2.lammpstrj");
            File.WriteAllLines(pa, lines.Take(starts[2]));
            File.WriteAllLines(pb, lines.Skip(starts[1]));
            vm.OpenJoined([pb, pa], Path.Combine(dir, "ps_melt.data"));
            Check(vm.Frames == 3 && vm.Status.Contains("2 files joined"), $"joined dumps: {vm.Frames} frames · {vm.Status}");
        }
        // File checks › Atom types: each type's element and mass, editable (the type kept); undo restores it
        {
            vm.Open(Path.Combine(dir, "ps_melt.data"));
            vm.OpenChecks();
            var typeCount = vm.TypeRows.Count;
            var h = vm.TypeRows.FirstOrDefault(r => r.Element == "H");
            var set = false;
            if (h != null)
            {
                h.Element = "Si";
                vm.ApplyTypeRow(h);
                set = vm.TypeRows.Any(r => r.Type == h.Type && r.Element == "Si" && r.Mismatch);   // the mass still points to H
                vm.UndoEdit(false);
            }
            vm.OpenChecks();
            var back = vm.TypeRows.Any(r => r.Element == "H");
            Check(typeCount >= 2 && set && back, $"atom types: {typeCount} types · set Si {set} · undone {back}");
            // Pack around it, filled to 0.6 g/cm³ with water: N = (0.6 × 35937 × 0.60221 − 8352.3) / 18.015 = 257
            var w = Path.Combine(outDir, "caps-selftest-one-water.xyz");
            File.WriteAllText(w, "3\nwater\nO 0 0 0\nH 0.9572 0 0\nH -0.2400 0.9266 0\n");
            vm.PackStart = 1;
            vm.NewPackInput();
            vm.AddPackStructure(w);
            vm.FillDensity = 0.6m;
            vm.FillToDensity();
            Check(vm.PackText.Contains("count 257") && vm.FillText.StartsWith("257 ×"), $"fill to density: {vm.FillText}");
            vm.PackStart = 0;
            // the Python console: doc = the open structure; a line runs; doc comes back as a new structure
            {
                vm.StartConsole();
                bool Wait(string what) { for (var t = 0; t < 150; ++t) { if (vm.ConsoleText.Contains(what)) return true; Thread.Sleep(100); } return false; }
                var started = Wait("1300 atoms)");
                vm.ConsoleInput = "print('elements', sorted(set(doc.atom_labels('element'))))";
                vm.ConsoleRun();
                var ran = Wait("elements ['C', 'H']");
                var itemsBefore = vm.ProjectItems.Count;
                vm.ConsoleOpenDoc().GetAwaiter().GetResult();
                var docBack = vm.ProjectItems.Count == itemsBefore + 1;
                vm.StopConsole();
                Check(started && ran && docBack, $"python console: started {started} · ran {ran} · doc back {docBack} · {vm.ConsoleText.Trim().Split('\n').LastOrDefault()}");
            }
            vm.SetModule(8);
        }

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
            // labels: every atom and bond kind answers; the first four still switch their kinds; the style follows the settings
            var bad = vm.AtomLabelKinds.Where(k => { try { return System.Text.Json.Nodes.JsonNode.Parse(adoc.AtomLabels(k.Id))!.AsArray().Count != adoc.Summary().Atoms; } catch { return true; } })
                .Concat(vm.BondLabelKinds.Where(k => { try { return System.Text.Json.Nodes.JsonNode.Parse(adoc.BondLabels(k.Id))!["labels"]!.AsArray().Count != adoc.Summary().Bonds; } catch { return true; } }))
                .Select(k => k.Id).ToList();
            vm.LabelCharge = true;
            var chargeOn = vm.AtomLabelKinds.First(k => k.Id == "charge").On && vm.AnyLabels;
            vm.LabelFont = 1; vm.LabelSize = 14; vm.AtomLabelColour = Array.IndexOf(MainViewModel.LabelColours, "Custom"); vm.AtomLabelHex = "#123456";
            var lk = vm.LabelLook;
            var bondOrders = System.Text.Json.Nodes.JsonNode.Parse(adoc.BondLabels("chemical"))!["labels"]!.AsArray().Select(x => x!.GetValue<string>()).Distinct().ToList();
            Check(vm.AtomLabelKinds.Count >= 30 && vm.BondLabelKinds.Count >= 12 && bad.Count == 0 && chargeOn && lk.Font == "IBM Plex Sans" && lk.Size == 14
                  && lk.AtomArgb == 0xFF123456 && bondOrders.Contains("C=O") && bondOrders.Contains("N–H"),
                  $"labels: {vm.AtomLabelKinds.Count} atom and {vm.BondLabelKinds.Count} bond kinds, failing [{string.Join(",", bad)}] · style {lk} · bonds {string.Join(" ", bondOrders.Take(8))}");
            vm.LabelCharge = false; vm.LabelFont = 0; vm.LabelSize = 10.5m; vm.AtomLabelColour = 0;
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
            // auto-clean (A): a placed atom is followed by a UFF clean-up of its neighbourhood, its own undo step
            if (!vm.AutoCleanOn) vm.ToggleAutoClean();
            vm.BuildElement = "O";
            vm.EditTool = 1;
            var carbonAt = Enumerable.Range(0, (int)vm.Document!.Summary().Atoms).First(i => vm.Document!.Atom(i).ElementSymbol == "C");
            var hist0 = vm.EditHistory.Count;
            vm.ToolClick(carbonAt);
            var cleaned = vm.Status.Contains("auto-cleaned") && vm.EditHistory.Count == hist0 + 2;
            vm.EditTool = 0;
            vm.UndoEdit(false);
            vm.UndoEdit(false);
            vm.ToggleAutoClean();
            Check(cleaned && !vm.AutoCleanOn && vm.Document!.Summary().Atoms == n0, $"auto-clean: {vm.EditHistory.Count} history · {(cleaned ? "cleaned after placing O" : "not cleaned")}");
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
            // exact geometry from the live monitor: a C–H bond picked, set to 1.200 Å (the hydrogen moves), undone
            {
                var hc = vm.Document!.Neighbours(carbon, 4).First(nb => vm.Document!.Atom(nb.Item1).ElementSymbol == "H").Item1;
                vm.Pick(carbon);
                vm.Pick(hc, true);
                var chBefore = vm.Document!.Measure([carbon, hc]);
                vm.MeasureTarget = "1.200";
                vm.SetMeasured();
                var after = vm.Document!.Measure([carbon, hc]);
                var kept = vm.CanSetMeasure && vm.MeasureText.Contains("1.200");
                vm.UndoEdit(false);
                var chUndone = vm.Document!.Measure([carbon, hc]);
                Check(Math.Abs(after - 1.2) < 1e-6 && Math.Abs(chUndone - chBefore) < 1e-9 && kept, $"set exact bond length: {chBefore:F3} → {after:F3} Å · undo {chUndone:F3} · monitor '{vm.MeasureText}'");
            }
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

        // Timeline: the time series of a pipeline over the three frames, with a centred running mean over it
        {
            vm.SetModule(20);
            vm.ClearPipeline();
            vm.AddStep("coordination");
            vm.ComputeSeries().GetAwaiter().GetResult();
            vm.SeriesWindow = 3;
            var spark = vm.Sparkline();
            var mean = vm.SparklineMean();
            var inRange = spark.Length > 0 && mean.All(m => m.Y >= spark.Min(p => p.Y) - 1e-9 && m.Y <= spark.Max(p => p.Y) + 1e-9);
            Check(vm.HasTimeline && spark.Length == 3 && mean.Length == 3 && inRange && Math.Abs(mean[1].Y - spark.Average(p => p.Y)) < 1e-9,
                  $"timeline running mean: {spark.Length} frames of {vm.SparkAttribute}, mean {mean.Length}");
            vm.SeriesWindow = 1;
            // manual selection: the lassoed atoms go into the step as ranges, and select the same atoms in the result
            vm.LassoSelect([3, 4, 5, 9], false);
            vm.AddStep("manual_selection");
            var manual = vm.PipelineRows.FirstOrDefault(r => r.Type == "manual_selection");
            Check(manual != null && (string?)manual.Params["atoms"] == "3-5 9" && manual.Summary.Contains("4 selected"),
                  $"manual selection: '{manual?.Params["atoms"]}' · {manual?.Summary}");
            if (manual != null) vm.PipelineRows.Remove(manual);
            vm.ClearDocSelection();
            // the per-frame cache: back on a frame, its result comes from the cache; a new setting empties it
            vm.ApplyPipeline();
            vm.Frame = 1;
            var computed = vm.PipeCacheText;
            vm.Frame = 0;
            var hit = vm.PipeCacheText;
            Check(computed.StartsWith("computed") && hit.StartsWith("from the cache") && hit.Contains("2 frames"), $"pipeline cache: '{computed}' then '{hit}'");
            // branches: a compute step only in "Numbers"; the shared coordination step runs once, switching back is cached
            vm.AddStep("compute_property");
            var numbers = vm.PipelineRows.First(r => r.Type == "compute_property");
            numbers.Params["branch"] = "Numbers";
            vm.ApplyPipeline();
            var shownNumbers = vm.HasPipeBranches && vm.PipeBranchIndex == 1 && !numbers.Summary.Contains("not shown");
            vm.PipeBranchIndex = 0;
            var hidden = numbers.Summary.Contains("branch Numbers · not shown");
            vm.PipeBranchIndex = 1;
            Check(shownNumbers && hidden && vm.PipeCacheText.StartsWith("from the cache"), $"pipeline branches: {string.Join(" | ", vm.PipeBranches)} · '{numbers.Summary}' · {vm.PipeCacheText}");
            vm.PipelineRows.Remove(numbers);
            vm.ApplyPipeline();
            // outputs: the rdf table as CSV and a plot, in the saved YAML and written now
            vm.OpenSavePipeline();
            vm.AddPipelineOutput("table");
            vm.AddPipelineOutput("plot");
            var outDirP = Path.Combine(outDir, "caps-selftest-outputs");
            vm.WritePipelineOutputs(outDirP);
            var wroteCsv = vm.PipelineOutputs.Count == 2 && File.Exists(Path.Combine(outDirP, vm.PipelineOutputs[0].Path));
            var wroteSvg = vm.PipelineOutputs.Count == 2 && File.Exists(Path.Combine(outDirP, vm.PipelineOutputs[1].Path));
            Check(vm.PipelineYaml.Contains("outputs:") && vm.PipelineYaml.Contains("- plot: rdf -> rdf.svg") && wroteCsv && wroteSvg,
                  $"pipeline outputs: {string.Join(", ", vm.PipelineOutputs.Select(o => $"{o.Kind} {o.What} -> {o.Path}"))} · {vm.Status}");
            while (vm.PipelineOutputs.Count > 0) vm.RemovePipelineOutput(vm.PipelineOutputs[0]);
            vm.SetModule(20);
            // the inspector as CSV: every particle matching the filter (not only the page), and a data table
            vm.InspectorTab = 0;
            vm.InspectorFilter = "Molecule == 1";
            var csv = Path.Combine(outDir, "caps-selftest-particles.csv");
            vm.ExportInspectorCsv(csv);
            var lines = File.Exists(csv) ? File.ReadAllLines(csv) : [];
            vm.InspectorTab = 3;
            var tcsv = Path.Combine(outDir, "caps-selftest-table.csv");
            vm.ExportInspectorCsv(tcsv);
            var tlines = File.Exists(tcsv) ? File.ReadAllLines(tcsv) : [];
            Check(lines.Length == 131 && lines[0].Contains("Molecule") && tlines.Length > 2, $"inspector CSV: {lines.Length - 1} particles of molecule 1 · table {tlines.Length - 1} rows · {vm.Status}");
            vm.InspectorFilter = "";
            vm.InspectorTab = 0;
            vm.ClearPipeline();
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
            vm.WaitAppearance();
            var previewOk = vm.ChargePreviewing && vm.Document!.AppearanceInfo().Contains("\"preview\":true") && Math.Abs(vm.Document!.Atom(0).Charge - q0) < 1e-12;
            vm.ApplyCharges();
            previewOk &= !vm.ChargePreviewing;
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
            vm.WaitAppearance();
            Check(previewOk && !vm.ChargePreviewing && !vm.Document!.AppearanceInfo().Contains("\"preview\":true"),
                  $"charges preview: the view shows the computed QEq charges before Apply, the structure's after ({previewOk})");
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
            // over the frames (when this structure has them): one total per frame, mean ± sd
            if (vm.HasFrames)
            {
                vm.SaOverFrames = true;
                vm.RunSurfaceArea().GetAwaiter().GetResult();
                sasa = sasa && vm.SaSeries.Length == vm.Frames && vm.SaSeriesText.Contains("±");
                vm.SaOverFrames = false;
            }
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
            // log-normal and a measured histogram
            vm.PdDist = 4;
            var logNormal = vm.PdK.StartsWith("ln N normal") && vm.PdLengths.Length == vm.GrowChains;
            vm.PdDist = 5;
            vm.PdHistogram = "10:1, 30:1";
            var hist = vm.PdLengths.All(x => x is "10" or "30") && vm.PdK.Contains("Nₙ 20");
            Check(logNormal && hist, $"polydispersity log-normal {logNormal} · histogram {hist}: {string.Join(",", vm.PdLengths)} · {vm.PdK}");
            vm.PdDist = 0;
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

        // A1 (design/boards/SelectionBar): one selection whatever made it; the bar's count, what and how; hide, ghost, show
        // only and show all as view states (the structure keeps every atom); commands act on the whole selection
        {
            vm.Open(Path.Combine(dir, "ps_melt.data"));
            var d = vm.Document!;
            var n = (int)d.Summary().Atoms;
            var mol = d.MoleculeIds();
            var a0 = 5;
            var molAtoms = mol.Count(m => m == mol[a0]);
            vm.SelectLike(a0, false, false);
            vm.RefreshSelBar();
            var molOk = vm.SelBarCount == molAtoms && vm.SelBarWhat == $"molecule {mol[a0]}" && vm.SelBarHow == "double-click";
            // a click adds a second molecule's atom: the bar counts both, made by double-click + click
            var other = Enumerable.Range(0, n).First(i => mol[i] != mol[a0]);
            vm.Pick(other, true);
            vm.RefreshSelBar();
            var unionOk = vm.SelBarCount == molAtoms + 1 && vm.SelBarHow == "double-click + click" && vm.SelectionAtoms().Contains(other);
            // hide: out of the view (and the selection), still in the structure
            vm.HideSelection();
            var st = d.AtomStates();
            var hideOk = vm.HiddenCount == molAtoms + 1 && st[a0] == 2 && st[other] == 2 && d.Summary().Atoms == n && vm.SelBarCount == 0 && vm.HasHiddenAtoms;
            // a lasso over hidden atoms does not take them
            vm.LassoSelect(new[] { a0, other, (a0 + 200) % n }, false);
            vm.RefreshSelBar();
            var lassoOk = !vm.SelectionAtoms().Contains(a0) && vm.SelBarHow == "lasso";
            vm.ShowAllAtoms();
            var showOk = vm.HiddenCount == 0 && d.AtomStates().All(v => v == 0);
            // same type (⌥ double-click), ghost, show only
            var t0 = d.Atom(a0).Type;
            vm.SelectLike(a0, true, false);
            vm.RefreshSelBar();
            var typeCount = Enumerable.Range(0, n).Count(i => d.Atom(i).Type == t0);
            var typeOk = vm.SelBarCount == typeCount && vm.SelBarHow.StartsWith($"type {t0}", StringComparison.Ordinal);
            vm.GhostSelection();
            var ghostOk = vm.GhostCount == typeCount;
            vm.ShowAllAtoms();
            vm.SelectLike(a0, false, false);
            vm.ShowOnlySelection();
            var onlyOk = vm.HiddenCount == n - molAtoms;
            vm.ShowAllAtoms();
            // delete acts on the whole selection (a molecule, not the four picks), undoably
            vm.SelectLike(a0, false, false);
            vm.DeleteSelection();
            var delOk = vm.Document!.Summary().Atoms == n - molAtoms;
            vm.UndoEdit(false);
            delOk &= vm.Document!.Summary().Atoms == n;
            var xyzOk = true;
            vm.SelectLike(a0, false, false);
            vm.RefreshSelBar();
            var xyz = vm.SelectionXyz().Split('\n', StringSplitOptions.RemoveEmptyEntries);
            xyzOk = xyz.Length == molAtoms + 2 && xyz[0] == molAtoms.ToString(System.Globalization.CultureInfo.InvariantCulture);
            vm.ClearAllSelection();
            vm.RefreshSelBar();
            var clearOk = vm.SelBarCount == 0 && !vm.HasSelBar;
            Check(molOk && unionOk && hideOk && lassoOk && showOk && typeOk && ghostOk && onlyOk && delOk && xyzOk && clearOk,
                  $"A1 selection bar: [{molOk} {unionOk} {hideOk} {lassoOk} {showOk} {typeOk} {ghostOk} {onlyOk} {delOk} {xyzOk} {clearOk}] molecule {molAtoms} atoms · type {t0}: {typeCount} · '{vm.SelBarWhat}'");
        }

        // A2/A9 layers (design/boards/Layers): kinds and their molecules, the eye cycle, the lock, rows lit by the selection
        {
            vm.Open(Path.Combine(dir, "ps_melt.data"));
            vm.RefreshLayers();
            var kind = vm.LayerRows.FirstOrDefault();
            var listOk = kind != null && kind.IsKind && kind.HasChildren && kind.Name == "10 × C64H66" && kind.Atoms == 1300 && vm.LayerRows.Count == 11
                         && kind.Z.Max() == 1.0 && vm.LayerRows.Skip(1).All(r => r.Atoms == 130);
            var m3 = vm.LayerRows[3];
            vm.CycleLayer(m3);
            var ghostOk = m3.State == "ghost" && kind!.State == "mixed" && vm.GhostCount == 130 && vm.LayerChip.Contains("ghosted", StringComparison.Ordinal);
            vm.CycleLayer(m3);
            var hiddenOk = m3.State == "hidden" && vm.HiddenCount == 130;
            vm.CycleLayer(m3);
            var backOk = m3.State == "shown" && vm.HiddenCount + vm.GhostCount == 0;
            var m1 = vm.LayerRows[1];
            vm.LockLayer(m1);
            var d = vm.Document!;
            var lockOk = m1.Locked && d.AtomStates().Take(1300).Count(v => (v & 4) != 0) == 130;
            // a locked layer is not selected by its row, the whole kind selects the rest
            vm.SelectLayer(m1, false);
            var lockSelOk = vm.SelectedCount == 0;
            vm.SelectLayer(kind!, false);
            vm.RefreshSelBar();
            var selOk = vm.SelectedCount == 1170 && kind.IsLit && !m1.IsLit && vm.LayerRows[2].IsFullyLit;
            vm.LockLayer(m1);
            vm.ToggleLayer(kind);
            var foldOk = vm.LayerRows.Count == 1;
            vm.ToggleLayer(kind);
            vm.ClearAllSelection();
            Check(listOk && ghostOk && hiddenOk && backOk && lockOk && lockSelOk && selOk && foldOk && !m1.Locked,
                  $"A2 layers: [{listOk} {ghostOk} {hiddenOk} {backOk} {lockOk} {lockSelOk} {selOk} {foldOk}] {kind?.Name} · {vm.LayerRows.Count} rows · chip '{vm.LayerChip}'");
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
            // runs are recorded too: a short MD becomes doc.md(steps=40, …) (the call alone; the macro opens and saves)
            vm.NewMacro();
            vm.Recording = true;
            var mdSteps0 = vm.MdStepsD;
            vm.MdStepsD = 40;
            vm.RunMd().GetAwaiter().GetResult();
            vm.MdStepsD = mdSteps0;
            vm.Recording = false;
            var mdLine = vm.RecordedCommands.FirstOrDefault(c => c.Text.StartsWith("doc.md("))?.Text ?? "";
            Check(mdLine.Contains("steps=40") && !vm.RecordedCommands.Any(c => c.Text.StartsWith("import") || c.Text.StartsWith("doc = caps.open")),
                  $"macro records MD: {mdLine}");
            // an equilibration records doc.equilibrate with the stages as run (40 steps)
            var eqText0 = vm.EqText;
            var eqUntil0 = vm.EqUntilConverged;
            vm.EqUntilConverged = false;
            vm.EqText = "nvt 0.04 ps T 300 # hold";
            vm.Recording = true;
            vm.RunEquilibrate().GetAwaiter().GetResult();
            vm.Recording = false;
            vm.EqText = eqText0;
            vm.EqUntilConverged = eqUntil0;
            var eqLine = vm.RecordedCommands.FirstOrDefault(c => c.Text.StartsWith("doc.equilibrate("))?.Text ?? "";
            Check(eqLine.Contains("\"nvt 0.04 ps T 300 # hold\"") && eqLine.Contains("dt="), $"macro records equilibration: {eqLine}");
            // Target: the open structure — the script takes it with caps.current() and hands a result back, which opens
            vm.NewMacro();
            vm.MacroText = "import caps\n\ndoc = caps.current()\nprint(\"atoms\", doc.atoms)\ncaps.hand_back(doc)\n";
            vm.MacroTarget = 1;
            var macroAtoms = vm.Document!.Summary().Atoms;
            vm.RunMacro().GetAwaiter().GetResult();
            for (var i = 0; i < 20; i++) { Avalonia.Threading.Dispatcher.UIThread.RunJobs(); Thread.Sleep(25); }
            Check(vm.MacroOutput.Contains($"atoms {macroAtoms}") && vm.MacroOutput.Contains("the result is open") && vm.Title.Contains("result"),
                  $"macro on the open structure: {vm.MacroOutput.Replace('\n', ' ').Trim()} · {vm.Title}");
            vm.MacroTarget = 0;
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
            // electron diffraction: Peng 1996 scattering factors, a 200 kV source for the angles
            var xrayPeak = vm.ScatterPeaks.FirstOrDefault()?.Key ?? "";
            vm.FormFactorMode = 1;
            vm.RunScattering().GetAwaiter().GetResult();
            var electron = vm.Analyze.Results.FirstOrDefault(r => r.Id == "electron");
            var electronCurve = vm.ScatterXrayCurve.Length;
            vm.FormFactorMode = 0;
            vm.RunScattering().GetAwaiter().GetResult();   // the X-ray curve back for the overlay below
            Check(electron != null && electronCurve > 10 && vm.ScatterXrayCurve.Length > 10,
                  $"electron scattering: {electronCurve} pts (Peng 1996 factors, TEM 200 kV) · X-ray back {vm.ScatterXrayCurve.Length} pts{(xrayPeak.Length > 0 ? " · X-ray peak " + xrayPeak : "")}");
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
            // the same slit with TraPPE-UA methane: one site per molecule; the walls UFF, the fluid TraPPE-UA, by group
            vm.PoreFluidIndex = MainViewModel.PoreFluids.ToList().FindIndex(f => f.UnitedAtom && f.Smiles == "C");
            vm.BuildNano().GetAwaiter().GetResult();
            var ua = vm.Document!.Summary().Atoms;
            var uffPore = vm.Field.Library.ToList().FindIndex(x => x.Id == "uff");
            var trPore = vm.Field.Library.ToList().FindIndex(x => x.Id == "trappe-ua");
            vm.Field.GroupMode = true;   // suggests groups of its own: replaced by these two
            vm.Field.Groups.Clear();
            vm.Field.Groups.Add(new FieldGroupRow(vm.Field) { Name = "walls", Molecules = "1", FfIndex = uffPore });
            vm.Field.Groups.Add(new FieldGroupRow(vm.Field) { Name = "methane", Molecules = "rest", FfIndex = trPore, Charges = 0 });
            vm.Field.AssignGroups().GetAwaiter().GetResult();
            Check(ua == 440 + 10 && uffPore >= 0 && trPore >= 0 && vm.Field.Complete && vm.NanoLog.Contains("united atoms"),
                  $"pore with TraPPE-UA methane: {ua} atoms · {vm.Field.ForceFieldName} · complete {vm.Field.Complete} · {vm.Field.Log}");
            vm.Field.Clear().GetAwaiter().GetResult();
            vm.Field.Groups.Clear();
            vm.Field.GroupMode = false;
            vm.HoldOn = false;
            vm.PoreFluidIndex = 0;
            vm.NanoKind = 1;
        }
        {
            // silica for rubber: a passivated quartz sphere, TESPT grafted on its silanols (four sulfurs each), undone with Undo
            vm.OpenNano();
            vm.NanoKind = 2;
            vm.ParticleCrystal = Math.Max(0, vm.Crystals.ToList().FindIndex(c => c.Name.Contains("quartz", StringComparison.OrdinalIgnoreCase)));
            vm.ParticleRadius = 10;
            vm.ParticlePassivate = true;
            vm.BuildNano().GetAwaiter().GetResult();
            var bare = vm.Document!.Summary().Atoms;
            {   // the same particle with its cut surface relaxed (UFF, core held)
                vm.ParticleRelax = true;
                vm.BuildNano().GetAwaiter().GetResult();
                var relaxed = vm.Document!.Summary().Atoms;
                Check(relaxed == bare && vm.NanoLog.Contains("surface relaxed with UFF"), $"particle surface relaxed: {vm.NanoLog.Split('\n').FirstOrDefault(l => l.Contains("relaxed"))}");
                vm.ParticleRelax = false;
                vm.BuildNano().GetAwaiter().GetResult();
            }
            vm.SilanePick = 0;
            vm.SilaneFractionD = 0.2m;
            vm.GraftSilane().GetAwaiter().GetResult();
            var grafted = vm.Document!.Summary().Atoms;
            var sulfur = Enumerable.Range(0, (int)grafted).Count(i => vm.Document.Atom(i).Element == 16);
            vm.UndoEdit(false);
            Check(grafted > bare && sulfur > 0 && sulfur % 4 == 0 && (grafted - bare) == sulfur / 4 * 63 && vm.Document!.Summary().Atoms == bare,
                  $"silane: TESPT on silica · {bare} → {grafted} atoms, {sulfur / 4} grafts · undone · {vm.Status}");
            // a copper particle: adaptive CNA finds its FCC core, surface atoms are other; centrosymmetry is large there
            vm.OpenNano();
            vm.NanoKind = 2;
            vm.ParticleCrystal = Math.Max(0, vm.Crystals.ToList().FindIndex(c => c.Name.Contains("copper", StringComparison.OrdinalIgnoreCase)));
            vm.ParticleRadius = 12;
            vm.ParticlePassivate = false;
            vm.ParticleThiolate = false;
            vm.BuildNano().GetAwaiter().GetResult();
            var cuAtoms = vm.Document!.Summary().Atoms;
            vm.OpenVisualize();
            vm.ClearPipeline();
            vm.AddStep("cna");
            string A(string k) => vm.PipeAttributes.FirstOrDefault(a => a.Key == k)?.Value ?? "0";
            var fcc = long.Parse(A("CommonNeighborAnalysis.counts.FCC"));
            var other = long.Parse(A("CommonNeighborAnalysis.counts.Other"));
            vm.ClearPipeline();
            vm.InspectorTab = 0;   // the step opened its table: back to the particles
            vm.SetModule(8);
            Check(fcc > cuAtoms / 2 && fcc + other == cuAtoms && other > 0, $"CNA on a copper particle: {fcc} FCC, {other} other of {cuAtoms}");
            // a gold particle capped with hexanethiolates; PTM finds the FCC core under the ligands
            vm.OpenNano();
            vm.NanoKind = 2;
            vm.ParticleCrystal = Math.Max(0, vm.Crystals.ToList().FindIndex(c => c.Id == "gold"));
            vm.ParticleRadius = 12;
            vm.ParticleThiolate = true;
            vm.ThiolatePick = 0;
            vm.BuildNano().GetAwaiter().GetResult();
            var auAll = vm.Document!.Summary().Atoms;
            var auS = Enumerable.Range(0, (int)auAll).Count(i => vm.Document.Atom(i).Element == 16);
            var auAu = Enumerable.Range(0, (int)auAll).Count(i => vm.Document.Atom(i).Element == 79);
            vm.ParticleThiolate = false;
            vm.OpenVisualize();
            vm.ClearPipeline();
            vm.AddStep("ptm");
            var ptmFcc = long.Parse(A("PolyhedralTemplateMatching.counts.FCC"));
            var ptmOther = long.Parse(A("PolyhedralTemplateMatching.counts.Other"));
            vm.ClearPipeline();
            vm.InspectorTab = 0;
            vm.SetModule(8);
            Check(auAu == 429 && auS > 30 && auAll == auAu + auS * 20 && ptmFcc > 150 && ptmFcc + ptmOther == auAll,
                  $"gold particle with thiolates: {auAu} Au, {auS} C6 thiolates · PTM {ptmFcc} FCC, {ptmOther} other of {auAll}");
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
            // the Fragments tab: the project drawer on its fragments tab; Monitors: the pinned measurements shown or hidden
            vm.Compact = true;
            vm.FragmentsDrawer = true;
            var frag = vm.ProjectPanelShown && vm.ShowFragmentsTab && vm.FragmentsDrawer && !vm.ProjectOnlyDrawer && !vm.InspectorShown;
            vm.MonitorsShown = false;
            var hidden = !vm.MonitorsVisible;
            vm.MonitorsShown = true;
            vm.Compact = false;
            Check(frag && hidden, $"compact side tabs: fragments drawer {frag} · monitors hidden {hidden}");
            vm.LeftTab = 0;
            vm.ProjectDrawer = false;
            vm.InspectorDrawer = true;
        }

        // Theory manual: every page's references resolve in CAPS's BibTeX table
        {
            vm.OpenManual("csvr");
            // this session's thermostatted runs are listed under "Used in this project"
            Check(vm.ManualUsedProject.Any(x => x.StartsWith("dynamics-") || x.StartsWith("equilibrate-")) && vm.ManualHasUsedProject,
                  $"manual used in the project: {string.Join(" | ", vm.ManualUsedProject.Take(3))}");
            var unresolved = new List<string>();
            foreach (var item in vm.ManualNav.Where(n => n.Page != null))
            {
                vm.ShowManualPage(item.Page);
                unresolved.AddRange(vm.ManualRefs.Where(r => !r.Contains(" (")));
                // every equation typesets (LaTeX parses), and every symbol row
                void Math(string what, string t, Views.MathView.MathMode m)
                {
                    try { Views.MathView.MathParser.Parse(t, m); } catch (Exception e) { unresolved.Add($"{item.Page!.Id} {what}: {e.Message}"); }
                }
                if (item.Page!.Tex.Length == 0) unresolved.Add(item.Page.Id + ": no LaTeX");
                Math("equation", item.Page.Tex, Views.MathView.MathMode.Tex);
                foreach (var r in vm.ManualSymbols)
                {
                    Math("symbol " + r.Key, r.Key, Views.MathView.MathMode.Symbol);
                    Math("meaning " + r.Key, r.Value, Views.MathView.MathMode.Prose);
                    Math("setting " + r.Key, r.Other, Views.MathView.MathMode.Prose);
                }
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
            // the results table with the literature row (atactic PS) as CSV
            vm.ProjectRefIndex = vm.Analyze.References.ToList().FindIndex(m => m.Id == "ps-atactic");
            var csvPath = Path.Combine(outDir, "caps-selftest-results.csv");
            vm.ExportProjectTable(csvPath);
            var csv = File.Exists(csvPath) ? File.ReadAllLines(csvPath) : [];
            Check(vm.HasProjectRef && csv.Length == 4 && csv[0].StartsWith("structure,atoms,density") && csv[3].Contains("literature") && csv[3].Contains("1.040–1.065"),
                  $"project table: {csv.Length} lines · {(csv.Length > 3 ? csv[3] : "")} · Tg {vm.ProjectRefTg}");
            vm.ProjectRefIndex = 0;
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
            // combining: 2 tacticities × 3 lengths × 1 seed is 6 runs as a grid, 4 one factor at a time, 2 paired
            vm.SweepDps = "3, 4, 5";
            vm.SweepSeeds = "1";
            var grid = vm.SweepSubtitle;
            vm.SweepCombine = 1;
            var ofat = vm.SweepSubtitle;
            vm.SweepCombine = 2;
            var paired = vm.SweepSubtitle;
            vm.SweepCombine = 0;
            Check(grid.StartsWith("6 runs") && ofat.StartsWith("4 runs") && paired.StartsWith("2 runs"), $"sweep combine: {grid} | {ofat} | {paired}");
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
            // charges kept through the reaction (a 'charges keep' line), and the save writes the reaction SMARTS and JSON beside it
            vm.OpenTemplateEditor("cc_crosslink");
            // drags in the pre-reaction pane: C1–H3 is bonded (break), C1–C2 already formed (the drag takes it back out)
            vm.ToggleTemplateBond(1, 3);
            var breakAdded = vm.TemplateText.Contains("break 1 3");   // (H3 is deleted: the drawing shows 1–3 gone either way)
            vm.ToggleTemplateBond(3, 1);
            vm.ToggleTemplateBond(2, 1);
            var formGone = !vm.TemplateText.Contains("break 1 3") && !vm.TemplateText.Contains("form 1 2") && vm.TemplateFormed.Count == 1;   // H3–H4 (the H2 byproduct) stays
            vm.ToggleTemplateBond(1, 2);
            Check(breakAdded && formGone && vm.TemplateText.Contains("form 1 2") && vm.TemplateFormed.Count == 2, $"template drags: break added {breakAdded} · form taken out {formGone}");
            vm.TemplateCharges = 1;
            var keptLine = vm.TemplateText.Contains("charges keep") && vm.TemplateCharges == 1;
            var tplSaved = vm.SaveTemplate();
            var smartsPath = Path.Combine(MainViewModel.TemplateFolder, "cc_crosslink.smarts");
            var smarts = File.Exists(smartsPath) ? File.ReadAllText(smartsPath) : "";
            vm.TemplateCharges = 0;
            Check(keptLine && !vm.TemplateText.Contains("charges keep") && smarts.StartsWith("[#6;X4;H2;A:1]~[#1:3].[#6;X4;H2;A:2]~[#1:4]>>[#6:1]-[#6:2]"),
                  $"template charges kept + SMARTS: kept {keptLine} · removed {!vm.TemplateText.Contains("charges keep")} · {smarts.Trim()} · {tplSaved} · {vm.TemplateError}");
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

        // a LAMMPS dump with a compute column and forces: the columns are offered for colouring and colour the atoms
        {
            var dumpPath = Path.Combine(outDir, "caps-selftest-columns.lammpstrj");
            var sb = new System.Text.StringBuilder();
            foreach (var step in new[] { 0, 100 })
            {
                sb.Append($"ITEM: TIMESTEP\n{step}\nITEM: NUMBER OF ATOMS\n3\nITEM: BOX BOUNDS pp pp pp\n0 10\n0 10\n0 10\nITEM: ATOMS id element x y z fx fy fz c_pe\n");
                sb.Append("1 C 1 1 1 3 4 0 -1\n2 C 2.5 1 1 0 0 0 -2\n3 C 4 1 1 0 0 4 -3\n");
            }
            File.WriteAllText(dumpPath, sb.ToString());
            var openBefore = vm.ActiveItem;
            vm.Open(dumpPath);
            var cols = vm.AppColumns.ToList();
            vm.AppColumn = cols.IndexOf("c_pe");
            var ok = vm.AppHasColumns && cols.Contains("c_pe") && cols.Contains("|f|") && vm.AppRampEnabled && vm.Status.IndexOf("fail", StringComparison.OrdinalIgnoreCase) < 0;
            vm.AppColumn = 0;
            Check(ok, $"dump columns: {string.Join(", ", cols)} · {vm.Status}");
            if (openBefore != null) vm.Activate(openBefore);   // the steps below work on the structure that was open
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
            // TIFF and PDF of the same figure: a baseline TIFF of the figure's pixels, a PDF page the figure's size in points
            {
                var tifPath = Path.Combine(outDir, "caps-selftest-figure.tiff");
                var pdfPath = Path.Combine(outDir, "caps-selftest-figure.pdf");
                vm.FigFormat = 2;
                vm.ExportFigure(tifPath, (_, _, _, _, _, _) => { }, (_, _) => { }, (px, _, _, _) => px).GetAwaiter().GetResult();   // no drawing platform here: the pixels as rendered
                vm.FigFormat = 3;
                vm.ExportFigure(pdfPath, (_, _, _, _, _, _) => { }, (_, _) => { }, (px, _, _, _) => px).GetAwaiter().GetResult();   // no drawing platform here: the pixels as rendered
                var tif = File.ReadAllBytes(tifPath);
                var tifOk = tif.Length > 3 * 2008 * 1130 && tif[0] == (byte)'I' && tif[1] == (byte)'I' && tif[2] == 42 && BitConverter.ToInt32(tif, 8 + 2 + 8) == 2008;
                var pdfText = System.Text.Encoding.ASCII.GetString(File.ReadAllBytes(pdfPath));
                var inches = 2008.0 / (double)vm.FigDpi;
                var widthPt = (inches * 72).ToString("0.###", System.Globalization.CultureInfo.InvariantCulture);
                var pdfOk = pdfText.StartsWith("%PDF") && pdfText.Contains("/MediaBox [0 0 " + widthPt + " ");
                Check(tifOk && pdfOk, $"figure TIFF {tif.Length} bytes ok {tifOk} · PDF ok {pdfOk} ({inches * 72:0} pt wide)");
            }
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
            // the Python overlay: the template script draws its text on the rendered image
            var ovPng = Path.Combine(outDir, "caps-selftest-overlay.png");
            var script = vm.NewOverlayScript();
            System.Text.Json.Nodes.JsonArray? drawn = null;
            vm.RenderOut(_ => ovPng, false, (rgba, w, h, sp, path) => { drawn = sp.Custom; }).GetAwaiter().GetResult();
            Check(File.Exists(script) && drawn is { Count: > 0 } && drawn.Any(c => ((string?)c?["s"] ?? "").Contains("ρ =")),
                  $"python overlay: {drawn?.Count ?? 0} commands · {vm.OvPythonNote} {vm.OvConsole}");
            vm.OvPython = false;
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
            // wrap by molecule: the editor offers the modes with a note; no bond is left across a face
            vm.OpenVisualize();
            vm.ClearPipeline();
            vm.AddStep("wrap");
            var modeField = vm.StepFields.FirstOrDefault(f => f.Key == "mode");
            if (modeField != null) modeField.Text = "molecules";
            var across = vm.PipeAttributes.FirstOrDefault(a => a.Key == "Wrap.bonds_across_faces")?.Value;
            var moved = vm.PipeAttributes.FirstOrDefault(a => a.Key == "Wrap.molecules_moved")?.Value;
            Check(modeField != null && vm.StepFields.Any(f => f.IsNote) && across == "0" && moved != null,
                  $"wrap molecules: {across} bonds across faces · {moved} molecules moved · {vm.PipelineRows.LastOrDefault()?.Summary}");
            // a Python step typed in the editor: the template runs; an edit runs only when sent
            vm.ClearPipeline();
            vm.AddStep("python");
            var codeField = vm.StepFields.FirstOrDefault(f => f.IsCode);
            var heavy = vm.PipeAttributes.FirstOrDefault(a => a.Key == "Heavy atoms")?.Value;
            if (codeField != null) codeField.Draft = codeField.Draft.Replace("sum(heavy)", "2 * sum(heavy)");
            var beforeRun = vm.PipeAttributes.FirstOrDefault(a => a.Key == "Heavy atoms")?.Value;
            var pending = codeField?.DraftChanged == true;
            codeField?.Commit();
            var afterRun = vm.PipeAttributes.FirstOrDefault(a => a.Key == "Heavy atoms")?.Value;
            var console = vm.StepFields.FirstOrDefault(f => f.IsCode)?.Output ?? "";
            Check(heavy == "640" && beforeRun == "640" && pending && afterRun == "1280" && codeField?.DraftChanged == false && console.Contains("1280 heavy atoms in 10 molecules"),
                  $"python step typed: {heavy} → {afterRun} after Run (pending {pending}) · console '{console}' · {vm.PipelineRows.LastOrDefault()?.Summary}");
            vm.ClearPipeline();
            // chain orientation per atom, then an affine strain of 10 % along x
            vm.AddStep("orientation");
            var orS = vm.PipeAttributes.FirstOrDefault(a => a.Key == "Orientation.S")?.Value;
            var hasProps = vm.PipeProperties().Contains("Orientation") && vm.PipeProperties().Contains("Crystalline");
            vm.ClearPipeline();
            vm.AddStep("affine_transform");
            var ratio = vm.PipeAttributes.FirstOrDefault(a => a.Key == "AffineTransformation.volume_ratio")?.Value;
            // freeze: frame-0 heights on the shown frame
            vm.ClearPipeline();
            vm.AddStep("freeze_property");
            var frozen = vm.PipeProperties().Contains("Position.Z frozen") && vm.PipelineRows.Last().Level == "ok";
            Check(frozen, $"freeze property: {vm.PipelineRows.LastOrDefault()?.Summary}");
            Check(orS != null && double.Parse(orS, System.Globalization.CultureInfo.InvariantCulture) is > 0 and < 0.5 && hasProps && ratio == "1.1",
                  $"orientation step: S {orS}, properties {hasProps} · affine: volume × {ratio}");
            vm.ClearPipeline();
            vm.SetModule(8);
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
            // on error: retried once then skipped; with "stop" the batch ends at the failure (the broken input first)
            {
                vm.BatchOnError = 1;
                vm.RunBatch(outBatch).GetAwaiter().GetResult();
                var badRow = vm.BatchInputs.First(b => b.Path == bad);
                var retried = badRow.State == "failed" && badRow.Error.EndsWith("(failed twice)") && vm.BatchInputs.Count(b => b.State == "done") == 2;
                vm.BatchOnError = 2;
                vm.BatchInputs.Move(vm.BatchInputs.IndexOf(badRow), 0);
                vm.RunBatch(outBatch).GetAwaiter().GetResult();
                var stopped = badRow.State == "failed" && vm.BatchInputs.Count(b => b.State == "done") == 0 && vm.BatchState.Contains("stopped");
                vm.BatchOnError = 0;
                Check(retried && stopped, $"batch on error: retried {retried} ({badRow.Error}) · stop {stopped} ({vm.BatchState})");
            }
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
            // a pipeline file with a Python step: loaded, but the step is switched off until turned on
            var pyPath = Path.Combine(outDir, "caps-selftest-python.json");
            File.WriteAllText(pyPath, "{\"steps\":[{\"type\":\"python\",\"code\":\"print('hi')\"},{\"type\":\"coordination\"}]}");
            vm.ClearPipeline();
            vm.LoadPipeline(pyPath);
            Check(vm.PipelineRows.Count == 2 && !vm.PipelineRows[0].Enabled && vm.PipelineRows[1].Enabled && vm.PipelineHeldPython == 1 && vm.Status.Contains("switched off"),
                  $"load pipeline with Python: held {vm.PipelineHeldPython} · {vm.Status}");
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
                // Tg from the potential energy per atom, two lines through the glassy and rubbery ranges, Berendsen at 2 atm
                vm.Analyze.TgProperty = 1;
                vm.Analyze.TgFit = 1;
                vm.Analyze.TgBarostat = 1;
                vm.Analyze.TgPressureD = 2;
                vm.RunGlass().GetAwaiter().GetResult();
                var tgCard = vm.Analyze.Results.FirstOrDefault(c => c.Id == "tg");
                // 2 ps holds are far too short for a break: the pooled result is a Tg or, when the lines cross outside the ranges, none
                var energyRun = vm.GtPoints.Length == 5 && vm.GtYLabel.StartsWith("potential energy") &&
                                (vm.GtTg.Value == "—" || vm.GtAlphaGlass.Label == "dE/dT glass") && tgCard != null && tgCard.Method.Contains("potential energy per atom") && tgCard.Method.Contains("glassy and the rubbery range") &&
                                tgCard.Method.Contains("2 atm") && tgCard.Method.Contains("Berendsen", StringComparison.OrdinalIgnoreCase);
                Check(energyRun, $"glass transition from the energy, two ranges: {vm.GtPoints.Length} points · {vm.GtTg.Value} · {tgCard?.Method}");
                // the saved recipe: the force field as assigned here (its file beside the recipe) and every cooling setting
                {
                    var gdir = System.IO.Path.Combine(System.IO.Path.GetTempPath(), "caps_glass_recipe");
                    System.IO.Directory.CreateDirectory(gdir);
                    var noFf = vm.Field.Assigned ? null : vm.SaveGlassForceField(System.IO.Path.Combine(gdir, "none.yaml"));
                    if (!vm.Field.Assigned) vm.Field.Assign().GetAwaiter().GetResult();
                    var ffName = vm.SaveGlassForceField(System.IO.Path.Combine(gdir, "x_tg.yaml"));
                    var grec = vm.GlassRecipe("/tmp/x.data", ffName);
                    Check(noFf == null && vm.GlassRecipe("/tmp/y.data").Contains("type: { forcefield: default }") && ffName == "x_tg.ff.json" && System.IO.File.ReadAllText(System.IO.Path.Combine(gdir, ffName)).Contains("caps-assigned-forcefield") &&
                          grec.Contains("type: { file: \"x_tg.ff.json\" }") && grec.Contains("barostat: berendsen") && grec.Contains("pressure: 2") &&
                          grec.Contains("property: energy") && grec.Contains("fit: ranges") && grec.Contains("glassy_max:"),
                          "glass recipe: " + grec.Split('\n').Last(l => l.Contains("tg:")).Trim());
                    // replicas run elsewhere (here: the recipe engine into job folders as a host would leave them) pooled back
                    var gjobs = new List<Job>();
                    for (var rep = 1; rep <= 3; ++rep)
                    {
                        var rd = System.IO.Path.Combine(gdir, $"rep{rep}");
                        System.IO.Directory.CreateDirectory(rd);
                        var job = new Job { Id = $"glass-t{rep}", Kind = "Glass", Module = 47, Title = $"replica {rep}", Document = "test", Status = "done",
                                            Remote = new RemoteRun { Host = "testhost", Local = rd, Stem = $"x_tg_r{rep}", Batch = "selftest-glass" } };
                        gjobs.Add(job);
                        if (rep == 3) continue;   // one that never came back
                        vm.Document!.Save(System.IO.Path.Combine(rd, "structure.caps.data"));
                        System.IO.File.Copy(System.IO.Path.Combine(gdir, ffName!), System.IO.Path.Combine(rd, "structure.ff.json"), true);
                        var (_, glassRep) = CapsDocument.RunRecipe(vm.GlassRecipe("structure.caps.data", "structure.ff.json", (ulong)rep, $"x_tg_r{rep}"),
                                                                $"{{\"base_dir\": \"{rd}\", \"out_dir\": \"{System.IO.Path.Combine(rd, "out")}\"}}", "replica", null);
                        _ = glassRep;
                    }
                    foreach (var job in gjobs) vm.Jobs.Add(job);
                    vm.CollectGlassReplicas("selftest-glass");
                    Check(vm.GtPoints.Length == 5 && vm.GtStatus.StartsWith("2 replicas") && vm.GtStatus.Contains("x_tg_r3") && vm.GtYLabel.StartsWith("potential energy") &&
                          (vm.GtTg.Value == "—" || vm.GtTg.Caption.Contains("two ranges")),
                          $"glass replicas pooled from job folders: {vm.GtPoints.Length} points · {vm.GtTg.Value} ({vm.GtTg.Caption}) · {vm.GtStatus}");
                    foreach (var job in gjobs) vm.Jobs.Remove(job);
                }
                vm.Analyze.TgProperty = 0; vm.Analyze.TgFit = 0; vm.Analyze.TgBarostat = 0; vm.Analyze.TgPressureD = 1;
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
                Check(templates == 5 && first != null && valid && flagged && files == 2 && vm.Document?.Summary().Atoms == 9 && prov.Contains("recipe.run") && prov.Contains(vm.RecipeSha),
                      $"recipes: {templates} · valid {valid} · flagged {flagged} · {files} files · {vm.Status}");
                vm.SelectedRecipe = vm.Recipes.First(r => r.Name.StartsWith("Sulfur-cured"));
                var cureOk = vm.RecipeOk && vm.RecipeStages.Any(st => st.Name == "React" && st.Ok);
                Check(cureOk, $"recipe template: sulfur-cured NR · {string.Join(", ", vm.RecipeStages.Select(st => st.Name + (st.Ok ? "" : "!")))} · {vm.RecipeValid}");

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
                psi.Environment["GMX_MAXBACKUP"] = "-1";   // repeated self-tests: no #system.tpr.N# backups (GROMACS stops at 99)
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
            // Nosé–Hoover NPT: MTK chooses the Nosé–Hoover thermostat; the decks name fix npt and Nose-Hoover/MTTK
            vm.MdEnsemble = 2;
            vm.MdBarostat = 2;
            vm.PreflightNow().GetAwaiter().GetResult();
            var nhDeck = vm.MdDeck;
            vm.MdGromacs = true;
            vm.PreflightNow().GetAwaiter().GetResult();
            var nhGmx = vm.MdDeck;
            vm.MdGromacs = false;
            Check(vm.MdThermostat == 2 && nhDeck.Contains("fix 1 all npt temp 300 300 100 iso 1 1 1000") && nhGmx.Contains("Nose-Hoover") && nhGmx.Contains("MTTK"),
                  $"Nosé–Hoover NPT: thermostat {vm.MdThermostat}, LAMMPS {nhDeck.Split('\n').FirstOrDefault(l => l.StartsWith("fix 1"))}, GROMACS MTTK {nhGmx.Contains("MTTK")}");
            vm.MdBarostat = 0;
            vm.MdThermostat = 0;
            vm.MdEnsemble = 1;
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
            for (int k = 0; k < 400 && vm.Field.Working; ++k) Thread.Sleep(25);
            // chains and toluene each a component: a group each with Grow's force field (their own types, grouped in LAMMPS)
            var gDir = Path.Combine(outDir, "grow-groups");
            Directory.CreateDirectory(gDir);
            vm.ExportFieldLammps(Path.Combine(gDir, "grown.data")).GetAwaiter().GetResult();
            var gIn = File.Exists(Path.Combine(gDir, "grown.in")) ? File.ReadAllText(Path.Combine(gDir, "grown.in")) : "";
            var gLines = gIn.Split('\n').Where(l => l.StartsWith("group ", StringComparison.Ordinal)).ToList();
            var growGroups = string.Join(" | ", gLines) + $" · assign {vm.GrowAssignField} · grouped {vm.Field.IsGrouped} · assigned {vm.Field.Assigned} {vm.Field.ForceFieldName} · ff {vm.GrowFfIndex} {vm.Field.Library[vm.GrowFfIndex].Id} · log {vm.Field.Log}";
            // a united-atom force field cannot type the components apart: the whole cell, assigned
            var ua = vm.Field.Library[vm.GrowFfIndex].Id.Contains("-ua", StringComparison.Ordinal);
            var grpOk = !vm.GrowAssignField || (ua ? vm.Field.Assigned && !vm.Field.IsGrouped
                                                   : vm.Field.IsGrouped && gLines.Count == 2 && gLines[0].Contains("chains") && gLines[1].Contains("Toluene"));
            // with an all-atom force field: by component
            if (ua)
            {
                vm.GrowFfIndex = vm.Field.Library.ToList().FindIndex(e => e.Id == "opls2005");
                vm.Grow().GetAwaiter().GetResult();
                for (int k = 0; k < 400 && vm.Field.Working; ++k) Thread.Sleep(25);
                vm.ExportFieldLammps(Path.Combine(gDir, "grown.data")).GetAwaiter().GetResult();
                gLines = File.ReadAllText(Path.Combine(gDir, "grown.in")).Split('\n').Where(l => l.StartsWith("group ", StringComparison.Ordinal)).ToList();
                grpOk &= vm.Field.IsGrouped && gLines.Count == 2 && gLines[0].Contains("chains") && gLines[1].Contains("Toluene");
                growGroups += " · OPLS: " + string.Join(" | ", gLines);
            }
            Check(mols == 3 + 5 && vm.GrowLog.Contains("5 × Toluene") && grpOk, $"grow with solvent: {mols} molecules · groups {growGroups} · {vm.GrowLog.Split('\n').FirstOrDefault(l => l.Contains("Toluene"))}");
            vm.Field.Clear().GetAwaiter().GetResult();
            vm.Field.Groups.Clear();
            vm.Field.GroupMode = false;
            vm.RemoveGrowSmall(vm.GrowSmall[0]);
            // oriented growth: chains drawn along z, the report gives their ⟨P₂⟩
            vm.GrowOrient = 3;
            vm.GrowOrientStrengthD = 4;
            vm.Grow().GetAwaiter().GetResult();
            var p2line = vm.GrowLog.Split('\n').FirstOrDefault(l => l.Contains("oriented growth")) ?? "";
            var p2 = double.TryParse(p2line.Split("⟨P₂⟩ = ").LastOrDefault()?.Split(" ")[0], System.Globalization.NumberStyles.Float, System.Globalization.CultureInfo.InvariantCulture, out var pv) ? pv : double.NaN;
            Check(p2 > 0.3 && vm.GrowRecipe().Contains("orientation: { axis: z"), $"oriented growth: {p2line}");
            vm.GrowOrient = 0;
        }

        // The pipeline in order (design/boards/PipelineGrow, ExportCenter): the force field chosen in Grow is assigned when
        // growing finishes, the strip says where the project is, and the Export center writes LAMMPS and GROMACS files
        {
            vm.UsePolystyreneInGrow();
            vm.GrowChainsD = 3;
            vm.GrowDpD = 6;
            vm.GrowDensityD = 0.3m;
            vm.GrowFfIndex = vm.Field.Library.ToList().FindIndex(e => e.Id == "gaff-amber25");   // Grow's own choice
            vm.GrowChargeMode = 0;
            vm.GrowAssignField = true;
            vm.Grow().GetAwaiter().GetResult();
            var steps = string.Join(" · ", vm.PipelineSteps.Select(s => $"{s.Name} {s.State}"));
            // the strip states what the active structure is and what was done to it, in no imposed order
            Check(vm.Field.Assigned && vm.Field.Complete && vm.PipelineSteps.Count == 7 && vm.PipelineSteps[1].Detail == "Polymer cell"
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
            // DL_POLY as well: FIELD, CONFIG, CONTROL in STEM_dlpoly
            vm.EngineDlpoly = true;
            vm.RefreshEnginesNow();
            var dlListed = vm.EngineGromacsFiles.Any(f => f.Name.EndsWith("_dlpoly/FIELD"));
            vm.WriteEngines().GetAwaiter().GetResult();
            var field = Path.Combine(pkg, "system_dlpoly", "FIELD");
            var fieldText = File.Exists(field) ? File.ReadAllText(field) : "";
            Check(dlListed && fieldText.Contains("units kcal") && fieldText.Contains("nummols") && fieldText.Contains("vdw") && File.Exists(Path.Combine(pkg, "system_dlpoly", "CONFIG")),
                  $"export center DL_POLY: listed {dlListed}, FIELD {fieldText.Length} bytes");
            vm.EngineDlpoly = false;
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
            // Analyze › Adsorption locator: two waters annealed on the PE cell (the built-in force field); the adsorbates are
            // added after the substrate, the lowest configuration shown, the energy negative, dE/dN for the component
            {
                var atoms0 = vm.Document!.Summary().Atoms;
                vm.OpenAdsorption();
                vm.AdsRows.Clear();
                vm.AdsRows.Add(new AdsorbateRow { Smiles = "O", Count = 2 });
                vm.AdsRegion = 0;
                vm.AdsCyclesD = 1;
                vm.AdsStepsD = 600;
                vm.RunAdsorption().GetAwaiter().GetResult();
                var added = vm.Document!.Summary().Atoms - atoms0;
                Check(vm.IsAdsorption && !vm.AdsHasError && added == 6 && vm.AdsEnergy.EndsWith("kcal/mol") && !vm.AdsEnergy.StartsWith("0") && vm.AdsComponents.Count == 1
                      && vm.AdsComponents[0].Name == "H2O" && vm.AdsConfigs.Count >= 1 && vm.Document!.Provenance().Contains("adsorption.locator"),
                      $"adsorption locator: {vm.AdsEnergy} · {vm.AdsSplit} · {added} atoms added · {vm.AdsError}");
                vm.UndoEdit(false);
                // Analyze › Sorption: methane in the PE cell (Widom, then two GCMC pressures); the structure is not changed
                vm.OpenSorption();
                vm.Sorbate = "C";
                vm.SorbInsertD = 20000;
                vm.SorbPressures = "100, 1000";
                vm.SorbStepsD = 4000;
                var atomsS = vm.Document!.Summary().Atoms;
                vm.RunSorption().GetAwaiter().GetResult();
                Check(vm.IsSorption && !vm.SorbHasError && vm.SorbS.Contains("cm³") && vm.SorbRows.Count == 2 && vm.Document!.Summary().Atoms == atomsS
                      && vm.Document!.Provenance().Contains("sorption.widom_gcmc"),
                      $"sorption: S {vm.SorbS} · K_H {vm.SorbHenry} · μex {vm.SorbMu} · {vm.SorbRows.Count} points · {vm.SorbError}");
                vm.SetModule(4);
            }
            // Equilibrate › Chain ends: CBMC regrowth with the built-in force field; a new frame, the provenance step
            vm.CbMovesD = 60;
            var framesBefore = vm.Frames;
            vm.RunCbmc().GetAwaiter().GetResult();
            Check(vm.CbNote.Contains("regrowths accepted") && vm.Frames > framesBefore && vm.Document!.Provenance().Contains("cbmc.regrow"),
                  $"CBMC: {vm.CbNote} · frames {framesBefore} → {vm.Frames}");
            vm.EqTarget = 1;
            var risTarget = vm.EqTargetNote.Contains("Flory");
            var tf = Path.Combine(outDir, "target.dat");
            File.WriteAllText(tf, "# n  C(n)\n1 1.0\n2 1.9\n3, 2.6\nbad line\n10 4.8\n");
            vm.LoadEqTarget(tf);
            var fileOk = vm.EqTarget == 2 && vm.ChainReference.Length == 4 && vm.ChainReference[^1] == (10.0, 4.8);
            vm.EqTarget = 1;
            vm.Open(Path.Combine(dir, "ps_melt.data"));
            Check(risOk && vm.RisCurve.Length == 0 && risTarget && fileOk && vm.EqTargetNote.Contains("no target"),
                  $"RIS reference: {vm.RisCurve.Length} on PS · alkane note ok {risOk} · RIS target {risTarget} · file target {fileOk} · {vm.EqTargetNote}");
            vm.EqTarget = 0;
        }

        // Field › Fill from OPLS 2005: natural rubber under OPLS-AA 2024 lacks the CM-CT-CT-CM torsion; one click borrows
        // it (and only what is missing) from OPLS 2005 by the shared OPLS classes
        {
            var (nrDoc, _) = CapsDocument.GrowChains("{\"units\":[{\"name\":\"isoprene\",\"smiles\":\"[*]C/C(C)=C\\\\C[*]\"}],\"dp\":4}",
                new CapsGrowOpts { Chains = 1, Dp = 0, Seed = 1, Density = 0.02, ContactScale = -0.8, Curve = 1 }, null, "nr");
            var nrFile = Path.Combine(outDir, "nr_fill.data");
            nrDoc.Save(nrFile);
            nrDoc.Dispose();
            vm.Open(nrFile);
            vm.Field.FfIndex = vm.Field.Library.ToList().FindIndex(x => x.Id == "oplsaa2024-moltemplate");
            vm.Field.Assign().GetAwaiter().GetResult();
            var fillBefore = (vm.Field.Complete, vm.Field.CanFillSuggested, vm.Field.FillSuggestedText);
            vm.Field.FillSuggested().GetAwaiter().GetResult();
            Check(!fillBefore.Complete && fillBefore.CanFillSuggested && fillBefore.FillSuggestedText == "Fill from OPLS 2005" && vm.Field.Complete && !vm.Field.CanFillSuggested,
                  $"fill from OPLS 2005: before complete {fillBefore.Complete}, offered '{fillBefore.FillSuggestedText}' · after complete {vm.Field.Complete} · {vm.Field.EstimatedText}");
            vm.Open(Path.Combine(dir, "ps_melt.data"));
        }

        // Field › Use OPLS 2005's charges: poly(2-vinylpyridine) under OPLS-AA 2024 — the substituted pyridine's fixed
        // charges do not balance; OPLS 2005's bond-increment charges do, with OPLS-AA 2024's types kept
        {
            var (pvp, _) = CapsDocument.GrowChains("{\"units\":[{\"name\":\"2-vinylpyridine\",\"smiles\":\"[*]CC([*])c1ccccn1\"}],\"dp\":4}",
                new CapsGrowOpts { Chains = 1, Dp = 0, Seed = 1, Density = 0.02, ContactScale = -0.8, Curve = 1 }, null, "p2vp");
            var pvpFile = Path.Combine(outDir, "p2vp_charges.data");
            pvp.Save(pvpFile);
            pvp.Dispose();
            vm.Open(pvpFile);
            vm.Field.ChargeMode = 1;
            vm.Field.FfIndex = vm.Field.Library.ToList().FindIndex(x => x.Id == "oplsaa2024-moltemplate");
            vm.Field.Assign().GetAwaiter().GetResult();
            var chBefore = (vm.Field.Complete, vm.Field.CanUseCompanionCharges, vm.Field.CompanionChargesText);
            vm.Field.UseCompanionCharges().GetAwaiter().GetResult();
            var pvpRep = vm.Document!.FieldReport();
            Check(!chBefore.Complete && chBefore.CanUseCompanionCharges && vm.Field.Complete && pvpRep.Contains("OPLS 2005's own, from its bond increments") && pvpRep.Contains("\"520_"),
                  $"OPLS 2005 charges with OPLS-AA 2024 types: before complete {chBefore.Complete}, offered '{chBefore.CompanionChargesText}' · after complete {vm.Field.Complete}");
            vm.Field.ChargeMode = 0;
            vm.Open(Path.Combine(dir, "ps_melt.data"));
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

        {
            // Pack on its own: 40 waters packed around the polystyrene melt, which stays fixed in its periodic cell; Pack's
            // own force field is assigned; the next step offered is the export (minimise / dynamics optional)
            vm.Open(Path.Combine(dir, "ps_melt.data"));
            vm.SetModule(5);
            vm.NewPackInput();
            vm.PackStart = 1;
            vm.PackCountD = 40;
            vm.PackAssignField = false;
            vm.AddPackStructure(Path.Combine(dir, "water.pdb"));
            vm.RunPack().GetAwaiter().GetResult();
            var packCell = vm.Document!.Summary();
            var export = vm.PipelineSteps.Any(p => p.Name == "Export");
            var optional = vm.PipelineSteps.First(p => p.Name == "Minimise").Detail == "optional";
            Check(packCell.Molecules == 50 && packCell.Atoms == 1300 + 120 && vm.PackDone && export && optional,
                  $"pack around the current structure: {packCell.Molecules} molecules, {packCell.Atoms} atoms · next: export (minimise {(optional ? "optional" : "?")}) · {vm.PackLog.Split('\n')[0]}");
            // the job tree: the Pack job's folder sits under the cell it made, with its result, report and chart
            var packJob = vm.ActiveItem?.Jobs.FirstOrDefault(j => j.Kind == "Pack");
            Check(packJob != null && packJob.Outputs.Any(o => o.Label.StartsWith("Result")) && packJob.Outputs.Any(o => o.Label.StartsWith("Report")),
                  $"job tree: {vm.ActiveItem?.Name} › {packJob?.FolderLabel} › {string.Join(" · ", packJob?.Outputs.Select(o => o.Label) ?? [])}");
            vm.PackStart = 0;
        }
        {
            // a dihedral restraint: a backbone torsion held at +60° (gauche; the sign as the measurement gives it) while the cell relaxes, measured the same way
            using var dr = CapsDocument.Open(Path.Combine(dir, "ps_melt.data"));
            int[] Heavy(int a) => System.Text.Json.Nodes.JsonNode.Parse(dr.AtomProperties(a))!["neighbours"]!.AsArray()
                .Where(n => (string?)n!["element"] == "C").Select(n => (int)n!["index"]!.GetValue<double>()).ToArray();
            int[]? quad = null;
            for (var j = 0; j < 200 && quad == null; j++)
                foreach (var k in Heavy(j))
                {
                    var i = Heavy(j).FirstOrDefault(x => x != k, -1);
                    var l = Heavy(k).FirstOrDefault(x => x != j && x != i, -1);
                    if (i >= 0 && l >= 0) { quad = [i, j, k, l]; break; }
                }
            var phiBefore = dr.Measure(quad!);
            dr.SetRestraints($"[{{\"i\":{quad![0]},\"j\":{quad[1]},\"k\":{quad[2]},\"l\":{quad[3]},\"phi0\":60,\"kphi\":200}}]");
            dr.Relax(new CapsRelaxOpts { Method = 2, Ftol = 1, MaxIterations = 2000, Pushoff = 1, Cutoff = 10, Coulomb = 1 }, null);
            var after = dr.Measure(quad);
            Check(Math.Abs(after - 60) < 6, $"dihedral restraint: torsion {string.Join("-", quad.Select(x => x + 1))} {phiBefore:F1}° → {after:F1}° (target +60°, the measured sign)");
        }
        {
            // the Properties explorer: the structure in numbers, then a picked atom with its bonded neighbours
            vm.Open(Path.Combine(dir, "ps_melt.data"));
            vm.RefreshProperties();
            var formula = vm.PropertyRows.FirstOrDefault(r => r.Name == "Formula")?.Value ?? "";
            var lattice = vm.PropertyRows.FirstOrDefault(r => r.Name.StartsWith("α"))?.Value ?? "";
            vm.PropertyFilter = "density";
            var filtered = vm.PropertyRows.Count(r => r.IsValue);
            vm.PropertyFilter = "";
            vm.Pick(0);
            vm.RefreshProperties();
            var bonded = vm.PropertyRows.FirstOrDefault(r => r.IsHeader && r.Name.StartsWith("Bonded to"))?.Name ?? "";
            Check(formula == "C640H660" && lattice == "90.00  90.00  90.00" && filtered == 1 && bonded.Length > 0 && vm.PropertyTitle.StartsWith("Atom 1"),
                  $"properties: {formula} · angles {lattice} · filter density → {filtered} row · {vm.PropertyTitle} · {bonded}");
        }
        {
            // a sketch with two attachment points becomes unit A of the polymer builder
            vm.MolSmiles = "*CC(*)C(=O)OC";
            var unitOk = vm.MolIsRepeatUnit;
            vm.UseSketchAsRepeatUnit();
            Check(unitOk && vm.PolyUnits.Count > 0 && vm.PolyUnits[0].Smiles == "*CC(*)C(=O)OC", $"sketch → repeat unit: {vm.PolyUnits.FirstOrDefault()?.Smiles} · {vm.Status}");
            vm.MolSmiles = "";
            // a repeat unit picked in 3D: ethylbenzene's CH3 (head) and CH2 (tail) → styrene's unit, C8H8
            vm.MolSmiles = "CCc1ccccc1";
            vm.BuildMolecule().GetAwaiter().GetResult();
            vm.OpenMoleculeInStudio();
            vm.Pick(0, false); vm.Pick(1, true);
            vm.UsePickedAsRepeatUnit();
            var picked = vm.PolyUnits.FirstOrDefault()?.Smiles ?? "";
            Check(picked.StartsWith("*") && picked.Count(ch => ch == '*') == 2 && vm.PolyUnits[0].Info.StartsWith("C₈H₈"), $"picked → repeat unit: {picked} · {vm.PolyUnits.FirstOrDefault()?.Info}");
            vm.MolSmiles = "";
        }
        {
            // Open in notebook: the notebook's code (all but the view) runs as written with the shipped caps package
            vm.Open(Path.Combine(dir, "ps_melt.data"));
            var nbPath = vm.WriteNotebook()!;
            var nbj = System.Text.Json.Nodes.JsonNode.Parse(File.ReadAllText(nbPath))!;
            var code = nbj["cells"]!.AsArray().Where(c => (string?)c!["cell_type"] == "code").Select(c => (string)c!["source"]!).ToList();
            var script = Path.Combine(outDir, "caps-notebook-check.py");
            File.WriteAllText(script, string.Join("\n\n", code.Where(c => !c.Contains("doc.view"))) + "\nprint('NOTEBOOK OK', doc.summary()['atoms'])\n");
            var psi = new System.Diagnostics.ProcessStartInfo("python3", $"\"{script}\"") { RedirectStandardOutput = true, RedirectStandardError = true };
            var run = System.Diagnostics.Process.Start(psi)!;
            var stdout = run.StandardOutput.ReadToEnd();
            var stderr = run.StandardError.ReadToEnd();
            run.WaitForExit(120000);
            Check(code.Count == 4 && run.ExitCode == 0 && stdout.Contains("NOTEBOOK OK 1300"),
                  $"notebook: {Path.GetFileName(nbPath)} · {code.Count} code cells · runs: {(run.ExitCode == 0 ? "yes" : stderr.Split('\n').LastOrDefault(l => l.Length > 0))}");
            try { File.Delete(nbPath); } catch { }
        }
        {
            // the keyboard map: ] grows the picked atom's selection one bond (a CH carbon and its four neighbours), ⌘I inverts, ⌘8 is Dynamics
            vm.SetModule(8);
            vm.Pick(0);
            vm.GrowSelectionKey();
            var grownSel = vm.SelectedCount;
            vm.InvertSelectionKey();
            var inverted = vm.SelectedCount;
            vm.ClearDocSelection();
            vm.GoRailPage(8);
            Check(grownSel == 5 && inverted == 1300 - 5 && vm.IsDynamics, $"keys: ] selects {grownSel}, ⌘I {inverted}, ⌘8 Dynamics {vm.IsDynamics}");
            vm.SetModule(8);
        }

        // Your shortcuts: a key for a palette command runs it; a key the Studio uses is flagged; Reset clears them
        {
            var cmdMod = OperatingSystem.IsMacOS() ? Avalonia.Input.KeyModifiers.Meta : Avalonia.Input.KeyModifiers.Control;
            vm.ShortcutFilter = "perspective orthographic";
            var row = vm.ShortcutRows.FirstOrDefault(r => r.Id == "view.projection");
            var persp = vm.Perspective;
            if (row != null)
            {
                vm.BeginRecordShortcut(row);
                vm.RecordShortcutKey(Avalonia.Input.Key.P, cmdMod | Avalonia.Input.KeyModifiers.Shift);
            }
            var ran = vm.TryUserShortcut(Avalonia.Input.Key.P, cmdMod | Avalonia.Input.KeyModifiers.Shift, false);
            var toggled = vm.Perspective != persp;
            var shown = vm.PaletteRowsFor("perspective").FirstOrDefault(r => r.Id == "view.projection")?.Shortcut ?? "";
            if (row != null) { vm.BeginRecordShortcut(row); vm.RecordShortcutKey(Avalonia.Input.Key.I, cmdMod); }
            var conflict = vm.ShortcutNote;
            vm.ResetShortcuts();
            var cleared = !vm.TryUserShortcut(Avalonia.Input.Key.I, cmdMod, false);
            if (toggled) vm.Perspective = persp;
            Check(row != null && ran && toggled && shown.Contains('P') && conflict.Contains("Invert the selection") && cleared,
                  $"shortcuts: ran {ran}, toggled {toggled}, palette shows '{shown}', conflict: {conflict}, cleared {cleared}");
            vm.ShortcutFilter = "";
        }

        // The session: the project tree kept on quitting and restored from Start, with its job folders
        var keptItems = vm.ProjectItems.Select(p => p.Name.Replace(" (unsaved)", "")).ToList();
        var keptJobs = vm.ProjectItems.Sum(p => p.Jobs.Count);
        var keptActive = vm.ActiveItem?.Name.Replace(" (unsaved)", "");
        vm.SaveSession();

        // Close goes back to Start
        vm.SetModule(1);
        vm.CloseAllStructures();
        Check(vm.NoDocument && vm.IsStudio && vm.Title == "" && vm.ProjectItems.Count == 0, "close all structures: back to Start with no document");
        vm.LoadLastSession();
        var offered = vm.HasLastSession;
        var offerText = vm.LastSessionText;
        vm.RestoreSession();
        var backItems = vm.ProjectItems.Select(p => p.Name).ToList();
        var backJobs = vm.ProjectItems.Sum(p => p.Jobs.Count);
        Check(offered && backItems.SequenceEqual(keptItems) && backJobs == keptJobs && vm.ActiveItem?.Name == keptActive && !vm.HasLastSession,
              $"session: offered {offered} ({offerText}) · back {backItems.Count}/{keptItems.Count} structures, {backJobs}/{keptJobs} jobs, active {vm.ActiveItem?.Name} · {vm.Status}");
        // Clear on the pipeline strip: asked once more, then every structure closed and the strip empty
        {
            var had = vm.ProjectItems.Count;
            vm.ClearAll();
            var asked = vm.ClearArmed && vm.ProjectItems.Count == had;
            vm.ClearAll();
            Check(had > 0 && asked && vm.ProjectItems.Count == 0 && !vm.ShowPipelineStrip && vm.PipelineSteps.Count == 0,
                  $"clear: {had} structures · asked {asked} · left {vm.ProjectItems.Count} · strip {vm.ShowPipelineStrip}");
        }
        // Pack › Add molecule › Sulfur (S8), built from its SMILES, 20 packed around the polymer cell
        {
            vm.Open(Path.Combine(dir, "ps_melt.data"));
            vm.PackStart = 1;
            vm.NewPackInput();
            var s8 = vm.PackAdditives.First(f => f.Name.StartsWith("Sulfur", StringComparison.Ordinal));
            vm.AddPackMolecule(s8.Smiles, s8.Name).GetAwaiter().GetResult();
            var hasRow = vm.PackItems.Any(r => r.Name.Contains("Sulfur", StringComparison.Ordinal));
            vm.PackText = System.Text.RegularExpressions.Regex.Replace(vm.PackText, @"count\s+\d+", "count   20");
            vm.RunPack().GetAwaiter().GetResult();
            var sum = vm.Document?.Summary();
            Check(hasRow && vm.PackAdditives.Count >= 10 && sum?.Molecules == 30 && sum?.Atoms == 1300 + 160,
                  $"pack curatives: {vm.PackAdditives.Count} additives · S8 row {hasRow} · {sum?.Molecules} molecules, {sum?.Atoms} atoms · d_min {vm.PackDmin}");
            vm.PackStart = 0;
        }

        // Pack rows with their own force fields: toluene by GAFF, water as TIP3P; the cell assigned by groups
        {
            vm.PackStart = 0;
            vm.PackAssignField = true;   // an earlier check turned it off
            vm.NewPackInput();
            vm.AddPackMolecule("Cc1ccccc1", "toluene").GetAwaiter().GetResult();
            vm.AddPackMolecule("O", "water").GetAwaiter().GetResult();
            vm.SetPackRowCount(0, 6);
            vm.SetPackRowCount(1, 6);
            var rowsCounted = vm.PackItems.Count == 2 && vm.PackItems.All(r => r.CountValue == 6);
            var gaffRow = vm.Field.Library.ToList().FindIndex(e => e.Id.StartsWith("gaff", StringComparison.Ordinal)) + 1;
            vm.SetPackRowForceField(0, gaffRow);
            vm.SetPackRowForceField(1, 1 + vm.Field.Library.Count + FieldViewModel.Waters.FindIndex(w => w.Id == "tip3p"));
            var rowsText = vm.PackText;
            vm.RunPack().GetAwaiter().GetResult();
            Check(rowsCounted && rowsText.Contains("forcefield gaff") && rowsText.Contains("water      tip3p") && vm.Field.ForceFieldName.Contains("TIP3P", StringComparison.Ordinal) && vm.Field.ForceFieldName.Contains("GAFF", StringComparison.Ordinal),
                  $"pack rows' force fields: {vm.Field.ForceFieldName} · {vm.Field.Log}");
        }

        // Water model on the Field page: a water with TIP4P/2005 gets its M site (4 atoms) and the model's force field
        {
            vm.Open(Path.Combine(dir, "water.pdb"));
            var gaffIx = vm.Field.Library.ToList().FindIndex(e => e.Id.StartsWith("gaff", StringComparison.Ordinal));
            vm.Field.FfIndex = Math.Max(0, gaffIx);
            vm.Field.WaterModelIndex = FieldViewModel.Waters.FindIndex(w => w.Id == "tip4p2005") + 1;
            vm.Field.Assign().GetAwaiter().GetResult();
            var atoms = vm.Document?.Summary().Atoms ?? 0;
            Check(atoms == 4 && vm.Field.ForceFieldName.Contains("TIP4P/2005", StringComparison.Ordinal) && FieldViewModel.Waters.Count == 11,
                  $"water model: {atoms} atoms · {vm.Field.ForceFieldName} · {FieldViewModel.Waters.Count} models · {vm.Field.Log}");
            vm.Field.WaterModelIndex = 0;
        }

        // Remote copy-back: the out folder listed on the host (run here with sh), the result and small files now, a large
        // trajectory left with copy actions (whole, every 10th / 100th frame thinned on the host by caps frames)
        {
            var od = Path.Combine(Path.GetTempPath(), "caps-remote-out-test");
            Directory.CreateDirectory(Path.Combine(od, "out"));
            File.WriteAllText(Path.Combine(od, "out", "cell.data"), "x");
            File.WriteAllText(Path.Combine(od, "out", "run log.txt"), "hello");
            var psi = new System.Diagnostics.ProcessStartInfo("sh") { RedirectStandardOutput = true, UseShellExecute = false };
            psi.ArgumentList.Add("-c");
            psi.ArgumentList.Add(MainViewModel.ListOutCommand(od));
            using var lsProc = System.Diagnostics.Process.Start(psi)!;
            var listing = MainViewModel.ParseListing(lsProc.StandardOutput.ReadToEnd());
            lsProc.WaitForExit();
            Directory.Delete(od, true);
            var files = listing.Concat([new RemoteFile("traj.lammpstrj", 5L << 30), new RemoteFile("big.xtc", 300L << 20)]).ToList();
            var (now, later) = MainViewModel.CopyPlan(files, "cell", MainViewModel.LargeRemoteBytes);
            var thin = MainViewModel.ThinCommand("/scratch/u/dyn-1", "big.xtc", 10, "cell.data");
            Check(listing.Count == 2 && listing.Any(f => f.Name == "run log.txt" && f.Bytes == 5) && now.Count == 2 && later.Count == 2
                  && thin == "cd '/scratch/u/dyn-1/out' && caps frames 'big.xtc' 'big.every10.trr' --stride 10 --topology 'cell.data'"
                  && MainViewModel.CopyTimeoutMs(5L << 30) == 5120 * 1000 && new RemoteFile("t", 5L << 30).Size == "5.0 GB",
                  $"remote copy-back: listed {string.Join(", ", listing.Select(f => $"{f.Name} {f.Bytes} B"))} · now {now.Count}, left {string.Join(", ", later.Select(f => f.Name + " " + f.Size))} · {thin}");
        }

        // Open page › read from / to / every: frames 1 to the end of the sample dump (2 of 3), recorded with first=1
        {
            var dump = Path.Combine(dir, "ps_melt.lammpstrj");
            vm.PreviewOpenNow(dump, Path.Combine(dir, "ps_melt.data"));
            vm.SetOpenFrames(1, -1, 1);
            var button = vm.OpenButton;
            vm.ConfirmOpen();
            vm.Recording = true;
            vm.Open(dump, Path.Combine(dir, "ps_melt.data"));
            vm.Recording = false;
            var line = vm.RecordedCommands.LastOrDefault(c => c.Text.StartsWith("doc = caps.open("))?.Text ?? "";
            Check(button == "Open 2 frames" && vm.Frames == 2 && line.Contains("first=1") && vm.Notes.Any(n => n.Contains("2 kept of 3 read")),
                  $"open a frame selection: {button} · {vm.Frames} frames · {line}");
        }
        vm.CloseAllStructures();

        Console.WriteLine(fails == 0 ? "all checks passed" : $"{fails} check(s) failed");
        return fails == 0 ? 0 : 1;
    }
}
