using Avalonia;
using Avalonia.Headless;
using Avalonia.Threading;
using CapsStudio.Views;

namespace CapsStudio;

/// <summary>Renders the real main window off-screen (Skia, headless platform) to a PNG, for design review and CI.</summary>
internal static class Screenshot
{
    public static int Run(string[] args)
    {
        var output = args.Length > 0 ? args[0] : "studio.png";
        var file = args.Length > 1 ? args[1] : null;
        var topo = args.Length > 2 ? args[2] : null;
        AppBuilder.Configure<App>()
            .UseSkia()
            .UseHeadless(new AvaloniaHeadlessPlatformOptions { UseHeadlessDrawing = false })
            .WithCapsFonts()
            .SetupWithoutStarting();

        // never touch the user's recent list: $CAPS_RECENT_DIR, else a scratch folder
        ViewModels.AppSettings.Override = Environment.GetEnvironmentVariable("CAPS_SETTINGS") is { Length: > 0 } sf ? sf
            : Path.Combine(Path.GetTempPath(), "caps-screenshot-settings.json");
        ViewModels.RecentFiles.Override = Environment.GetEnvironmentVariable("CAPS_RECENT_DIR") is { Length: > 0 } rd ? rd
            : Path.Combine(Path.GetTempPath(), "caps-screenshot-recent");
        // CAPS_SHOT_SIZE=WxH (default 1440x900): the compact layout starts below 1440 px
        var size = (Environment.GetEnvironmentVariable("CAPS_SHOT_SIZE") ?? "1440x900").Split('x');
        var w = new MainWindow { Width = double.Parse(size[0]), Height = double.Parse(size[1]) };
        w.Show();
        if (!string.IsNullOrEmpty(file)) w.OpenOnStart(file, string.IsNullOrEmpty(topo) ? null : topo);
        // Let layout, the async render and the bitmap swap complete.
        for (var i = 0; i < 60; i++)
        {
            Dispatcher.UIThread.RunJobs();
            AvaloniaHeadlessPlatform.ForceRenderTimerTick();
            Thread.Sleep(30);
        }
        // Extra args (Field: ff=ID charges=N field=1 fieldpick=N fieldset=TYPE): pick=0,1,2  colour=3  style=4  module=0|1|2|3  grow=1  relax=1  mdsteps=N ensemble=0|1|2  md=1  eqprotocol=0|1|2 eqscale=X eq=1  tab=0..4  view=0|1  wrap=1
        foreach (var kv in args.Skip(3).Select(a => a.Split('=', 2)).Where(p => p.Length == 2))
        {
            if (kv[0] == "pick") w.PickForTest(kv[1].Split(',').Select(int.Parse).ToArray());
            if (kv[0] == "colour") w.ViewModel.ColourIndex = int.Parse(kv[1]);
            if (kv[0] == "style") w.ViewModel.StyleIndex = int.Parse(kv[1]);
            if (kv[0] == "module") w.ViewModel.SetModule(int.Parse(kv[1]));
            if (kv[0] == "figure")   // figure=BG: Export › Figure with that background selected
            {
                w.ViewModel.OpenFigure();
                w.ViewModel.FigBackground = int.Parse(kv[1]);
                for (int k = 0; k < 40; ++k) { Dispatcher.UIThread.RunJobs(); AvaloniaHeadlessPlatform.ForceRenderTimerTick(); Thread.Sleep(25); }
            }
            if (kv[0] == "figexport")   // figexport=PATH: export the figure (PNG or SVG by extension)
            {
                if (kv[1].EndsWith(".svg")) w.ViewModel.FigFormat = 1;
                var t = w.ViewModel.ExportFigure(kv[1], CapsStudio.Views.FigureDrawing.SavePng, CapsStudio.Views.FigureDrawing.AddToSvg);
                while (!t.IsCompleted) { Dispatcher.UIThread.RunJobs(); AvaloniaHeadlessPlatform.ForceRenderTimerTick(); Thread.Sleep(20); }
                Console.WriteLine(t.IsFaulted ? "figure export failed: " + t.Exception?.InnerException?.Message : "figure: " + t.Result);
            }
            if (kv[0] == "render") { w.ViewModel.OpenRender(); for (int k = 0; k < 60; ++k) { Dispatcher.UIThread.RunJobs(); AvaloniaHeadlessPlatform.ForceRenderTimerTick(); Thread.Sleep(25); } }
            if (kv[0] == "renderout")   // renderout=PATH: render the image with its overlays
            {
                var t = w.ViewModel.RenderOut(_ => kv[1], false, CapsStudio.Views.FigureDrawing.SaveRenderPng);
                while (!t.IsCompleted) { Dispatcher.UIThread.RunJobs(); AvaloniaHeadlessPlatform.ForceRenderTimerTick(); Thread.Sleep(20); }
                Console.WriteLine(t.IsFaulted ? "render failed: " + t.Exception?.InnerException?.Message : $"render: {t.Result} image(s)");
            }
            if (kv[0] == "visualize")   // visualize=PIPELINE.json|1: Analyze › Visualize with those steps (1: the default)
            {
                if (kv[1] != "1") w.ViewModel.LoadPipeline(kv[1]);
                w.ViewModel.OpenVisualize();
                for (int k = 0; k < 40; ++k) { Dispatcher.UIThread.RunJobs(); AvaloniaHeadlessPlatform.ForceRenderTimerTick(); Thread.Sleep(25); }
            }
            if (kv[0] == "inspector") { w.ViewModel.InspectorTab = int.Parse(kv[1]); for (int k = 0; k < 10; ++k) { Dispatcher.UIThread.RunJobs(); Thread.Sleep(20); } }
            if (kv[0] == "filter") w.ViewModel.InspectorFilter = kv[1];
            if (kv[0] == "steplib") w.ViewModel.StepLibraryOpen = true;
            if (kv[0] == "savepipeline") { w.ViewModel.PipelineName = kv[1]; w.ViewModel.OpenSavePipeline(); for (int k = 0; k < 20; ++k) { Dispatcher.UIThread.RunJobs(); Thread.Sleep(20); } }
            if (kv[0] == "openpreview") { w.ViewModel.PreviewOpen(kv[1]); for (int k = 0; k < 80 && !w.ViewModel.OpenIdle || k < 10; ++k) { Dispatcher.UIThread.RunJobs(); AvaloniaHeadlessPlatform.ForceRenderTimerTick(); Thread.Sleep(25); } for (int k = 0; k < 10; ++k) { Dispatcher.UIThread.RunJobs(); Thread.Sleep(20); } }
            if (kv[0] == "bundle") { w.ViewModel.OpenBundle(); w.ViewModel.BundleInput = kv[1] == "1"; for (int k = 0; k < 80; ++k) { Dispatcher.UIThread.RunJobs(); AvaloniaHeadlessPlatform.ForceRenderTimerTick(); Thread.Sleep(25); } }
            if (kv[0] == "viewports") { w.ViewModel.OpenViewports(); w.ViewModel.ViewportLayout = int.Parse(kv[1]); for (int k = 0; k < 80; ++k) { Dispatcher.UIThread.RunJobs(); AvaloniaHeadlessPlatform.ForceRenderTimerTick(); Thread.Sleep(25); } }
            if (kv[0] == "colourby") { w.ViewModel.OpenColourBy(); for (int k = 0; k < 80; ++k) { Dispatcher.UIThread.RunJobs(); AvaloniaHeadlessPlatform.ForceRenderTimerTick(); Thread.Sleep(25); } }
            if (kv[0] == "compare")   // compare=A|B: Analyze › Compare with those two files
            {
                var ab = kv[1].Split('|');
                w.ViewModel.OpenCompare();
                foreach (var (path, isA) in new[] { (ab[0], true), (ab[1], false) })
                {
                    var t = w.ViewModel.SetCompareInput(isA, path);
                    while (!t.IsCompleted) { Dispatcher.UIThread.RunJobs(); Thread.Sleep(20); }
                }
                for (int k = 0; k < 40; ++k) { Dispatcher.UIThread.RunJobs(); AvaloniaHeadlessPlatform.ForceRenderTimerTick(); Thread.Sleep(25); }
            }
            if (kv[0] == "batch")   // batch=PATTERN: run the current pipeline over the matching files
            {
                w.ViewModel.OpenBatch();
                w.ViewModel.BatchPattern = kv[1];
                var t = w.ViewModel.RunBatch();
                while (!t.IsCompleted) { Dispatcher.UIThread.RunJobs(); AvaloniaHeadlessPlatform.ForceRenderTimerTick(); Thread.Sleep(20); }
                for (int k = 0; k < 20; ++k) { Dispatcher.UIThread.RunJobs(); AvaloniaHeadlessPlatform.ForceRenderTimerTick(); Thread.Sleep(20); }
                Console.WriteLine("batch: " + w.ViewModel.BatchState + " · " + w.ViewModel.BatchOut);
            }
            if (kv[0] == "export")   // export=FORMATINDEX: Export › Data with that format, after the preview is written
            {
                w.ViewModel.OpenExport();
                w.ViewModel.ExportFormatIndex = int.Parse(kv[1]);
                for (int k = 0; k < 80 && !w.ViewModel.ExportIdle || k < 10; ++k) { Dispatcher.UIThread.RunJobs(); AvaloniaHeadlessPlatform.ForceRenderTimerTick(); Thread.Sleep(25); }
                for (int k = 0; k < 10; ++k) { Dispatcher.UIThread.RunJobs(); Thread.Sleep(20); }
            }
            if (kv[0] == "series")
            {
                var t = w.ViewModel.ComputeSeries();
                while (!t.IsCompleted) { Dispatcher.UIThread.RunJobs(); AvaloniaHeadlessPlatform.ForceRenderTimerTick(); Thread.Sleep(20); }
                for (int k = 0; k < 20; ++k) { Dispatcher.UIThread.RunJobs(); Thread.Sleep(20); }
                if (int.TryParse(kv[1], out var col)) w.ViewModel.PipeYColumn = col;
            }
            if (kv[0] == "pipestep") w.ViewModel.PipeSelected = w.ViewModel.PipelineRows[int.Parse(kv[1])];
            if (kv[0] == "focus") w.ViewModel.FocusOn(int.Parse(kv[1]));   // focus=N: keyboard-walk focus on atom N
            if (kv[0] == "walk")   // walk=keys: d(own) u(p) b(ond) ](next molecule) s(elect) m(easure)
                foreach (var c in kv[1])
                {
                    if (c == 'd') w.ViewModel.FocusStep(1); else if (c == 'u') w.ViewModel.FocusStep(-1); else if (c == 'b') w.ViewModel.FocusBond();
                    else if (c == ']') w.ViewModel.FocusMolecule(1); else if (c == 's') w.ViewModel.FocusSelect(); else if (c == 'm') w.ViewModel.FocusMeasure();
                }
            if (kv[0] == "waitload")   // waitload=F: until the progressive open has read F of the file (CAPS_PROGRESSIVE_BYTES, CAPS_LOAD_DELAY_MS)
            {
                var f = double.Parse(kv[1], System.Globalization.CultureInfo.InvariantCulture);
                for (int k = 0; k < 2000 && w.ViewModel.IsLoading && w.ViewModel.LoadFraction < f; ++k) { Dispatcher.UIThread.RunJobs(); AvaloniaHeadlessPlatform.ForceRenderTimerTick(); Thread.Sleep(10); }
            }
            if (kv[0] == "quick") w.ViewModel.QuickText = kv[1];
            if (kv[0] == "settab") { w.ViewModel.SetModule(10); w.ViewModel.SettingsTab = int.Parse(kv[1]); }
            if (kv[0] == "colours") w.ViewModel.SetPalette = int.Parse(kv[1]);
            if (kv[0] == "theme") w.ViewModel.SetTheme = kv[1];
            if (kv[0] == "jobs") w.ViewModel.SetModule(11);
            if (kv[0] == "polymer")
            {
                w.ViewModel.SetModule(13);
                w.ViewModel.LoadPolymerLibrary();
                if (w.ViewModel.PolymerLibrary.FirstOrDefault(e => e.Id == kv[1]) is { } entry) w.ViewModel.UseLibrary(entry, null);
                var t = w.ViewModel.BuildPolyPreview();
                while (!t.IsCompleted) { Dispatcher.UIThread.RunJobs(); AvaloniaHeadlessPlatform.ForceRenderTimerTick(); Thread.Sleep(20); }
            }
            if (kv[0] == "surface")   // surface=<crystal id>: the Surface builder with that crystal's slab previewed
            {
                w.ViewModel.OpenSurface();
                var ix = w.ViewModel.Crystals.ToList().FindIndex(c => c.Id == kv[1]);
                if (ix >= 0) w.ViewModel.SurfCrystal = ix;
                for (int k = 0; k < 40; ++k) { Dispatcher.UIThread.RunJobs(); AvaloniaHeadlessPlatform.ForceRenderTimerTick(); Thread.Sleep(25); }
                Console.WriteLine($"surface: termination {w.ViewModel.SurfTermination} of {w.ViewModel.SurfTerminations.Count}");
            }
            if (kv[0] == "nano")   // nano=0|1|2 (sheet, tube, particle): the Nanostructure builder with its preview
            {
                w.ViewModel.OpenNano();
                w.ViewModel.NanoKind = int.Parse(kv[1]);
                for (int k = 0; k < 40; ++k) { Dispatcher.UIThread.RunJobs(); AvaloniaHeadlessPlatform.ForceRenderTimerTick(); Thread.Sleep(25); }
            }
            if (kv[0] == "crystal")   // crystal=1: the Crystal builder (polyethylene); crystal=ID: a library crystal imported through Find symmetry
            {
                w.ViewModel.OpenCrystal();
                if (kv[1] != "1" && w.ViewModel.Crystals.FirstOrDefault(c => c.Id == kv[1]) is { } cx)
                {
                    var t = w.ViewModel.ImportCrystalCif(cx.File);
                    while (!t.IsCompleted) { Dispatcher.UIThread.RunJobs(); Thread.Sleep(20); }
                }
                for (int k = 0; k < 40; ++k) { Dispatcher.UIThread.RunJobs(); AvaloniaHeadlessPlatform.ForceRenderTimerTick(); Thread.Sleep(25); }
            }
            if (kv[0] == "crystalsuper") { w.ViewModel.CrystalSupercell = kv[1].Replace('x', '×'); for (int k = 0; k < 40; ++k) { Dispatcher.UIThread.RunJobs(); AvaloniaHeadlessPlatform.ForceRenderTimerTick(); Thread.Sleep(25); } }
            if (kv[0] == "crystalquery") { w.ViewModel.CrystalQuery = kv[1]; Dispatcher.UIThread.RunJobs(); }
            if (kv[0] == "bio")   // bio=1: the Biomolecule builder (the board's peptide); bio=SEQUENCE sets the sequence
            {
                w.ViewModel.OpenBio();
                if (kv[1] != "1") w.ViewModel.BioSequence = kv[1];
                for (int k = 0; k < 40; ++k) { Dispatcher.UIThread.RunJobs(); AvaloniaHeadlessPlatform.ForceRenderTimerTick(); Thread.Sleep(25); }
            }
            if (kv[0] == "bioss")   // bioss=FROM-TO:K selects residues FROM..TO (1-based) and applies structure K (0 helix 1 strand 2 PPII 3 coil)
            {
                var parts = kv[1].Split(':');
                var range = parts[0].Split('-').Select(int.Parse).ToArray();
                w.ViewModel.SelectResidue(w.ViewModel.BioCells[range[0] - 1], false);
                w.ViewModel.SelectResidue(w.ViewModel.BioCells[range[^1] - 1], true);
                w.ViewModel.BioType = int.Parse(parts[1]);
                for (int k = 0; k < 40; ++k) { Dispatcher.UIThread.RunJobs(); AvaloniaHeadlessPlatform.ForceRenderTimerTick(); Thread.Sleep(25); }
            }
            if (kv[0] == "solvation")   // solvation=1: the Solvation builder around the open structure (or pure solvent); solvation=ID picks a solvent
            {
                w.ViewModel.OpenSolvation();
                if (kv[1] != "1" && w.ViewModel.Solvents.FirstOrDefault(sv => sv.Id == kv[1]) is { } svi) w.ViewModel.SolvSolvent = svi;
                for (int k = 0; k < 40; ++k) { Dispatcher.UIThread.RunJobs(); AvaloniaHeadlessPlatform.ForceRenderTimerTick(); Thread.Sleep(25); }
            }
            if (kv[0] == "solvpad") { w.ViewModel.SolvPadding = decimal.Parse(kv[1], System.Globalization.CultureInfo.InvariantCulture); for (int k = 0; k < 40; ++k) { Dispatcher.UIThread.RunJobs(); AvaloniaHeadlessPlatform.ForceRenderTimerTick(); Thread.Sleep(25); } }
            if (kv[0] == "appearance") { w.ViewModel.AppearanceOpen = kv[1] == "1"; Dispatcher.UIThread.RunJobs(); }
            if (kv[0] == "appexpr") { w.ViewModel.AppTarget = 1; w.ViewModel.AppExpression = kv[1]; }   // the target of the next appstyle
            if (kv[0] == "appstyle") w.ViewModel.AppStyle = int.Parse(kv[1]);
            if (kv[0] == "appcolour") w.ViewModel.AppColour = int.Parse(kv[1]);
            if (kv[0] == "appsurface") w.ViewModel.AppSurface = int.Parse(kv[1]);
            if (kv[0] == "appmap") w.ViewModel.AppSurfaceMap = int.Parse(kv[1]);
            if (kv[0] == "appsurfatoms") w.ViewModel.AppSurfaceAtoms = kv[1];
            if (kv[0] == "labels")   // labels=element,rs,type,charge
            {
                var ks = kv[1].Split(',');
                w.ViewModel.LabelElement = ks.Contains("element"); w.ViewModel.LabelRs = ks.Contains("rs");
                w.ViewModel.LabelType = ks.Contains("type"); w.ViewModel.LabelCharge = ks.Contains("charge");
            }
            if (kv[0] == "appwait") for (int k = 0; k < 60; ++k) { Dispatcher.UIThread.RunJobs(); AvaloniaHeadlessPlatform.ForceRenderTimerTick(); Thread.Sleep(30); }
            if (kv[0] == "player")   // player=1: the trajectory player (a log beside the file is picked up)
            {
                w.ViewModel.OpenTrajectory();
                for (int k = 0; k < 60; ++k) { Dispatcher.UIThread.RunJobs(); AvaloniaHeadlessPlatform.ForceRenderTimerTick(); Thread.Sleep(30); }
            }
            if (kv[0] == "frame") { w.ViewModel.Frame = int.Parse(kv[1]); for (int k = 0; k < 30; ++k) { Dispatcher.UIThread.RunJobs(); AvaloniaHeadlessPlatform.ForceRenderTimerTick(); Thread.Sleep(25); } }
            if (kv[0] == "torsion")   // torsion=1: the scan page with the default torsion; torsion=run also runs it
            {
                w.ViewModel.OpenTorsion();
                if (kv[1] == "run") { var t = w.ViewModel.RunTorsionScan(); while (!t.IsCompleted) { Dispatcher.UIThread.RunJobs(); Thread.Sleep(20); } }
                for (int k = 0; k < 40; ++k) { Dispatcher.UIThread.RunJobs(); AvaloniaHeadlessPlatform.ForceRenderTimerTick(); Thread.Sleep(25); }
            }
            if (kv[0] == "element") w.ViewModel.BuildElement = kv[1];
            if (kv[0] == "picker") { w.ViewModel.ElementPickerOpen = kv[1] == "1"; Dispatcher.UIThread.RunJobs(); w.MarkChosenForTest(); }
            if (kv[0] == "edit")   // edit=add_h | clean | tool=N
            {
                if (kv[1] == "add_h") w.ViewModel.AddHydrogensAll();
                else if (kv[1] == "clean") { var t = w.ViewModel.AutoClean(); while (!t.IsCompleted) { Dispatcher.UIThread.RunJobs(); Thread.Sleep(20); } }
                else if (kv[1].StartsWith("tool")) w.ViewModel.EditTool = int.Parse(kv[1][4..]);
                for (int k = 0; k < 20; ++k) { Dispatcher.UIThread.RunJobs(); AvaloniaHeadlessPlatform.ForceRenderTimerTick(); Thread.Sleep(25); }
            }
            if (kv[0] == "selection") { w.ViewModel.SelectionOpen = kv[1] == "1"; Dispatcher.UIThread.RunJobs(); }
            if (kv[0] == "selsmarts") { w.ViewModel.SelectMode = 0; w.ViewModel.SelectPattern = kv[1]; w.ViewModel.RunSelect("replace"); w.ViewModel.SaveSelectionAsSet(); }
            if (kv[0] == "seliso") { var t = w.ViewModel.MakeTactic(kv[1] == "1"); while (!t.IsCompleted) { Dispatcher.UIThread.RunJobs(); Thread.Sleep(20); } }
            if (kv[0] == "checks") { w.ViewModel.InteractionsOpen = kv[1] == "1"; for (int k = 0; k < 20; ++k) { Dispatcher.UIThread.RunJobs(); AvaloniaHeadlessPlatform.ForceRenderTimerTick(); Thread.Sleep(25); } }
            if (kv[0] == "split")   // split=PATH: the split view with PATH on the right
            {
                w.ViewModel.OpenSplit();
                var t = w.ViewModel.SetSplitB(kv[1]);
                while (!t.IsCompleted) { Dispatcher.UIThread.RunJobs(); Thread.Sleep(20); }
                for (int k = 0; k < 40; ++k) { Dispatcher.UIThread.RunJobs(); AvaloniaHeadlessPlatform.ForceRenderTimerTick(); Thread.Sleep(25); }
            }
            if (kv[0] == "theme") { w.ViewModel.SetTheme = kv[1]; for (int k = 0; k < 10; ++k) { Dispatcher.UIThread.RunJobs(); Thread.Sleep(20); } }
            if (kv[0] == "settingstab") w.ViewModel.SettingsTab = int.Parse(kv[1]);
            if (kv[0] == "a11y")   // a11y=1: list the visible buttons a screen reader would reach without a name
            {
                for (int k = 0; k < 10; ++k) { Dispatcher.UIThread.RunJobs(); AvaloniaHeadlessPlatform.ForceRenderTimerTick(); Thread.Sleep(20); }
                var un = Views.AccessibleNames.Unnamed(w);
                var all = Avalonia.VisualTree.VisualExtensions.GetVisualDescendants(w).OfType<Avalonia.Controls.Button>().Where(b => b.IsEffectivelyVisible).ToList();
                var named = all.Where(b => !string.IsNullOrEmpty(Avalonia.Automation.AutomationProperties.GetName(b))).ToList();
                Console.WriteLine($"visible buttons: {all.Count} · named from tooltips: {named.Count} · e.g. {string.Join(" | ", named.Take(6).Select(b => Avalonia.Automation.AutomationProperties.GetName(b)))}");
                Console.WriteLine($"unnamed buttons: {un.Count}");
                foreach (var b in un.Take(40)) Console.WriteLine($"  {b.GetType().Name} {b.Name} {(b.Content as Avalonia.Controls.Control)?.GetType().Name} classes={string.Join(".", b.Classes)}");
            }
            if (kv[0] == "update")   // update=1: check for updates (set CAPS_UPDATE_FEED to a local release JSON)
            {
                var t = w.ViewModel.CheckForUpdates();
                while (!t.IsCompleted) { Dispatcher.UIThread.RunJobs(); Thread.Sleep(20); }
                for (int k = 0; k < 20; ++k) { Dispatcher.UIThread.RunJobs(); AvaloniaHeadlessPlatform.ForceRenderTimerTick(); Thread.Sleep(25); }
            }
            if (kv[0] == "history")   // history=1: the History panel after a snapshot, three edits and an undo
            {
                var vm = w.ViewModel;
                vm.HistoryOpen = true;
                vm.TakeSnapshot("as opened");
                vm.BuildElement = "O"; vm.EditTool = 1;
                vm.ToolClick(-1); vm.ToolClick(-1);
                vm.TakeSnapshot("two waters' oxygens");
                vm.ToolClick(-1);
                vm.EditTool = 0;
                vm.UndoEdit(false);
                for (int k = 0; k < 20; ++k) { Dispatcher.UIThread.RunJobs(); AvaloniaHeadlessPlatform.ForceRenderTimerTick(); Thread.Sleep(25); }
            }
            if (kv[0] == "analysis")   // analysis=glass|interface|diffusion[:run] — the focus page, optionally run (glass: a short scan)
            {
                var vm = w.ViewModel;
                var parts = kv[1].Split(':');
                System.Threading.Tasks.Task? t = null;
                switch (parts[0])
                {
                    case "glass":
                        vm.Analyze.TgFromD = 450; vm.Analyze.TgToD = 250; vm.Analyze.TgStepD = 25; vm.Analyze.TgPsD = 5; vm.Analyze.EqPsD = 2; vm.GtReplicas = 3;
                        vm.OpenGlass();
                        if (parts.Length > 1) t = vm.RunGlass();
                        break;
                    case "interface": vm.OpenInterface(); if (parts.Length > 1) t = vm.RunInterface(); break;
                    case "diffusion": vm.OpenDiffusion(); if (parts.Length > 1) t = vm.RunDiffusion(); break;
                }
                while (t is { IsCompleted: false }) { Dispatcher.UIThread.RunJobs(); Thread.Sleep(20); }
                for (int k = 0; k < 30; ++k) { Dispatcher.UIThread.RunJobs(); AvaloniaHeadlessPlatform.ForceRenderTimerTick(); Thread.Sleep(25); }
            }
            if (kv[0] == "query")   // query=TEXT: the ⌘F query bar with that query, applied
            {
                w.ViewModel.QueryOpen = true;
                w.ViewModel.QueryText = kv[1];
                w.ViewModel.ApplyQuery();
                for (int k = 0; k < 20; ++k) { Dispatcher.UIThread.RunJobs(); AvaloniaHeadlessPlatform.ForceRenderTimerTick(); Thread.Sleep(25); }
            }
            if (kv[0] == "vision") w.ViewModel.OpenColourVision();                    // vision=1: Settings › Colour vision
            if (kv[0] == "cvd") w.ViewModel.VisionPreview = int.Parse(kv[1]);          // cvd=0..3: the view as seen with a deficiency
            if (kv[0] == "host")   // host=NAME,HOSTNAME,USER: a remote host added and tested
            {
                var f = kv[1].Split(',');
                w.ViewModel.AddHost();
                w.ViewModel.HostName = f[0];
                if (f.Length > 1) w.ViewModel.HostHostname = f[1];
                if (f.Length > 2) w.ViewModel.HostUser = f[2];
                var t = w.ViewModel.TestHost();
                while (!t.IsCompleted) { Dispatcher.UIThread.RunJobs(); Thread.Sleep(50); }
            }
            if (kv[0] == "macro")   // macro=1: a recorded session (open, place, field, relax, save), a parameter, the run
            {
                var vm = w.ViewModel;
                CapsStudio.ViewModels.MainViewModel.MacroFolderOverride = Path.Combine(Path.GetTempPath(), "caps-shot-macros");
                if (Directory.Exists(CapsStudio.ViewModels.MainViewModel.MacroFolderOverride)) Directory.Delete(CapsStudio.ViewModels.MainViewModel.MacroFolderOverride, true);
                vm.OpenMacro();
                vm.Recording = true;
                if (vm.Document?.Path is { } path) vm.RecordOpen(path, null);
                vm.RecordEdit("{\"op\":\"place\",\"smiles\":\"CO\",\"name\":\"Methanol\",\"resname\":\"MOH\"}");
                vm.Record("doc.field.assign(\"uff\", charges=\"qeq\")");
                vm.Record("doc.relax(ftol=0.5, method=\"lbfgs\", max_iterations=5000)");
                vm.Record("print(\"atoms\", doc.atoms, \"·\", doc.summary()[\"density\"], \"g/cm3\")");
                vm.PromoteToParameter("0.5", "ftol");
                var t = vm.RunMacro();
                while (!t.IsCompleted) { Dispatcher.UIThread.RunJobs(); Thread.Sleep(50); }
                for (int k = 0; k < 60; ++k) { Dispatcher.UIThread.RunJobs(); AvaloniaHeadlessPlatform.ForceRenderTimerTick(); Thread.Sleep(25); }
            }
            if (kv[0] == "import")   // import=PATH: the import dialog on that file
            {
                w.ViewModel.ShowImport(kv[1]);
                w.ViewModel.WaitImport();
                for (int k = 0; k < 60; ++k) { Dispatcher.UIThread.RunJobs(); AvaloniaHeadlessPlatform.ForceRenderTimerTick(); Thread.Sleep(25); }
            }
            if (kv[0] == "importff") { w.ViewModel.ImportForceField = int.Parse(kv[1]); Dispatcher.UIThread.RunJobs(); }
            if (kv[0] == "importok")   // importok=1: confirm the import dialog opened by import=PATH
            {
                w.ViewModel.ImportChecks = false;
                var t = w.ViewModel.ConfirmImport();
                while (!t.IsCompleted) { Dispatcher.UIThread.RunJobs(); Thread.Sleep(20); }
                for (int k = 0; k < 40; ++k) { Dispatcher.UIThread.RunJobs(); AvaloniaHeadlessPlatform.ForceRenderTimerTick(); Thread.Sleep(25); }
            }
            if (kv[0] == "exportdlg")   // exportdlg=image|movie: the export dialog
            {
                w.ViewModel.OpenExportDialog(kv[1] == "movie" ? 1 : 0);
                for (int k = 0; k < 20; ++k) { Dispatcher.UIThread.RunJobs(); AvaloniaHeadlessPlatform.ForceRenderTimerTick(); Thread.Sleep(25); }
                Avalonia.Controls.ControlExtensions.FindControl<CapsStudio.Views.Pages.ExportDialog>(w, "ExportDlg")?.RefreshPreview();
                for (int k = 0; k < 40; ++k) { Dispatcher.UIThread.RunJobs(); AvaloniaHeadlessPlatform.ForceRenderTimerTick(); Thread.Sleep(25); }
            }
            if (kv[0] == "exportbg") w.ViewModel.ExportDlgBackground = int.Parse(kv[1]);
            if (kv[0] == "provenance")   // provenance=OTHER|-: the provenance page, compared with OTHER's
            {
                w.ViewModel.OpenProvenance();
                if (kv[1] != "-") w.ViewModel.CompareProvenanceWith(kv[1]);
                for (int k = 0; k < 40; ++k) { Dispatcher.UIThread.RunJobs(); AvaloniaHeadlessPlatform.ForceRenderTimerTick(); Thread.Sleep(25); }
            }
            if (kv[0] == "mechanics")   // mechanics=run|show: the Mechanics page (run: a short strain + tensile run)
            {
                var vm = w.ViewModel;
                vm.OpenMechanics();
                if (kv[1] == "run")
                {
                    vm.MechRate = 2;
                    vm.Analyze.EqPsD = 2;
                    vm.Analyze.TensMaxD = 0.1m;
                    var t = vm.RunMechanics();
                    while (!t.IsCompleted) { Dispatcher.UIThread.RunJobs(); Thread.Sleep(50); }
                }
                for (int k = 0; k < 40; ++k) { Dispatcher.UIThread.RunJobs(); AvaloniaHeadlessPlatform.ForceRenderTimerTick(); Thread.Sleep(25); }
            }
            if (kv[0] == "scattering")   // scattering=run|show[,EXPERIMENT]: the Scattering page
            {
                var vm = w.ViewModel;
                var f = kv[1].Split(',');
                vm.OpenScattering();
                if (f.Length > 1) vm.LoadExperiment(f[1]);
                if (f[0] == "run")
                {
                    var t = vm.RunScattering();
                    while (!t.IsCompleted) { Dispatcher.UIThread.RunJobs(); Thread.Sleep(50); }
                }
                for (int k = 0; k < 40; ++k) { Dispatcher.UIThread.RunJobs(); AvaloniaHeadlessPlatform.ForceRenderTimerTick(); Thread.Sleep(25); }
            }
            if (kv[0] == "isotope") w.ViewModel.IsotopePattern = int.Parse(kv[1]);
            if (kv[0] == "freevolume")   // freevolume=run|show: the Free volume page
            {
                var vm = w.ViewModel;
                vm.OpenFreeVolume();
                if (kv[1] == "run")
                {
                    var t = vm.RunFreeVolume();
                    while (!t.IsCompleted) { Dispatcher.UIThread.RunJobs(); Thread.Sleep(50); }
                }
                for (int k = 0; k < 60; ++k) { Dispatcher.UIThread.RunJobs(); AvaloniaHeadlessPlatform.ForceRenderTimerTick(); Thread.Sleep(25); }
            }
            if (kv[0] == "pore")   // pore=slit|cylinder|framework: the nanostructure builder on a pore
            {
                var vm = w.ViewModel;
                vm.OpenNano();
                vm.NanoKind = 3;
                vm.PoreType = kv[1] switch { "cylinder" => 1, "framework" => 2, _ => 0 };
                if (kv[1] == "cylinder") { vm.PoreCrystal = Math.Max(0, vm.Crystals.ToList().FindIndex(c => c.Name.Contains("quartz", StringComparison.OrdinalIgnoreCase))); vm.PoreWidth = 14; vm.PoreFluidIndex = 2; vm.PoreCount = 20; }
                for (int k = 0; k < 80; ++k) { Dispatcher.UIThread.RunJobs(); AvaloniaHeadlessPlatform.ForceRenderTimerTick(); Thread.Sleep(25); }
            }
            if (kv[0] == "tour")   // tour=STEP (1-based): the first-run tour on that step
            {
                w.ViewModel.StartTour();
                for (int k = 1; k < int.Parse(kv[1]); ++k) w.ViewModel.TourNext();
                for (int k = 0; k < 40; ++k) { Dispatcher.UIThread.RunJobs(); AvaloniaHeadlessPlatform.ForceRenderTimerTick(); Thread.Sleep(25); }
            }
            if (kv[0] == "manual")   // manual=PAGE-ID: the theory manual on that page
            {
                w.ViewModel.OpenManual(kv[1]);
                for (int k = 0; k < 40; ++k) { Dispatcher.UIThread.RunJobs(); AvaloniaHeadlessPlatform.ForceRenderTimerTick(); Thread.Sleep(25); }
            }
            if (kv[0] == "project")   // project=FOLDER|-: the project home (- = the open structure's folder)
            {
                w.ViewModel.OpenProject(kv[1] == "-" ? null : kv[1]);
                for (int k = 0; k < 120; ++k) { Dispatcher.UIThread.RunJobs(); AvaloniaHeadlessPlatform.ForceRenderTimerTick(); Thread.Sleep(25); }
            }
            if (kv[0] == "sweep")   // sweep=plan|run: the parameter sweep (run: DP 3,5 · seeds 1,2 · 3 chains, into a temp folder)
            {
                var vm = w.ViewModel;
                vm.SweepFolder = Path.Combine(Path.GetTempPath(), "caps-shot-sweep");
                vm.OpenSweep();
                if (kv[1] == "run")
                {
                    vm.SweepDps = "3, 5"; vm.SweepSeeds = "1, 2"; vm.SweepChains = 3; vm.SweepSyn = false;
                    var t = vm.RunSweep();
                    while (!t.IsCompleted) { Dispatcher.UIThread.RunJobs(); Thread.Sleep(50); }
                }
                for (int k = 0; k < 40; ++k) { Dispatcher.UIThread.RunJobs(); AvaloniaHeadlessPlatform.ForceRenderTimerTick(); Thread.Sleep(25); }
            }
            if (kv[0] == "cg")   // cg=1: the coarse-grained builder
            {
                w.ViewModel.OpenCg();
                for (int k = 0; k < 60; ++k) { Dispatcher.UIThread.RunJobs(); AvaloniaHeadlessPlatform.ForceRenderTimerTick(); Thread.Sleep(25); }
            }
            if (kv[0] == "template")   // template=NAME: the reaction template editor
            {
                w.ViewModel.OpenTemplateEditor(kv[1]);
                w.ViewModel.TestTemplate();
                for (int k = 0; k < 40; ++k) { Dispatcher.UIThread.RunJobs(); AvaloniaHeadlessPlatform.ForceRenderTimerTick(); Thread.Sleep(25); }
            }
            if (kv[0] == "lod")   // lod=1: level of detail on, the Performance panel and HUD
            {
                w.ViewModel.LodOn = true; w.ViewModel.PerfHud = true; w.ViewModel.LodOpen = true;
                for (int k = 0; k < 80; ++k) { Dispatcher.UIThread.RunJobs(); AvaloniaHeadlessPlatform.ForceRenderTimerTick(); Thread.Sleep(25); }
            }
            if (kv[0] == "fragments")   // fragments=CATEGORY: the library on that category
            {
                w.ViewModel.OpenFragments();
                w.ViewModel.FragmentCategory = kv[1];
                for (int k = 0; k < 120; ++k) { Dispatcher.UIThread.RunJobs(); AvaloniaHeadlessPlatform.ForceRenderTimerTick(); Thread.Sleep(25); }
            }
            if (kv[0] == "fragment") { w.ViewModel.SelectedFragment = w.ViewModel.FragmentTiles.FirstOrDefault(f => f.Name == kv[1]); for (int k = 0; k < 40; ++k) { Dispatcher.UIThread.RunJobs(); AvaloniaHeadlessPlatform.ForceRenderTimerTick(); Thread.Sleep(25); } }
            if (kv[0] == "nanomatrix") w.ViewModel.NanoMatrix = kv[1] == "1";
            if (kv[0] == "blend")
            {
                w.ViewModel.OpenBlend();
                for (int k = 0; k < 20; ++k) { Dispatcher.UIThread.RunJobs(); AvaloniaHeadlessPlatform.ForceRenderTimerTick(); Thread.Sleep(25); }
            }
            if (kv[0] == "nanobuild")
            {
                var t = w.ViewModel.BuildNano();
                while (!t.IsCompleted) { Dispatcher.UIThread.RunJobs(); AvaloniaHeadlessPlatform.ForceRenderTimerTick(); Thread.Sleep(20); }
            }
            if (kv[0] == "surfbuild")   // grows the film (or builds the slab) and opens it in the Studio
            {
                var t = w.ViewModel.BuildSurface();
                while (!t.IsCompleted) { Dispatcher.UIThread.RunJobs(); AvaloniaHeadlessPlatform.ForceRenderTimerTick(); Thread.Sleep(20); }
            }
            if (kv[0] == "polygrow")
            {
                w.ViewModel.SendPolymerToGrow();
                var t = w.ViewModel.Grow();
                while (!t.IsCompleted) { Dispatcher.UIThread.RunJobs(); AvaloniaHeadlessPlatform.ForceRenderTimerTick(); Thread.Sleep(20); }
            }
            if (kv[0] == "bench")
            {
                w.ViewModel.SetModule(12);
                w.ViewModel.LoadBench();
                w.ViewModel.BenchQuick = true;
                w.ViewModel.BenchRepeats = 1;
                if (kv[1] != "0")
                {
                    var t = w.ViewModel.RunBench(true);
                    while (!t.IsCompleted) { Dispatcher.UIThread.RunJobs(); AvaloniaHeadlessPlatform.ForceRenderTimerTick(); Thread.Sleep(20); }
                }
                w.ViewModel.SelectedBench = w.ViewModel.BenchItems.FirstOrDefault(b => b.Id == "T5");
            }
            if (kv[0] == "palette") { w.ViewModel.PaletteOpen = true; w.ViewModel.PaletteQuery = kv[1]; }
            if (kv[0] == "molecule")
            {
                w.ViewModel.OpenBuilder(kv[1]);
                var t = w.ViewModel.BuildMolecule();
                while (!t.IsCompleted) { Dispatcher.UIThread.RunJobs(); AvaloniaHeadlessPlatform.ForceRenderTimerTick(); Thread.Sleep(20); }
            }
            // pointer input in window coordinates: click=X,Y  drag=X1,Y1,X2,Y2 (after any building so the layout is settled)
            if (kv[0] is "click" or "drag")
            {
                var v = kv[1].Split(',').Select(x => double.Parse(x, System.Globalization.CultureInfo.InvariantCulture)).ToArray();
                var a = new Point(v[0], v[1]);
                w.MouseMove(a);
                w.MouseDown(a, Avalonia.Input.MouseButton.Left);
                if (kv[0] == "drag")
                {
                    var b = new Point(v[2], v[3]);
                    for (var k = 1; k <= 8; k++) w.MouseMove(new Point(a.X + (b.X - a.X) * k / 8, a.Y + (b.Y - a.Y) * k / 8), Avalonia.Input.RawInputModifiers.LeftMouseButton);
                    w.MouseUp(b, Avalonia.Input.MouseButton.Left);
                }
                else w.MouseUp(a, Avalonia.Input.MouseButton.Left);
                for (var i = 0; i < 40; i++) { Dispatcher.UIThread.RunJobs(); AvaloniaHeadlessPlatform.ForceRenderTimerTick(); Thread.Sleep(20); }
            }
            if (kv[0] == "wait")
                for (var i = 0; i < int.Parse(kv[1]) / 20; i++) { Dispatcher.UIThread.RunJobs(); AvaloniaHeadlessPlatform.ForceRenderTimerTick(); Thread.Sleep(20); }
            if (kv[0] == "conformers") w.ViewModel.MolConfCount = decimal.Parse(kv[1], System.Globalization.CultureInfo.InvariantCulture);
            if (kv[0] == "grow")
            {
                w.ViewModel.SetModule(0);
                var t = w.ViewModel.Grow();
                while (!t.IsCompleted) { Dispatcher.UIThread.RunJobs(); Thread.Sleep(20); }
            }
            if (kv[0] == "relax")
            {
                w.ViewModel.SetModule(2);
                var t = w.ViewModel.Relax();
                while (!t.IsCompleted) { Dispatcher.UIThread.RunJobs(); AvaloniaHeadlessPlatform.ForceRenderTimerTick(); Thread.Sleep(20); }
            }
            if (kv[0] == "mdsteps") w.ViewModel.MdStepsD = decimal.Parse(kv[1], System.Globalization.CultureInfo.InvariantCulture);
            if (kv[0] == "ensemble") w.ViewModel.MdEnsemble = int.Parse(kv[1]);
            if (kv[0] == "md")
            {
                w.ViewModel.SetModule(3);
                var t = w.ViewModel.RunMd();
                while (!t.IsCompleted) { Dispatcher.UIThread.RunJobs(); AvaloniaHeadlessPlatform.ForceRenderTimerTick(); Thread.Sleep(20); }
            }
            if (kv[0] == "eqprotocol") w.ViewModel.EqProtocol = int.Parse(kv[1]);
            if (kv[0] == "equntil") w.ViewModel.EqUntilConverged = kv[1] == "1";
            if (kv[0] == "eqblock") w.ViewModel.EqBlockD = decimal.Parse(kv[1], System.Globalization.CultureInfo.InvariantCulture);
            if (kv[0] == "eqmax") w.ViewModel.EqMaxBlocksD = decimal.Parse(kv[1], System.Globalization.CultureInfo.InvariantCulture);
            if (kv[0] == "eqscale") w.ViewModel.EqScaleD = decimal.Parse(kv[1], System.Globalization.CultureInfo.InvariantCulture);
            if (kv[0] == "eq")
            {
                w.ViewModel.SetModule(4);
                var t = w.ViewModel.RunEquilibrate();
                while (!t.IsCompleted) { Dispatcher.UIThread.RunJobs(); AvaloniaHeadlessPlatform.ForceRenderTimerTick(); Thread.Sleep(20); }
            }
            if (kv[0] == "packexample") { w.ViewModel.SetModule(5); w.ViewModel.PackCountD = decimal.Parse(kv[1], System.Globalization.CultureInfo.InvariantCulture); w.RunPackExampleForTest(); }
            if (kv[0] == "pack")
            {
                w.ViewModel.SetModule(5);
                var t = w.ViewModel.RunPack();
                while (!t.IsCompleted) { Dispatcher.UIThread.RunJobs(); AvaloniaHeadlessPlatform.ForceRenderTimerTick(); Thread.Sleep(20); }
            }
            if (kv[0] == "rxcycles") w.ViewModel.RxCyclesD = decimal.Parse(kv[1], System.Globalization.CultureInfo.InvariantCulture);
            if (kv[0] == "rxmd") w.ViewModel.RxMdPsD = decimal.Parse(kv[1], System.Globalization.CultureInfo.InvariantCulture);
            if (kv[0] == "react")
            {
                w.ViewModel.SetModule(6);
                var t = w.ViewModel.RunReact();
                while (!t.IsCompleted) { Dispatcher.UIThread.RunJobs(); AvaloniaHeadlessPlatform.ForceRenderTimerTick(); Thread.Sleep(20); }
            }
            // Field: ff=<catalogue id>  charges=0|1|2  field=1 (assign)  fieldpick=<atom index>  fieldset=<type> (override the picked atom)
            if (kv[0] == "ff")
            {
                var k = w.ViewModel.Field.Library.ToList().FindIndex(x => x.Id == kv[1]);
                if (k >= 0) w.ViewModel.Field.FfIndex = k;
            }
            if (kv[0] == "charges") w.ViewModel.Field.ChargeMode = int.Parse(kv[1]);
            if (kv[0] == "compress") w.ViewModel.RelaxCompress = kv[1] == "1";
            if (kv[0] == "field")
            {
                w.ViewModel.SetModule(7);
                var t = w.ViewModel.Field.Assign();
                while (!t.IsCompleted) { Dispatcher.UIThread.RunJobs(); AvaloniaHeadlessPlatform.ForceRenderTimerTick(); Thread.Sleep(20); }
            }
            if (kv[0] == "fieldpick") w.ViewModel.Field.SelectAtom(int.Parse(kv[1]));
            if (kv[0] == "fieldmissing") w.ViewModel.Field.SelectedMissing = w.ViewModel.Field.Missing.ElementAtOrDefault(int.Parse(kv[1]));
            if (kv[0] == "fieldset")
            {
                w.ViewModel.Field.OverrideType = kv[1];
                var t = w.ViewModel.Field.ApplyOverride();
                while (!t.IsCompleted) { Dispatcher.UIThread.RunJobs(); AvaloniaHeadlessPlatform.ForceRenderTimerTick(); Thread.Sleep(20); }
            }
            // Analyze › Properties: props=density,rdf,...  frameps=X  ref=<material id>  analyze=1 (run)  curve=N
            if (kv[0] == "props")
            {
                w.ViewModel.SetModule(1);
                w.ViewModel.AnalyzeProperties = true;
                var ids = kv[1].Split(',');
                foreach (var c in w.ViewModel.Analyze.Groups.SelectMany(g => g.Chips)) c.IsOn = ids.Contains(c.Id);
            }
            if (kv[0] == "frameps") w.ViewModel.Analyze.FramePsD = decimal.Parse(kv[1], System.Globalization.CultureInfo.InvariantCulture);
            if (kv[0] == "ref") w.ViewModel.Analyze.RefIndex = w.ViewModel.Analyze.References.ToList().FindIndex(r => r.Id == kv[1]);
            if (kv[0] == "analyze")
            {
                var t = w.ViewModel.Analyze.Run();
                while (!t.IsCompleted) { Dispatcher.UIThread.RunJobs(); AvaloniaHeadlessPlatform.ForceRenderTimerTick(); Thread.Sleep(20); }
            }
            if (kv[0] == "curve") w.ViewModel.Analyze.CurveIndex = int.Parse(kv[1]);
            // protocol settings: tgset=from:to:step:ps   tensset=rate:max:T
            if (kv[0] == "tgset")
            {
                var v = kv[1].Split(':').Select(x => decimal.Parse(x, System.Globalization.CultureInfo.InvariantCulture)).ToArray();
                (w.ViewModel.Analyze.TgFromD, w.ViewModel.Analyze.TgToD, w.ViewModel.Analyze.TgStepD, w.ViewModel.Analyze.TgPsD) = (v[0], v[1], v[2], v[3]);
            }
            if (kv[0] == "tensset")
            {
                var v = kv[1].Split(':').Select(x => decimal.Parse(x, System.Globalization.CultureInfo.InvariantCulture)).ToArray();
                (w.ViewModel.Analyze.TensRateD, w.ViewModel.Analyze.TensMaxD, w.ViewModel.Analyze.TensTD) = (v[0], v[1], v[2]);
            }
            if (kv[0] == "wrap") w.ViewModel.Wrap = kv[1] == "1";
            if (kv[0] == "tab") w.SelectAnalysisTab(int.Parse(kv[1]));
            if (kv[0] == "view") w.ViewModel.ViewBackground = int.Parse(kv[1]);
        }
        Pump();
        var frame = w.CaptureRenderedFrame();
        frame?.Save(output);
        Console.WriteLine(frame == null ? "no frame captured" : $"wrote {output}");
        return frame == null ? 1 : 0;
    }

    private static void Pump()
    {
        for (var i = 0; i < 40; i++)
        {
            Dispatcher.UIThread.RunJobs();
            AvaloniaHeadlessPlatform.ForceRenderTimerTick();
            Thread.Sleep(30);
        }
    }
}
