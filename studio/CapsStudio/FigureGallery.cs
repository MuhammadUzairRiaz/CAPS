using Avalonia;
using Avalonia.Headless;
using Avalonia.Media.Imaging;
using CapsStudio.ViewModels;
using CapsStudio.Views;

namespace CapsStudio;

/// <summary>CapsStudio --figures OUTDIR [SAMPLES]: short real runs on the samples (each MD run at most 5,000 steps) and every
/// figure they give written the way the Studio writes them for print (PDF, SVG and a 300 dpi PNG, one journal column
/// wide) and as the screen draws them (name_screen.png), to look at side by side. Nothing is made up: each figure holds
/// the rows the run produced.</summary>
public static class FigureGallery
{
    public static int Run(string[] args)
    {
        var outDir = Path.GetFullPath(args.Length > 0 ? args[0] : "figures");
        var dir = args.Length > 1 ? args[1] : "samples";
        Directory.CreateDirectory(outDir);
        Environment.SetEnvironmentVariable("CAPS_NO_TOUR", "1");
        var settings = Path.Combine(outDir, "settings.json");
        AppSettings.Override = settings;
        RecentFiles.Override = Path.Combine(outDir, "recent");
        var vm = new MainViewModel();
        var written = new List<string>();
        var figs = new List<(string Name, SciFigure? Fig)>();
        // the runs first (their view model posts to a dispatcher that has no window yet, as in the self-test), the drawing after
        void Write(string name, SciFigure? fig) { figs.Add((name, fig)); Console.WriteLine($"ran {name}"); }
        // the run's rows reach the view model through the dispatcher: let it deliver them
        void Pump() { for (var i = 0; i < 5; ++i) { Avalonia.Threading.Dispatcher.UIThread.RunJobs(); Thread.Sleep(20); } }
        void Draw(string name, SciFigure? fig, double aspect = 0.75)
        {
            if (fig == null || fig.Series.Length == 0 && fig.Heat == null) { Console.WriteLine($"skip {name}: no data"); return; }
            var h = SciExport.Single * aspect;
            foreach (var ext in new[] { ".pdf", ".svg", ".png" }) SciExport.Save(fig, Path.Combine(outDir, name + ext), SciExport.Single, h);
            // the screen version (the theme the Studio starts in), 560 × 320 px at 2× pixels
            using (var rtb = new RenderTargetBitmap(new PixelSize(1120, 640), new Vector(96, 96)))
            {
                using (var ctx = rtb.CreateDrawingContext())
                {
                    ctx.FillRectangle(Tokens.Brush("Bg1B"), new Rect(0, 0, 1120, 640));
                    using var t = ctx.PushTransform(Matrix.CreateScale(2, 2));
                    var st = SciStyle.Screen();
                    SciPlot.DrawTo(new ScreenSink(ctx, st.Font), fig, new Size(560, 320), st);
                }
                rtb.Save(Path.Combine(outDir, name + "_screen.png"));
            }
            var n = fig.Series.Sum(s => s.X.Length);
            written.Add($"{name}: {fig.Series.Length} series, {n} points, kinds {string.Join("/", fig.Series.Select(s => s.Kind).Distinct())}");
            Console.WriteLine("wrote " + written[^1]);
        }

        // Analyze on the sample trajectory: density, g(r), Rg and MSD over its frames
        vm.Open(Path.Combine(dir, "ps_melt.lammpstrj"), Path.Combine(dir, "ps_melt.data"));
        foreach (var c in vm.Analyze.Groups.SelectMany(g => g.Chips)) c.IsOn = c.Id is "rdf" or "rg" or "msd";
        vm.Analyze.FramePsD = 1;
        vm.Analyze.Run().GetAwaiter().GetResult();
        foreach (var c in vm.Analyze.Curves.Where(c => c.Figure != null).GroupBy(c => c.Property).Select(g => g.First()))
            Write("analyze_" + c.Property, c.Figure);

        // Minimise the melt (energy and largest force per reported iteration)
        vm.Open(Path.Combine(dir, "ps_melt.data"));
        vm.RelaxFtolD = 2;
        vm.Relax().GetAwaiter().GetResult(); Pump();
        Write("relax_energy", PageFigures.RelaxEnergy(vm.RelaxEnergyCurve));
        Write("relax_force", PageFigures.RelaxForce(vm.RelaxForceCurve, 2));

        // Dynamics: 4,000 steps of NVT at 1 fs (Studio sampling: every steps/400, at least 10)
        vm.MdStepsD = 4000;
        vm.RunMd().GetAwaiter().GetResult(); Pump();
        Console.WriteLine($"md: {vm.Thermo.Count} thermo rows · {vm.MdLog.Split('\n')[0]}");
        Write("md_temperature", PageFigures.ThermoTemperature(vm.Thermo, (double)(vm.MdTempD ?? 300m)));
        Write("md_energy", PageFigures.ThermoEnergy(vm.Thermo));

        // a short tensile protocol (0.02 strain at 0.01 /ps: 2,000 steps) and its stress–strain figure, as Mechanics shows it
        vm.Open(Path.Combine(dir, "ps_melt.data"));
        foreach (var c in vm.Analyze.Groups.SelectMany(g => g.Chips)) c.IsOn = c.Id is "tensile";
        vm.Analyze.TensRateD = 0.01m;
        vm.Analyze.TensMaxD = 0.02m;
        vm.Analyze.TensTD = 300;
        vm.Analyze.EqPsD = 1;
        vm.Analyze.Run().GetAwaiter().GetResult();
        if (vm.Analyze.Curves.FirstOrDefault(c => c.Label.Contains("smoothed")) is { } ss)
        {
            Write("analyze_tensile", ss.Figure);
            var curve = ss.X.Zip(ss.Y, (x, y) => (x * 100, y / 1000)).ToArray();
            var smooth = ss.OverlayX is { } ox && ss.OverlayY is { } oy ? ox.Zip(oy, (x, y) => (x * 100, y / 1000)).ToArray() : [];
            var e = vm.Analyze.Results.FirstOrDefault(r => r.Id == "tensile_modulus")?.Value ?? double.NaN;
            Write("mechanics", PageFigures.Mechanics(curve, double.IsFinite(e) ? [(0, 0), (3, e * 0.03)] : [], 2, smooth));
        }

        // DPD: an A5B5 diblock melt, 2,000 steps
        vm.OpenDpd();
        vm.DpdSpecies.Clear();
        vm.DpdSpecies.Add(new DpdSpeciesRow { Name = "diblock", Sequence = "A5B5", Count = 80 });
        vm.DpdSteps = 2000;
        vm.RunDpd().GetAwaiter().GetResult(); Pump();
        var o = vm.DpdOrderSeries;
        Write("dpd_order", new SciFigure
        {
            XLabel = "step", YLabel = "*ψ*", YMin = 0, YZero = true, Empty = "",
            Series = [new SciSeries("", o.Select(p => p.X).ToArray(), o.Select(p => p.Y).ToArray()) { Kind = PageFigures.Auto(o.Length), Width = 1.3 }],
        });

        // Pack: four maleic anhydride molecules into the melt
        vm.Open(Path.Combine(dir, "ps_melt.data"));
        vm.SetModule(5);
        vm.PackStart = 1;
        vm.NewPackInput();
        vm.AddPackMolecule("O=C1OC(=O)C=C1", "maleic anhydride").GetAwaiter().GetResult();
        vm.SetPackRowCount(vm.PackItems.Count - 1, 300);
        vm.PackAssignField = false;
        vm.RunPack().GetAwaiter().GetResult(); Pump();
        Console.WriteLine($"pack: {vm.PackLog} · rounds {vm.PackLoop} · {vm.PackCurve.Count} points");
        Write("pack_penalty", PageFigures.PackPenalty(vm.PackCurve));

        // thermal conductivity: Green–Kubo on a small PE lamella, 0.5 ps (500 steps) after 0.05 ps
        vm.OpenSemicrystalline();
        vm.ScNx = 3; vm.ScNy = 5; vm.ScCrystal = 1.5m; vm.ScAmorphous = 2.5m;
        vm.BuildSemicrystalline().GetAwaiter().GetResult();
        vm.OpenThermal();
        vm.ThDim = 0; vm.ThEq = 0.05m; vm.ThPs = 0.5m; vm.ThCorr = 0.1m;
        vm.RunThermal().GetAwaiter().GetResult(); Pump();
        Write("thermal_gk", PageFigures.Thermal(vm.ThCurves, true));

        AppBuilder.Configure<App>().UseSkia().UseHeadless(new AvaloniaHeadlessPlatformOptions { UseHeadlessDrawing = false }).WithCapsFonts().SetupWithoutStarting();
        foreach (var (name, fig) in figs) Draw(name, fig);
        // a two-column page of four of them, as the figure composer writes it (vector PDF and SVG, a 300 dpi PNG)
        var page = new FigureComposition { WidthIn = 7.0, HeightIn = 4.6, Rows = 2, Cols = 2 };
        foreach (var name in new[] { "md_temperature", "relax_force", "thermal_gk", "mechanics" })
            if (figs.FirstOrDefault(f => f.Name == name).Fig is { } pf)
                page.Panels.Add(new FigurePanel { Kind = "curve", Figure = pf, Source = name, Colour = Avalonia.Media.Color.Parse(SciStyle.PaperPalette[0]) });
        File.WriteAllBytes(Path.Combine(outDir, "composer_page.pdf"), page.ToPdf());
        File.WriteAllText(Path.Combine(outDir, "composer_page.svg"), page.ToSvg(_ => null));
        using (var rtb = new RenderTargetBitmap(new PixelSize(2100, 1380), new Vector(96, 96)))
        {
            using (var ctx = rtb.CreateDrawingContext()) page.Draw(ctx, 300 / 72.0);
            rtb.Save(Path.Combine(outDir, "composer_page.png"));
        }
        written.Add("composer_page: 4 panels, 7 × 4.6 in");
        File.WriteAllLines(Path.Combine(outDir, "figures.txt"), written);
        Console.WriteLine($"{written.Count} figures in {outDir}");
        return 0;
    }
}
