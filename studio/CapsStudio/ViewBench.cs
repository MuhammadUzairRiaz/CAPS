using Avalonia;
using Avalonia.Threading;
using CapsStudio.Views;

namespace CapsStudio;

/// <summary>--viewbench FILE [SECONDS]: the real Studio (on screen, OpenGL where available) opens FILE, turns the main
/// view for SECONDS, then opens the full view and turns it as long, and prints each view's frame times and the frames a
/// second it delivered (a third argument "cpu": the full view on the CPU renderer, for comparison). A check of view
/// speed on the user's own machine; nothing is saved.</summary>
internal static class ViewBench
{
    public static double Seconds { get; private set; }
    public static bool Active => Seconds > 0;
    private static bool _cpu;

    public static int Run(string[] args)
    {
        if (args.Length < 1) { Console.Error.WriteLine("usage: CapsStudio --viewbench FILE [SECONDS]"); return 2; }
        _cpu = args.Length > 2 && args[2] == "cpu";
        Seconds = args.Length > 1 ? double.Parse(args[1], System.Globalization.CultureInfo.InvariantCulture) : 4;
        // the user's settings and recent list are left alone
        ViewModels.AppSettings.Override = Path.Combine(Path.GetTempPath(), "caps-viewbench-settings.json");
        ViewModels.RecentFiles.Override = Path.Combine(Path.GetTempPath(), "caps-viewbench-recent");
        return Program.BuildAvaloniaApp().StartWithClassicDesktopLifetime([args[0]]);
    }

    public static async void Start(MainWindow w, Action quit)
    {
        try
        {
            var t0 = DateTime.UtcNow;
            while (!w.BenchReady(true) && (DateTime.UtcNow - t0).TotalSeconds < 60) await Task.Delay(100);
            if (!w.BenchReady(false)) { Console.WriteLine("viewbench: the structure did not open"); quit(); return; }
            await Task.Delay(1000);
            var gl = w.BenchGl;
            Console.WriteLine($"viewbench: {w.ViewModel.Document!.Summary().Atoms} atoms · main view GPU {(gl.Ready ? "on · " + gl.Status : "off")}");

            // the main view
            var main = new List<double>();
            w.BenchFrames = main;
            var gl0 = gl.FramesDrawn;
            var secs = await Drive(() => w.BenchTurn(0.02));
            var glFrames = gl.FramesDrawn - gl0;
            w.BenchFrames = null;
            Report("main view", w.BenchPixels, main, gl.Ready ? glFrames : main.Count, secs);

            // the full view
            var fv = w.OpenFullView(false)!;
            await Task.Delay(1000);
            fv.Activate();
            fv.WindowState = Avalonia.Controls.WindowState.FullScreen;
            await Task.Delay(3000);   // full screen (macOS animates it)
            if (fv.Screens.ScreenFromWindow(fv) is { } scr && fv.Bounds.Width < scr.Bounds.Width / scr.Scaling - 1)
            {
                // not in full screen (the app was not in front): the window over the whole screen instead
                fv.WindowState = Avalonia.Controls.WindowState.Normal;
                fv.Position = scr.WorkingArea.Position;
                fv.Width = scr.WorkingArea.Width / scr.Scaling;
                fv.Height = scr.WorkingArea.Height / scr.Scaling;
                await Task.Delay(1500);
            }
            var view = fv.BenchView;
            Console.WriteLine($"viewbench: full view window {fv.WindowState} {fv.Bounds.Width:0}×{fv.Bounds.Height:0} pt, view {view.Bounds.Width:0}×{view.Bounds.Height:0} pt");
            if (_cpu) { view.UseGpu = false; await Task.Delay(1000); }
            var full = new List<double>();
            view.FrameDone += full.Add;
            var f0 = view.Frames;
            var g0 = view.GpuFramesDrawn;
            secs = await Drive(() => { var c = view.Camera; c.Yaw += 0.02; view.Camera = c; });
            await Task.Delay(300);
            var delivered = view.UsingGpu ? view.GpuFramesDrawn - g0 : view.Frames - f0;
            Report($"full view ({(view.UsingGpu ? "GPU" : "CPU")})", view.FramePixels, full, delivered, secs);
            fv.Close();
        }
        catch (Exception e) { Console.WriteLine("viewbench failed: " + e); }
        quit();
    }

    /// <summary>Turns the camera at 125 Hz (faster than any display) for the benchmark's seconds.</summary>
    private static async Task<double> Drive(Action turn)
    {
        var sw = System.Diagnostics.Stopwatch.StartNew();
        var tcs = new TaskCompletionSource();
        var timer = new DispatcherTimer(TimeSpan.FromMilliseconds(8), DispatcherPriority.Input, (_, _) => { });
        timer.Tick += (_, _) =>
        {
            if (sw.Elapsed.TotalSeconds >= Seconds) { timer.Stop(); tcs.TrySetResult(); return; }
            turn();
        };
        timer.Start();
        await tcs.Task;
        return sw.Elapsed.TotalSeconds;
    }

    private static void Report(string what, (int W, int H) px, List<double> ms, long frames, double secs)
    {
        if (ms.Count == 0) { Console.WriteLine($"viewbench: {what}: no frames"); return; }
        var s = ms.OrderBy(x => x).ToArray();
        double Q(double q) => s[Math.Min(s.Length - 1, (int)(q * s.Length))];
        Console.WriteLine($"viewbench: {what}: {px.W}×{px.H} px · frame median {Q(0.5):0.0} ms · mean {ms.Average():0.0} ms · p95 {Q(0.95):0.0} ms" +
                          $" · {frames / secs:0.0} frames/s delivered ({frames} in {secs:0.0} s, {ms.Count} camera updates)");
    }
}
