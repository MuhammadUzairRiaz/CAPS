using Avalonia;
using Avalonia.Media;

namespace CapsStudio;

internal static partial class Program
{
    [STAThread]
    public static int Main(string[] args)
    {
        // Scientific numbers use a decimal point everywhere (inputs, tables, reports), whatever the system locale.
        System.Globalization.CultureInfo.DefaultThreadCurrentCulture = System.Globalization.CultureInfo.InvariantCulture;
        System.Globalization.CultureInfo.DefaultThreadCurrentUICulture = System.Globalization.CultureInfo.InvariantCulture;
        System.Globalization.CultureInfo.CurrentCulture = System.Globalization.CultureInfo.InvariantCulture;
        if (args.Length > 0 && args[0] is "--selftest" or "--screenshot" or "--shelftest" or "--viewbench" or "--figures")
        {
            Environment.SetEnvironmentVariable("CAPS_NO_TOUR", "1");   // started explicitly there
            ViewModels.UpdateRules.NetworkAllowed = false;            // and nothing there reaches the network (update checks read local files only)
        }
        if (args.Length > 0 && args[0] == "--selftest") return SelfTest.Run(args.Skip(1).ToArray());
        if (args.Length > 0 && args[0] == "--screenshot") return Screenshot.Run(args.Skip(1).ToArray());
        if (args.Length > 0 && args[0] == "--shelftest") return ShelfTest.Run(args.Skip(1).ToArray());
        if (args.Length > 0 && args[0] == "--viewbench") return ViewBench.Run(args.Skip(1).ToArray());
        if (args.Length > 0 && args[0] == "--figures") return FigureGallery.Run(args.Skip(1).ToArray());
        // CAPS_SETTINGS: another settings file (a second copy of the Studio for testing leaves the user's own untouched)
        if (Environment.GetEnvironmentVariable("CAPS_SETTINGS") is { Length: > 0 } sf) ViewModels.AppSettings.Override = sf;
        return AxPress.KeepRunning(() => BuildAvaloniaApp().StartWithClassicDesktopLifetime(args),
            token => Avalonia.Threading.Dispatcher.UIThread.MainLoop(token),
            done => ((Avalonia.Controls.ApplicationLifetimes.IClassicDesktopStyleApplicationLifetime)Application.Current!.ApplicationLifetime!).Exit += (_, e) => done(e.ApplicationExitCode));
    }

    public static AppBuilder BuildAvaloniaApp() =>
        AppBuilder.Configure<App>().UsePlatformDetect().WithCapsFonts().LogToTrace()
            // macOS: the compositor on OpenGL, so the GPU 3D view (GlMolView) shares its context
            .With(new AvaloniaNativePlatformOptions { RenderingMode = [AvaloniaNativeRenderingMode.OpenGl, AvaloniaNativeRenderingMode.Software] });

    /// <summary>IBM Plex Sans (and Mono) from the app's assets as the default family: the design system's type.</summary>
    public static AppBuilder WithCapsFonts(this AppBuilder b) =>
        b.With(new FontManagerOptions { DefaultFamilyName = "avares://CapsStudio/Assets/Fonts#IBM Plex Sans" });
}

/// <summary>An accessibility client (VoiceOver, AX automation) can press a button that is disabled; Avalonia
/// answers with an ElementNotEnabledException. On macOS that exception comes through a native callback, which Avalonia.Native
/// rethrows out of the run loop itself, past the dispatcher's UnhandledException, so it ended the app. A mouse click on a
/// disabled button does nothing, and so does such a press now: the UI loop is entered again until the session exits.</summary>
internal static class AxPress
{
    public static bool Ignorable(Exception x) => x is Avalonia.Automation.ElementNotEnabledException;

    /// <param name="start">starts the app and runs its loop to the end (returns the exit code)</param>
    /// <param name="loop">runs the UI loop again until the token is cancelled</param>
    /// <param name="onExit">registers what to call, with the exit code, when the session exits</param>
    public static int KeepRunning(Func<int> start, Action<CancellationToken> loop, Action<Action<int>> onExit)
    {
        try { return start(); }
        catch (Exception x) when (Ignorable(x)) { }
        var code = 0;
        using var exited = new CancellationTokenSource();
        onExit(c => { code = c; exited.Cancel(); });
        while (!exited.IsCancellationRequested)
        {
            try { loop(exited.Token); break; }   // a loop that returns has ended: the session exited or the dispatcher shut down
            catch (Exception x) when (Ignorable(x)) { }
        }
        Environment.ExitCode = code;
        return code;
    }
}
