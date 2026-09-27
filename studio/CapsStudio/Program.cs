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
        if (args.Length > 0 && args[0] is "--selftest" or "--screenshot") Environment.SetEnvironmentVariable("CAPS_NO_TOUR", "1");   // started explicitly there
        if (args.Length > 0 && args[0] == "--selftest") return SelfTest.Run(args.Skip(1).ToArray());
        if (args.Length > 0 && args[0] == "--screenshot") return Screenshot.Run(args.Skip(1).ToArray());
        return BuildAvaloniaApp().StartWithClassicDesktopLifetime(args);
    }

    public static AppBuilder BuildAvaloniaApp() =>
        AppBuilder.Configure<App>().UsePlatformDetect().WithCapsFonts().LogToTrace()
            // macOS: the compositor on OpenGL, so the GPU 3D view (GlMolView) shares its context
            .With(new AvaloniaNativePlatformOptions { RenderingMode = [AvaloniaNativeRenderingMode.OpenGl, AvaloniaNativeRenderingMode.Software] });

    /// <summary>IBM Plex Sans (and Mono) from the app's assets as the default family: the design system's type.</summary>
    public static AppBuilder WithCapsFonts(this AppBuilder b) =>
        b.With(new FontManagerOptions { DefaultFamilyName = "avares://CapsStudio/Assets/Fonts#IBM Plex Sans" });
}
