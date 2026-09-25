using System;
using System.Diagnostics;

namespace CapsStudio;

/// <summary>Motion tokens (design/boards/Motion): short, purposeful, and off when the system asks for reduced motion.
/// Durations: instant 0 (selection, hover) · fast 120 ms (chips, toggles, tooltips) · base 200 ms (drawers, menus,
/// toasts) · slow 320 ms (dialogs, page changes) · camera 450 ms (fly-to, reset view). Easings as CSS cubic-béziers:
/// standard (0.2, 0, 0, 1), enter (0, 0, 0, 1), exit (0.3, 0, 1, 1). Data — plots, energies, trajectories — is never
/// tweened between values.</summary>
public static class Motion
{
    public static readonly TimeSpan Fast = TimeSpan.FromMilliseconds(120), Base = TimeSpan.FromMilliseconds(200),
                                    Slow = TimeSpan.FromMilliseconds(320), Camera = TimeSpan.FromMilliseconds(450);

    /// <summary>"system" follows the operating system, "on" always reduces, "off" never.</summary>
    public static string Mode { get; set; } = "system";
    private static bool? _system;
    public static bool SystemReduces => _system ??= ReadSystem();
    public static bool Reduced => Mode == "on" || (Mode == "system" && SystemReduces);
    public static string SystemText => SystemReduces ? "the system asks for reduced motion" : "the system allows motion";

    /// <summary>Scaled duration: zero when motion is reduced (moves become cuts).</summary>
    public static TimeSpan Of(TimeSpan d) => Reduced ? TimeSpan.Zero : d;

    public static double Standard(double t) => Bezier(0.2, 0, 0, 1, t);
    public static double Enter(double t) => Bezier(0, 0, 0, 1, t);
    public static double Exit(double t) => Bezier(0.3, 0, 1, 1, t);

    /// <summary>y at x = t of the CSS cubic-bézier (0,0) (x1,y1) (x2,y2) (1,1), by Newton on x(s) with bisection fallback.</summary>
    public static double Bezier(double x1, double y1, double x2, double y2, double t)
    {
        t = Math.Clamp(t, 0, 1);
        static double B(double a, double b, double s) => 3 * a * s * (1 - s) * (1 - s) + 3 * b * s * s * (1 - s) + s * s * s;
        static double dB(double a, double b, double s) => 3 * a * (1 - s) * (1 - s) + 6 * (b - a) * s * (1 - s) + 3 * (1 - b) * s * s;
        double u = t;
        for (int i = 0; i < 8; ++i)
        {
            var d = dB(x1, x2, u);
            if (Math.Abs(d) < 1e-6) break;
            u = Math.Clamp(u - (B(x1, x2, u) - t) / d, 0, 1);
        }
        if (Math.Abs(B(x1, x2, u) - t) > 1e-4)
        {
            double lo = 0, hi = 1;
            for (int i = 0; i < 40; ++i) { u = (lo + hi) / 2; if (B(x1, x2, u) < t) lo = u; else hi = u; }
        }
        return B(y1, y2, u);
    }

    private static bool ReadSystem()
    {
        if (Environment.GetEnvironmentVariable("CAPS_REDUCE_MOTION") is { Length: > 0 } env) return env != "0";
        try
        {
            if (OperatingSystem.IsMacOS()) return Run("defaults", "read com.apple.universalaccess reduceMotion") == "1";
            if (OperatingSystem.IsLinux()) return Run("gsettings", "get org.gnome.desktop.interface enable-animations") == "false";
            if (OperatingSystem.IsWindows())
                return Microsoft.Win32.Registry.GetValue(@"HKEY_CURRENT_USER\Control Panel\Desktop\WindowMetrics", "MinAnimate", "1") as string == "0";
        }
        catch { /* no such setting: motion allowed */ }
        return false;
    }

    private static string Run(string exe, string args)
    {
        using var p = Process.Start(new ProcessStartInfo(exe, args) { RedirectStandardOutput = true, RedirectStandardError = true, UseShellExecute = false, CreateNoWindow = true });
        if (p == null) return "";
        var o = p.StandardOutput.ReadToEnd().Trim();
        p.WaitForExit(2000);
        return o;
    }
}
