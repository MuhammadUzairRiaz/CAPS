using System.Diagnostics;
using System.Runtime.InteropServices;

namespace CapsStudio;

/// <summary>A short system sound for a selection (Settings › Accessibility "Sound on selection", design/boards/VisAccess):
/// the platform's own sound, played by its own tool and never waited for — macOS afplay with a system sound, Windows the
/// system's asterisk sound, Linux canberra-gtk-play (or nothing when it is not there). At most one every 120 ms.</summary>
public static partial class SoundCue
{
    private static long _last;

    [LibraryImport("user32.dll")]
    [return: MarshalAs(UnmanagedType.Bool)]
    private static partial bool MessageBeep(uint type);   // 0x40: MB_ICONASTERISK, the system's asterisk sound

    public static void Selection()
    {
        var now = Environment.TickCount64;
        if (now - Interlocked.Read(ref _last) < 120) return;
        Interlocked.Exchange(ref _last, now);
        try
        {
            if (OperatingSystem.IsWindows()) { MessageBeep(0x40); return; }
            ProcessStartInfo psi;
            if (OperatingSystem.IsMacOS())
            {
                const string tink = "/System/Library/Sounds/Tink.aiff";
                if (!File.Exists(tink)) return;
                psi = new ProcessStartInfo("afplay") { ArgumentList = { "-v", "0.4", tink } };
            }
            else psi = new ProcessStartInfo("canberra-gtk-play") { ArgumentList = { "--id", "message" } };
            psi.UseShellExecute = false;
            psi.CreateNoWindow = true;
            psi.RedirectStandardError = true;
            psi.RedirectStandardOutput = true;
            using var p = Process.Start(psi);
        }
        catch { /* no player on this system: silence */ }
    }
}
