using System.Diagnostics;
using System.Runtime.Versioning;

namespace CapsStudio;

/// <summary>A Python 3 interpreter for macros, overlays, notebooks and Python pipeline steps when none is set (Settings ›
/// Python &amp; scripting, or CAPS_PYTHON): python3 / python on the PATH — but never Windows' Microsoft Store placeholders
/// (WindowsApps), which only print "Python was not found" — then the py launcher and the registry (PEP 514) on Windows,
/// then the usual conda, Miniforge, python.org and Homebrew folders (an app started from Finder or the Start menu does
/// not see the shell's PATH). Each candidate must run and report Python 3. Searched once per session.</summary>
public static class PythonLocator
{
    private static readonly object Gate = new();
    private static bool _searched;
    private static string? _found;

    /// <summary>The interpreter found, or null when this machine has none.</summary>
    public static string? Detect()
    {
        lock (Gate)
        {
            if (!_searched) { _found = Search(); _searched = true; }
            return _found;
        }
    }

    private static IEnumerable<string> Candidates()
    {
        var win = OperatingSystem.IsWindows();
        var names = win ? new[] { "python.exe", "python3.exe" } : new[] { "python3", "python" };
        foreach (var dir in (Environment.GetEnvironmentVariable("PATH") ?? "").Split(Path.PathSeparator, StringSplitOptions.RemoveEmptyEntries))
        {
            if (win && dir.Contains("WindowsApps", StringComparison.OrdinalIgnoreCase)) continue;   // the Store placeholders
            foreach (var n in names) yield return Path.Combine(dir, n);
        }
        var home = Environment.GetFolderPath(Environment.SpecialFolder.UserProfile);
        var conda = new[] { "anaconda3", "miniconda3", "miniforge3", "mambaforge", "Anaconda3", "Miniconda3" };
        if (win)
        {
            if (PyLauncher() is { } viaPy) yield return viaPy;
            foreach (var r in Registry()) yield return r;
            var local = Environment.GetFolderPath(Environment.SpecialFolder.LocalApplicationData);
            var common = Environment.GetFolderPath(Environment.SpecialFolder.CommonApplicationData);
            foreach (var root in new[] { home, local, common, @"C:\", @"D:\" })
                foreach (var c in conda) yield return Path.Combine(root, c, "python.exe");
            foreach (var root in new[] { Path.Combine(local, "Programs", "Python"), @"C:\" })
                if (Directory.Exists(root))
                    foreach (var d in Directory.EnumerateDirectories(root, "Python3*").OrderByDescending(x => x))
                        yield return Path.Combine(d, "python.exe");
        }
        else
        {
            foreach (var c in conda) yield return Path.Combine(home, c, "bin", "python3");
            foreach (var c in conda) yield return Path.Combine("/opt", c, "bin", "python3");
            yield return "/opt/homebrew/bin/python3";
            yield return "/usr/local/bin/python3";
            yield return "/usr/bin/python3";
        }
        // conda's own list of environments (the base first)
        var envs = Path.Combine(home, ".conda", "environments.txt");
        if (File.Exists(envs))
            foreach (var e in File.ReadAllLines(envs).Where(l => l.Trim().Length > 0))
                yield return win ? Path.Combine(e.Trim(), "python.exe") : Path.Combine(e.Trim(), "bin", "python3");
    }

    private static string? Search()
    {
        var seen = new HashSet<string>(StringComparer.OrdinalIgnoreCase);
        foreach (var c in Candidates())
        {
            try
            {
                if (!seen.Add(c) || !File.Exists(c)) continue;
                // macOS: /usr/bin/python3 is a stub that offers to install the developer tools when they are missing
                if (OperatingSystem.IsMacOS() && c == "/usr/bin/python3" && !Directory.Exists("/Library/Developer/CommandLineTools") && !Directory.Exists("/Applications/Xcode.app")) continue;
                if (Run(c, "-c \"import sys; print(sys.version_info[0])\"")?.Trim() == "3") return c;
            }
            catch { /* the next one */ }
        }
        return null;
    }

    private static string? PyLauncher()
    {
        try { return Run("py", "-3 -c \"import sys; print(sys.executable)\"")?.Trim() is { Length: > 0 } p && File.Exists(p) ? p : null; }
        catch { return null; }
    }

    [SupportedOSPlatform("windows")]
    private static IEnumerable<string> Registry()
    {
        var out_ = new List<string>();
        foreach (var hive in new[] { Microsoft.Win32.Registry.CurrentUser, Microsoft.Win32.Registry.LocalMachine })
        {
            try
            {
                using var py = hive.OpenSubKey(@"Software\Python");
                if (py == null) continue;
                foreach (var company in py.GetSubKeyNames())
                {
                    using var ck = py.OpenSubKey(company);
                    if (ck == null) continue;
                    foreach (var tag in ck.GetSubKeyNames().OrderByDescending(x => x))
                    {
                        using var ip = ck.OpenSubKey(tag + @"\InstallPath");
                        if (ip == null) continue;
                        if (ip.GetValue("ExecutablePath") is string exe && exe.Length > 0) out_.Add(exe);
                        else if (ip.GetValue(null) is string dir && dir.Length > 0) out_.Add(Path.Combine(dir, "python.exe"));
                    }
                }
            }
            catch { /* no access */ }
        }
        return out_;
    }

    private static string? Run(string exe, string args)
    {
        var psi = new ProcessStartInfo(exe, args) { RedirectStandardOutput = true, RedirectStandardError = true, UseShellExecute = false, CreateNoWindow = true };
        using var p = Process.Start(psi);
        if (p == null) return null;
        var text = p.StandardOutput.ReadToEndAsync();
        if (!p.WaitForExit(8000)) { try { p.Kill(true); } catch { } return null; }
        return p.ExitCode == 0 ? text.Result : null;
    }
}
