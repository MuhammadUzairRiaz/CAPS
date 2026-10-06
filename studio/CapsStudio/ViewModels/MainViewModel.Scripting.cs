using System.Collections.ObjectModel;
using System.Diagnostics;
using CapsStudio.Interop;

namespace CapsStudio.ViewModels;

/// <summary>A user shortcut that shadows one of the Studio's own keys (Settings › Accessibility).</summary>
public sealed record ShortcutConflict(string Id, string Gesture, string Text);

/// <summary>Settings › Accessibility and Settings › Python &amp; scripting (design/boards/Settings).</summary>
public sealed partial class MainViewModel
{
    // ---------------------------------------------------------------- Python & scripting
    /// <summary>The interpreter macros, render overlays and Python pipeline steps run with: the one set here, else
    /// CAPS_PYTHON, else python3 on the PATH.</summary>
    public static string PythonExe => _pythonSetting is { Length: > 0 } p ? p : Environment.GetEnvironmentVariable("CAPS_PYTHON") is { Length: > 0 } e ? e : "python3";
    private static string _pythonSetting = "";
    public string SetPythonExe
    {
        get => _settings.PythonExe;
        set
        {
            var v = (value ?? "").Trim();
            if (_settings.PythonExe == v) return;
            _settings.PythonExe = v;
            Raise(); Raise(nameof(PythonSource)); ApplyPython(); Changed("Python interpreter");
            PythonCheck = "";
        }
    }
    public string PythonSource => _settings.PythonExe.Length > 0 ? "set here" : Environment.GetEnvironmentVariable("CAPS_PYTHON") is { Length: > 0 } ? "from CAPS_PYTHON" : "python3 on the PATH";
    /// <summary>The core runs Python pipeline steps with the same interpreter.</summary>
    private void ApplyPython()
    {
        _pythonSetting = _settings.PythonExe;
        try { if (Paths.Python is { } py) Native.SetPython(py, PythonExe); } catch { /* an older core */ }
    }
    private string _pythonCheck = "";
    public string PythonCheck { get => _pythonCheck; private set => Set(ref _pythonCheck, value); }
    private bool _pythonChecking;
    public bool PythonChecking { get => _pythonChecking; private set => Set(ref _pythonChecking, value); }

    /// <summary>Runs the interpreter: its version, NumPy (optional), and the caps package with this Studio's library.</summary>
    public async Task CheckPython()
    {
        PythonChecking = true;
        try
        {
            var exe = PythonExe;
            var code = "import sys, importlib.util as u\nprint('Python ' + sys.version.split()[0] + ' · ' + sys.executable)\n" +
                       "print('NumPy ' + (__import__('numpy').__version__ if u.find_spec('numpy') else 'not installed (optional: arrays from caps.Table)'))\n" +
                       "import caps\nprint('caps package: ABI %d · %s' % (caps.abi_version(), caps.__file__))\n";
            var (ok, text) = await Task.Run(() =>
            {
                var psi = PythonProcess.Utf8Io(new ProcessStartInfo(exe) { RedirectStandardOutput = true, RedirectStandardError = true, UseShellExecute = false });
                psi.ArgumentList.Add("-c");
                psi.ArgumentList.Add(code);
                if (Paths.Python is { } pkg) psi.Environment["PYTHONPATH"] = pkg + (Environment.GetEnvironmentVariable("PYTHONPATH") is { Length: > 0 } pp ? Path.PathSeparator + pp : "");
                if (Native.LibraryPath is { } lib) psi.Environment["CAPS_LIB"] = lib;
                try
                {
                    using var p = Process.Start(psi)!;
                    var err = p.StandardError.ReadToEndAsync();
                    var output = p.StandardOutput.ReadToEnd();
                    if (!p.WaitForExit(20000)) { try { p.Kill(true); } catch { } return (false, "the interpreter did not answer within 20 s"); }
                    return (p.ExitCode == 0, (output + (p.ExitCode != 0 ? err.Result.Trim().Split('\n').LastOrDefault() ?? "" : "")).Trim());
                }
                catch (Exception e) { return (false, $"cannot run {exe}: {e.Message}"); }
            });
            PythonCheck = (ok ? "✓ " : "✗ ") + text.Replace("\n", "\n   ");
        }
        finally { PythonChecking = false; }
    }

    /// <summary>The shell lines that give your own scripts and notebooks this Studio's caps package and library.</summary>
    public string PythonShellLines
    {
        get
        {
            var pkg = Paths.Python ?? "<the Studio's data/python>";
            var lib = Native.LibraryPath ?? "<the Studio's libcaps>";
            return OperatingSystem.IsWindows()
                ? $"set PYTHONPATH={pkg};%PYTHONPATH%\nset CAPS_LIB={lib}"
                : $"export PYTHONPATH=\"{pkg}:$PYTHONPATH\"\nexport CAPS_LIB=\"{lib}\"";
        }
    }

    // ---------------------------------------------------------------- Accessibility
    public bool SetHighContrast
    {
        get => _settings.HighContrast;
        set { if (_settings.HighContrast == value) return; _settings.HighContrast = value; Tokens.UseHighContrast(value); Raise(); RenderRequested?.Invoke(); MolViewChanged?.Invoke(); Changed("High-contrast outlines"); }
    }
    public bool SetFocusRings
    {
        get => _settings.FocusRingsAlways;
        set { if (_settings.FocusRingsAlways == value) return; _settings.FocusRingsAlways = value; Tokens.UseFocusRings(value); Raise(); Changed("Focus rings"); }
    }
    public bool SetAnnounceProgress
    {
        get => _settings.AnnounceProgress;
        set { if (_settings.AnnounceProgress == value) return; _settings.AnnounceProgress = value; Raise(); Changed("Announce job progress"); }
    }
    /// <summary>The view's outline level for the renderer: 0 none, 1 light, 2 high contrast.</summary>
    private int OutlineLevel => _outlines ? (_settings.HighContrast ? 2 : 1) : 0;

    private readonly Dictionary<Job, int> _announcedQuarter = new();
    /// <summary>With Announce job progress on, a running job is read out when it passes a quarter of its run.</summary>
    private void AnnounceProgress(Job? j)
    {
        if (j == null || !_settings.AnnounceProgress || j.Status is not ("running" or "")) return;
        var q = (int)Math.Floor(Math.Clamp(j.Progress, 0, 1) * 4);
        if (q <= 0 || q >= 4) return;
        if (_announcedQuarter.TryGetValue(j, out var had) && had >= q) return;
        _announcedQuarter[j] = q;
        Announcement = $"{j.Title}: {q * 25} % done.";
    }

    public ObservableCollection<ShortcutConflict> ShortcutConflicts { get; } = new();
    public bool HasShortcutConflicts => ShortcutConflicts.Count > 0;
    /// <summary>Your shortcuts that replace one of the Studio's own keys.</summary>
    public void FillShortcutConflicts()
    {
        if (!_modelCommands) { _modelCommands = true; AddModelCommands(); }
        ShortcutConflicts.Clear();
        foreach (var (id, g) in _settings.Shortcuts)
            foreach (var (bg, what) in BuiltInKeys)
                if ((OperatingSystem.IsMacOS() ? bg : bg.Replace("Meta", "Ctrl")) == g)
                    ShortcutConflicts.Add(new ShortcutConflict(id, g, $"{ShowGesture(g)} runs “{_commands.FirstOrDefault(c => c.Id == id)?.Title ?? id}” and no longer {what}"));
        Raise(nameof(HasShortcutConflicts));
    }
    /// <summary>Resolve: your key goes, the Studio's own comes back.</summary>
    public void ResolveShortcutConflict(ShortcutConflict c)
    {
        if (_settings.Shortcuts.Remove(c.Id)) _settings.Save();
        FillShortcutConflicts();
        FillShortcutRows();
        Status = $"{ShowGesture(c.Gesture)} is the Studio's own key again";
    }
}
