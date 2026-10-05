using Avalonia;
using Avalonia.Controls;
using Avalonia.Controls.ApplicationLifetimes;
using Avalonia.Markup.Xaml;
using Avalonia.Platform.Storage;
using CapsStudio.Views;

namespace CapsStudio;

public partial class App : Application
{
    public override void Initialize()
    {
        AvaloniaXamlLoader.Load(this);
        GuardNumberBoxes();
    }

    private static bool _numberGuard;
    /// <summary>Every number box in CAPS: a value typed beyond its range is brought to the nearest end, and a box emptied
    /// (or holding text that is no number) goes back to its last value — the page never receives "nothing" for a number.</summary>
    public static void GuardNumberBoxes()
    {
        if (_numberGuard) return;
        _numberGuard = true;
        NumericUpDown.ValueProperty.Changed.AddClassHandler<NumericUpDown>((box, e) =>
        {
            if (e.NewValue is not null || e.OldValue is not decimal old) return;
            var back = Math.Clamp(old, box.Minimum, Math.Max(box.Minimum, box.Maximum));
            Avalonia.Threading.Dispatcher.UIThread.Post(() => { if (box.Value is null) box.Value = back; });
        });
    }

    private MainWindow? _main;
    private void OnAbout(object? s, System.EventArgs e) => _main?.ShowSettings();
    private void OnSettings(object? s, System.EventArgs e) => _main?.ShowSettings();

    public override void OnFrameworkInitializationCompleted()
    {
        if (ApplicationLifetime is IClassicDesktopStyleApplicationLifetime desktop)
        {
            // an error in the interface is written down and reported, not the end of the session (the user's structures stay)
            Avalonia.Threading.Dispatcher.UIThread.UnhandledException += (_, e) =>
            {
                var where = CrashLog.Write(e.Exception);
                e.Handled = true;
                if (_main is MainWindow mw) mw.ViewModel.Status = $"Something went wrong ({e.Exception.GetType().Name}: {e.Exception.Message}) · CAPS kept running · details in {where}";
            };
            AppDomain.CurrentDomain.UnhandledException += (_, e) => { if (e.ExceptionObject is Exception x) CrashLog.Write(x); };
            TaskScheduler.UnobservedTaskException += (_, e) => { CrashLog.Write(e.Exception); e.SetObserved(); };
            var w = new MainWindow();
            desktop.MainWindow = w;
            _main = w;
            var args = desktop.Args ?? [];
            if (args.Length > 0) w.OpenOnStart(args[0], args.Length > 1 ? args[1] : null);
            // macOS hands a double-clicked file (a .capsproj, a structure) to the running app, not as an argument
            if (TryGetFeature(typeof(IActivatableLifetime)) is IActivatableLifetime act)
                act.Activated += (_, e) =>
                {
                    if (e is FileActivatedEventArgs f)
                        foreach (var p in f.Files.Select(x => x.TryGetLocalPath()).OfType<string>().Take(1)) w.OpenOnStart(p, null);
                };
        }
        base.OnFrameworkInitializationCompleted();
    }
}

/// <summary>Errors CAPS could not handle, with their stack, in the user's application data (macOS ~/Library/Application
/// Support/CAPS/logs, Windows %LOCALAPPDATA%\CAPS\logs, Linux ~/.local/share/CAPS/logs): crash.log, the newest last.</summary>
public static class CrashLog
{
    public static string Write(Exception x)
    {
        try
        {
            var dir = Path.Combine(Environment.GetFolderPath(Environment.SpecialFolder.LocalApplicationData), "CAPS", "logs");
            Directory.CreateDirectory(dir);
            var path = Path.Combine(dir, "crash.log");
            if (File.Exists(path) && new FileInfo(path).Length > 2_000_000) File.Delete(path);
            File.AppendAllText(path, $"---- {DateTime.Now:yyyy-MM-dd HH:mm:ss} · CAPS {typeof(CrashLog).Assembly.GetName().Version}\n{x}\n\n");
            return path;
        }
        catch { return "(no log: the folder is not writable)"; }
    }
}
