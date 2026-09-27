using Avalonia;
using Avalonia.Controls;
using Avalonia.Controls.ApplicationLifetimes;
using Avalonia.Markup.Xaml;
using CapsStudio.Views;

namespace CapsStudio;

public partial class App : Application
{
    public override void Initialize() => AvaloniaXamlLoader.Load(this);

    public override void OnFrameworkInitializationCompleted()
    {
        if (ApplicationLifetime is IClassicDesktopStyleApplicationLifetime desktop)
        {
            var w = new MainWindow();
            desktop.MainWindow = w;
            // the application menu (macOS): about, settings; Hide and Quit are the system's
            var app = new NativeMenu();
            var about = new NativeMenuItem("About CAPS Studio");
            about.Click += (_, _) => w.ShowSettings();
            var settings = new NativeMenuItem("Settings…") { Gesture = Avalonia.Input.KeyGesture.Parse("Meta+OemComma") };
            settings.Click += (_, _) => w.ShowSettings();
            app.Add(about);
            app.Add(new NativeMenuItemSeparator());
            app.Add(settings);
            NativeMenu.SetMenu(this, app);
            var args = desktop.Args ?? [];
            if (args.Length > 0) w.OpenOnStart(args[0], args.Length > 1 ? args[1] : null);
        }
        base.OnFrameworkInitializationCompleted();
    }
}
