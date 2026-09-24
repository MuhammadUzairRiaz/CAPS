using Avalonia;
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
            var args = desktop.Args ?? [];
            if (args.Length > 0) w.OpenOnStart(args[0], args.Length > 1 ? args[1] : null);
        }
        base.OnFrameworkInitializationCompleted();
    }
}
