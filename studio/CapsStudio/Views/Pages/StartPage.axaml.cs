using Avalonia.Controls;
using Avalonia.Input;
using Avalonia.Interactivity;
using Avalonia.Markup.Xaml;
using Avalonia.Platform.Storage;
using CapsStudio.ViewModels;

namespace CapsStudio.Views.Pages;

public partial class StartPage : PageBase
{
    public StartPage()
    {
        AvaloniaXamlLoader.Load(this);
        AddHandler(DragDrop.DropEvent, (_, e) =>
        {
            if (e.Data.GetFiles() is { } files) Window?.OpenMany(files.Select(f => f.TryGetLocalPath()).OfType<string>().ToList());
        });
    }

    private void OnQuickKey(object? s, KeyEventArgs e)
    {
        if (e.Key == Key.Enter) { OnQuickGo(s, e); e.Handled = true; }
    }

    private void OnQuickGo(object? s, RoutedEventArgs e)
    {
        if (Vm.QuickGo() is { } path) Window?.OpenMany([path]);
    }

    private void OnBuilder(object? s, RoutedEventArgs e)
    {
        if ((s as Control)?.Tag is StartBuilder b && b.Available)
        {
            if (b.Module == 9) Vm.OpenBuilder();
            else if (b.Module == 14) Vm.OpenSurface();
            else if (b.Module == 29) Vm.OpenCrystal();
            else if (b.Module == 31) Vm.OpenSolvation();
            else Vm.SetModule(b.Module);
        }
    }

    private void OnRestoreSession(object? s, RoutedEventArgs e) => Vm.RestoreSession();
    private void OnForgetSession(object? s, RoutedEventArgs e) => Vm.ForgetSession();

    private void OnRecent(object? s, RoutedEventArgs e)
    {
        if ((s as Control)?.Tag is not RecentItem r) return;
        if (!File.Exists(r.Path)) { Vm.Status = $"{r.Name} is no longer at {r.Path}"; RecentFiles.Forget(r.Path); Vm.LoadRecent(); return; }
        Window?.OpenMany(r.Topology != null && File.Exists(r.Topology) ? [r.Path, r.Topology] : [r.Path]);
    }

    private void OnGuide(object? s, RoutedEventArgs e)
    {
        if ((s as Control)?.Tag is not StartGuide g) return;
        if (g.Sample && !Vm.HasDocument) Window?.OnOpenSample(s, e);
        Vm.StartGuideAction(g);
    }

    private async void OnRecipe(object? s, RoutedEventArgs e)
    {
        if (Window == null) return;
        await Window.RecipeOpenAsync();
    }

    private async void OnOpen(object? s, RoutedEventArgs e) { if (Window != null) await Window.OpenDialog(); }
    private void OnSample(object? s, RoutedEventArgs e) => Window?.OnOpenSample(s, e);
    private void OnSettings(object? s, RoutedEventArgs e) => Window?.ShowSettings();
}
