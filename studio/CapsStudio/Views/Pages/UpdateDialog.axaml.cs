using Avalonia.Controls;
using Avalonia.Input;
using Avalonia.Interactivity;
using Avalonia.Markup.Xaml;

namespace CapsStudio.Views.Pages;

public partial class UpdateDialog : PageBase
{
    public UpdateDialog() { AvaloniaXamlLoader.Load(this); }

    private void OnBackdrop(object? s, PointerPressedEventArgs e) => Vm.CloseUpdate();
    private void OnLater(object? s, RoutedEventArgs e) => Vm.CloseUpdate();
    private async void OnChangelog(object? s, RoutedEventArgs e)
    {
        if (TopLevel.GetTopLevel(this)?.Launcher is { } l) await l.LaunchUriAsync(new System.Uri("https://github.com/MuhammadUzairRiaz/CAPS/releases"));
    }
    private async void OnDownload(object? s, RoutedEventArgs e)
    {
        if (Vm.UpdateUrl.Length > 0 && TopLevel.GetTopLevel(this)?.Launcher is { } l) await l.LaunchUriAsync(new System.Uri(Vm.UpdateUrl));
        Vm.CloseUpdate();
    }
}
