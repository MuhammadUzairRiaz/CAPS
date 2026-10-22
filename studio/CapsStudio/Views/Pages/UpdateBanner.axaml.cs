using Avalonia.Controls;
using Avalonia.Interactivity;
using Avalonia.Markup.Xaml;

namespace CapsStudio.Views.Pages;

public partial class UpdateBanner : PageBase
{
    public UpdateBanner() { AvaloniaXamlLoader.Load(this); }

    private void OnWhatsNew(object? s, RoutedEventArgs e) => Vm.OpenRemoteWhatsNew();
    private async void OnDownload(object? s, RoutedEventArgs e)
    {
        if (TopLevel.GetTopLevel(this)?.Launcher is { } l) await l.LaunchUriAsync(new System.Uri(Vm.UpdateDownloadUrl));
        Vm.Status = "Downloading the new CAPS in your browser; install it over this one (projects and settings stay)";
    }
    private void OnLater(object? s, RoutedEventArgs e) => Vm.UpdateLater();
    private void OnSkip(object? s, RoutedEventArgs e) => Vm.SkipThisUpdate();
}
