using Avalonia.Controls;
using Avalonia.Input;
using Avalonia.Interactivity;
using Avalonia.Markup.Xaml;
using CapsStudio.ViewModels;

namespace CapsStudio.Views.Pages;

public partial class WhatsNewPanel : PageBase
{
    public WhatsNewPanel() { AvaloniaXamlLoader.Load(this); }

    private void OnBackdrop(object? s, PointerPressedEventArgs e) => Vm.CloseWhatsNew();
    private void OnClose(object? s, RoutedEventArgs e) => Vm.CloseWhatsNew();
    private void OnSkip(object? s, RoutedEventArgs e) => Vm.SkipThisUpdate();
    private void OnOpenItem(object? s, RoutedEventArgs e) { if ((s as Control)?.Tag is WhatsNewItem item) Vm.OpenWhatsNewItem(item); }
    private async void OnNotes(object? s, RoutedEventArgs e)
    {
        if (TopLevel.GetTopLevel(this)?.Launcher is { } l) await l.LaunchUriAsync(new System.Uri(Vm.UpdateNotesUrl));
    }
    private async void OnDownload(object? s, RoutedEventArgs e)
    {
        if (TopLevel.GetTopLevel(this)?.Launcher is { } l) await l.LaunchUriAsync(new System.Uri(Vm.UpdateDownloadUrl));
        Vm.CloseWhatsNew();
        Vm.Status = "Downloading the new CAPS in your browser; install it over this one (projects and settings stay)";
    }
}
