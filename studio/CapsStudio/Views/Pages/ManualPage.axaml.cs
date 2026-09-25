using Avalonia.Controls;
using Avalonia.Input.Platform;
using Avalonia.Interactivity;
using Avalonia.Markup.Xaml;
using CapsStudio.ViewModels;

namespace CapsStudio.Views.Pages;

public partial class ManualPageView : PageBase
{
    public ManualPageView() { AvaloniaXamlLoader.Load(this); }

    private void OnPage(object? s, RoutedEventArgs e) { if (s is Button { Tag: ManualPage p }) Vm.ShowManualPage(p); }

    private async void OnCopyBibtex(object? s, RoutedEventArgs e)
    {
        if (TopLevel.GetTopLevel(this)?.Clipboard is not { } cb) return;
        await cb.SetTextAsync(Vm.ManualBibtex());
        Vm.Status = $"Copied the BibTeX of {Vm.ManualCurrent?.Title}";
    }
}
