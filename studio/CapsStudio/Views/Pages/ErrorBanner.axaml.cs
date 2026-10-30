using Avalonia.Interactivity;
using Avalonia.Markup.Xaml;

namespace CapsStudio.Views.Pages;

public partial class ErrorBanner : PageBase
{
    public ErrorBanner() { AvaloniaXamlLoader.Load(this); }

    private void OnSend(object? s, RoutedEventArgs e) => Vm.RunCommand("app.report");
    private void OnDismiss(object? s, RoutedEventArgs e) => Vm.ErrorBannerOpen = false;
}
