using Avalonia.Controls;
using Avalonia.Interactivity;
using Avalonia.Markup.Xaml;
using CapsStudio.ViewModels;

namespace CapsStudio.Views.Pages;

public partial class SweepPage : PageBase
{
    public SweepPage() { AvaloniaXamlLoader.Load(this); }

    private async void OnRun(object? s, RoutedEventArgs e) => await Vm.RunSweep();
    private void OnPause(object? s, RoutedEventArgs e) => Vm.PauseSweep();
    private void OnStop(object? s, RoutedEventArgs e) => Vm.StopSweep();
    private void OnAddSeeds(object? s, RoutedEventArgs e) => Vm.AddSweepSeeds();
    private void OnRunSquare(object? s, RoutedEventArgs e) { if (s is Button { Tag: SweepRun r }) Vm.OpenSweepRun(r); }
}
