using Avalonia;
using Avalonia.Controls;
using Avalonia.Interactivity;
using Avalonia.Markup.Xaml;
using Avalonia.Platform.Storage;
using CapsStudio.Interop;
using CapsStudio.ViewModels;

namespace CapsStudio.Views.Pages;

public partial class DensityCalcPage : PageBase
{
    public DensityCalcPage() { AvaloniaXamlLoader.Load(this); }
    private void OnSend(object? s, RoutedEventArgs e) => Vm.SendDensityToPack();
    private void OnAdd(object? s, RoutedEventArgs e) => Vm.AddSpeciesFromSmiles();
}
