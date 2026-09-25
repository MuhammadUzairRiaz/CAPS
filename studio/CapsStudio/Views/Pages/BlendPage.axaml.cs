using Avalonia.Controls;
using Avalonia.Interactivity;
using Avalonia.Markup.Xaml;
using CapsStudio.ViewModels;

namespace CapsStudio.Views.Pages;

public partial class BlendPage : PageBase
{
    public BlendPage() => AvaloniaXamlLoader.Load(this);

    private void OnPolymer(object? s, RoutedEventArgs e) => Vm.SetModule(13);
    private void OnSurface(object? s, RoutedEventArgs e) => Vm.OpenSurface();
    private void OnNano(object? s, RoutedEventArgs e) => Vm.OpenNano();
    private void OnSolvation(object? s, RoutedEventArgs e) => Vm.SetModule(5);
    private void OnCrystal(object? s, RoutedEventArgs e) => Vm.OpenCrystal();
    private void OnBio(object? s, RoutedEventArgs e) => Vm.OpenBio();
    private void OnCancel(object? s, RoutedEventArgs e) => Vm.SetModule(8);
    private void OnAdd(object? s, RoutedEventArgs e) => Vm.AddBlendRow();
    private void OnRemove(object? s, RoutedEventArgs e) { if ((s as Control)?.Tag is BlendRow r) Vm.RemoveBlendRow(r); }
    private async void OnBuild(object? s, RoutedEventArgs e) => await Vm.BuildBlend();
}
