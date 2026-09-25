using Avalonia.Controls;
using Avalonia.Interactivity;
using Avalonia.Markup.Xaml;
using CapsStudio.ViewModels;

namespace CapsStudio.Views.Pages;

public partial class NanoPage : PageBase
{
    public NanoPage()
    {
        AvaloniaXamlLoader.Load(this);
        var preview = this.FindControl<MolView>("Preview")!;
        DataContextChanged += (_, _) =>
        {
            if (DataContext is not MainViewModel vm) return;
            vm.NanoViewChanged += () => { preview.Document = vm.NanoDoc; preview.Reset(); };
        };
    }

    private void OnPolymer(object? s, RoutedEventArgs e) => Vm.SetModule(13);
    private void OnSurface(object? s, RoutedEventArgs e) => Vm.OpenSurface();
    private void OnSolvation(object? s, RoutedEventArgs e) => Vm.SetModule(5);
    private void OnCrystal(object? s, RoutedEventArgs e) => Vm.OpenCrystal();
    private void OnBio(object? s, RoutedEventArgs e) => Vm.OpenBio();
    private void OnCancel(object? s, RoutedEventArgs e) => Vm.SetModule(8);
    private async void OnBuild(object? s, RoutedEventArgs e) => await Vm.BuildNano();
}
