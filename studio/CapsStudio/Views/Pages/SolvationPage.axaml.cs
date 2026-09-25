using Avalonia.Controls;
using Avalonia.Interactivity;
using Avalonia.Markup.Xaml;
using CapsStudio.ViewModels;

namespace CapsStudio.Views.Pages;

public partial class SolvationPage : PageBase
{
    public SolvationPage()
    {
        AvaloniaXamlLoader.Load(this);
        var preview = this.FindControl<MolView>("Preview")!;
        DataContextChanged += (_, _) =>
        {
            if (DataContext is not MainViewModel vm) return;
            vm.SolvViewChanged += () => { var fresh = preview.Document == null; preview.Document = vm.SolvDoc; if (fresh) preview.Reset(); };
        };
    }

    private void OnPolymer(object? s, RoutedEventArgs e) => Vm.SetModule(13);
    private void OnCrystal(object? s, RoutedEventArgs e) => Vm.OpenCrystal();
    private void OnSurface(object? s, RoutedEventArgs e) => Vm.OpenSurface();
    private void OnNano(object? s, RoutedEventArgs e) => Vm.OpenNano();
    private void OnBio(object? s, RoutedEventArgs e) => Vm.OpenBio();
    private void OnCancel(object? s, RoutedEventArgs e) => Vm.SetModule(8);
    private async void OnSolvate(object? s, RoutedEventArgs e) => await Vm.Solvate();
    private void OnCg(object? s, Avalonia.Interactivity.RoutedEventArgs e) => Vm.OpenCg();
}
