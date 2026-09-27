using Avalonia.Controls;
using Avalonia.Interactivity;
using Avalonia.Markup.Xaml;
using CapsStudio.Interop;
using CapsStudio.ViewModels;

namespace CapsStudio.Views.Pages;

public partial class NanoPage : PageBase
{
    private void OnGraftSilane(object? s, Avalonia.Interactivity.RoutedEventArgs e) => Vm.GraftSilane();
    public NanoPage()
    {
        AvaloniaXamlLoader.Load(this);
        var preview = this.FindControl<MolView>("Preview")!;
        DataContextChanged += (_, _) =>
        {
            if (DataContext is not MainViewModel vm) return;
            vm.NanoViewChanged += () =>
            {
                preview.Document = vm.NanoDoc;
                preview.Reset();
                // a pore from the side: the walls horizontal, z up (screen up is y, so pitch towards −90°)
                if (vm.NanoIsPore && vm.PoreIsSlit) preview.Camera = new CapsCamera { Yaw = 0.35, Pitch = -1.22, Zoom = 1.0 };
            };
        };
    }

    private void OnPolymer(object? s, RoutedEventArgs e) => Vm.SetModule(13);
    private void OnSurface(object? s, RoutedEventArgs e) => Vm.OpenSurface();
    private void OnSolvation(object? s, RoutedEventArgs e) => Vm.OpenSolvation();
    private void OnCrystal(object? s, RoutedEventArgs e) => Vm.OpenCrystal();
    private void OnBio(object? s, RoutedEventArgs e) => Vm.OpenBio();
    private void OnCancel(object? s, RoutedEventArgs e) => Vm.SetModule(8);
    private async void OnBuild(object? s, RoutedEventArgs e) => await Vm.BuildNano();
    private void OnCg(object? s, Avalonia.Interactivity.RoutedEventArgs e) => Vm.OpenCg();
}
