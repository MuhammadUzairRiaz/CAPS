using Avalonia;
using Avalonia.Controls;
using Avalonia.Interactivity;
using Avalonia.Markup.Xaml;
using CapsStudio.Interop;
using CapsStudio.ViewModels;

namespace CapsStudio.Views.Pages;

public partial class SemicrystallinePage : PageBase
{
    private MainViewModel? _hooked;

    public SemicrystallinePage()
    {
        AvaloniaXamlLoader.Load(this);
        DataContextChanged += (_, _) =>
        {
            if (DataContext is not MainViewModel vm || vm == _hooked) return;
            _hooked = vm;
            vm.ScChanged += () =>
            {
                ShowView(vm);
                foreach (var (name, data) in new[] { ("Rho", vm.ScRho), ("P2", vm.ScP2), ("Bridges", vm.ScBridgeHist) })
                {
                    var plot = this.FindControl<LinePlot>(name)!;
                    plot.RefY = null;
                    plot.SetData(data);
                }
            };
        };
    }

    private void ShowView(MainViewModel vm)
    {
        var view = this.FindControl<MolView>("View")!;
        view.ShowCell = true;
        if (view.Document != vm.Document)
        {
            view.Document = vm.Document;
            view.Reset();
            view.Camera = new CapsCamera { Yaw = 0.0, Pitch = -1.35, Zoom = 1.0 };   // from the side: the stack normal up
        }
    }

    protected override void OnPropertyChanged(AvaloniaPropertyChangedEventArgs change)
    {
        base.OnPropertyChanged(change);
        if (change.Property == IsVisibleProperty && IsVisible && DataContext is MainViewModel vm) ShowView(vm);
    }

    private async void OnBuild(object? s, RoutedEventArgs e) => await Vm.BuildSemicrystalline();
    private void OnAnalyse(object? s, RoutedEventArgs e) => Vm.AnalyseSemicrystalline();
    private void OnColour(object? s, RoutedEventArgs e) => Vm.ColourByPhase();
}
