using Avalonia.Controls;
using Avalonia.Interactivity;
using Avalonia.Markup.Xaml;
using CapsStudio.ViewModels;

namespace CapsStudio.Views.Pages;

public partial class MechanicsPage : PageBase
{
    private MainViewModel? _hooked;

    public MechanicsPage()
    {
        AvaloniaXamlLoader.Load(this);
        DataContextChanged += (_, _) =>
        {
            if (DataContext is not MainViewModel vm || vm == _hooked) return;
            _hooked = vm;
            vm.MechCurveChanged += () =>
            {
                var plot = this.FindControl<LinePlot>("Curve")!;
                plot.RefY = null;
                plot.Band = vm.MechHasCurve ? (0, vm.MechFitTo) : null;
                plot.SetWithFit(vm.MechCurve, vm.MechFit);
            };
        };
    }

    private async void OnRun(object? s, RoutedEventArgs e) => await Vm.RunMechanics();
    private void OnCancel(object? s, RoutedEventArgs e) => Vm.Analyze.Cancel();
    private void OnExport(object? s, RoutedEventArgs e) => Window?.ExportAnalysis();
}
