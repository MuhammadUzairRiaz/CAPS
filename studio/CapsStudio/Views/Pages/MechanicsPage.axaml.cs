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
                this.FindControl<SciPlot>("Curve")!.Figure = vm.MechHasCurve ? ViewModels.PageFigures.Mechanics(vm.MechCurve, vm.MechFit, vm.MechFitTo, vm.MechSmooth) : null;
            };
        };
    }

    private async void OnRun(object? s, RoutedEventArgs e) { if (Vm.RunsRemote) await Vm.SubmitMechanicsRemote(); else await Vm.RunMechanics(); }
    private void OnCancel(object? s, RoutedEventArgs e) => Vm.Analyze.Cancel();
    private void OnExport(object? s, RoutedEventArgs e) => Window?.ExportAnalysis();
}
