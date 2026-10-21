using Avalonia.Controls;
using Avalonia.Interactivity;
using Avalonia.Markup.Xaml;
using CapsStudio.ViewModels;

namespace CapsStudio.Views.Pages;

public partial class OrientationPage : PageBase
{
    private MainViewModel? _hooked;

    public OrientationPage()
    {
        AvaloniaXamlLoader.Load(this);
        DataContextChanged += (_, _) =>
        {
            if (DataContext is not MainViewModel vm || vm == _hooked) return;
            _hooked = vm;
            vm.OrChanged += () =>
            {
                this.FindControl<SciPlot>("PerFrame")!.Figure = ViewModels.PageFigures.OrderPerFrame(vm.OrPerFrame);
                this.FindControl<SciPlot>("AlongZ")!.Figure = ViewModels.PageFigures.OrderAlong(vm.Analyze.AxisName, vm.OrAlongZ);
            };
        };
    }

    private async void OnRun(object? s, RoutedEventArgs e) => await Vm.RunOrientation();
    private void OnCancel(object? s, RoutedEventArgs e) => Vm.Analyze.Cancel();
    private void OnExport(object? s, RoutedEventArgs e) => Window?.ExportAnalysis();
}
