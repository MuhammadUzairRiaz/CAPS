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
                var f = this.FindControl<LinePlot>("PerFrame")!;
                f.RefY = null;
                f.SetData(vm.OrPerFrame);
                var z = this.FindControl<LinePlot>("AlongZ")!;
                z.XLabel = $"{vm.Analyze.AxisName} (Å)";
                z.RefY = 0;
                z.SetData(vm.OrAlongZ);
            };
        };
    }

    private async void OnRun(object? s, RoutedEventArgs e) => await Vm.RunOrientation();
    private void OnCancel(object? s, RoutedEventArgs e) => Vm.Analyze.Cancel();
    private void OnExport(object? s, RoutedEventArgs e) => Window?.ExportAnalysis();
}
