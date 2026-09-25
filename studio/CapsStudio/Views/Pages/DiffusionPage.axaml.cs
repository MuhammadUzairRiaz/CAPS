using Avalonia.Controls;
using Avalonia.Interactivity;
using Avalonia.Markup.Xaml;
using CapsStudio.ViewModels;

namespace CapsStudio.Views.Pages;

public partial class DiffusionPage : PageBase
{
    private MainViewModel? _hooked;

    public DiffusionPage()
    {
        AvaloniaXamlLoader.Load(this);
        DataContextChanged += (_, _) =>
        {
            if (DataContext is not MainViewModel vm || vm == _hooked) return;
            _hooked = vm;
            vm.DfChanged += () =>
            {
                var msd = this.FindControl<LinePlot>("Msd")!;
                msd.RefY = null;
                msd.Band = vm.DfWindow;
                msd.SetData(vm.DfMsd);
                var slope = this.FindControl<LinePlot>("Slope")!;
                slope.RefY = 1.0;
                slope.Band = vm.DfWindow;
                slope.SetData(vm.DfSlopeCurve);
            };
        };
    }

    private async void OnRun(object? s, RoutedEventArgs e) => await Vm.RunDiffusion();
    private void OnCancel(object? s, RoutedEventArgs e) => Vm.Analyze.Cancel();
    private void OnExport(object? s, RoutedEventArgs e) => Window?.ExportAnalysis();
}
