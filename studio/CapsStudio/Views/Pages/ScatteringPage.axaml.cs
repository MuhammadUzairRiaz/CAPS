using Avalonia.Controls;
using Avalonia.Interactivity;
using Avalonia.Markup.Xaml;
using Avalonia.Platform.Storage;
using CapsStudio.ViewModels;

namespace CapsStudio.Views.Pages;

public partial class ScatteringPage : PageBase
{
    private MainViewModel? _hooked;

    public ScatteringPage()
    {
        AvaloniaXamlLoader.Load(this);
        DataContextChanged += (_, _) =>
        {
            if (DataContext is not MainViewModel vm || vm == _hooked) return;
            _hooked = vm;
            vm.ScatterChanged += () =>
            {
                this.FindControl<SciPlot>("Plot")!.Figure = ViewModels.PageFigures.Scattering(vm.ScatterXray ? vm.ScatterXrayCurve : [], vm.ScatterNeutron ? vm.ScatterNeutronCurve : [], vm.ScatterExperiment);
            };
        };
    }

    private async void OnRun(object? s, RoutedEventArgs e) => await Vm.RunScattering();
    private void OnCancel(object? s, RoutedEventArgs e) => Vm.Analyze.Cancel();
    private void OnClearExperiment(object? s, RoutedEventArgs e) => Vm.ClearExperiment();

    private async void OnLoadExperiment(object? s, RoutedEventArgs e)
    {
        var top = TopLevel.GetTopLevel(this);
        if (top == null) return;
        var files = await top.StorageProvider.OpenFilePickerAsync(new FilePickerOpenOptions { Title = "Measured pattern (two columns)", AllowMultiple = false });
        if (files.Count > 0 && files[0].TryGetLocalPath() is { } p)
        {
            var why = Vm.LoadExperiment(p);
            if (why != null) Vm.Status = why;
        }
    }
}
