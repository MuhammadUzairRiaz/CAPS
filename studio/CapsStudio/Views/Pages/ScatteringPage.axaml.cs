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
                this.FindControl<SciPlot>("Plot")!.Figure = vm.PowderView
                    ? vm.RefineHas ? ViewModels.PageFigures.Refinement(vm.RefineFit.Obs, vm.RefineFit.Calc, vm.RefineFit.Diff) : ViewModels.PageFigures.Powder(vm.PowderCurve, vm.PowderExperiment)
                    : ViewModels.PageFigures.Scattering(vm.ScatterXray ? vm.ScatterXrayCurve : [], vm.ScatterNeutron ? vm.ScatterNeutronCurve : [], vm.ScatterExperiment);
            };
        };
    }

    private async void OnRun(object? s, RoutedEventArgs e) => await Vm.RunScattering();
    private async void OnPowder(object? s, RoutedEventArgs e) => await Vm.RunPowder();
    private async void OnRefine(object? s, RoutedEventArgs e) => await Vm.RunRefine();
    private async void OnIndex(object? s, RoutedEventArgs e) => await Vm.RunIndex();
    private async void OnAdp(object? s, RoutedEventArgs e) => await Vm.RunAdp();
    private async void OnAdpColour(object? s, RoutedEventArgs e) => await Vm.ColourByMobility();
    private async void OnAdpCif(object? s, RoutedEventArgs e)
    {
        var top = TopLevel.GetTopLevel(this);
        if (top?.StorageProvider is not { } sp) return;
        var f = await sp.SaveFilePickerAsync(new Avalonia.Platform.Storage.FilePickerSaveOptions { Title = "CIF with displacement parameters", SuggestedFileName = "adp.cif", DefaultExtension = "cif" });
        if (f?.TryGetLocalPath() is { } path) await Vm.RunAdp(path);
    }
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
