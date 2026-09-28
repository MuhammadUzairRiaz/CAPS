using Avalonia.Controls;
using Avalonia.Interactivity;
using Avalonia.Markup.Xaml;
using Avalonia.Platform.Storage;

namespace CapsStudio.Views.Pages;

public partial class EquilibratePage : PageBase
{
    public EquilibratePage()
    {
        AvaloniaXamlLoader.Load(this);
        var rho = this.FindControl<LinePlot>("DensityPlot")!;
        var e = this.FindControl<LinePlot>("EnergyPlot")!;
        var rg = this.FindControl<LinePlot>("RgPlot")!;
        var chain = this.FindControl<LinePlot>("ChainPlot")!;
        rho.RefY = null; e.RefY = null; rg.RefY = null;
        DataContextChanged += (_, _) =>
        {
            if (DataContext is not ViewModels.MainViewModel vm) return;
            vm.ThermoChanged += () =>
            {
                if (!vm.IsEquilibrate && !vm.EqRunning) return;
                rho.Stages = vm.EqStageBands; e.Stages = vm.EqStageBands;
                rho.SetData(vm.Thermo.Select(r => (r.TimePs, r.Density)).ToArray());
                e.SetData(vm.Thermo.Select(r => (r.TimePs, r.Potential)).ToArray());
            };
            vm.EqChecksChanged += () =>
            {
                rg.SetData(vm.EqRgBlocks);
                rho.SetThird(vm.EqDensityBlocks);   // block means over the time series
                e.SetThird(vm.EqEnergyBlocks);
            };
            void Draw() { if (vm.ChainReference.Length > 0) chain.SetCompare(vm.ChainCurve, vm.ChainReference); else chain.SetData(vm.ChainCurve); }
            vm.PropertyChanged += (_, a) => { if (a.PropertyName is nameof(vm.ChainCurve) or nameof(vm.ChainReference)) Draw(); };
            Draw();
        };
    }

    private async void OnRun(object? s, RoutedEventArgs e) => await Vm.RunEquilibrate();
    private void OnQueue(object? s, RoutedEventArgs e) => Vm.QueueEquilibrate();
    private void OnPause(object? s, RoutedEventArgs e) => Vm.TogglePause();
    private void OnCancel(object? s, RoutedEventArgs e) => Vm.CancelEquilibrate();
    private async void OnCbmc(object? s, RoutedEventArgs e) => await Vm.RunCbmc();
    private async void OnExtend(object? s, RoutedEventArgs e) => await Vm.ExtendEquilibrate();
    private void OnAccept(object? s, RoutedEventArgs e) => Vm.AcceptEquilibration();
    private async void OnLoadTarget(object? s, RoutedEventArgs e)
    {
        if (Window == null) return;
        var files = await Window.StorageProvider.OpenFilePickerAsync(new FilePickerOpenOptions
        {
            Title = "Target curve: n and ⟨R²(n)⟩/(n⟨b²⟩)",
            AllowMultiple = false,
            FileTypeFilter = [new FilePickerFileType("Two columns") { Patterns = ["*.csv", "*.dat", "*.txt", "*.tsv", "*.xvg"] }, FilePickerFileTypes.All],
        });
        if (files.Count > 0 && files[0].TryGetLocalPath() is { } path) Vm.LoadEqTarget(path);
    }
    private async void OnSaveTrajectory(object? s, RoutedEventArgs e) { if (Window != null) await Window.SaveTrajectoryAsync(); }
}
