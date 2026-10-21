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
        var rho = this.FindControl<SciPlot>("DensityPlot")!;
        var e = this.FindControl<SciPlot>("EnergyPlot")!;
        var rg = this.FindControl<SciPlot>("RgPlot")!;
        var chain = this.FindControl<SciPlot>("ChainPlot")!;
        DataContextChanged += (_, _) =>
        {
            if (DataContext is not ViewModels.MainViewModel vm) return;
            // the time series with the protocol's stages and the block means over them
            void Series()
            {
                rho.Figure = ViewModels.PageFigures.ThermoDensity(vm.Thermo, vm.EqStageBands, vm.EqDensityBlocks);
                e.Figure = ViewModels.PageFigures.ThermoEnergy(vm.Thermo, vm.EqStageBands, vm.EqEnergyBlocks);
            }
            vm.ThermoChanged += () => { if (vm.IsEquilibrate || vm.EqRunning) Series(); };
            vm.EqChecksChanged += () =>
            {
                rg.Figure = ViewModels.PageFigures.RgBlocks(vm.EqRgBlocks);
                Series();
            };
            void Draw() => chain.Figure = ViewModels.PageFigures.ChainCn(vm.ChainCurve, vm.ChainReference, vm.ChainReferenceName);
            vm.PropertyChanged += (_, a) => { if (a.PropertyName is nameof(vm.ChainCurve) or nameof(vm.ChainReference)) Draw(); };
            Draw();
        };
    }

    private void OnCopyCommand(object? s, Avalonia.Interactivity.RoutedEventArgs e) => CopyCommand(Vm.EqCommand, "equilibration");
    private async void OnRun(object? s, RoutedEventArgs e) { if (Vm.RunsRemote) await Vm.SubmitRemote("Equilibrate"); else await Vm.RunEquilibrate(); }
    private void OnQueue(object? s, RoutedEventArgs e) => Vm.QueueEquilibrate();
    private void OnPause(object? s, RoutedEventArgs e) => Vm.TogglePause();
    private void OnCancel(object? s, RoutedEventArgs e) => Vm.CancelEquilibrate();
    private async void OnCbmc(object? s, RoutedEventArgs e) { if (Vm.RunsRemote) await Vm.SubmitCbmcRemote(); else await Vm.RunCbmc(); }
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
