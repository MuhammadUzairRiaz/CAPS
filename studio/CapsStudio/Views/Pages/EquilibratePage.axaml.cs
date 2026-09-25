using Avalonia.Controls;
using Avalonia.Interactivity;
using Avalonia.Markup.Xaml;

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
                rho.SetData(vm.Thermo.Select(r => (r.TimePs, r.Density)).ToArray());
                e.SetData(vm.Thermo.Select(r => (r.TimePs, r.Potential)).ToArray());
            };
            vm.EqChecksChanged += () => rg.SetData(vm.EqRgBlocks);
            void Draw() { if (vm.RisCurve.Length > 0) chain.SetCompare(vm.ChainCurve, vm.RisCurve); else chain.SetData(vm.ChainCurve); }
            vm.PropertyChanged += (_, a) => { if (a.PropertyName == nameof(vm.ChainCurve)) Draw(); };
            Draw();
        };
    }

    private async void OnRun(object? s, RoutedEventArgs e) => await Vm.RunEquilibrate();
    private void OnCancel(object? s, RoutedEventArgs e) => Vm.CancelEquilibrate();
    private async void OnSaveTrajectory(object? s, RoutedEventArgs e) { if (Window != null) await Window.SaveTrajectoryAsync(); }
}
