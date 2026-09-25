using Avalonia.Controls;
using Avalonia.Interactivity;
using Avalonia.Markup.Xaml;

namespace CapsStudio.Views.Pages;

public partial class RelaxPage : PageBase
{
    public RelaxPage()
    {
        AvaloniaXamlLoader.Load(this);
        var energy = this.FindControl<LinePlot>("EnergyPlot")!;
        var force = this.FindControl<LinePlot>("ForcePlot")!;
        energy.RefY = null;
        DataContextChanged += (_, _) =>
        {
            if (DataContext is not ViewModels.MainViewModel vm) return;
            vm.RelaxCurvesChanged += () =>
            {
                force.RefY = System.Math.Log10((double)(vm.RelaxFtolD ?? 0.5m));
                energy.SetData(vm.RelaxEnergyCurve);
                force.SetData(vm.RelaxForceCurve);
            };
        };
    }

    private async void OnRelax(object? s, RoutedEventArgs e) => await Vm.Relax();
    private void OnCancel(object? s, RoutedEventArgs e) => Vm.CancelRelax();
    private void OnAddRestraint(object? s, RoutedEventArgs e) => Vm.AddMeasuredRestraint();
    private void OnRemoveRestraint(object? s, RoutedEventArgs e) { if (s is Control { Tag: CapsStudio.ViewModels.RestraintRow r }) Vm.RemoveRestraint(r); }
    private async void OnSave(object? s, RoutedEventArgs e) { if (Window != null) await Window.SaveAsAsync("data", "LAMMPS data"); }
}
