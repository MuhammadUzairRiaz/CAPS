using Avalonia.Controls;
using Avalonia.Interactivity;
using Avalonia.Markup.Xaml;
using CapsStudio.ViewModels;

namespace CapsStudio.Views.Pages;

public partial class DpdPage : PageBase
{
    private MainViewModel? _hooked;

    public DpdPage()
    {
        AvaloniaXamlLoader.Load(this);
        DataContextChanged += (_, _) =>
        {
            if (DataContext is not MainViewModel vm || vm == _hooked) return;
            _hooked = vm;
            vm.DpdChanged += () =>
            {
                var o = this.FindControl<LinePlot>("Order")!;
                o.RefY = null;
                o.SetData(vm.DpdOrderSeries);
                var q = this.FindControl<LinePlot>("Sq")!;
                q.RefY = null;
                q.SetData(vm.DpdSq);
            };
        };
    }

    private async void OnRun(object? s, RoutedEventArgs e) => await Vm.RunDpd();
    private void OnCancel(object? s, RoutedEventArgs e) => Vm.CancelDpd();
    private void OnAddSpecies(object? s, RoutedEventArgs e) => Vm.AddDpdSpecies();
    private void OnAddChi(object? s, RoutedEventArgs e) => Vm.AddDpdChi();
    private void OnRemoveSpecies(object? s, RoutedEventArgs e) { if ((s as Control)?.DataContext is DpdSpeciesRow r) Vm.RemoveDpdSpecies(r); }
    private void OnRemoveChi(object? s, RoutedEventArgs e) { if ((s as Control)?.DataContext is DpdChiRow r) Vm.RemoveDpdChi(r); }
}
