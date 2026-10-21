using Avalonia;
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
                ShowView(vm);
            };
            // another frame, or another structure made active: the view follows while the page is shown
            vm.PropertyChanged += (_, e) =>
            {
                if (IsVisible && e.PropertyName is nameof(MainViewModel.Frame) or nameof(MainViewModel.DpdHasBeads)) ShowView(vm);
            };
        };
    }

    /// <summary>The run's beads in the page's own view (only while they are the open structure).</summary>
    private void ShowView(MainViewModel vm)
    {
        var view = this.FindControl<MolView>("View")!;
        var doc = vm.DpdHasBeads ? vm.Document : null;
        if (view.Document != doc) { view.Document = doc; if (doc != null) view.Reset(); }
        else if (doc != null) view.Refresh();
    }

    protected override void OnPropertyChanged(AvaloniaPropertyChangedEventArgs change)
    {
        base.OnPropertyChanged(change);
        if (change.Property == IsVisibleProperty && IsVisible && DataContext is MainViewModel vm) ShowView(vm);
    }

    private async void OnRun(object? s, RoutedEventArgs e) => await Vm.RunDpd();
    private void OnCancel(object? s, RoutedEventArgs e) => Vm.CancelDpd();
    private void OnAddSpecies(object? s, RoutedEventArgs e) => Vm.AddDpdSpecies();
    private void OnAddChi(object? s, RoutedEventArgs e) => Vm.AddDpdChi();
    private void OnRemoveSpecies(object? s, RoutedEventArgs e) { if ((s as Control)?.DataContext is DpdSpeciesRow r) Vm.RemoveDpdSpecies(r); }
    private void OnRemoveChi(object? s, RoutedEventArgs e) { if ((s as Control)?.DataContext is DpdChiRow r) Vm.RemoveDpdChi(r); }
}
