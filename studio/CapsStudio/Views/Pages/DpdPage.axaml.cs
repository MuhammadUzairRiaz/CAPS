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
                var o = vm.DpdOrderSeries;
                this.FindControl<SciPlot>("Order")!.Figure = o.Length < 2 ? null : new SciFigure
                {
                    XLabel = "step", YLabel = "*ψ*", YMin = 0, YZero = true, Empty = "",
                    Series = [new SciSeries("", o.Select(p => p.X).ToArray(), o.Select(p => p.Y).ToArray()) { Kind = o.Length <= 40 ? SciKind.LineMarkers : SciKind.Line, Width = 1.3 }],
                };
                var q = vm.DpdSq;
                this.FindControl<SciPlot>("Sq")!.Figure = q.Length < 2 ? null : new SciFigure
                {
                    XLabel = "*q* / *r*_{c}⁻¹", YLabel = "*S*(*q*)", YMin = 0, YZero = true, Empty = "",
                    Series = [new SciSeries("", q.Select(p => p.X).ToArray(), q.Select(p => Math.Max(0, p.Y)).ToArray()) { Kind = SciKind.LineMarkers, Width = 1.2 }],
                };
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
