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
            vm.FieldsChanged += () =>
            {
                var (x, y, z, ty) = vm.FieldSlice;
                this.FindControl<SciPlot>("FieldMap")!.Figure = x.Length == 0 ? null : ViewModels.PageFigures.Heat("x (grid)", "y (grid)", "θ_" + ty, x, y, z);
                var o = vm.FieldOrder;
                this.FindControl<SciPlot>("FieldOrder")!.Figure = o.Length < 2 ? null : new SciFigure
                {
                    XLabel = "*τ*", YLabel = "⟨(θ − θ̄)²⟩", YZero = true, Empty = "",
                    Series = [new SciSeries("", o.Select(p => p.X).ToArray(), o.Select(p => p.Y).ToArray()) { Kind = SciKind.Line, Width = 1.2 }],
                };
            };
            vm.DpdChanged += () =>
            {
                var o = vm.DpdOrderSeries;
                this.FindControl<SciPlot>("Order")!.Figure = o.Length < 2 ? null : new SciFigure
                {
                    XLabel = "step", YLabel = "*ψ*", YMin = 0, YZero = true, Empty = "",
                    Series = [new SciSeries("", o.Select(p => p.X).ToArray(), o.Select(p => p.Y).ToArray()) { Kind = ViewModels.PageFigures.Auto(o.Length), Width = 1.3 }],
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
    private void OnAddRegion(object? s, RoutedEventArgs e) => Vm.AddDpdRegion();
    private void OnRemoveRegion(object? s, RoutedEventArgs e) { if ((s as Control)?.DataContext is DpdRegionRow r) Vm.RemoveDpdRegion(r); }
    private void OnRemoveSpecies(object? s, RoutedEventArgs e) { if ((s as Control)?.DataContext is DpdSpeciesRow r) Vm.RemoveDpdSpecies(r); }
    private void OnRemoveChi(object? s, RoutedEventArgs e) { if ((s as Control)?.DataContext is DpdChiRow r) Vm.RemoveDpdChi(r); }
    private async void OnFields(object? s, Avalonia.Interactivity.RoutedEventArgs e) => await Vm.RunFields();
}
