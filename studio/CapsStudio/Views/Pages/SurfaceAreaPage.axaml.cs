using System.Linq;
using Avalonia;
using Avalonia.Controls;
using Avalonia.Interactivity;
using Avalonia.Markup.Xaml;
using Avalonia.Platform.Storage;
using CapsStudio.Interop;
using CapsStudio.ViewModels;

namespace CapsStudio.Views.Pages;

public partial class SurfaceAreaPage : PageBase
{
    private MainViewModel? _hooked;
    private const string ConvHint = "Press Run to see the SASA converge with the point count";
    private const string SeriesHint = "Press Run to compute the SASA over the frames";
    private readonly Control? _legend;

    public SurfaceAreaPage()
    {
        AvaloniaXamlLoader.Load(this);
        var conv = this.FindControl<SciPlot>("Conv")!;
        var series = this.FindControl<SciPlot>("Series")!;
        conv.EmptyText = ConvHint;
        series.EmptyText = SeriesHint;
        // The buried-to-exposed ramp under the view: the DockPanel child docked at the bottom of the view's card.
        if (this.FindControl<MolView>("View")?.Parent is Control viewPanel && viewPanel.Parent is DockPanel card)
            _legend = card.Children.OfType<Border>().FirstOrDefault(b => DockPanel.GetDock(b) == Dock.Bottom);
        // The empty-state hint for the per-group table is the XAML TextBlock bound to !SaGroups.Count; no second one here.
        DataContextChanged += (_, _) =>
        {
            if (DataContext is not MainViewModel vm || vm == _hooked) return;
            _hooked = vm;
            vm.SaChanged += () => Update(vm);
            // OpenSurfaceArea may have raised SaChanged before this page existed: fill the panels now.
            Update(vm);
        };
    }
    private void Update(MainViewModel vm)
    {
        Show(vm);
        var conv = this.FindControl<SciPlot>("Conv")!;
        var series = this.FindControl<SciPlot>("Series")!;
        conv.Figure = vm.SaCurve.Length > 1 ? ViewModels.PageFigures.SasaConvergence(vm.SaCurve) : null;
        series.Figure = vm.SaSeries.Length > 1 ? ViewModels.PageFigures.SasaSeries(vm.SaSeries) : null;
    }
    private void Show(MainViewModel vm)
    {
        var v = this.FindControl<MolView>("View")!;
        // Exposure colours only once a run has put per-atom values on this document. Before that, colour by property would fall
        // back to each atom's distance from its molecule's centre on the same ramp, which reads as a SASA result.
        var coloured = vm.Document != null && vm.SaColoured == vm.Document;
        v.ColourMode = coloured ? 3 : 0;   // 3: the per-atom values (exposure) on the ramp; 0: by element
        if (_legend != null) _legend.IsVisible = coloured;
        if (v.Document != vm.Document) { v.Document = vm.Document; v.Reset(); } else v.Refresh();
    }
    protected override void OnPropertyChanged(AvaloniaPropertyChangedEventArgs change)
    {
        base.OnPropertyChanged(change);
        if (change.Property == IsVisibleProperty && IsVisible && DataContext is MainViewModel vm) Show(vm);
    }
    private async void OnRun(object? s, RoutedEventArgs e) => await Vm.RunSurfaceArea();
    private async void OnExport(object? s, RoutedEventArgs e)
    {
        var top = TopLevel.GetTopLevel(this);
        if (top == null) return;
        var file = await top.StorageProvider.SaveFilePickerAsync(new FilePickerSaveOptions
        {
            Title = "Per-atom SASA", SuggestedFileName = "sasa.csv", DefaultExtension = "csv",
            FileTypeChoices = [new FilePickerFileType("CSV") { Patterns = ["*.csv"] }],
        });
        if (file?.TryGetLocalPath() is { } p) Vm.ExportSasaCsv(p);
    }
}
