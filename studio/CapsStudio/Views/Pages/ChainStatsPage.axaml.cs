using Avalonia;
using Avalonia.Controls;
using Avalonia.Interactivity;
using Avalonia.Markup.Xaml;
using Avalonia.Platform.Storage;
using CapsStudio.Interop;
using CapsStudio.ViewModels;

namespace CapsStudio.Views.Pages;

public partial class ChainStatsPage : PageBase
{
    private MainViewModel? _hooked;
    private const string CnHint = "Press Run to compute the characteristic ratio";
    private const string ReeHint = "Press Run to see the end-to-end distance distribution";
    public ChainStatsPage()
    {
        AvaloniaXamlLoader.Load(this);
        // Until a run has filled them, the two charts say what to do instead of showing a blank card.
        var cnPlot = this.FindControl<SciPlot>("Cn")!;
        var reePlot = this.FindControl<SciPlot>("Ree")!;
        cnPlot.EmptyText = CnHint;
        reePlot.EmptyText = ReeHint;
        DataContextChanged += (_, _) =>
        {
            if (DataContext is not MainViewModel vm || vm == _hooked) return;
            _hooked = vm;
            vm.CsChanged += () =>
            {
                Show(vm);
                // The shared figure builders carry an empty hint of "", which would draw nothing at all:
                // an empty curve gets the page's own hint instead.
                // After a run that gave no curve (no chain long enough), say so rather than asking for another run.
                var ran = vm.CsStatus.Length > 0 && !vm.Analyze.Working;
                cnPlot.EmptyText = ran ? "No sub-chain curve from this run" : CnHint;
                reePlot.EmptyText = ran ? "No end-to-end distances from this run" : ReeHint;
                cnPlot.Figure = vm.CsCnCurve.Length > 1 ? ViewModels.PageFigures.ChainCnStats(vm.CsCnCurve) : null;
                reePlot.Figure = vm.CsReeHist.Length > 1 ? ViewModels.PageFigures.EndToEnd(vm.CsReeHist) : null;
            };
        };
    }
    private void Show(MainViewModel vm)
    {
        var v = this.FindControl<MolView>("View")!;
        v.DrawStyle = 4; v.ColourMode = 1; v.ShowCell = true;   // backbone, by molecule
        if (v.Document != vm.Document) { v.Document = vm.Document; v.Reset(); } else v.Refresh();
    }
    protected override void OnPropertyChanged(AvaloniaPropertyChangedEventArgs change)
    {
        base.OnPropertyChanged(change);
        if (change.Property == IsVisibleProperty && IsVisible && DataContext is MainViewModel vm) Show(vm);
    }
    private async void OnRun(object? s, RoutedEventArgs e) => await Vm.RunChainStats();
    private void OnCancel(object? s, RoutedEventArgs e) => Vm.Analyze.Cancel();
    private void OnExport(object? s, RoutedEventArgs e) => Window?.ExportAnalysis();
}
