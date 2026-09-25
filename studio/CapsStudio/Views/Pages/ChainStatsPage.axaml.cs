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
    public ChainStatsPage()
    {
        AvaloniaXamlLoader.Load(this);
        DataContextChanged += (_, _) =>
        {
            if (DataContext is not MainViewModel vm || vm == _hooked) return;
            _hooked = vm;
            vm.CsChanged += () =>
            {
                Show(vm);
                var cn = this.FindControl<LinePlot>("Cn")!; cn.RefY = null; cn.SetData(vm.CsCnCurve);
                var ree = this.FindControl<LinePlot>("Ree")!; ree.RefY = null; ree.SetData(vm.CsReeHist);
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
