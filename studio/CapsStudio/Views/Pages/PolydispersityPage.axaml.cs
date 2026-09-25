using System.Linq;
using Avalonia;
using Avalonia.Controls;
using Avalonia.Interactivity;
using Avalonia.Markup.Xaml;
using CapsStudio.ViewModels;

namespace CapsStudio.Views.Pages;

public partial class PolydispersityPage : PageBase
{
    private MainViewModel? _hooked;
    private int _style = 4;
    public PolydispersityPage()
    {
        AvaloniaXamlLoader.Load(this);
        DataContextChanged += (_, _) =>
        {
            if (DataContext is not MainViewModel vm || vm == _hooked) return;
            _hooked = vm;
            vm.PdChanged += () => { Chart(vm); if (IsVisible) Show(vm); };
        };
    }
    private void Chart(MainViewModel vm)
    {
        var drawn = vm.PdLengths.Select(s => (double.Parse(s, System.Globalization.CultureInfo.InvariantCulture), double.NaN)).ToArray();
        this.FindControl<XyChart>("Dist")!.Set(
        [
            new ChartSeries("number fraction", vm.PdNumber, "line", "SelB"),
            new ChartSeries("weight fraction", vm.PdWeight, "line", "AccB"),
            new ChartSeries($"the {drawn.Length} drawn chains", drawn, "dots", "TextB", 2.6),
        ]);
    }
    /// <summary>The grown cell with each chain coloured by its length (only while this page shows).</summary>
    private void Show(MainViewModel vm)
    {
        var v = this.FindControl<MolView>("View")!;
        var doc = vm.PdGrown != null && vm.PdGrown == vm.Document ? vm.Document : null;
        this.FindControl<Border>("Hint")!.IsVisible = doc == null;
        this.FindControl<Border>("Styles")!.IsVisible = doc != null;
        if (doc == null) { v.Document = null; return; }
        var lengths = vm.PdGrownLengths;
        var (mol, _) = doc.MoleculeIndex((int)doc.Summary().Atoms);
        double lo = lengths.Min(), hi = lengths.Max();
        doc.SetAtomValues(mol.Select(m => m >= 0 && m < lengths.Length ? (hi > lo ? (lengths[m] - lo) / (hi - lo) : 0.5) : 0.5).ToArray(), 1);
        v.DrawStyle = _style; v.ShowCell = true;
        if (v.Document != doc) { v.Document = doc; v.Reset(); } else v.Refresh();
    }
    protected override void OnPropertyChanged(AvaloniaPropertyChangedEventArgs change)
    {
        base.OnPropertyChanged(change);
        if (change.Property != IsVisibleProperty || DataContext is not MainViewModel vm) return;
        if (IsVisible) { Chart(vm); Show(vm); }
        else if (vm.PdGrown != null && vm.PdGrown == vm.Document) vm.Document?.SetAtomValues(null);   // the main view keeps its own colours
    }
    private void OnRedraw(object? s, RoutedEventArgs e) => Vm.PdRedraw();
    private void OnUse(object? s, RoutedEventArgs e) => Vm.UsePdLengths();
    private void OnStyle(object? s, RoutedEventArgs e)
    {
        if (s is Control { Tag: string t } && int.TryParse(t, out var k)) { _style = k; Show(Vm); }
    }
}
