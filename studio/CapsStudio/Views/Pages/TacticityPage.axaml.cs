using System.Linq;
using Avalonia;
using Avalonia.Controls;
using Avalonia.Interactivity;
using Avalonia.Markup.Xaml;
using CapsStudio.ViewModels;

namespace CapsStudio.Views.Pages;

public partial class TacticityPage : PageBase
{
    private MainViewModel? _hooked;
    public TacticityPage()
    {
        AvaloniaXamlLoader.Load(this);
        this.FindControl<SequenceStrip>("Dyads")!.Palette = ["#F0A83C", "#2271DB"];
        DataContextChanged += (_, _) =>
        {
            if (DataContext is not MainViewModel vm || vm == _hooked) return;
            _hooked = vm;
            vm.TsChanged += () => Update(vm);
        };
    }
    private void Update(MainViewModel vm)
    {
        this.FindControl<SequenceStrip>("Dyads")!.SetUnits(vm.TsDyads);
        var chart = this.FindControl<XyChart>("Pentads")!;
        chart.Categories = vm.TsPentadNames;
        chart.Set(
        [
            new ChartSeries(vm.TsMarkov ? "Markov model" : "Bernoulli model", vm.TsModelPentads.Select((y, i) => ((double)i, y)).ToArray(), "bars", "AccB"),
            new ChartSeries("this chain", vm.TsChainPentads.Select((y, i) => ((double)i, y)).ToArray(), "bars", "#2271DB"),
        ]);
        var v = this.FindControl<MolView>("View")!;
        v.DrawStyle = 0;
        if (v.Document != vm.TsDoc) { v.Document = vm.TsDoc; v.Reset(); }
        if (vm.TsDoc != null && vm.TsFirstUnits.Length > 0)
            try { v.Camera = vm.TsDoc.Focus(v.Camera, vm.TsFirstUnits, 0.85); } catch { v.Refresh(); }
        else v.Refresh();
    }
    protected override void OnPropertyChanged(AvaloniaPropertyChangedEventArgs change)
    {
        base.OnPropertyChanged(change);
        if (change.Property == IsVisibleProperty && IsVisible && DataContext is MainViewModel vm) Update(vm);
    }
    private void OnRedraw(object? s, RoutedEventArgs e) => Vm.TsRedraw();
    private void OnApply(object? s, RoutedEventArgs e) => Vm.TsApply();
    private void OnFit(object? s, RoutedEventArgs e) => Vm.TsFitMeasured();
}
