using Avalonia;
using Avalonia.Controls;
using Avalonia.Interactivity;
using Avalonia.Markup.Xaml;
using CapsStudio.ViewModels;

namespace CapsStudio.Views.Pages;

public partial class CopolymerPage : PageBase
{
    private MainViewModel? _hooked;
    public CopolymerPage()
    {
        AvaloniaXamlLoader.Load(this);
        this.FindControl<SequenceStrip>("Strip")!.Palette = ["#F0A83C", "#2271DB"];
        DataContextChanged += (_, _) =>
        {
            if (DataContext is not MainViewModel vm || vm == _hooked) return;
            _hooked = vm;
            vm.CoChanged += () => Update(vm);
            vm.PropertyChanged += (_, e) => { if (e.PropertyName == nameof(MainViewModel.CoDoc)) Show(vm); };
        };
    }
    private void Update(MainViewModel vm)
    {
        var strip = this.FindControl<SequenceStrip>("Strip")!;
        strip.SetUnits(vm.CoSequence);
        var f1 = (double)vm.CoF1;
        var markers = new System.Collections.Generic.List<ChartMarker> { new(f1, vm.CoF1Model, "AccB") };
        if (vm.CoAzeotrope is double z) markers.Insert(0, new(z, z, "OkB", "azeotrope", 3.5));
        this.FindControl<XyChart>("Comp")!.Set(
            [new ChartSeries("f₁ = F₁", [(0, 0), (1, 1)], "dash", "DimB", 1), new ChartSeries("Mayo–Lewis", vm.CoCurve, "line", "AccB")],
            markers.ToArray());
    }
    private void Show(MainViewModel vm)
    {
        var v = this.FindControl<MolView>("View")!;
        v.DrawStyle = 4;   // backbone, coloured by unit (atom values)
        if (v.Document == vm.CoDoc) { v.Refresh(); return; }
        v.Document = vm.CoDoc;
        v.Reset();
        if (vm.CoDoc is { } d)   // frame the chain, not its (dilute) growth box
            try { v.Camera = d.Focus(v.Camera, System.Linq.Enumerable.Range(0, (int)d.Summary().Atoms).ToArray(), 0.8); } catch { }
    }
    protected override void OnPropertyChanged(AvaloniaPropertyChangedEventArgs change)
    {
        base.OnPropertyChanged(change);
        if (change.Property == IsVisibleProperty && IsVisible && DataContext is MainViewModel vm) { Update(vm); Show(vm); }
    }
    private void OnRedraw(object? s, RoutedEventArgs e) => Vm.CoRedraw();
    private async void OnBuild(object? s, RoutedEventArgs e) => await Vm.BuildCopolymerChains();
}
