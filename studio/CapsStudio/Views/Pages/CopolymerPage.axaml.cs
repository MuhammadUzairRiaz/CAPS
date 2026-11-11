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
        this.FindControl<SequenceStrip>("Strip")!.Palette = ["#F0A83C", "#2271DB", "#7DC884"];
        // the view's empty state: said in words, so a black viewport is never left unexplained
        _empty = new TextBlock
        {
            Text = "No chain to show yet. Choose two different monomers; a chain is drawn from the model as you change it, and Redraw draws another.",
            TextWrapping = Avalonia.Media.TextWrapping.Wrap, TextAlignment = Avalonia.Media.TextAlignment.Center, MaxWidth = 380,
            HorizontalAlignment = Avalonia.Layout.HorizontalAlignment.Center, VerticalAlignment = Avalonia.Layout.VerticalAlignment.Center,
            IsHitTestVisible = false,
        };
        _empty.Classes.Add("dim");
        if (this.FindControl<MolView>("View")!.Parent is Panel host) host.Children.Add(_empty);
        DataContextChanged += (_, _) =>
        {
            if (DataContext is not MainViewModel vm || vm == _hooked) return;
            _hooked = vm;
            vm.CoChanged += () => Update(vm);
            vm.PropertyChanged += (_, e) =>
            {
                if (e.PropertyName == nameof(MainViewModel.CoDoc)) Show(vm);
                else if (e.PropertyName == nameof(MainViewModel.CoBuilding)) ShowEmpty(vm);
            };
            ShowEmpty(vm);
        };
    }
    private readonly TextBlock _empty;
    private void ShowEmpty(MainViewModel vm) => _empty.IsVisible = vm.CoDoc == null && !vm.CoBuilding;
    private void Update(MainViewModel vm)
    {
        var strip = this.FindControl<SequenceStrip>("Strip")!;
        strip.SetUnits(vm.CoSequence);
        var f1 = (double)vm.CoF1;
        var points = new System.Collections.Generic.List<SciPoint> { new(f1, vm.CoF1Model, "this feed", 1) };
        if (vm.CoAzeotrope is double z) points.Insert(0, new(z, z, "azeotrope", 2));
        this.FindControl<SciPlot>("Comp")!.Figure = ViewModels.PageFigures.Copolymer(vm.CoCurve, points.ToArray());
    }
    private void Show(MainViewModel vm)
    {
        ShowEmpty(vm);
        var v = this.FindControl<MolView>("View")!;
        v.DrawStyle = 4;    // backbone
        v.ColourMode = 3;   // by the per-atom values: amber for M₁, blue for M₂ (set on the chain's atoms by the view model)
        // The chain grows in a dilute box about 90 Å wide. With the cell hidden the renderer frames the shown atoms
        // alone, so a plain reset centres and fits the chain. A camera focus would add a pan and zoom worked out
        // against the cell centre on top of that, and push the chain out of the view.
        v.ShowCell = false;
        if (v.Document == vm.CoDoc) { v.Refresh(); return; }
        v.Document = vm.CoDoc;
        v.Reset();
        // The fit fills 90 % of the viewport with the chain, which puts its top under the title chip (top left,
        // about 40 px deep). Stepping back a little leaves the chip a clear band above the chain.
        var cam = v.Camera;
        cam.Zoom = 0.82;
        v.Camera = cam;
    }
    protected override void OnPropertyChanged(AvaloniaPropertyChangedEventArgs change)
    {
        base.OnPropertyChanged(change);
        if (change.Property == IsVisibleProperty && IsVisible && DataContext is MainViewModel vm) { Update(vm); Show(vm); }
    }
    private void OnRedraw(object? s, RoutedEventArgs e) => Vm.CoRedraw();
    private async void OnBuild(object? s, RoutedEventArgs e) => await Vm.BuildCopolymerChains();
}
