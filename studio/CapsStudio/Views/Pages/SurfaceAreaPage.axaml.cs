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
    public SurfaceAreaPage()
    {
        AvaloniaXamlLoader.Load(this);
        DataContextChanged += (_, _) =>
        {
            if (DataContext is not MainViewModel vm || vm == _hooked) return;
            _hooked = vm;
            vm.SaChanged += () =>
            {
                Show(vm);
                var c = this.FindControl<LinePlot>("Conv")!; c.RefY = null; c.SetData(vm.SaCurve);
                var sp = this.FindControl<LinePlot>("Series")!; sp.RefY = null; sp.SetData(vm.SaSeries);
            };
        };
    }
    private void Show(MainViewModel vm)
    {
        var v = this.FindControl<MolView>("View")!;
        v.ColourMode = 3;   // the per-atom values (exposure) on the ramp
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
