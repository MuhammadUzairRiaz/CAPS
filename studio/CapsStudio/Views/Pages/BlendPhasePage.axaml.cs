using System.IO;
using Avalonia;
using Avalonia.Controls;
using Avalonia.Interactivity;
using Avalonia.Markup.Xaml;
using Avalonia.Platform.Storage;
using CapsStudio.ViewModels;

namespace CapsStudio.Views.Pages;

public partial class BlendPhasePage : PageBase
{
    private MainViewModel? _hooked;
    public BlendPhasePage()
    {
        AvaloniaXamlLoader.Load(this);
        DataContextChanged += (_, _) =>
        {
            if (DataContext is not MainViewModel vm || vm == _hooked) return;
            _hooked = vm;
            vm.BpChanged += () => Update(vm);
        };
    }
    private void Update(MainViewModel vm)
    {
        var markers = vm.BpCritical is { } c ? new[] { new ChartMarker(c.X, c.Y, "TextB", "") } : [];
        var chart = this.FindControl<XyChart>("Diagram")!;
        chart.Set(
        [
            new ChartSeries("binodal (coexistence)", vm.BpBinodal, "line", "AccB"),
            new ChartSeries("spinodal", vm.BpSpinodal, "dash", "SelB", 1.6),
            new ChartSeries("critical point", [], "dots", "TextB"),
        ], markers, [new ChartGuide(vm.BpTNow, false, "DimB", $"{vm.BpTNow:0} K")]);
    }
    protected override void OnPropertyChanged(AvaloniaPropertyChangedEventArgs change)
    {
        base.OnPropertyChanged(change);
        if (change.Property == IsVisibleProperty && IsVisible && DataContext is MainViewModel vm) Update(vm);
    }
    private async void OnExport(object? s, RoutedEventArgs e)
    {
        if (Window?.StorageProvider is not { } sp) return;
        var f = await sp.SaveFilePickerAsync(new FilePickerSaveOptions { Title = "Export the phase diagram", SuggestedFileName = "blend_phase.csv", DefaultExtension = "csv" });
        if (f?.TryGetLocalPath() is { } path) { File.WriteAllText(path, Vm.BlendCsv()); Vm.Status = "Saved " + Path.GetFileName(path); }
    }
}
