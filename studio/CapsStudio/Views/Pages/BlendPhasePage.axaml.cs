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
        this.FindControl<SciPlot>("Diagram")!.Figure = ViewModels.PageFigures.BlendPhase(vm.BpBinodal, vm.BpSpinodal, vm.BpCritical, vm.BpTNow);
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
    private async void OnFit(object? s, Avalonia.Interactivity.RoutedEventArgs e) { if (DataContext is ViewModels.MainViewModel vm) await vm.BpFitContacts(); }
}
