using System.IO;
using Avalonia.Controls;
using Avalonia.Interactivity;
using Avalonia.Markup.Xaml;
using Avalonia.Platform.Storage;
using CapsStudio.ViewModels;

namespace CapsStudio.Views.Pages;

public partial class ThermalPage : PageBase
{
    private MainViewModel? _hooked;

    public ThermalPage()
    {
        AvaloniaXamlLoader.Load(this);
        DataContextChanged += (_, _) =>
        {
            if (DataContext is not MainViewModel vm || vm == _hooked) return;
            _hooked = vm;
            vm.ThChanged += () =>
            {
                var plot = this.FindControl<LinePlot>("Curve")!;
                plot.RefY = null;
                var c = vm.ThCurves;
                plot.SetCompare(c.Length > 0 ? c[0] : [], c.Length > 1 ? c[1] : []);
                plot.SetThird(c.Length > 2 ? c[2] : []);
            };
        };
    }

    private async void OnRun(object? s, RoutedEventArgs e) => await Vm.RunThermal();

    private async void OnLammps(object? s, RoutedEventArgs e)
    {
        var top = TopLevel.GetTopLevel(this);
        if (top == null) return;
        var file = await top.StorageProvider.SaveFilePickerAsync(new FilePickerSaveOptions
        {
            Title = "LAMMPS input for the thermal conductivity (the data file is written beside it)", SuggestedFileName = "in.kappa.lmp", DefaultExtension = "lmp",
            FileTypeChoices = [new FilePickerFileType("LAMMPS input") { Patterns = ["*.lmp", "*.in"] }],
        });
        if (file?.TryGetLocalPath() is not { } p) return;
        Vm.Status = Vm.WriteThermalLammps(p);
    }
}
