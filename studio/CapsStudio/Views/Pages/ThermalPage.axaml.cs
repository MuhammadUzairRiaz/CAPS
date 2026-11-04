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
            vm.PhononsChanged += () =>
            {
                var (x, bands, ticks) = vm.PhDispersion;
                this.FindControl<SciPlot>("Bands")!.Figure = bands.Length == 0 ? null : new SciFigure
                {
                    XLabel = "wave vector along the path / Å⁻¹", YLabel = "*ν̃* / cm⁻¹", Empty = "",
                    Series = bands.Select(b => new SciSeries("", x, b) { Kind = SciKind.Line, Colour = 0, Width = 1.0 }).ToArray(),
                    Marks = ticks.Select(tk => new SciMark(tk.X, tk.Label)).ToArray(), XMin = x.FirstOrDefault(), XMax = x.LastOrDefault(),
                };
                var (w, dos, parts) = vm.PhDos;
                var series = new List<SciSeries> { new("total", w, dos) { Kind = SciKind.Line, Colour = -1, Width = 1.2 } };
                if (parts.Length > 1) series.AddRange(parts.Select((p, i) => new SciSeries(p.Label, w, p.Y) { Kind = SciKind.Line, Colour = i, Width = 1.0 }));
                this.FindControl<SciPlot>("Dos")!.Figure = w.Length == 0 ? null : new SciFigure
                {
                    XLabel = "*ν̃* / cm⁻¹", YLabel = "states / cm⁻¹ per cell", YZero = true, YMin = 0, Series = [.. series], Empty = "",
                };
            };
            vm.DielectricChanged += () =>
            {
                var (w, refl) = vm.DlReflectivity;
                this.FindControl<SciPlot>("Reflect")!.Figure = w.Length == 0 ? null : new SciFigure
                {
                    XLabel = "*ν̃* / cm⁻¹", YLabel = "reflectivity *R*", YMin = 0, Empty = "",
                    Series = refl.Select((r, i) => new SciSeries($"along principal axis {i + 1}", w, r) { Kind = SciKind.Line, Colour = i, Width = 1.2 }).ToArray(),
                };
            };
            vm.ThChanged += () =>
            {
                this.FindControl<SciPlot>("Curve")!.Figure = ViewModels.PageFigures.Thermal(vm.ThCurves, vm.ThGreenKubo);
            };
        };
    }

    private async void OnRun(object? s, RoutedEventArgs e) => await Vm.RunThermal();
    private async void OnDielectric(object? s, RoutedEventArgs e) => await Vm.RunDielectric();
    private async void OnPhonons(object? s, RoutedEventArgs e) => await Vm.RunPhonons();

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
