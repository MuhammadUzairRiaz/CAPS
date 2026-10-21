using Avalonia.Controls;
using Avalonia.Interactivity;
using Avalonia.Markup.Xaml;
using CapsStudio.ViewModels;

namespace CapsStudio.Views.Pages;

public partial class SorptionPage : PageBase
{
    private MainViewModel? _hooked;

    public SorptionPage()
    {
        AvaloniaXamlLoader.Load(this);
        DataContextChanged += (_, _) =>
        {
            if (DataContext is not MainViewModel vm || vm == _hooked) return;
            _hooked = vm;
            vm.SorptionChanged += () =>
            {
                var sp = vm.SorbSpeciesIsotherms;
                var iso = vm.SorbIsotherm;
                SciSeries[] series = sp.Count >= 2   // a mixture: the gases' own isotherms
                    ? sp.Select((s, i) => new SciSeries(i < vm.SorbSpeciesRows.Count ? vm.SorbSpeciesRows[i].Smiles : $"gas {i + 1}", s.Select(q => q.X).ToArray(), s.Select(q => q.Y).ToArray())
                        { Kind = SciKind.LineMarkers, Colour = i, Width = 1.2 }).ToArray()
                    : [new SciSeries("GCMC, ± block error", iso.Select(q => q.X).ToArray(), iso.Select(q => q.Y).ToArray())
                        { Kind = SciKind.LineMarkers, Width = 1.2, Err = vm.SorbIsothermErr.Length == iso.Length ? vm.SorbIsothermErr : null }];
                this.FindControl<SciPlot>("Iso")!.Figure = iso.Length == 0 ? null : new SciFigure
                {
                    XLabel = "*p* / kPa", YLabel = "uptake / mol kg⁻¹", Series = series, YMin = 0, YZero = true, XMin = 0, Empty = "",
                };
            };
            vm.SorptionMapChanged += () =>
            {
                var m = this.FindControl<LinePlot>("Map")!;
                var (x, y, z, xl, yl) = vm.SorbMapProjection();
                m.XLabel = xl;
                m.YLabel = yl;
                m.SetHeat(x, y, z);
            };
        };
    }

    private async void OnRun(object? s, RoutedEventArgs e) => await Vm.RunSorption();
    private void OnCancel(object? s, RoutedEventArgs e) => Vm.CancelSorption();
}
