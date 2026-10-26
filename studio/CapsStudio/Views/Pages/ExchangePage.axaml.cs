using Avalonia;
using Avalonia.Controls;
using Avalonia.Interactivity;
using Avalonia.Markup.Xaml;
using Avalonia.Platform.Storage;
using CapsStudio.ViewModels;

namespace CapsStudio.Views.Pages;

public partial class ExchangePage : PageBase
{
    private MainViewModel? _hooked;

    public ExchangePage()
    {
        AvaloniaXamlLoader.Load(this);
        DataContextChanged += (_, _) =>
        {
            if (DataContext is not MainViewModel vm || vm == _hooked) return;
            _hooked = vm;
            vm.ExChanged += () => Draw(vm);
        };
    }

    private void Draw(MainViewModel vm)
    {
        var live = vm.ExLive.ToArray();
        this.FindControl<SciPlot>("Live")!.Figure = live.Length < 2 ? null : new SciFigure
        {
            XLabel = "*t* / ps", YLabel = "accepted moves", YMin = 0, YZero = true, Empty = "",
            Series = [new SciSeries("", live.Select(p => p.X).ToArray(), live.Select(p => p.Y).ToArray()) { Kind = SciKind.Steps, Width = 1.4 }],
        };
        var c = vm.ExCurve;
        var c2 = vm.ExCurve2;
        var title = this.FindControl<TextBlock>("CurveTitle")!;
        SciFigure? fig = null;
        switch (vm.ExCurveKind)
        {
            case "run":
                title.Text = "Accepted moves against time";
                if (c.Length > 1) fig = new SciFigure { XLabel = "*t* / ps", YLabel = "accepted moves", YMin = 0, YZero = true, Series = [new SciSeries("", X(c), Y(c)) { Kind = SciKind.Steps }] };
                break;
            case "step_strain":
                title.Text = "Relaxation modulus after the step strain";
                var pos = c.Where(p => p.Y > 0).ToArray();
                if (pos.Length > 1)
                    fig = new SciFigure
                    {
                        XLabel = "*t* / ps", YLabel = "*G*(*t*) = *σ*(*t*)/*ε* / MPa", XLog = true, YLog = true,
                        Series = [new SciSeries("", X(pos), Y(pos)) { Kind = SciKind.LineMarkers, Width = 1.2 }],
                        Notes = pos.Length < c.Length ? ["points where the stress is not positive are left out of the log plot"] : [],
                    };
                break;
            case "tv":
                title.Text = "Specific volume and exchange rate against temperature";
                if (c.Length > 1)
                    fig = new SciFigure
                    {
                        XLabel = "*T* / K", YLabel = "*v* / cm³ g⁻¹", Y2Label = "accepted / ps",
                        Series = [new SciSeries("specific volume", X(c), Y(c)) { Kind = SciKind.LineMarkers },
                                  new SciSeries("exchange rate", X(c2), Y(c2)) { Kind = SciKind.Dashed, Right = true, Colour = 1 }],
                    };
                break;
            case "heal_run":
                title.Text = "Bonds across the cut while it heals";
                if (c.Length > 0)
                    fig = new SciFigure
                    {
                        XLabel = "*t* / ps", YLabel = "bonds", YMin = 0, YZero = true,
                        Series = [new SciSeries("bridging the two sides", X(c), Y(c)) { Kind = SciKind.LineMarkers },
                                  new SciSeries("crossing the plane", X(c2), Y(c2)) { Kind = SciKind.Dashed, Colour = 1 }],
                    };
                break;
        }
        this.FindControl<SciPlot>("Curve")!.Figure = fig;
        var h = vm.ExHealCurve;
        var i = vm.ExIntactCurve;
        this.FindControl<SciPlot>("Heal")!.Figure = h.Length < 2 || i.Length < 2 ? null : new SciFigure
        {
            XLabel = "strain", YLabel = "*σ* / MPa", YZero = true,
            Series = [new SciSeries("intact", X(i), Y(i)) { Colour = -1 }, new SciSeries("healed", X(h), Y(h)) { Colour = 0 }],
        };
    }

    // shown by any route (the menu, the palette, a project reopened): the library is read once
    protected override void OnPropertyChanged(AvaloniaPropertyChangedEventArgs change)
    {
        base.OnPropertyChanged(change);
        if (change.Property == IsVisibleProperty && IsVisible && DataContext is MainViewModel vm && vm.ExChemistries.Count == 0) vm.LoadExchangeLibrary();
    }

    private static double[] X((double X, double Y)[] c) => c.Select(p => p.X).ToArray();
    private static double[] Y((double X, double Y)[] c) => c.Select(p => p.Y).ToArray();

    private async Task<string?> Pick(string title)
    {
        if (TopLevel.GetTopLevel(this) is not { } top) return null;
        var f = await top.StorageProvider.OpenFilePickerAsync(new FilePickerOpenOptions { Title = title, AllowMultiple = false });
        return f.Count > 0 ? f[0].TryGetLocalPath() : null;
    }

    private async void OnRun(object? s, RoutedEventArgs e) => await Vm.RunExchange();
    private void OnCancel(object? s, RoutedEventArgs e) => Vm.CancelExchange();
    private void OnCut(object? s, RoutedEventArgs e) => Vm.BuildExchangeStart();
    private void OnCompare(object? s, RoutedEventArgs e) => Vm.CompareHealing();
    private void OnCopyCommand(object? s, RoutedEventArgs e) => CopyCommand(Vm.ExchangeCommand, "exchange run");
    private async void OnHealed(object? s, RoutedEventArgs e) { if (await Pick("The healed sample's stress–strain table") is { } p) Vm.ExHealed = p; }
    private async void OnIntact(object? s, RoutedEventArgs e) { if (await Pick("The intact sample's stress–strain table") is { } p) Vm.ExIntact = p; }
}
