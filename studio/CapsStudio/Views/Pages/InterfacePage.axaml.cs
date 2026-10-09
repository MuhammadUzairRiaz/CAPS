using System.Linq;
using Avalonia;
using Avalonia.Controls;
using Avalonia.Interactivity;
using Avalonia.Markup.Xaml;
using CapsStudio.Interop;
using CapsStudio.ViewModels;

namespace CapsStudio.Views.Pages;

public partial class InterfacePage : PageBase
{
    private MainViewModel? _hooked;

    public InterfacePage()
    {
        AvaloniaXamlLoader.Load(this);
        DataContextChanged += (_, _) =>
        {
            if (DataContext is not MainViewModel vm || vm == _hooked) return;
            _hooked = vm;
            vm.IfChanged += () =>
            {
                ShowView(vm);
                var plot = this.FindControl<LinePlot>("Profile")!;
                plot.RefY = null;
                plot.XLabel = $"{vm.Analyze.AxisName} (Å) · normal to the surface";
                plot.Band = vm.IfGapBand;
                var film = vm.IfSeries.FirstOrDefault(s => s.Label == "film").Data ?? [];
                var surf = vm.IfSeries.Where(s => s.Label != "film").ToArray();
                plot.SetCompare(film, surf.Length > 0 ? surf[0].Data : []);
                plot.SetThird(surf.Length > 1 ? surf[1].Data : []);
                this.FindControl<TextBlock>("Legend")!.Text = vm.IfSeries.Length == 0 ? "Run to compute the profile."
                    : "film (solid)" + (surf.Length > 0 ? $" · {surf[0].Label} (dashed)" : "") + (surf.Length > 1 ? $" · {surf[1].Label} (dotted)" : "") +
                      (vm.IfGapBand != null ? " · shaded: the gap" : "");
            };
        };
    }

    private void ShowView(MainViewModel vm)
    {
        var view = this.FindControl<MolView>("View")!;
        view.ShowCell = true;
        if (view.Document != vm.Document)
        {
            view.Document = vm.Document;
            view.Reset();
            view.Camera = new CapsCamera { Yaw = 0.0, Pitch = -1.35, Zoom = 1.0 };   // from the side: the surface at the bottom, z up
        }
    }

    protected override void OnPropertyChanged(AvaloniaPropertyChangedEventArgs change)
    {
        base.OnPropertyChanged(change);
        if (change.Property == IsVisibleProperty && IsVisible && DataContext is MainViewModel vm) ShowView(vm);
    }

    private async void OnRun(object? s, RoutedEventArgs e) { if (Vm.RunsRemote && Vm.Analyze.PullOn) await Vm.SubmitPullRemote(); else await Vm.RunInterface(); }
    private void OnCancel(object? s, RoutedEventArgs e) => Vm.Analyze.Cancel();
    private void OnUseHeld(object? s, RoutedEventArgs e) => Vm.IfUseHeld();
    private void OnExport(object? s, RoutedEventArgs e) => Window?.ExportAnalysis();
}
