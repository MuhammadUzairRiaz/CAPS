using Avalonia.Controls;
using Avalonia.Interactivity;
using Avalonia.Markup.Xaml;
using Avalonia.Platform.Storage;
using Avalonia.VisualTree;
using CapsStudio.ViewModels;

namespace CapsStudio.Views.Pages;

public partial class TrajectoryPage : PageBase
{
    public TrajectoryPage()
    {
        AvaloniaXamlLoader.Load(this);
        var preview = this.FindControl<MolView>("Preview")!;
        var markers = this.FindControl<MarkerStrip>("Markers")!;
        DataContextChanged += (_, _) =>
        {
            if (DataContext is not MainViewModel vm) return;
            vm.TrajectoryChanged += () =>
            {
                if (!vm.IsTrajectory) return;
                if (!ReferenceEquals(preview.Document, vm.Document)) { preview.Document = vm.Document; preview.Reset(); }
                preview.LineA = vm.TrajEnds.Length == 2 ? vm.TrajEnds[0] : -1;
                preview.LineB = vm.TrajEnds.Length == 2 ? vm.TrajEnds[1] : -1;
                preview.LineLabel = vm.TrajReeText;
                preview.Refresh();
                markers.SetMarkers(vm.TrajRunFrames.Select(f => vm.FrameMax > 0 ? f / (double)vm.FrameMax : 0).ToList());
                // the plots: their series, colour and the cursor at this frame's time
                var plots = this.GetVisualDescendants().OfType<LinePlot>().ToList();
                for (var k = 0; k < plots.Count && k < vm.TrajPlots.Length; k++)
                {
                    plots[k].LineBrush = vm.TrajPlots[k].Colour;
                    plots[k].RefY = null;
                    plots[k].YLabel = vm.TrajPlots[k].Series ?? "";
                    plots[k].CursorX = vm.TrajCursor;
                    plots[k].SetData(vm.TrajPlots[k].Points);
                }
            };
        };
    }

    private void OnBack(object? s, RoutedEventArgs e) { if (Vm.IsPlaying) Window?.TogglePlayback(); Vm.SetModule(8); }
    private void OnPlay(object? s, RoutedEventArgs e) => Window?.TogglePlayback();
    private void OnPrev(object? s, RoutedEventArgs e) => Vm.StepFrame(-1);
    private void OnNext(object? s, RoutedEventArgs e) => Vm.StepFrame(1);
    private void OnMovie(object? s, RoutedEventArgs e) => Vm.OpenRender();

    private async void OnLog(object? s, RoutedEventArgs e)
    {
        var top = TopLevel.GetTopLevel(this);
        if (top == null) return;
        var files = await top.StorageProvider.OpenFilePickerAsync(new FilePickerOpenOptions { Title = "A LAMMPS log for the thermo plots", AllowMultiple = false });
        if (files.Count > 0 && files[0].TryGetLocalPath() is { } path) Vm.TrajLogPath = path;
    }
}
