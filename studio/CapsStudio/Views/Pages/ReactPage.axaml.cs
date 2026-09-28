using Avalonia.Controls;
using Avalonia.Interactivity;
using Avalonia.Markup.Xaml;

namespace CapsStudio.Views.Pages;

public partial class ReactPage : PageBase
{
    public ReactPage()
    {
        AvaloniaXamlLoader.Load(this);
        var conv = this.FindControl<LinePlot>("ConvPlot")!;
        var gel = this.FindControl<LinePlot>("GelPlot")!;
        conv.RefY = null;
        gel.RefY = null;
        DataContextChanged += (_, _) =>
        {
            if (DataContext is not ViewModels.MainViewModel vm) return;
            vm.ReactChanged += () =>
            {
                var rows = vm.ReactRows;
                conv.SetData(rows.Select(r => ((double)r.Cycle, r.Conversion)).ToArray());
                gel.SetData(rows.Select(r => (r.Conversion, r.LargestFraction)).ToArray());
                vm.RaiseGel();
            };
        };
    }

    private async void OnRun(object? s, RoutedEventArgs e) { if (Vm.RunsRemote) await Vm.SubmitRemote("React"); else await Vm.RunReact(); }
    private void OnCancel(object? s, RoutedEventArgs e) => Vm.CancelReact();
    private async void OnSave(object? s, RoutedEventArgs e) { if (Window != null) await Window.SaveAsAsync("data", "LAMMPS data"); }
    private async void OnInsertCurative(object? s, RoutedEventArgs e) => await Vm.InsertCurative();
    private void OnTemplateEditor(object? s, Avalonia.Interactivity.RoutedEventArgs e) => Vm.OpenTemplateEditor();

    private async void OnCopyPython(object? s, RoutedEventArgs e)
    {
        var clip = TopLevel.GetTopLevel(this)?.Clipboard;
        if (clip != null) await clip.SetTextAsync(Vm.ReactPython());
        Vm.Status = "Copied the crosslinking as Python (import caps; doc = caps.open(…))";
    }
    private void OnQueue(object? s, RoutedEventArgs e) => Vm.QueueReact();
}
