using Avalonia.Controls;
using Avalonia.Interactivity;
using Avalonia.Markup.Xaml;

namespace CapsStudio.Views.Pages;

public partial class PackPage : PageBase
{
    public PackPage()
    {
        AvaloniaXamlLoader.Load(this);
        DataContextChanged += (_, _) =>
        {
            if (DataContext is ViewModels.MainViewModel vm)
                vm.PackCurveChanged += () => this.FindControl<LinePlot>("ConvPlot")!.SetData(vm.PackCurve.ToArray());
        };
        this.FindControl<LinePlot>("ConvPlot")!.RefY = null;
    }
    public Decorator Slot => this.FindControl<Decorator>("ViewSlot")!;

    private async void OnPack(object? s, RoutedEventArgs e) => await Vm.RunPack();
    private void OnCancel(object? s, RoutedEventArgs e) => Vm.CancelPack();
    private void OnNew(object? s, RoutedEventArgs e) => Vm.NewPackInput();
    private void OnExample(object? s, RoutedEventArgs e) => Window?.PackExample();
    private async void OnAddMolecule(object? s, RoutedEventArgs e) { if (Window != null) await Window.PackAddAsync(); }
    private async void OnOpenInput(object? s, RoutedEventArgs e) { if (Window != null) await Window.PackOpenAsync(); }
    private void OnExportEngines(object? s, RoutedEventArgs e) => Vm.PackExport();
    private void OnMinimise(object? s, RoutedEventArgs e) => Vm.SetModule(2);
    private void OnDynamics(object? s, RoutedEventArgs e) => Vm.SetModule(3);
    private async void OnSave(object? s, RoutedEventArgs e) { if (Window != null) await Window.SaveAsAsync("data", "LAMMPS data"); }
    private void OnQueue(object? s, RoutedEventArgs e) => Vm.QueuePack();
    private async void OnCopyPython(object? s, RoutedEventArgs e)
    {
        var clip = TopLevel.GetTopLevel(this)?.Clipboard;
        if (clip != null) await clip.SetTextAsync(Vm.PackPython());
        Vm.Status = "Copied the packing as Python (caps.pack with this input)";
    }
}
