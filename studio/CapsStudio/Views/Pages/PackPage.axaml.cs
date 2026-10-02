using CapsStudio.ViewModels;
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
    private async void OnAddFragment(object? s, RoutedEventArgs e) { if ((s as Control)?.Tag is FragmentItem f) await Vm.AddPackMolecule(f.Smiles, f.Name); }
    private void OnSavePackMolecule(object? s, RoutedEventArgs e) => Vm.Status = Vm.SavePackMolecule();
    private void OnRemovePackMolecule(object? s, RoutedEventArgs e) { if ((s as Control)?.Tag is FragmentItem f) Vm.Status = Vm.RemovePackMolecule(f); }
    private async void OnAddSmiles(object? s, RoutedEventArgs e) => await Vm.AddPackMolecule(Vm.PackSmiles);
    private void OnRowCount(object? s, RoutedEventArgs e) => ApplyRowCount(s);
    private void OnRowCountKey(object? s, Avalonia.Input.KeyEventArgs e) { if (e.Key == Avalonia.Input.Key.Enter) ApplyRowCount(s); }
    private void ApplyRowCount(object? s)
    {
        if (s is NumericUpDown { Tag: PackItem row, Value: { } v } && (int)v != (int)(row.CountValue ?? 0)) Vm.SetPackRowCount(row.Row, (int)v);
    }
    private void OnRowForceField(object? s, SelectionChangedEventArgs e)
    {
        if (s is ComboBox { Tag: PackItem row } c && c.IsLoaded && c.SelectedIndex >= 0 && c.SelectedIndex != row.FfChoice)
            Vm.SetPackRowForceField(row.Row, c.SelectedIndex);
    }
    private async void OnSmilesKey(object? s, Avalonia.Input.KeyEventArgs e) { if (e.Key == Avalonia.Input.Key.Enter) { e.Handled = true; await Vm.AddPackMolecule(Vm.PackSmiles); } }
    private async void OnOpenInput(object? s, RoutedEventArgs e) { if (Window != null) await Window.PackOpenAsync(); }
    private async void OnSaveInput(object? s, RoutedEventArgs e) { if (Window != null) await Window.PackSaveAsync(); }
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
    private void OnFillToDensity(object? s, Avalonia.Interactivity.RoutedEventArgs e) => Vm.FillToDensity();
}
