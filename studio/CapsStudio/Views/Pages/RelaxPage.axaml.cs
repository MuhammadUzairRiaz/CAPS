using Avalonia.Controls;
using Avalonia.Interactivity;
using Avalonia.Markup.Xaml;

namespace CapsStudio.Views.Pages;

public partial class RelaxPage : PageBase
{
    public RelaxPage()
    {
        AvaloniaXamlLoader.Load(this);
        var energy = this.FindControl<SciPlot>("EnergyPlot")!;
        var force = this.FindControl<SciPlot>("ForcePlot")!;
        DataContextChanged += (_, _) =>
        {
            if (DataContext is not ViewModels.MainViewModel vm) return;
            vm.RelaxCurvesChanged += () =>
            {
                energy.Figure = ViewModels.PageFigures.RelaxEnergy(vm.RelaxEnergyCurve);
                force.Figure = ViewModels.PageFigures.RelaxForce(vm.RelaxForceCurve, (double)(vm.RelaxFtolD ?? 0.5m));
            };
        };
    }

    private async void OnRelax(object? s, RoutedEventArgs e) { if (Vm.RunsRemote) await Vm.SubmitRemote("Relax"); else await Vm.Relax(); }
    private void OnCancel(object? s, RoutedEventArgs e) => Vm.CancelRelax();
    private void OnAddRestraint(object? s, RoutedEventArgs e) => Vm.AddMeasuredRestraint();
    private void OnRemoveRestraint(object? s, RoutedEventArgs e) { if (s is Control { Tag: CapsStudio.ViewModels.RestraintRow r }) Vm.RemoveRestraint(r); }
    private async void OnSave(object? s, RoutedEventArgs e) { if (Window != null) await Window.SaveAsAsync("data", "LAMMPS data"); }

    private void OnCopyCommand(object? s, Avalonia.Interactivity.RoutedEventArgs e) => CopyCommand(Vm.RelaxCommand, "minimisation");
    private async void OnCopyPython(object? s, RoutedEventArgs e)
    {
        var clip = TopLevel.GetTopLevel(this)?.Clipboard;
        if (clip != null) await clip.SetTextAsync(Vm.RelaxPython());
        Vm.Status = "Copied the minimisation as Python (import caps; doc = caps.open(…))";
    }
    private void OnQueue(object? s, RoutedEventArgs e) => Vm.QueueRelax();
    private void OnHoldSelection(object? s, RoutedEventArgs e) => Vm.HoldSelection();
    private void OnRestrainForRuns(object? s, RoutedEventArgs e) => Vm.RestrainSelectionForRuns();
    private void OnApplyRunRestraints(object? s, RoutedEventArgs e) => Vm.ApplyRunRestraints();
    private void OnFreeAtoms(object? s, RoutedEventArgs e) => Vm.FreeFixedAtoms();

    private async void OnChoosePathB(object? s, Avalonia.Interactivity.RoutedEventArgs e)
    {
        var top = Avalonia.Controls.TopLevel.GetTopLevel(this);
        if (top == null) return;
        var files = await top.StorageProvider.OpenFilePickerAsync(new Avalonia.Platform.Storage.FilePickerOpenOptions { Title = "End B of the path", AllowMultiple = false });
        if (files.Count > 0 && Avalonia.Platform.Storage.StorageProviderExtensions.TryGetLocalPath(files[0]) is { } p && DataContext is CapsStudio.ViewModels.MainViewModel vm) vm.ChoosePathB(p);
    }
    private async void OnRunPath(object? s, Avalonia.Interactivity.RoutedEventArgs e) { if (DataContext is CapsStudio.ViewModels.MainViewModel vm) await vm.RunReactionPath(); }
}
