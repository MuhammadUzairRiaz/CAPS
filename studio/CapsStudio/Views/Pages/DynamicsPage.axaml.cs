using Avalonia.Controls;
using Avalonia.Interactivity;
using Avalonia.Markup.Xaml;
using Avalonia.Platform.Storage;

namespace CapsStudio.Views.Pages;

public partial class DynamicsPage : PageBase
{
    public DynamicsPage()
    {
        AvaloniaXamlLoader.Load(this);
        var t = this.FindControl<SciPlot>("TPlot")!;
        var rho = this.FindControl<SciPlot>("RhoPlot")!;
        DataContextChanged += (_, _) =>
        {
            if (DataContext is not ViewModels.MainViewModel vm) return;
            vm.ThermoChanged += () =>
            {
                t.Figure = ViewModels.PageFigures.ThermoTemperature(vm.Thermo, (double)(vm.MdTempD ?? 300m));
                rho.Figure = ViewModels.PageFigures.ThermoDensity(vm.Thermo);
            };
            vm.PropertyChanged += (_, e) =>
            {
                if (vm.IsDynamics && e.PropertyName is nameof(vm.IsDynamics) or nameof(vm.Title) or "MdEnsembleText" || (vm.IsDynamics && e.PropertyName != null &&
                    (e.PropertyName.StartsWith("Md", StringComparison.Ordinal) && e.PropertyName is not ("MdDeck" or "MdPreflightSummary" or "MdLog" or "MdRunning"))))
                    vm.RefreshPreflight();
            };
        };
    }

    private async void OnRun(object? s, RoutedEventArgs e) { if (Vm.RunsRemote) await Vm.SubmitRemote("Dynamics"); else await Vm.RunMd(); }
    private void OnQueue(object? s, RoutedEventArgs e) => Vm.QueueMd();
    private void OnPause(object? s, RoutedEventArgs e) => Vm.TogglePause();
    private void OnCancel(object? s, RoutedEventArgs e) => Vm.CancelMd();
    private async void OnContinue(object? s, RoutedEventArgs e) => await Vm.ContinueMd();
    private async void OnSaveTrajectory(object? s, RoutedEventArgs e) { if (Window != null) await Window.SaveTrajectoryAsync(); }
    /// <summary>The deck is written from the structure: not while a job is still changing it.</summary>
    private bool DeckBlocked()
    {
        if (!Vm.Busy) return false;
        Vm.Status = "Wait for the running job to finish: the input is written from the structure, which the job is still changing";
        return true;
    }

    private async void OnCopyDeck(object? s, RoutedEventArgs e)
    {
        if (Vm.Document == null || DeckBlocked()) return;
        var gromacs = Vm.MdGromacs;
        string deck;
        try { deck = await Vm.MdDeckNow(); }
        catch (Exception ex) { Vm.Status = (gromacs ? "Could not write the GROMACS run parameters: " : "Could not write the LAMMPS input: ") + ex.Message; return; }
        var clip = TopLevel.GetTopLevel(this)?.Clipboard;
        if (clip != null) await clip.SetTextAsync(deck);
        Vm.Status = gromacs ? "Copied the GROMACS run parameters" : "Copied the LAMMPS input";
    }

    /// <summary>Writes the deck into a folder: system.data and system.in (LAMMPS), or system.top, system.gro and system.mdp (GROMACS).</summary>
    private async void OnSaveDeck(object? s, RoutedEventArgs e)
    {
        var top = TopLevel.GetTopLevel(this);
        if (top == null || DeckBlocked()) return;
        var dirs = await top.StorageProvider.OpenFolderPickerAsync(new FolderPickerOpenOptions { Title = "Folder for " + Vm.MdDeckFiles.Split(" (")[0] });
        var dir = dirs.Count > 0 ? dirs[0].TryGetLocalPath() : null;
        if (dir == null || Vm.Document == null || DeckBlocked()) return;
        try
        {
            if (Vm.MdGromacs || Vm.MdBoth) Vm.SaveGromacs(dir);
            if (!Vm.MdGromacs)
            {
                Vm.Document!.Save(System.IO.Path.Combine(dir, "system.data"));   // before the await: no job can start in between
                var deck = await Vm.MdDeckNow();   // for the engine and settings as they are now, not the preview
                await System.IO.File.WriteAllTextAsync(System.IO.Path.Combine(dir, "system.in"), deck);
            }
            Vm.Status = $"Wrote {Vm.MdDeckFiles} to {dir}";
        }
        catch (Exception ex) { Vm.Status = "Could not write the deck: " + ex.Message; }
    }

    private void OnCopyCommand(object? s, Avalonia.Interactivity.RoutedEventArgs e) => CopyCommand(Vm.MdCommand, "dynamics run");
    private async void OnCopyPython(object? s, RoutedEventArgs e)
    {
        var clip = TopLevel.GetTopLevel(this)?.Clipboard;
        if (clip != null) await clip.SetTextAsync(Vm.MdPython());
        Vm.Status = "Copied the dynamics run as Python (import caps; doc = caps.open(…))";
    }
    private async void OnCompareEnergies(object? s, RoutedEventArgs e) => await Vm.CompareEnergies();
    private async void OnNveProbe(object? s, RoutedEventArgs e) => await Vm.RunNveProbe();
}
