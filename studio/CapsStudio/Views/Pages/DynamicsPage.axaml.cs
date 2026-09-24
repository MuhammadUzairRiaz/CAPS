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
        var t = this.FindControl<LinePlot>("TPlot")!;
        var rho = this.FindControl<LinePlot>("RhoPlot")!;
        rho.RefY = null;
        DataContextChanged += (_, _) =>
        {
            if (DataContext is not ViewModels.MainViewModel vm) return;
            vm.ThermoChanged += () =>
            {
                t.RefY = (double)(vm.MdTempD ?? 300m);
                t.SetData(vm.Thermo.Select(r => (r.TimePs, r.Temperature)).ToArray());
                rho.SetData(vm.Thermo.Select(r => (r.TimePs, r.Density)).ToArray());
            };
            vm.PropertyChanged += (_, e) =>
            {
                if (vm.IsDynamics && e.PropertyName is nameof(vm.IsDynamics) or nameof(vm.Title) or "MdEnsembleText" || (vm.IsDynamics && e.PropertyName != null &&
                    (e.PropertyName.StartsWith("Md", StringComparison.Ordinal) && e.PropertyName is not ("MdDeck" or "MdPreflightSummary" or "MdLog" or "MdRunning"))))
                    vm.RefreshPreflight();
            };
        };
    }

    private async void OnRun(object? s, RoutedEventArgs e) => await Vm.RunMd();
    private void OnCancel(object? s, RoutedEventArgs e) => Vm.CancelMd();
    private async void OnSaveTrajectory(object? s, RoutedEventArgs e) { if (Window != null) await Window.SaveTrajectoryAsync(); }
    private async void OnCopyDeck(object? s, RoutedEventArgs e)
    {
        var clip = TopLevel.GetTopLevel(this)?.Clipboard;
        if (clip != null) await clip.SetTextAsync(Vm.MdDeck);
        Vm.Status = "Copied the LAMMPS input";
    }

    /// <summary>Writes system.data (with the force field) and system.in into a folder.</summary>
    private async void OnSaveDeck(object? s, RoutedEventArgs e)
    {
        var top = TopLevel.GetTopLevel(this);
        if (top == null) return;
        var dirs = await top.StorageProvider.OpenFolderPickerAsync(new FolderPickerOpenOptions { Title = "Folder for system.data and system.in" });
        var dir = dirs.Count > 0 ? dirs[0].TryGetLocalPath() : null;
        if (dir == null) return;
        try
        {
            Vm.Document!.Save(System.IO.Path.Combine(dir, "system.data"));
            await System.IO.File.WriteAllTextAsync(System.IO.Path.Combine(dir, "system.in"), Vm.MdDeck);
            Vm.Status = $"Wrote system.data and system.in to {dir} (lmp -in system.in)";
        }
        catch (Exception ex) { Vm.Status = "Could not write the deck: " + ex.Message; }
    }
}
