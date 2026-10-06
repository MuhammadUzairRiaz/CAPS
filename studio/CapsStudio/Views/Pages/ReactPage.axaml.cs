using System.Linq;
using Avalonia.Controls;
using Avalonia.Platform.Storage;
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
                conv.SetData(rows.Count == 0 ? [] : new[] { (0.0, 0.0) }.Concat(rows.Select(r => ((double)r.Cycle, (double)r.Crosslinks))).ToArray());
                conv.RefY = rows.Count > 0 && rows[^1].Target > 0 ? rows[^1].Target : null;
                gel.SetData(rows.Select(r => (r.Conversion, r.LargestFraction)).ToArray());
                this.FindControl<LinkMap>("LinkMap")!.SetLinks(vm.RxLinks.ToList(), vm.Reacting ? 0 : vm.RxChainCount, vm.RxUnitsPerChain);
                vm.RaiseGel();
            };
        };
    }

    private async void OnRun(object? s, RoutedEventArgs e) { if (Vm.RunsRemote) await Vm.SubmitRemote("React"); else await Vm.RunReact(); }
    private void OnCancel(object? s, RoutedEventArgs e) => Vm.CancelReact();
    private async void OnSave(object? s, RoutedEventArgs e) { if (Window != null) await Window.SaveAsAsync("data", "LAMMPS data"); }
    private async void OnInsertCurative(object? s, RoutedEventArgs e) => await Vm.InsertCurative();
    private async void OnInsertSulfur(object? s, RoutedEventArgs e) => await Vm.InsertSulfurDonors();
    private void OnAddReaction(object? s, RoutedEventArgs e) => Vm.AddReaction();
    private void OnLibUse(object? s, RoutedEventArgs e) => Vm.UseRxLib(false);
    private void OnLibAdd(object? s, RoutedEventArgs e) => Vm.UseRxLib(true);
    private async void OnExportBondReact(object? s, RoutedEventArgs e)
    {
        var top = TopLevel.GetTopLevel(this);
        if (top == null) return;
        var dirs = await top.StorageProvider.OpenFolderPickerAsync(new Avalonia.Platform.Storage.FolderPickerOpenOptions { Title = "Folder for the fix bond/react files", AllowMultiple = false });
        if (dirs.Count > 0 && dirs[0].TryGetLocalPath() is { } path) await Vm.ExportBondReact(path);
    }
    private async void OnRunBondReact(object? s, RoutedEventArgs e)
    {
        var top = TopLevel.GetTopLevel(this);
        if (top == null) return;
        var dirs = await top.StorageProvider.OpenFolderPickerAsync(new Avalonia.Platform.Storage.FolderPickerOpenOptions { Title = "Folder for the fix bond/react run", AllowMultiple = false });
        if (dirs.Count > 0 && dirs[0].TryGetLocalPath() is { } path) await Vm.RunBondReact(path);
    }
    private void OnStopBondReact(object? s, RoutedEventArgs e) => Vm.StopBondReact();
    private async void OnImportBondReact(object? s, RoutedEventArgs e)
    {
        var top = TopLevel.GetTopLevel(this);
        if (top == null) return;
        // the three files of one reaction (and the data file for the elements, when the molecule files number their types)
        var files = await top.StorageProvider.OpenFilePickerAsync(new Avalonia.Platform.Storage.FilePickerOpenOptions
        {
            Title = "The pre- and post-reaction templates, the map file (and the data file)", AllowMultiple = true,
        });
        var paths = files.Select(f => f.TryGetLocalPath()).Where(p => p != null).Cast<string>().ToList();
        string? Find(params string[] keys) => paths.FirstOrDefault(p => keys.Any(k => System.IO.Path.GetFileName(p).Contains(k, System.StringComparison.OrdinalIgnoreCase)));
        var map = Find("map");
        var pre = Find("pre");
        var post = Find("post");
        var data = Find(".data", ".lmp", "data.");
        if (map == null || pre == null || post == null)
        {
            Vm.Status = "Choose the three files of a reaction: names with pre, post and map (and the data file)";
            return;
        }
        Vm.ImportBondReact(pre, post, map, data ?? "");
    }
    private void OnTemplateEditor(object? s, Avalonia.Interactivity.RoutedEventArgs e) => Vm.OpenTemplateEditor();

    private async void OnCopyPython(object? s, RoutedEventArgs e)
    {
        var clip = TopLevel.GetTopLevel(this)?.Clipboard;
        if (clip != null) await clip.SetTextAsync(Vm.ReactPython());
        Vm.Status = "Copied the crosslinking as Python (import caps; doc = caps.open(…))";
    }
    private void OnQueue(object? s, RoutedEventArgs e) => Vm.QueueReact();
    private async void OnInsertCrosslinker(object? sender, Avalonia.Interactivity.RoutedEventArgs e) => await Vm.InsertRxCrosslinker();
}
