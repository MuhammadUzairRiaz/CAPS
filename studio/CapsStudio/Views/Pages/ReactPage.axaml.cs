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
        var conv = this.FindControl<SciPlot>("ConvPlot")!;
        var gel = this.FindControl<SciPlot>("GelPlot")!;
        DataContextChanged += (_, _) =>
        {
            if (DataContext is not ViewModels.MainViewModel vm) return;
            vm.ReactChanged += () =>
            {
                var rows = vm.ReactRows;
                conv.Figure = ViewModels.PageFigures.ReactLinks(rows);
                gel.Figure = ViewModels.PageFigures.ReactGel(rows);
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
        // each role takes a file out of the list, so one file never fills two roles; a key that stands as its own
        // word in the name (rxn1_pre.lmp, pre.data_template) wins over one buried in a longer word (prepolymer.data)
        static string Name(string p) => System.IO.Path.GetFileName(p);
        static bool Word(string name, string key) => System.Text.RegularExpressions.Regex.IsMatch(name,
            "(^|[^a-z])" + System.Text.RegularExpressions.Regex.Escape(key) + "([^a-z]|$)", System.Text.RegularExpressions.RegexOptions.IgnoreCase);
        string? Take(System.Func<string, bool> test)
        {
            var hit = paths.FirstOrDefault(p => test(Name(p)));
            if (hit != null) paths.Remove(hit);
            return hit;
        }
        string? TakeKey(string key) => Take(n => Word(n, key)) ?? Take(n => n.Contains(key, System.StringComparison.OrdinalIgnoreCase));
        var map = TakeKey("map");
        var pre = TakeKey("pre");
        var post = TakeKey("post");
        // the data file: matched on its extension (or a name starting with data.), so .data_template does not count
        var data = Take(n =>
        {
            var ext = System.IO.Path.GetExtension(n);
            return ext.Equals(".data", System.StringComparison.OrdinalIgnoreCase) || ext.Equals(".lmp", System.StringComparison.OrdinalIgnoreCase)
                || n.StartsWith("data.", System.StringComparison.OrdinalIgnoreCase);
        });
        if (map == null || pre == null || post == null)
        {
            Vm.Status = "Choose the three files of a reaction: names with pre, post and map (and the data file)";
            return;
        }
        Vm.ImportBondReact(pre, post, map, data ?? "");
    }
    private void OnTemplateEditor(object? s, Avalonia.Interactivity.RoutedEventArgs e) => Vm.OpenTemplateEditor();

    private void OnCopyCommand(object? s, Avalonia.Interactivity.RoutedEventArgs e) => CopyCommand(Vm.ReactCommand, "crosslinking");
    private async void OnCopyPython(object? s, RoutedEventArgs e)
    {
        var clip = TopLevel.GetTopLevel(this)?.Clipboard;
        if (clip != null) await clip.SetTextAsync(Vm.ReactPython());
        Vm.Status = "Copied the crosslinking as Python (import caps; doc = caps.open(…))";
    }
    private void OnQueue(object? s, RoutedEventArgs e) => Vm.QueueReact();
    private async void OnInsertCrosslinker(object? sender, Avalonia.Interactivity.RoutedEventArgs e) => await Vm.InsertRxCrosslinker();
}
