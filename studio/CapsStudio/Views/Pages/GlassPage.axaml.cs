using System.Linq;
using Avalonia.Controls;
using Avalonia.Interactivity;
using Avalonia.Markup.Xaml;
using Avalonia.Platform.Storage;
using CapsStudio.ViewModels;

namespace CapsStudio.Views.Pages;

public partial class GlassPage : PageBase
{
    private void OnCopyCommand(object? s, Avalonia.Interactivity.RoutedEventArgs e) => CopyCommand(Vm.TgCommand, "cooling scan");
    private MainViewModel? _hooked;

    public GlassPage()
    {
        AvaloniaXamlLoader.Load(this);
        DataContextChanged += (_, _) =>
        {
            if (DataContext is not MainViewModel vm || vm == _hooked) return;
            _hooked = vm;
            vm.GtChanged += () =>
            {
                this.FindControl<SciPlot>("Schedule")!.Figure = ViewModels.PageFigures.GlassSchedule(vm.GtSchedule);
                this.FindControl<SciPlot>("Vt")!.Figure = vm.GtPoints.Length > 1
                    ? AnalyzeFigures.Cooling(vm.GtPoints, vm.GtErrors, vm.GtLines, vm.GtYLabel.StartsWith("potential", StringComparison.Ordinal)) : null;
            };
        };
    }

    private void OnPlan(object? s, NumericUpDownValueChangedEventArgs e) => Vm.RaiseGlassPlan();
    private async void OnRun(object? s, RoutedEventArgs e) { if (Vm.RunsRemote) await Vm.SubmitGlassRemote(); else await Vm.RunGlass(); }
    private void OnCancel(object? s, RoutedEventArgs e) => Vm.Analyze.Cancel();

    private async void OnRecipe(object? s, RoutedEventArgs e)
    {
        var top = TopLevel.GetTopLevel(this);
        var src = Vm.Document?.Path;
        if (top == null || string.IsNullOrEmpty(src)) { Vm.Status = "Save the structure first: the recipe reads it from its file"; return; }
        var file = await top.StorageProvider.SaveFilePickerAsync(new FilePickerSaveOptions
        {
            Title = "Save the cooling scan as a recipe", SuggestedFileName = System.IO.Path.GetFileNameWithoutExtension(src) + "_tg.yaml", DefaultExtension = "yaml",
            FileTypeChoices = [new FilePickerFileType("CAPS recipe") { Patterns = ["*.yaml", "*.yml"] }],
        });
        if (file?.TryGetLocalPath() is not { } p) return;
        string? ff;
        try { ff = Vm.SaveGlassForceField(p); }
        catch (Exception ex) { Vm.Status = "The force field was not saved: " + ex.Message; return; }
        System.IO.File.WriteAllText(p, Vm.GlassRecipe(src, ff));
        Vm.Status = $"Saved {System.IO.Path.GetFileName(p)}{(ff != null ? " and " + ff : "")} · caps run {System.IO.Path.GetFileName(p)}";
    }
    private async void OnOpenFinishedFolders(object? s, RoutedEventArgs e)
    {
        var top = TopLevel.GetTopLevel(this);
        if (top == null) return;
        var dirs = await top.StorageProvider.OpenFolderPickerAsync(new FolderPickerOpenOptions { Title = "Finished cooling runs: one folder per replica (caps tg -o OUT.data --csv csv)", AllowMultiple = true });
        var paths = dirs.Select(d => d.TryGetLocalPath()).OfType<string>().ToList();
        if (paths.Count > 0) Vm.OpenFinishedGlass(paths);
    }
    private async void OnOpenFinishedFiles(object? s, RoutedEventArgs e)
    {
        var top = TopLevel.GetTopLevel(this);
        if (top == null) return;
        var files = await top.StorageProvider.OpenFilePickerAsync(new FilePickerOpenOptions
        {
            Title = "Finished cooling runs: their OUT.tg-points.json, properties files or tables (one per replica)", AllowMultiple = true,
            FileTypeFilter = [new FilePickerFileType("Cooling points") { Patterns = ["*.json", "*.csv", "*.txt", "*.dat"] }, FilePickerFileTypes.All],
        });
        var paths = files.Select(f => f.TryGetLocalPath()).OfType<string>().ToList();
        if (paths.Count > 0) Vm.OpenFinishedGlass(paths);
    }
    private void OnRefit(object? s, RoutedEventArgs e) => Vm.RefitFinishedGlass();
    private void OnQueue(object? s, Avalonia.Interactivity.RoutedEventArgs e) => Vm.QueueGlass();
}
