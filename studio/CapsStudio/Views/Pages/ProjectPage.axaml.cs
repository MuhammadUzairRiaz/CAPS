using Avalonia.Controls;
using Avalonia.Input;
using Avalonia.Input.Platform;
using Avalonia.Interactivity;
using Avalonia.Markup.Xaml;
using Avalonia.Platform.Storage;
using CapsStudio.ViewModels;

namespace CapsStudio.Views.Pages;

public partial class ProjectPage : PageBase
{
    public ProjectPage() { AvaloniaXamlLoader.Load(this); }

    private void OnDoc(object? s, RoutedEventArgs e) { if (s is Button { Tag: ProjectDoc d }) Vm.SelectProjectDoc(d); }
    private void OnDocOpen(object? s, TappedEventArgs e)
    {
        if (s is Button { Tag: ProjectDoc d }) { Vm.SetModule(8); Vm.Open(d.Path, null); }
    }
    private void OnNew(object? s, RoutedEventArgs e) => Window?.OpenCommand.Execute(null);

    private async void OnFolder(object? s, RoutedEventArgs e)
    {
        var top = TopLevel.GetTopLevel(this);
        if (top == null) return;
        var dirs = await top.StorageProvider.OpenFolderPickerAsync(new FolderPickerOpenOptions { Title = "Project folder" });
        if (dirs.Count > 0 && dirs[0].TryGetLocalPath() is { } p) Vm.OpenProject(p);
    }

    private async void OnShare(object? s, RoutedEventArgs e)
    {
        var top = TopLevel.GetTopLevel(this);
        if (top == null) return;
        var file = await top.StorageProvider.SaveFilePickerAsync(new FilePickerSaveOptions
        {
            Title = "Share project", SuggestedFileName = Vm.ProjectName + ".zip", DefaultExtension = "zip",
            FileTypeChoices = [new FilePickerFileType("Zip") { Patterns = ["*.zip"] }],
        });
        if (file?.TryGetLocalPath() is { } p)
            try { Vm.Status = Vm.ShareProject(p); } catch (Exception ex) { Vm.Status = "Could not write the zip: " + ex.Message; }
    }

    private async void OnCopyMethods(object? s, RoutedEventArgs e)
    {
        if (TopLevel.GetTopLevel(this)?.Clipboard is not { } cb) return;
        await cb.SetTextAsync(Vm.ProjectMethodsWithRefs);
        Vm.Status = "Copied the methods section and its references";
    }

    private async void OnExportTable(object? s, RoutedEventArgs e)
    {
        if (DataContext is not MainViewModel vm || TopLevel.GetTopLevel(this) is not { } top) return;
        var f = await top.StorageProvider.SaveFilePickerAsync(new Avalonia.Platform.Storage.FilePickerSaveOptions
        {
            Title = "Export the results table", SuggestedFileName = "results.csv", DefaultExtension = "csv",
            FileTypeChoices = [new Avalonia.Platform.Storage.FilePickerFileType("CSV") { Patterns = ["*.csv"] }],
        });
        if (f?.TryGetLocalPath() is { } path) { try { vm.ExportProjectTable(path); } catch (Exception ex) { vm.Status = "Could not write the table: " + ex.Message; } }
    }

    private async void OnBibtex(object? s, RoutedEventArgs e)
    {
        var top = TopLevel.GetTopLevel(this);
        if (top == null) return;
        var file = await top.StorageProvider.SaveFilePickerAsync(new FilePickerSaveOptions
        {
            Title = "BibTeX", SuggestedFileName = "references.bib", DefaultExtension = "bib",
            FileTypeChoices = [new FilePickerFileType("BibTeX") { Patterns = ["*.bib"] }],
        });
        if (file?.TryGetLocalPath() is { } p) { await File.WriteAllTextAsync(p, Vm.ProjectBibtex()); Vm.Status = $"Wrote {Path.GetFileName(p)}"; }
    }
}
