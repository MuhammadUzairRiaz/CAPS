using Avalonia.Controls;
using Avalonia.Interactivity;
using Avalonia.Markup.Xaml;
using Avalonia.Platform.Storage;

namespace CapsStudio.Views.Pages;

public partial class ProvenancePage : PageBase
{
    public ProvenancePage() { AvaloniaXamlLoader.Load(this); }

    private async void OnCompare(object? s, RoutedEventArgs e)
    {
        var top = TopLevel.GetTopLevel(this);
        if (top == null) return;
        var files = await top.StorageProvider.OpenFilePickerAsync(new FilePickerOpenOptions { Title = "Compare with the provenance of", AllowMultiple = false });
        if (files.Count > 0 && files[0].TryGetLocalPath() is { } p)
        {
            // the sidecar itself or the file it describes
            if (p.EndsWith(".provenance.json", StringComparison.Ordinal)) p = p[..^".provenance.json".Length];
            Vm.CompareProvenanceWith(p);
        }
    }

    private void OnClearCompare(object? s, RoutedEventArgs e) => Vm.ClearProvenanceCompare();

    private async void OnBibtex(object? s, RoutedEventArgs e)
    {
        var top = TopLevel.GetTopLevel(this);
        if (top == null) return;
        var file = await top.StorageProvider.SaveFilePickerAsync(new FilePickerSaveOptions
        {
            Title = "Export BibTeX", SuggestedFileName = "caps-methods.bib", DefaultExtension = "bib",
            FileTypeChoices = [new FilePickerFileType("BibTeX") { Patterns = ["*.bib"] }],
        });
        if (file?.TryGetLocalPath() is not { } path) return;
        await File.WriteAllTextAsync(path, Vm.ProvenanceBibtex());
        Vm.Status = $"Wrote {Path.GetFileName(path)}";
    }
}
