using Avalonia.Controls;
using Avalonia.Interactivity;
using Avalonia.Markup.Xaml;
using Avalonia.Platform.Storage;
using CapsStudio.ViewModels;

namespace CapsStudio.Views.Pages;

public partial class ChecksPage : PageBase
{
    public ChecksPage() => AvaloniaXamlLoader.Load(this);

    private void OnAction(object? s, RoutedEventArgs e) { if ((s as Control)?.Tag is FileCheckRow r) Vm.RunCheckAction(r); }
    private void OnContinue(object? s, RoutedEventArgs e) => Vm.SetModule(8);

    private async void OnExport(object? s, RoutedEventArgs e)
    {
        var top = TopLevel.GetTopLevel(this);
        if (top == null) return;
        var f = await top.StorageProvider.SaveFilePickerAsync(new FilePickerSaveOptions
        {
            Title = "Export the check report", SuggestedFileName = "file-checks.md", DefaultExtension = "md",
            FileTypeChoices = [new FilePickerFileType("Markdown") { Patterns = ["*.md"] }],
        });
        if (f?.TryGetLocalPath() is { } path) Vm.ExportChecks(path);
    }
    private void OnApplyType(object? s, Avalonia.Interactivity.RoutedEventArgs e) { if ((s as Control)?.Tag is ViewModels.AtomTypeRow r) Vm.ApplyTypeRow(r); }
}
