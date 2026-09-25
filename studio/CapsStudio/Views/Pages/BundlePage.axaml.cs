using Avalonia.Controls;
using Avalonia.Interactivity;
using Avalonia.Markup.Xaml;
using Avalonia.Platform.Storage;

namespace CapsStudio.Views.Pages;

public partial class BundlePage : PageBase
{
    public BundlePage() => AvaloniaXamlLoader.Load(this);

    private void OnBack(object? s, RoutedEventArgs e) => Vm.SetModule(Vm.PipelineRows.Count > 0 ? 20 : 8);

    private async void OnExport(object? s, RoutedEventArgs e)
    {
        var top = TopLevel.GetTopLevel(this);
        if (top == null) return;
        var f = await top.StorageProvider.SaveFilePickerAsync(new FilePickerSaveOptions
        {
            Title = "Export figure bundle", SuggestedFileName = Vm.BundleFileName, DefaultExtension = "zip",
            FileTypeChoices = [new FilePickerFileType("CAPS bundle (zip)") { Patterns = ["*.zip"] }],
        });
        if (f?.TryGetLocalPath() is { } path) await Vm.ExportBundle(path);
    }
}
