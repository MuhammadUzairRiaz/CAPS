using Avalonia.Controls;
using Avalonia.Interactivity;
using Avalonia.Markup.Xaml;
using Avalonia.Platform.Storage;

namespace CapsStudio.Views.Pages;

public partial class SavePipelinePage : PageBase
{
    public SavePipelinePage() => AvaloniaXamlLoader.Load(this);

    private void OnBack(object? s, RoutedEventArgs e) => Vm.SetModule(20);
    private void OnSaveScope(object? s, RoutedEventArgs e) { if (Vm.SavePipelineToScope() != null) Vm.SetModule(20); }
    private void OnAddOutput(object? s, RoutedEventArgs e) { if ((s as Control)?.Tag is string kind) Vm.AddPipelineOutput(kind); }
    private void OnRemoveOutput(object? s, RoutedEventArgs e) { if ((s as Control)?.Tag is CapsStudio.ViewModels.PipelineOutputRow row) Vm.RemovePipelineOutput(row); }
    private async void OnWriteOutputs(object? s, RoutedEventArgs e)
    {
        var top = TopLevel.GetTopLevel(this);
        if (top == null) return;
        var dirs = await top.StorageProvider.OpenFolderPickerAsync(new FolderPickerOpenOptions { Title = "Write the outputs into", AllowMultiple = false });
        if (dirs.Count > 0 && dirs[0].TryGetLocalPath() is { } dir) Vm.WritePipelineOutputs(dir);
    }

    private async void OnSave(object? s, RoutedEventArgs e)
    {
        var top = TopLevel.GetTopLevel(this);
        if (top == null) return;
        var f = await top.StorageProvider.SaveFilePickerAsync(new FilePickerSaveOptions
        {
            Title = "Save pipeline", SuggestedFileName = Vm.PipelineFileName, DefaultExtension = "yaml",
            FileTypeChoices = [new FilePickerFileType("CAPS pipeline (YAML)") { Patterns = ["*.yaml", "*.yml"] }],
        });
        if (f?.TryGetLocalPath() is { } path) { Vm.SavePipelineYaml(path); Vm.SetModule(20); }
    }
}
