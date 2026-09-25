using Avalonia.Controls;
using Avalonia.Interactivity;
using Avalonia.Markup.Xaml;
using Avalonia.Platform.Storage;

namespace CapsStudio.Views.Pages;

public partial class SavePipelinePage : PageBase
{
    public SavePipelinePage() => AvaloniaXamlLoader.Load(this);

    private void OnBack(object? s, RoutedEventArgs e) => Vm.SetModule(20);

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
