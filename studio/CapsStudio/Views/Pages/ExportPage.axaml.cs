using Avalonia.Controls;
using Avalonia.Interactivity;
using Avalonia.Markup.Xaml;
using Avalonia.Platform.Storage;
using CapsStudio.ViewModels;

namespace CapsStudio.Views.Pages;

public partial class ExportPage : PageBase
{
    public ExportPage() => AvaloniaXamlLoader.Load(this);

    private void OnBack(object? s, RoutedEventArgs e) => Vm.SetModule(Vm.ExportPipeline ? 20 : 8);

    private async void OnExport(object? s, RoutedEventArgs e)
    {
        var top = TopLevel.GetTopLevel(this);
        if (top == null) return;
        var fmt = MainViewModel.ExportFormats[Vm.ExportFormatIndex];
        var f = await top.StorageProvider.SaveFilePickerAsync(new FilePickerSaveOptions
        {
            Title = "Export " + fmt.Name, SuggestedFileName = Vm.ExportSuggestedName, DefaultExtension = fmt.Extension,
            FileTypeChoices = [new FilePickerFileType(fmt.Name) { Patterns = [$"*.{fmt.Extension}"] }],
        });
        if (f?.TryGetLocalPath() is { } path) await Vm.ExportNow(path);
    }
}
