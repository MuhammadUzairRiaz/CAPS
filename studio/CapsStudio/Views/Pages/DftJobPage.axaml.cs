using Avalonia.Controls;
using Avalonia.Interactivity;
using Avalonia.Markup.Xaml;
using Avalonia.Platform.Storage;

namespace CapsStudio.Views.Pages;

public partial class DftJobPage : PageBase
{
    public DftJobPage() => AvaloniaXamlLoader.Load(this);

    private async void OnPreview(object? s, RoutedEventArgs e) => await Vm.PreviewJob();
    private async void OnWrite(object? s, RoutedEventArgs e) => await Vm.WriteJob();
    private async void OnCopyCli(object? s, RoutedEventArgs e) { if (TopLevel.GetTopLevel(this)?.Clipboard is { } c) await c.SetTextAsync(Vm.DftCli); }
    private async void OnStructure(object? s, RoutedEventArgs e)
    {
        if (TopLevel.GetTopLevel(this) is not { } top) return;
        var f = await top.StorageProvider.OpenFilePickerAsync(new FilePickerOpenOptions { Title = "The structure (POSCAR, CIF …)", AllowMultiple = false });
        if (f.Count > 0 && f[0].TryGetLocalPath() is { } p) Vm.JobStructure = p;
    }
    private async Task<string?> Folder(string title)
    {
        if (TopLevel.GetTopLevel(this) is not { } top) return null;
        var d = await top.StorageProvider.OpenFolderPickerAsync(new FolderPickerOpenOptions { Title = title });
        return d.Count > 0 ? d[0].TryGetLocalPath() : null;
    }
    private async void OnPaw(object? s, RoutedEventArgs e) { if (await Folder("Your licensed PAW folder (potpaw_PBE)") is { } p) Vm.PotcarDir = p; }
    private async void OnStore(object? s, RoutedEventArgs e) { if (await Folder("The project folder that gets the .store file") is { } p) Vm.SetStore(p); }
}
