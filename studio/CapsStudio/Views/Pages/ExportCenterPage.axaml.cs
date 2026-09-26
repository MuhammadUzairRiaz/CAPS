using Avalonia.Controls;
using Avalonia.Interactivity;
using Avalonia.Markup.Xaml;
using Avalonia.Platform.Storage;

namespace CapsStudio.Views.Pages;

public partial class ExportCenterPage : PageBase
{
    public ExportCenterPage() => AvaloniaXamlLoader.Load(this);

    private void OnForceField(object? s, RoutedEventArgs e) => Vm.SetModule(7);
    private void OnExportData(object? s, RoutedEventArgs e) => Vm.OpenExport();
    private void OnExportFigure(object? s, RoutedEventArgs e) => Vm.SetModule(18);
    private async void OnWrite(object? s, RoutedEventArgs e) => await Vm.WriteEngines();

    private async void OnPickFolder(object? s, RoutedEventArgs e)
    {
        var top = TopLevel.GetTopLevel(this);
        if (top == null) return;
        var dirs = await top.StorageProvider.OpenFolderPickerAsync(new FolderPickerOpenOptions { Title = "Folder for the LAMMPS and GROMACS files", AllowMultiple = false });
        if (dirs.Count > 0 && dirs[0].TryGetLocalPath() is { } path) Vm.EngineFolder = path;
    }
}
