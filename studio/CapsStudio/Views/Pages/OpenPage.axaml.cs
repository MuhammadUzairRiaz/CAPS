using Avalonia.Controls;
using Avalonia.Interactivity;
using Avalonia.Markup.Xaml;
using Avalonia.Platform.Storage;

namespace CapsStudio.Views.Pages;

public partial class OpenPage : PageBase
{
    public OpenPage() => AvaloniaXamlLoader.Load(this);

    private void OnCancel(object? s, RoutedEventArgs e) => Vm.CancelOpen();
    private void OnOpen(object? s, RoutedEventArgs e) => Vm.ConfirmOpen();

    private async void OnTopology(object? s, RoutedEventArgs e)
    {
        var top = TopLevel.GetTopLevel(this);
        if (top == null) return;
        var files = await top.StorageProvider.OpenFilePickerAsync(new FilePickerOpenOptions
        {
            Title = "Topology (LAMMPS data)", AllowMultiple = false,
            FileTypeFilter = [new FilePickerFileType("LAMMPS data") { Patterns = ["*.data", "*.lmp"] }],
        });
        if (files.Count > 0 && files[0].TryGetLocalPath() is { } p) Vm.SetOpenTopology(p);
    }
}
