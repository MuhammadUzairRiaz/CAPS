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
    private void OnOpenGmxSet(object? s, RoutedEventArgs e) => Vm.OpenGmxSet();

    private async void OnTopology(object? s, RoutedEventArgs e)
    {
        var top = TopLevel.GetTopLevel(this);
        if (top == null) return;
        var files = await top.StorageProvider.OpenFilePickerAsync(new FilePickerOpenOptions
        {
            Title = "Topology or structure (LAMMPS data, GROMACS .top/.gro, PDB …)", AllowMultiple = false,
            FileTypeFilter =
            [
                new FilePickerFileType("Topologies and structures") { Patterns = ["*.data", "*.lmp", "*.top", "*.itp", "*.prmtop", "*.parm7", "*.gro", "*.pdb", "*.mol2", "*.xyz"] },
                new FilePickerFileType("All files") { Patterns = ["*"] },
            ],
        });
        if (files.Count > 0 && files[0].TryGetLocalPath() is { } p) Vm.SetOpenTopology(p);
    }
}
