using Avalonia.Controls;
using Avalonia.Controls.Primitives;
using Avalonia.Interactivity;
using Avalonia.Markup.Xaml;
using Avalonia.Platform.Storage;
using CapsStudio.ViewModels;

namespace CapsStudio.Views.Pages;

public partial class MoleculePage : PageBase
{
    private void OnUseAsUnit(object? s, Avalonia.Interactivity.RoutedEventArgs e) => Vm.UseSketchAsRepeatUnit();

    public MoleculePage()
    {
        AvaloniaXamlLoader.Load(this);
    }

    private async void OnSave(object? s, RoutedEventArgs e)
    {
        if (Window == null) return;
        var file = await Window.StorageProvider.SaveFilePickerAsync(new FilePickerSaveOptions
        {
            Title = "Save the molecule",
            SuggestedFileName = "molecule.mol2",
            DefaultExtension = "mol2",
            FileTypeChoices =
            [
                new FilePickerFileType("Tripos mol2") { Patterns = ["*.mol2"] },
                new FilePickerFileType("PDB") { Patterns = ["*.pdb"] },
                new FilePickerFileType("XYZ") { Patterns = ["*.xyz"] },
                new FilePickerFileType("LAMMPS data") { Patterns = ["*.data"] },
            ],
        });
        if (file?.TryGetLocalPath() is not { } path) return;
        try { Vm.SaveMolecule(path); }
        catch (Exception ex) { Vm.Status = "Could not save: " + ex.Message; }
    }

    private async void OnPack(object? s, RoutedEventArgs e)
    {
        if (Window == null) return;
        var file = await Window.StorageProvider.SaveFilePickerAsync(new FilePickerSaveOptions
        {
            Title = "Save the molecule for Pack",
            SuggestedFileName = "molecule.pdb",
            DefaultExtension = "pdb",
            FileTypeChoices = [new FilePickerFileType("PDB") { Patterns = ["*.pdb"] }],
        });
        if (file?.TryGetLocalPath() is not { } path) return;
        try
        {
            Vm.SaveMolecule(path);
            Vm.AddPackStructure(path);
            Vm.SetModule(5);
        }
        catch (Exception ex) { Vm.Status = "Could not save: " + ex.Message; }
    }

    private async void OnBuildBeads(object? sender, Avalonia.Interactivity.RoutedEventArgs e) { if (DataContext is MainViewModel vm) await vm.BuildBeadsMolecule(); }
    void OnInsertIntoDocument(object? s, RoutedEventArgs e) => Vm.InsertMoleculeIntoDocument();
    void OnOpenInStudio(object? s, RoutedEventArgs e) => Vm.OpenMoleculeInStudio();
    private async void OnAnalogs(object? s, RoutedEventArgs e) => await Vm.BuildAnalogs();
}
