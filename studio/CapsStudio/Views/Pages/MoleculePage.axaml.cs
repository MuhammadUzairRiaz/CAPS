using Avalonia.Controls;
using Avalonia.Controls.Primitives;
using Avalonia.Interactivity;
using Avalonia.Markup.Xaml;
using Avalonia.Platform.Storage;
using CapsStudio.ViewModels;

namespace CapsStudio.Views.Pages;

public partial class MoleculePage : PageBase
{
    private readonly SketchCanvas _sketch;
    private readonly MolView _preview;

    public MoleculePage()
    {
        AvaloniaXamlLoader.Load(this);
        _sketch = this.FindControl<SketchCanvas>("Sketch")!;
        _preview = this.FindControl<MolView>("Preview")!;
        _sketch.Edited += smiles => Vm.SmilesFromSketch(smiles);
        DataContextChanged += (_, _) =>
        {
            if (DataContext is not MainViewModel vm) return;
            vm.SketchGraphChanged += json => _sketch.Load(json);
            vm.MolViewChanged += () => _preview.Document = vm.MolDoc;
            vm.PropertyChanged += (_, e) =>
            {
                if (e.PropertyName == nameof(MainViewModel.MolDoc)) { _preview.Document = vm.MolDoc; _preview.Reset(); }
            };
        };
    }

    private void OnTool(object? s, RoutedEventArgs e)
    {
        if (s is not ToggleButton b) return;
        foreach (var n in new[] { "ToolDraw", "ToolRing", "ToolCharge", "ToolErase" })
            if (this.FindControl<ToggleButton>(n) is { } t) t.IsChecked = t == b;
        _sketch.Mode = (b.Tag as string) switch
        {
            "ring" => SketchCanvas.Tool.Ring,
            "erase" => SketchCanvas.Tool.Erase,
            "charge" => SketchCanvas.Tool.Charge,
            _ => SketchCanvas.Tool.Draw,
        };
    }

    private void OnElement(object? s, RoutedEventArgs e)
    {
        if (s is not ToggleButton b || !int.TryParse(b.Tag as string, out var z)) return;
        foreach (var t in this.FindControl<StackPanel>("Elements")!.Children.OfType<ToggleButton>()) t.IsChecked = t == b;
        _sketch.Element = z;
        // choosing an element means drawing with it
        if (this.FindControl<ToggleButton>("ToolDraw") is { } d && d.IsChecked != true) OnTool(d, e);
    }

    private void OnClear(object? s, RoutedEventArgs e)
    {
        _sketch.Load("");
        Vm.SmilesFromSketch("");
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

    private void OnOpenInStudio(object? s, RoutedEventArgs e)
    {
        _preview.Document = null;
        Vm.OpenMoleculeInStudio();
    }
}
