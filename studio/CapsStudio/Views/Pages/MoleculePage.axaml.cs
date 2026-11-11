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
    private readonly SketchCanvas _sketch;
    private readonly MolView _preview;

    public MoleculePage()
    {
        AvaloniaXamlLoader.Load(this);
        _sketch = this.FindControl<SketchCanvas>("Sketch")!;
        _preview = this.FindControl<MolView>("Preview")!;
        _sketch.Edited += smiles => Vm.SmilesFromSketch(smiles);
        _sketch.HistoryChanged += () =>
        {
            this.FindControl<Button>("UndoButton")!.IsEnabled = _sketch.CanUndo;
            this.FindControl<Button>("RedoButton")!.IsEnabled = _sketch.CanRedo;
        };
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

    private static readonly string[] ToolNames = ["ToolDraw", "ToolBond", "ToolRing", "ToolCharge", "ToolMove", "ToolCut", "ToolErase"];

    private void OnTool(object? s, RoutedEventArgs e)
    {
        if (s is not ToggleButton b) return;
        SelectTool(b.Tag as string ?? "draw");
    }

    /// <summary>One tool at a time: its button checked, the canvas in its mode.</summary>
    public void SelectTool(string tag)
    {
        foreach (var n in ToolNames)
            if (this.FindControl<ToggleButton>(n) is { } t) t.IsChecked = (t.Tag as string) == tag;
        _sketch.Mode = tag switch
        {
            "ring" => SketchCanvas.Tool.Ring,
            "erase" => SketchCanvas.Tool.Erase,
            "charge" => SketchCanvas.Tool.Charge,
            "bond" => SketchCanvas.Tool.Bond,
            "cut" => SketchCanvas.Tool.Cut,
            "move" => SketchCanvas.Tool.Move,
            _ => SketchCanvas.Tool.Draw,
        };
    }

    private void OnOrder(object? s, RoutedEventArgs e)
    {
        if (s is not ToggleButton b || !int.TryParse(b.Tag as string, out var n)) return;
        foreach (var t in this.FindControl<StackPanel>("Orders")!.Children.OfType<ToggleButton>()) t.IsChecked = t == b;
        _sketch.BondOrder = n;
        // choosing an order means drawing bonds with it
        if (_sketch.Mode != SketchCanvas.Tool.Bond) SelectTool("bond");
    }

    private void OnRingKind(object? s, SelectionChangedEventArgs e)
    {
        if (_sketch == null || s is not ComboBox c) return;
        // Benzene, then 3- to 8-rings (the 6-ring saturated)
        _sketch.RingAromatic = c.SelectedIndex <= 0;
        _sketch.RingSize = c.SelectedIndex <= 0 ? 6 : c.SelectedIndex + 2;
        if (_sketch.Mode != SketchCanvas.Tool.Ring) SelectTool("ring");
    }

    private void OnElement(object? s, RoutedEventArgs e)
    {
        if (s is not ToggleButton b || !int.TryParse(b.Tag as string, out var z) || z <= 0) return;
        ChooseElement(z);
    }

    /// <summary>The element to draw with: its key checked (an element from the table gets the spare key).</summary>
    public void ChooseElement(int z)
    {
        var keys = this.FindControl<StackPanel>("Elements")!.Children.OfType<ToggleButton>().ToList();
        var other = this.FindControl<ToggleButton>("OtherElement")!;
        var key = keys.FirstOrDefault(k => k != other && k.Tag as string == z.ToString(System.Globalization.CultureInfo.InvariantCulture));
        if (key == null)
        {
            other.Tag = z.ToString(System.Globalization.CultureInfo.InvariantCulture);
            other.Content = z <= MainViewModel.Syms.Length ? MainViewModel.Syms[z - 1] : "?";
            other.IsVisible = true;
            key = other;
        }
        foreach (var t in keys) t.IsChecked = t == key;
        _sketch.Element = z;
        // choosing an element means placing it (the Bond tool places it too)
        if (_sketch.Mode is not (SketchCanvas.Tool.Draw or SketchCanvas.Tool.Bond)) SelectTool("draw");
    }

    /// <summary>Any element: the periodic table in a flyout (the Studio's own layout and colours).</summary>
    private void OnPeriodicTable(object? s, RoutedEventArgs e)
    {
        if (s is not Control anchor || DataContext is not MainViewModel vm) return;
        var flyout = new Flyout { Placement = PlacementMode.Center };
        flyout.FlyoutPresenterClasses.Add("wide");
        var canvas = new Canvas { Width = 738, Height = 420 };
        foreach (var c in vm.PeriodicTable)
        {
            var b = new Button
            {
                Width = 38, Height = 38, MinHeight = 38, Padding = new Avalonia.Thickness(0), CornerRadius = new Avalonia.CornerRadius(5),
                Classes = { "pt" }, Tag = c.Z, [ToolTip.TipProperty] = c.Tip,
                HorizontalContentAlignment = Avalonia.Layout.HorizontalAlignment.Stretch, VerticalContentAlignment = Avalonia.Layout.VerticalAlignment.Stretch,
                Content = new Panel
                {
                    Children =
                    {
                        new Avalonia.Controls.Shapes.Rectangle { Width = 3, HorizontalAlignment = Avalonia.Layout.HorizontalAlignment.Left, Fill = Avalonia.Media.Brush.Parse(c.Colour) },
                        new TextBlock { Text = c.ZText, FontSize = 8, Opacity = 0.75, Margin = new Avalonia.Thickness(6, 2, 0, 0) },
                        new TextBlock { Text = c.Symbol, FontSize = 13.5, FontWeight = Avalonia.Media.FontWeight.SemiBold, HorizontalAlignment = Avalonia.Layout.HorizontalAlignment.Center,
                                        VerticalAlignment = Avalonia.Layout.VerticalAlignment.Bottom, Margin = new Avalonia.Thickness(2, 0, 0, 4) },
                    },
                },
            };
            if (c.Common) b.Classes.Add("common");
            b.Click += (_, _) => { ChooseElement(c.Z); flyout.Hide(); };
            Canvas.SetLeft(b, c.X);
            Canvas.SetTop(b, c.Y);
            canvas.Children.Add(b);
        }
        flyout.Content = new StackPanel
        {
            Spacing = 10,
            Children = { new TextBlock { Text = "Choose an element to draw with", Classes = { "dim" } }, new ScrollViewer { Content = canvas, HorizontalScrollBarVisibility = Avalonia.Controls.Primitives.ScrollBarVisibility.Auto, MaxHeight = 520, MaxWidth = 780 } },
        };
        // over the middle of the sketch: from the key at the toolbar's right end the table ran past the window's edge
        flyout.ShowAt(_sketch.Bounds.Width > 400 ? _sketch : anchor);
    }

    private void OnUndo(object? s, RoutedEventArgs e) => _sketch.Undo();
    private void OnRedo(object? s, RoutedEventArgs e) => _sketch.Undo(redo: true);

    private void OnTidy(object? s, RoutedEventArgs e)
    {
        var smiles = Vm.MolSmiles ?? "";
        if (smiles.Trim().Length == 0) return;
        try { _sketch.Tidy(Interop.CapsDocument.SmilesDepict(smiles)); }
        catch (Exception ex) { Vm.Status = "Could not lay the drawing out: " + ex.Message; }
    }

    private void OnClear(object? s, RoutedEventArgs e) => _sketch.Clear();

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
    void OnOpenInStudio(object? s, RoutedEventArgs e)
    {
        _preview.Document = null;
        Vm.OpenMoleculeInStudio();
    }
    private async void OnAnalogs(object? s, RoutedEventArgs e) => await Vm.BuildAnalogs();
}
