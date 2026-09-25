using Avalonia.Controls;
using Avalonia.Interactivity;
using Avalonia.Markup.Xaml;
using Avalonia.Platform.Storage;

namespace CapsStudio.Views.Pages;

public partial class GrowPage : PageBase
{
    public GrowPage() => AvaloniaXamlLoader.Load(this);
    /// <summary>Where the 3D view goes while this page shows.</summary>
    public Decorator Slot => this.FindControl<Decorator>("ViewSlot")!;

    private async void OnGrow(object? s, RoutedEventArgs e) => await Vm.Grow();
    private void OnCancel(object? s, RoutedEventArgs e) => Vm.CancelGrow();
    private void OnDensityMode(object? s, RoutedEventArgs e) => Vm.GrowUseBox = false;
    private void OnDisplay(object? s, RoutedEventArgs e) { if (s is Control { Tag: string t }) Vm.StyleIndex = int.Parse(t); }
    private async void OnSave(object? s, RoutedEventArgs e) { if (Window != null) await Window.SaveAsAsync("data", "LAMMPS data"); }
    private async void OnCopyCommand(object? s, RoutedEventArgs e)
    {
        var clip = TopLevel.GetTopLevel(this)?.Clipboard;
        if (clip != null) await clip.SetTextAsync(Vm.GrowCommand);
        Vm.Status = "Copied: " + Vm.GrowCommand;
    }
    private async void OnCopyPython(object? s, RoutedEventArgs e)
    {
        var clip = TopLevel.GetTopLevel(this)?.Clipboard;
        if (clip != null) await clip.SetTextAsync(Vm.GrowPython());
        Vm.Status = "Copied the build as Python (import caps; caps.polymer(...))";
    }

    private async void OnSaveRecipe(object? s, RoutedEventArgs e)
    {
        var top = TopLevel.GetTopLevel(this);
        if (top == null) return;
        var file = await top.StorageProvider.SaveFilePickerAsync(new Avalonia.Platform.Storage.FilePickerSaveOptions
        {
            Title = "Save the build as a recipe",
            SuggestedFileName = "grow_recipe.yaml",
            DefaultExtension = "yaml",
            FileTypeChoices = [new Avalonia.Platform.Storage.FilePickerFileType("CAPS recipe") { Patterns = ["*.yaml", "*.yml"] }],
        });
        if (file?.TryGetLocalPath() is not string path) return;
        try
        {
            System.IO.File.WriteAllText(path, Vm.GrowRecipe());
            Vm.Status = $"Saved {System.IO.Path.GetFileName(path)} · run it with caps run {System.IO.Path.GetFileName(path)} or Start › From a recipe";
        }
        catch (System.Exception ex) { Vm.Status = "Could not save the recipe: " + ex.Message; }
    }

    private void OnChoosePolymer(object? s, Avalonia.Interactivity.RoutedEventArgs e) => Vm.SetModule(13);
    private void OnUsePs(object? s, Avalonia.Interactivity.RoutedEventArgs e) => Vm.UsePolystyreneInGrow();
    private void OnDispersity(object? s, Avalonia.Interactivity.RoutedEventArgs e) => Vm.OpenPolydispersity();
    private void OnClearLengths(object? s, Avalonia.Interactivity.RoutedEventArgs e) => Vm.ClearPdLengths();
    private void OnClearStereo(object? s, Avalonia.Interactivity.RoutedEventArgs e) => Vm.ClearGrowStereo();
}
