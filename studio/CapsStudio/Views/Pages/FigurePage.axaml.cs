using Avalonia.Controls;
using Avalonia.Interactivity;
using Avalonia.Markup.Xaml;
using Avalonia.Platform.Storage;
using CapsStudio.ViewModels;

namespace CapsStudio.Views.Pages;

public partial class FigurePage : PageBase
{
    public FigurePage() => AvaloniaXamlLoader.Load(this);

    private void OnTile(object? s, RoutedEventArgs e) { if ((s as Control)?.Tag is FigureTile t) Vm.FigBackground = t.Background; }
    private void OnBack(object? s, RoutedEventArgs e) => Vm.SetModule(8);

    private async void OnExport(object? s, RoutedEventArgs e)
    {
        var top = TopLevel.GetTopLevel(this);
        if (top == null || Vm.Document == null) return;
        var ext = Vm.FigExtension;
        var bg = MainViewModel.Backgrounds[Vm.FigBackground].ToLowerInvariant();
        var f = await top.StorageProvider.SaveFilePickerAsync(new FilePickerSaveOptions
        {
            Title = "Export figure", DefaultExtension = ext,
            SuggestedFileName = $"{Path.GetFileNameWithoutExtension(Vm.Document.Path)}_{bg}.{ext}",
            FileTypeChoices = [new FilePickerFileType(ext.ToUpperInvariant()) { Patterns = [$"*.{ext}"] }],
        });
        if (f?.TryGetLocalPath() is not { } path) return;
        try { await Vm.ExportFigure(path, FigureDrawing.SavePng, FigureDrawing.AddToSvg, FigureDrawing.Compose); }
        catch (Exception ex) { Vm.Status = "Export failed: " + ex.Message; }
    }
}
