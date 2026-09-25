using Avalonia.Controls;
using Avalonia.Interactivity;
using Avalonia.Markup.Xaml;
using Avalonia.Platform.Storage;
using CapsStudio.ViewModels;

namespace CapsStudio.Views.Pages;

public partial class ComposerPage : PageBase
{
    private MainViewModel? _hooked;

    public ComposerPage()
    {
        AvaloniaXamlLoader.Load(this);
        var canvas = this.FindControl<FigureCanvas>("Canvas")!;
        canvas.PanelClicked += k => { if (DataContext is MainViewModel vm) vm.ComposerPanel = k; };
        DataContextChanged += (_, _) =>
        {
            if (DataContext is not MainViewModel vm || vm == _hooked) return;
            _hooked = vm;
            vm.ComposerChanged += () =>
            {
                canvas.Composition = vm.Composition;
                canvas.Selected = vm.ComposerPanel;
                canvas.InvalidateVisual();
                this.FindControl<TextBlock>("ScaleText")!.Text = canvas.ScaleText;
            };
        };
    }

    private async System.Threading.Tasks.Task Export(string ext, string kind)
    {
        var top = TopLevel.GetTopLevel(this);
        if (top == null) return;
        var file = await top.StorageProvider.SaveFilePickerAsync(new FilePickerSaveOptions
        {
            Title = "Export the figure", SuggestedFileName = "figure." + ext, DefaultExtension = ext,
            FileTypeChoices = [new FilePickerFileType(kind) { Patterns = ["*." + ext] }],
        });
        if (file?.TryGetLocalPath() is not { } p) return;
        try { Vm.Status = ext == "svg" ? Vm.ExportComposerSvg(p) : Vm.ExportComposerRaster(p); }
        catch (System.Exception ex) { Vm.Status = "Could not export: " + ex.Message; }
    }

    private async void OnSvg(object? s, RoutedEventArgs e) => await Export("svg", "SVG");
    private async void OnTiff(object? s, RoutedEventArgs e) => await Export("tiff", "TIFF");
    private async void OnPng(object? s, RoutedEventArgs e) => await Export("png", "PNG");
    private async void OnPdf(object? s, RoutedEventArgs e) => await Export("pdf", "PDF");
}
