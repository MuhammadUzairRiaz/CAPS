using Avalonia;
using Avalonia.Controls;
using Avalonia.Input;
using Avalonia.Input.Platform;
using Avalonia.Interactivity;
using Avalonia.Markup.Xaml;
using Avalonia.Media;
using Avalonia.Platform.Storage;
using Avalonia.Threading;
using Avalonia.VisualTree;
using CapsStudio.ViewModels;

namespace CapsStudio.Views.Pages;

public partial class ExportDialog : PageBase
{
    private MainViewModel? _hooked;
    private readonly DispatcherTimer _debounce = new() { Interval = TimeSpan.FromMilliseconds(120) };

    public ExportDialog()
    {
        AvaloniaXamlLoader.Load(this);
        _debounce.Tick += (_, _) => { _debounce.Stop(); RefreshPreview(); };
        DataContextChanged += (_, _) =>
        {
            if (DataContext is not MainViewModel vm || vm == _hooked) return;
            _hooked = vm;
            vm.ExportPreviewChanged += () => { _debounce.Stop(); _debounce.Start(); };
        };
    }

    /// <summary>Renders the preview now (also used by screenshots).</summary>
    public void RefreshPreview()
    {
        if (DataContext is not MainViewModel vm || !vm.ExportDialogOpen) return;
        var frame = this.FindControl<Border>("PreviewFrame")!;
        var fw = Math.Max(200, (int)(frame.Bounds.Width > 0 ? frame.Bounds.Width : 598) - 2);
        var fh = Math.Max(120, (int)(frame.Bounds.Height > 0 ? frame.Bounds.Height : 340) - 2);
        var r = vm.ExportPreview(fw, fh);
        if (r is not { } p) return;
        var rgba = p.Rgba;
        if (vm.ExportLabels && vm.ExportIsImage)
        {
            var (labels, line, measure) = vm.ExportMarks(p.W, p.H);
            if (labels.Count > 0 || line != null)
            {
                var layer = FigureDrawing.MarksLayer(p.W, p.H, labels, line, measure, vm.ExportDlgBackground == 0);
                Blend(rgba, layer);
            }
        }
        var img = this.FindControl<Image>("Preview")!;
        (img.Source as IDisposable)?.Dispose();
        img.Source = MainViewModel.ToBitmap(rgba, p.W, p.H);
        this.FindControl<Border>("CheckerLayer")!.Background = vm.ExportDlgBackground == 2 ? Checker() : Brushes.Transparent;
    }

    private static void Blend(byte[] dst, byte[] src)
    {
        for (int k = 0; k < dst.Length; k += 4)
        {
            var a = src[k + 3] / 255.0;
            if (a == 0) continue;
            var da = dst[k + 3] / 255.0;
            var oa = a + da * (1 - a);
            for (int c = 0; c < 3; c++) dst[k + c] = (byte)Math.Round((src[k + c] * a + dst[k + c] * da * (1 - a)) / oa);
            dst[k + 3] = (byte)Math.Round(oa * 255);
        }
    }

    private static IBrush Checker()
    {
        var tile = new DrawingGroup();
        tile.Children.Add(new GeometryDrawing { Brush = Brushes.White, Geometry = new RectangleGeometry(new Rect(0, 0, 16, 16)) });
        var grey = new SolidColorBrush(Color.FromRgb(0xE9, 0xE9, 0xE9));
        tile.Children.Add(new GeometryDrawing { Brush = grey, Geometry = new RectangleGeometry(new Rect(0, 0, 8, 8)) });
        tile.Children.Add(new GeometryDrawing { Brush = grey, Geometry = new RectangleGeometry(new Rect(8, 8, 8, 8)) });
        return new DrawingBrush(tile) { TileMode = TileMode.Tile, DestinationRect = new RelativeRect(0, 0, 16, 16, RelativeUnit.Absolute), Stretch = Stretch.None };
    }

    private void OnBackdrop(object? s, PointerPressedEventArgs e) { if (Vm.ExportDlgIdle) Vm.ExportDialogOpen = false; }
    private void OnClose(object? s, RoutedEventArgs e) { Vm.CancelExport(); Vm.ExportDialogOpen = false; }
    private void OnStop(object? s, RoutedEventArgs e) => Vm.CancelExport();
    private void OnPreset(object? s, RoutedEventArgs e)
    {
        if (s is Button { Tag: string name }) Vm.ExportPreset = Array.FindIndex(MainViewModel.ExportPresets, p => p.Name == name);
        UpdatePresetClasses();
    }

    protected override void OnPropertyChanged(AvaloniaPropertyChangedEventArgs change)
    {
        base.OnPropertyChanged(change);
        if (change.Property == IsVisibleProperty && IsVisible) Dispatcher.UIThread.Post(() => { UpdatePresetClasses(); RefreshPreview(); });
    }

    private void UpdatePresetClasses()
    {
        if (DataContext is not MainViewModel vm) return;
        foreach (var b in this.GetVisualDescendants().OfType<Button>().Where(b => b.Classes.Contains("preset")))
            b.Classes.Set("on", b.Tag is string n && vm.ExportPreset >= 0 && MainViewModel.ExportPresets[vm.ExportPreset].Name == n);
    }

    private byte[]? Marks(int w, int h)
    {
        var (labels, line, measure) = Vm.ExportMarks(w, h);
        return labels.Count > 0 || line != null ? FigureDrawing.MarksLayer(w, h, labels, line, measure, Vm.ExportDlgBackground == 0) : null;
    }

    private async void OnExport(object? s, RoutedEventArgs e)
    {
        var top = TopLevel.GetTopLevel(this);
        if (top == null) return;
        string? path;
        if (Vm.ExportIsMovie && Vm.MovieFormat == 1)
        {
            var dirs = await top.StorageProvider.OpenFolderPickerAsync(new FolderPickerOpenOptions { Title = "Folder for the PNG sequence" });
            path = dirs.Count > 0 ? dirs[0].TryGetLocalPath() : null;
            if (path != null) path = Path.Combine(path, Vm.ExportDefaultName);
        }
        else
        {
            var ext = Path.GetExtension(Vm.ExportDefaultName).TrimStart('.');
            var file = await top.StorageProvider.SaveFilePickerAsync(new FilePickerSaveOptions
            {
                Title = Vm.ExportButton, SuggestedFileName = Vm.ExportDefaultName, DefaultExtension = ext,
                FileTypeChoices = [new FilePickerFileType(ext.ToUpperInvariant()) { Patterns = ["*." + ext] }],
            });
            path = file?.TryGetLocalPath();
        }
        if (path == null) return;
        if (Vm.ExportIsMovie) await Vm.ExportDialogMovie(path);
        else await Vm.ExportDialogImage(path, Marks);
    }

    /// <summary>The image at full size on the clipboard.</summary>
    private async void OnCopy(object? s, RoutedEventArgs e)
    {
        var top = TopLevel.GetTopLevel(this);
        if (top?.Clipboard == null) return;
        var path = Path.Combine(Path.GetTempPath(), "caps-export-" + Guid.NewGuid().ToString("N")[..8] + ".png");
        if (await Vm.ExportDialogImage(path, Marks) == null) return;
        using var bmp = new Avalonia.Media.Imaging.Bitmap(path);
        await top.Clipboard.SetBitmapAsync(bmp);
        Vm.Status = "Copied the image to the clipboard";
    }
}
