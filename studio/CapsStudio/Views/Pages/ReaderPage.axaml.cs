using Avalonia.Controls;
using Avalonia.Input;
using Avalonia.Interactivity;
using Avalonia.Markup.Xaml;
using Avalonia.Media;
using Avalonia.Platform.Storage;
using CapsStudio.ViewModels;

namespace CapsStudio.Views.Pages;

public partial class ReaderPage : PageBase
{
    private MainViewModel? _hooked;

    public ReaderPage()
    {
        AvaloniaXamlLoader.Load(this);
        DataContextChanged += (_, _) =>
        {
            if (DataContext is not MainViewModel vm || vm == _hooked) return;
            _hooked = vm;
            vm.ReaderChanged += () =>
            {
                var p = this.FindControl<LinePlot>("Plot")!;
                p.XLabel = vm.ReaderXLabel;
                p.YLabel = vm.ReaderYLabel;
                p.RefY = null;
                if (vm.ReaderFitLine.Length > 0) p.SetWithFit(vm.ReaderCurve, vm.ReaderFitLine);
                else p.SetData(vm.ReaderCurve);
            };
            vm.ReaderSelect += (start, length) =>
            {
                var b = this.FindControl<TextBox>("Body")!;
                var n = b.Text?.Length ?? 0;
                if (start < 0) { b.CaretIndex = n; return; }
                b.SelectionStart = Math.Min(start, n);
                b.SelectionEnd = Math.Min(start + length, n);
                b.CaretIndex = Math.Min(start + length, n);
            };
            vm.PropertyChanged += (_, e) =>
            {
                if (e.PropertyName == nameof(MainViewModel.ReaderWrap))
                    this.FindControl<TextBox>("Body")!.TextWrapping = vm.ReaderWrap ? TextWrapping.Wrap : TextWrapping.NoWrap;
            };
        };
    }

    private void OnReload(object? s, RoutedEventArgs e) => Vm.ReaderReload();
    private void OnNext(object? s, RoutedEventArgs e) => Vm.ReaderFindNext();
    private void OnPrev(object? s, RoutedEventArgs e) => Vm.ReaderFindPrev();
    private void OnGoTo(object? s, RoutedEventArgs e) => Vm.ReaderGoToTable();
    private void OnFindKey(object? s, KeyEventArgs e)
    {
        if (e.Key != Key.Enter) return;
        if (e.KeyModifiers.HasFlag(KeyModifiers.Shift)) Vm.ReaderFindPrev(); else Vm.ReaderFindNext();
        e.Handled = true;
    }
    private async void OnCopyCsv(object? s, RoutedEventArgs e)
    {
        var clip = TopLevel.GetTopLevel(this)?.Clipboard;
        if (clip != null) await clip.SetTextAsync(Vm.ReaderTableCsv());
        Vm.Status = "Copied the table as CSV";
    }
    private async void OnOpen(object? s, RoutedEventArgs e)
    {
        var top = TopLevel.GetTopLevel(this);
        if (top == null) return;
        var files = await top.StorageProvider.OpenFilePickerAsync(new FilePickerOpenOptions
        {
            Title = "Read a text output",
            AllowMultiple = false,
            FileTypeFilter =
            [
                new FilePickerFileType("Logs, inputs and tables") { Patterns = ["*.log", "log.*", "*.in", "in.*", "*.mdp", "*.xvg", "*.csv", "*.tsv", "*.txt", "*.md", "*.out"] },
                new FilePickerFileType("All files") { Patterns = ["*"] },
            ],
        });
        if (files.Count > 0 && files[0].TryGetLocalPath() is { } p) Vm.OpenReader(p);
    }
}
