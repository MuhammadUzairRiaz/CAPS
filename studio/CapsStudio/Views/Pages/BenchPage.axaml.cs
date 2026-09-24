using System.ComponentModel;
using Avalonia;
using Avalonia.Controls;
using Avalonia.Interactivity;
using Avalonia.Layout;
using Avalonia.Markup.Xaml;
using Avalonia.Media;
using Avalonia.Platform.Storage;
using Avalonia.Threading;
using CapsStudio.ViewModels;

namespace CapsStudio.Views.Pages;

public partial class BenchPage : PageBase
{
    private BenchItem? _shown;

    public BenchPage()
    {
        AvaloniaXamlLoader.Load(this);
        DataContextChanged += (_, _) =>
        {
            if (DataContext is not MainViewModel vm) return;
            vm.LoadBench();
            vm.BenchSelected += () => Dispatcher.UIThread.Post(ShowTable);
            ShowTable();
        };
    }

    private void OnItem(object? s, PropertyChangedEventArgs e) => Dispatcher.UIThread.Post(ShowTable);

    /// <summary>The chosen table's rows as a grid (its columns vary), its note, its time.</summary>
    private void ShowTable()
    {
        var t = (DataContext as MainViewModel)?.SelectedBench;
        if (_shown != t)
        {
            if (_shown != null) _shown.PropertyChanged -= OnItem;
            _shown = t;
            if (t != null) t.PropertyChanged += OnItem;
        }
        var detail = this.FindControl<Grid>("Detail")!;
        detail.Children.Clear();
        detail.RowDefinitions.Clear();
        detail.ColumnDefinitions.Clear();
        this.FindControl<TextBlock>("DetailTitle")!.Text = t == null ? "" : $"Table {t.Id[1..]} · {t.Title}";
        this.FindControl<TextBlock>("Seconds")!.Text = t?.SecondsText ?? "";
        this.FindControl<Button>("RunOne")!.IsVisible = t?.Runnable == true;
        var note = t?.Note ?? "";
        this.FindControl<TextBlock>("Note")!.Text = note;
        this.FindControl<Border>("NoteBox")!.IsVisible = note.Length > 0;
        if (t == null) return;
        var cols = t.Columns;
        var rows = t.Rows;
        if (cols.Length == 0 || rows.Count == 0)
        {
            detail.Children.Add(new TextBlock
            {
                Text = t.IsRunning ? "Running…" : t.Runnable ? "Not run yet. Run the table, or the whole suite." : "Not part of the built-in suite.",
                Classes = { "dim" }, FontSize = 12.5, Margin = new Thickness(0, 6),
            });
            return;
        }
        for (var c = 0; c < cols.Length; c++) detail.ColumnDefinitions.Add(new ColumnDefinition(c == 0 ? GridLength.Star : GridLength.Auto));
        detail.RowDefinitions.Add(new RowDefinition(GridLength.Auto));
        for (var c = 0; c < cols.Length; c++)
        {
            var h = new TextBlock { Text = cols[c].ToUpperInvariant(), Classes = { "h3" }, FontSize = 10, Margin = new Thickness(c == 0 ? 0 : 16, 0, 0, 8) };
            if (c > 0) h.HorizontalAlignment = HorizontalAlignment.Right;
            Grid.SetColumn(h, c);
            detail.Children.Add(h);
        }
        for (var r = 0; r < rows.Count; r++)
        {
            detail.RowDefinitions.Add(new RowDefinition(new GridLength(32)));
            var band = new Border { Background = Tokens.Brush("Bg2B"), CornerRadius = new CornerRadius(5), Opacity = r % 2 == 0 ? 0.55 : 0 };
            Grid.SetRow(band, r + 1);
            Grid.SetColumnSpan(band, cols.Length);
            detail.Children.Add(band);
            var (cells, status) = rows[r];
            for (var c = 0; c < Math.Min(cols.Length, cells.Length); c++)
            {
                var text = cells[c];
                var isStatus = cols[c] == "Status" || cols[c] == "Success";
                var tb = new TextBlock
                {
                    Text = text, VerticalAlignment = VerticalAlignment.Center, FontSize = c == 0 ? 12.5 : 12,
                    Margin = new Thickness(c == 0 ? 8 : 16, 0, c == cols.Length - 1 ? 8 : 0, 0),
                    TextTrimming = TextTrimming.CharacterEllipsis,
                };
                if (c > 0) { tb.FontFamily = Tokens.Mono; tb.HorizontalAlignment = HorizontalAlignment.Right; }
                if (isStatus) tb.Foreground = Tokens.Brush(status == "pass" ? "OkB" : status == "fail" ? "ErrB" : "MutedB");
                Grid.SetRow(tb, r + 1);
                Grid.SetColumn(tb, c);
                detail.Children.Add(tb);
            }
        }
    }

    private async void OnRunAll(object? s, RoutedEventArgs e) => await Vm.RunBench(true);
    private async void OnRunOne(object? s, RoutedEventArgs e) => await Vm.RunBench(false);
    private void OnCancel(object? s, RoutedEventArgs e) => Vm.CancelBench();

    private async void OnExport(object? s, RoutedEventArgs e)
    {
        if (Window == null) return;
        var dirs = await Window.StorageProvider.OpenFolderPickerAsync(new FolderPickerOpenOptions { Title = "Folder for results.md, results.tex and the CSV tables" });
        if (dirs.Count == 0 || dirs[0].TryGetLocalPath() is not { } dir) return;
        try { Vm.ExportBench(dir); }
        catch (Exception ex) { Vm.Status = "Could not export: " + ex.Message; }
    }
}
