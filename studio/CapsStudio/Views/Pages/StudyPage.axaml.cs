using Avalonia.Controls;
using Avalonia.Interactivity;
using Avalonia.Markup.Xaml;
using Avalonia.Platform.Storage;
using CapsStudio.ViewModels;

namespace CapsStudio.Views.Pages;

public partial class StudyPage : PageBase
{
    private MainViewModel? _hooked;

    public StudyPage()
    {
        AvaloniaXamlLoader.Load(this);
        DataContextChanged += (_, _) =>
        {
            if (DataContext is not MainViewModel vm || vm == _hooked) return;
            _hooked = vm;
            vm.StudyChanged += () =>
            {
                var p = this.FindControl<LinePlot>("Plot")!;
                p.XLabel = vm.StudyXLabel;
                p.YLabel = vm.StudyYLabel;
                p.RefY = null;
                if (vm.StudyFitLine.Length > 0) p.SetWithFit(vm.StudyCurve, vm.StudyFitLine);
                else p.SetData(vm.StudyCurve);
            };
        };
    }

    private void OnAddColumn(object? s, RoutedEventArgs e) => Vm.Status = Vm.StudyAddColumn();
    private void OnAddRow(object? s, RoutedEventArgs e) => Vm.StudyAddEmptyRow();
    private void OnAddOpen(object? s, RoutedEventArgs e) => Vm.Status = Vm.StudyAddOpen();
    private void OnAddBatch(object? s, RoutedEventArgs e) => Vm.Status = Vm.StudyAddBatch();
    private void OnAddReader(object? s, RoutedEventArgs e) => Vm.Status = Vm.StudyAddReader();
    private void OnScores(object? s, RoutedEventArgs e) => Vm.Status = Vm.StudyAddScores();
    private void OnRemoveColumn(object? s, RoutedEventArgs e) { if ((s as Control)?.DataContext is StudyColumn c) Vm.StudyRemoveColumn(c); }
    private void OnRemoveRow(object? s, RoutedEventArgs e) { if ((s as Control)?.DataContext is StudyRow r) Vm.StudyRemoveRow(r); }
    private void OnOpenRow(object? s, RoutedEventArgs e) { if ((s as Control)?.DataContext is StudyRow r) Vm.StudyOpenRow(r); }

    private async void OnPaste(object? s, RoutedEventArgs e)
    {
        var clip = TopLevel.GetTopLevel(this)?.Clipboard;
        var text = clip == null ? null : await clip.GetTextAsync();
        Vm.Status = Vm.StudyPasteCsv(text ?? "");
    }

    private async void OnImport(object? s, RoutedEventArgs e)
    {
        var top = TopLevel.GetTopLevel(this);
        if (top == null) return;
        var files = await top.StorageProvider.OpenFilePickerAsync(new FilePickerOpenOptions
        {
            Title = "Rows from a CSV file", AllowMultiple = true,
            FileTypeFilter = [new FilePickerFileType("CSV or TSV") { Patterns = ["*.csv", "*.tsv", "*.txt"] }, new FilePickerFileType("All files") { Patterns = ["*"] }],
        });
        foreach (var f in files) if (f.TryGetLocalPath() is { } p) Vm.Status = Vm.StudyImportCsv(p);
    }

    private async void OnLoad(object? s, RoutedEventArgs e)
    {
        var top = TopLevel.GetTopLevel(this);
        if (top == null) return;
        var files = await top.StorageProvider.OpenFilePickerAsync(new FilePickerOpenOptions
        {
            Title = "Open a study table", AllowMultiple = false,
            FileTypeFilter = [new FilePickerFileType("CAPS study table") { Patterns = ["*.capstable"] }, new FilePickerFileType("CSV") { Patterns = ["*.csv", "*.tsv"] }],
        });
        if (files.Count > 0 && files[0].TryGetLocalPath() is { } p) Vm.Status = Vm.StudyLoad(p);
    }

    private async Task SaveAs(string title, string ext, string kind)
    {
        var top = TopLevel.GetTopLevel(this);
        if (top == null) return;
        var file = await top.StorageProvider.SaveFilePickerAsync(new FilePickerSaveOptions
        {
            Title = title, SuggestedFileName = "study" + ext, DefaultExtension = ext,
            FileTypeChoices = [new FilePickerFileType(kind) { Patterns = ["*" + ext] }],
        });
        if (file?.TryGetLocalPath() is { } p) Vm.Status = Vm.StudySave(p);
    }
    private async void OnSave(object? s, RoutedEventArgs e) => await SaveAs("Save the study table", ".capstable", "CAPS study table");
    private async void OnExport(object? s, RoutedEventArgs e) => await SaveAs("Export the table as CSV", ".csv", "CSV");
}
