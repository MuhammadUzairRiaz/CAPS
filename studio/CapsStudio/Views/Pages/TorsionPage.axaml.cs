using Avalonia.Controls;
using Avalonia.Interactivity;
using Avalonia.Markup.Xaml;
using Avalonia.Platform.Storage;
using CapsStudio.ViewModels;

namespace CapsStudio.Views.Pages;

public partial class TorsionPage : PageBase
{
    public TorsionPage()
    {
        AvaloniaXamlLoader.Load(this);
        var preview = this.FindControl<MolView>("Preview")!;
        var plot = this.FindControl<LinePlot>("Plot")!;
        DataContextChanged += (_, _) =>
        {
            if (DataContext is not MainViewModel vm) return;
            vm.TorsionChanged += () =>
            {
                if (!vm.IsTorsion) return;
                if (!ReferenceEquals(preview.Document, vm.Document)) { preview.Document = vm.Document; preview.Reset(); }
                var at = vm.TorsionAtoms;
                preview.Highlights = at.Where(a => a >= 0).ToArray();
                preview.LineA = at.Length == 4 ? at[1] : -1;
                preview.LineB = at.Length == 4 ? at[2] : -1;
                preview.LineLabel = vm.TorsionPhiText;
                preview.Refresh();
                var curve = vm.TorsionCurve();
                plot.RefY = null;
                plot.SetData(curve, curve);
            };
        };
    }

    private void OnBack(object? s, RoutedEventArgs e) => Vm.SetModule(8);
    private async void OnRun(object? s, RoutedEventArgs e) => await Vm.RunTorsionScan();
    private void OnRow(object? s, RoutedEventArgs e) { if ((s as Control)?.Tag is TorsionRow r) Vm.ShowTorsionRow(r); }

    private async void OnCsv(object? s, RoutedEventArgs e)
    {
        var top = TopLevel.GetTopLevel(this);
        if (top == null) return;
        var file = await top.StorageProvider.SaveFilePickerAsync(new FilePickerSaveOptions
        {
            Title = "Save the torsion scan", SuggestedFileName = "torsion_scan.csv",
            FileTypeChoices = [new FilePickerFileType("CSV") { Patterns = ["*.csv"] }],
        });
        if (file?.TryGetLocalPath() is { } path) { File.WriteAllText(path, Vm.TorsionCsv()); Vm.Status = "Saved " + path; }
    }
}
