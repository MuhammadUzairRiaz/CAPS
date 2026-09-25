using Avalonia;
using Avalonia.Controls;
using Avalonia.Interactivity;
using Avalonia.Markup.Xaml;
using Avalonia.Platform.Storage;
using CapsStudio.ViewModels;

namespace CapsStudio.Views.Pages;

public partial class ChargesPage : PageBase
{
    private MainViewModel? _hooked;

    public ChargesPage()
    {
        AvaloniaXamlLoader.Load(this);
        DataContextChanged += (_, _) =>
        {
            if (DataContext is not MainViewModel vm || vm == _hooked) return;
            _hooked = vm;
            vm.ChargesChanged += () =>
            {
                ShowView(vm);
                var h = this.FindControl<LinePlot>("Hist")!;
                h.RefY = null;
                h.SetData(vm.ChargeHistogram);
            };
        };
    }

    private void ShowView(MainViewModel vm)
    {
        var view = this.FindControl<MolView>("View")!;
        if (view.Document != vm.Document) { view.Document = vm.Document; view.Reset(); }
        else view.Refresh();
    }

    protected override void OnPropertyChanged(AvaloniaPropertyChangedEventArgs change)
    {
        base.OnPropertyChanged(change);
        if (change.Property == IsVisibleProperty && IsVisible && DataContext is MainViewModel vm) ShowView(vm);
    }

    private void OnApply(object? s, RoutedEventArgs e) => Vm.ApplyCharges();

    private async void OnImport(object? s, RoutedEventArgs e)
    {
        var top = TopLevel.GetTopLevel(this);
        if (top == null) return;
        var files = await top.StorageProvider.OpenFilePickerAsync(new FilePickerOpenOptions
        {
            Title = "Charges from RESP, AM1-BCC or another program",
            FileTypeFilter = [new FilePickerFileType("Charges") { Patterns = ["*.chg", "*.txt", "*.dat"] }, new FilePickerFileType("All files") { Patterns = ["*"] }],
        });
        if (files.Count > 0 && files[0].TryGetLocalPath() is { } p) Vm.ChooseChargeFile(p);
    }

    private async void OnExport(object? s, RoutedEventArgs e)
    {
        var top = TopLevel.GetTopLevel(this);
        if (top == null) return;
        var file = await top.StorageProvider.SaveFilePickerAsync(new FilePickerSaveOptions
        {
            Title = "Export with these charges", SuggestedFileName = "charged.mol2", DefaultExtension = "mol2",
            FileTypeChoices = [new FilePickerFileType("Tripos mol2") { Patterns = ["*.mol2"] }],
        });
        if (file?.TryGetLocalPath() is { } p) Vm.ExportChargesMol2(p);
    }
}
