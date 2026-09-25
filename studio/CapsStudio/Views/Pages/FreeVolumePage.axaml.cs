using Avalonia;
using Avalonia.Controls;
using Avalonia.Interactivity;
using Avalonia.Markup.Xaml;
using Avalonia.Platform.Storage;
using CapsStudio.ViewModels;

namespace CapsStudio.Views.Pages;

public partial class FreeVolumePage : PageBase
{
    private MainViewModel? _hooked;

    public FreeVolumePage()
    {
        AvaloniaXamlLoader.Load(this);
        DataContextChanged += (_, _) =>
        {
            if (DataContext is not MainViewModel vm || vm == _hooked) return;
            _hooked = vm;
            vm.FvChanged += () =>
            {
                var view = this.FindControl<MolView>("View")!;
                view.DrawStyle = 2;       // sticks, so the voids read
                view.ShowCell = true;
                if (view.Document != vm.Document) { view.Document = vm.Document; view.Reset(); }
                else view.Document = vm.Document;
                var psd = this.FindControl<LinePlot>("Psd")!;
                psd.RefY = null;
                psd.SetData(vm.FvPsd);
            };
        };
    }

    protected override void OnPropertyChanged(AvaloniaPropertyChangedEventArgs change)
    {
        base.OnPropertyChanged(change);
        if (change.Property == IsVisibleProperty && IsVisible && DataContext is MainViewModel vm)
        {
            var view = this.FindControl<MolView>("View")!;
            if (view.Document != vm.Document) { view.DrawStyle = 2; view.ShowCell = true; view.Document = vm.Document; view.Reset(); }
        }
    }

    private async void OnRun(object? s, RoutedEventArgs e) => await Vm.RunFreeVolume();
    private void OnCancel(object? s, RoutedEventArgs e) => Vm.Analyze.Cancel();

    private async void OnExport(object? s, RoutedEventArgs e)
    {
        var top = TopLevel.GetTopLevel(this);
        if (top == null) return;
        var file = await top.StorageProvider.SaveFilePickerAsync(new FilePickerSaveOptions
        {
            Title = "Export voids", SuggestedFileName = "voids.pdb", DefaultExtension = "pdb",
            FileTypeChoices = [new FilePickerFileType("PDB") { Patterns = ["*.pdb"] }],
        });
        if (file?.TryGetLocalPath() is { } p) Vm.ExportVoids(p);
    }
}
