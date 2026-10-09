using Avalonia;
using Avalonia.Controls;
using Avalonia.Interactivity;
using Avalonia.Markup.Xaml;
using Avalonia.Platform.Storage;
using CapsStudio.ViewModels;

namespace CapsStudio.Views.Pages;

public partial class SheetsPage : PageBase
{
    private MainViewModel? _hooked;

    public SheetsPage()
    {
        AvaloniaXamlLoader.Load(this);
        DataContextChanged += (_, _) =>
        {
            if (DataContext is not MainViewModel vm || vm == _hooked) return;
            _hooked = vm;
            vm.DftChanged += () => { if (vm.IsSheets) Show(vm); };
        };
    }

    private void Show(MainViewModel vm)
    {
        var view = this.FindControl<MolView>("View")!;
        view.ShowCell = true;
        view.Highlights = [];
        if (view.Document != vm.SlabDoc) { view.Document = vm.SlabDoc; view.Reset(); } else view.Refresh();
        this.FindControl<SiteMap>("Sites")!.Set(vm.SheetSites, vm.SheetCell);
    }

    protected override void OnPropertyChanged(AvaloniaPropertyChangedEventArgs change)
    {
        base.OnPropertyChanged(change);
        if (change.Property == IsVisibleProperty && IsVisible && DataContext is MainViewModel vm) Show(vm);
    }

    private async void OnBuild(object? s, RoutedEventArgs e) => await Vm.BuildSlab();
    private void OnAddTop(object? s, RoutedEventArgs e) => Vm.AddTerm(true);
    private void OnAddBottom(object? s, RoutedEventArgs e) => Vm.AddTerm(false);
    private void OnRemoveTerm(object? s, RoutedEventArgs e) { if ((s as Control)?.DataContext is TermRow r) Vm.RemoveTerm(r); }
    private void OnFromClear(object? s, RoutedEventArgs e) => Vm.SheetFrom = "";
    private void OnFinding(object? s, RoutedEventArgs e)
    {
        if ((s as Control)?.Tag is not DftFinding f) return;
        var view = this.FindControl<MolView>("View")!;
        view.Highlights = f.Atoms;
        view.Refresh();
    }
    private async void OnCopyCli(object? s, RoutedEventArgs e) { if (TopLevel.GetTopLevel(this)?.Clipboard is { } c) await c.SetTextAsync(Vm.DftCli); }
    private void OnOpenInStudio(object? s, RoutedEventArgs e) { if (Vm.HasSlab) { Vm.Open(Vm.SlabPath); Vm.SetModule(8); } }

    private async void OnFrom(object? s, RoutedEventArgs e)
    {
        if (TopLevel.GetTopLevel(this) is not { } top) return;
        var f = await top.StorageProvider.OpenFilePickerAsync(new FilePickerOpenOptions { Title = "A structure to cut one layer out of", AllowMultiple = false });
        if (f.Count > 0 && f[0].TryGetLocalPath() is { } p) Vm.SheetFrom = p;
    }
    private async Task<string?> SaveAs(string title, string name, string ext)
    {
        if (TopLevel.GetTopLevel(this) is not { } top) return null;
        var f = await top.StorageProvider.SaveFilePickerAsync(new FilePickerSaveOptions { Title = title, SuggestedFileName = name, DefaultExtension = ext });
        return f?.TryGetLocalPath();
    }
    private async void OnExportPoscar(object? s, RoutedEventArgs e) { if (await SaveAs("Slab as POSCAR", System.IO.Path.GetFileNameWithoutExtension(Vm.SlabPath) + ".vasp", "vasp") is { } p) await Vm.ExportSlab(p, false); }
    private async void OnExportCif(object? s, RoutedEventArgs e) { if (await SaveAs("Slab as CIF", System.IO.Path.GetFileNameWithoutExtension(Vm.SlabPath) + ".cif", "cif") is { } p) await Vm.ExportSlab(p, false); }
    private async void OnExportSet(object? s, RoutedEventArgs e)
    {
        if (TopLevel.GetTopLevel(this) is not { } top) return;
        var d = await top.StorageProvider.OpenFolderPickerAsync(new FolderPickerOpenOptions { Title = "VASP set: the case folder (created inside the one you pick)" });
        if (d.Count > 0 && d[0].TryGetLocalPath() is { } p) await Vm.ExportSlab(System.IO.Path.Combine(p, System.IO.Path.GetFileNameWithoutExtension(Vm.SlabPath)), true);
    }
}
