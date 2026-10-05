using Avalonia.Controls;
using Avalonia.Interactivity;
using Avalonia.Markup.Xaml;
using Avalonia.Platform.Storage;
using CapsStudio.ViewModels;

namespace CapsStudio.Views.Pages;

public partial class CrystalPage : PageBase
{
    public CrystalPage()
    {
        AvaloniaXamlLoader.Load(this);
        var preview = this.FindControl<MolView>("Preview")!;
        DataContextChanged += (_, _) =>
        {
            if (DataContext is not MainViewModel vm) return;
            vm.CrystalViewChanged += () => { var fresh = preview.Document == null; preview.Document = vm.CrystalDoc; if (fresh) preview.Reset(); };
        };
    }

    private void OnPolymer(object? s, RoutedEventArgs e) => Vm.SetModule(13);
    private void OnSurface(object? s, RoutedEventArgs e) => Vm.OpenSurface();
    private void OnNano(object? s, RoutedEventArgs e) => Vm.OpenNano();
    private void OnSolvation(object? s, RoutedEventArgs e) => Vm.OpenSolvation();
    private void OnBio(object? s, RoutedEventArgs e) => Vm.OpenBio();
    private void OnCancel(object? s, RoutedEventArgs e) => Vm.SetModule(8);
    private async void OnBuild(object? s, RoutedEventArgs e) => await Vm.BuildCrystal();
    private void OnAddSite(object? s, RoutedEventArgs e) => Vm.AddCrystalSite();
    private void OnSitesFromMolecule(object? s, RoutedEventArgs e) => Vm.SitesFromOpenMolecule();
    private void OnRemoveSite(object? s, RoutedEventArgs e) { if ((s as Control)?.Tag is CrystalSiteRow r) Vm.RemoveCrystalSite(r); }
    private void OnApplySymmetry(object? s, RoutedEventArgs e) => Vm.CrystalApplySymmetry();
    private async void OnFindSymmetry(object? s, RoutedEventArgs e) => await Vm.CrystalFindSymmetry();

    private async void OnImportCif(object? s, RoutedEventArgs e)
    {
        var top = TopLevel.GetTopLevel(this);
        if (top == null) return;
        var files = await top.StorageProvider.OpenFilePickerAsync(new FilePickerOpenOptions
        {
            Title = "Import a crystal (CIF)",
            AllowMultiple = false,
            FileTypeFilter = [new FilePickerFileType("Crystallographic Information File") { Patterns = ["*.cif"] }],
        });
        if (files.Count > 0 && files[0].TryGetLocalPath() is { } path) await Vm.ImportCrystalCif(path);
    }

    private async void OnLibrary(object? s, SelectionChangedEventArgs e)
    {
        if (s is not ComboBox box || box.SelectedItem is not CrystalEntry c) return;
        await Vm.ImportCrystalCif(c.File);
        box.SelectedItem = null;
    }
    private void OnCg(object? s, Avalonia.Interactivity.RoutedEventArgs e) => Vm.OpenCg();
    private void OnFindPrimitive(object? s, RoutedEventArgs e) => Vm.FindPrimitiveCell();
    private void OnNiggli(object? s, RoutedEventArgs e) => Vm.NiggliCell();
    private void OnConventional(object? s, RoutedEventArgs e) => Vm.ConventionalCell();
    private void OnCellEditor(object? s, RoutedEventArgs e) => Vm.OpenCellEditor();
    private void OnRedefine(object? s, RoutedEventArgs e) => Vm.RedefineLattice();
    private void OnVacuumSlab(object? s, RoutedEventArgs e) => Vm.MakeVacuumSlab();
    private void OnCluster(object? s, RoutedEventArgs e) => Vm.CutCluster();
    private void OnDefects(object? s, RoutedEventArgs e) => Vm.MakeDefects();
    private void OnNanowire(object? s, RoutedEventArgs e) => Vm.MakeNanowire();
}
