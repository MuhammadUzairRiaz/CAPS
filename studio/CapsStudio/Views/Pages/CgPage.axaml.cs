using Avalonia.Controls;
using Avalonia.Interactivity;
using Avalonia.Markup.Xaml;
using Avalonia.Platform.Storage;
using CapsStudio.ViewModels;

namespace CapsStudio.Views.Pages;

public partial class CgPage : PageBase
{
    public CgPage()
    {
        AvaloniaXamlLoader.Load(this);
        var preview = this.FindControl<MolView>("Preview")!;
        DataContextChanged += (_, _) =>
        {
            if (DataContext is not MainViewModel vm) return;
            vm.CgViewChanged += () => { var fresh = preview.Document == null; preview.Document = vm.CgDoc; if (fresh) preview.Reset(); };
        };
    }

    private void OnPolymer(object? s, RoutedEventArgs e) => Vm.SetModule(13);
    private void OnModelResolution(object? s, RoutedEventArgs e) => Vm.OpenModelResolution();
    private async void OnCgBackmap(object? s, RoutedEventArgs e) => await Vm.BackmapCg();
    private void OnCrystal(object? s, RoutedEventArgs e) => Vm.OpenCrystal();
    private void OnSurface(object? s, RoutedEventArgs e) => Vm.OpenSurface();
    private void OnNano(object? s, RoutedEventArgs e) => Vm.OpenNano();
    private void OnBio(object? s, RoutedEventArgs e) => Vm.OpenBio();
    private void OnSolvation(object? s, RoutedEventArgs e) => Vm.OpenSolvation();
    private async void OnBuild(object? s, RoutedEventArgs e) { if (Vm.CgIsMartini) await Vm.BuildMartini(); else Vm.BuildCg(); }

    private async void OnExport(object? s, RoutedEventArgs e)
    {
        var top = TopLevel.GetTopLevel(this);
        if (top == null) return;
        var file = await top.StorageProvider.SaveFilePickerAsync(new FilePickerSaveOptions
        {
            Title = "LAMMPS deck (writes STEM.data and STEM.in)", SuggestedFileName = $"KG_melt_{Vm.CgChains:0}x{Vm.CgBeads:0}", DefaultExtension = "in",
            FileTypeChoices = [new FilePickerFileType("LAMMPS input") { Patterns = ["*.in"] }],
        });
        if (file?.TryGetLocalPath() is not { } p) return;
        var stem = Path.Combine(Path.GetDirectoryName(p)!, Path.GetFileNameWithoutExtension(p));
        try { Vm.Status = Vm.ExportCg(stem); } catch (Exception ex) { Vm.Status = "Could not write: " + ex.Message; }
    }
}
