using Avalonia.Controls;
using Avalonia.Input;
using Avalonia.Interactivity;
using Avalonia.Markup.Xaml;
using Avalonia.Platform.Storage;
using CapsStudio.ViewModels;

namespace CapsStudio.Views.Pages;

public partial class BioPage : PageBase
{
    public BioPage()
    {
        AvaloniaXamlLoader.Load(this);
        var preview = this.FindControl<MolView>("Preview")!;
        DataContextChanged += (_, _) =>
        {
            if (DataContext is not MainViewModel vm) return;
            vm.BioViewChanged += () => { var fresh = preview.Document == null; preview.Document = vm.BioDoc; if (fresh) preview.Reset(); };
        };
    }

    private void OnPolymer(object? s, RoutedEventArgs e) => Vm.SetModule(13);
    private void OnCrystal(object? s, RoutedEventArgs e) => Vm.OpenCrystal();
    private void OnSurface(object? s, RoutedEventArgs e) => Vm.OpenSurface();
    private void OnNano(object? s, RoutedEventArgs e) => Vm.OpenNano();
    private void OnSolvation(object? s, RoutedEventArgs e) => Vm.OpenSolvation();
    private void OnCancel(object? s, RoutedEventArgs e) => Vm.SetModule(8);
    private async void OnBuild(object? s, RoutedEventArgs e) => await Vm.BuildPeptide();

    private void OnResidue(object? s, RoutedEventArgs e)
    {
        if ((s as Control)?.Tag is not ResidueCell c) return;
        var shift = (TopLevel.GetTopLevel(this) as Window) is { } && _shift;
        Vm.SelectResidue(c, shift);
    }

    // Shift held while clicking a residue extends the selection
    private bool _shift;
    protected override void OnKeyDown(KeyEventArgs e) { _shift = e.KeyModifiers.HasFlag(KeyModifiers.Shift); base.OnKeyDown(e); }
    protected override void OnKeyUp(KeyEventArgs e) { _shift = e.KeyModifiers.HasFlag(KeyModifiers.Shift); base.OnKeyUp(e); }
    protected override void OnPointerPressed(PointerPressedEventArgs e) { _shift = e.KeyModifiers.HasFlag(KeyModifiers.Shift); base.OnPointerPressed(e); }

    private void OnSs(object? s, RoutedEventArgs e)
    {
        if ((s as Control)?.Tag is string t && int.TryParse(t, out var k)) Vm.BioType = k;
    }

    private async void OnFasta(object? s, RoutedEventArgs e)
    {
        var top = TopLevel.GetTopLevel(this);
        if (top == null) return;
        var files = await top.StorageProvider.OpenFilePickerAsync(new FilePickerOpenOptions
        {
            Title = "Import a sequence (FASTA)",
            AllowMultiple = false,
            FileTypeFilter = [new FilePickerFileType("FASTA") { Patterns = ["*.fasta", "*.fa", "*.faa", "*.txt"] }],
        });
        if (files.Count > 0 && files[0].TryGetLocalPath() is { } path) Vm.ImportFasta(path);
    }
}
