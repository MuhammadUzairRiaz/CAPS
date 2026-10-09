using Avalonia;
using Avalonia.Controls;
using Avalonia.Interactivity;
using Avalonia.Markup.Xaml;
using Avalonia.Platform.Storage;
using CapsStudio.ViewModels;

namespace CapsStudio.Views.Pages;

public partial class AdsorbDftPage : PageBase
{
    private MainViewModel? _hooked;

    public AdsorbDftPage()
    {
        AvaloniaXamlLoader.Load(this);
        DataContextChanged += (_, _) =>
        {
            if (DataContext is not MainViewModel vm || vm == _hooked) return;
            _hooked = vm;
            vm.DftChanged += () => { if (vm.IsAdsorbDft) Show(vm); };
        };
    }

    private void Show(MainViewModel vm)
    {
        var view = this.FindControl<MolView>("View")!;
        view.ShowCell = true;
        if (view.Document != vm.AdsPreviewDoc) { view.Document = vm.AdsPreviewDoc; view.Reset(); } else view.Refresh();
    }

    protected override void OnPropertyChanged(AvaloniaPropertyChangedEventArgs change)
    {
        base.OnPropertyChanged(change);
        if (change.Property == IsVisibleProperty && IsVisible && DataContext is MainViewModel vm) Show(vm);
    }

    private async void OnBuild(object? s, RoutedEventArgs e) => await Vm.BuildAdsorptionSet();
    private void OnComplex(object? s, RoutedEventArgs e) { if ((s as Control)?.Tag is DftComplexRow c) Vm.PreviewComplex(c); }
    private void OnToJobs(object? s, RoutedEventArgs e) => Vm.OpenDftJob();
    private void OnMolClear(object? s, RoutedEventArgs e) => Vm.AdsMolFile = "";
    private async void OnCopyCli(object? s, RoutedEventArgs e) { if (TopLevel.GetTopLevel(this)?.Clipboard is { } c) await c.SetTextAsync(Vm.DftCli); }
    private async Task<string?> Pick(string title)
    {
        if (TopLevel.GetTopLevel(this) is not { } top) return null;
        var f = await top.StorageProvider.OpenFilePickerAsync(new FilePickerOpenOptions { Title = title, AllowMultiple = false });
        return f.Count > 0 ? f[0].TryGetLocalPath() : null;
    }
    private async void OnSlab(object? s, RoutedEventArgs e) { if (await Pick("The slab (POSCAR, or a relaxed CONTCAR)") is { } p) Vm.AdsSlab = p; }
    private async void OnMolFile(object? s, RoutedEventArgs e) { if (await Pick("The molecule (any structure file)") is { } p) Vm.AdsMolFile = p; }
}
