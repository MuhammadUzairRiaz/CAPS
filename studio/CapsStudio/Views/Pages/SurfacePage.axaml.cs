using Avalonia.Controls;
using Avalonia.Interactivity;
using Avalonia.Markup.Xaml;
using Avalonia.Platform.Storage;
using CapsStudio.ViewModels;

namespace CapsStudio.Views.Pages;

public partial class SurfacePage : PageBase
{
    public SurfacePage()
    {
        AvaloniaXamlLoader.Load(this);
        var preview = this.FindControl<MolView>("Preview")!;
        DataContextChanged += (_, _) =>
        {
            if (DataContext is not MainViewModel vm) return;
            vm.SurfViewChanged += () => { preview.Document = vm.SurfDoc; };
            // a new termination list resets the box's selection: put the model's back once the items are in
            var terms = this.FindControl<ComboBox>("Terms")!;
            vm.PropertyChanged += (_, e) =>
            {
                if (e.PropertyName == nameof(MainViewModel.SurfTerminations))
                    Avalonia.Threading.Dispatcher.UIThread.Post(() => terms.SelectedIndex = vm.SurfTermination, Avalonia.Threading.DispatcherPriority.Background);
            };
        };
    }

    private void OnPolymer(object? s, RoutedEventArgs e) => Vm.SetModule(13);
    private void OnSolvation(object? s, RoutedEventArgs e) => Vm.OpenSolvation();
    private void OnNano(object? s, RoutedEventArgs e) => Vm.OpenNano();
    private void OnCrystal(object? s, RoutedEventArgs e) => Vm.OpenCrystal();
    private void OnBio(object? s, RoutedEventArgs e) => Vm.OpenBio();
    private void OnCancel(object? s, RoutedEventArgs e) => Vm.SetModule(8);
    private async void OnBuild(object? s, RoutedEventArgs e) => await Vm.BuildSurface();

    private async void OnOpenCif(object? s, RoutedEventArgs e)
    {
        var top = TopLevel.GetTopLevel(this);
        if (top == null) return;
        var files = await top.StorageProvider.OpenFilePickerAsync(new FilePickerOpenOptions
        {
            Title = "Open a crystal (CIF)",
            AllowMultiple = false,
            FileTypeFilter = [new FilePickerFileType("Crystallographic Information File") { Patterns = ["*.cif"] }],
        });
        if (files.Count > 0 && files[0].TryGetLocalPath() is { } path) Vm.UseCif(path);
    }
    private void OnCg(object? s, Avalonia.Interactivity.RoutedEventArgs e) => Vm.OpenCg();
    private void OnRemoveLayer(object? s, RoutedEventArgs e) { if ((s as Control)?.Tag is int k) Vm.RemoveSurfLayer(k); }
    private void OnAddSlabLayer(object? s, RoutedEventArgs e) => Vm.AddSlabLayer();
    private void OnLayerFromGrow(object? s, RoutedEventArgs e) => Vm.AddLayerFromOpen();
    private async void OnAddLayerFile(object? s, RoutedEventArgs e)
    {
        var top = TopLevel.GetTopLevel(this);
        if (top == null) return;
        var files = await top.StorageProvider.OpenFilePickerAsync(new FilePickerOpenOptions { Title = "A layer: a structure with a rectangular periodic cell", AllowMultiple = false });
        if (files.Count > 0 && files[0].TryGetLocalPath() is { } path) Vm.AddLayerFile(path);
    }
}
