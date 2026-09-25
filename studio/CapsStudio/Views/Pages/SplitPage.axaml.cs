using Avalonia.Controls;
using Avalonia.Interactivity;
using Avalonia.Markup.Xaml;
using Avalonia.Platform.Storage;
using CapsStudio.ViewModels;

namespace CapsStudio.Views.Pages;

public partial class SplitPage : PageBase
{
    public SplitPage()
    {
        AvaloniaXamlLoader.Load(this);
        var a = this.FindControl<MolView>("PaneA")!;
        var b = this.FindControl<MolView>("PaneB")!;
        DataContextChanged += (_, _) =>
        {
            if (DataContext is not MainViewModel vm) return;
            vm.SplitChanged += () =>
            {
                if (!vm.IsSplit) return;
                if (!ReferenceEquals(a.Document, vm.Document)) { a.Document = vm.Document; a.Reset(); }
                if (!ReferenceEquals(b.Document, vm.SplitDocB)) { b.Document = vm.SplitDocB; if (vm.SyncCamera) b.Camera = a.Camera; else b.Reset(); }
                var picked = vm.PickedAtoms;
                a.Highlights = picked;
                b.Highlights = vm.SyncSelection ? picked : [];
                a.Refresh();
                b.Refresh();
            };
            // one camera: a drag in either pane moves both
            a.CameraChanged += cam => { if (vm.SyncCamera) b.Camera = cam; };
            b.CameraChanged += cam => { if (vm.SyncCamera) a.Camera = cam; };
            a.AtomClicked += i => { vm.Pick(i); vm.RefreshSplitSelection(); };
            b.AtomClicked += i => { if (vm.SyncSelection) { vm.Pick(i); vm.RefreshSplitSelection(); } };
        };
    }

    private void OnBack(object? s, RoutedEventArgs e) => Vm.SetModule(8);
    private void OnTheme(object? s, RoutedEventArgs e) => Vm.ToggleSplitTheme();
    private void OnCloseB(object? s, RoutedEventArgs e) => Vm.CloseSplitB();
    private async void OnRecentB(object? s, SelectionChangedEventArgs e)
    {
        if (s is not ComboBox box || box.SelectedItem is not RecentItem r) return;
        await Vm.SetSplitB(r.Path, r.Topology);
        box.SelectedItem = null;
    }
    private async void OnOpenB(object? s, RoutedEventArgs e)
    {
        var top = TopLevel.GetTopLevel(this);
        if (top == null) return;
        var files = await top.StorageProvider.OpenFilePickerAsync(new FilePickerOpenOptions { Title = "A structure for the right pane", AllowMultiple = false });
        if (files.Count > 0 && files[0].TryGetLocalPath() is { } path) await Vm.SetSplitB(path);
    }
}
