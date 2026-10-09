using Avalonia.Controls;
using Avalonia.Interactivity;
using Avalonia.Markup.Xaml;
using Avalonia.Platform.Storage;
using CapsStudio.ViewModels;

namespace CapsStudio.Views.Pages;

public partial class DftRunsPage : PageBase
{
    private MainViewModel? _hooked;

    public DftRunsPage()
    {
        AvaloniaXamlLoader.Load(this);
        DataContextChanged += (_, _) =>
        {
            if (DataContext is not MainViewModel vm || vm == _hooked) return;
            _hooked = vm;
            vm.DftChanged += () =>
            {
                if (!vm.IsDftRuns) return;
                var p = this.FindControl<LinePlot>("Force")!;
                p.RefY = 0.01;
                p.SetData(vm.SelectedDftRun?.Fmax.Select((f, i) => ((double)(i + 1), f)).ToArray() ?? []);
            };
        };
    }

    private async void OnRefresh(object? s, RoutedEventArgs e) => await Vm.RefreshRuns();
    private async void OnResubmit(object? s, RoutedEventArgs e) => await Vm.ResubmitRun();
    private async void OnReset(object? s, RoutedEventArgs e) => await Vm.ResetRunStage();
    private async void OnDerived(object? s, RoutedEventArgs e) { if ((s as Control)?.Tag is string k) await Vm.MakeDerived(k); }
    private void OnOpen(object? s, RoutedEventArgs e) { Vm.OpenRunStructure(); Vm.SetModule(8); }
    private async void OnCleanup(object? s, RoutedEventArgs e) => await Vm.CleanupRuns(this.FindControl<CheckBox>("KeepBest")!.IsChecked == true);
    private async void OnCopyCli(object? s, RoutedEventArgs e) { if (TopLevel.GetTopLevel(this)?.Clipboard is { } c) await c.SetTextAsync(Vm.DftCli); }
    private async void OnFolder(object? s, RoutedEventArgs e)
    {
        if (TopLevel.GetTopLevel(this) is not { } top) return;
        var d = await top.StorageProvider.OpenFolderPickerAsync(new FolderPickerOpenOptions { Title = "The project folder (structures/, adsorption/ …)" });
        if (d.Count > 0 && d[0].TryGetLocalPath() is { } p) { Vm.RunsRoot = p; await Vm.RefreshRuns(); }
    }
}
