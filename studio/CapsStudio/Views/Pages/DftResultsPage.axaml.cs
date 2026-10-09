using Avalonia.Controls;
using Avalonia.Interactivity;
using Avalonia.Markup.Xaml;
using Avalonia.Platform.Storage;
using CapsStudio.ViewModels;

namespace CapsStudio.Views.Pages;

public partial class DftResultsPage : PageBase
{
    private MainViewModel? _hooked;

    public DftResultsPage()
    {
        AvaloniaXamlLoader.Load(this);
        DataContextChanged += (_, _) =>
        {
            if (DataContext is not MainViewModel vm || vm == _hooked) return;
            _hooked = vm;
            vm.DftChanged += () => { if (vm.IsDftResults) Draw(vm); };
        };
    }

    private void Draw(MainViewModel vm)
    {
        var b = this.FindControl<LinePlot>("Bind")!;
        b.RefY = 0;
        b.SetData(vm.BindBars);
        var d = this.FindControl<LinePlot>("Dos")!;
        d.RefY = null;
        d.CursorX = 0;
        d.SetData(vm.DosGroups.TryGetValue(vm.DosGroup, out var g) ? g : []);
        var w = this.FindControl<LinePlot>("Wf")!;
        w.RefY = 0;
        w.SetData(vm.WfProfile);
        var c = this.FindControl<LinePlot>("Cdd")!;
        c.RefY = 0;
        c.SetCompare(vm.CddProfile, vm.CddCumulative);
        var m = this.FindControl<LinePlot>("Md")!;
        m.RefY = null;
        m.SetData(vm.MdHeight);
    }

    private async void OnRefresh(object? s, RoutedEventArgs e) => await Vm.RefreshResults();
    private async Task<string?> Folder(string title)
    {
        if (TopLevel.GetTopLevel(this) is not { } top) return null;
        var d = await top.StorageProvider.OpenFolderPickerAsync(new FolderPickerOpenOptions { Title = title });
        return d.Count > 0 ? d[0].TryGetLocalPath() : null;
    }
    private async void OnSet(object? s, RoutedEventArgs e) { if (await Folder("The adsorption set (slab/, molecule/, complex_*)") is { } p) Vm.ResultsSet = p; }
    private async void OnAimd(object? s, RoutedEventArgs e) { if (await Folder("The AIMD folder (seg_001 …)") is { } p) Vm.ResultsAimd = p; }
    private async void OnMain(object? s, RoutedEventArgs e) { if (await Folder("The project (structures/, adsorption/, aimd/)") is { } p) Vm.ResultsMain = p; }
}
