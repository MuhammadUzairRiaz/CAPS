using Avalonia.Controls;
using Avalonia.Interactivity;
using Avalonia.Markup.Xaml;
using Avalonia.Platform.Storage;
using CapsStudio.ViewModels;

namespace CapsStudio.Views.Pages;

public partial class GlassPage : PageBase
{
    private MainViewModel? _hooked;

    public GlassPage()
    {
        AvaloniaXamlLoader.Load(this);
        DataContextChanged += (_, _) =>
        {
            if (DataContext is not MainViewModel vm || vm == _hooked) return;
            _hooked = vm;
            vm.GtChanged += () =>
            {
                var sch = this.FindControl<LinePlot>("Schedule")!;
                sch.RefY = null;
                sch.XLabel = vm.GtScheduleAxis;
                sch.SetData(vm.GtSchedule);
                var vt = this.FindControl<LinePlot>("Vt")!;
                vt.YLabel = vm.GtYLabel;
                vt.RefY = null;
                vt.Errors = vm.GtErrors.Length == vm.GtPoints.Length ? vm.GtErrors : null;
                vt.SetData(vm.GtPoints, vm.GtFit);
            };
        };
    }

    private void OnPlan(object? s, NumericUpDownValueChangedEventArgs e) => Vm.RaiseGlassPlan();
    private async void OnRun(object? s, RoutedEventArgs e) => await Vm.RunGlass();
    private void OnCancel(object? s, RoutedEventArgs e) => Vm.Analyze.Cancel();

    private async void OnRecipe(object? s, RoutedEventArgs e)
    {
        var top = TopLevel.GetTopLevel(this);
        var src = Vm.Document?.Path;
        if (top == null || string.IsNullOrEmpty(src)) { Vm.Status = "Save the structure first: the recipe reads it from its file"; return; }
        var file = await top.StorageProvider.SaveFilePickerAsync(new FilePickerSaveOptions
        {
            Title = "Save the cooling scan as a recipe", SuggestedFileName = System.IO.Path.GetFileNameWithoutExtension(src) + "_tg.yaml", DefaultExtension = "yaml",
            FileTypeChoices = [new FilePickerFileType("CAPS recipe") { Patterns = ["*.yaml", "*.yml"] }],
        });
        if (file?.TryGetLocalPath() is not { } p) return;
        System.IO.File.WriteAllText(p, Vm.GlassRecipe(src));
        Vm.Status = $"Saved {System.IO.Path.GetFileName(p)} · caps run {System.IO.Path.GetFileName(p)}";
    }
}
