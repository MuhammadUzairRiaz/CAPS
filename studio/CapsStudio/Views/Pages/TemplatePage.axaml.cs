using Avalonia.Controls;
using Avalonia.Interactivity;
using Avalonia.Markup.Xaml;
using Avalonia.Platform.Storage;
using CapsStudio.ViewModels;

namespace CapsStudio.Views.Pages;

public partial class TemplatePage : PageBase
{
    public TemplatePage()
    {
        AvaloniaXamlLoader.Load(this);
        DataContextChanged += (_, _) =>
        {
            if (DataContext is not MainViewModel vm) return;
            // pre-reaction pane: a drag from one atom to another forms that bond (or breaks it when the pattern has it)
            this.FindControl<TemplateDrawing>("Pre")!.Linked += (a, b) => vm.ToggleTemplateBond(a, b);
            var pre = this.FindControl<TemplateDrawing>("Pre")!;
            var post = this.FindControl<TemplateDrawing>("Post")!;
            pre.Hovered += m => post.Lit = m;
            post.Hovered += m => pre.Lit = m;
            vm.TemplateDrawn += () =>
            {
                this.FindControl<TemplateDrawing>("Pre")!.Set(vm.TemplatePre, [], vm.TemplateBroken);
                this.FindControl<TemplateDrawing>("Post")!.Set(vm.TemplatePost, vm.TemplateFormed, []);
            };
        };
    }

    private void OnTest(object? s, RoutedEventArgs e) => Vm.TestTemplate();
    private void OnUse(object? s, RoutedEventArgs e) => Vm.UseTemplateInReact();
    private void OnSave(object? s, RoutedEventArgs e) { try { Vm.Status = Vm.SaveTemplate(); } catch (Exception ex) { Vm.Status = ex.Message; } }

    private async void OnImport(object? s, RoutedEventArgs e)
    {
        var top = TopLevel.GetTopLevel(this);
        if (top == null) return;
        var files = await top.StorageProvider.OpenFilePickerAsync(new FilePickerOpenOptions { Title = "Reaction template (CAPS template text)", AllowMultiple = false });
        if (files.Count > 0 && files[0].TryGetLocalPath() is { } p) Vm.ImportTemplate(p);
    }
}
