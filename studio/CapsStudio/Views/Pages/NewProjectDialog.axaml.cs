using Avalonia.Controls;
using Avalonia.Input;
using Avalonia.Interactivity;
using Avalonia.Markup.Xaml;
using Avalonia.Platform.Storage;

namespace CapsStudio.Views.Pages;

public partial class NewProjectDialog : PageBase
{
    public NewProjectDialog()
    {
        AvaloniaXamlLoader.Load(this);
        // the name is ready to type over when the sheet opens
        PropertyChanged += (_, e) =>
        {
            if (e.Property == IsVisibleProperty && IsVisible && this.FindControl<TextBox>("NameBox") is { } b)
                Avalonia.Threading.Dispatcher.UIThread.Post(() => { b.Focus(); b.SelectAll(); }, Avalonia.Threading.DispatcherPriority.Input);
        };
    }

    private void OnBackdrop(object? s, PointerPressedEventArgs e) => Vm.NewProjectOpen = false;
    private void OnCancel(object? s, RoutedEventArgs e) => Vm.NewProjectOpen = false;
    private void OnCreate(object? s, RoutedEventArgs e) => Vm.CreateNewProject();
    private void OnNameKey(object? s, KeyEventArgs e)
    {
        if (e.Key == Key.Enter) { e.Handled = true; Vm.CreateNewProject(); }
        else if (e.Key == Key.Escape) { e.Handled = true; Vm.NewProjectOpen = false; }
    }

    private async void OnChoose(object? s, RoutedEventArgs e)
    {
        if (TopLevel.GetTopLevel(this)?.StorageProvider is not { } sp) return;
        var start = Vm.NewProjectParent.Length > 0 ? await sp.TryGetFolderFromPathAsync(Vm.NewProjectParent) : null;
        var picked = await sp.OpenFolderPickerAsync(new FolderPickerOpenOptions { Title = "Where the project folder goes", AllowMultiple = false, SuggestedStartLocation = start });
        if (picked.Count > 0 && picked[0].TryGetLocalPath() is { } p) Vm.NewProjectParent = p;
    }
}
