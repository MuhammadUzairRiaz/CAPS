using Avalonia.Controls;
using Avalonia.VisualTree;
using CapsStudio.ViewModels;

namespace CapsStudio.Views.Pages;

/// <summary>A module page (one design board): reads the main view model; actions that need the window (file
/// dialogs, the 3D view) go through it.</summary>
public abstract class PageBase : UserControl
{
    protected MainViewModel Vm => (MainViewModel)DataContext!;
    protected MainWindow? Window => this.FindAncestorOfType<MainWindow>();

    /// <summary>Copy as command: the page's run as the caps command line runs it (the same engine and options).</summary>
    protected async void CopyCommand(Func<string> command, string what)
    {
        var clip = TopLevel.GetTopLevel(this)?.Clipboard;
        if (clip != null) await clip.SetTextAsync(command());
        Vm.Status = $"Copied the {what} as a caps command (the same options and defaults as this page)";
    }
}
