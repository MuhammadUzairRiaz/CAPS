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
}
