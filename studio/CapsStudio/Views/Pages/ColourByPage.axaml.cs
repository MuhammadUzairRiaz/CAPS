using Avalonia.Controls;
using Avalonia.Interactivity;
using Avalonia.Markup.Xaml;
using CapsStudio.ViewModels;

namespace CapsStudio.Views.Pages;

public partial class ColourByPage : PageBase
{
    public ColourByPage() => AvaloniaXamlLoader.Load(this);

    private void OnTile(object? s, RoutedEventArgs e) { if ((s as Control)?.Tag is ColourTile t) Vm.ApplyColourTile(t); }
    private void OnBack(object? s, RoutedEventArgs e) => Vm.SetModule(20);
}
