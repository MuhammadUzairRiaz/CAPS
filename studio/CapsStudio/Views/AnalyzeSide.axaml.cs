using Avalonia.Controls;
using Avalonia.Interactivity;
using Avalonia.Markup.Xaml;
using CapsStudio.ViewModels;

namespace CapsStudio.Views;

public partial class AnalyzeSide : UserControl
{
    public AnalyzeSide() { AvaloniaXamlLoader.Load(this); }
    private void OnChip(object? s, RoutedEventArgs e) { if (s is Button { Tag: CalcChip c } && DataContext is MainViewModel vm) vm.FocusChip(c); }
    private void OnAll(object? s, RoutedEventArgs e) { if (DataContext is MainViewModel vm) vm.SetModule(1); }
}
