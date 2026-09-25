using Avalonia;
using Avalonia.Controls;
using Avalonia.Interactivity;
using Avalonia.Markup.Xaml;
using Avalonia.Platform.Storage;
using CapsStudio.Interop;
using CapsStudio.ViewModels;

namespace CapsStudio.Views.Pages;

public partial class CellEditorPage : PageBase
{
    public CellEditorPage() { AvaloniaXamlLoader.Load(this); }
    private void Show()
    {
        if (DataContext is not MainViewModel vm) return;
        var v = this.FindControl<MolView>("View")!;
        v.ShowCell = true;
        if (v.Document != vm.Document) { v.Document = vm.Document; v.Reset(); } else v.Refresh();
    }
    protected override void OnPropertyChanged(AvaloniaPropertyChangedEventArgs change)
    {
        base.OnPropertyChanged(change);
        if (change.Property == IsVisibleProperty && IsVisible) Show();
    }
    private void OnApply(object? s, RoutedEventArgs e) { Vm.ApplyCell(); Show(); }
    private void OnSupercell(object? s, RoutedEventArgs e) { Vm.MakeSupercell(); Show(); }
    private void OnPreset(object? s, RoutedEventArgs e) { if ((s as Control)?.Tag is CellPresetRow p) Vm.UseCellPreset(p); }
}
