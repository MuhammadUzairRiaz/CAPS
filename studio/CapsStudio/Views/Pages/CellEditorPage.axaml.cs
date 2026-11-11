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
    private MainViewModel? _hooked;
    public CellEditorPage()
    {
        AvaloniaXamlLoader.Load(this);
        // the page draws the structure in its own view: every change the main view is told about (edits, undo,
        // symmetrise, space group, lattice planes) redraws it too while the page is open
        DataContextChanged += (_, _) =>
        {
            if (_hooked != null) _hooked.RenderRequested -= OnRenderRequested;
            _hooked = DataContext as MainViewModel;
            if (_hooked != null) _hooked.RenderRequested += OnRenderRequested;
        };
    }
    private void OnRenderRequested() { if (IsVisible) Show(); }
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
    private void OnFindPointGroup(object? s, Avalonia.Interactivity.RoutedEventArgs e) => Vm.FindPointGroup();
    private void OnSymmetrizeMolecule(object? s, Avalonia.Interactivity.RoutedEventArgs e) => Vm.SymmetrizeMolecule();
    private void OnShowMiller(object? s, Avalonia.Interactivity.RoutedEventArgs e) => Vm.ShowMiller();
    private void OnViewMiller(object? s, Avalonia.Interactivity.RoutedEventArgs e)
    {
        Vm.ViewOntoMiller();
        if (!Vm.MillerShown || Vm.MillerHasError) return;
        // turn the page's own view (the main view is hidden behind the page) straight down the planes' normal
        var v = this.FindControl<MolView>("View")!;
        var (yaw, pitch) = Vm.MillerView;
        var c = v.Camera;
        v.Camera = new CapsCamera { Yaw = yaw, Pitch = pitch, Zoom = c.Zoom, Perspective = c.Perspective, FovDeg = c.FovDeg };
    }
    private void OnHideMiller(object? s, Avalonia.Interactivity.RoutedEventArgs e) => Vm.HideMiller();
    private void OnFindGroup(object? s, Avalonia.Interactivity.RoutedEventArgs e) => Vm.FindSpaceGroup();
    private void OnDropGroup(object? s, Avalonia.Interactivity.RoutedEventArgs e) => Vm.DropSpaceGroup();
    private void OnSupercell(object? s, RoutedEventArgs e) { Vm.MakeSupercell(); Show(); }
    private void OnPreset(object? s, RoutedEventArgs e) { if ((s as Control)?.Tag is CellPresetRow p) Vm.UseCellPreset(p); }
}
