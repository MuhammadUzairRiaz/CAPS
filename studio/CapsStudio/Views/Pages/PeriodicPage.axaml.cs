using Avalonia;
using Avalonia.Controls;
using Avalonia.Interactivity;
using Avalonia.Markup.Xaml;
using CapsStudio.ViewModels;

namespace CapsStudio.Views.Pages;

public partial class PeriodicPage : PageBase
{
    private MainViewModel? _hooked;

    public PeriodicPage()
    {
        AvaloniaXamlLoader.Load(this);
        DataContextChanged += (_, _) =>
        {
            if (DataContext is not MainViewModel vm || vm == _hooked) return;
            _hooked = vm;
            vm.PeriodicChanged += () => ShowView(vm);
        };
    }

    private void ShowView(MainViewModel vm)
    {
        var view = this.FindControl<MolView>("View")!;
        view.ShowCell = true;
        if (view.Document != vm.Document) { view.Document = vm.Document; view.Reset(); }
        else view.Refresh();
    }

    protected override void OnPropertyChanged(AvaloniaPropertyChangedEventArgs change)
    {
        base.OnPropertyChanged(change);
        if (change.Property == IsVisibleProperty && IsVisible && DataContext is MainViewModel vm) ShowView(vm);
    }

    private void OnWhole(object? s, RoutedEventArgs e) => Vm.PbShow = 1;
    private void OnCentre(object? s, RoutedEventArgs e) => Vm.CentreOnSelection();
}
