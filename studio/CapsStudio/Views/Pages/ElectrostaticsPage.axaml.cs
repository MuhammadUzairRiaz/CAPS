using Avalonia;
using Avalonia.Controls;
using Avalonia.Interactivity;
using Avalonia.Markup.Xaml;
using CapsStudio.ViewModels;

namespace CapsStudio.Views.Pages;

public partial class ElectrostaticsPage : PageBase
{
    private MainViewModel? _hooked;
    public ElectrostaticsPage()
    {
        AvaloniaXamlLoader.Load(this);
        DataContextChanged += (_, _) =>
        {
            if (DataContext is not MainViewModel vm || vm == _hooked) return;
            _hooked = vm;
            vm.EsChanged += () => Update(vm);
        };
    }
    private void Update(MainViewModel vm)
    {
        var rc = (double)vm.EsCutoff;
        var tol = vm.EsTolValue;
        this.FindControl<SciPlot>("Curve")!.Figure = ViewModels.PageFigures.Electrostatics(vm.EsCurve, rc, tol, vm.EsPme);
    }
    protected override void OnPropertyChanged(AvaloniaPropertyChangedEventArgs change)
    {
        base.OnPropertyChanged(change);
        if (change.Property == IsVisibleProperty && IsVisible && DataContext is MainViewModel vm) Update(vm);
    }
    private void OnApply(object? s, RoutedEventArgs e) => Vm.EsApply();
    private void OnReset(object? s, RoutedEventArgs e) => Vm.EsReset();
}
