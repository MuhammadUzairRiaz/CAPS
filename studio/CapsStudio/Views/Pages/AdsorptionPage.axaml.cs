using Avalonia;
using Avalonia.Controls;
using Avalonia.Interactivity;
using Avalonia.Markup.Xaml;
using CapsStudio.ViewModels;

namespace CapsStudio.Views.Pages;

public partial class AdsorptionPage : PageBase
{
    private MainViewModel? _hooked;

    public AdsorptionPage()
    {
        AvaloniaXamlLoader.Load(this);
        DataContextChanged += (_, _) =>
        {
            if (DataContext is not MainViewModel vm || vm == _hooked) return;
            _hooked = vm;
            vm.AdsorptionChanged += () =>
            {
                ShowView(vm);
                var h = this.FindControl<LinePlot>("Hist")!;
                h.RefY = null;
                h.SetData(vm.AdsHistogram);
            };
        };
    }

    private void ShowView(MainViewModel vm)
    {
        var view = this.FindControl<MolView>("View")!;
        if (view.Document != vm.Document) { view.Document = vm.Document; view.Reset(); }
        else view.Refresh();
    }

    protected override void OnPropertyChanged(AvaloniaPropertyChangedEventArgs change)
    {
        base.OnPropertyChanged(change);
        if (change.Property == IsVisibleProperty && IsVisible && DataContext is MainViewModel vm) ShowView(vm);
    }

    private async void OnRun(object? s, RoutedEventArgs e) => await Vm.RunAdsorption();
    private void OnCancel(object? s, RoutedEventArgs e) => Vm.CancelAdsorption();
    private void OnAdd(object? s, RoutedEventArgs e) => Vm.AddAdsorbate();
    private void OnRemove(object? s, RoutedEventArgs e) { if ((s as Control)?.DataContext is AdsorbateRow r) Vm.RemoveAdsorbate(r); }
}
