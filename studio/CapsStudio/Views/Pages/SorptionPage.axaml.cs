using Avalonia.Controls;
using Avalonia.Interactivity;
using Avalonia.Markup.Xaml;
using CapsStudio.ViewModels;

namespace CapsStudio.Views.Pages;

public partial class SorptionPage : PageBase
{
    private MainViewModel? _hooked;

    public SorptionPage()
    {
        AvaloniaXamlLoader.Load(this);
        DataContextChanged += (_, _) =>
        {
            if (DataContext is not MainViewModel vm || vm == _hooked) return;
            _hooked = vm;
            vm.SorptionChanged += () =>
            {
                var p = this.FindControl<LinePlot>("Iso")!;
                p.RefY = null;
                p.SetData(vm.SorbIsotherm);
            };
        };
    }

    private async void OnRun(object? s, RoutedEventArgs e) => await Vm.RunSorption();
    private void OnCancel(object? s, RoutedEventArgs e) => Vm.CancelSorption();
}
