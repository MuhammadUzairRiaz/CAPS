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
                var sp = vm.SorbSpeciesIsotherms;
                if (sp.Count >= 2)   // a mixture: the gases' own isotherms (up to three)
                {
                    p.SetCompare(sp[0], sp[1]);
                    if (sp.Count >= 3) p.SetThird(sp[2]);
                }
                else p.SetData(vm.SorbIsotherm);
            };
            vm.SorptionMapChanged += () =>
            {
                var m = this.FindControl<LinePlot>("Map")!;
                var (x, y, z, xl, yl) = vm.SorbMapProjection();
                m.XLabel = xl;
                m.YLabel = yl;
                m.SetHeat(x, y, z);
            };
        };
    }

    private async void OnRun(object? s, RoutedEventArgs e) => await Vm.RunSorption();
    private void OnCancel(object? s, RoutedEventArgs e) => Vm.CancelSorption();
}
