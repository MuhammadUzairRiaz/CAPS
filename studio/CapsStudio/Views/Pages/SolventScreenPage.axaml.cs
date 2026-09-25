using System.Linq;
using Avalonia;
using Avalonia.Controls;
using Avalonia.Markup.Xaml;
using CapsStudio.ViewModels;

namespace CapsStudio.Views.Pages;

public partial class SolventScreenPage : PageBase
{
    private MainViewModel? _hooked;
    public SolventScreenPage()
    {
        AvaloniaXamlLoader.Load(this);
        DataContextChanged += (_, _) =>
        {
            if (DataContext is not MainViewModel vm || vm == _hooked) return;
            _hooked = vm;
            vm.SsRows.CollectionChanged += (_, _) => Update(vm);
        };
    }
    private void Update(MainViewModel vm)
    {
        var bars = vm.SsRows.Select(r => (r.Name, double.Parse(r.Chi, System.Globalization.CultureInfo.InvariantCulture), r.Good ? "OkB" : r.Bad ? "ErrB" : "WarnB", r.Chi)).ToArray();
        this.FindControl<HBarChart>("Bars")!.Set(bars, 0.5, "χ = 0.5");
    }
    protected override void OnPropertyChanged(AvaloniaPropertyChangedEventArgs change)
    {
        base.OnPropertyChanged(change);
        if (change.Property == IsVisibleProperty && IsVisible && DataContext is MainViewModel vm) Update(vm);
    }
    private async void OnContacts(object? s, Avalonia.Interactivity.RoutedEventArgs e) { if (DataContext is ViewModels.MainViewModel vm) await vm.SsComputeContacts(); }
    private void OnStopContacts(object? s, Avalonia.Interactivity.RoutedEventArgs e) { if (DataContext is ViewModels.MainViewModel vm) vm.SsStopContacts(); }
}
