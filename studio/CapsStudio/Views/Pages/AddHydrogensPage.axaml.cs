using Avalonia;
using Avalonia.Controls;
using Avalonia.Interactivity;
using Avalonia.Markup.Xaml;
using CapsStudio.ViewModels;

namespace CapsStudio.Views.Pages;

public partial class AddHydrogensPage : PageBase
{
    private MainViewModel? _hooked;
    public AddHydrogensPage()
    {
        AvaloniaXamlLoader.Load(this);
        DataContextChanged += (_, _) =>
        {
            if (DataContext is not MainViewModel vm || vm == _hooked) return;
            _hooked = vm;
            vm.PropertyChanged += (_, e) => { if (e.PropertyName is nameof(MainViewModel.AhAfterDoc)) Show(vm); };
        };
    }
    private void Show(MainViewModel vm)
    {
        foreach (var (name, doc) in new[] { ("Before", vm.AhBeforeDoc), ("After", vm.AhAfterDoc) })
        {
            var v = this.FindControl<MolView>(name)!;
            v.DrawStyle = 0; v.ColourMode = 0; v.ShowCell = false;
            if (v.Document != doc) { v.Document = doc; v.Reset(); } else v.Refresh();
        }
    }
    protected override void OnPropertyChanged(AvaloniaPropertyChangedEventArgs change)
    {
        base.OnPropertyChanged(change);
        if (change.Property == IsVisibleProperty && IsVisible && DataContext is MainViewModel vm) Show(vm);
    }
    private async void OnApply(object? s, RoutedEventArgs e) => await Vm.AhApply();
    private async void OnUndo(object? s, RoutedEventArgs e) => await Vm.AhUndo();
}
