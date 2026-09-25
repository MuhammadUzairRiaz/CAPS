using Avalonia;
using Avalonia.Controls;
using Avalonia.Interactivity;
using Avalonia.Markup.Xaml;
using CapsStudio.ViewModels;

namespace CapsStudio.Views.Pages;

public partial class DisplayStylesPage : PageBase
{
    private MainViewModel? _hooked;
    public DisplayStylesPage()
    {
        AvaloniaXamlLoader.Load(this);
        DataContextChanged += (_, _) =>
        {
            if (DataContext is not MainViewModel vm || vm == _hooked) return;
            _hooked = vm;
            vm.PropertyChanged += (_, e) =>
            {
                if (IsVisible && e.PropertyName is nameof(MainViewModel.DsStyle) or nameof(MainViewModel.ColourIndex) or nameof(MainViewModel.DisplayStatus)) Show(vm);
            };
        };
    }
    private void Show(MainViewModel vm)
    {
        foreach (var (name, style, card, k) in new[] { ("All", vm.DsAtoms, "CardAll", 0), ("NoH", 3, "CardNoH", 1), ("Bb", 4, "CardBb", 2) })
        {
            var v = this.FindControl<MolView>(name)!;
            v.DrawStyle = style; v.ColourMode = vm.ColourIndex; v.ShowCell = true;
            if (v.Document != vm.Document) { v.Document = vm.Document; v.Reset(); } else v.Refresh();
            this.FindControl<Border>(card)!.Classes.Set("on", vm.DsStyle == k);
        }
    }
    protected override void OnPropertyChanged(AvaloniaPropertyChangedEventArgs change)
    {
        base.OnPropertyChanged(change);
        if (change.Property == IsVisibleProperty && IsVisible && DataContext is MainViewModel vm) Show(vm);
    }
    private void OnStyle(object? s, RoutedEventArgs e) { if (s is Control { Tag: string t }) Vm.DsStyle = int.Parse(t); }
    private void OnReset(object? s, RoutedEventArgs e) => Vm.ResetDisplay();
    private void OnSaveDefault(object? s, RoutedEventArgs e) => Vm.SaveDisplayDefault();
}
