using Avalonia.Controls;
using Avalonia.Input;
using Avalonia.Interactivity;
using Avalonia.Markup.Xaml;
using CapsStudio.ViewModels;

namespace CapsStudio.Views.Pages;

public partial class UnitTypingPage : PageBase
{
    private MainViewModel? _hooked;
    public UnitTypingPage()
    {
        AvaloniaXamlLoader.Load(this);
        var view = this.FindControl<MolView>("View")!;
        view.AtomClicked += i => Vm.PickUnitAtom(i);
        DataContextChanged += (_, _) =>
        {
            if (DataContext is not MainViewModel vm || vm == _hooked) return;
            _hooked = vm;
            vm.UtChanged += () => Show(vm);
            vm.PropertyChanged += (_, e) =>
            {
                if (e.PropertyName is nameof(MainViewModel.UtDoc)) { view.Document = vm.UtDoc; view.Reset(); Show(vm); }
                if (e.PropertyName is nameof(MainViewModel.UtHighlights)) Show(vm);
                if (e.PropertyName is nameof(MainViewModel.UtSelected) && vm.UtSelected != null) this.FindControl<ListBox>("Atoms")?.ScrollIntoView(vm.UtSelected);
            };
        };
    }

    private void Show(MainViewModel vm)
    {
        var view = this.FindControl<MolView>("View")!;
        if (view.Document != vm.UtDoc) { view.Document = vm.UtDoc; view.Reset(); }
        view.ColourMode = 2;   // by force-field type
        view.Highlights = vm.UtHighlights;
        view.Refresh();
    }

    private void OnAuto(object? s, RoutedEventArgs e) => Vm.AutoTypeUnitExample();
    private void OnRebuild(object? s, RoutedEventArgs e) => Vm.BuildUnitExample();
    private void OnApply(object? s, RoutedEventArgs e) => Vm.ApplyUnitTyping();
    private void OnAssign(object? s, RoutedEventArgs e) => Vm.AssignUnitType();
    private void OnClear(object? s, RoutedEventArgs e) => Vm.ClearUnitType();
    private void OnTypeDouble(object? s, TappedEventArgs e) => Vm.AssignUnitType();
}
