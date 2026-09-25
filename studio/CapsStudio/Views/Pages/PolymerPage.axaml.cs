using Avalonia.Controls;
using Avalonia.Input;
using Avalonia.Interactivity;
using Avalonia.Markup.Xaml;
using CapsStudio.ViewModels;

namespace CapsStudio.Views.Pages;

public partial class PolymerPage : PageBase
{
    private PolyUnit? _slot;

    public PolymerPage()
    {
        AvaloniaXamlLoader.Load(this);
        var strip = this.FindControl<SequenceStrip>("Strip")!;
        var preview = this.FindControl<MolView>("Preview")!;
        DataContextChanged += (_, _) =>
        {
            if (DataContext is not MainViewModel vm) return;
            vm.LoadPolymerLibrary();
            vm.PolyStripChanged += () => strip.SetUnits(vm.PolyStripUnits);
            strip.SetUnits(vm.PolyStripUnits);
            vm.PropertyChanged += (_, e) =>
            {
                if (e.PropertyName == nameof(MainViewModel.PolyDoc)) { preview.Document = vm.PolyDoc; preview.Reset(); }
            };
        };
    }

    private void OnLibraryUse(object? s, TappedEventArgs e)
    {
        if (this.FindControl<ListBox>("Library")!.SelectedItem is LibraryEntry entry) Vm.UseLibrary(entry, _slot);
    }

    private void OnUnitPressed(object? s, PointerPressedEventArgs e) => _slot = (s as Control)?.Tag as PolyUnit;
    private void OnAddUnit(object? s, RoutedEventArgs e) { Vm.AddPolyUnit(); _slot = Vm.PolyUnits.LastOrDefault(); }
    private void OnRemoveUnit(object? s, RoutedEventArgs e) { if ((s as Control)?.Tag is PolyUnit u) { Vm.RemovePolyUnit(u); if (_slot == u) _slot = null; } }
    private async void OnPreview(object? s, RoutedEventArgs e) => await Vm.BuildPolyPreview();
    private void OnSendGrow(object? s, RoutedEventArgs e) => Vm.SendPolymerToGrow();
    private void OnArch(object? s, RoutedEventArgs e) { if (s is Control { Tag: string t }) Vm.PolyArch = int.Parse(t); }
    private void OnNetwork(object? s, RoutedEventArgs e) { Vm.SetModule(6); Vm.Status = "Networks: grow the cell, then crosslink it here (C–C, sulfur or peroxide cures, epoxy–amine)"; }
    private void OnBlend(object? s, RoutedEventArgs e) => Vm.OpenBlend();
    private async void OnBuildOne(object? s, RoutedEventArgs e) => await Vm.BuildPolymerInStudio();
}
