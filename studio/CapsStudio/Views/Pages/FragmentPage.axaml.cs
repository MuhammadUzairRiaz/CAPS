using Avalonia.Controls;
using Avalonia.Input;
using Avalonia.Interactivity;
using Avalonia.Markup.Xaml;
using CapsStudio.ViewModels;

namespace CapsStudio.Views.Pages;

public partial class FragmentPage : PageBase
{
    public FragmentPage()
    {
        AvaloniaXamlLoader.Load(this);
        var detail = this.FindControl<MolView>("Detail")!;
        DataContextChanged += (_, _) =>
        {
            if (DataContext is not MainViewModel vm) return;
            vm.FragmentChanged += () =>
            {
                var f = vm.SelectedFragment;
                detail.Smiles = f?.SmilesH;
                detail.Highlights = f?.StarIndices ?? [];
                var cats = this.FindControl<ListBox>("Cats")!;
                var want = vm.FragmentCategories.FirstOrDefault(c => c.Name == vm.FragmentCategory);
                if (want != null && !Equals(cats.SelectedItem, want)) cats.SelectedItem = want;
            };
        };
    }

    private void OnBack(object? s, RoutedEventArgs e) => Vm.SetModule(8);
    private void OnCategory(object? s, SelectionChangedEventArgs e) { if ((s as ListBox)?.SelectedItem is FragmentCategory c) { Vm.FragmentQuery = ""; Vm.FragmentCategory = c.Name; } }
    private async void OnUse(object? s, RoutedEventArgs e) => await Vm.UseFragment();
    private async void OnTileDouble(object? s, TappedEventArgs e) => await Vm.UseFragment();
    private void OnSaveCopy(object? s, RoutedEventArgs e) => Vm.SaveFragmentCopy();
    private void OnAddMine(object? s, RoutedEventArgs e)
    {
        var name = this.FindControl<TextBox>("NewName")!.Text ?? "";
        var smiles = this.FindControl<TextBox>("NewSmiles")!.Text ?? "";
        Vm.AddMyFragment(name, smiles);
    }
    private void OnFromSelection(object? s, RoutedEventArgs e) => Vm.SaveSelectionAsFragment();
}
