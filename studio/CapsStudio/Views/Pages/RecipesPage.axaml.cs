using Avalonia.Controls;
using Avalonia.Interactivity;
using Avalonia.Markup.Xaml;
using CapsStudio.ViewModels;

namespace CapsStudio.Views.Pages;

public partial class RecipesPage : PageBase
{
    private MainViewModel? _hooked;

    public RecipesPage()
    {
        AvaloniaXamlLoader.Load(this);
        DataContextChanged += (_, _) =>
        {
            if (DataContext is not MainViewModel vm || vm == _hooked) return;
            _hooked = vm;
            vm.RecipeChanged += () =>
            {
                var t = this.FindControl<LinePlot>("Temp")!;
                t.RefY = null;
                t.SetData(vm.RecipeT);
                var p = this.FindControl<LinePlot>("Pres")!;
                p.RefY = null;
                p.SetData(vm.RecipeP);
            };
        };
    }

    private void OnSave(object? s, RoutedEventArgs e) => Vm.SaveRecipe();
    private void OnDuplicate(object? s, RoutedEventArgs e) => Vm.DuplicateRecipe();
    private async void OnRun(object? s, RoutedEventArgs e) => await Vm.RunSelectedRecipe();
    private async void OnRunOnDoc(object? s, RoutedEventArgs e) => await Vm.RunRecipeOnDocument();
}
