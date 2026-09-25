using Avalonia.Controls;
using Avalonia.Interactivity;
using Avalonia.Markup.Xaml;
using CapsStudio.ViewModels;

namespace CapsStudio.Views.Pages;

public partial class MacroPage : PageBase
{
    public MacroPage()
    {
        AvaloniaXamlLoader.Load(this);
        DataContextChanged += (_, _) =>
        {
            if (DataContext is not MainViewModel vm) return;
            vm.PropertyChanged += (_, e) =>
            {
                if (e.PropertyName == nameof(MainViewModel.MacroOutput)) this.FindControl<ScrollViewer>("OutScroll")?.ScrollToEnd();
            };
        };
    }

    private void OnBack(object? s, RoutedEventArgs e) => Vm.SetModule(8);
    private void OnNew(object? s, RoutedEventArgs e) => Vm.NewMacro();
    private void OnSave(object? s, RoutedEventArgs e) => Vm.SaveMacro();
    private async void OnRun(object? s, RoutedEventArgs e) => await Vm.RunMacro();
    private void OnStop(object? s, RoutedEventArgs e) => Vm.StopMacro();
    private void OnPromote(object? s, RoutedEventArgs e)
    {
        var code = this.FindControl<TextBox>("Code")!;
        var note = this.FindControl<TextBlock>("PromoteNote")!;
        var why = Vm.PromoteToParameter(code.SelectedText ?? "");
        note.Text = why ?? "Promoted: it is now an argument of macro(…) with its old value as the default.";
    }
}
