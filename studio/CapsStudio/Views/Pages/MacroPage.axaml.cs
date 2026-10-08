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
            vm.MacroGoToLine += line =>
            {
                // the caret at the start of that line, the editor focused
                var code = this.FindControl<TextBox>("Code")!;
                var text = code.Text ?? "";
                var at = 0;
                for (var l = 1; l < line && at >= 0; ++l) { at = text.IndexOf('\n', at); if (at >= 0) ++at; }
                if (at < 0) at = text.Length;
                code.Focus();
                code.CaretIndex = at;
                var end = text.IndexOf('\n', at);
                code.SelectionStart = at;
                code.SelectionEnd = end < 0 ? text.Length : end;
            };
        };
    }

    private void OnBack(object? s, RoutedEventArgs e) => Vm.SetModule(8);
    private void OnNew(object? s, RoutedEventArgs e) => Vm.NewMacro();
    private void OnSave(object? s, RoutedEventArgs e) => Vm.SaveMacro();
    private async void OnRun(object? s, RoutedEventArgs e) => await Vm.RunMacro();
    private void OnStop(object? s, RoutedEventArgs e) => Vm.StopMacro();
    private async void OnCheck(object? s, RoutedEventArgs e) => await Vm.CheckMacro();
    private void OnOutTab(object? s, RoutedEventArgs e) => Vm.MacroOutTab = (s as Control)?.Tag as string == "1" ? 1 : 0;
    private void OnProblem(object? s, RoutedEventArgs e) { if ((s as Control)?.Tag is MacroProblem p) Vm.GoToProblem(p); }
    private void OnPromote(object? s, RoutedEventArgs e)
    {
        var code = this.FindControl<TextBox>("Code")!;
        var note = this.FindControl<TextBlock>("PromoteNote")!;
        var why = Vm.PromoteToParameter(code.SelectedText ?? "");
        note.Text = why ?? "Promoted: it is now an argument of macro(…) with its old value as the default.";
    }
}
