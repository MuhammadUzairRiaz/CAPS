using Avalonia.Controls;
using Avalonia.Input;
using Avalonia.Interactivity;
using Avalonia.Markup.Xaml;
using CapsStudio.ViewModels;

namespace CapsStudio.Views.Pages;

public partial class ImportDialog : PageBase
{
    private MainViewModel? _hooked;

    public ImportDialog()
    {
        AvaloniaXamlLoader.Load(this);
        DataContextChanged += (_, _) =>
        {
            if (DataContext is not MainViewModel vm || vm == _hooked) return;
            _hooked = vm;
            vm.ImportFragmentChanged += () =>
            {
                var view = this.FindControl<MolView>("Frag")!;
                view.Document = vm.ImportFragmentDoc;
                view.Reset();
            };
        };
    }

    private void OnBackdrop(object? s, PointerPressedEventArgs e) => Vm.CloseImport();
    private void OnCancel(object? s, RoutedEventArgs e) => Vm.CloseImport();
    private async void OnImport(object? s, RoutedEventArgs e) => await Vm.ConfirmImport();
}
