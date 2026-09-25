using Avalonia;
using Avalonia.Controls;
using Avalonia.Interactivity;
using Avalonia.Markup.Xaml;
using CapsStudio.ViewModels;

namespace CapsStudio.Views.Pages;

public partial class ModelResolutionPage : PageBase
{
    private MainViewModel? _hooked;
    public ModelResolutionPage()
    {
        AvaloniaXamlLoader.Load(this);
        DataContextChanged += (_, _) =>
        {
            if (DataContext is not MainViewModel vm || vm == _hooked) return;
            _hooked = vm;
            vm.PropertyChanged += (_, e) => { if (IsVisible && e.PropertyName is nameof(MainViewModel.MrCgDoc) or nameof(MainViewModel.Document)) Show(vm); };
        };
    }
    private void Show(MainViewModel vm)
    {
        foreach (var (name, doc, style) in new[] { ("Aa", vm.Document, 0), ("Ua", vm.MrUaDoc, 0), ("Cg", vm.MrCgDoc, 1) })
        {
            var v = this.FindControl<MolView>(name)!;
            v.DrawStyle = style; v.ColourMode = name == "Cg" ? 1 : 0; v.ShowCell = false;
            if (v.Document == doc) { v.Refresh(); continue; }
            v.Document = doc;
            v.Reset();
            if (doc != null)   // frame the atoms, not a roomy cell
                try { v.Camera = doc.Focus(v.Camera, System.Linq.Enumerable.Range(0, (int)doc.Summary().Atoms).ToArray(), 0.8); } catch { }
        }
    }
    protected override void OnPropertyChanged(AvaloniaPropertyChangedEventArgs change)
    {
        base.OnPropertyChanged(change);
        if (change.Property == IsVisibleProperty && IsVisible && DataContext is MainViewModel vm) Show(vm);
    }
    private void OnBackmap(object? s, RoutedEventArgs e) => Vm.SetModule(44);
}
