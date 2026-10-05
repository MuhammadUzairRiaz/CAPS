using Avalonia.Controls;
using Avalonia.Input;
using Avalonia.Interactivity;
using CapsStudio.ViewModels;

namespace CapsStudio.Views;

/// <summary>Brush to select (design/boards/BrushSelect).</summary>
public partial class MainWindow
{
    private void OnBrushClose(object? s, RoutedEventArgs e) => _vm.BrushOpen = false;
    // Bond rules (design/boards/BondRules)
    private void OnRulesClose(object? s, RoutedEventArgs e) => _vm.RulesOpen = false;
    private void OnRulesApply(object? s, RoutedEventArgs e) { _vm.ApplyBondRules(); RequestRender(); }
    private void OnRulesSave(object? s, RoutedEventArgs e) => _vm.SaveRuleSet();
    private void OnBrushTag(object? s, RoutedEventArgs e) => _vm.TagSelection();
    private void OnBrushClear(object? s, RoutedEventArgs e) => _vm.ClearBrushes();
    private void OnBrushRemove(object? s, RoutedEventArgs e) { if ((s as Control)?.Tag is BrushHist h) _vm.RemoveBrushHist(h); }
    private void OnBrushQueryKey(object? s, KeyEventArgs e) { if (e.Key == Key.Enter) { e.Handled = true; _vm.RunBrushQuery(); } }

    private void OnBrushAdd(object? s, RoutedEventArgs e)
    {
        var items = MainViewModel.BrushChoices.Where(c => _vm.BrushHists.All(h => h.Key != c.Key))
            .Select(c => { var m = new MenuItem { Header = c.Title + (c.Unit.Length > 0 ? $" ({c.Unit})" : "") }; m.Click += (_, _) => _vm.AddBrushHist(c.Key); return m; }).ToList();
        if (items.Count == 0) return;
        new ContextMenu { ItemsSource = items }.Open(s as Control);
    }

    private void OnBrushDistance(object? s, RoutedEventArgs e)
    {
        var items = new List<MenuItem>();
        var atom = new MenuItem { Header = "From the picked atom" };
        atom.Click += (_, _) => _vm.AddBrushDistance(null);
        items.Add(atom);
        foreach (var t in _vm.TagChips)
        {
            var m = new MenuItem { Header = "From the tag " + t.Name };
            var name = t.Name;
            m.Click += (_, _) => _vm.AddBrushDistance(name);
            items.Add(m);
        }
        new ContextMenu { ItemsSource = items }.Open(s as Control);
    }
}
