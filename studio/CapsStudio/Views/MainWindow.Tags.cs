using Avalonia.Controls;
using Avalonia.Input;
using Avalonia.Interactivity;
using CapsStudio.ViewModels;

namespace CapsStudio.Views;

/// <summary>The tag strip over the view (design/boards/Tags).</summary>
public partial class MainWindow
{
    private static TagChip? ChipOf(object? s) => (s as Control)?.Tag as TagChip;
    private KeyModifiers _lastKeyModifiers;   // the keys held at the last press (a click carries none)

    private void InitTags() => AddHandler(PointerPressedEvent, (_, e) => _lastKeyModifiers = e.KeyModifiers, RoutingStrategies.Tunnel, handledEventsToo: true);

    private void OnTagChip(object? s, RoutedEventArgs e)
    {
        if (ChipOf(s) is not { } chip) return;
        var mods = _lastKeyModifiers;
        _vm.UseTag(chip, mods.HasFlag(KeyModifiers.Shift), mods.HasFlag(KeyModifiers.Alt));
        RequestRender();
    }
    private void OnTagNew(object? s, RoutedEventArgs e) => _vm.TagSelection();
    private void OnTagRename(object? s, RoutedEventArgs e) { if (ChipOf(s) is { } c) c.Renaming = true; }
    private void OnTagAddSel(object? s, RoutedEventArgs e) { if (ChipOf(s) is { } c) _vm.TagWithSelection(c, false); }
    private void OnTagRemoveSel(object? s, RoutedEventArgs e) { if (ChipOf(s) is { } c) _vm.TagWithSelection(c, true); }
    private void OnTagDelete(object? s, RoutedEventArgs e) { if (ChipOf(s) is { } c) _vm.DeleteTag(c); RequestRender(); }

    private void OnTagColourMenu(object? s, RoutedEventArgs e)
    {
        if (s is not MenuItem m || m.Tag is not TagChip chip || m.Items.Count > 0) return;
        foreach (var col in MainViewModel.TagColours)
        {
            var it = new MenuItem { Header = new Avalonia.Controls.Shapes.Ellipse { Width = 12, Height = 12, Fill = Avalonia.Media.SolidColorBrush.Parse(col) } };
            it.Click += (_, _) => _vm.RecolourTag(chip, col);
            m.Items.Add(it);
        }
    }

    private void OnTagNameShown(object? s, Avalonia.VisualTreeAttachmentEventArgs e)
    {
        if (s is TextBox b && b.IsVisible) Avalonia.Threading.Dispatcher.UIThread.Post(() => { b.Focus(); b.SelectAll(); });
    }
    private void OnTagNameKey(object? s, KeyEventArgs e)
    {
        if (s is not TextBox b || b.Tag is not TagChip chip) return;
        if (e.Key == Key.Enter) { e.Handled = true; _vm.RenameTag(chip, b.Text ?? ""); ViewHost.Focus(); }
        else if (e.Key == Key.Escape) { e.Handled = true; chip.Renaming = false; ViewHost.Focus(); }
    }
    private void OnTagNameLost(object? s, RoutedEventArgs e)
    {
        if (s is TextBox b && b.Tag is TagChip { Renaming: true } chip) _vm.RenameTag(chip, b.Text ?? "");
    }
}
