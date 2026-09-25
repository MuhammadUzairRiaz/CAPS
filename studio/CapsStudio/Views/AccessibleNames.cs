using System;
using System.Collections.Generic;
using System.Linq;
using Avalonia;
using Avalonia.Automation;
using Avalonia.Controls;
using Avalonia.LogicalTree;
using Avalonia.VisualTree;

namespace CapsStudio.Views;

/// <summary>Accessibility map (design/boards/AccessibilityMap): every button a screen reader reaches has a name. Buttons
/// with text are read by their text; icon-only buttons take the first clause of their tooltip ("Export image or movie
/// (⌘E)" → "Export image or movie"), set once when they load, unless a name was given in XAML.</summary>
public static class AccessibleNames
{
    public static void Install() => Control.LoadedEvent.AddClassHandler<Button>((b, _) => Name(b));

    public static void Name(Button b)
    {
        if (!string.IsNullOrEmpty(AutomationProperties.GetName(b)) || HasText(b)) return;
        if (ToolTip.GetTip(b) is string tip && From(tip) is { Length: > 0 } n) AutomationProperties.SetName(b, n);
    }

    /// <summary>The spoken name from a tooltip: its first clause, without shortcuts in brackets, at most 60 characters.</summary>
    public static string From(string tip)
    {
        var t = tip.Trim();
        foreach (var cut in new[] { " — ", ": ", " · ", " (", "; ", ". " })
        {
            var i = t.IndexOf(cut, StringComparison.Ordinal);
            if (i > 2) t = t[..i];
        }
        t = t.TrimEnd('.', ' ');
        return t.Length > 60 ? t[..60].TrimEnd() : t;
    }

    private static bool HasText(Button b) => b.Content switch
    {
        string s => s.Trim().Length > 0,
        TextBlock tb => !string.IsNullOrWhiteSpace(tb.Text),
        Control c => c.GetLogicalDescendants().OfType<TextBlock>().Any(t => !string.IsNullOrWhiteSpace(t.Text) && t.IsVisible),
        _ => false,
    };

    /// <summary>Visible buttons a screen reader would reach without a name (for checks).</summary>
    public static List<Button> Unnamed(Visual root) =>
        root.GetVisualDescendants().OfType<Button>()
            .Where(b => b.IsEffectivelyVisible && b.IsEffectivelyEnabled && string.IsNullOrEmpty(AutomationProperties.GetName(b)) && !HasText(b)).ToList();
}
