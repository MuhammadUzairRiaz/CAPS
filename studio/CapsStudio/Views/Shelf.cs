using Avalonia;
using Avalonia.Controls;
using Avalonia.Controls.Primitives;
using Avalonia.Controls.Shapes;
using Avalonia.Input;
using Avalonia.Layout;
using Avalonia.Media;
using Avalonia.LogicalTree;
using Avalonia.VisualTree;

namespace CapsStudio.Views;

/// <summary>A tool shelf (design/boards/Shelves): an amber tab, then its tools. Docked in a row of the toolbar it is flat,
/// like the row it replaced; floating, popped out or peeking from its puck it is a capsule with a shadow. The tools are
/// the window's own controls, moved in (their bindings and handlers go with them), or proxies of them on a shelf the user
/// made. Drag the tab to move the shelf, click it to fold the shelf into a puck, double-click to pop it out.</summary>
public sealed class Shelf : Border
{
    protected override Type StyleKeyOverride => typeof(Border);
    public string Id { get; }
    public string Title { get; }
    public string Glyph { get; }
    /// <summary>The panel holding the tools (the original row).</summary>
    public Panel Tools { get; }
    public Button Tab { get; }
    private readonly DockPanel _inner;
    private bool _vertical, _floating;

    public Shelf(string id, string title, string glyph, Panel tools)
    {
        Id = id; Title = title; Glyph = glyph; Tools = tools;
        Classes.Add("shelf");
        Tab = new Button
        {
            Classes = { "shelftab" },
            Content = new Icon { Kind = glyph, Size = 12, StrokeWidth = 2 },
            [ToolTip.TipProperty] = $"{title} · drag to move or dock · click to fold · double-click to pop out",
            [Avalonia.Automation.AutomationProperties.NameProperty] = $"Move the {title} shelf",
        };
        _inner = new DockPanel { LastChildFill = true };
        DockPanel.SetDock(Tab, Dock.Left);
        _inner.Children.Add(Tab);
        _inner.Children.Add(tools);
        Child = _inner;
    }

    public bool Vertical
    {
        get => _vertical;
        set
        {
            _vertical = value;
            Classes.Set("vertical", value);
            DockPanel.SetDock(Tab, value ? Dock.Top : Dock.Left);
            switch (Tools)
            {
                case WrapPanel w: w.Orientation = value ? Orientation.Vertical : Orientation.Horizontal; break;
                case StackPanel s: s.Orientation = value ? Orientation.Vertical : Orientation.Horizontal; break;
            }
            // inner rows (the element buttons of Modify) turn with it
            foreach (var p in Tools.GetVisualDescendants().OfType<Panel>().Where(p => p.Classes.Contains("shelfflow")))
                if (p is StackPanel sp) sp.Orientation = value ? Orientation.Vertical : Orientation.Horizontal;
                else if (p is WrapPanel wp) wp.Orientation = value ? Orientation.Vertical : Orientation.Horizontal;
            Tools.Classes.Set("compact", value);
            // a label beside an icon goes on a vertical shelf (the icon says it, the tooltip names it); a text-only tool keeps its text
            foreach (var b in Tools.GetLogicalDescendants().OfType<Control>().Where(c => c is Button || c is ToggleButton))
                if (b.GetLogicalDescendants().OfType<Icon>().Any())
                    foreach (var t in b.GetLogicalDescendants().OfType<TextBlock>()) t.Classes.Set("shelflabel", true);
        }
    }

    public bool Floating { get => _floating; set { _floating = value; Classes.Set("floating", value); } }

    /// <summary>How many tools it holds (the puck's badge).</summary>
    public int ToolCount
    {
        get
        {
            var n = Tools.GetVisualDescendants().OfType<Control>().Count(c => (c is Button || c is ToggleButton) && c.IsVisible);
            return n > 0 ? n : Tools.Children.Count;
        }
    }
}

/// <summary>A folded shelf: a round puck on the edge it was docked to. Hovering it opens the shelf beside it (a peek);
/// clicking unfolds the shelf where it was.</summary>
public sealed class Puck : Button
{
    protected override Type StyleKeyOverride => typeof(Button);
    public Shelf Shelf { get; }
    public Puck(Shelf shelf, int count)
    {
        Shelf = shelf;
        Classes.Add("puck");
        var badge = new Border
        {
            Classes = { "puckbadge" },
            Child = new TextBlock { Text = count.ToString(System.Globalization.CultureInfo.InvariantCulture) },
            HorizontalAlignment = HorizontalAlignment.Right, VerticalAlignment = VerticalAlignment.Top, Margin = new Thickness(0, -6, -6, 0),
        };
        Content = new Panel { Children = { new Icon { Kind = shelf.Glyph, Size = 15, HorizontalAlignment = HorizontalAlignment.Center, VerticalAlignment = VerticalAlignment.Center }, badge } };
        ToolTip.SetTip(this, $"{shelf.Title} ({count} tools) · hover to peek · click to unfold");
        Avalonia.Automation.AutomationProperties.SetName(this, $"Open the {shelf.Title} shelf");
    }
}

/// <summary>A tool on a shelf the user made: it shows the original tool's icon or text and tooltip, is enabled and lit
/// when the original is, and clicking it clicks the original (toggles a toggle).</summary>
public sealed class ProxyTool : Button
{
    protected override Type StyleKeyOverride => typeof(Button);
    public string ToolId { get; }
    public ProxyTool(string id, Control original)
    {
        ToolId = id;
        Classes.Add("tool");
        Content = ShelfTools.Face(original);
        ToolTip.SetTip(this, ToolTip.GetTip(original));
        Avalonia.Automation.AutomationProperties.SetName(this, ShelfTools.Label(original));
        IsEnabled = original.IsEnabled;
        original.PropertyChanged += (_, e) =>
        {
            if (e.Property == IsEnabledProperty || e.Property == IsEffectivelyEnabledProperty) IsEnabled = original.IsEffectivelyEnabled;
            if (e.Property == ToggleButton.IsCheckedProperty) Classes.Set("on", original is ToggleButton { IsChecked: true });
        };
        original.Classes.CollectionChanged += (_, _) => Classes.Set("on", original.Classes.Contains("on") || original is ToggleButton { IsChecked: true });
        Classes.Set("on", original.Classes.Contains("on") || original is ToggleButton { IsChecked: true });
        Click += (_, _) =>
        {
            if (original is ToggleButton t) t.IsChecked = !(t.IsChecked ?? false);
            else if (original is Button b)
            {
                if (b.Flyout is { } f) f.ShowAt(this);
                else b.RaiseEvent(new Avalonia.Interactivity.RoutedEventArgs(Button.ClickEvent, b));
            }
        };
    }
}

/// <summary>The tools of the built-in shelves, by a stable id (their tooltip's first words), for shelves the user makes.</summary>
public static class ShelfTools
{
    public static string Label(Control c)
    {
        var tip = ToolTip.GetTip(c) as string ?? "";
        var cut = tip.IndexOfAny([':', '(', '·', '—', ',']);
        var s = (cut > 0 ? tip[..cut] : tip).Trim();
        if (s.Length == 0) s = c.GetVisualDescendants().OfType<TextBlock>().Select(t => t.Text).FirstOrDefault(t => !string.IsNullOrWhiteSpace(t)) ?? "tool";
        return s.Length > 40 ? s[..40] : s;
    }
    public static string Id(Control c) => new string(Label(c).ToLowerInvariant().Select(ch => char.IsLetterOrDigit(ch) ? ch : '_').ToArray()).Trim('_');

    /// <summary>A copy of the tool's face: its icon, else its text.</summary>
    public static Control Face(Control c)
    {
        var icon = c.GetVisualDescendants().OfType<Icon>().FirstOrDefault() ?? c.GetLogicalDescendants().OfType<Icon>().FirstOrDefault();
        if (icon != null) return new Icon { Kind = icon.Kind, Size = 18 };
        var text = c.GetLogicalDescendants().OfType<TextBlock>().Select(t => t.Text).FirstOrDefault(t => !string.IsNullOrWhiteSpace(t)) ?? "?";
        return new TextBlock { Text = text, FontWeight = FontWeight.SemiBold, FontSize = 12, VerticalAlignment = VerticalAlignment.Center };
    }

    /// <summary>The tool controls of a shelf: buttons and toggles, not the ones inside a flyout or a template.</summary>
    public static IEnumerable<Control> Of(Panel tools) =>
        tools.GetLogicalDescendants().OfType<Control>().Where(c => c is Button or ToggleButton && c is not ProxyTool && !c.Classes.Contains("shelftab")
                                                                 && ToolTip.GetTip(c) is string { Length: > 0 });
}

