using System.Globalization;
using System.Runtime.CompilerServices;
using Avalonia;
using Avalonia.Controls;
using Avalonia.Controls.Primitives;
using Avalonia.Controls.Shapes;
using Avalonia.Input;
using Avalonia.Layout;
using Avalonia.Media;
using Avalonia.VisualTree;

namespace CapsStudio.Views;

/// <summary>Resizable, collapsible panes. Panes.Split="True" on a Grid puts a divider on each boundary between its
/// columns (Panes.SplitRows between its rows): drag the line to resize; hovering it shows ◀ ▶ (▲ ▼) buttons that push
/// the divider that way — hiding the pane on that side, or bringing back the one hidden on the other; a double click
/// puts both panes back to their first sizes. Sizes are kept in the settings (Panes.Store) under the grid's key.
/// A child with Panes.Fixed="True" gets no divider beside it. Plain Avalonia: the same on macOS, Windows and Linux.</summary>
public static class Panes
{
    public static readonly AttachedProperty<bool> SplitProperty = AvaloniaProperty.RegisterAttached<Grid, bool>("Split", typeof(Panes));
    public static readonly AttachedProperty<bool> SplitRowsProperty = AvaloniaProperty.RegisterAttached<Grid, bool>("SplitRows", typeof(Panes));
    public static readonly AttachedProperty<string?> KeyProperty = AvaloniaProperty.RegisterAttached<Grid, string?>("Key", typeof(Panes));
    public static readonly AttachedProperty<bool> FixedProperty = AvaloniaProperty.RegisterAttached<Control, bool>("Fixed", typeof(Panes));
    public static bool GetSplit(Grid g) => g.GetValue(SplitProperty);
    public static void SetSplit(Grid g, bool v) => g.SetValue(SplitProperty, v);
    public static bool GetSplitRows(Grid g) => g.GetValue(SplitRowsProperty);
    public static void SetSplitRows(Grid g, bool v) => g.SetValue(SplitRowsProperty, v);
    public static string? GetKey(Grid g) => g.GetValue(KeyProperty);
    public static void SetKey(Grid g, string? v) => g.SetValue(KeyProperty, v);
    /// <summary>The size a hidden pane folds to (a dock keeps its tab row); 0 by default.</summary>
    public static readonly AttachedProperty<double> FoldedSizeProperty = AvaloniaProperty.RegisterAttached<Control, double>("FoldedSize", typeof(Panes));
    public static double GetFoldedSize(Control c) => c.GetValue(FoldedSizeProperty);
    public static void SetFoldedSize(Control c, double v) => c.SetValue(FoldedSizeProperty, v);
    public static bool GetFixed(Control c) => c.GetValue(FixedProperty);
    public static void SetFixed(Control c, bool v) => c.SetValue(FixedProperty, v);

    /// <summary>Where sizes are kept (the settings' PaneSizes) and how they are saved.</summary>
    public static Dictionary<string, string>? Store { get; set; }
    public static Action? Save { get; set; }
    /// <summary>A pane was hidden or shown by its divider: (grid key, index, shown).</summary>
    public static event Action<string, int, bool>? Toggled;

    private sealed class State
    {
        public bool Rows;
        public string Key = "";
        public string[] Defaults = [];
        public readonly Dictionary<int, string> Before = new();   // a hidden pane's size, to bring it back
        public readonly Dictionary<int, (int J, Control C, double Size)> Filled = new();   // hiding k let set-size pane j fill the space
        public readonly List<PaneDivider> Dividers = new();
    }
    private static readonly ConditionalWeakTable<Grid, State> States = new();

    static Panes()
    {
        SplitProperty.Changed.AddClassHandler<Grid>((g, e) => { if (e.NewValue is true) Hook(g, false); });
        SplitRowsProperty.Changed.AddClassHandler<Grid>((g, e) => { if (e.NewValue is true) Hook(g, true); });
    }

    private static void Hook(Grid g, bool rows)
    {
        if (g.IsAttachedToVisualTree()) { Install(g, rows); return; }
        void Once(object? s, VisualTreeAttachmentEventArgs e) { g.AttachedToVisualTree -= Once; Install(g, rows); }
        g.AttachedToVisualTree += Once;
    }

    // ---------------------------------------------------------------- sizes
    private static int Count(Grid g, bool rows) => rows ? g.RowDefinitions.Count : g.ColumnDefinitions.Count;
    private static GridLength Def(Grid g, bool rows, int k) => rows ? g.RowDefinitions[k].Height : g.ColumnDefinitions[k].Width;
    private static void SetDef(Grid g, bool rows, int k, GridLength v)
    {
        if (rows) g.RowDefinitions[k].Height = v; else g.ColumnDefinitions[k].Width = v;
    }
    private static double Actual(Grid g, bool rows, int k) => rows ? g.RowDefinitions[k].ActualHeight : g.ColumnDefinitions[k].ActualWidth;
    private static int Index(Control c, bool rows) => rows ? Grid.GetRow(c) : Grid.GetColumn(c);
    private static int Span(Control c, bool rows) => Math.Max(1, rows ? Grid.GetRowSpan(c) : Grid.GetColumnSpan(c));
    private static IEnumerable<Control> Panes_(Grid g) => g.Children.Where(c => c is not PaneDivider);

    /// <summary>An Auto track sized by its one child's own Width / Height (a panel with a set size): that size moves.</summary>
    private static Control? SizedChild(Grid g, bool rows, int k)
    {
        if (!Def(g, rows, k).IsAuto) return null;
        var kids = Panes_(g).Where(c => Index(c, rows) == k && Span(c, rows) == 1).ToList();
        // the one shown (another panel may share the track, hidden until its moment)
        bool Sized(Control c) => !double.IsNaN(rows ? c.Height : c.Width);
        return kids.FirstOrDefault(c => GetFoldedSize(c) > 0 && Sized(c)) ?? kids.FirstOrDefault(c => c.IsVisible && Sized(c)) ?? kids.FirstOrDefault(Sized);
    }

    private static double Size(Grid g, bool rows, int k) =>
        SizedChild(g, rows, k) is { } c ? (c.IsVisible ? (rows ? c.Height : c.Width) : 0) : Actual(g, rows, k);

    private static string Describe(Grid g, bool rows, int k) =>
        SizedChild(g, rows, k) is { } c ? "child:" + (rows ? c.Height : c.Width).ToString("0", CultureInfo.InvariantCulture) : Def(g, rows, k).ToString();

    private static void Apply(Grid g, bool rows, int k, string v)
    {
        if (v.StartsWith("child:", StringComparison.Ordinal))
        {
            if (SizedChild(g, rows, k) is { } c && double.TryParse(v[6..], NumberStyles.Float, CultureInfo.InvariantCulture, out var px))
            {
                if (rows) { c.MinHeight = Math.Min(c.MinHeight, px); c.Height = px; } else { c.MinWidth = Math.Min(c.MinWidth, px); c.Width = px; }
            }
            return;
        }
        try { SetDef(g, rows, k, GridLength.Parse(v)); } catch (FormatException) { }
    }

    private static void SetPixels(Grid g, bool rows, int k, double px)
    {
        if (SizedChild(g, rows, k) is { } c)
        {
            if (rows) { c.MinHeight = Math.Min(c.MinHeight, px); c.Height = px; } else { c.MinWidth = Math.Min(c.MinWidth, px); c.Width = px; }
        }
        else SetDef(g, rows, k, new GridLength(px, GridUnitType.Pixel));
    }

    private static void Remember(Grid g, State st)
    {
        if (Store == null) return;
        Store[st.Key] = string.Join("|", Enumerable.Range(0, Count(g, st.Rows)).Select(k => Describe(g, st.Rows, k)));
        Save?.Invoke();
    }

    // ---------------------------------------------------------------- install
    private static void Install(Grid g, bool rows)
    {
        if (States.TryGetValue(g, out _)) return;
        var n = Count(g, rows);
        if (n < 2) return;
        var st = new State { Rows = rows };
        States.Add(g, st);
        var host = g.FindAncestorOfType<UserControl>() as Control ?? g.FindAncestorOfType<Window>();
        var index = host?.GetVisualDescendants().OfType<Grid>().Where(x => GetSplit(x) || GetSplitRows(x)).ToList().IndexOf(g) ?? 0;
        st.Key = (GetKey(g) ?? $"{host?.GetType().Name ?? "grid"}.{g.Name ?? index.ToString(CultureInfo.InvariantCulture)}") + (rows ? ".rows" : ".cols");
        st.Defaults = Enumerable.Range(0, n).Select(k => Describe(g, rows, k)).ToArray();
        if (Store != null && Store.TryGetValue(st.Key, out var saved) && saved.Split('|') is { } parts && parts.Length == n)
            for (int k = 0; k < n; ++k) Apply(g, rows, k, parts[k]);
        for (int k = 0; k + 1 < n; ++k)
        {
            var left = Panes_(g).Where(c => Index(c, rows) + Span(c, rows) - 1 == k && !GetFixed(c)).ToList();
            var right = Panes_(g).Where(c => Index(c, rows) == k + 1 && !GetFixed(c)).ToList();
            if (left.Count == 0 || right.Count == 0) continue;
            // along the other axis: where panes meet on both sides
            int From(Control c) => rows ? Grid.GetColumn(c) : Grid.GetRow(c);
            int To(Control c) => From(c) + Math.Max(1, rows ? Grid.GetColumnSpan(c) : Grid.GetRowSpan(c));
            var a = Math.Max(left.Min(From), right.Min(From));
            var b = Math.Min(left.Max(To), right.Max(To));
            if (b <= a) continue;
            var d = new PaneDivider(g, st, k, left, right);
            if (rows) { Grid.SetRow(d, k + 1); Grid.SetColumn(d, a); Grid.SetColumnSpan(d, b - a); }
            else { Grid.SetColumn(d, k + 1); Grid.SetRow(d, a); Grid.SetRowSpan(d, b - a); }
            st.Dividers.Add(d);
            g.Children.Add(d);
        }
    }

    // ---------------------------------------------------------------- hide and show
    public static bool IsHidden(Grid g, int k) => States.TryGetValue(g, out var st) && st.Before.ContainsKey(k);

    /// <summary>Hides pane k of a split grid (its size kept to bring it back).</summary>
    public static void Hide(Grid g, int k, bool remember = true)
    {
        if (!States.TryGetValue(g, out var st) || st.Before.ContainsKey(k) || k >= Count(g, st.Rows)) return;
        st.Before[k] = Describe(g, st.Rows, k);
        var wasStar = Def(g, st.Rows, k).IsStar;
        SetPixels(g, st.Rows, k, SizedChild(g, st.Rows, k) is { } fold ? GetFoldedSize(fold) : 0);
        // a proportional pane (the 3D view) hidden: a set-size neighbour (the dock) takes its space
        if (wasStar)
            foreach (var j in new[] { k + 1, k - 1 })
                if (j >= 0 && j < Count(g, st.Rows) && !st.Before.ContainsKey(j) && SizedChild(g, st.Rows, j) is { } c)
                {
                    st.Filled[k] = (j, c, st.Rows ? c.Height : c.Width);
                    SetDef(g, st.Rows, j, new GridLength(1, GridUnitType.Star));
                    if (st.Rows) c.Height = double.NaN; else c.Width = double.NaN;
                    break;
                }
        foreach (var d in st.Dividers) d.Refresh();
        if (remember) Remember(g, st);
        Toggled?.Invoke(st.Key, k, false);
    }

    /// <summary>Shows pane k again at the size it had (or its first size).</summary>
    public static void Show(Grid g, int k, bool remember = true)
    {
        if (!States.TryGetValue(g, out var st) || k >= Count(g, st.Rows)) return;
        if (st.Filled.Remove(k, out var f))
        {
            SetDef(g, st.Rows, f.J, GridLength.Auto);
            if (st.Rows) f.C.Height = f.Size; else f.C.Width = f.Size;
        }
        var v = st.Before.Remove(k, out var b) ? b : null;
        var folded = SizedChild(g, st.Rows, k) is { } fc ? GetFoldedSize(fc) : 0;
        if (v == null && Size(g, st.Rows, k) > folded + 1) return;   // already shown
        v ??= st.Defaults[k];
        if (v == "0" || v == "child:0") v = st.Defaults[k];
        Apply(g, st.Rows, k, v);
        foreach (var d in st.Dividers) d.Refresh();
        if (remember) Remember(g, st);
        Toggled?.Invoke(st.Key, k, true);
    }

    public static void Toggle(Grid g, int k)
    {
        if (!States.TryGetValue(g, out var st)) return;
        var folded = SizedChild(g, st.Rows, k) is { } fc ? GetFoldedSize(fc) : 0;
        if (st.Before.ContainsKey(k) || Size(g, st.Rows, k) < folded + 1) Show(g, k); else Hide(g, k);
    }

    /// <summary>Both sides of a divider back to their first sizes.</summary>
    private static void Reset(Grid g, State st, int k)
    {
        foreach (var i in new[] { k, k + 1 }) { st.Before.Remove(i); Apply(g, st.Rows, i, st.Defaults[i]); }
        foreach (var d in st.Dividers) d.Refresh();
        Remember(g, st);
    }

    /// <summary>The screenshots: the n-th divider shown under root as if the pointer were on it.</summary>
    internal static void HoverForShot(Visual root, int n) =>
        root.GetVisualDescendants().OfType<PaneDivider>().Where(d => d.IsEffectivelyVisible).Skip(n).FirstOrDefault()?.ShowHover();

    // ---------------------------------------------------------------- the divider
    private sealed class PaneDivider : Panel
    {
        private readonly Grid _g;
        private readonly State _st;
        private readonly int _k;
        private readonly Border _line = new();
        private readonly Border _pill;
        private readonly Button _back, _forward;
        private bool _dragging;
        private double _a, _b;   // the two panes' sizes at the start of the drag
        private double _moved;

        public PaneDivider(Grid g, State st, int k, IEnumerable<Control> left, IEnumerable<Control> right)
        {
            _g = g; _st = st; _k = k;
            var rows = st.Rows;
            ZIndex = 60;
            ClipToBounds = false;
            Background = Brushes.Transparent;
            Cursor = new Cursor(rows ? StandardCursorType.SizeNorthSouth : StandardCursorType.SizeWestEast);
            if (rows) { Height = 8; Margin = new Thickness(0, -4, 0, 0); VerticalAlignment = VerticalAlignment.Top; }
            else { Width = 8; Margin = new Thickness(-4, 0, 0, 0); HorizontalAlignment = HorizontalAlignment.Left; }

            _line.Background = Brushes.Transparent;
            if (rows) { _line.Height = 2; _line.VerticalAlignment = VerticalAlignment.Center; }
            else { _line.Width = 2; _line.HorizontalAlignment = HorizontalAlignment.Center; }
            _line.IsHitTestVisible = false;
            Children.Add(_line);

            var thumb = new Thumb { Cursor = Cursor, Template = new Avalonia.Controls.Templates.FuncControlTemplate<Thumb>((_, _) => new Border { Background = Brushes.Transparent }) };
            thumb.DragStarted += (_, _) => { _dragging = true; _moved = 0; _a = Size(_g, rows, _k); _b = Size(_g, rows, _k + 1); };
            thumb.DragDelta += (_, e) => Drag(rows ? e.Vector.Y : e.Vector.X);
            thumb.DragCompleted += (_, _) => { _dragging = false; Snap(); Remember(_g, _st); Refresh(); };
            thumb.DoubleTapped += (_, _) => Reset(_g, _st, _k);
            Children.Add(thumb);

            _back = Arrow(rows ? "chevu" : "chevl");
            _forward = Arrow(rows ? "chev" : "chevr");
            _back.Click += (_, _) => { if (IsHidden(_g, _k + 1) || Size(_g, rows, _k + 1) < 1) Show(_g, _k + 1); else Hide(_g, _k); };
            _forward.Click += (_, _) => { if (IsHidden(_g, _k) || Size(_g, rows, _k) < 1) Show(_g, _k); else Hide(_g, _k + 1); };
            var stack = new StackPanel { Orientation = rows ? Orientation.Horizontal : Orientation.Vertical, Spacing = 2, HorizontalAlignment = HorizontalAlignment.Center, VerticalAlignment = VerticalAlignment.Center };
            stack.Children.Add(_back);
            stack.Children.Add(_forward);
            _pill = new Border
            {
                Child = stack, Padding = new Thickness(2), CornerRadius = new CornerRadius(8), BorderThickness = new Thickness(1),
                HorizontalAlignment = HorizontalAlignment.Center, VerticalAlignment = VerticalAlignment.Center, IsVisible = false,
                // wider than the 8 px divider: it reaches over both panes' edges, centred on the line
                Margin = rows ? new Thickness(0, -9, 0, -9) : new Thickness(-9, 0, -9, 0),
                Width = rows ? double.NaN : 26, Height = rows ? 26 : double.NaN,
                BoxShadow = BoxShadows.Parse("0 4 14 0 #50000000"),
            };
            _pill.Bind(Border.BackgroundProperty, _pill.GetResourceObservable("Bg2B"));
            _pill.Bind(Border.BorderBrushProperty, _pill.GetResourceObservable("AccB"));
            Children.Add(_pill);

            PointerEntered += (_, _) => Hover(true);
            PointerExited += (_, _) => { if (!_dragging) Hover(false); };
            // no divider beside a pane that is not shown (a dock with nothing open)
            foreach (var c in left.Concat(right)) c.PropertyChanged += (_, e) => { if (e.Property == IsVisibleProperty) Refresh(); };
            _left = left.ToList();
            _right = right.ToList();
            Refresh();
        }
        private readonly List<Control> _left, _right;

        private static Button Arrow(string icon)
        {
            var b = new Button { Width = 20, Height = 20, MinWidth = 20, MinHeight = 20, Padding = new Thickness(0), Cursor = new Cursor(StandardCursorType.Hand) };
            b.Classes.Add("tool");
            b.Content = new Icon { Kind = icon, Size = 12, StrokeWidth = 2, HorizontalAlignment = HorizontalAlignment.Center, VerticalAlignment = VerticalAlignment.Center };
            return b;
        }

        public void ShowHover() => Hover(true);
        private void Hover(bool on)
        {
            _pill.IsVisible = on;
            _line.Background = on && this.TryFindResource("AccB", ActualThemeVariant, out var acc) && acc is IBrush br ? br : Brushes.Transparent;
            if (on) Refresh();
        }

        public void Refresh()
        {
            IsVisible = _left.Any(c => c.IsVisible) && _right.Any(c => c.IsVisible);
            var rows = _st.Rows;
            var aHidden = IsHidden(_g, _k) || Size(_g, rows, _k) < 1;
            var bHidden = IsHidden(_g, _k + 1) || Size(_g, rows, _k + 1) < 1;
            string side(bool first) => rows ? (first ? "upper" : "lower") : (first ? "left" : "right");
            ToolTip.SetTip(_back, bHidden ? $"Show the {side(false)} pane" : $"Hide the {side(true)} pane");
            ToolTip.SetTip(_forward, aHidden ? $"Show the {side(true)} pane" : $"Hide the {side(false)} pane");
            ToolTip.SetTip(this, "Drag to resize · double-click for the first sizes");
        }

        private void Drag(double delta)
        {
            var rows = _st.Rows;
            _moved += delta;
            var total = _a + _b;
            var a = Math.Clamp(_a + _moved, 0, total);
            var b = total - a;
            var da = Def(_g, rows, _k);
            var db = Def(_g, rows, _k + 1);
            if (da.IsStar && db.IsStar)
            {
                // two proportional panes: every proportional pane in pixels as its star value, then these two move
                for (int i = 0; i < Count(_g, rows); ++i)
                    if (Def(_g, rows, i).IsStar) SetDef(_g, rows, i, new GridLength(Math.Max(1, Actual(_g, rows, i)), GridUnitType.Star));
                SetDef(_g, rows, _k, new GridLength(Math.Max(1, a), GridUnitType.Star));
                SetDef(_g, rows, _k + 1, new GridLength(Math.Max(1, b), GridUnitType.Star));
            }
            else
            {
                if (!da.IsStar) SetPixels(_g, rows, _k, a);
                if (!db.IsStar) SetPixels(_g, rows, _k + 1, b);
            }
            _st.Before.Remove(_k);
            _st.Before.Remove(_k + 1);
        }

        /// <summary>A pane dragged thinner than 48 px is hidden (its first size comes back when shown).</summary>
        private void Snap()
        {
            var rows = _st.Rows;
            foreach (var i in new[] { _k, _k + 1 })
            {
                var s = Size(_g, rows, i);
                if (s > 0.5 && s < 48 && !Def(_g, rows, i).IsStar) { SetPixels(_g, rows, i, 0); _st.Before[i] = _st.Defaults[i]; Toggled?.Invoke(_st.Key, i, false); }
            }
        }
    }
}
