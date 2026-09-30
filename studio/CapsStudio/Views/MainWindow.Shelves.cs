using Avalonia;
using Avalonia.Controls;
using Avalonia.Controls.Primitives;
using Avalonia.Input;
using Avalonia.Layout;
using Avalonia.Media;
using Avalonia.Threading;
using Avalonia.VisualTree;
using CapsStudio.ViewModels;

namespace CapsStudio.Views;

/// <summary>Tool shelves (design/boards/Shelves, ShelfEditor): the toolbar rows as shelves that dock in the rows above the
/// view, on its left, right or bottom edge, float over the window, pop out into a window of their own, or fold into a
/// puck; workspaces (Sketch, Build, Analyse, Present, Mine) each keep a layout; shelves the user makes hold proxies of
/// the built-in tools. Layouts are saved in the settings.</summary>
public partial class MainWindow
{
    public static readonly string[] Workspaces = ["Sketch", "Assemble", "Analyse", "Present", "Mine"];
    private readonly Dictionary<string, Shelf> _shelves = new();
    private readonly Dictionary<string, ShelfState> _shelfState = new();
    private readonly Dictionary<string, Puck> _pucks = new();
    private readonly Dictionary<string, Window> _shelfWindows = new();
    private Dictionary<string, Panel> _slots = new();
    private Popup? _peek;
    private Border? _guide;
    // a drag in progress: the shelf, the pointer's offset in it, where it started, whether it has moved
    private Shelf? _dragShelf;
    private Point _dragOffset, _dragStart;
    private bool _dragMoved, _dragBlocked;
    private string? _dragTarget;

    public IReadOnlyDictionary<string, Shelf> Shelves => _shelves;
    public ShelfState? ShelfStateOf(string id) => _shelfState.GetValueOrDefault(id);

    private void InitShelves()
    {
        _slots = new Dictionary<string, Panel>
        {
            ["top-left"] = SlotTopLeft, ["top-right"] = SlotTopRight, ["top-2"] = SlotTop2,
            ["left"] = SlotLeft, ["right"] = SlotRight, ["bottom"] = SlotBottom,
        };
        ((Panel)ToolbarLeft.Parent!).Children.Remove(ToolbarLeft);
        ((Panel)ToolbarRight.Parent!).Children.Remove(ToolbarRight);
        SlotTop2.Children.Remove(ModifyRow);
        Add(new Shelf("tools", "Sketch & edit", "cursor", ToolbarLeft));
        Add(new Shelf("view", "View & panels", "eye", ToolbarRight));
        var modify = new Shelf("modify", "Modify", "atom", ModifyRow);
        modify.Bind(IsVisibleProperty, new Avalonia.Data.Binding(nameof(MainViewModel.IsStudio)));
        Add(modify);
        foreach (var c in _vm.Settings.CustomShelves) AddCustom(c);
        _guide = new Border { Classes = { "shelfguide" }, IsVisible = false, IsHitTestVisible = false };
        ShelfFloat.Children.Add(_guide);
        ShelfFloat.PointerMoved += OnShelfDragMove;
        ShelfFloat.PointerReleased += OnShelfDragEnd;
        ShelfFloat.PointerCaptureLost += (_, _) => EndShelfDrag(commit: true);
        foreach (var w in Workspaces)
        {
            var b = new Button { Classes = { "ws" }, Content = w, Tag = w, [ToolTip.TipProperty] = $"{w} workspace: its own shelves and places" };
            b.Click += (_, _) => UseWorkspace(w);
            WorkspaceButtons.Children.Add(b);
        }
        _vm.PropertyChanged += (_, e) => { if (e.PropertyName == nameof(MainViewModel.IsStudio)) UpdateModifyBar(); };
        _vm.SettingsReplaced += () => { Panes.Store = _vm.Settings.PaneSizes; ReloadShelves(); };
        // the first structure opened after shelves arrived: say where they are, once
        _vm.PropertyChanged += (_, e) =>
        {
            if (e.PropertyName != nameof(MainViewModel.HasDocument) || !_vm.HasDocument || _vm.Settings.ShelfHintShown) return;
            _vm.Settings.ShelfHintShown = true;
            _vm.Settings.Save();
            DispatcherTimer.RunOnce(() => _vm.Status = "Tool shelves: drag the small tab at the left of a toolbar to move it (dock on any edge of the view or float), click it to fold · workspaces at the top", TimeSpan.FromSeconds(1.5));
        };
        _vm.AddCommand(new PaletteCommand { Title = "Make a shelf…", Id = "shelf.make", Icon = "pin", Section = "View", Keywords = "toolbar custom tools", Run = () => _ = MakeShelfDialog() });
        _vm.AddCommand(new PaletteCommand { Title = "Reset this workspace's shelves", Id = "shelf.reset", Icon = "layers", Section = "View", Keywords = "toolbar layout default", Run = ResetWorkspace });
        foreach (var (sid, title) in new[] { ("tools", "Sketch & edit"), ("view", "View & panels"), ("modify", "Modify") })
        {
            _vm.AddCommand(new PaletteCommand { Title = $"Shelf {title}: pop out into a window", Id = $"shelf.{sid}.window", Icon = "expand", Section = "View", Keywords = "toolbar float detach", Run = () => MoveShelf(sid, "window") });
            _vm.AddCommand(new PaletteCommand { Title = $"Shelf {title}: float over the view", Id = $"shelf.{sid}.float", Icon = "layers", Section = "View", Keywords = "toolbar move", Run = () => MoveShelf(sid, "float") });
            _vm.AddCommand(new PaletteCommand { Title = $"Shelf {title}: fold / unfold", Id = $"shelf.{sid}.fold", Icon = "chev", Section = "View", Keywords = "toolbar collapse puck", Run = () => FoldShelf(sid, !(ShelfStateOf(sid)?.Folded ?? false)) });
            _vm.AddCommand(new PaletteCommand { Title = $"Shelf {title}: show / hide", Id = $"shelf.{sid}.toggle", Icon = "eye", Section = "View", Keywords = "toolbar", Run = () => ToggleShelf(sid) });
        }
        foreach (var ws in Workspaces)
            _vm.AddCommand(new PaletteCommand { Title = $"Workspace: {ws}", Id = "workspace." + ws.ToLowerInvariant(), Icon = "layers", Section = "View", Keywords = "shelves toolbar layout", Run = () => UseWorkspace(ws) });
        // the workspace called Build before (it shared its name with the Build module) is Assemble
        if (_vm.Settings.Workspace == "Build") _vm.Settings.Workspace = "Assemble";
        if (_vm.Settings.ShelfLayouts.Remove("Build", out var old)) _vm.Settings.ShelfLayouts["Assemble"] = old;
        UseWorkspace(Workspaces.Contains(_vm.Settings.Workspace) ? _vm.Settings.Workspace : "Sketch", save: false);
    }

    /// <summary>The settings were reset or imported: the shelves the user made follow them, and the workspace's layout.</summary>
    private void ReloadShelves()
    {
        foreach (var id in _shelves.Keys.Where(k => k.StartsWith("my-", StringComparison.Ordinal) && _vm.Settings.CustomShelves.All(c => c.Id != k)).ToList())
        {
            Detach(_shelves[id]);
            if (_pucks.Remove(id, out var p)) Detach(p);
            _shelves.Remove(id);
        }
        foreach (var c in _vm.Settings.CustomShelves.Where(c => !_shelves.ContainsKey(c.Id))) AddCustom(c);
        UseWorkspace(Workspaces.Contains(_vm.Settings.Workspace) ? _vm.Settings.Workspace : "Sketch", save: false);
    }

    private void Add(Shelf s)
    {
        s.DataContext = _vm;
        _shelves[s.Id] = s;
        s.Tab.AddHandler(PointerPressedEvent, OnShelfTabPressed, Avalonia.Interactivity.RoutingStrategies.Tunnel);
        s.Tab.ContextMenu = ShelfMenu(s);
        // keyboard (Space, Enter) and assistive technology press the tab: it folds, as a mouse click does
        s.Tab.Click += (_, _) => { if (_dragShelf == null) FoldShelf(s.Id, true); };
    }

    private ContextMenu ShelfMenu(Shelf s)
    {
        MenuItem Item(string h, Action a) { var m = new MenuItem { Header = h }; m.Click += (_, _) => a(); return m; }
        var dock = new MenuItem { Header = "Dock" };
        foreach (var (h, d) in new[] { ("In the rows above the view", "top-left"), ("Second row", "top-2"), ("Left of the view", "left"), ("Right of the view", "right"), ("Below the view", "bottom") })
            dock.Items.Add(Item(h, () => MoveShelf(s.Id, d)));
        var items = new List<object>
        {
            Item("Fold into a puck", () => FoldShelf(s.Id, true)), Item("Float over the view", () => MoveShelf(s.Id, "float")),
            Item("Pop out into a window", () => MoveShelf(s.Id, "window")), dock, new Separator(),
            Item("Hide", () => MoveShelf(s.Id, "hidden")), Item("Reset this workspace", ResetWorkspace),
        };
        if (s.Id.StartsWith("my-", StringComparison.Ordinal)) items.Add(Item("Delete this shelf", () => DeleteCustomShelf(s.Id)));
        return new ContextMenu { ItemsSource = items };
    }

    // ---------------------------------------------------------------- layouts
    private static List<ShelfState> DefaultLayout(string ws)
    {
        ShelfState S(string id, string dock, int order, bool folded = false) => new() { Id = id, Dock = dock, Order = order, Folded = folded };
        return ws switch
        {
            "Assemble" => [S("tools", "top-left", 0), S("view", "top-right", 0), S("modify", "top-2", 0, true)],
            "Analyse" => [S("tools", "top-left", 0, true), S("view", "top-right", 0), S("modify", "hidden", 0)],
            "Present" => [S("tools", "top-left", 0, true), S("view", "top-right", 0, true), S("modify", "hidden", 0)],
            _ => [S("tools", "top-left", 0), S("view", "top-right", 0), S("modify", "top-2", 0)],
        };
    }

    public string Workspace => _vm.Settings.Workspace;

    public void UseWorkspace(string ws, bool save = true)
    {
        if (save) SaveShelfLayout();
        _vm.Settings.Workspace = ws;
        var layout = _vm.Settings.ShelfLayouts.TryGetValue(ws, out var l) && l.Count > 0 ? l : DefaultLayout(ws == "Mine" ? "Sketch" : ws);
        ApplyShelfLayout(layout);
        foreach (var b in WorkspaceButtons.Children.OfType<Button>()) b.Classes.Set("on", (string?)b.Tag == ws);
        if (save) { SaveShelfLayout(); _vm.Status = $"{ws} workspace"; }
    }

    public void ResetWorkspace()
    {
        _vm.Settings.ShelfLayouts.Remove(Workspace);
        ApplyShelfLayout(DefaultLayout(Workspace == "Mine" ? "Sketch" : Workspace));
        SaveShelfLayout();
        _vm.Status = $"{Workspace} workspace reset";
    }

    private void ApplyShelfLayout(List<ShelfState> layout)
    {
        foreach (var s in _shelves.Values) Detach(s);
        foreach (var p in _pucks.Values) Detach(p);
        _pucks.Clear();
        foreach (var w in _shelfWindows.Values.ToList()) { w.Tag = "closing"; w.Close(); }
        _shelfWindows.Clear();
        _shelfState.Clear();
        foreach (var id in _shelves.Keys)
        {
            var st = layout.FirstOrDefault(x => x.Id == id) ?? new ShelfState { Id = id, Dock = id.StartsWith("my-", StringComparison.Ordinal) ? "float" : "hidden", X = 360, Y = 150 };
            _shelfState[id] = new ShelfState { Id = id, Dock = st.Dock, Order = st.Order, X = st.X, Y = st.Y, Folded = st.Folded };
        }
        foreach (var st in _shelfState.Values.OrderBy(x => x.Order)) Place(st.Id);
        UpdateModifyBar();
        Dispatcher.UIThread.Post(FitToolRow, DispatcherPriority.Background);
    }

    /// <summary>The layout as it is now, kept for the workspace in the settings.</summary>
    public void SaveShelfLayout()
    {
        foreach (var (id, st) in _shelfState)
        {
            Control shown = st.Folded && _pucks.TryGetValue(id, out var p) ? p : _shelves[id];
            if (_slots.TryGetValue(st.Dock, out var slot)) st.Order = Math.Max(0, slot.Children.IndexOf(shown));
        }
        _vm.Settings.ShelfLayouts[Workspace] = _shelfState.Values.Select(s => new ShelfState { Id = s.Id, Dock = s.Dock, Order = s.Order, X = s.X, Y = s.Y, Folded = s.Folded }).ToList();
        _vm.Settings.Save();
    }

    private static void Detach(Control c)
    {
        switch (c.Parent)
        {
            case Panel p: p.Children.Remove(c); break;
            case Popup pp: pp.Child = null; break;
            case ContentControl cc: cc.Content = null; break;
            case Decorator d: d.Child = null; break;
        }
    }

    /// <summary>Puts a shelf (or its puck) where its state says.</summary>
    private void Place(string id, int? index = null)
    {
        var s = _shelves[id];
        var st = _shelfState[id];
        Detach(s);
        if (_pucks.Remove(id, out var old)) Detach(old);
        if (_shelfWindows.Remove(id, out var win)) { win.Tag = "closing"; win.Content = null; win.Close(); }
        if (st.Dock == "hidden") { UpdateModifyBar(); return; }
        var vertical = st.Dock is "left" or "right";
        s.Floating = st.Dock is "float" or "window";
        s.Vertical = vertical;
        if (st.Dock == "window") { PopOut(s, st); return; }
        Control shown = s;
        if (st.Folded)
        {
            var count = CountTools(s);
            var puck = new Puck(s, count) { DataContext = _vm };
            puck.Click += (_, _) => FoldShelf(id, false);
            puck.PointerEntered += (_, _) => PeekShelf(id, puck);
            _pucks[id] = puck;
            shown = puck;
        }
        s.MaxWidth = double.PositiveInfinity;
        if (st.Dock == "float")
        {
            var (x, y) = ClampToWindow(st.X, st.Y);
            Canvas.SetLeft(shown, x);
            Canvas.SetTop(shown, y);
            // a floating shelf wraps rather than run past the window's edge
            if (ShelfFloat.Bounds.Width > 0) s.MaxWidth = Math.Max(260, ShelfFloat.Bounds.Width - x - 12);
            ShelfFloat.Children.Add(shown);
        }
        else if (_slots.TryGetValue(st.Dock, out var slot))
        {
            var k = Math.Clamp(index ?? st.Order, 0, slot.Children.Count);
            slot.Children.Insert(k, shown);
        }
        UpdateModifyBar();
        Dispatcher.UIThread.Post(FitToolRow, DispatcherPriority.Background);   // the rows fit again (labels, overflow)
    }

    private int CountTools(Shelf s)
    {
        var n = ShelfTools.Of(s.Tools).Count(c => c.IsVisible);
        return n > 0 ? n : s.Tools.Children.Count;
    }

    private (double, double) ClampToWindow(double x, double y)
    {
        var w = Math.Max(200, ShelfFloat.Bounds.Width);
        var h = Math.Max(200, ShelfFloat.Bounds.Height);
        return (Math.Clamp(x, 0, w - 60), Math.Clamp(y, 44, h - 60));
    }

    public void MoveShelf(string id, string dock, int? index = null)
    {
        if (!_shelfState.TryGetValue(id, out var st)) return;
        if (dock == "float" && st.Dock != "float" && st.X == 0 && st.Y == 0) (st.X, st.Y) = (360, 160);
        st.Dock = dock;
        if (dock == "hidden") st.Folded = false;
        Place(id, index);
        SaveShelfLayout();
        Dispatcher.UIThread.Post(FitToolRow, DispatcherPriority.Background);
    }

    public void FoldShelf(string id, bool fold)
    {
        if (!_shelfState.TryGetValue(id, out var st)) return;
        ClosePeek();
        st.Folded = fold;
        if (st.Dock is "window" or "hidden") st.Dock = "top-2";
        Place(id);
        SaveShelfLayout();
        Dispatcher.UIThread.Post(FitToolRow, DispatcherPriority.Background);
    }

    public void ShowShelf(string id, bool show)
    {
        if (!_shelfState.TryGetValue(id, out var st)) return;
        if (show && st.Dock == "hidden") MoveShelf(id, DefaultLayout("Sketch").FirstOrDefault(x => x.Id == id)?.Dock ?? "float");
        else if (!show) MoveShelf(id, "hidden");
    }

    /// <summary>The Modify row keeps its line only while a shelf sits in it (and the Studio shows).</summary>
    private void UpdateModifyBar() => ModifyBar.IsVisible = _vm.IsStudio && SlotTop2.Children.OfType<Control>().Any(c => c.IsVisible || c is Puck);

    // ---------------------------------------------------------------- peek and pop-out
    private void PeekShelf(string id, Puck puck)
    {
        if (_dragShelf != null) return;
        ClosePeek();
        var s = _shelves[id];
        var st = _shelfState[id];
        Detach(s);
        s.Floating = true;
        s.Vertical = st.Dock is "left" or "right";
        _peek = new Popup
        {
            PlacementTarget = puck, Child = s, IsLightDismissEnabled = false,   // it closes when the pointer leaves the puck and the shelf, so a click on the puck still unfolds
            Placement = st.Dock switch { "left" => PlacementMode.Right, "right" => PlacementMode.Left, "bottom" => PlacementMode.Top, _ => PlacementMode.Bottom },
        };
        ((ISetLogicalParent)_peek).SetParent(this);
        _peek.Closed += (_, _) => { if (_peek?.Child == s) _peek.Child = null; };
        s.PointerExited += PeekLeft;
        puck.PointerExited += PeekLeft;
        _peek.Open();
        void PeekLeft(object? o, PointerEventArgs e) =>
            DispatcherTimer.RunOnce(() =>
            {
                if (s.IsPointerOver || puck.IsPointerOver || _peek?.Child != s) return;
                s.PointerExited -= PeekLeft;
                puck.PointerExited -= PeekLeft;
                ClosePeek();
            }, TimeSpan.FromMilliseconds(300));
    }

    private void ClosePeek()
    {
        if (_peek == null) return;
        var p = _peek;
        _peek = null;
        p.Close();
        p.Child = null;
        ((ISetLogicalParent)p).SetParent(null);
    }

    private void PopOut(Shelf s, ShelfState st)
    {
        s.Vertical = false;
        var w = new Window
        {
            Title = s.Title, SizeToContent = SizeToContent.WidthAndHeight, CanResize = false, ShowInTaskbar = false, Content = s, DataContext = _vm,
            Background = (IBrush?)this.FindResource("Bg1B") ?? Brushes.Black, Padding = new Thickness(6), Topmost = true,
        };
        if (st.X != 0 || st.Y != 0) { w.WindowStartupLocation = WindowStartupLocation.Manual; w.Position = new PixelPoint((int)st.X, (int)st.Y); }
        else w.WindowStartupLocation = WindowStartupLocation.CenterOwner;
        w.PositionChanged += (_, e) => { st.X = e.Point.X; st.Y = e.Point.Y; };
        w.Closing += (_, _) =>
        {
            if ((string?)w.Tag == "closing") return;   // closed by a layout change
            w.Content = null;
            _shelfWindows.Remove(s.Id);
            st.Dock = DefaultLayout("Sketch").FirstOrDefault(x => x.Id == s.Id)?.Dock ?? "float";
            Place(s.Id);
            SaveShelfLayout();
        };
        _shelfWindows[s.Id] = w;
        w.Show(this);
    }

    // ---------------------------------------------------------------- drag and dock
    private void OnShelfTabPressed(object? sender, PointerPressedEventArgs e)
    {
        if (sender is not Button tab || !e.GetCurrentPoint(tab).Properties.IsLeftButtonPressed) return;
        var s = _shelves.Values.FirstOrDefault(x => x.Tab == tab);
        if (s == null) return;
        e.Handled = true;
        ClosePeek();
        _dragShelf = s;
        _dragMoved = false;
        _dragBlocked = false;
        _dragStart = e.GetPosition(ShelfFloat);
        var origin = s.TranslatePoint(new Point(0, 0), ShelfFloat) ?? _dragStart;
        _dragOffset = _dragStart - origin;
        e.Pointer.Capture(ShelfFloat);
    }

    private void OnShelfDragMove(object? sender, PointerEventArgs e)
    {
        if (_dragShelf is not { } s) return;
        var p = e.GetPosition(ShelfFloat);
        if (!_dragMoved)
        {
            if (Math.Abs(p.X - _dragStart.X) + Math.Abs(p.Y - _dragStart.Y) < 5) return;
            if (_vm.Settings.ShelvesLocked) { _dragBlocked = true; _vm.Status = "Shelves are locked (View › Shelves › Lock shelves in place)"; return; }
            _dragMoved = true;
            if (_shelfWindows.Remove(s.Id, out var w)) { w.Tag = "closing"; w.Content = null; w.Close(); }
            Detach(s);
            s.Floating = true;
            s.Vertical = false;
            ShelfFloat.Children.Add(s);
        }
        Canvas.SetLeft(s, p.X - _dragOffset.X);
        Canvas.SetTop(s, p.Y - _dragOffset.Y);
        _dragTarget = DockAt(p, out var rect);
        if (_guide != null)
        {
            _guide.IsVisible = _dragTarget != null;
            if (_dragTarget != null)
            {
                Canvas.SetLeft(_guide, rect.X); Canvas.SetTop(_guide, rect.Y);
                _guide.Width = rect.Width; _guide.Height = rect.Height;
                ToolTip.SetTip(_guide, null);
            }
            _guide.ZIndex = -1;
        }
    }

    private void OnShelfDragEnd(object? sender, PointerReleasedEventArgs e)
    {
        if (_dragShelf == null) return;
        e.Pointer.Capture(null);
        EndShelfDrag(commit: true, e.GetPosition(ShelfFloat));
    }

    private void EndShelfDrag(bool commit, Point? at = null)
    {
        if (_dragShelf is not { } s) return;
        _dragShelf = null;
        if (_guide != null) _guide.IsVisible = false;
        var st = _shelfState[s.Id];
        if (!_dragMoved)
        {
            if (!_dragBlocked) FoldShelf(s.Id, true);   // a click on the tab folds the shelf (a drag on a locked one does nothing)
            return;
        }
        var target = _dragTarget;
        _dragTarget = null;
        if (target != null && _slots.TryGetValue(target, out var slot))
        {
            var index = InsertIndex(slot, at ?? default);
            st.Dock = target;
            st.Folded = false;
            Place(s.Id, index);
        }
        else
        {
            st.Dock = "float";
            st.X = Canvas.GetLeft(s);
            st.Y = Canvas.GetTop(s);
            Place(s.Id);
        }
        SaveShelfLayout();
        Dispatcher.UIThread.Post(FitToolRow, DispatcherPriority.Background);
    }

    /// <summary>Which dock the pointer is over: the rows above the view, or near an edge of the view (and its guide).</summary>
    private string? DockAt(Point p, out Rect guide)
    {
        guide = default;
        Rect R(Control c) => c.TranslatePoint(new Point(0, 0), ShelfFloat) is { } o ? new Rect(o, c.Bounds.Size) : default;
        var row = R(ToolRow);
        var mod = ModifyBar.IsVisible ? R(ModifyBar) : new Rect(row.X, row.Bottom, row.Width, 34);
        var view = R(ViewDock);
        const double edge = 56;
        if (row.Contains(p)) { guide = new Rect(row.X, row.Bottom - 3, row.Width, 3); return p.X < row.X + row.Width * 0.62 ? "top-left" : "top-right"; }
        if (mod.Contains(p) || (p.Y >= row.Bottom && p.Y < view.Y && p.X >= row.X)) { guide = new Rect(mod.X, mod.Bottom - 3, mod.Width, 3); return "top-2"; }
        if (!view.Contains(p)) return null;
        if (p.X < view.X + edge) { guide = new Rect(view.X, view.Y, 4, view.Height); return "left"; }
        if (p.X > view.Right - edge) { guide = new Rect(view.Right - 4, view.Y, 4, view.Height); return "right"; }
        if (p.Y > view.Bottom - edge) { guide = new Rect(view.X, view.Bottom - 4, view.Width, 4); return "bottom"; }
        return null;
    }

    private int InsertIndex(Panel slot, Point p)
    {
        var k = 0;
        foreach (var c in slot.Children)
        {
            if (c.TranslatePoint(new Point(c.Bounds.Width / 2, c.Bounds.Height / 2), ShelfFloat) is { } mid)
            {
                var vertical = slot is StackPanel { Orientation: Orientation.Vertical };
                if (vertical ? p.Y < mid.Y : p.X < mid.X) return k;
            }
            ++k;
        }
        return slot.Children.Count;
    }

    // ---------------------------------------------------------------- shelves the user makes
    /// <summary>Every tool of the built-in shelves, by id, for Make a shelf.</summary>
    public List<(string Id, string Label, Control Tool)> AllTools()
    {
        var seen = new HashSet<string>();
        var list = new List<(string, string, Control)>();
        foreach (var id in new[] { "tools", "view", "modify" })
            foreach (var c in ShelfTools.Of(_shelves[id].Tools))
            {
                var tid = ShelfTools.Id(c);
                if (tid.Length > 0 && seen.Add(tid)) list.Add((tid, ShelfTools.Label(c), c));
            }
        return list;
    }

    private void AddCustom(CustomShelf c)
    {
        var tools = AllTools().ToDictionary(t => t.Id, t => t.Tool);
        var panel = new WrapPanel { Orientation = Orientation.Horizontal, ItemSpacing = 2, LineSpacing = 2, MinHeight = 40, VerticalAlignment = VerticalAlignment.Center };
        foreach (var id in c.Tools) if (tools.TryGetValue(id, out var t)) panel.Children.Add(new ProxyTool(id, t));
        Add(new Shelf(c.Id, c.Name, c.Glyph, panel));
    }

    /// <summary>A new shelf with these tools, floating over the view.</summary>
    public string MakeCustomShelf(string name, string glyph, IEnumerable<string> toolIds)
    {
        var id = "my-" + Guid.NewGuid().ToString("N")[..8];
        var c = new CustomShelf { Id = id, Name = name.Trim().Length > 0 ? name.Trim() : "My shelf", Glyph = glyph, Tools = toolIds.ToList() };
        _vm.Settings.CustomShelves.Add(c);
        AddCustom(c);
        _shelfState[id] = new ShelfState { Id = id, Dock = "float", X = 360, Y = 160 };
        Place(id);
        SaveShelfLayout();
        return id;
    }

    public void ToggleShelf(string id)
    {
        if (!_shelfState.TryGetValue(id, out var st)) return;
        ShowShelf(id, st.Dock == "hidden");
    }

    public static readonly string[] ShelfGlyphs = ["pin", "cursor", "atom", "bond", "ring", "ruler", "tag", "eye", "filter", "wand", "layers", "chart", "flask", "cube", "hex"];

    /// <summary>Make a shelf (design/boards/ShelfEditor): click tools to put them on the shelf (click one on the shelf to take
    /// it off), name it, pick its glyph; it appears floating over the view.</summary>
    public async Task MakeShelfDialog()
    {
        var chosen = new List<string>();
        var tools = AllTools();
        var preview = new WrapPanel { Orientation = Orientation.Horizontal, ItemSpacing = 2, LineSpacing = 2, MinHeight = 40 };
        var name = new TextBox { Text = "My shelf", Width = 240 };
        var glyph = new ComboBox { ItemsSource = ShelfGlyphs, SelectedIndex = 0, Width = 120 };
        var hint = new TextBlock { Classes = { "dim" }, Text = "drop tools here: click them above", VerticalAlignment = VerticalAlignment.Center };
        void Refresh()
        {
            preview.Children.Clear();
            var byId = tools.ToDictionary(t => t.Id, t => t.Tool);
            foreach (var id in chosen.ToList())
            {
                var b = new ProxyTool(id, byId[id]);
                b.Click += (_, _) => { chosen.Remove(id); Refresh(); };
                ToolTip.SetTip(b, "Click to take it off the shelf");
                preview.Children.Add(b);
            }
            hint.IsVisible = chosen.Count == 0;
        }
        var grid = new WrapPanel { Orientation = Orientation.Horizontal, ItemSpacing = 4, LineSpacing = 8, MaxWidth = 640 };
        foreach (var (id, label, tool) in tools)
        {
            var face = new StackPanel { Spacing = 4, Width = 74, Children = { new Border { Child = ShelfTools.Face(tool), Height = 22, HorizontalAlignment = HorizontalAlignment.Center },
                new TextBlock { Text = label, FontSize = 10.5, TextWrapping = TextWrapping.Wrap, TextAlignment = TextAlignment.Center, MaxHeight = 44, TextTrimming = TextTrimming.CharacterEllipsis } } };
            var b = new Button { Classes = { "ghost" }, Content = face, Padding = new Thickness(2, 6), [ToolTip.TipProperty] = ToolTip.GetTip(tool) };
            b.Click += (_, _) => { if (!chosen.Contains(id)) { chosen.Add(id); Refresh(); } };
            grid.Children.Add(b);
        }
        var done = new Button { Classes = { "primary" }, Content = "Make the shelf" };
        var cancel = new Button { Content = "Cancel" };
        var dlg = new Window
        {
            Title = "Make a shelf", SizeToContent = SizeToContent.WidthAndHeight, CanResize = false, ShowInTaskbar = false, WindowStartupLocation = WindowStartupLocation.CenterOwner,
            Background = (IBrush?)this.FindResource("Bg1B") ?? Brushes.Black,
            Content = new StackPanel
            {
                Margin = new Thickness(18), Spacing = 12, Width = 660,
                Children =
                {
                    new TextBlock { Text = "Make a shelf", Classes = { "h2" } },
                    new TextBlock { Text = "Click tools to put them on the shelf below; click one on the shelf to take it off. It appears floating over the view: drag its tab to dock it.", TextWrapping = TextWrapping.Wrap, Classes = { "muted" } },
                    new ScrollViewer { MaxHeight = 330, Content = grid },
                    new Border { Classes = { "dropzone" }, Padding = new Thickness(10), Child = new Panel { Children = { hint, preview } } },
                    new StackPanel { Orientation = Orientation.Horizontal, Spacing = 10, Children = { new TextBlock { Text = "Name", VerticalAlignment = VerticalAlignment.Center }, name,
                        new TextBlock { Text = "Glyph", VerticalAlignment = VerticalAlignment.Center, Margin = new Thickness(8, 0, 0, 0) }, glyph } },
                    new StackPanel { Orientation = Orientation.Horizontal, Spacing = 8, HorizontalAlignment = HorizontalAlignment.Right, Children = { cancel, done } },
                },
            },
        };
        var ok = false;
        done.Click += (_, _) => { if (chosen.Count == 0) { hint.Text = "Put at least one tool on the shelf"; return; } ok = true; dlg.Close(); };
        cancel.Click += (_, _) => dlg.Close();
        await dlg.ShowDialog(this);
        if (ok) { MakeCustomShelf(name.Text ?? "", (string?)glyph.SelectedItem ?? "pin", chosen); _vm.Status = $"Shelf '{name.Text}' made · drag its tab to dock it"; }
    }

    public void DeleteCustomShelf(string id)
    {
        if (!_shelves.Remove(id, out var s)) return;
        Detach(s);
        if (_pucks.Remove(id, out var p)) Detach(p);
        _shelfState.Remove(id);
        _vm.Settings.CustomShelves.RemoveAll(c => c.Id == id);
        foreach (var l in _vm.Settings.ShelfLayouts.Values) l.RemoveAll(x => x.Id == id);
        SaveShelfLayout();
    }
}
