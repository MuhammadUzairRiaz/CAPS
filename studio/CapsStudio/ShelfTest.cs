using Avalonia;
using Avalonia.Headless;
using Avalonia.Input;
using Avalonia.Threading;
using Avalonia.VisualTree;
using CapsStudio.Views;

namespace CapsStudio;

/// <summary>--shelftest: the tool shelves in the real main window on the headless platform — pointer drags onto each dock,
/// fold and unfold, workspaces keeping their own layouts, a shelf the user makes, lock, and the layout coming back when
/// the Studio starts again. Exit code 0 when every check passes.</summary>
internal static class ShelfTest
{
    private static int _fail;

    private static void Check(bool ok, string what)
    {
        Console.WriteLine((ok ? "ok   " : "FAIL ") + what);
        if (!ok) ++_fail;
    }

    private static void Pump(int n = 20)
    {
        for (var i = 0; i < n; i++) { Dispatcher.UIThread.RunJobs(); AvaloniaHeadlessPlatform.ForceRenderTimerTick(); Thread.Sleep(15); }
    }

    private static void Drag(MainWindow w, Point a, Point b)
    {
        w.MouseDown(a, MouseButton.Left);
        for (var k = 1; k <= 8; ++k) { w.MouseMove(a + (b - a) * (k / 8.0)); Dispatcher.UIThread.RunJobs(); }
        w.MouseUp(b, MouseButton.Left);
        Pump();
    }

    private static Point TabOf(MainWindow w, string id)
    {
        var tab = w.Shelves[id].Tab;
        return tab.TranslatePoint(new Point(tab.Bounds.Width / 2, tab.Bounds.Height / 2), w) ?? default;
    }

    private static string Where(MainWindow w, string id) => w.ShelfStateOf(id) is { } s ? s.Dock + (s.Folded ? "·folded" : "") : "none";

    public static int Run(string[] args)
    {
        AppBuilder.Configure<App>().UseSkia().UseHeadless(new AvaloniaHeadlessPlatformOptions { UseHeadlessDrawing = false }).WithCapsFonts().SetupWithoutStarting();
        var settings = Path.Combine(Path.GetTempPath(), $"caps-shelftest-{Environment.ProcessId}.json");
        File.Delete(settings);
        ViewModels.AppSettings.Override = settings;
        ViewModels.RecentFiles.Override = Path.Combine(Path.GetTempPath(), "caps-shelftest-recent");
        var file = args.Length > 0 ? args[0] : Path.Combine("samples", "ps_melt.data");

        var w = new MainWindow { Width = 1440, Height = 900 };
        w.Show();
        w.OpenOnStart(file, null);
        Pump(60);
        Check(Where(w, "tools") == "top-left" && Where(w, "view") == "top-right" && Where(w, "modify") == "top-2" && w.Workspace == "Sketch",
              $"Sketch layout: tools {Where(w, "tools")}, view {Where(w, "view")}, modify {Where(w, "modify")}");

        // a drag from the Modify tab to the view's left edge docks it there, vertical
        var t = TabOf(w, "modify");
        Drag(w, t, new Point(t.X + 12, 420));
        Check(Where(w, "modify") == "left" && w.Shelves["modify"].Vertical, $"drag to the view's left edge: {Where(w, "modify")}, vertical {w.Shelves["modify"].Vertical}");
        // to the middle of the view: it floats, a capsule
        t = TabOf(w, "modify");
        Drag(w, t, new Point(700, 430));
        Check(Where(w, "modify") == "float" && w.Shelves["modify"].Floating && !w.Shelves["modify"].Vertical, $"drag into the view: {Where(w, "modify")}");
        // back into the second row
        t = TabOf(w, "modify");
        Drag(w, t, new Point(600, 192));
        Check(Where(w, "modify") == "top-2" && !w.Shelves["modify"].Floating, $"drag back to the second row: {Where(w, "modify")}");
        // the right edge and the bottom
        t = TabOf(w, "view");
        var vd = Avalonia.Controls.NameScopeExtensions.Find<Avalonia.Controls.DockPanel>(Avalonia.Controls.NameScope.GetNameScope(w)!, "ViewDock")!;
        var vr = new Rect(vd.TranslatePoint(new Point(0, 0), w)!.Value, vd.Bounds.Size);
        Drag(w, t, new Point(vr.Right - 20, vr.Y + vr.Height / 2));
        Check(Where(w, "view") == "right", $"view shelf to the right edge: {Where(w, "view")}");
        t = TabOf(w, "tools");
        Drag(w, t, new Point(vr.X + vr.Width / 2 - 100, vr.Bottom - 20));
        Check(Where(w, "tools") == "bottom", $"tools shelf below the view: {Where(w, "tools")}");

        // a click on a tab folds; a click on the puck unfolds
        t = TabOf(w, "modify");
        w.MouseDown(t, MouseButton.Left); w.MouseUp(t, MouseButton.Left); Pump();
        Check(Where(w, "modify") == "top-2·folded", $"click on the tab folds: {Where(w, "modify")}");
        // hovering the puck peeks the shelf beside it (a popup); leaving it closes the peek; clicking the puck unfolds
        var puck = w.GetVisualDescendants().OfType<Puck>().FirstOrDefault(p => p.Shelf.Id == "modify");
        var pc = puck?.TranslatePoint(new Point(puck.Bounds.Width / 2, puck.Bounds.Height / 2), w) ?? default;
        w.MouseMove(pc - new Vector(80, 0)); Pump(3);
        w.MouseMove(pc); Pump();
        Check(puck != null && w.Shelves["modify"].Parent is Avalonia.Controls.Primitives.Popup { IsOpen: true }, $"hover on the puck peeks: {w.Shelves["modify"].Parent?.GetType().Name}");
        w.MouseDown(pc, MouseButton.Left); w.MouseUp(pc, MouseButton.Left); Pump();
        Check(Where(w, "modify") == "top-2" && w.Shelves["modify"].Parent is Avalonia.Controls.Panel, $"click on the puck unfolds: {Where(w, "modify")}");
        // Space / Enter or a screen reader presses the tab: it folds too
        var tab = w.Shelves["modify"].Tab;
        tab.RaiseEvent(new Avalonia.Interactivity.RoutedEventArgs(Avalonia.Controls.Button.ClickEvent, tab)); Pump();
        Check(Where(w, "modify") == "top-2·folded", $"a keyboard press on the tab folds: {Where(w, "modify")}");
        w.FoldShelf("modify", false); Pump();
        Check(Where(w, "modify") == "top-2", $"unfold: {Where(w, "modify")}");

        // workspaces keep their own layouts
        var mine = Where(w, "tools");
        w.UseWorkspace("Present"); Pump();
        Check(Where(w, "tools") == "top-left·folded" && Where(w, "view") == "top-right·folded" && Where(w, "modify") == "hidden", $"Present: tools {Where(w, "tools")}, view {Where(w, "view")}, modify {Where(w, "modify")}");
        w.UseWorkspace("Sketch"); Pump();
        Check(Where(w, "tools") == mine && Where(w, "view") == "right", $"back to Sketch keeps its layout: tools {Where(w, "tools")}, view {Where(w, "view")}");

        // a shelf the user makes: proxies of the built-in tools, floating; clicking a proxy works the original
        var ids = w.AllTools().Select(x => x.Id).ToList();
        var pick = ids.Where(x => x is "lasso_select" or "measure" or "auto_clean").ToList();
        var id = w.MakeCustomShelf("Interface", "pin", pick);
        Pump();
        var proxies = w.Shelves[id].Tools.Children.OfType<ProxyTool>().ToList();
        Check(pick.Count == 3 && proxies.Count == 3 && Where(w, id) == "float", $"made a shelf: {proxies.Count} tools ({string.Join(", ", pick)}) of {ids.Count}, {Where(w, id)}");
        var lasso = proxies.FirstOrDefault(p => p.ToolId == "lasso_select");
        lasso?.RaiseEvent(new Avalonia.Interactivity.RoutedEventArgs(Avalonia.Controls.Button.ClickEvent, lasso));
        Pump();
        Check(w.ViewModel.IsLassoTool, "the shelf's Lasso turns the lasso tool on");

        // locked shelves do not move
        w.ViewModel.Settings.ShelvesLocked = true;
        t = TabOf(w, "modify");
        Drag(w, t, new Point(700, 430));
        Check(Where(w, "modify") == "top-2", $"locked: the drag leaves it {Where(w, "modify")}");
        w.ViewModel.Settings.ShelvesLocked = false;
        w.SaveShelfLayout();

        // the next start reads the file: the same layout comes back
        var saved = ViewModels.AppSettings.Load(settings);
        var sk = saved.ShelfLayouts.GetValueOrDefault("Sketch") ?? new();
        string D(string sid) => sk.FirstOrDefault(x => x.Id == sid) is { } x ? x.Dock + (x.Folded ? "·folded" : "") : "none";
        Check(saved.Workspace == "Sketch" && D("view") == "right" && D("tools") == mine && D(id) == "float" && saved.CustomShelves.Any(c => c.Id == id && c.Tools.Count == 3)
              && saved.ShelfLayouts.GetValueOrDefault("Present")?.FirstOrDefault(x => x.Id == "modify")?.Dock == "hidden",
              $"saved for the next start: view {D("view")}, tools {D("tools")}, my shelf {D(id)} with {saved.CustomShelves.FirstOrDefault(c => c.Id == id)?.Tools.Count} tools");
        w.DeleteCustomShelf(id);
        w.ResetWorkspace();
        Pump();
        Check(!w.Shelves.ContainsKey(id) && Where(w, "view") == "top-right" && Where(w, "tools") == "top-left", $"delete my shelf, reset: view {Where(w, "view")}, tools {Where(w, "tools")}");
        File.Delete(settings);
        Console.WriteLine(_fail == 0 ? "all shelf checks passed" : $"{_fail} shelf checks failed");
        Environment.Exit(_fail == 0 ? 0 : 1);   // no window closing: the main window's quit path is not part of this test
        return _fail == 0 ? 0 : 1;
    }
}
