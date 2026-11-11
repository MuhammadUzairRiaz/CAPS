using Avalonia;
using Avalonia.Controls;
using Avalonia.Controls.Primitives;
using Avalonia.Headless;
using Avalonia.Input;
using Avalonia.Threading;
using Avalonia.VisualTree;
using CapsStudio.ViewModels;
using CapsStudio.Views;

namespace CapsStudio;

/// <summary>--guisweep [SAMPLE.data]: the real main window on the headless platform, a structure open, and on every module
/// page every visible, enabled button, toggle, check box and drop-down clicked with the pointer; every menu command run;
/// every pane divider dragged. A click that throws, leaves an error, or closes the window is reported with its page and
/// its label. Runs a click starts are cancelled at once (the runs themselves are the self-test's). CAPS_SWEEP_MODULES=8,9
/// limits it to those pages. Exit code 0 when nothing failed.</summary>
internal static class GuiSweep
{
    private static readonly List<string> Failures = new();
    private static readonly List<string> Errors = new();
    private static int _clicks, _menus, _dividers;

    // labels that leave the program, change the user's files or settings for good, or reach other machines: not clicked
    private static readonly string[] Skip =
    [
        "quit", "exit", "delete", "remove project", "forget", "discard", "reset settings", "reset all", "uninstall", "submit",
        "cluster", "ssh", "upload", "send report", "report a problem", "check for updates", "install", "download",
        "sign in", "log in", "shut down", "restart", "erase", "move to trash", "overwrite",
    ];

    private static void Pump(int n = 6)
    {
        for (var i = 0; i < n; i++) { Dispatcher.UIThread.RunJobs(); AvaloniaHeadlessPlatform.ForceRenderTimerTick(); Thread.Sleep(8); }
    }

    private static string Describe(Control c)
    {
        string? text = c switch
        {
            ContentControl cc when cc.Content is string s => s,
            ContentControl cc when cc.Content is TextBlock tb => tb.Text,
            ContentControl cc when cc.Content is Control inner => inner.GetVisualDescendants().OfType<TextBlock>().Select(t => t.Text).FirstOrDefault(t => !string.IsNullOrWhiteSpace(t)),
            _ => null,
        };
        var tip = ToolTip.GetTip(c) as string;
        var name = string.IsNullOrEmpty(c.Name) ? "" : "#" + c.Name;
        var label = !string.IsNullOrWhiteSpace(text) ? text : tip ?? c.GetType().Name;
        if (label.Length > 70) label = label[..70] + "…";
        return $"{c.GetType().Name}{name} '{label.Replace('\n', ' ')}'";
    }

    private static bool Dangerous(Control c)
    {
        var all = (Describe(c) + " " + (ToolTip.GetTip(c) as string ?? "")).ToLowerInvariant();
        return Skip.Any(all.Contains);
    }

    private static void StopEverything(MainViewModel vm)
    {
        foreach (var j in vm.Jobs.Where(j => j.IsRunning).ToList()) try { vm.CancelJob(j); } catch { }
        try { vm.CancelAdsorption(); vm.CancelDpd(); vm.CancelExchange(); vm.CancelBench(); vm.StopBatch(); vm.Analyze.Cancel(); vm.CancelExport(); } catch { }
        for (var k = 0; k < 400 && (!vm.Idle || vm.Busy); k++) Pump(1);
    }

    private static void CloseExtras(MainWindow w)
    {
        if (Application.Current?.ApplicationLifetime is Avalonia.Controls.ApplicationLifetimes.IClassicDesktopStyleApplicationLifetime life)
            foreach (var x in life.Windows.Where(x => x != w).ToList()) try { x.Close(); } catch { }
        foreach (var x in w.OwnedWindows.ToList()) try { x.Close(); } catch { }
        // popups, flyouts and overlays: Escape, as the user would
        w.KeyPress(Key.Escape, RawInputModifiers.None, PhysicalKey.Escape, null);
        Pump(3);
    }

    private static void Click(Control c)
    {
        if (TopLevel.GetTopLevel(c) is not { } top) return;
        var p = c.TranslatePoint(new Point(c.Bounds.Width / 2, c.Bounds.Height / 2), top);
        if (p is not { } pt) return;
        top.MouseMove(pt);
        top.MouseDown(pt, MouseButton.Left);
        top.MouseUp(pt, MouseButton.Left);
    }

    public static int Run(string[] args)
    {
        AppBuilder.Configure<App>().UseSkia().UseHeadless(new AvaloniaHeadlessPlatformOptions { UseHeadlessDrawing = false }).WithCapsFonts().SetupWithoutStarting();
        var tmp = Path.Combine(Path.GetTempPath(), $"caps-guisweep-{Environment.ProcessId}");
        Directory.CreateDirectory(tmp);
        AppSettings.Override = Path.Combine(tmp, "settings.json");
        RecentFiles.Override = Path.Combine(tmp, "recent");
        var sample = args.Length > 0 ? args[0] : Path.Combine("samples", "ps_melt.data");
        var only = Environment.GetEnvironmentVariable("CAPS_SWEEP_MODULES") is { Length: > 0 } list ? list.Split(',').Select(int.Parse).ToHashSet() : null;

        var current = "";
        Dispatcher.UIThread.UnhandledException += (_, e) => { Errors.Add($"{current}: {e.Exception.GetType().Name}: {e.Exception.Message}"); e.Handled = true; };
        TaskScheduler.UnobservedTaskException += (_, e) => { Errors.Add($"{current}: (task) {e.Exception.InnerException?.Message ?? e.Exception.Message}"); e.SetObserved(); };

        var w = new MainWindow { Width = 1440, Height = 900 };
        w.Show();
        var vm = w.ViewModel;
        void Reopen()
        {
            w.OpenOnStart(Path.GetFullPath(sample), null);
            for (var k = 0; k < 400 && vm.Document == null; k++) Pump(1);
            Pump(20);
        }
        Reopen();
        Console.WriteLine($"info sample {sample}: {(vm.Document != null ? "open" : "NOT OPEN")}");

        // ---------------------------------------------------------------- every page
        for (var m = 0; m < MainViewModel.ModuleCount; m++)
        {
            if (only != null && !only.Contains(m)) continue;
            var crumb = MainViewModel.ModuleCrumb(m);
            if (crumb.Length == 0 && m != 8) continue;
            current = $"page {m} {(m == 8 ? "Studio" : crumb)}";
            if (vm.Document == null && m != 9) Reopen();
            vm.SetModule(m);
            Pump(30);
            StopEverything(vm);
            var seen = new HashSet<string>();
            var errorsBefore = Errors.Count;
            for (var pass = 0; pass < 2; pass++)   // a click can reveal more controls (a panel that opens): a second pass takes them
            {
                var controls = w.GetVisualDescendants().OfType<Control>()
                    .Where(c => c is Button or CheckBox or ComboBox or RadioButton && c.IsEffectivelyVisible && c.IsEffectivelyEnabled && c.Bounds.Width > 2 && c.Bounds.Height > 2)
                    .ToList();
                foreach (var c in controls)
                {
                    var key = Describe(c) + "@" + string.Join("/", c.GetVisualAncestors().Take(4).Select(a => a.GetType().Name));
                    if (!seen.Add(key) || !c.IsEffectivelyVisible || !c.IsEffectivelyEnabled || c.GetVisualRoot() == null) continue;
                    if (Dangerous(c)) continue;
                    var what = $"{current} · {Describe(c)}";
                    var before = Errors.Count;
                    try
                    {
                        Click(c);
                        _clicks++;
                        Pump(8);
                        if (c is ComboBox cb && cb.IsDropDownOpen)
                        {
                            // a drop-down: the next item, as a user picks one
                            if (cb.ItemCount > 1) cb.SelectedIndex = (Math.Max(0, cb.SelectedIndex) + 1) % cb.ItemCount;
                            cb.IsDropDownOpen = false;
                            Pump(6);
                        }
                    }
                    catch (Exception e) { Failures.Add($"{what}: threw {e.GetType().Name}: {e.Message}"); }
                    if (Errors.Count > before) Failures.Add($"{what}: {Errors[^1]}");
                    StopEverything(vm);
                    CloseExtras(w);
                    if (vm.Module != m) { vm.SetModule(m); Pump(10); }
                    if (vm.Document == null && m != 9) Reopen();
                    if (!w.IsVisible) { Failures.Add($"{what}: the main window closed"); w.Show(); }
                }
            }
            Console.WriteLine($"page {m,2} {(m == 8 ? "Studio" : crumb),-34} {seen.Count,4} controls{(Errors.Count > errorsBefore ? $" · {Errors.Count - errorsBefore} error(s)" : "")}");
        }

        // ---------------------------------------------------------------- the pane dividers: each dragged both ways
        foreach (var m in new[] { 8, 9, 0, 1, 11, 12 })
        {
            vm.SetModule(m);
            Pump(15);
            current = $"dividers on page {m}";
            foreach (var d in w.GetVisualDescendants().OfType<Control>().Where(c => c.GetType().Name == "PaneDivider" && c.IsEffectivelyVisible && c.Bounds.Width > 0 && c.Bounds.Height > 0).ToList())
            {
                var thumb = d.GetVisualDescendants().OfType<Thumb>().FirstOrDefault();
                if (thumb == null || TopLevel.GetTopLevel(thumb) is not { } top) continue;
                var rows = d.Bounds.Height < d.Bounds.Width;
                // grabbed a third of the way along, as a hand would (the middle holds the pill's hide and show arrows)
                if (thumb.TranslatePoint(rows ? new Point(thumb.Bounds.Width * 0.33, thumb.Bounds.Height / 2) : new Point(thumb.Bounds.Width / 2, thumb.Bounds.Height * 0.33), top) is not { } a) continue;
                var grid = d.GetVisualAncestors().OfType<Grid>().FirstOrDefault();
                string Sizes() => grid == null ? "" : string.Join(",", rows ? grid.RowDefinitions.Select(r => Math.Round(r.ActualHeight)) : grid.ColumnDefinitions.Select(cd => Math.Round(cd.ActualWidth)));
                // a divider under a page (the Studio's, beneath another module) is not the user's to drag there
                if (Avalonia.Input.InputExtensions.InputHitTest(top, a) is not Visual under || !under.GetVisualAncestors().Contains(d) && !ReferenceEquals(under, d)) continue;
                var s0 = Sizes();
                var before = Errors.Count;
                int started = 0, deltas = 0;
                thumb.DragStarted += (_, _) => started++;
                thumb.DragDelta += (_, _) => deltas++;
                var presses = "";
                void WinPress(object? s, PointerPressedEventArgs e) => presses += $"[window got press: source {(e.Source is Visual sv ? string.Join(" < ", new[] { sv }.Concat(sv.GetVisualAncestors()).Take(6).Select(x => x.GetType().Name + (x is Control xc && !string.IsNullOrEmpty(xc.Name) ? "#" + xc.Name : ""))) : "?")} · captured {(e.Pointer.Captured is Visual cv ? string.Join(" < ", new[] { cv }.Concat(cv.GetVisualAncestors()).Take(6).Select(x => x.GetType().Name + (x is Control xc && !string.IsNullOrEmpty(xc.Name) ? "#" + xc.Name : ""))) : "none")}]";
                top.AddHandler(InputElement.PointerPressedEvent, WinPress, Avalonia.Interactivity.RoutingStrategies.Tunnel, handledEventsToo: true);
                thumb.AddHandler(InputElement.PointerPressedEvent, (_, e) => presses += $"[tunnel handled={e.Handled} by {e.Source?.GetType().Name}]", Avalonia.Interactivity.RoutingStrategies.Tunnel, handledEventsToo: true);
                thumb.AddHandler(InputElement.PointerPressedEvent, (_, e) => presses += $"[bubble handled={e.Handled} left={e.GetCurrentPoint(thumb).Properties.IsLeftButtonPressed}]", Avalonia.Interactivity.RoutingStrategies.Bubble, handledEventsToo: true);
                foreach (var delta in new[] { 40.0, -40.0 })
                {
                    var b = rows ? new Point(a.X, a.Y + delta) : new Point(a.X + delta, a.Y);
                    top.MouseMove(a);
                    top.MouseDown(a, MouseButton.Left);
                    for (var k = 1; k <= 6; k++) { top.MouseMove(a + (b - a) * (k / 6.0), RawInputModifiers.LeftMouseButton); Pump(1); }
                    top.MouseUp(b, MouseButton.Left);
                    Pump(6);
                    if (delta > 0 && Sizes() == s0)
                    {
                        var hit = Avalonia.Input.InputExtensions.InputHitTest(top, a);
                        var chain = hit is Visual hv ? string.Join(" < ", new[] { hv }.Concat(hv.GetVisualAncestors()).Take(5).Select(x => x.GetType().Name)) : "nothing";
                        Failures.Add($"{current}: a {(rows ? "row" : "column")} divider did not move (sizes {s0} → {Sizes()}; drag started {started}, moves {deltas}, press {presses}; under the pointer: {chain}; grid {grid?.Name})");
                    }
                    a = b;
                }
                _dividers++;
                if (Errors.Count > before) Failures.Add($"{current}: {Errors[^1]}");
            }
        }
        // the Molecule module's wall between the sketch and the view
        vm.SetModule(9);
        Pump(15);
        if (w.GetVisualDescendants().OfType<Thumb>().FirstOrDefault(x => x.Name == "MolSketchWall" && x.IsEffectivelyVisible) is { } wall
            && TopLevel.GetTopLevel(wall) is { } wtop && wall.TranslatePoint(new Point(wall.Bounds.Width / 2, wall.Bounds.Height * 0.4), wtop) is { } wa)
        {
            var host = wall.GetVisualAncestors().OfType<Border>().First(x => x.Name == "MolSketchHost");
            var w0 = host.Bounds.Width;
            var wb = new Point(wa.X + 80, wa.Y);
            wtop.MouseMove(wa);
            wtop.MouseDown(wa, MouseButton.Left);
            for (var k = 1; k <= 6; k++) { wtop.MouseMove(wa + (wb - wa) * (k / 6.0), RawInputModifiers.LeftMouseButton); Pump(1); }
            wtop.MouseUp(wb, MouseButton.Left);
            Pump(6);
            _dividers++;
            if (Math.Abs(host.Bounds.Width - w0 - 80) > 12) Failures.Add($"Molecule: the wall between the sketch and the view moved {host.Bounds.Width - w0:F0} px for an 80 px drag");
        }
        else Failures.Add("Molecule: the wall between the sketch and the view is missing");
        Console.WriteLine($"dividers: {_dividers} dragged");

        // ---------------------------------------------------------------- the menu bar: every command
        vm.SetModule(8);
        Pump(10);
        if (NativeMenu.GetMenu(w) is { } menu)
        {
            var raise = typeof(NativeMenuItem).GetMethods(System.Reflection.BindingFlags.Instance | System.Reflection.BindingFlags.Public | System.Reflection.BindingFlags.NonPublic)
                .FirstOrDefault(x => x.Name.EndsWith("RaiseClicked", StringComparison.Ordinal) && x.GetParameters().Length == 0);
            void Walk(NativeMenu mm, string path)
            {
                foreach (var item in mm.Items.OfType<NativeMenuItem>().ToList())
                {
                    var p = path.Length == 0 ? item.Header ?? "" : path + " › " + item.Header;
                    if (item.Menu is { } sub) { Walk(sub, p); continue; }
                    if (Skip.Any(s => p.ToLowerInvariant().Contains(s)) || raise == null) continue;
                    current = "menu " + p;
                    var before = Errors.Count;
                    try { raise.Invoke(item, null); _menus++; Pump(10); }
                    catch (Exception e) { Failures.Add($"{current}: threw {(e.InnerException ?? e).Message}"); }
                    if (Errors.Count > before) Failures.Add($"{current}: {Errors[^1]}");
                    StopEverything(vm);
                    CloseExtras(w);
                    if (vm.Document == null) Reopen();
                }
            }
            Walk(menu, "");
            Console.WriteLine($"menu: {_menus} commands run{(raise == null ? " (could not invoke the items)" : "")}");
        }

        foreach (var f in Failures.Distinct()) Console.WriteLine("FAIL " + f);
        var rest = Errors.Where(e => Failures.All(f => !f.Contains(e))).Distinct().ToList();
        foreach (var e in rest) Console.WriteLine("ERROR " + e);
        Console.WriteLine($"{_clicks} controls clicked · {_menus} menu commands · {_dividers} dividers · {Failures.Distinct().Count()} failure(s) · {rest.Count} other error(s)");
        w.Close();
        try { Directory.Delete(tmp, true); } catch { }
        return Failures.Count == 0 && rest.Count == 0 ? 0 : 1;
    }
}
