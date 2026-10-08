using Avalonia;
using Avalonia.Controls;
using Avalonia.Input;
using Avalonia.Interactivity;
using Avalonia.Layout;
using Avalonia.Markup.Xaml;
using Avalonia.Media;
using CapsStudio.ViewModels;

namespace CapsStudio.Views.Pages;

public partial class ViewportsPage : PageBase
{
    private readonly Canvas _stage;
    private readonly Dictionary<ViewportTile, Border> _tiles = new();

    public ViewportsPage()
    {
        AvaloniaXamlLoader.Load(this);
        _stage = this.FindControl<Canvas>("Stage")!;
        _stage.SizeChanged += (_, _) => Layout();
        DataContextChanged += (_, _) =>
        {
            if (DataContext is not MainViewModel vm) return;
            foreach (var t in vm.ViewportTiles) _tiles[t] = MakeTile(t);
            vm.ViewportsArranged += Layout;
        };
        KeyDown += (_, e) =>
        {
            if (DataContext is not MainViewModel vm || !vm.IsViewports) return;
            if (e.Key == Key.Enter) { vm.MaximiseViewport(); e.Handled = true; }
            else if (e.Key == Key.Tab) { vm.CycleViewport(); e.Handled = true; }
        };
    }

    private Border MakeTile(ViewportTile t)
    {
        var image = new Image { Stretch = Stretch.Fill };
        image.Bind(Image.SourceProperty, new Avalonia.Data.Binding(nameof(ViewportTile.Image)) { Source = t });
        var title = new Border
        {
            Classes = { "chip" }, Margin = new Thickness(10), HorizontalAlignment = HorizontalAlignment.Left, VerticalAlignment = VerticalAlignment.Top,
            Cursor = new Cursor(StandardCursorType.Hand),
        };
        var label = new TextBlock();
        label.Bind(TextBlock.TextProperty, new Avalonia.Data.Binding(nameof(ViewportTile.Label)) { Source = t });
        title.Child = label;
        // right-click the title: the view this pane shows
        var menu = new ContextMenu();
        for (var k = 0; k < MainViewModel.ViewportViews.Length; ++k)
        {
            var view = k;
            var item = new MenuItem { Header = MainViewModel.ViewportViews[k].Name + (MainViewModel.ViewportViews[k].Ortho ? " · ortho" : "") };
            item.Click += (_, _) => Vm.SwapViewport(t, view);
            menu.Items.Add(item);
        }
        title.ContextMenu = menu;
        ToolTip.SetTip(title, "Click: make active · right-click: show another view");
        title.Bind(Border.BackgroundProperty, new Avalonia.Data.Binding(nameof(ViewportTile.Active))
            { Source = t, Converter = new Avalonia.Data.Converters.FuncValueConverter<bool, IBrush>(a => Tokens.Brush(a ? "Bg3B" : "HudB")) });
        var tripod = new TripodView { Width = 56, Height = 56, HorizontalAlignment = HorizontalAlignment.Right, VerticalAlignment = VerticalAlignment.Bottom, Margin = new Thickness(6) };
        tripod.Bind(IsVisibleProperty, new Avalonia.Data.Binding(nameof(MainViewModel.ViewportTripod)) { Source = DataContext });
        var panel = new Panel { Children = { image, title, tripod } };
        var border = new Border { Child = panel, BorderThickness = new Thickness(1), CornerRadius = new CornerRadius(4), ClipToBounds = true, Background = Tokens.Brush("Bg0B") };
        border.Bind(Border.BorderBrushProperty, new Avalonia.Data.Binding(nameof(ViewportTile.Active))
            { Source = t, Converter = new Avalonia.Data.Converters.FuncValueConverter<bool, IBrush>(a => Tokens.Brush(a ? "AccB" : "LineB")) });
        border.Bind(IsVisibleProperty, new Avalonia.Data.Binding(nameof(ViewportTile.Shown)) { Source = t });
        Point? press = null, last = null;
        bool moved = false;
        border.PointerPressed += (_, e) =>
        {
            Focus();
            Vm.ActivateViewport(t);
            press = last = e.GetPosition(border);
            moved = false;
            e.Pointer.Capture(border);
        };
        border.PointerMoved += (_, e) =>
        {
            if (last is not { } l) return;
            var p = e.GetPosition(border);
            if (press is { } p0 && Math.Abs(p.X - p0.X) + Math.Abs(p.Y - p0.Y) > 3) moved = true;
            if (!moved) return;
            Vm.DragViewport(t, p.X - l.X, p.Y - l.Y, border.Bounds.Height);
            tripod.Set(t.Ortho ? t.Yaw : Vm.Camera.Yaw, t.Ortho ? t.Pitch : Vm.Camera.Pitch);
            last = p;
        };
        border.PointerReleased += (_, e) =>
        {
            if (!moved && press is { } p0)
            {
                var k = VisualRoot?.RenderScaling ?? 1;
                var sx = t.PixelW / Math.Max(1, border.Bounds.Width);
                Vm.PickInViewport(t, (int)(p0.X * sx), (int)(p0.Y * sx), e.KeyModifiers.HasFlag(KeyModifiers.Shift));
            }
            press = last = null;
            e.Pointer.Capture(null);
        };
        border.PointerWheelChanged += (_, e) => { Vm.ZoomViewport(t, e.Delta.Y); e.Handled = true; };
        tripod.Set(t.Yaw, t.Pitch);
        _stage.Children.Add(border);
        return border;
    }

    private void Layout()
    {
        if (DataContext is not MainViewModel vm) return;
        var W = _stage.Bounds.Width;
        var H = _stage.Bounds.Height;
        if (W < 20 || H < 20) return;
        const double gap = 6;
        foreach (var (t, b) in _tiles)
        {
            var x = t.NX * W + gap / 2;
            var y = t.NY * H + gap / 2;
            var w = Math.Max(10, t.NW * W - gap);
            var h = Math.Max(10, t.NH * H - gap);
            Canvas.SetLeft(b, x);
            Canvas.SetTop(b, y);
            b.Width = w;
            b.Height = h;
            if (t.Shown) vm.SetViewportSize(t, (int)w, (int)h);
            if (!t.Ortho) (((Panel)b.Child!).Children[2] as TripodView)?.Set(vm.Camera.Yaw, vm.Camera.Pitch);
        }
    }

    private void OnBack(object? s, RoutedEventArgs e) => Vm.SetModule(8);
}
