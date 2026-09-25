using Avalonia;
using Avalonia.Controls;
using Avalonia.Controls.Shapes;
using Avalonia.Interactivity;
using Avalonia.Layout;
using Avalonia.Markup.Xaml;
using Avalonia.Media;
using CapsStudio.ViewModels;

namespace CapsStudio.Views;

public partial class TourOverlay : UserControl
{
    private Rect _target;

    public TourOverlay()
    {
        AvaloniaXamlLoader.Load(this);
        SizeChanged += (_, _) => Place();
    }

    private MainViewModel? Vm => DataContext as MainViewModel;

    /// <summary>Lights this rectangle (in the overlay's coordinates) and places the callout beside it.</summary>
    public void SetTarget(Rect r)
    {
        _target = r.Inflate(3);
        Place();
    }

    private void Place()
    {
        var W = Bounds.Width;
        var H = Bounds.Height;
        if (W <= 0 || H <= 0 || Vm is not { } vm) return;
        var t = _target;
        var dim = this.FindControl<Avalonia.Controls.Shapes.Path>("Dim")!;
        var geo = new PathGeometry { FillRule = FillRule.EvenOdd };
        geo.Figures!.Add(Figure(new Rect(0, 0, W, H)));
        geo.Figures.Add(Figure(t));
        dim.Data = geo;
        var spot = this.FindControl<Border>("Spot")!;
        Canvas.SetLeft(spot, t.X);
        Canvas.SetTop(spot, t.Y);
        spot.Width = t.Width;
        spot.Height = t.Height;

        var call = this.FindControl<Border>("Callout")!;
        call.Measure(new Size(420, double.PositiveInfinity));
        var ch = call.DesiredSize.Height;
        double x, y;
        if (t.Height < 120 && t.Y + t.Height + ch + 24 < H) { x = t.X + 60; y = t.Bottom + 16; }            // a bar at the top: below it
        else if (t.Height < 320 && t.Y - ch - 16 > 0) { x = t.X + 40; y = t.Y - ch - 16; }                 // a dock at the bottom: above it
        else if (t.Width < 420 && t.Right + 440 < W) { x = t.Right + 16; y = t.Y + 60; }                   // a tall column on the left: beside it
        else if (t.Width < 420 && t.X - 436 > 0) { x = t.X - 436; y = t.Y + 60; }                          // a tall column on the right
        else { x = t.X + 24; y = t.Y + 24; }                                                               // a large region: inside it
        Canvas.SetLeft(call, Math.Clamp(x, 8, Math.Max(8, W - 428)));
        Canvas.SetTop(call, Math.Clamp(y, 8, Math.Max(8, H - ch - 8)));

        var outline = this.FindControl<Border>("Outline")!;
        outline.Measure(new Size(290, double.PositiveInfinity));
        var oh = outline.DesiredSize.Height;
        // a corner clear of the lit region: bottom right, bottom left (past the rail), top right
        var corners = new[] { new Point(W - 290 - 24, H - oh - 40), new Point(96, H - oh - 40), new Point(W - 290 - 24, 100) };
        var at = corners.FirstOrDefault(c => !new Rect(c, new Size(290, oh)).Intersects(t), corners[0]);
        Canvas.SetLeft(outline, at.X);
        Canvas.SetTop(outline, at.Y);
        FillSteps(vm);
    }

    private static PathFigure Figure(Rect r) => new()
    {
        StartPoint = r.TopLeft, IsClosed = true, IsFilled = true,
        Segments = [new LineSegment { Point = r.TopRight }, new LineSegment { Point = r.BottomRight }, new LineSegment { Point = r.BottomLeft }],
    };

    private void FillSteps(MainViewModel vm)
    {
        var acc = Tokens.Brush("AccB");
        var dots = this.FindControl<ItemsControl>("Dots")!;
        dots.ItemsSource = MainViewModel.TourSteps.Select((_, k) => new Border
        {
            Width = k == vm.TourStep ? 18 : 7, Height = 7, CornerRadius = new CornerRadius(4),
            Background = k == vm.TourStep ? acc : Tokens.Brush("Bg3B"),
        }).ToList();
        var steps = this.FindControl<ItemsControl>("Steps")!;
        steps.ItemsSource = MainViewModel.TourSteps.Select((s, k) =>
        {
            var on = k == vm.TourStep;
            var badge = new Border
            {
                Width = 22, Height = 22, CornerRadius = new CornerRadius(11), Background = on ? acc : Tokens.Brush("Bg3B"),
                Child = new TextBlock { Text = (k + 1).ToString(), FontFamily = Tokens.Mono, FontSize = 11, HorizontalAlignment = HorizontalAlignment.Center, VerticalAlignment = VerticalAlignment.Center,
                                        Foreground = on ? Tokens.Brush("AccInkB") : Tokens.Brush("MutedB") },
            };
            var text = new StackPanel { Spacing = 2 };
            text.Children.Add(new TextBlock { Text = s.Name, FontSize = 12.5, FontWeight = FontWeight.SemiBold });
            text.Children.Add(new TextBlock { Text = s.Outline, FontSize = 11.5, Foreground = Tokens.Brush("MutedB"), TextWrapping = TextWrapping.Wrap, MaxWidth = 220 });
            var row = new StackPanel { Orientation = Orientation.Horizontal, Spacing = 10, Opacity = on ? 1 : 0.75 };
            row.Children.Add(badge);
            row.Children.Add(text);
            return row;
        }).ToList();
    }

    private void OnNext(object? s, RoutedEventArgs e) => Vm?.TourNext();
    private void OnBack(object? s, RoutedEventArgs e) => Vm?.TourBack();
    private void OnSkip(object? s, RoutedEventArgs e) => Vm?.EndTour();
}
