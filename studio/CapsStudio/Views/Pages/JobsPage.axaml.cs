using System.ComponentModel;
using Avalonia;
using Avalonia.Controls;
using Avalonia.Controls.Primitives;
using Avalonia.Controls.Shapes;
using Avalonia.Interactivity;
using Avalonia.Layout;
using Avalonia.Markup.Xaml;
using Avalonia.Media;
using Avalonia.Threading;
using CapsStudio.ViewModels;

namespace CapsStudio.Views.Pages;

public partial class JobsPage : PageBase
{
    private readonly LinePlot _a, _b;
    private readonly UniformGrid _stages;
    private readonly StackPanel _legend;
    private readonly ProgressBar _bar;
    private readonly Grid _curves;
    private readonly TextBlock _progressTitle;
    private Job? _watched;
    private readonly SelectableTextBlock _liveText;
    private readonly ScrollViewer _liveScroll;
    private readonly CheckBox _follow;
    private readonly TextBox _logSearch;
    private readonly TextBlock _cancelHostText;
    private bool _cancelArmed;
    private readonly DispatcherTimer _clock = new() { Interval = TimeSpan.FromSeconds(1) };

    public JobsPage()
    {
        AvaloniaXamlLoader.Load(this);
        _a = this.FindControl<LinePlot>("PlotA")!;
        _b = this.FindControl<LinePlot>("PlotB")!;
        _a.AutoRange = _b.AutoRange = true;
        _a.RefY = _b.RefY = null;
        _stages = this.FindControl<UniformGrid>("StageBar")!;
        _legend = this.FindControl<StackPanel>("StageLegend")!;
        _bar = this.FindControl<ProgressBar>("Bar")!;
        _curves = this.FindControl<Grid>("Curves")!;
        _progressTitle = this.FindControl<TextBlock>("ProgressTitle")!;
        _liveText = this.FindControl<SelectableTextBlock>("LiveText")!;
        _liveScroll = this.FindControl<ScrollViewer>("LiveScroll")!;
        _follow = this.FindControl<CheckBox>("FollowLog")!;
        _logSearch = this.FindControl<TextBox>("LogSearch")!;
        _cancelHostText = this.FindControl<TextBlock>("CancelHostText")!;
        _logSearch.TextChanged += (_, _) => ShowLive();
        _clock.Tick += (_, _) => _watched?.Tick();
        DataContextChanged += (_, _) =>
        {
            if (DataContext is not MainViewModel vm) return;
            vm.HookJobs();
            vm.JobSelected += Show;
            vm.JobCurvesChanged += ShowCurves;
            Show();
        };
    }

    private void Show()
    {
        if (_watched != null) _watched.PropertyChanged -= OnJob;
        _watched = (DataContext as MainViewModel)?.SelectedJob;
        if (_watched != null) _watched.PropertyChanged += OnJob;
        ShowStages();
        ShowCurves();
        ShowLive();
        DisarmCancel();
        _clock.IsEnabled = _watched?.IsRunning == true;
    }

    /// <summary>The run's own log: every line, or those with the search text; errors and warnings in their colours;
    /// kept at its end while Follow is on.</summary>
    private void ShowLive()
    {
        var text = _watched?.LiveLog ?? "";
        var q = _logSearch.Text?.Trim() ?? "";
        var lines = text.Split('\n');
        if (lines.Length > 4000) lines = lines[^4000..];
        _liveText.Inlines?.Clear();
        _liveText.Inlines ??= new Avalonia.Controls.Documents.InlineCollection();
        foreach (var line in lines)
        {
            if (q.Length > 0 && !line.Contains(q, StringComparison.OrdinalIgnoreCase)) continue;
            var run = new Avalonia.Controls.Documents.Run(line + "\n");
            var low = line.ToLowerInvariant();
            if (low.Contains("error") || low.Contains("failed") || low.Contains("caps job: ") || low.StartsWith("caps ", StringComparison.Ordinal) && low.Contains(": "))
                run.Foreground = this.FindResource("ErrB") as IBrush;
            else if (low.Contains("warning") || low.Contains("note:") || low.Contains("not converged"))
                run.Foreground = this.FindResource("WarnB") as IBrush ?? this.FindResource("AccB") as IBrush;
            _liveText.Inlines.Add(run);
        }
        if (_follow.IsChecked == true) Dispatcher.UIThread.Post(() => _liveScroll.ScrollToEnd(), DispatcherPriority.Background);
    }

    private void DisarmCancel() { _cancelArmed = false; _cancelHostText.Text = "Cancel on host"; }
    private async void OnCancelRemote(object? s, RoutedEventArgs e)
    {
        if (!_cancelArmed) { _cancelArmed = true; _cancelHostText.Text = "Confirm cancel"; return; }   // a second press confirms
        DisarmCancel();
        await Vm.CancelRemote(Vm.SelectedJob);
    }
    private async void OnLatestFrame(object? s, RoutedEventArgs e) => await Vm.LatestFrame(Vm.SelectedJob);
    private async void OnResumeRemote(object? s, RoutedEventArgs e) => await Vm.ResumeRemote(Vm.SelectedJob);
    private async void OnCopyTerminal(object? s, RoutedEventArgs e)
    {
        var text = Vm.TerminalCommands(Vm.SelectedJob);
        if (text.Length == 0 || Window?.Clipboard is not { } cb) return;
        try { await cb.SetTextAsync(text); Vm.Status = "Copied: ssh, caps job status, caps job tail -f"; }
        catch { Vm.Status = text.Replace("\n", "  ·  "); }
    }

    private void OnJob(object? s, PropertyChangedEventArgs e)
    {
        if (e.PropertyName is nameof(Job.Stage) or nameof(Job.Stages) or nameof(Job.Status)) Dispatcher.UIThread.Post(ShowStages);
        if (e.PropertyName == nameof(Job.Status)) _clock.IsEnabled = _watched?.IsRunning == true;
        if (e.PropertyName == nameof(Job.LiveLog)) Dispatcher.UIThread.Post(ShowLive);
        if (e.PropertyName == nameof(Job.HasCurves)) Dispatcher.UIThread.Post(ShowCurves);
    }

    private void ShowCurves()
    {
        var j = _watched;
        _curves.IsVisible = j?.HasCurves == true;
        if (j == null) return;
        _a.XLabel = _b.XLabel = j.AxisX;
        _a.YLabel = j.AxisA;
        _b.YLabel = j.AxisB;
        _a.SetData(j.A.ToArray());
        _b.SetData(j.B.ToArray());
    }

    /// <summary>The protocol as a row of stages (done, running, remaining), or a plain bar for other runs.</summary>
    private void ShowStages()
    {
        _stages.Children.Clear();
        _legend.Children.Clear();
        var j = _watched;
        var n = j?.Stages ?? 0;
        _stages.IsVisible = _legend.IsVisible = n > 0;
        _bar.IsVisible = n == 0;
        _progressTitle.Text = n > 0 ? "Protocol" : "Progress";
        if (j == null || n == 0) return;
        var done = j.IsDone ? n : Math.Max(0, j.Stage - 1);
        var running = j.IsRunning ? 1 : 0;
        for (var k = 1; k <= n; k++)
        {
            IBrush fill = k <= done ? Tokens.Brush("OkB") : k == j.Stage && j.IsRunning ? Tokens.Brush("AccB") : Tokens.Brush("Bg3B");
            _stages.Children.Add(new Border { Background = fill, Opacity = k <= done ? 0.75 : 1, CornerRadius = new CornerRadius(3), Margin = new Thickness(0, 0, 3, 0) });
        }
        void Chip(string text, string brush)
        {
            _legend.Children.Add(new Border
            {
                Classes = { "chip" },
                Child = new StackPanel
                {
                    Orientation = Orientation.Horizontal, Spacing = 6,
                    Children =
                    {
                        new Ellipse { Width = 7, Height = 7, Fill = Tokens.Brush(brush), VerticalAlignment = VerticalAlignment.Center },
                        new TextBlock { Text = text, Foreground = Tokens.Brush("MutedB") },
                    },
                },
            });
        }
        Chip($"done {done}", "OkB");
        if (running > 0) Chip("running 1", "AccB");
        Chip($"remaining {Math.Max(0, n - done - running)}", "DimB");
    }

    private void OnOpenStudio(object? s, RoutedEventArgs e) => Vm.SetModule(8);
    private void OnGoModule(object? s, RoutedEventArgs e) { if (Vm.SelectedJob is { } j) Vm.SetModule(j.Module); }
    private async void OnRestartFix(object? s, RoutedEventArgs e) { if (Vm.SelectedJob is { } j) await Vm.RestartWithFix(j); }
    private void OnSuggest(object? s, RoutedEventArgs e) { if (Vm.SelectedJob is { SuggestModule: >= 0 } j) Vm.SetModule(j.SuggestModule); }
    private void OnCancel(object? s, RoutedEventArgs e) => Vm.CancelJob(Vm.SelectedJob);
    private void OnPause(object? s, RoutedEventArgs e) => Vm.TogglePause();
    private async void OnCheckRemote(object? s, RoutedEventArgs e) => await Vm.RefreshJob(Vm.SelectedJob);
    private void OnOpenRemote(object? s, RoutedEventArgs e) => Vm.OpenRemoteResult(Vm.SelectedJob);
    private void OnClear(object? s, RoutedEventArgs e) => Vm.ClearJobs();
    private void OnNew(object? s, RoutedEventArgs e) { if ((s as MenuItem)?.Tag is string t && int.TryParse(t, out var m)) Vm.SetModule(m); }

    private async void OnCopyJson(object? s, RoutedEventArgs e)
    {
        if (Vm.SelectedJob is not { } j || Window?.Clipboard is not { } cb) return;
        await cb.SetTextAsync(Vm.JobJson(j));
        Vm.Status = $"Copied {j.Id} (log and provenance) as JSON";
    }

    private void OnFullProvenance(object? s, RoutedEventArgs e) => Vm.OpenProvenance();
    private void OnSweep(object? s, RoutedEventArgs e) => Vm.OpenSweep();
}
