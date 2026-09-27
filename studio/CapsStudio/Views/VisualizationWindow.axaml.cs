using System.Text.Json.Nodes;
using Avalonia.Controls;
using Avalonia.Input;
using Avalonia.Interactivity;
using Avalonia.Markup.Xaml;
using Avalonia.Platform.Storage;
using CapsStudio.Interop;
using CapsStudio.ViewModels;

namespace CapsStudio.Views;

/// <summary>The structure alone, full screen: a copy of the open structure (edits and runs in the main window do not
/// change it), or, during a run, the run's live snapshots as they come. Style, colour, background, the cell and a
/// "show only" expression are this window's own; the main view is untouched.</summary>
public partial class VisualizationWindow : Window
{
    private readonly MainViewModel _vm;
    private readonly bool _live;
    private CapsDocument? _copy;     // the structure's copy (static view), disposed on close
    private CapsDocument? _shown;
    private string _expr = "";

    public VisualizationWindow() : this(null!, false) { }

    public VisualizationWindow(MainViewModel vm, bool live)
    {
        _vm = vm;
        _live = live;
        AvaloniaXamlLoader.Load(this);
        if (vm == null) return;
        var view = this.FindControl<MolView>("View")!;
        var style = this.FindControl<ComboBox>("StyleBox")!;
        var colour = this.FindControl<ComboBox>("ColourBox")!;
        var back = this.FindControl<ComboBox>("BackBox")!;
        var cell = this.FindControl<CheckBox>("CellBox")!;
        var expr = this.FindControl<TextBox>("Expr")!;
        style.SelectedIndex = Math.Clamp(vm.StyleIndex, 0, 4);
        style.SelectionChanged += (_, _) => Apply();
        colour.SelectionChanged += (_, _) => Apply();
        back.SelectionChanged += (_, _) => Apply();
        cell.IsCheckedChanged += (_, _) => Apply();
        expr.KeyDown += (_, e) => { if (e.Key == Key.Enter) { _expr = expr.Text ?? ""; Apply(); e.Handled = true; } };
        expr.LostFocus += (_, _) => { if ((expr.Text ?? "") != _expr) { _expr = expr.Text ?? ""; Apply(); } };
        this.FindControl<Button>("ZoomIn")!.Click += (_, _) => view.ZoomBy(1.25);
        this.FindControl<Button>("ZoomOut")!.Click += (_, _) => view.ZoomBy(0.8);
        this.FindControl<Button>("ResetBtn")!.Click += (_, _) => view.Reset();
        this.FindControl<Button>("HideBtn")!.Click += (_, _) => ToggleBar();
        this.FindControl<Button>("CloseBtn")!.Click += (_, _) => Close();
        this.FindControl<Button>("SaveBtn")!.Click += async (_, _) => await SavePng();
        KeyDown += (_, e) =>
        {
            if (e.Source is TextBox) return;
            if (e.Key == Key.Escape) { Close(); e.Handled = true; }
            else if (e.Key == Key.H) { ToggleBar(); e.Handled = true; }
            else if (e.Key is Key.OemPlus or Key.Add) { view.ZoomBy(1.25); e.Handled = true; }
            else if (e.Key is Key.OemMinus or Key.Subtract) { view.ZoomBy(0.8); e.Handled = true; }
            else if (e.Key == Key.R) { view.Reset(); e.Handled = true; }
            else if (e.Key == Key.F) { WindowState = WindowState == WindowState.FullScreen ? WindowState.Normal : WindowState.FullScreen; e.Handled = true; }
        };
        if (live) vm.PropertyChanged += OnVm;
        Closed += (_, _) =>
        {
            if (_live) _vm.PropertyChanged -= OnVm;
            view.Document = null;
            _copy?.Dispose();
        };
        Show(Source());
        Opened += (_, _) => { Activate(); WindowState = WindowState.FullScreen; };   // after it opens (macOS ignores it before)
    }

    private void OnVm(object? s, System.ComponentModel.PropertyChangedEventArgs e)
    {
        if (e.PropertyName == nameof(MainViewModel.RunLiveDoc) && _vm.RunLiveDoc is { } d) Show(d);
    }

    /// <summary>What to show: the run's latest snapshot, else a copy of the open structure.</summary>
    private CapsDocument? Source()
    {
        if (_live && _vm.RunLiveDoc is { } d) return d;
        if (_vm.Document is not { } doc) return null;
        try { _copy = doc.Copy("view"); } catch { return null; }
        return _copy;
    }

    private void Show(CapsDocument? doc)
    {
        var view = this.FindControl<MolView>("View")!;
        var first = view.Document == null;
        _shown = doc;
        ApplyLook(doc);
        view.Document = doc;
        if (first && doc != null) view.Reset();
        UpdateNote();
    }

    private void Apply()
    {
        ApplyLook(_shown);
        this.FindControl<MolView>("View")!.Refresh();
        UpdateNote();
    }

    private string _error = "";
    private void ApplyLook(CapsDocument? doc)
    {
        var view = this.FindControl<MolView>("View")!;
        var style = this.FindControl<ComboBox>("StyleBox")!.SelectedIndex;
        view.DrawStyle = Math.Clamp(style, 0, 4);
        view.ColourMode = Math.Clamp(this.FindControl<ComboBox>("ColourBox")!.SelectedIndex, 0, 3);
        view.Backdrop = this.FindControl<ComboBox>("BackBox")!.SelectedIndex - 1;
        view.ShowCell = this.FindControl<CheckBox>("CellBox")!.IsChecked == true;
        if (doc == null) return;
        // "show only": every atom hidden, then those of the expression in the chosen style
        _error = "";
        try
        {
            if (_expr.Trim().Length == 0) doc.SetAppearance("{\"active\":false}");
            else
            {
                var names = new[] { "ball_and_stick", "space_filling", "sticks", "no_hydrogens", "ribbon" };
                var look = new JsonObject
                {
                    ["active"] = true,
                    ["layers"] = new JsonArray(new JsonObject { ["expression"] = "", ["style"] = "hidden" },
                                               new JsonObject { ["expression"] = _expr.Trim(), ["style"] = names[Math.Clamp(style, 0, 4)] }),
                };
                doc.SetAppearance(look.ToJsonString());
            }
        }
        catch (Exception e) { _error = e.Message; }
    }

    private void UpdateNote()
    {
        var t = this.FindControl<TextBlock>("NoteText")!;
        var what = _live ? (_vm.RunLiveText.Length > 0 ? _vm.RunLiveText : "live view") : _vm.Title.Replace(" (unsaved)", "");
        t.Text = _error.Length > 0 ? "Expression: " + _error : what + (_expr.Trim().Length > 0 ? " · showing " + _expr.Trim() : "") +
                                                          " · drag rotates · ⇧ drag pans · wheel zooms · H hides the controls · F full screen";
    }

    private void ToggleBar()
    {
        var bar = this.FindControl<Border>("Bar")!;
        bar.IsVisible = !bar.IsVisible;
        this.FindControl<Border>("Note")!.IsVisible = bar.IsVisible;
    }

    private async Task SavePng()
    {
        var file = await StorageProvider.SaveFilePickerAsync(new FilePickerSaveOptions
        {
            Title = "Save the view",
            SuggestedFileName = "view.png",
            FileTypeChoices = [new FilePickerFileType("PNG image") { Patterns = ["*.png"] }],
        });
        if (file?.TryGetLocalPath() is not { } path) return;
        var ok = await this.FindControl<MolView>("View")!.SavePng(path, 2);
        this.FindControl<TextBlock>("NoteText")!.Text = ok ? "Saved " + System.IO.Path.GetFileName(path) : "Could not save the image";
    }
}
