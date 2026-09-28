using Avalonia;
using Avalonia.Controls;
using Avalonia.Interactivity;
using Avalonia.Layout;
using Avalonia.Markup.Xaml;
using Avalonia.Platform.Storage;
using Avalonia.VisualTree;
using CapsStudio.ViewModels;

namespace CapsStudio.Views.Pages;

public partial class SettingsPage : PageBase
{
    private static readonly (string Keys, string What)[] Shortcuts =
    [
        ("⌘ K", "Command palette"), ("⌘ O", "Open a structure"), ("⌘ S", "Save as LAMMPS data"), ("⌘ W", "Close the document"),
        ("R", "Reset the view"), ("← →", "Previous / next frame"), ("Space", "Play the trajectory"), ("⇧ click", "Measure up to four atoms"),
        ("⌘ + / ⌘ −", "Interface scale"), ("Esc", "Close the palette"),
        ("F", "Frame the selection"), ("1 · 2 · 3", "Front · top · side view"), ("5", "Perspective / orthographic"),
        ("⌘ Z / ⇧ ⌘ Z", "Undo / redo"), ("⌘ F", "Select by query"), ("⌘ I", "Invert the selection"), ("]", "Grow the selection one bond"),
        ("⌘ H", "Add hydrogens"), ("⇧ ⌘ C", "Clean the geometry (UFF)"), ("⌥ I", "Invert the picked stereocentre"), ("⇧ E", "Periodic table"),
        ("⌫", "Delete the picked atoms"), ("⌘ E", "Export an image"), ("⌘ 1 … ⌘ 9", "Studio · Build · Polymer cell · Force field · Packing · Minimise · Equilibrate · Dynamics · Analyze"),
        ("⌃ Tab", "Next structure"), ("F1", "Theory manual"),
    ];

    private readonly ItemsControl _nav;

    public SettingsPage()
    {
        AvaloniaXamlLoader.Load(this);
        _nav = this.FindControl<ItemsControl>("Nav")!;
        for (var k = 0; k < MainViewModel.SettingsTabs.Length; k++)
        {
            var tab = k;
            var b = new Button { Classes = { "nav" }, Tag = tab };
            b.Content = new StackPanel
            {
                Orientation = Orientation.Horizontal, Spacing = 10,
                Children =
                {
                    new Icon { Kind = MainViewModel.SettingsIcons[k], Size = 16, VerticalAlignment = VerticalAlignment.Center },
                    new TextBlock { Text = MainViewModel.SettingsTabs[k], VerticalAlignment = VerticalAlignment.Center },
                },
            };
            b.Click += (_, _) => { Vm.SettingsTab = tab; MarkNav(); };
            _nav.Items.Add(b);
        }
        var keys = this.FindControl<ItemsControl>("Keys")!;
        foreach (var (k, what) in Shortcuts)
            keys.Items.Add(new Grid
            {
                ColumnDefinitions = new ColumnDefinitions("*,Auto"), Height = 30,
                Children =
                {
                    new TextBlock { Text = what, FontSize = 12.5, VerticalAlignment = VerticalAlignment.Center },
                    new Border
                    {
                        [Grid.ColumnProperty] = 1, BorderBrush = Tokens.Brush("LineB"), BorderThickness = new Thickness(1), CornerRadius = new CornerRadius(4),
                        Padding = new Thickness(6, 1), VerticalAlignment = VerticalAlignment.Center,
                        Child = new TextBlock { Text = k, FontFamily = Tokens.Mono, FontSize = 11, Foreground = Tokens.Brush("MutedB") },
                    },
                },
            });
        DataContextChanged += (_, _) =>
        {
            if (DataContext is not MainViewModel vm) return;
            vm.PropertyChanged += (_, e) =>
            {
                if (e.PropertyName is nameof(MainViewModel.SetPalette)) MarkPalettes();
                if (e.PropertyName is nameof(MainViewModel.SettingsTab)) MarkNav();
            };
            MarkNav();
        };
        AttachedToVisualTree += (_, _) => { MarkNav(); MarkPalettes(); Vm.FillShortcutRows(); };
        LayoutUpdated += (_, _) => { if (!_marked) MarkPalettes(); };
    }

    private void OnRecordShortcut(object? s, RoutedEventArgs e) { if ((s as Control)?.Tag is ShortcutRow r) Vm.BeginRecordShortcut(r); }
    private void OnClearShortcut(object? s, RoutedEventArgs e) { if ((s as Control)?.Tag is ShortcutRow r) Vm.ClearShortcut(r); }
    private void OnResetShortcuts(object? s, RoutedEventArgs e) => Vm.ResetShortcuts();
    private void OnResolveConflict(object? s, RoutedEventArgs e) { if ((s as Control)?.Tag is ShortcutConflict c) Vm.ResolveShortcutConflict(c); }
    private void OnOpenShortcuts(object? s, RoutedEventArgs e) => Vm.SettingsTab = 4;   // Files: your shortcuts
    private async void OnCheckPython(object? s, RoutedEventArgs e) => await Vm.CheckPython();
    private async void OnCopyPythonLines(object? s, RoutedEventArgs e)
    {
        if (TopLevel.GetTopLevel(this)?.Clipboard is { } cb) { await cb.SetTextAsync(Vm.PythonShellLines); Vm.Status = "Copied the two lines for your shell"; }
    }

    private bool _marked;

    private void MarkNav()
    {
        if (DataContext is not MainViewModel vm) return;
        foreach (var b in _nav.Items.OfType<Button>()) b.Classes.Set("on", (int)b.Tag! == vm.SettingsTab);
    }

    private void MarkPalettes()
    {
        if (DataContext is not MainViewModel vm) return;
        var buttons = this.GetVisualDescendants().OfType<Button>().Where(b => b.Tag is int && b.Classes.Contains("tile")).ToList();
        foreach (var b in buttons) b.Classes.Set("sel", (int)b.Tag! == vm.SetPalette);
        _marked = buttons.Count > 0;
    }

    private void OnThemeDark(object? s, RoutedEventArgs e) => Vm.ThemeDark = true;
    private void OnThemeLight(object? s, RoutedEventArgs e) => Vm.ThemeLight = true;
    private void OnThemeSystem(object? s, RoutedEventArgs e) => Vm.ThemeSystem = true;
    private void OnPalette(object? s, RoutedEventArgs e) { if ((s as Control)?.Tag is int p) Vm.SetPalette = p; }
    private void OnClearRecent(object? s, RoutedEventArgs e) => Vm.ClearRecent();
    private void OnReset(object? s, RoutedEventArgs e) => Vm.ResetSettings();

    private async void OnExport(object? s, RoutedEventArgs e)
    {
        if (Window == null) return;
        var f = await Window.StorageProvider.SaveFilePickerAsync(new FilePickerSaveOptions { Title = "Export settings", SuggestedFileName = "caps-settings.json", DefaultExtension = "json" });
        if (f?.TryGetLocalPath() is { } p) Vm.ExportSettings(p);
    }

    private async void OnImport(object? s, RoutedEventArgs e)
    {
        if (Window == null) return;
        var f = await Window.StorageProvider.OpenFilePickerAsync(new FilePickerOpenOptions
        {
            Title = "Import settings", AllowMultiple = false, FileTypeFilter = [new FilePickerFileType("Settings") { Patterns = ["*.json"] }],
        });
        if (f.Count > 0 && f[0].TryGetLocalPath() is { } p) Vm.ImportSettings(p);
    }

    private void OnAddHost(object? s, RoutedEventArgs e) => Vm.AddHost();
    private void OnRemoveHost(object? s, RoutedEventArgs e) => Vm.RemoveHost();
    private async void OnTestHost(object? s, RoutedEventArgs e) => await Vm.TestHost();
    private void OnResetTemplate(object? s, RoutedEventArgs e) => Vm.ResetJobTemplate();
    private void OnColourVision(object? s, RoutedEventArgs e) => Vm.OpenColourVision();
    private async void OnCheckUpdates(object? s, RoutedEventArgs e) => await Vm.CheckForUpdates();
}
