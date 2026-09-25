using Avalonia.Controls;
using Avalonia.Input;
using Avalonia.Input.Platform;
using Avalonia.Interactivity;
using Avalonia.Markup.Xaml;
using Avalonia.Platform.Storage;
using CapsStudio.ViewModels;

namespace CapsStudio.Views.Pages;

public partial class EmptyPage : PageBase
{
    public EmptyPage()
    {
        AvaloniaXamlLoader.Load(this);
        AddHandler(DragDrop.DropEvent, (_, e) =>
        {
            if (e.Data.GetFiles() is { } files) Window?.OpenMany(files.Select(f => f.TryGetLocalPath()).OfType<string>().ToList());
        });
    }

    private async void OnOpen(object? s, RoutedEventArgs e) { if (Window != null) await Window.OpenWithPreview(); }
    private void OnGrow(object? s, RoutedEventArgs e) => Vm.SetModule(0);

    private async void OnPaste(object? s, RoutedEventArgs e)
    {
        var clip = TopLevel.GetTopLevel(this)?.Clipboard;
        var text = clip == null ? null : await clip.TryGetTextAsync();
        Vm.EmptyError = "";
        var paths = (text ?? "").Split('\n', StringSplitOptions.RemoveEmptyEntries | StringSplitOptions.TrimEntries)
            .Select(p => p.Trim('"', '\'')).Select(p => p.StartsWith('~') ? Path.Combine(Environment.GetFolderPath(Environment.SpecialFolder.UserProfile), p.TrimStart('~', '/')) : p)
            .ToList();
        var found = paths.Where(File.Exists).ToList();
        if (found.Count == 0)
        {
            Vm.EmptyError = paths.Count == 0 ? "The clipboard holds no path" : $"No file at {paths[0]}";
            return;
        }
        Window?.OpenMany(found);
    }

    private void OnRecent(object? s, RoutedEventArgs e)
    {
        if ((s as Control)?.Tag is not RecentItem r) return;
        if (!File.Exists(r.Path)) { Vm.EmptyError = $"{r.Name} is no longer at {r.Path}"; RecentFiles.Forget(r.Path); Vm.LoadRecent(); return; }
        Window?.OpenMany(r.Topology != null && File.Exists(r.Topology) ? [r.Path, r.Topology] : [r.Path]);
    }

    private void OnSample(object? s, RoutedEventArgs e) => Window?.OpenSample((s as Control)?.Tag as string ?? "ps");
}
