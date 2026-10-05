using Avalonia.VisualTree;
using Avalonia.Controls;
using Avalonia.Interactivity;
using Avalonia.Platform.Storage;

namespace CapsStudio.Views;

/// <summary>Projects (design/boards/ProjectsStart, ProjectOpen): the file pickers and folder views the project list,
/// the project switcher and Start share.</summary>
public partial class MainWindow
{
    /// <summary>Open project…: the project's folder (its .capsproj is found inside; no need to go into the folder).</summary>
    public async Task OpenProjectFileDialog()
    {
        var picked = await StorageProvider.OpenFolderPickerAsync(new FolderPickerOpenOptions { Title = "Open a CAPS project: choose its folder", AllowMultiple = false });
        if (picked.Count > 0 && picked[0].TryGetLocalPath() is { } p) _vm.Status = _vm.OpenCapsProject(p);
    }

    /// <summary>A project whose folder moved: its .capsproj in the new place.</summary>
    public async Task LocateProjectDialog(ViewModels.KnownProject k)
    {
        var picked = await StorageProvider.OpenFolderPickerAsync(new FolderPickerOpenOptions { Title = $"Where is {k.Name}? Choose its folder", AllowMultiple = false });
        if (picked.Count > 0 && picked[0].TryGetLocalPath() is { } p) _vm.Status = _vm.LocateProject(k, p);
    }

    /// <summary>The folder in Finder, Explorer or the Linux file manager.</summary>
    public async Task ShowFolder(string folder)
    {
        if (!Directory.Exists(folder)) { _vm.Status = $"{folder} is not there"; return; }
        if (Launcher is { } l) await l.LaunchDirectoryInfoAsync(new DirectoryInfo(folder));
    }

    private void OnNewProjectButton(object? s, RoutedEventArgs e) => _vm.OpenNewProject();
    private async void OnOpenProjectButton(object? s, RoutedEventArgs e) => await OpenProjectFileDialog();

    // ---------------------------------------------------------------- the project switcher on the top bar

    /// <summary>For the screenshot tool: the project switcher open under its chip.</summary>
    public void OpenProjectMenuForShot()
    {
        if (this.GetVisualDescendants().OfType<Button>().FirstOrDefault(b => b.Classes.Contains("projchip")) is { } chip) OnProjectChip(chip, new RoutedEventArgs());
    }

    private void OnProjectChip(object? s, RoutedEventArgs e)
    {
        if (s is not Control anchor) return;
        var items = new List<Control>();
        foreach (var k in _vm.KnownProjects.Take(8))
        {
            var kk = k;
            var m = new MenuItem
            {
                Header = new StackPanel
                {
                    Spacing = 1,
                    Children =
                    {
                        new TextBlock { Text = k.Name + (k.Missing ? "  (not found)" : ""), FontSize = 12.5 },
                        new TextBlock { Text = k.FolderText, FontSize = 11, Classes = { "mono", "dim" } },
                    },
                },
                Icon = new Icon { Kind = k.Current ? "check" : k.Missing ? "alert" : "folder", Size = 14 },
            };
            m.Click += async (_, _) => { if (kk.Missing) await LocateProjectDialog(kk); else _vm.Status = _vm.OpenCapsProject(kk.File); };
            items.Add(m);
        }
        if (items.Count > 0) items.Add(new Separator());
        MenuItem M(string header, string icon, Func<Task> run, string gesture = "", bool enabled = true)
        {
            var m = new MenuItem { Header = header, Icon = new Icon { Kind = icon, Size = 14 }, IsEnabled = enabled,
                                   InputGesture = gesture.Length > 0 ? Avalonia.Input.KeyGesture.Parse(gesture) : null };
            m.Click += async (_, _) => await run();
            return m;
        }
        var has = _vm.HasCapsProject;
        items.Add(M(has ? $"Save · {_vm.CapsProjectSavedText}" : "Save", "save", () => { _vm.Status = _vm.SaveCapsProject(); return Task.CompletedTask; }, OperatingSystem.IsMacOS() ? "Cmd+S" : "Ctrl+S", has));
        items.Add(M("New project…", "plus", () => { _vm.OpenNewProject(); return Task.CompletedTask; }));
        items.Add(M("Open project…", "folder", OpenProjectFileDialog));
        items.Add(M("Show the folder", "file", () => ShowFolder(Path.GetDirectoryName(_vm.CapsProjectPath) ?? ""), "", has));
        if (has)
        {
            items.Add(new Separator());
            items.Add(M("Close project", "close", () => { _vm.Status = _vm.CloseCapsProject(); return Task.CompletedTask; }));
        }
        new ContextMenu { ItemsSource = items }.Open(anchor);
    }
}
