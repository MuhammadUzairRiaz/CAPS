using System.Collections.ObjectModel;
using Avalonia.Input;

namespace CapsStudio.ViewModels;

/// <summary>A command of the palette with the key the user gave it (Settings › Your shortcuts).</summary>
public sealed class ShortcutRow : ObservableObject
{
    public required string Id { get; init; }
    public required string Title { get; init; }
    public string Section { get; init; } = "";
    public string BuiltIn { get; init; } = "";
    private string _gesture = "";
    /// <summary>Stored form: "Meta+Shift+G" (Avalonia's KeyGesture syntax); "" none.</summary>
    public string Gesture { get => _gesture; set { if (Set(ref _gesture, value)) { Raise(nameof(Display)); Raise(nameof(HasGesture)); } } }
    public bool HasGesture => _gesture.Length > 0;
    public string Display => _gesture.Length > 0 ? MainViewModel.ShowGesture(_gesture) : BuiltIn.Length > 0 ? BuiltIn : "—";
    private bool _recording;
    public bool Recording { get => _recording; set { if (Set(ref _recording, value)) Raise(nameof(ButtonText)); } }
    public string ButtonText => _recording ? "Press keys…" : _gesture.Length > 0 ? "Change" : "Set";
}

/// <summary>User keyboard shortcuts: any palette command can have a key, kept in the settings and checked before the
/// Studio's own keys (a key the Studio uses is flagged when chosen, as are keys the text fields need).</summary>
public sealed partial class MainViewModel
{
    public ObservableCollection<ShortcutRow> ShortcutRows { get; } = new();
    private string _shortcutFilter = "", _shortcutNote = "";
    public string ShortcutFilter { get => _shortcutFilter; set { if (Set(ref _shortcutFilter, value)) FillShortcutRows(); } }
    public string ShortcutNote { get => _shortcutNote; private set { if (Set(ref _shortcutNote, value)) Raise(nameof(HasShortcutNote)); } }
    public bool HasShortcutNote => _shortcutNote.Length > 0;
    private ShortcutRow? _recordingRow;
    public bool RecordingShortcut => _recordingRow != null;

    /// <summary>The Studio's own keys (Settings › Keyboard), as gestures, for conflict warnings.</summary>
    private static readonly (string Gesture, string What)[] BuiltInKeys =
    [
        ("Meta+K", "the command palette"), ("Meta+O", "Open"), ("Meta+S", "Save"), ("Meta+W", "Close"), ("Meta+Shift+O", "the most recent file"),
        ("Meta+Z", "Undo"), ("Meta+Shift+Z", "Redo"), ("Meta+F", "Select by query"), ("Meta+I", "Invert the selection"), ("Meta+H", "Add hydrogens"),
        ("Meta+Shift+C", "Clean the geometry"), ("A", "Auto-clean while editing"), ("Alt+I", "Invert the picked stereocentre"), ("Shift+E", "the periodic table"), ("Meta+E", "Export an image"),
        ("Meta+Shift+A", "Announce the selection"), ("Ctrl+Tab", "Next structure"), ("F1", "Theory manual"), ("R", "Reset the view"), ("F", "Frame the selection"),
        ("D1", "Front view"), ("D2", "Top view"), ("D3", "Side view"), ("D5", "Perspective"), ("Space", "Play"), ("Left", "Previous frame"), ("Right", "Next frame"),
        ("Escape", "Clear the selection"), ("OemCloseBrackets", "Grow the selection"), ("Delete", "Delete the picked atoms"), ("Back", "Delete the picked atoms"),
        ("Meta+D1", "Studio"), ("Meta+D2", "Build"), ("Meta+D3", "Polymer cell"), ("Meta+D4", "Force field"), ("Meta+D5", "Packing"), ("Meta+D6", "Minimise"),
        ("Meta+D7", "Equilibrate"), ("Meta+D8", "Dynamics"), ("Meta+D9", "Analyze"), ("Meta+OemPlus", "Interface scale"), ("Meta+OemMinus", "Interface scale"),
        ("Meta+OemComma", "Settings"), ("Meta+C", "Copy"), ("Meta+V", "Paste"), ("Meta+X", "Cut"), ("Meta+A", "Select all in text"), ("Meta+Q", "Quit"),
    ];

    /// <summary>Settings › Your shortcuts: the commands matching the filter (those with a key first when it is empty).</summary>
    public void FillShortcutRows()
    {
        if (!_modelCommands) { _modelCommands = true; AddModelCommands(); }
        var q = _shortcutFilter.Trim();
        var mine = _settings.Shortcuts;
        IEnumerable<PaletteCommand> list = _commands.Where(c => c.Id.Length > 0);
        list = q.Length == 0 ? list.Where(c => mine.ContainsKey(c.Id)).Concat(list.Where(c => !mine.ContainsKey(c.Id)).Take(8))
                             : list.Select(c => (c, s: Score(q, c.Title + " " + c.Keywords, c.Id))).Where(x => x.s > 0).OrderByDescending(x => x.s).Select(x => x.c).Take(14);
        ShortcutRows.Clear();
        foreach (var c in list.DistinctBy(c => c.Id))
            ShortcutRows.Add(new ShortcutRow { Id = c.Id, Title = c.Title, Section = c.Section, BuiltIn = c.Shortcut, Gesture = mine.GetValueOrDefault(c.Id, "") });
    }

    /// <summary>Set / Change: the next key pressed (with its modifiers) becomes the row's shortcut; Esc cancels.</summary>
    public void BeginRecordShortcut(ShortcutRow row)
    {
        if (_recordingRow != null) _recordingRow.Recording = false;
        _recordingRow = row;
        row.Recording = true;
        ShortcutNote = $"Press the keys for “{row.Title}” (Esc cancels)";
        Raise(nameof(RecordingShortcut));
    }

    /// <summary>A key while recording: stored unless it is only a modifier; true when the key was used.</summary>
    public bool RecordShortcutKey(Key key, KeyModifiers mods)
    {
        if (_recordingRow is not { } row) return false;
        if (key is Key.LeftShift or Key.RightShift or Key.LeftCtrl or Key.RightCtrl or Key.LeftAlt or Key.RightAlt or Key.LWin or Key.RWin) return true;
        row.Recording = false;
        _recordingRow = null;
        Raise(nameof(RecordingShortcut));
        if (key == Key.Escape && mods == KeyModifiers.None) { ShortcutNote = ""; return true; }
        var g = GestureText(key, mods);
        var notes = new List<string>();
        if (mods is KeyModifiers.None or KeyModifiers.Shift && key is >= Key.A and <= Key.Z or >= Key.D0 and <= Key.D9)
            notes.Add("a plain letter or digit does nothing while a text field has the focus");
        if (BuiltInKeys.FirstOrDefault(b => (OperatingSystem.IsMacOS() ? b.Gesture : b.Gesture.Replace("Meta", "Ctrl")) == g) is { What: { } what }) notes.Add($"it replaces the Studio's key for {what}");
        foreach (var (id, other) in _settings.Shortcuts.ToList())
            if (other == g && id != row.Id)
            {
                _settings.Shortcuts.Remove(id);
                notes.Add($"taken from “{_commands.FirstOrDefault(c => c.Id == id)?.Title ?? id}”");
                if (ShortcutRows.FirstOrDefault(r => r.Id == id) is { } r) r.Gesture = "";
            }
        _settings.Shortcuts[row.Id] = g;
        row.Gesture = g;
        _settings.Save();
        ShortcutNote = $"{ShowGesture(g)} runs “{row.Title}”" + (notes.Count > 0 ? " · " + string.Join(" · ", notes) : "");
        return true;
    }

    public void ClearShortcut(ShortcutRow row)
    {
        if (_settings.Shortcuts.Remove(row.Id)) _settings.Save();
        row.Gesture = "";
        ShortcutNote = "";
    }

    public void ResetShortcuts()
    {
        _settings.Shortcuts.Clear();
        _settings.Save();
        ShortcutNote = "Your shortcuts are cleared; the Studio's own keys are back";
        FillShortcutRows();
    }

    /// <summary>A key the user assigned: runs its command (true), unless it is not available now.</summary>
    public bool TryUserShortcut(Key key, KeyModifiers mods, bool inText)
    {
        if (_settings.Shortcuts.Count == 0) return false;
        if (inText && !mods.HasFlag(KeyModifiers.Meta) && !mods.HasFlag(KeyModifiers.Control) && key is not (>= Key.F1 and <= Key.F12)) return false;
        var g = GestureText(key, mods);
        var id = _settings.Shortcuts.FirstOrDefault(kv => kv.Value == g).Key;
        if (id == null) return false;
        if (!_modelCommands) { _modelCommands = true; AddModelCommands(); }
        if (RunCommand(id)) return true;
        Status = $"“{_commands.FirstOrDefault(c => c.Id == id)?.Title ?? id}” is not available now";
        return true;
    }

    /// <summary>The palette's rows for a query (tests).</summary>
    public List<PaletteRow> PaletteRowsFor(string query)
    {
        _paletteQuery = query;
        FilterPalette();
        return PaletteRows.Where(r => r.IsCommand).ToList();
    }

    /// <summary>The shortcut shown for a command in the palette: the user's key, else the Studio's.</summary>
    private string ShortcutShown(PaletteCommand c) => _settings.Shortcuts.TryGetValue(c.Id, out var g) ? ShowGesture(g) : c.Shortcut;

    public static string GestureText(Key key, KeyModifiers mods)
    {
        var parts = new List<string>();
        if (mods.HasFlag(KeyModifiers.Control)) parts.Add("Ctrl");
        if (mods.HasFlag(KeyModifiers.Alt)) parts.Add("Alt");
        if (mods.HasFlag(KeyModifiers.Shift)) parts.Add("Shift");
        if (mods.HasFlag(KeyModifiers.Meta)) parts.Add("Meta");
        parts.Add(key.ToString());
        return string.Join("+", parts);
    }

    /// <summary>"Meta+Shift+G" → "⇧ ⌘ G" on macOS, "Ctrl+Shift+G" elsewhere.</summary>
    public static string ShowGesture(string g)
    {
        var parts = g.Split('+');
        var key = parts[^1] switch
        {
            var d when d.Length == 2 && d[0] == 'D' && char.IsDigit(d[1]) => d[1..],
            "OemCloseBrackets" => "]", "OemOpenBrackets" => "[", "OemComma" => ",", "OemPeriod" => ".", "OemPlus" => "+", "OemMinus" => "−",
            "OemQuestion" => "/", "OemSemicolon" => ";", "OemQuotes" => "'", "Escape" => "Esc", "Back" => "⌫", "Return" or "Enter" => "↵",
            var k => k,
        };
        if (!OperatingSystem.IsMacOS()) return string.Join("+", parts[..^1].Select(m => m == "Meta" ? "Win" : m).Append(key));
        var sym = parts[..^1].Select(m => m switch { "Ctrl" => "⌃", "Alt" => "⌥", "Shift" => "⇧", "Meta" => "⌘", _ => m });
        return string.Join(" ", sym.Append(key));
    }
}
