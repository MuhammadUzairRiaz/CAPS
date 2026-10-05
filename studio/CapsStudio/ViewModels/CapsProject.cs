using System.Globalization;
using System.Text.Json;
using System.Text.Json.Nodes;

namespace CapsStudio.ViewModels;

/// <summary>A CAPS project on disk (design/boards/ProjectsStart, NewProject, ProjectOpen): a folder the user puts anywhere
/// holding NAME.capsproj — CAPS's own file, readable JSON naming the project's structures and its sessions — beside
/// structures/ (each structure as LAMMPS data with its force field and provenance), runs/ and exports/.</summary>
public static class CapsProjectFile
{
    public const string Extension = ".capsproj";
    public const string Format = "caps-project";
    private static readonly JsonSerializerOptions Indented = new() { WriteIndented = true };

    /// <summary>A name a disk can hold: path separators and characters Windows refuses become "-"; never empty.</summary>
    public static string SafeName(string name)
    {
        var bad = Path.GetInvalidFileNameChars().Concat(['<', '>', ':', '"', '/', '\\', '|', '?', '*']).ToHashSet();
        var s = new string(name.Trim().Select(c => bad.Contains(c) || char.IsControl(c) ? '-' : c).ToArray()).Trim(' ', '.');
        return s.Length == 0 ? "CAPS project" : s.Length > 80 ? s[..80].TrimEnd(' ', '.') : s;
    }

    /// <summary>Makes PARENT/NAME/NAME.capsproj with structures/, runs/ and exports/; an existing folder is never used
    /// (NAME 2, NAME 3 … instead). Returns the project file's path.</summary>
    public static string Create(string parent, string name)
    {
        var safe = SafeName(name);
        var dir = Path.Combine(parent, safe);
        for (var k = 2; Directory.Exists(dir) || File.Exists(dir); ++k) dir = Path.Combine(parent, $"{safe} {k}");
        Directory.CreateDirectory(dir);
        foreach (var sub in new[] { "structures", "runs", "exports" }) Directory.CreateDirectory(Path.Combine(dir, sub));
        var file = Path.Combine(dir, Path.GetFileName(dir) + Extension);
        var now = DateTime.Now.ToString("o", CultureInfo.InvariantCulture);
        Write(file, new JsonObject
        {
            ["format"] = Format, ["version"] = 1, ["name"] = name.Trim().Length > 0 ? name.Trim() : safe, ["created"] = now, ["saved"] = now,
            ["structures"] = new JsonArray(), ["sessions"] = new JsonArray(),
        });
        return file;
    }

    public static JsonObject Read(string file)
    {
        var o = JsonNode.Parse(File.ReadAllText(file)) as JsonObject ?? throw new InvalidDataException("not a CAPS project file");
        if ((string?)o["format"] != Format) throw new InvalidDataException("not a CAPS project file (format is not caps-project)");
        o["structures"] ??= new JsonArray();
        o["sessions"] ??= new JsonArray();
        return o;
    }

    /// <summary>Written whole to a side file first, then moved over: a crash mid-write never leaves half a project file.</summary>
    public static void Write(string file, JsonObject o)
    {
        var tmp = file + ".saving";
        File.WriteAllText(tmp, o.ToJsonString(Indented));
        File.Move(tmp, file, true);
    }

    public static string FolderOf(string file) => Path.GetDirectoryName(Path.GetFullPath(file))!;
}

/// <summary>A project CAPS knows (made or opened here), wherever its folder is.</summary>
public sealed class KnownProject : ObservableObject
{
    public string File { get; init; } = "";
    private string _name = "", _counts = "", _when = "";
    private bool _missing, _current;
    private Avalonia.Media.Imaging.Bitmap? _thumb;
    public string Name { get => _name; set => Set(ref _name, value); }
    public string Folder => System.IO.Path.GetDirectoryName(File) ?? "";
    public string FolderText => RecentFiles.Tilde(Folder);
    public string Counts { get => _counts; set => Set(ref _counts, value); }
    public string When { get => _when; set => Set(ref _when, value); }
    public bool Missing { get => _missing; set { if (Set(ref _missing, value)) Raise(nameof(Found)); } }
    public bool Found => !_missing;
    public bool Current { get => _current; set => Set(ref _current, value); }
    public Avalonia.Media.Imaging.Bitmap? Thumbnail { get => _thumb; set { if (Set(ref _thumb, value)) Raise(nameof(HasThumbnail)); } }
    public bool HasThumbnail => _thumb != null;
    public DateTime LastOpened { get; set; }
}

/// <summary>The list of every project CAPS has made or opened: ~/.caps/projects.json (beside the settings).</summary>
public static class ProjectRegistry
{
    public static string File => AppSettings.Override != null
        ? Path.Combine(Path.GetDirectoryName(AppSettings.Override)!, "caps-projects.json")
        : Path.Combine(AppSettings.Folder, "projects.json");

    public sealed record Entry(string File, string Name, DateTime Opened);

    public static List<Entry> Load()
    {
        try
        {
            if (!System.IO.File.Exists(File)) return [];
            var j = JsonNode.Parse(System.IO.File.ReadAllText(File));
            return (j?["projects"] as JsonArray ?? []).OfType<JsonObject>()
                .Select(o => new Entry((string?)o["file"] ?? "", (string?)o["name"] ?? "",
                                       DateTime.TryParse((string?)o["opened"], CultureInfo.InvariantCulture, DateTimeStyles.RoundtripKind, out var t) ? t : DateTime.MinValue))
                .Where(e => e.File.Length > 0).ToList();
        }
        catch { return []; }
    }

    private static void Save(List<Entry> list)
    {
        try
        {
            Directory.CreateDirectory(Path.GetDirectoryName(File)!);
            var a = new JsonArray(list.Select(e => (JsonNode)new JsonObject
            {
                ["file"] = e.File, ["name"] = e.Name, ["opened"] = e.Opened.ToString("o", CultureInfo.InvariantCulture),
            }).ToArray());
            System.IO.File.WriteAllText(File, new JsonObject { ["format"] = "caps-projects", ["projects"] = a }.ToJsonString(new JsonSerializerOptions { WriteIndented = true }));
        }
        catch { }
    }

    private static bool Same(string a, string b) => string.Equals(Path.GetFullPath(a), Path.GetFullPath(b),
        OperatingSystem.IsWindows() || OperatingSystem.IsMacOS() ? StringComparison.OrdinalIgnoreCase : StringComparison.Ordinal);

    /// <summary>The project at FILE was made or opened now: first in the list.</summary>
    public static void Touch(string file, string name)
    {
        var list = Load().Where(e => !Same(e.File, file)).ToList();
        list.Insert(0, new Entry(Path.GetFullPath(file), name, DateTime.Now));
        Save(list);
    }

    /// <summary>A moved project found again: its entry points at the new file.</summary>
    public static void Relocate(string oldFile, string newFile, string name)
    {
        var list = Load();
        var k = list.FindIndex(e => Same(e.File, oldFile));
        var e = new Entry(Path.GetFullPath(newFile), name, k >= 0 ? list[k].Opened : DateTime.Now);
        list = list.Where(x => !Same(x.File, oldFile) && !Same(x.File, newFile)).ToList();
        list.Insert(Math.Max(0, Math.Min(k, list.Count)), e);
        Save(list);
    }

    /// <summary>Dropped from the list; the files stay where they are.</summary>
    public static void Forget(string file) => Save(Load().Where(e => !Same(e.File, file)).ToList());
}
