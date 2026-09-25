using System.Globalization;
using System.Security.Cryptography;
using System.Text;
using System.Text.Json;
using Avalonia.Media.Imaging;

namespace CapsStudio.ViewModels;

/// <summary>One structure opened or saved recently (the Start board's Recent row).</summary>
public sealed class RecentItem
{
    public string Path { get; set; } = "";
    public string? Topology { get; set; }
    public string Detail { get; set; } = "";
    public DateTime When { get; set; }
    public string Name => System.IO.Path.GetFileName(Path);
    public string Folder => System.IO.Path.GetDirectoryName(Path) is { } d ? RecentFiles.Tilde(d) : "";
    public string Subtitle => $"{Detail} · {RecentFiles.Ago(When)}";
    [System.Text.Json.Serialization.JsonIgnore] public Bitmap? Thumb { get; set; }
    [System.Text.Json.Serialization.JsonIgnore] public bool HasThumb => Thumb != null;
}

/// <summary>The recent-files list, kept next to the settings (~/.caps/recent.json) with a thumbnail of each structure
/// (~/.caps/thumbs). Nothing here is required: a missing or unreadable list is an empty one.</summary>
public static class RecentFiles
{
    public const int Max = 12;
    private static string Dir => AppSettings.Folder;
    /// <summary>Tests and screenshots point the list elsewhere so they never touch the user's own.</summary>
    public static string? Override { get; set; }
    private static string Root => Override ?? Dir;

    public static string ThumbPath(string path)
    {
        var h = Convert.ToHexString(SHA1.HashData(Encoding.UTF8.GetBytes(path)))[..16].ToLowerInvariant();
        return System.IO.Path.Combine(Root, "thumbs", h + ".t.png");   // transparent: sits on either theme's card
    }

    /// <summary>The thumbnail to show: the transparent one, else one from before (drawn on the dark view background).</summary>
    private static string? ShownThumb(string path)
    {
        var t = ThumbPath(path);
        if (File.Exists(t)) return t;
        var old = t[..^6] + ".png";
        return File.Exists(old) ? old : null;
    }

    public static List<RecentItem> Load()
    {
        try
        {
            var p = System.IO.Path.Combine(Root, "recent.json");
            if (!File.Exists(p)) return [];
            var list = JsonSerializer.Deserialize<List<RecentItem>>(File.ReadAllText(p)) ?? [];
            list = list.Where(x => File.Exists(x.Path)).Take(Max).ToList();
            foreach (var x in list)
            {
                var t = ShownThumb(x.Path);
                if (t != null) try { x.Thumb = new Bitmap(t); } catch { /* a broken thumbnail is left out */ }
            }
            return list;
        }
        catch { return []; }
    }

    /// <summary>Puts `path` first; `thumb` (optional) writes the thumbnail PNG to the given file.</summary>
    public static void Touch(string path, string? topology, string detail, Action<string>? thumb)
    {
        try
        {
            Directory.CreateDirectory(System.IO.Path.Combine(Root, "thumbs"));
            var list = Load();
            list.RemoveAll(x => string.Equals(x.Path, path, StringComparison.Ordinal));
            list.Insert(0, new RecentItem { Path = path, Topology = topology, Detail = detail, When = DateTime.Now });
            if (list.Count > Max) list.RemoveRange(Max, list.Count - Max);
            File.WriteAllText(System.IO.Path.Combine(Root, "recent.json"), JsonSerializer.Serialize(list, new JsonSerializerOptions { WriteIndented = true }));
            thumb?.Invoke(ThumbPath(path));
        }
        catch { /* the list is a convenience; failing to write it never stops the work */ }
    }

    public static void Clear()
    {
        try
        {
            var p = System.IO.Path.Combine(Root, "recent.json");
            if (File.Exists(p)) File.Delete(p);
            var t = System.IO.Path.Combine(Root, "thumbs");
            if (Directory.Exists(t)) Directory.Delete(t, true);
        }
        catch { }
    }

    public static void Forget(string path)
    {
        try
        {
            var list = Load();
            list.RemoveAll(x => x.Path == path);
            File.WriteAllText(System.IO.Path.Combine(Root, "recent.json"), JsonSerializer.Serialize(list, new JsonSerializerOptions { WriteIndented = true }));
        }
        catch { }
    }

    public static string Tilde(string dir)
    {
        var home = Environment.GetFolderPath(Environment.SpecialFolder.UserProfile);
        return home.Length > 0 && dir.StartsWith(home, StringComparison.Ordinal) ? "~" + dir[home.Length..] : dir;
    }

    public static string Ago(DateTime t)
    {
        var d = DateTime.Now - t;
        if (d.TotalMinutes < 1) return "just now";
        if (d.TotalHours < 1) return $"{(int)d.TotalMinutes} min ago";
        if (d.TotalHours < 24 && t.Date == DateTime.Today) return $"{(int)d.TotalHours} h ago";
        if (t.Date == DateTime.Today.AddDays(-1)) return "yesterday";
        return t.ToString("MMM d", CultureInfo.InvariantCulture);
    }
}
