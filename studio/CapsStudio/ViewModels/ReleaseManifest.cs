using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Runtime.InteropServices;
using System.Text.Json.Nodes;

namespace CapsStudio.ViewModels;

public sealed record ReleaseHighlight(string Title, string Detail, string Page);

/// <summary>release.json (format caps-release/1): a CAPS release's version, build commit and its date, title, what is
/// new, and its installers per platform. The same file is uploaded beside the installers and bundled in each app as
/// data/release.json (tools/release/make_release_json.py writes both). The rules match the core's caps/update.hpp, which
/// caps --check-update uses.</summary>
public sealed class ReleaseManifest
{
    public string Version { get; init; } = "";
    public string Commit { get; init; } = "";
    public string Date { get; init; } = "";
    public string Title { get; init; } = "";
    public string NotesUrl { get; init; } = "";
    public List<ReleaseHighlight> Highlights { get; init; } = new();
    public Dictionary<string, (string File, string Sha256)> Assets { get; init; } = new();

    /// <summary>Throws InvalidDataException on a file that is not caps-release/1 or lacks version, commit or date.</summary>
    public static ReleaseManifest Parse(string json)
    {
        JsonNode? j;
        try { j = JsonNode.Parse(json); }
        catch (Exception e) { throw new InvalidDataException("release manifest: not JSON (" + e.Message + ")"); }
        if (j is not JsonObject o) throw new InvalidDataException("release manifest: not a JSON object");
        static string S(JsonNode? n) => n is JsonValue v && v.TryGetValue<string>(out var s) ? s : "";
        if (!S(o["format"]).StartsWith("caps-release/1", StringComparison.Ordinal))
            throw new InvalidDataException($"release manifest: format is '{S(o["format"])}', not caps-release/1");
        var m = new ReleaseManifest
        {
            Version = S(o["version"]).Trim(), Commit = S(o["commit"]).Trim().ToLowerInvariant(), Date = S(o["date"]).Trim(),
            Title = S(o["title"]), NotesUrl = S(o["notes_url"]),
        };
        if (m.Version.Length == 0 || m.Commit.Length == 0 || m.Date.Length == 0)
            throw new InvalidDataException("release manifest: version, commit and date are required");
        if (o["highlights"] is JsonArray hs)
            foreach (var h in hs)
            {
                if (h is JsonValue) { if (S(h) is { Length: > 0 } t) m.Highlights.Add(new(t, "", "")); continue; }
                if (h is JsonObject ho && S(ho["title"]) is { Length: > 0 } title) m.Highlights.Add(new(title, S(ho["detail"]), S(ho["page"])));
            }
        if (o["assets"] is JsonObject assets)
            foreach (var (key, a) in assets)
            {
                if (a is JsonValue) { if (S(a) is { Length: > 0 } f) m.Assets[key] = (f, ""); }
                else if (a is JsonObject ao && S(ao["file"]) is { Length: > 0 } file) m.Assets[key] = (file, S(ao["sha256"]).ToLowerInvariant());
            }
        return m;
    }

    public static ReleaseManifest? TryLoad(string? path)
    {
        try { return path != null && File.Exists(path) ? Parse(File.ReadAllText(path)) : null; }
        catch { return null; }
    }

    public DateTimeOffset? When => UpdateRules.ParseDate(Date);
    /// <summary>"22 Oct 2026", or the date as written.</summary>
    public string DateText => When is { } t ? t.ToUniversalTime().ToString("d MMM yyyy", CultureInfo.InvariantCulture) : Date;
    public string ShortCommit => Commit.Length > 10 ? Commit[..10] : Commit;
}

public enum UpdateVerdict { Update, Current, LocalNewer, Unknown }

/// <summary>This build: version, commit and the commit's date (the native core's caps_build_info).</summary>
public sealed record BuildIdentity(string Version, string Commit, string Date);

public static class UpdateRules
{
    public const string ManifestUrl = "https://github.com/MuhammadUzairRiaz/CAPS/releases/latest/download/release.json";
    public const string ReleasesPage = "https://github.com/MuhammadUzairRiaz/CAPS/releases/latest";
    public const string DownloadBase = "https://github.com/MuhammadUzairRiaz/CAPS/releases/latest/download/";

    /// <summary>False in the self-test, --screenshot, --figures and the other tool modes: nothing there may reach the
    /// network (a local CAPS_UPDATE_FEED file still works, for tests).</summary>
    public static bool NetworkAllowed { get; set; } = true;

    /// <summary>CAPS_NO_UPDATE_CHECK set to anything but "", 0, false, no or off.</summary>
    public static bool Disabled(string? env = null)
    {
        var v = (env ?? Environment.GetEnvironmentVariable("CAPS_NO_UPDATE_CHECK"))?.Trim().ToLowerInvariant();
        return v != null && v is not ("" or "0" or "false" or "no" or "off");
    }

    /// <summary>-1, 0, 1 for a &lt; = &gt; b ("0.1.0" &lt; "0.2.0"; "0.2.0-rc1" &lt; "0.2.0"; a missing part is 0).</summary>
    public static int CompareVersions(string a, string b)
    {
        static (long[] Nums, string Pre) Split(string v)
        {
            var s = v.Trim().TrimStart('v', 'V');
            var pre = "";
            var cut = s.IndexOfAny(['-', '+']);
            if (cut >= 0)
            {
                if (s[cut] == '-') { pre = s[(cut + 1)..]; var plus = pre.IndexOf('+'); if (plus >= 0) pre = pre[..plus]; }
                s = s[..cut];
            }
            var parts = s.Split('.');
            var nums = new long[parts.Length];
            for (var k = 0; k < parts.Length; ++k)
            {
                long n = 0;
                foreach (var c in parts[k]) { if (c is < '0' or > '9') break; n = n * 10 + (c - '0'); }
                nums[k] = n;
            }
            return (nums, pre);
        }
        var (na, pa) = Split(a);
        var (nb, pb) = Split(b);
        for (var k = 0; k < Math.Max(na.Length, nb.Length); ++k)
        {
            var x = k < na.Length ? na[k] : 0;
            var y = k < nb.Length ? nb[k] : 0;
            if (x != y) return x < y ? -1 : 1;
        }
        if (pa == pb) return 0;
        if (pa.Length == 0) return 1;
        if (pb.Length == 0) return -1;
        return string.CompareOrdinal(pa, pb) < 0 ? -1 : 1;
    }

    public static DateTimeOffset? ParseDate(string? s) =>
        DateTimeOffset.TryParse(s?.Trim(), CultureInfo.InvariantCulture, DateTimeStyles.AssumeUniversal, out var t) ? t : null;

    /// <summary>One is a prefix of the other (a short and a full SHA), at least 7 characters, case ignored.</summary>
    public static bool SameCommit(string? a, string? b)
    {
        var x = (a ?? "").Trim().ToLowerInvariant();
        var y = (b ?? "").Trim().ToLowerInvariant();
        if (x.Length < 7 || y.Length < 7 || x == "unknown" || y == "unknown") return false;
        var n = Math.Min(x.Length, y.Length);
        return string.CompareOrdinal(x, 0, y, 0, n) == 0;
    }

    /// <summary>The version decides first; with equal versions, another commit with a later date is an update.</summary>
    public static UpdateVerdict Decide(BuildIdentity local, ReleaseManifest remote)
    {
        var v = CompareVersions(remote.Version, local.Version);
        if (v > 0) return UpdateVerdict.Update;
        if (v < 0) return UpdateVerdict.LocalNewer;
        if (SameCommit(local.Commit, remote.Commit)) return UpdateVerdict.Current;
        if (ParseDate(remote.Date) is not { } tr || ParseDate(local.Date) is not { } tl) return UpdateVerdict.Unknown;
        return tr > tl ? UpdateVerdict.Update : UpdateVerdict.LocalNewer;
    }

    /// <summary>The release file for this machine and the way CAPS was installed: macos-arm64; windows-x64 (the
    /// installer) or windows-x64-portable (a folder without the uninstaller); linux-x86_64 (AppImage), -deb (/opt/caps
    /// installed by dpkg), -pacman (/opt/caps otherwise: the PKGBUILD, so the release page) or -tar.</summary>
    public static string PlatformKey(string? baseDir = null)
    {
        baseDir ??= AppContext.BaseDirectory;
        var arm = RuntimeInformation.OSArchitecture == Architecture.Arm64;
        if (OperatingSystem.IsMacOS()) return arm ? "macos-arm64" : "macos-x86_64";
        if (OperatingSystem.IsWindows()) return File.Exists(Path.Combine(baseDir, "unins000.exe")) ? "windows-x64" : "windows-x64-portable";
        var linux = arm ? "linux-aarch64" : "linux-x86_64";
        if (Environment.GetEnvironmentVariable("APPIMAGE") is { Length: > 0 }) return linux;
        if (!baseDir.StartsWith("/opt/caps", StringComparison.Ordinal)) return linux + "-tar";
        return File.Exists("/var/lib/dpkg/info/caps-studio.list") ? linux + "-deb" : linux + "-pacman";
    }

    public static string DownloadUrl(ReleaseManifest m, string key)
    {
        if (!m.Assets.TryGetValue(key, out var a))
        {
            // an installer kind the release lacks: the platform's main file (a pacman build is rebuilt from the PKGBUILD instead)
            var dash = key.LastIndexOf('-');
            if (dash > key.IndexOf('-') && key[(dash + 1)..] is "deb" or "tar" or "portable") m.Assets.TryGetValue(key[..dash], out a);
        }
        if (!string.IsNullOrEmpty(a.File)) return DownloadBase + a.File;
        return m.NotesUrl.Length > 0 ? m.NotesUrl : ReleasesPage;
    }
}
