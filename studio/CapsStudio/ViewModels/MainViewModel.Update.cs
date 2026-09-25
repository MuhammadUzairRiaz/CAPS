using System;
using System.Collections.ObjectModel;
using System.IO;
using System.Linq;
using System.Net.Http;
using System.Text.Json.Nodes;
using System.Threading.Tasks;

namespace CapsStudio.ViewModels;

public sealed record UpdateSection(string Title, string[] Items, bool AltersResults);

/// <summary>Update available (design/boards/UpdateDialog). Checks only when asked (Settings › Files, or the palette):
/// the latest release of the CAPS repository on GitHub, or CAPS_UPDATE_FEED (a file or URL with the same JSON) for tests.
/// Release notes are split at their "## " headings; a heading that mentions results is shown first, as changes that can
/// alter published numbers. Nothing installs itself: Download opens the release page.</summary>
public partial class MainViewModel
{
    public const string ReleasesApi = "https://api.github.com/repos/MuhammadUzairRiaz/CAPS/releases/latest";
    public static string CurrentVersion => typeof(MainViewModel).Assembly.GetName().Version?.ToString(3) ?? "0.1.0";

    private bool _updateOpen, _updateChecking;
    private string _updateVersion = "", _updateUrl = "", _updateText = "";
    public bool UpdateOpen { get => _updateOpen; set => Set(ref _updateOpen, value); }
    public bool UpdateChecking { get => _updateChecking; private set => Set(ref _updateChecking, value); }
    public string UpdateVersion { get => _updateVersion; private set { Set(ref _updateVersion, value); Raise(nameof(UpdateTitle)); } }
    public string UpdateTitle => $"CAPS {_updateVersion} is available";
    public string UpdateSubtitle => $"You have {CurrentVersion}";
    public string UpdateUrl { get => _updateUrl; private set => Set(ref _updateUrl, value); }
    /// <summary>The result of the last check, for Settings.</summary>
    public string UpdateText { get => _updateText; private set => Set(ref _updateText, value); }
    public ObservableCollection<UpdateSection> UpdateSections { get; } = new();

    public async Task CheckForUpdates(bool quietWhenCurrent = false)
    {
        if (_updateChecking) return;
        UpdateChecking = true;
        UpdateText = "Checking…";
        try
        {
            var json = await Task.Run(FetchRelease);
            var j = JsonNode.Parse(json)!;
            var tag = ((string?)j["tag_name"] ?? "").TrimStart('v', 'V');
            if (!Version.TryParse(tag, out var latest)) throw new InvalidDataException("the release has no version tag");
            UpdateUrl = (string?)j["html_url"] ?? "";
            if (latest <= Version.Parse(CurrentVersion))
            {
                UpdateText = $"CAPS {CurrentVersion} is the latest release · checked {DateTime.Now:HH:mm}";
                if (!quietWhenCurrent) Status = UpdateText;
                return;
            }
            UpdateVersion = tag;
            UpdateSections.Clear();
            foreach (var s in ParseReleaseNotes((string?)j["body"] ?? "")) UpdateSections.Add(s);
            UpdateText = $"CAPS {tag} is available · you have {CurrentVersion}";
            UpdateOpen = true;
        }
        catch (HttpRequestException e)
        {
            UpdateText = "Offline · " + e.Message;
            Notify(new Notice
            {
                Key = "offline", Severity = "info", Icon = "server", Title = "Offline",
                Body = "The theory manual, fragment and force-field libraries are on this computer. Remote hosts and update checks wait until you are back online.",
                Primary = "OK",
            });
        }
        catch (Exception e) { UpdateText = "Could not check: " + e.Message; }
        finally { UpdateChecking = false; }
    }

    private static string FetchRelease()
    {
        var feed = Environment.GetEnvironmentVariable("CAPS_UPDATE_FEED");
        if (feed is { Length: > 0 } && !feed.StartsWith("http", StringComparison.OrdinalIgnoreCase)) return File.ReadAllText(feed);
        using var http = new HttpClient { Timeout = TimeSpan.FromSeconds(10) };
        http.DefaultRequestHeaders.UserAgent.ParseAdd($"CAPS-Studio/{CurrentVersion}");
        http.DefaultRequestHeaders.Accept.ParseAdd("application/vnd.github+json");
        return http.GetStringAsync(feed is { Length: > 0 } ? feed : ReleasesApi).GetAwaiter().GetResult();
    }

    /// <summary>"## Heading" sections with "- item" lines; results-altering sections first.</summary>
    public static UpdateSection[] ParseReleaseNotes(string md)
    {
        var sections = new System.Collections.Generic.List<(string Title, System.Collections.Generic.List<string> Items)>();
        foreach (var raw in md.Replace("\r", "").Split('\n'))
        {
            var line = raw.Trim();
            if (line.StartsWith("#"))
            {
                sections.Add((line.TrimStart('#').Trim(), new()));
                continue;
            }
            if (line.Length == 0) continue;
            if (sections.Count == 0) sections.Add(("Changes", new()));
            var item = line.StartsWith("- ") || line.StartsWith("* ") ? line[2..].Trim() : line;
            sections[^1].Items.Add(item.Replace("**", "").Replace("`", ""));
        }
        static bool Alters(string t) => t.Contains("result", StringComparison.OrdinalIgnoreCase) || t.Contains("numer", StringComparison.OrdinalIgnoreCase);
        return sections.Where(s => s.Items.Count > 0).Select(s => new UpdateSection(s.Title, s.Items.ToArray(), Alters(s.Title)))
                       .OrderByDescending(s => s.AltersResults).ToArray();
    }

    public void CloseUpdate() => UpdateOpen = false;
}
