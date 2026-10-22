using System;
using System.Collections.ObjectModel;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Net.Http;
using System.Text.Json.Nodes;
using System.Threading.Tasks;

namespace CapsStudio.ViewModels;

/// <summary>A line of What's new: what changed, and the page that shows it (Open).</summary>
public sealed record WhatsNewItem(string Title, string Detail, string Command, bool CanOpen)
{
    public bool HasDetail => Detail.Length > 0;
}

/// <summary>Update notices and What's new.
/// Once a day at most, after the window is up, the Studio fetches the latest release's release.json (one GET of a public
/// file, no cookies, no identifiers, nothing else sent; a few seconds at most, silent when offline). When that release is
/// newer than this build (a higher version; or, with equal versions, another commit with a later date) a banner at the
/// top offers it: What's new, Download (the installer for this machine, in the browser), Later, Skip this update.
/// After an update, the first start of a build whose bundled data/release.json names another commit than the one last
/// seen opens What's new with that build's highlights, each with Open for its page. CAPS_NO_UPDATE_CHECK=1 turns the
/// check off; CAPS_UPDATE_FEED points it at a local file (tests) or another URL.</summary>
public partial class MainViewModel
{
    private bool _updateChecking, _bannerOpen, _laterThisSession, _whatsNewOpen, _whatsNewRemote;
    private string _updateText = "", _bannerText = "", _whatsNewTitle = "", _whatsNewSubtitle = "", _whatsNewNote = "";
    private ReleaseManifest? _remoteRelease;
    private BuildIdentity? _thisBuild;

    public bool UpdateChecking { get => _updateChecking; private set => Set(ref _updateChecking, value); }
    /// <summary>The result of the last check, for Settings.</summary>
    public string UpdateText { get => _updateText; private set => Set(ref _updateText, value); }
    public bool UpdateBannerOpen { get => _bannerOpen; private set => Set(ref _bannerOpen, value); }
    public string UpdateBannerText { get => _bannerText; private set { Set(ref _bannerText, value); Raise(nameof(UpdateBannerLead)); Raise(nameof(UpdateBannerDetail)); } }
    /// <summary>The banner's text in two parts: "CAPS update available" and " · 5 Nov 2026 · November update".</summary>
    public string UpdateBannerLead => _bannerText.Split(" · ")[0];
    public string UpdateBannerDetail => _bannerText.Length > UpdateBannerLead.Length ? _bannerText[UpdateBannerLead.Length..] : "";
    public ReleaseManifest? RemoteRelease => _remoteRelease;
    public string UpdateDownloadUrl => _remoteRelease == null ? UpdateRules.ReleasesPage : UpdateRules.DownloadUrl(_remoteRelease, UpdateRules.PlatformKey());
    public string UpdateNotesUrl => (_whatsNewRemote ? _remoteRelease : BundledRelease ?? _remoteRelease)?.NotesUrl is { Length: > 0 } u ? u : UpdateRules.ReleasesPage;

    public bool WhatsNewOpen { get => _whatsNewOpen; set => Set(ref _whatsNewOpen, value); }
    /// <summary>The panel shows a release not installed yet (from the banner): Download and Skip are offered.</summary>
    public bool WhatsNewRemote { get => _whatsNewRemote; private set { Set(ref _whatsNewRemote, value); Raise(nameof(UpdateNotesUrl)); } }
    public string WhatsNewTitle { get => _whatsNewTitle; private set => Set(ref _whatsNewTitle, value); }
    public string WhatsNewSubtitle { get => _whatsNewSubtitle; private set => Set(ref _whatsNewSubtitle, value); }
    public string WhatsNewNote { get => _whatsNewNote; private set { Set(ref _whatsNewNote, value); Raise(nameof(HasWhatsNewNote)); } }
    public bool HasWhatsNewNote => _whatsNewNote.Length > 0;
    public ObservableCollection<WhatsNewItem> WhatsNewItems { get; } = new();

    public bool SetCheckUpdates
    {
        get => _settings.CheckUpdates;
        set { if (_settings.CheckUpdates == value) return; _settings.CheckUpdates = value; Raise(); Changed(value ? "Update checks on" : "Update checks off"); }
    }
    public bool UpdateChecksDisabledByEnv => UpdateRules.Disabled();

    /// <summary>This build as the native core reports it (version, commit, the commit's date).</summary>
    public BuildIdentity ThisBuild
    {
        get
        {
            if (_thisBuild != null) return _thisBuild;
            try
            {
                var j = JsonNode.Parse(Interop.CapsDocument.BuildInfo())!;
                _thisBuild = new((string?)j["version"] ?? "0.1.0", (string?)j["commit"] ?? "unknown", (string?)j["date"] ?? "unknown");
            }
            catch { _thisBuild = new("0.1.0", "unknown", "unknown"); }
            return _thisBuild;
        }
    }
    /// <summary>For tests: what this build reports itself as.</summary>
    public void SetThisBuildForTest(BuildIdentity? b) => _thisBuild = b;

    /// <summary>The release.json bundled with this app (absent in a build from source).</summary>
    public ReleaseManifest? BundledRelease => _bundledOverride ?? ReleaseManifest.TryLoad(Paths.Data is { } d ? Path.Combine(d, "release.json") : null);
    private ReleaseManifest? _bundledOverride;
    public void SetBundledReleaseForTest(ReleaseManifest? m) => _bundledOverride = m;

    /// <summary>After the window is up: What's new after an update, then the daily check (or yesterday's answer).</summary>
    public void StartupUpdateTasks()
    {
        if (!UpdateRules.NetworkAllowed) return;   // the tool modes (self-test, screenshots, figures)
        MaybeShowWhatsNewAfterUpdate();
        if (!_settings.CheckUpdates || UpdateRules.Disabled()) return;
        var last = UpdateRules.ParseDate(_settings.LastUpdateCheck);
        if (last is { } t && DateTimeOffset.UtcNow - t < TimeSpan.FromHours(24) && _settings.LatestRelease.Length > 0)
        {
            try { ApplyRelease(ReleaseManifest.Parse(_settings.LatestRelease), manual: false); } catch { /* stale cache: next check replaces it */ }
            return;
        }
        _ = CheckForUpdates(manual: false);
    }

    /// <summary>The first start of a build whose bundled manifest names another commit than the last one seen opens What's
    /// new (not on a brand-new install, where the tour introduces the Studio instead).</summary>
    public bool MaybeShowWhatsNewAfterUpdate()
    {
        if (BundledRelease is not { } m || UpdateRules.SameCommit(m.Commit, _settings.LastSeenCommit)) return false;
        var firstEver = _settings.LastSeenCommit.Length == 0 && !_settings.TourDone;
        _settings.LastSeenCommit = m.Commit;
        _settings.Save();
        if (firstEver) return false;
        ShowWhatsNew(m, remote: false);
        return true;
    }

    /// <summary>Help › What's new: this build's highlights (or, in a build from source, the latest release's if known).</summary>
    public void OpenWhatsNew()
    {
        if (BundledRelease is { } m) ShowWhatsNew(m, remote: false);
        else if (_remoteRelease is { } r) ShowWhatsNew(r, remote: UpdateRules.Decide(ThisBuild, r) == UpdateVerdict.Update);
        else
        {
            WhatsNewItems.Clear();
            WhatsNewRemote = false;
            WhatsNewTitle = "What's new in CAPS";
            WhatsNewSubtitle = $"CAPS {ThisBuild.Version} · build {ThisBuild.Commit}";
            WhatsNewNote = "This build was made from source and carries no release notes. The release notes on GitHub list every change.";
            WhatsNewOpen = true;
        }
    }

    /// <summary>The banner's What's new: the highlights of the release on offer, before updating.</summary>
    public void OpenRemoteWhatsNew() { if (_remoteRelease is { } r) ShowWhatsNew(r, remote: true); }

    public void ShowWhatsNew(ReleaseManifest m, bool remote)
    {
        WhatsNewItems.Clear();
        foreach (var h in m.Highlights)
        {
            var cmd = PageCommand(h.Page);
            WhatsNewItems.Add(new WhatsNewItem(h.Title, h.Detail, cmd, !remote && cmd.Length > 0 && HasCommand(cmd)));
        }
        WhatsNewRemote = remote;
        var version = UpdateRules.CompareVersions(m.Version, ThisBuild.Version) != 0
            ? (remote ? $" · version {m.Version} (you have {ThisBuild.Version})" : $" · version {m.Version}") : "";
        WhatsNewTitle = remote ? "A newer CAPS is ready" : "What's new in CAPS";
        WhatsNewSubtitle = $"{(m.Title.Length > 0 ? m.Title : "CAPS " + m.Version)} · {m.DateText}{version} · build {m.ShortCommit}";
        WhatsNewNote = remote
            ? "Your projects keep the CAPS version they were made with. After updating, re-run a project and compare the two provenance files before quoting new numbers."
            : m.Highlights.Count == 0 ? "This release lists no highlights; the release notes on GitHub list every change." : "";
        WhatsNewOpen = true;
    }

    /// <summary>"semicrystalline" → the palette's "module.open semicrystalline"; a full command id stays as it is.</summary>
    public static string PageCommand(string page)
    {
        page = page.Trim();
        if (page.Length == 0) return "";
        return page.Contains('.') || page.Contains(' ') ? page : "module.open " + page.ToLowerInvariant();
    }
    public bool HasCommand(string id)
    {
        if (!_modelCommands) { _modelCommands = true; AddModelCommands(); }
        return _commands.Any(c => c.Id == id);
    }
    public void OpenWhatsNewItem(WhatsNewItem item)
    {
        WhatsNewOpen = false;
        if (!RunCommand(item.Command)) Status = $"{item.Title}: that page is not available now";
    }
    public void CloseWhatsNew() => WhatsNewOpen = false;

    /// <summary>Fetch the latest release's manifest and compare. manual: Help › Check for updates now (says the result,
    /// shows a skipped update again); otherwise silent on any failure.</summary>
    public async Task CheckForUpdates(bool manual = true)
    {
        if (_updateChecking) return;
        if (UpdateRules.Disabled())
        {
            UpdateText = "Update checks are off on this computer (CAPS_NO_UPDATE_CHECK)";
            if (manual) Status = UpdateText;
            return;
        }
        UpdateChecking = true;
        if (manual) UpdateText = "Checking…";
        try
        {
            var json = await Task.Run(FetchManifest);
            var m = ReleaseManifest.Parse(json);
            _settings.LastUpdateCheck = DateTimeOffset.UtcNow.ToString("o", CultureInfo.InvariantCulture);
            _settings.LatestRelease = json.Length < 64_000 ? json : "";
            _settings.Save();
            ApplyRelease(m, manual);
        }
        catch (Exception e) when (manual)
        {
            UpdateText = e is HttpRequestException or TaskCanceledException ? "Could not reach GitHub: offline, or no release manifest yet"
                                                                            : "Could not check: " + e.Message;
            Status = UpdateText;
        }
        catch { /* automatic: silent */ }
        finally { UpdateChecking = false; }
    }

    /// <summary>What a manifest means for this build: the banner when it is an update (unless skipped, or Later this
    /// session, for automatic checks).</summary>
    public UpdateVerdict ApplyRelease(ReleaseManifest m, bool manual)
    {
        _remoteRelease = m;
        Raise(nameof(UpdateDownloadUrl));
        var verdict = UpdateRules.Decide(ThisBuild, m);
        var when = DateTime.Now.ToString("HH:mm", CultureInfo.InvariantCulture);
        if (verdict != UpdateVerdict.Update)
        {
            UpdateBannerOpen = false;
            UpdateText = verdict == UpdateVerdict.Unknown ? $"Cannot tell whether the latest release is newer (this build has no date) · checked {when}"
                                                          : $"You have the latest CAPS · checked {when}";
            if (manual)
            {
                Status = UpdateText;
                Notify(new Notice { Key = "update", Severity = "check", Icon = "check", Title = verdict == UpdateVerdict.Unknown ? "Update check" : "You have the latest CAPS",
                                    Body = $"Latest release: {(m.Title.Length > 0 ? m.Title + " · " : "")}{m.DateText} · build {m.ShortCommit}. This build: {ThisBuild.Commit}.", Primary = "OK" });
            }
            return verdict;
        }
        var version = UpdateRules.CompareVersions(m.Version, ThisBuild.Version) != 0 ? $" · version {m.Version}" : "";
        UpdateBannerText = $"CAPS update available{version} · {m.DateText}{(m.Title.Length > 0 ? " · " + m.Title : "")}";
        UpdateText = $"An update is available{version} · {m.DateText} · checked {when}";
        var skipped = UpdateRules.SameCommit(m.Commit, _settings.SkippedUpdateCommit);
        if (manual) { _laterThisSession = false; UpdateBannerOpen = true; }
        else UpdateBannerOpen = !skipped && !_laterThisSession;
        return verdict;
    }

    /// <summary>Later: hidden until the next start.</summary>
    public void UpdateLater() { _laterThisSession = true; UpdateBannerOpen = false; }
    /// <summary>Skip this update: this release's commit is remembered and its banner not shown again.</summary>
    public void SkipThisUpdate()
    {
        if (_remoteRelease is { } m) { _settings.SkippedUpdateCommit = m.Commit; _settings.Save(); }
        UpdateBannerOpen = false;
        WhatsNewOpen = false;
        Status = "This update will not be offered again; Help › Check for updates now still shows it";
    }
    public bool IsSkipped(string commit) => UpdateRules.SameCommit(commit, _settings.SkippedUpdateCommit);

    private static string FetchManifest()
    {
        var feed = Environment.GetEnvironmentVariable("CAPS_UPDATE_FEED");
        if (feed is { Length: > 0 } && !feed.StartsWith("http", StringComparison.OrdinalIgnoreCase)) return File.ReadAllText(feed);
        if (!UpdateRules.NetworkAllowed) throw new InvalidOperationException("no network in this mode");
        var url = feed is { Length: > 0 } ? feed : UpdateRules.ManifestUrl;
        using var handler = new HttpClientHandler { UseCookies = false, AllowAutoRedirect = true, MaxAutomaticRedirections = 5, UseDefaultCredentials = false };
        using var http = new HttpClient(handler) { Timeout = TimeSpan.FromSeconds(6), MaxResponseContentBufferSize = 1 << 20 };
        http.DefaultRequestHeaders.UserAgent.ParseAdd("CAPS");   // GitHub refuses requests without one; nothing else identifies the user
        return http.GetStringAsync(url).GetAwaiter().GetResult();
    }
}
