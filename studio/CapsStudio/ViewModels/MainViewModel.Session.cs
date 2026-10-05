using System.Globalization;
using System.Text.Json.Nodes;
using CapsStudio.Interop;

namespace CapsStudio.ViewModels;

/// <summary>The project tree kept across sessions (as a Materials Studio project reopens with its documents and job
/// folders): on quitting, each structure is written to the session folder (its current frame as a LAMMPS data file with
/// the force field and provenance) or, for a trajectory, remembered by its path; Start offers to bring them back, each
/// with its job folders from the job history.</summary>
public sealed partial class MainViewModel
{
    public static string SessionFolder => AppSettings.Override != null
        ? Path.Combine(Path.GetDirectoryName(AppSettings.Override)!, "caps-session")
        : Path.Combine(AppSettings.Folder, "session");
    private static string SessionFile => Path.Combine(SessionFolder, "session.json");

    private JsonArray? _lastSession;
    public bool HasLastSession => _lastSession is { Count: > 0 } && ProjectItems.Count == 0;
    public string LastSessionText
    {
        get
        {
            if (_lastSession is not { Count: > 0 } s) return "";
            var names = s.Select(n => (string?)n?["name"] ?? "").Where(n => n.Length > 0).ToList();
            var jobs = s.Sum(n => (n?["jobs"] as JsonArray)?.Count ?? 0);
            return $"{names.Count} structure{(names.Count == 1 ? "" : "s")}" + (jobs > 0 ? $" · {jobs} job{(jobs == 1 ? "" : "s")}" : "") + " · " +
                   string.Join(", ", names.Take(3)) + (names.Count > 3 ? $" and {names.Count - 3} more" : "");
        }
    }
    public string LastSessionWhen => _lastSessionTime is { } t ? "closed " + t.ToString("d MMM, HH:mm", CultureInfo.InvariantCulture) : "";
    private DateTime? _lastSessionTime;

    /// <summary>Reads the last session's list (Start shows it while nothing is open).</summary>
    public void LoadLastSession()
    {
        _lastSession = null;
        _lastSessionTime = null;
        try
        {
            if (File.Exists(SessionFile) && JsonNode.Parse(File.ReadAllText(SessionFile)) is JsonObject o)
            {
                _lastSession = o["items"] as JsonArray;
                if (DateTime.TryParse((string?)o["saved"], CultureInfo.InvariantCulture, DateTimeStyles.RoundtripKind, out var t)) _lastSessionTime = t;
            }
        }
        catch { _lastSession = null; }
        RaiseSession();
    }

    private void RaiseSession()
    {
        Raise(nameof(HasLastSession));
        Raise(nameof(LastSessionText));
        Raise(nameof(LastSessionWhen));
    }

    /// <summary>Writes the project tree to the session folder (on quitting). Nothing open: the last session is kept.</summary>
    public void SaveSession()
    {
        if (ProjectItems.Count == 0) return;
        try
        {
            StashActive();
            Directory.CreateDirectory(SessionFolder);
            foreach (var f in Directory.EnumerateFiles(SessionFolder, "item*")) File.Delete(f);
            var items = new JsonArray();
            var k = 0;
            foreach (var it in ProjectItems.ToList())
            {
                if (it.Doc.IsDisposed) continue;
                ++k;
                string file;
                string? topology = null;
                var s = it.Doc.Summary();
                if (s.Frames > 1 && it.Doc.Path.Length > 0 && File.Exists(it.Doc.Path))
                {
                    // a trajectory stays where it is (the frames are the file's)
                    file = it.Doc.Path;
                    topology = RecentFiles.Load().FirstOrDefault(r => r.Path == file)?.Topology ?? TopologyFor(file);
                }
                else
                {
                    file = Path.Combine(SessionFolder, $"item{k}.data");
                    it.Doc.Save(file);
                }
                items.Add(new JsonObject
                {
                    ["name"] = it.Name.Replace(" (unsaved)", ""), ["origin"] = it.Origin, ["history"] = it.History, ["force_field"] = it.ForceField,
                    ["file"] = file, ["topology"] = topology, ["frame"] = it.Frame, ["active"] = it == _activeItem, ["builder"] = it.BuildSettings?.DeepClone(),
                    ["jobs"] = new JsonArray(it.Jobs.Select(j => (JsonNode)j.Id).ToArray()),
                });
            }
            File.WriteAllText(SessionFile, new JsonObject { ["saved"] = DateTime.Now.ToString("o"), ["items"] = items }.ToJsonString(new System.Text.Json.JsonSerializerOptions { WriteIndented = true }));
            SaveJobs();
        }
        catch (Exception e) { Status = "Could not keep the session: " + e.Message; }
    }

    /// <summary>Start › Restore: the last session's structures back in the project, with their job folders.</summary>
    public void RestoreSession()
    {
        if (_lastSession is not { Count: > 0 } items) return;
        if (Busy) { Status = "Wait for the run to finish before restoring the session"; return; }
        HookJobs();
        ProjectItem? active = null;
        var missing = new List<string>();
        foreach (var n in items.OfType<JsonObject>())
        {
            var file = (string?)n["file"] ?? "";
            var name = (string?)n["name"] ?? Path.GetFileName(file);
            if (!File.Exists(file)) { missing.Add(name); continue; }
            var topo = (string?)n["topology"] is { Length: > 0 } t && File.Exists(t) ? t : null;
            CapsDocument doc;
            try { doc = CapsDocument.Open(file, topo); }
            catch { missing.Add(name); continue; }
            Show(doc, name);
            var it = _activeItem;
            if (it == null) continue;
            it.Origin = (string?)n["origin"] ?? "";
            it.History = (string?)n["history"] ?? "";
            if (n["builder"] is JsonObject bs) it.BuildSettings = (JsonObject)bs.DeepClone();
            if (it.ForceField.Length == 0) it.ForceField = (string?)n["force_field"] ?? "";
            var frame = (int?)n["frame"] ?? 0;
            if (frame > 0 && frame < Frames) { it.Frame = frame; Frame = frame; }
            if (n["jobs"] is JsonArray ids)
                foreach (var id in ids.Select(x => (string?)x ?? "").Reverse())
                    if (Jobs.FirstOrDefault(j => j.Id == id) is { } job && (job.Item == null || !ProjectItems.Contains(job.Item)))
                    {
                        job.Item = null;   // its structure was closed: the restored one holds it now
                        AttachJob(job, it);
                        RestoredOutputs(job);
                    }
            if ((bool?)n["active"] == true) active = it;
        }
        if (active != null) Activate(active);
        _lastSession = null;
        RaiseSession();
        Status = missing.Count == 0 ? $"Restored {ProjectItems.Count} structure{(ProjectItems.Count == 1 ? "" : "s")} from the last session"
                                    : $"Restored the last session; not found: {string.Join(", ", missing)}";
    }

    /// <summary>A job from an earlier session: its structure and its report (the curves were not kept).</summary>
    private void RestoredOutputs(Job job)
    {
        job.Outputs.Clear();
        if (job.Item is { } item && job.Status is "done" or "stopped" or "cancelled" && job.Kind != "Analyze")
        {
            job.Outputs.Add(new JobOutput($"Structure · {item.Name}", "cube", () => { Activate(item); SetModule(8); }));
            job.Outputs.Add(new JobOutput("Provenance", "history", () => { Activate(item); OpenProvenance(); }));
        }
        job.Outputs.Add(new JobOutput(job.IsFailed ? "Error report" : "Report & log", job.IsFailed ? "alert" : "file", () => ShowJob(job)));
    }

    /// <summary>Start › Forget: the last session is not offered again.</summary>
    public void ForgetSession()
    {
        try { if (File.Exists(SessionFile)) File.Delete(SessionFile); } catch { }
        _lastSession = null;
        RaiseSession();
    }
}
