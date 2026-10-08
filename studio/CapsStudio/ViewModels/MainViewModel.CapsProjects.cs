using System.Collections.ObjectModel;
using System.Globalization;
using System.Text.Json.Nodes;
using CapsStudio.Interop;

namespace CapsStudio.ViewModels;

/// <summary>One session of a project: when it was open and what was built or run in it.</summary>
public sealed record ProjectSession(string When, string Span, string What, bool Now);

/// <summary>CAPS projects (design/boards/ProjectsStart, NewProject, ProjectOpen): Start lists every project CAPS knows,
/// wherever its folder is; a project opens with all its structures (force fields, workflow state, job folders) and the
/// history of its sessions, and is saved as the user works — after each change, on ⌘S, and when CAPS closes.</summary>
public sealed partial class MainViewModel
{
    public ObservableCollection<KnownProject> KnownProjects { get; } = new();
    public bool HasKnownProjects => KnownProjects.Count > 0;
    public string KnownProjectsText => KnownProjects.Count == 1 ? "1 known · wherever it is saved" : $"{KnownProjects.Count} known · wherever they are saved";
    public ObservableCollection<ProjectSession> ProjectSessions { get; } = new();
    public string ProjectSessionsTitle => $"Sessions · {ProjectSessions.Count}";

    private string? _projFile;
    private JsonObject? _projManifest;
    private DateTime _projOpened, _projSavedAt;
    private int _projSession = -1;   // the index of this session in the manifest's sessions
    // the structures as they were when the project opened (name → history): what this session did is measured from it
    private Dictionary<ProjectItem, (string Name, string History)> _projStart = new();
    private readonly Dictionary<ProjectItem, string> _projFiles = new();   // each structure's file under structures/
    private bool _projClosing;

    public bool HasCapsProject => _projFile != null;
    public bool NoCapsProject => _projFile == null;
    public string CapsProjectPath => _projFile ?? "";
    public string CapsProjectName => (string?)_projManifest?["name"] ?? "";
    public string CapsProjectFolderText => _projFile == null ? "" : RecentFiles.Tilde(CapsProjectPath);
    public string CapsProjectSavedText => _projFile == null ? "" : _projSavedAt == default ? "opened " + _projOpened.ToString("HH:mm", CultureInfo.InvariantCulture)
                                                                  : "saved " + _projSavedAt.ToString("HH:mm", CultureInfo.InvariantCulture);

    public string CapsProjectChipTip => $"{CapsProjectFolderText}\n{CapsProjectSavedText} · click to switch, save or close the project";
    public static string SaveKey => OperatingSystem.IsMacOS() ? "⌘S" : "Ctrl+S";
    public string NewProjectSaveNote => $"Saved as you work: after every change, on {SaveKey}, and when CAPS closes.";

    private void RaiseCapsProject()
    {
        foreach (var n in new[] { nameof(HasCapsProject), nameof(NoCapsProject), nameof(CapsProjectPath), nameof(CapsProjectName), nameof(CapsProjectFolderText),
                                  nameof(CapsProjectSavedText), nameof(CapsProjectChipTip), nameof(ExplorerProjectName), nameof(ProjectSessionsTitle) }) Raise(n);
        foreach (var k in KnownProjects) k.Current = _projFile != null && SameFile(k.File, _projFile);
    }

    private static bool SameFile(string a, string b) => string.Equals(Path.GetFullPath(a), Path.GetFullPath(b),
        OperatingSystem.IsLinux() ? StringComparison.Ordinal : StringComparison.OrdinalIgnoreCase);

    // ---------------------------------------------------------------- the list on Start

    /// <summary>Every project CAPS knows, newest first: its name, folder, counts, when last opened; one whose file is not
    /// where it was is shown as not found (Locate / Forget).</summary>
    public void LoadKnownProjects()
    {
        KnownProjects.Clear();
        foreach (var e in ProjectRegistry.Load())
        {
            var k = new KnownProject { File = e.File, Name = e.Name, LastOpened = e.Opened };
            if (!File.Exists(e.File)) { k.Missing = true; k.Counts = e.Opened == DateTime.MinValue ? "" : "last seen " + RecentFiles.Ago(e.Opened); }
            else
            {
                try
                {
                    var m = CapsProjectFile_Read(e.File);
                    k.Name = (string?)m["name"] ?? e.Name;
                    var ns = (m["structures"] as JsonArray)?.Count ?? 0;
                    var nx = (m["sessions"] as JsonArray)?.Count ?? 0;
                    k.Counts = $"{ns} structure{(ns == 1 ? "" : "s")} · {nx} session{(nx == 1 ? "" : "s")}";
                    var thumb = Path.Combine(CapsProjectFile_Folder(e.File), "thumbnail.png");
                    if (File.Exists(thumb)) try { k.Thumbnail = new Avalonia.Media.Imaging.Bitmap(thumb); } catch { }
                }
                catch { k.Counts = "cannot be read"; }
            }
            k.Current = _projFile != null && SameFile(e.File, _projFile);
            k.When = k.Current ? "open now" : e.Opened == DateTime.MinValue ? "" : "opened " + RecentFiles.Ago(e.Opened);
            KnownProjects.Add(k);
        }
        Raise(nameof(HasKnownProjects));
        Raise(nameof(KnownProjectsText));
    }
    private static JsonObject CapsProjectFile_Read(string f) => CapsProjectFile.Read(f);
    private static string CapsProjectFile_Folder(string f) => CapsProjectFile.FolderOf(f);

    /// <summary>A project whose folder moved: its .capsproj in the new place.</summary>
    public string LocateProject(KnownProject k, string newFile)
    {
        if (CapsProjectFile.Resolve(newFile) is not { } resolved) return $"{Path.GetFileName(newFile)} holds no CAPS project (no .capsproj file in it)";
        newFile = resolved;
        try
        {
            var m = CapsProjectFile.Read(newFile);
            ProjectRegistry.Relocate(k.File, newFile, (string?)m["name"] ?? k.Name);
            LoadKnownProjects();
            return $"Found {(string?)m["name"] ?? k.Name} at {RecentFiles.Tilde(CapsProjectFile.FolderOf(newFile))}";
        }
        catch (Exception e) { return $"{Path.GetFileName(newFile)} is not a CAPS project: {e.Message}"; }
    }

    /// <summary>Dropped from the list (its files stay where they are).</summary>
    public string ForgetProject(KnownProject k)
    {
        if (_projFile != null && SameFile(k.File, _projFile)) return "The open project stays in the list: close it first";
        ProjectRegistry.Forget(k.File);
        LoadKnownProjects();
        return $"{k.Name} is no longer listed (its folder was left as it is)";
    }

    // ---------------------------------------------------------------- the New project sheet

    private bool _newProjOpen, _newProjBring = true;
    private string _newProjName = "", _newProjParent = "";
    public bool NewProjectOpen { get => _newProjOpen; set => Set(ref _newProjOpen, value); }
    public string NewProjectName { get => _newProjName; set { if (Set(ref _newProjName, value ?? "")) RaiseNewProject(); } }
    public string NewProjectParent { get => _newProjParent; set { if (Set(ref _newProjParent, value ?? "")) RaiseNewProject(); } }
    public bool NewProjectBring { get => _newProjBring; set => Set(ref _newProjBring, value); }
    public bool NewProjectCanBring => ProjectItems.Count > 0;
    public string NewProjectBringText => ProjectItems.Count == 1 ? "Bring the structure open now into it (with its force field and runs)"
                                       : $"Bring the {ProjectItems.Count} structures open now into it (with their force fields and runs)";
    public string NewProjectParentText => RecentFiles.Tilde(_newProjParent);
    public string NewProjectFolderName => CapsProjectFile.SafeName(_newProjName) + "/";
    public string NewProjectFileName => "  " + CapsProjectFile.SafeName(_newProjName) + CapsProjectFile.Extension;
    public bool NewProjectReady => _newProjName.Trim().Length > 0 && Directory.Exists(_newProjParent);
    public string NewProjectNote => !Directory.Exists(_newProjParent) ? "Choose the folder the project goes in"
        : Directory.Exists(Path.Combine(_newProjParent, CapsProjectFile.SafeName(_newProjName))) ? $"{CapsProjectFile.SafeName(_newProjName)} is there already: the new project gets a folder of its own beside it"
        : "The folder can be moved or copied whole; CAPS finds it again.";
    private void RaiseNewProject()
    {
        foreach (var n in new[] { nameof(NewProjectParentText), nameof(NewProjectFolderName), nameof(NewProjectFileName), nameof(NewProjectReady), nameof(NewProjectNote) }) Raise(n);
    }

    /// <summary>New project…: the sheet, named after the open structure, in the folder the last project went to.</summary>
    public void OpenNewProject()
    {
        var docs = Environment.GetFolderPath(Environment.SpecialFolder.MyDocuments);
        NewProjectParent = Directory.Exists(_settings.ProjectParent) ? _settings.ProjectParent : Directory.Exists(docs) ? docs : Environment.GetFolderPath(Environment.SpecialFolder.UserProfile);
        NewProjectName = _activeItem != null ? CleanName(_activeItem.Name) is var n && n.Contains('.') ? Path.GetFileNameWithoutExtension(n) : CleanName(_activeItem.Name) : "New project";
        NewProjectBring = ProjectItems.Count > 0;
        Raise(nameof(NewProjectCanBring));
        Raise(nameof(NewProjectBringText));
        NewProjectOpen = true;
    }

    public void CreateNewProject()
    {
        if (!NewProjectReady) { Status = NewProjectNote; return; }
        _settings.ProjectParent = _newProjParent;
        _settings.Save();
        NewProjectOpen = false;
        Status = NewCapsProject(_newProjName, _newProjParent, _newProjBring && ProjectItems.Count > 0);
    }

    // ---------------------------------------------------------------- new, open, close

    /// <summary>A new project: PARENT/NAME/NAME.capsproj. The structures open now go into it (bring), or stay in this
    /// session's list (Start › last session) while the new project starts empty.</summary>
    public string NewCapsProject(string name, string parent, bool bring)
    {
        if (Busy) return "Wait for the run to finish before making a project";
        if (string.IsNullOrWhiteSpace(parent) || !Directory.Exists(parent)) return "Choose the folder the project goes in";
        string file;
        try { file = CapsProjectFile.Create(parent, name); }
        catch (Exception e) { return "Could not make the project: " + e.Message; }
        if (_projFile != null)
        {
            SaveCapsProject();
            if (!bring) CloseProjectStructures();
            _projFiles.Clear();   // brought along: written afresh into the new folder
        }
        else if (!bring && ProjectItems.Count > 0) { SaveSession(); CloseProjectStructures(); }
        Attach(file, CapsProjectFile.Read(file));
        SaveCapsProject();
        LoadKnownProjects();
        return $"Made the project {CapsProjectName} in {RecentFiles.Tilde(CapsProjectFile.FolderOf(file))}" +
               (ProjectItems.Count > 0 ? $" with {ProjectItems.Count} structure{(ProjectItems.Count == 1 ? "" : "s")}" : "");
    }

    /// <summary>Opens a project file: every structure comes back (force field, workflow state, job folders, the frame it
    /// was on) and a new session starts. The project open before is saved first; structures open outside a project are
    /// kept as the last session.</summary>
    public string OpenCapsProject(string file)
    {
        if (CapsProjectFile.Resolve(file) is not { } resolved) return $"{Path.GetFileName(file)} holds no CAPS project (no .capsproj file in it)";
        file = resolved;
        if (_projFile != null && SameFile(file, _projFile)) { SetModule(8); return $"{CapsProjectName} is open"; }
        if (Busy) return "Wait for the run to finish before opening another project";
        JsonObject m;
        try { m = CapsProjectFile.Read(file); }
        catch (Exception e) { return $"{Path.GetFileName(file)} could not be opened: {e.Message}"; }
        if (_projFile != null) { SaveCapsProject(closing: true); CloseProjectStructures(); }
        else if (ProjectItems.Count > 0) { SaveSession(); CloseProjectStructures(); }
        _projFiles.Clear();
        HookJobs();
        var folder = CapsProjectFile.FolderOf(file);
        ProjectItem? active = null;
        var missing = new List<string>();
        foreach (var n in (m["structures"] as JsonArray ?? []).OfType<JsonObject>())
        {
            var rel = (string?)n["file"] ?? "";
            var path = Path.IsPathRooted(rel) ? rel : Path.Combine(folder, rel);
            var name = (string?)n["name"] ?? Path.GetFileNameWithoutExtension(path);
            if (!File.Exists(path)) { missing.Add(name); continue; }
            var topo = (string?)n["topology"] is { Length: > 0 } t ? (Path.IsPathRooted(t) ? t : Path.Combine(folder, t)) : null;
            CapsDocument doc;
            try { doc = CapsDocument.Open(path, topo != null && File.Exists(topo) ? topo : null); }
            catch { missing.Add(name); continue; }
            Show(doc, name);
            var it = _activeItem;
            if (it == null) continue;
            var ff = Path.ChangeExtension(path, ".ff.json");
            if (File.Exists(ff)) try { if (doc.FieldLoad(ff)) Field.LoadReport(doc); } catch { }
            it.Origin = (string?)n["origin"] ?? "";
            it.History = (string?)n["history"] ?? "";
            if (n["builder"] is JsonObject bs) it.BuildSettings = (JsonObject)bs.DeepClone();
            if (n["reaction"] is JsonObject rs) it.ReactSettings = (JsonObject)rs.DeepClone();
            it.ForceField = (string?)n["force_field"] ?? it.ForceField;
            if (n["done"] is JsonArray done) { it.Done = new HashSet<string>(done.Select(x => (string?)x ?? "")); _pipeDone.Clear(); foreach (var d in it.Done) _pipeDone.Add(d); RefreshSteps(); }
            var frame = (int?)n["frame"] ?? 0;
            if (frame > 0 && frame < Frames) { it.Frame = frame; Frame = frame; }
            if (!Path.IsPathRooted(rel)) _projFiles[it] = rel;
            if (n["jobs"] is JsonArray ids)
                foreach (var id in ids.Select(x => (string?)x ?? "").Reverse())
                    if (Jobs.FirstOrDefault(j => j.Id == id) is { } job && (job.Item == null || !ProjectItems.Contains(job.Item)))
                    {
                        job.Item = null;
                        AttachJob(job, it);
                        RestoredOutputs(job);
                    }
            if ((bool?)n["active"] == true) active = it;
        }
        if (active != null) Activate(active);
        Attach(file, m);
        LoadKnownProjects();
        if (ProjectItems.Count == 0) SetModule(8);
        return $"Opened {CapsProjectName}: {ProjectItems.Count} structure{(ProjectItems.Count == 1 ? "" : "s")}" +
               (missing.Count > 0 ? $"; not found: {string.Join(", ", missing)}" : "") + _engineNote;
    }

    /// <summary>The project file becomes the open project and a new session starts.</summary>
    private string _engineNote = "";
    public bool EngineNotice { get => _settings.EngineNotice; set { if (_settings.EngineNotice == value) return; _settings.EngineNotice = value; Raise(); _settings.Save(); } }

    private void Attach(string file, JsonObject m)
    {
        _projFile = Path.GetFullPath(file);
        _projManifest = m;
        _projOpened = DateTime.Now;
        _projSavedAt = default;
        var sessions = m["sessions"] as JsonArray ?? new JsonArray();
        m["sessions"] = sessions;
        // the engine each session ran with: a project last worked on with another CAPS version is said so on opening
        var lastEngine = sessions.OfType<JsonObject>().Select(x => (string?)x["engine"]).LastOrDefault(e => e is { Length: > 0 });
        _engineNote = lastEngine != null && lastEngine != ProvGenerator && _settings.EngineNotice
            ? $" · last saved with {lastEngine}, this is {ProvGenerator}: re-run and compare the provenance manifests (Jobs › Provenance) before quoting new numbers"
            : "";
        sessions.Add(new JsonObject { ["opened"] = _projOpened.ToString("o", CultureInfo.InvariantCulture), ["closed"] = null, ["did"] = new JsonArray(), ["engine"] = ProvGenerator });
        _projSession = sessions.Count - 1;
        _projStart = ProjectItems.ToDictionary(it => it, it => (CleanName(it.Name), it.History));   // by structure: two may share a name
        ProjectFolder = Path.Combine(CapsProjectFile.FolderOf(file), "structures");   // the Project home page shows its structures
        ProjectRegistry.Touch(_projFile, CapsProjectName);
        RaiseCapsProject();
        RefreshSessions();
    }

    /// <summary>Saves and closes the project: its structures leave the window and Start shows again.</summary>
    public string CloseCapsProject()
    {
        if (_projFile == null) return "";
        if (Busy) return "Wait for the run to finish before closing the project";
        var name = CapsProjectName;
        SaveCapsProject(closing: true);
        CloseProjectStructures();
        _projFile = null;
        _projManifest = null;
        _projFiles.Clear();
        ProjectFolder = "";
        ProjectSessions.Clear();
        RaiseCapsProject();
        LoadKnownProjects();
        SetModule(8);
        return $"Saved and closed {name}";
    }

    /// <summary>Every structure out of the window, none of them copied aside (the project or the session holds them).</summary>
    private void CloseProjectStructures()
    {
        _projClosing = true;
        try { while (_doc != null) CloseDocument(); }
        finally { _projClosing = false; }
    }

    // ---------------------------------------------------------------- saving

    private static string CleanName(string n) => n.Replace(" (unsaved)", "").Trim();

    /// <summary>Writes every structure into structures/ (a trajectory stays where it is, by its path) with its force field
    /// and provenance, then the project file with this session's line. Nothing open: only the project file.</summary>
    public string SaveCapsProject() => SaveCapsProject(false);

    /// <summary>Save a copy as… (design/boards/ProjectOpen, ⇧⌘S): the project saved, then its whole folder (structures,
    /// runs, exports, the session history) copied to PARENT/NAME with the project file renamed after it; the project open
    /// now stays open, and the copy joins the projects Start lists.</summary>
    public string SaveCapsProjectCopy(string parent, string name)
    {
        if (_projFile == null) return "No project is open";
        if (Busy) return "Wait for the run to finish before copying the project";
        if (string.IsNullOrWhiteSpace(parent) || !Directory.Exists(parent)) return "Choose the folder the copy goes in";
        SaveCapsProject();
        var src = CapsProjectFile.FolderOf(_projFile);
        var safe = CapsProjectFile.SafeName(name);
        var dest = Path.Combine(parent, safe);
        if (Path.GetFullPath(dest).StartsWith(Path.GetFullPath(src) + Path.DirectorySeparatorChar, StringComparison.Ordinal) || SameFile(dest, src))
            return "The copy cannot go inside the project itself";
        if (Directory.Exists(dest) || File.Exists(dest)) return $"{RecentFiles.Tilde(dest)} exists already: choose another name";
        try
        {
            foreach (var dir in Directory.EnumerateDirectories(src, "*", SearchOption.AllDirectories))
                Directory.CreateDirectory(Path.Combine(dest, Path.GetRelativePath(src, dir)));
            Directory.CreateDirectory(dest);
            foreach (var f in Directory.EnumerateFiles(src, "*", SearchOption.AllDirectories))
            {
                if (f.EndsWith(".saving", StringComparison.Ordinal)) continue;
                File.Copy(f, Path.Combine(dest, Path.GetRelativePath(src, f)));
            }
            var oldFile = Path.Combine(dest, Path.GetFileName(_projFile));
            var newFile = Path.Combine(dest, safe + CapsProjectFile.Extension);
            var m = CapsProjectFile.Read(oldFile);
            m["name"] = name.Trim().Length > 0 ? name.Trim() : safe;
            CapsProjectFile.Write(newFile, m);
            if (!SameFile(oldFile, newFile)) File.Delete(oldFile);
            ProjectRegistry.Touch(newFile, (string?)m["name"] ?? safe);
            LoadKnownProjects();
            return $"Saved a copy of {CapsProjectName} as {RecentFiles.Tilde(newFile)} (this project stays open)";
        }
        catch (Exception e) { return "Could not copy the project: " + e.Message; }
    }

    /// <summary>closing: a session in which nothing was built, run or changed is left out of the history.</summary>
    private string SaveCapsProject(bool closing)
    {
        if (_projFile == null || _projManifest == null) return "No project is open";
        _saveQueued = false;
        var folder = CapsProjectFile.FolderOf(_projFile);
        var sdir = Path.Combine(folder, "structures");
        try
        {
            Directory.CreateDirectory(sdir);
            StashActive();
            var before = (_projManifest["structures"] as JsonArray ?? []).OfType<JsonObject>().Select(o => (string?)o["file"] ?? "").Where(f => f.Length > 0 && !Path.IsPathRooted(f)).ToHashSet();
            var items = new JsonArray();
            var used = new HashSet<string>(StringComparer.OrdinalIgnoreCase);
            foreach (var it in ProjectItems.ToList())
            {
                if (it.Doc.IsDisposed) continue;
                var name = CleanName(it.Name);
                if (it.Name != name) it.Name = name;
                if (it == _activeItem && Title.Contains(" (unsaved)", StringComparison.Ordinal)) Title = name;
                string rel;
                string? topology = null;
                var s = it.Doc.Summary();
                if (s.Frames > 1 && it.Doc.Path.Length > 0 && File.Exists(it.Doc.Path) && !it.Doc.Path.StartsWith(sdir, StringComparison.Ordinal))
                {
                    rel = Path.GetFullPath(it.Doc.Path);   // a trajectory stays where it is (the frames are the file's)
                    topology = RecentFiles.Load().FirstOrDefault(r => r.Path == rel)?.Topology ?? TopologyFor(rel);
                }
                else
                {
                    if (!_projFiles.TryGetValue(it, out var mine) || used.Contains(mine))
                    {
                        // a structure named after its file (ps_melt.data) keeps that stem, not ps_melt.data.data
                        var stem = CapsProjectFile.SafeName(StructureExtensions.Contains(Path.GetExtension(name).ToLowerInvariant()) ? Path.GetFileNameWithoutExtension(name) : name);
                        mine = $"structures/{stem}.data";
                        for (var k = 2; used.Contains(mine) || (!before.Contains(mine) && File.Exists(Path.Combine(folder, mine))); ++k) mine = $"structures/{stem} {k}.data";
                        _projFiles[it] = mine;
                    }
                    rel = mine;
                    var path = Path.Combine(folder, rel);
                    it.Doc.Save(path);
                    var ff = Path.ChangeExtension(path, ".ff.json");
                    try { it.Doc.FieldSave(ff); } catch { if (File.Exists(ff)) File.Delete(ff); }   // no force field assigned
                }
                used.Add(rel);
                items.Add(new JsonObject
                {
                    ["name"] = name, ["file"] = rel, ["topology"] = topology, ["origin"] = it.Origin, ["history"] = it.History, ["force_field"] = it.ForceField,
                    ["done"] = new JsonArray((it == _activeItem ? _pipeDone : it.Done).Select(d => (JsonNode)d).ToArray()),
                    ["frame"] = it.Frame, ["active"] = it == _activeItem, ["atoms"] = (double)s.Atoms,
                    ["builder"] = it.BuildSettings?.DeepClone(), ["reaction"] = it.ReactSettings?.DeepClone(),
                    ["jobs"] = new JsonArray(it.Jobs.Select(j => (JsonNode)j.Id).ToArray()),
                });
            }
            // structures taken out of the project: their files go too
            foreach (var gone in before.Where(f => !used.Contains(f)))
                foreach (var f in new[] { gone, Path.ChangeExtension(gone, ".ff.json"), gone + ".provenance.json", gone + ".tags.json" })
                    try { var p = Path.Combine(folder, f); if (File.Exists(p)) File.Delete(p); } catch { }
            _projManifest["structures"] = items;
            _projManifest["saved"] = DateTime.Now.ToString("o", CultureInfo.InvariantCulture);
            if (_projManifest["sessions"] is JsonArray ss && _projSession >= 0 && _projSession < ss.Count && ss[_projSession] is JsonObject cur)
            {
                var did = SessionDid();
                cur["closed"] = DateTime.Now.ToString("o", CultureInfo.InvariantCulture);
                cur["did"] = new JsonArray(did.Select(x => (JsonNode)x).ToArray());
                if (closing && did.Count == 0 && ss.Count > 1) { ss.RemoveAt(_projSession); _projSession = -1; }
            }
            CapsProjectFile.Write(_projFile, _projManifest);
            SaveJobs();
            ProjectThumbnail(folder);
            _projSavedAt = DateTime.Now;
            Raise(nameof(CapsProjectSavedText));
            Raise(nameof(CapsProjectChipTip));
            RefreshSessions();
            return $"Saved {CapsProjectName} ({ProjectItems.Count} structure{(ProjectItems.Count == 1 ? "" : "s")})";
        }
        catch (Exception e) { return $"Could not save {CapsProjectName}: {e.Message}"; }
    }

    /// <summary>What this session did, from the structures now and as they were when it began: new structures, what was
    /// done to the ones there (their history grew), the ones taken out.</summary>
    private List<string> SessionDid()
    {
        var did = new List<string>();
        var now = ProjectItems.Where(it => !it.Doc.IsDisposed).ToList();
        foreach (var it in now)
        {
            var name = CleanName(it.Name);
            var how = it.Origin.Length == 0 || it.Origin.Equals("file", StringComparison.OrdinalIgnoreCase) ? "opened" : it.Origin.ToLowerInvariant();
            if (!_projStart.TryGetValue(it, out var was))
                did.Add(it.History.Length > 0 ? $"{name}: {how} · {it.History}" : $"{name}: {how}");
            else
            {
                if (was.Name != name) did.Add($"{was.Name} renamed {name}");
                if (it.History != was.History && it.History.Length > 0)
                {
                    var h0 = was.History;
                    var grew = it.History.StartsWith(h0, StringComparison.Ordinal) && h0.Length > 0 ? it.History[h0.Length..].TrimStart(' ', '·') : it.History;
                    did.Add($"{name}: {grew}");
                }
            }
        }
        foreach (var (it, was) in _projStart) if (!now.Contains(it)) did.Add($"{was.Name}: taken out");
        return did;
    }

    /// <summary>The active structure's picture as the project's card on Start.</summary>
    private void ProjectThumbnail(string folder)
    {
        var doc = _activeItem?.Doc;
        if (doc == null || doc.IsDisposed) return;
        var opt = new CapsRenderOpts { Width = 440, Height = 240, Supersample = 2, Background = 2, Style = 0, ColourBy = 1, Outlines = 1, DepthCue = 1, ShowCell = 1,
                                       Highlight0 = -1, Highlight1 = -1, Highlight2 = -1, Highlight3 = -1 };
        try { doc.ExportPng(new CapsCamera { Yaw = 0.55, Pitch = 0.40, Zoom = 1.0 }, opt, Path.Combine(folder, "thumbnail.png")); } catch { }
    }

    /// <summary>The sessions, newest first, as the project panel lists them.</summary>
    private void RefreshSessions()
    {
        ProjectSessions.Clear();
        if (_projManifest?["sessions"] is not JsonArray ss) { Raise(nameof(ProjectSessionsTitle)); return; }
        var inv = CultureInfo.InvariantCulture;
        for (var k = ss.Count - 1; k >= 0; --k)
        {
            if (ss[k] is not JsonObject o) continue;
            var now = k == _projSession;
            DateTime.TryParse((string?)o["opened"], inv, DateTimeStyles.RoundtripKind, out var t0);
            var closed = DateTime.TryParse((string?)o["closed"], inv, DateTimeStyles.RoundtripKind, out var t1);
            var did = now ? SessionDid() : (o["did"] as JsonArray ?? []).Select(x => (string?)x ?? "").ToList();
            var when = t0.Date == DateTime.Today ? "Today" : t0.Date == DateTime.Today.AddDays(-1) ? "Yesterday" : t0.ToString("d MMM", inv);
            var span = t0.ToString("HH:mm", inv) + " – " + (now ? "now" : closed ? t1.ToString(t1.Date == t0.Date ? "HH:mm" : "d MMM HH:mm", inv) : "?");
            ProjectSessions.Add(new ProjectSession(when, span, did.Count == 0 ? (now ? "nothing changed yet" : "looked, changed nothing") : string.Join(" · ", did), now));
        }
        Raise(nameof(ProjectSessionsTitle));
    }

    // ---------------------------------------------------------------- saving as the user works

    private bool _saveQueued;
    private Avalonia.Threading.DispatcherTimer? _saveTimer;

    /// <summary>A change to the project (a structure built, opened, renamed or taken out, a run finished): saved a moment
    /// later, once the work is still; never during a run.</summary>
    public void QueueProjectSave()
    {
        if (_projFile == null || _projClosing) return;
        _saveQueued = true;
        _saveTimer ??= new Avalonia.Threading.DispatcherTimer(TimeSpan.FromSeconds(2), Avalonia.Threading.DispatcherPriority.Background, (_, _) =>
        {
            if (Busy) return;   // tried again on the next tick
            _saveTimer!.Stop();
            if (_saveQueued && _projFile != null) Status = SaveCapsProject();
        });
        _saveTimer.Stop();
        _saveTimer.Start();
    }

    private void OnProjectItemsChanged(object? s, System.Collections.Specialized.NotifyCollectionChangedEventArgs e)
    {
        foreach (var it in e.NewItems?.OfType<ProjectItem>() ?? []) it.PropertyChanged += OnProjectItemChanged;
        foreach (var it in e.OldItems?.OfType<ProjectItem>() ?? []) { it.PropertyChanged -= OnProjectItemChanged; _projFiles.Remove(it); }
        QueueProjectSave();
    }
    private void OnProjectItemChanged(object? s, System.ComponentModel.PropertyChangedEventArgs e)
    {
        if (e.PropertyName is nameof(ProjectItem.Name) or nameof(ProjectItem.History) or nameof(ProjectItem.ForceField)) QueueProjectSave();
    }

    /// <summary>The project's own name changed (its folder and file keep theirs).</summary>
    private void RenameCapsProject(string name)
    {
        if (_projFile == null || _projManifest == null || name.Length == 0) return;
        _projManifest["name"] = name;
        ProjectRegistry.Touch(_projFile, name);
        Status = SaveCapsProject();
        RaiseCapsProject();
        LoadKnownProjects();
    }

    /// <summary>On quitting: the open project saved (else the session, as before).</summary>
    public void SaveOnClose()
    {
        if (_projFile != null) SaveCapsProject(closing: true);
        else SaveSession();
    }
}
