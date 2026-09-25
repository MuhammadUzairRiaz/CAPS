using System.Collections.ObjectModel;
using System.Globalization;
using System.IO.Compression;
using System.Text.Json.Nodes;
using Avalonia.Media.Imaging;
using CapsStudio.Interop;

namespace CapsStudio.ViewModels;

/// <summary>A structure file of the project folder.</summary>
public sealed class ProjectDoc : ObservableObject
{
    public string Path { get; init; } = "";
    public string Name => System.IO.Path.GetFileName(Path);
    private string _subtitle = "reading…", _density = "—", _atoms = "—", _status = "", _last = "";
    private Bitmap? _thumb;
    private bool _selected;
    public string Subtitle { get => _subtitle; set => Set(ref _subtitle, value); }
    public string Density { get => _density; set => Set(ref _density, value); }
    public string Atoms { get => _atoms; set => Set(ref _atoms, value); }
    public string Status { get => _status; set => Set(ref _status, value); }
    public string LastStep { get => _last; set => Set(ref _last, value); }
    public Bitmap? Thumbnail { get => _thumb; set => Set(ref _thumb, value); }
    public bool Selected { get => _selected; set => Set(ref _selected, value); }
    public JsonNode? Manifest { get; set; }
    public bool HasProvenance => Manifest != null;
    public int StatusLevel { get; set; }   // 0 none, 1 in progress, 2 done
}

/// <summary>Project home (design/boards/ProjectHome): the structures of a folder with thumbnails, a table of what each
/// holds and how far it went, and a methods section written from the selected structure's provenance — replicas that
/// differ only in their seeds counted — with numbered references, ready to copy.</summary>
public sealed partial class MainViewModel
{
    public bool IsProject => _module == 42;
    public ObservableCollection<ProjectDoc> ProjectDocs { get; } = new();
    public ObservableCollection<string> ProjectRefs { get; } = new();
    private string _projectFolder = "", _projectMethods = "", _projectMethodsFor = "";
    public string ProjectFolder { get => _projectFolder; private set { if (Set(ref _projectFolder, value)) { Raise(nameof(ProjectName)); Raise(nameof(ProjectFolderText)); } } }
    public string ProjectName => _projectFolder.Length == 0 ? "No project" : System.IO.Path.GetFileName(_projectFolder.TrimEnd('/', '\\'));
    public string ProjectFolderText => _projectFolder.Length == 0 ? "" : RecentFiles.Tilde(_projectFolder);
    public string ProjectMethods { get => _projectMethods; private set { if (Set(ref _projectMethods, value)) Raise(nameof(ProjectHasMethods)); } }
    public bool ProjectHasMethods => _projectMethods.Length > 0;
    public string ProjectMethodsFor { get => _projectMethodsFor; private set => Set(ref _projectMethodsFor, value); }
    public string ProjectSubtitle => $"Project · {ProjectDocs.Count} structure{(ProjectDocs.Count == 1 ? "" : "s")} · {ProjectDocs.Count(d => d.HasProvenance)} with provenance · {Jobs.Count} jobs this session";
    private static readonly string[] StructureExtensions = [".data", ".pdb", ".gro", ".xyz", ".mol2", ".cif", ".lammpstrj", ".dump", ".extxyz"];
    private int _projectTicket;

    /// <summary>The folder of the open structure (or the given one) as a project.</summary>
    public void OpenProject(string? folder = null)
    {
        folder ??= _doc?.Path is { Length: > 0 } p && File.Exists(p) ? System.IO.Path.GetDirectoryName(System.IO.Path.GetFullPath(p)) : _projectFolder;
        SetModule(42);
        if (string.IsNullOrEmpty(folder) || !Directory.Exists(folder)) { ProjectFolder = ""; ProjectDocs.Clear(); Raise(nameof(ProjectSubtitle)); return; }
        ProjectFolder = folder;
        LoadProject();
    }

    private void LoadProject()
    {
        var ticket = ++_projectTicket;
        ProjectDocs.Clear();
        var files = Directory.EnumerateFiles(_projectFolder)
            .Where(f => StructureExtensions.Contains(System.IO.Path.GetExtension(f).ToLowerInvariant()))
            .OrderBy(f => f, StringComparer.OrdinalIgnoreCase).Take(24).ToList();
        foreach (var f in files) ProjectDocs.Add(new ProjectDoc { Path = f });
        Raise(nameof(ProjectSubtitle));
        var inv = CultureInfo.InvariantCulture;
        foreach (var d in ProjectDocs.ToList())
        {
            // the provenance first (cheap), then the structure for its numbers and thumbnail
            try
            {
                var j = JsonNode.Parse(CapsDocument.ProvenanceFile(d.Path));
                if (j?["ok"]?.GetValue<bool>() == true) d.Manifest = j;
            }
            catch { }
            var steps = d.Manifest?["steps"] as JsonArray;
            var last = steps?.LastOrDefault()?["engine"]?.GetValue<string>() ?? "";
            d.LastStep = last;
            (d.Status, d.StatusLevel) = last switch
            {
                "" => ("no provenance", 0),
                _ when last.StartsWith("equilibrate.", StringComparison.Ordinal) => ("equilibrated", 2),
                _ when last.StartsWith("dynamics.", StringComparison.Ordinal) => ("simulated", 2),
                _ when last.StartsWith("relax.", StringComparison.Ordinal) => ("minimised", 1),
                "field.assign" => ("typed", 1),
                _ when last.StartsWith("io.", StringComparison.Ordinal) => ("read", 0),
                _ => ("built", 1),
            };
            var path = d.Path;
            if (new FileInfo(path).Length > 60_000_000) { d.Subtitle = "large file · open it to see it"; continue; }
            Task.Run(() =>
            {
                try
                {
                    using var doc = CapsDocument.Open(path, TopologyFor(path) is { Length: > 0 } t ? t : null);
                    var s = doc.Summary();
                    const int w = 420, h = 240;
                    var rgba = new byte[w * h * 4];
                    var opt = ViewOptions(w, h, 2);
                    opt.Style = 3;   // no hydrogens: a thumbnail reads better
                    opt.Background = _viewBackground;
                    opt.Highlight0 = opt.Highlight1 = opt.Highlight2 = opt.Highlight3 = -1;
                    opt.Focus = 0;
                    doc.Render(new CapsCamera { Yaw = 0.5, Pitch = 0.4, Zoom = 1.0 }, opt, rgba);
                    return (Summary: (CapsSummary?)s, Rgba: rgba, W: w, H: h);
                }
                catch { return (Summary: (CapsSummary?)null, Rgba: Array.Empty<byte>(), W: 0, H: 0); }
            }).ContinueWith(t => Avalonia.Threading.Dispatcher.UIThread.Post(() =>
            {
                if (ticket != _projectTicket) return;
                var (s, rgba, w, h) = t.Result;
                if (s is not { } sum) { d.Subtitle = "cannot read"; return; }
                d.Atoms = sum.Atoms.ToString("N0", inv);
                d.Density = sum.CellValid != 0 ? sum.Density.ToString("0.000", inv) : "—";
                d.Subtitle = $"{sum.Atoms:N0} atoms · {(d.HasProvenance ? d.Status : "no provenance")}";
                if (rgba.Length > 0) d.Thumbnail = ToBitmap(rgba, w, h);
            }));
        }
        Raise(nameof(ProjectSubtitle));
        SelectProjectDoc(ProjectDocs.FirstOrDefault(d => _doc != null && d.Path == System.IO.Path.GetFullPath(_doc.Path)) ?? ProjectDocs.FirstOrDefault(d => d.HasProvenance) ?? ProjectDocs.FirstOrDefault());
    }

    /// <summary>The methods section from this structure's provenance; the others identical but for their seeds are replicas.</summary>
    public void SelectProjectDoc(ProjectDoc? d)
    {
        foreach (var x in ProjectDocs) x.Selected = x == d;
        ProjectRefs.Clear();
        if (d?.Manifest == null) { ProjectMethods = ""; ProjectMethodsFor = d == null ? "" : d.Name + " has no provenance (FILE.provenance.json) yet"; return; }
        var a = d.Manifest.ToJsonString();
        var replicas = new JsonArray();
        foreach (var o in ProjectDocs.Where(o => o != d && o.Manifest != null))
        {
            try
            {
                var c = JsonNode.Parse(CapsDocument.ProvenanceCompare(a, o.Manifest!.ToJsonString()))!;
                var rows = (JsonArray)c["rows"]!;
                var notes = (JsonArray)c["notes"]!;
                if (rows.Count > 0 && notes.Count == 0 && rows.All(r => r!["key"]!.GetValue<string>() == "rng")) replicas.Add(JsonNode.Parse(o.Manifest.ToJsonString()));
            }
            catch { }
        }
        var m = JsonNode.Parse(CapsDocument.MethodsText(a, replicas.Count > 0 ? replicas.ToJsonString() : null))!;
        ProjectMethods = m["text"]?.GetValue<string>() ?? "";
        foreach (var (r, k) in ((JsonArray)m["refs"]!).Select((r, k) => (r, k))) ProjectRefs.Add($"[{k + 1}] {r!.GetValue<string>()}");
        ProjectMethodsFor = d.Name + (replicas.Count > 0 ? $" · {replicas.Count + 1} replicas" : "");
    }

    public string ProjectMethodsWithRefs => _projectMethods + "\n\n" + string.Join("\n", ProjectRefs);

    public string ProjectBibtex()
    {
        var d = ProjectDocs.FirstOrDefault(x => x.Selected);
        return d?.Manifest == null ? "" : CapsDocument.ProvenanceBibtex(d.Manifest.ToJsonString());
    }

    /// <summary>The project's structures, their provenance, the methods text and its BibTeX in one zip.</summary>
    public string ShareProject(string zipPath)
    {
        if (File.Exists(zipPath)) File.Delete(zipPath);
        using var zip = System.IO.Compression.ZipFile.Open(zipPath, System.IO.Compression.ZipArchiveMode.Create);
        var n = 0;
        foreach (var d in ProjectDocs)
        {
            zip.CreateEntryFromFile(d.Path, d.Name);
            n++;
            var side = d.Path + ".provenance.json";
            if (File.Exists(side)) zip.CreateEntryFromFile(side, d.Name + ".provenance.json");
        }
        if (_projectMethods.Length > 0)
        {
            using (var w = new StreamWriter(zip.CreateEntry("methods.txt").Open())) w.Write(ProjectMethodsWithRefs);
            using (var w = new StreamWriter(zip.CreateEntry("references.bib").Open())) w.Write(ProjectBibtex());
        }
        return $"Wrote {System.IO.Path.GetFileName(zipPath)} · {n} structures with their provenance" + (_projectMethods.Length > 0 ? ", methods.txt and references.bib" : "");
    }
}
