using System.Collections.ObjectModel;
using System.Globalization;
using System.Linq;
using CapsStudio.Interop;

namespace CapsStudio.ViewModels;

/// <summary>One structure of the project (as a Materials Studio document): a molecule, a grown polymer cell, a packed
/// box, a crystal, a file. Every module (Force field, Minimise, Equilibrate, Dynamics, React, Analyze, Export) works on
/// the active one; the others keep their state (camera, frame, what was done to them) until they are picked again.</summary>
public sealed class ProjectItem : ObservableObject
{
    public ProjectItem(CapsDocument doc, string name)
    {
        Doc = doc;
        _name = name;
        Jobs.CollectionChanged += (_, _) => Raise(nameof(HasJobs));
    }
    public bool HasJobs => Jobs.Count > 0;
    public CapsDocument Doc { get; }

    private string _name;
    public string Name { get => _name; set => Set(ref _name, value); }
    private string _origin = "";
    /// <summary>How it was made: "Polymer cell" (Grow), "Packing" (Pack), a builder, a file.</summary>
    public string Origin { get => _origin; set { if (Set(ref _origin, value)) { Raise(nameof(IconKind)); Raise(nameof(EditTip)); } } }
    private bool _active;
    public bool Active { get => _active; set => Set(ref _active, value); }
    private string _detail = "";
    public string Detail { get => _detail; set => Set(ref _detail, value); }
    private string _forceField = "";
    public string ForceField { get => _forceField; set { if (Set(ref _forceField, value)) Raise(nameof(HasForceField)); } }
    public bool HasForceField => _forceField.Length > 0;
    private string _history = "";
    /// <summary>What has been done to it: "minimised · equilibrated · MD".</summary>
    public string History { get => _history; set { if (Set(ref _history, value)) Raise(nameof(HasHistory)); } }
    public bool HasHistory => _history.Length > 0;
    /// <summary>Its job folders, newest first (Materials Studio's project tree).</summary>
    public ObservableCollection<Job> Jobs { get; } = new();
    private bool _expanded = true;
    public bool Expanded { get => _expanded; set => Set(ref _expanded, value); }
    public string IconKind => _origin switch { "Polymer cell" or "Blend" => "grow", "Packing" or "Solvation" => "pack", "Molecule builder" => "atom", "File" => "file", _ => "cube" };
    // Delete asks twice: the first click arms it (the row says so), a second within a few seconds removes the structure
    private bool _deleteArmed;
    public bool DeleteArmed { get => _deleteArmed; set { if (Set(ref _deleteArmed, value)) Raise(nameof(DeleteTip)); } }
    public string DeleteTip => _deleteArmed ? "Click again to remove it from the project (a saved file on disk stays)" : "Delete: remove from the project (asks once more)";
    public string EditTip => _origin.Length > 0 && _origin != "File" ? $"Edit in the {_origin}" : "Edit with the Studio's builder tools";

    // the state kept while another structure is active
    internal CapsCamera Camera = new() { Yaw = 0.55, Pitch = 0.40, Zoom = 1.0 };
    internal int Frame;
    internal HashSet<string> Done = new();
    internal string Build = "";
    internal bool GrownUnsaved, AutoFf;
}

public sealed partial class MainViewModel
{
    public ObservableCollection<ProjectItem> ProjectItems { get; } = new();
    private ProjectItem? _activeItem;

    /// <summary>The structure the modules work on; picking another one (tab, project list or a module's "Structure"
    /// box) makes it active.</summary>
    public ProjectItem? ActiveItem
    {
        get => _activeItem;
        set { if (value != null && value != _activeItem) Activate(value); }
    }
    public bool HasSeveralStructures => ProjectItems.Count > 1;
    public string ProjectCountText => ProjectItems.Count == 1 ? "1 structure" : $"{ProjectItems.Count} structures";

    /// <summary>A new structure (built, grown, packed, opened): added to the project and made active; the one that was
    /// active stays in the project.</summary>
    private ProjectItem AddProjectItem(CapsDocument doc, string title)
    {
        StashActive();
        var it = new ProjectItem(doc, title)
        {
            Origin = title.StartsWith("packed_", StringComparison.Ordinal) ? "Packing" : MakerOf(doc).Label,
        };
        ProjectItems.Add(it);
        foreach (var p in ProjectItems) p.Active = p == it;
        _activeItem = it;
        RaiseProject();
        return it;
    }

    /// <summary>What made a structure, from the first step of its provenance: the builder's name, its page, and that step.</summary>
    private static (string Label, int Module, System.Text.Json.Nodes.JsonNode? Step) MakerOf(CapsDocument doc)
    {
        System.Text.Json.Nodes.JsonNode? step = null;
        try
        {
            if (System.Text.Json.Nodes.JsonNode.Parse(doc.Provenance())?["steps"] is System.Text.Json.Nodes.JsonArray st)
                step = st.FirstOrDefault(x => x?["engine"]?.GetValue<string>() is string e && e != "scratch");
        }
        catch { }
        var engine = step?["engine"]?.GetValue<string>() ?? "";
        (string, int) m = engine switch
        {
            "grow.trials" => ("Polymer cell", 0),
            "grow.blend" => ("Blend", 16),
            "pack.lbfgs" or "pack.insert" => ("Packing", 5),
            "chem.build" => ("Molecule builder", 9),
            "crystal.build" => ("Crystal builder", 29),
            "nano.build" or "nano.embed" or "nano.pore" => ("Nanostructure builder", 15),
            "surface.build" or "interface.build" or "build.stack" => ("Surface builder", 14),
            "bio.peptide" => ("Biomolecule builder", 30),
            "solvate.pack" => ("Solvation", 31),
            "cg.build" or "cg.kremer_grest" or "cg.martini_melt" => ("Coarse-grained builder", 44),
            "io.read" or "io.import" => ("File", 8),
            _ => ("", 8),
        };
        return (m.Item1, m.Item2, step);
    }

    /// <summary>Edit (the project tree): the structure made active and the builder that made it opened — the molecule
    /// builder with its SMILES; a file, a copy or an edited structure: the Studio with its builder tools.</summary>
    public void EditProjectItem(ProjectItem it)
    {
        if (!Idle) { Status = "Wait for the run to finish before editing"; return; }
        Activate(it);
        var (label, module, step) = MakerOf(it.Doc);
        if (label.Length > 0 && it.Origin.Length == 0) it.Origin = label;
        switch (module)
        {
            case 9:
                var smi = (step?["params"] as System.Text.Json.Nodes.JsonObject)?["smiles"]?.GetValue<string>()
                          ?? (step?["params"] as System.Text.Json.Nodes.JsonArray)?.FirstOrDefault(p => p?[0]?.GetValue<string>() == "smiles")?[1]?.GetValue<string>();
                OpenBuilder(smi);
                break;
            case 0: SetModule(0); break;
            case 5: SetModule(5); break;
            case 16: OpenBlend(); break;
            case 29: OpenCrystal(); break;
            case 15: OpenNano(); break;
            case 14: OpenSurface(); break;
            case 30: OpenBio(); break;
            case 31: OpenSolvation(); break;
            case 44: OpenCg(); break;
            default:
                SetModule(8);
                Status = $"Editing {it.Name} with the builder tools (place, bond, delete, +H, fragments; ⌘Z undoes)";
                return;
        }
        Status = $"{label}: change the settings and build again — the new structure joins the project beside {it.Name}";
    }

    /// <summary>Delete (the project tree), asked twice: the structure leaves the project (a saved file stays on disk).</summary>
    public void DeleteProjectItem(ProjectItem it)
    {
        if (!it.DeleteArmed)
        {
            foreach (var p in ProjectItems) p.DeleteArmed = p == it;
            Status = $"Delete {it.Name}? Click Delete again to remove it from the project";
            Avalonia.Threading.DispatcherTimer.RunOnce(() => it.DeleteArmed = false, TimeSpan.FromSeconds(4));
            return;
        }
        it.DeleteArmed = false;
        if (!Idle || it.Doc.LongRunning) { Status = "A run is using this structure: delete it when the run finishes"; return; }
        if (it == _activeItem) { CloseDocument(); return; }
        ProjectItems.Remove(it);
        it.Doc.Dispose();
        RaiseProject();
        Status = $"Removed {it.Name} from the project";
    }

    private void RaiseProject()
    {
        Raise(nameof(ActiveItem));
        Raise(nameof(HasSeveralStructures));
        Raise(nameof(ProjectCountText));
        Raise(nameof(HasLastSession));
    }

    /// <summary>Keeps the active structure's view and workflow state in its item.</summary>
    private void StashActive()
    {
        var it = _activeItem;
        if (it == null) return;
        it.Camera = Camera;
        it.Frame = _frame;
        it.Done = new HashSet<string>(_pipeDone);
        it.Build = _pipeBuild;
        it.GrownUnsaved = GrownUnsaved;
        it.AutoFf = _pipeAutoFf;
        if (Title.Length > 0) it.Name = Title;
        UpdateItemInfo(it);
    }

    /// <summary>The project list's line for a structure: atoms, cell, force field and what was done to it.</summary>
    private void UpdateItemInfo(ProjectItem it)
    {
        try
        {
            var s = it.Doc.Summary();
            var inv = CultureInfo.InvariantCulture;
            it.Detail = s.Atoms.ToString("N0", inv) + " atoms" + (s.CellValid != 0 ? " · " + s.Density.ToString("0.000", inv) + " g/cm³" : "") +
                        (s.Frames > 1 ? " · " + s.Frames.ToString("N0", inv) + " frames" : "");
        }
        catch (ObjectDisposedException) { return; }
        if (it == _activeItem)
        {
            it.ForceField = Field.Assigned ? Field.ForceFieldName + (Field.Complete ? "" : " (incomplete)") : "";
            var done = new List<string>();
            if (_pipeDone.Contains("Relax")) done.Add("minimised");
            if (_pipeDone.Contains("Equilibrate")) done.Add("equilibrated");
            if (_pipeDone.Contains("Dynamics")) done.Add("MD");
            if (_pipeDone.Contains("React")) done.Add("reacted");
            it.History = string.Join(" · ", done);
            if (_pipeDone.Contains("Grow") && it.Origin.Length == 0) it.Origin = "Polymer cell";
            if (_pipeDone.Contains("Pack") && it.Origin.Length == 0) it.Origin = "Packing";
        }
    }

    /// <summary>Makes a structure of the project the active one: its view, force field and workflow state come back.</summary>
    public void Activate(ProjectItem it)
    {
        if (it == _activeItem || !ProjectItems.Contains(it)) return;
        if (it.Doc.IsDisposed)   // closed elsewhere: it is no longer part of the project
        {
            ProjectItems.Remove(it);
            RaiseProject();
            return;
        }
        if (Busy) { Status = "Wait for the run to finish (or cancel it) before switching structures"; Raise(nameof(ActiveItem)); return; }
        StashActive();
        foreach (var p in ProjectItems) p.Active = p == it;
        _activeItem = it;
        var doc = it.Doc;
        ClearFocus();
        Document = doc;
        RestraintsFollow(doc);
        Field.Reset();
        try { Field.LoadReport(doc); } catch { Field.Reset(); }
        Analyze.Load("");
        SyncHeld(); SyncFixed();
        _pipeDone.Clear();
        foreach (var d in it.Done) _pipeDone.Add(d);
        _pipeBuild = it.Build;
        _pipeAutoFf = it.AutoFf;
        GrownUnsaved = it.GrownUnsaved;
        Title = it.Name;
        FieldInfoText = "";
        var s = doc.Summary();
        Frames = (int)Math.Max(1, s.Frames);
        Raise(nameof(FrameMax));
        _frame = Math.Clamp(it.Frame, 0, Frames - 1);
        Raise(nameof(Frame));
        Raise(nameof(FrameLabel));
        Camera = it.Camera;
        _selection.Clear();
        Raise(nameof(CanAddRestraint));
        RefreshSelection();
        RefreshSummary();
        RefreshRdf();
        RefreshLegend();
        RefreshMolecules();
        IsPlaying = false;
        Notes.Clear();
        foreach (var n in doc.Notes()) Notes.Add(n);
        if (IsProvenance) LoadProvenance();
        LoadFileChecks();
        RefreshSteps();
        UpdateItemInfo(it);
        RaiseProject();
        if (IsExportCenter) RefreshEngines();
        Status = $"Working on {it.Name}";
        RenderRequested?.Invoke();
    }

    /// <summary>Closes every structure of the project and goes back to Start.</summary>
    public void CloseAllStructures()
    {
        if (Busy) { Status = "Wait for the run to finish (or cancel it) before closing"; return; }
        while (_doc != null) CloseDocument();
    }

    // Clear (the pipeline strip): every structure closed and a clean Start — asked twice (saved files stay on disk)
    private bool _clearArmed;
    public bool ClearArmed { get => _clearArmed; private set { if (Set(ref _clearArmed, value)) Raise(nameof(ClearText)); } }
    public string ClearText => _clearArmed ? "Click again to clear" : "Clear";
    public void ClearAll()
    {
        if (Busy) { Status = "Wait for the run to finish (or cancel it) before clearing"; return; }
        if (!_clearArmed)
        {
            ClearArmed = true;
            Status = $"Clear closes all {ProjectItems.Count} structure{(ProjectItems.Count == 1 ? "" : "s")} and starts clean (saved files stay): click Clear again";
            Avalonia.Threading.DispatcherTimer.RunOnce(() => ClearArmed = false, TimeSpan.FromSeconds(4));
            return;
        }
        ClearArmed = false;
        var n = ProjectItems.Count;
        CloseAllStructures();
        SetModule(8);
        Status = $"Cleared {n} structure{(n == 1 ? "" : "s")} · a clean start";
    }

    /// <summary>A copy of the active structure (its force field kept) as a new structure of the project.</summary>
    public void DuplicateStructure(string? suffix = null)
    {
        if (_doc == null || _activeItem == null) return;
        if (Busy) { Status = "Wait for the run to finish before copying the structure"; return; }
        var name = Title.Replace(" (unsaved)", "") + (suffix ?? " copy");
        var src = _activeItem;
        var copy = _doc.Copy(name);
        Show(copy, name + " (unsaved)");
        // the copy keeps what was done to the original and its force field
        foreach (var d in src.Done) _pipeDone.Add(d);
        _pipeBuild = src.Build;
        if (_activeItem != null) _activeItem.Origin = src.Origin;
        try { Field.LoadReport(copy); } catch { }
        RefreshSteps();
        if (_activeItem != null) UpdateItemInfo(_activeItem);
    }

    // Minimise, Equilibrate and Dynamics can keep the structure they start from (as Materials Studio writes a job's result
    // to a new document): the run works on a copy that joins the project
    private bool _resultAsNew;
    public bool ResultAsNew { get => _resultAsNew; set => Set(ref _resultAsNew, value); }

    /// <summary>Before a run: when asked, the run's own copy of the structure becomes the active one.</summary>
    private void PrepareRunTarget(string what)
    {
        if (!_resultAsNew || _doc == null) return;
        DuplicateStructure("");   // the run adds its own " · relaxed", " · MD" … to the copy's name
    }
}
