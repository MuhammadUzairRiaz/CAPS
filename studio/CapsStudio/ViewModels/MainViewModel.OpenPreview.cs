using System.Collections.ObjectModel;
using System.Globalization;
using System.Text.Json.Nodes;
using CapsStudio.Interop;

namespace CapsStudio.ViewModels;

public sealed record ColumnRow(string Column, string MapsTo, string Kind, bool Used);
public sealed record TypeRow(string Type, string Name, string Mass, string Element);
public sealed record FormatRow(string Name, string About, string Reads, string Writes);

/// <summary>Open file (design/boards/OpenLammps, OpenGromacs): the format found from the content, the first lines, how
/// dump columns map to properties, types and masses, where bonds come from and how many frames there are — then open.</summary>
public sealed partial class MainViewModel
{
    public bool IsOpenPreview => _module == 27;
    public ObservableCollection<PreviewLine> OpenHead { get; } = new();
    public ObservableCollection<ColumnRow> OpenColumns { get; } = new();
    public ObservableCollection<TypeRow> OpenTypes { get; } = new();
    public ObservableCollection<string> OpenNotes { get; } = new();
    public static readonly FormatRow[] FormatMatrix =
    [
        new("LAMMPS data", "atom_style full · molecular · charge · atomic; Masses, Bonds", "yes", "yes, with force-field sections"),
        new("LAMMPS dump", "text custom columns: id mol type element q x/xu/xs ix", "yes", "yes"),
        new("GROMACS .gro", "positions in nm, residues, box", "yes", "yes"),
        new("PDB", "ATOM/HETATM, CONECT, CRYST1", "yes", "yes"),
        new("XYZ · extended XYZ", "lattice in the comment line", "yes", "yes"),
        new("Tripos mol2", "bond orders, charges, atom types", "yes", "yes"),
        new("CIF", "crystals, symmetry expanded", "yes", "—"),
    ];
    private string _openPath = "", _openTopo = "", _openFormat = "", _openTitle = "", _openAtoms = "", _openFrames = "", _openBonds = "", _openUnits = "", _openMapped = "";
    private bool _openBusy, _openIsDump;
    private int _openGen;
    private int _returnModule = 8;

    public string OpenPath => _openPath;
    public string OpenTopology => _openTopo;
    public string OpenFileName => System.IO.Path.GetFileName(_openPath);
    public string OpenFormat { get => _openFormat; private set => Set(ref _openFormat, value); }
    public string OpenTitle { get => _openTitle; private set => Set(ref _openTitle, value); }
    public string OpenAtoms { get => _openAtoms; private set => Set(ref _openAtoms, value); }
    public string OpenFrames { get => _openFrames; private set => Set(ref _openFrames, value); }
    public string OpenBonds { get => _openBonds; private set => Set(ref _openBonds, value); }
    public string OpenUnits { get => _openUnits; private set => Set(ref _openUnits, value); }
    public string OpenMapped { get => _openMapped; private set => Set(ref _openMapped, value); }
    public bool OpenIsDump { get => _openIsDump; private set => Set(ref _openIsDump, value); }
    public bool OpenBusy { get => _openBusy; private set { if (Set(ref _openBusy, value)) Raise(nameof(OpenIdle)); } }
    public bool OpenIdle => !_openBusy;
    public string OpenButton => OpenManyFrames ? $"Open {OpenSelection().Kept(OpenFrameCount):N0} frames" : "Open";

    // Which frames to read (a long trajectory held at the size of what is kept; the others are passed over as it is read)
    private double _openFrameFrom, _openFrameTo, _openFrameEvery = 1;
    private FrameSelection? _pendingFrames;   // taken by the next Open
    private long OpenFrameCount => long.TryParse(_openFrames, out var n) ? n : 0;
    public bool OpenManyFrames => OpenFrameCount > 1;
    public double OpenFrameMax => Math.Max(0, OpenFrameCount - 1);
    public double OpenFrameFrom { get => _openFrameFrom; set { if (Set(ref _openFrameFrom, Math.Clamp(Math.Round(value), 0, OpenFrameMax))) FramesChanged(); } }
    public double OpenFrameTo { get => _openFrameTo; set { if (Set(ref _openFrameTo, Math.Clamp(Math.Round(value), 0, OpenFrameMax))) FramesChanged(); } }
    public double OpenFrameEvery { get => _openFrameEvery; set { if (Set(ref _openFrameEvery, Math.Max(1, Math.Round(value)))) FramesChanged(); } }
    public string OpenFramesKept
    {
        get
        {
            var sel = OpenSelection();
            var kept = sel.Kept(OpenFrameCount);
            return sel.All ? "all frames are read" : kept == 0 ? "no frames in this range"
                 : string.Format(CultureInfo.InvariantCulture, "{0:N0} of {1:N0} frames kept · the others are passed over as the file is read", kept, OpenFrameCount);
        }
    }
    private FrameSelection OpenSelection() =>
        new((long)_openFrameFrom, _openFrameTo >= OpenFrameMax ? -1 : (long)_openFrameTo, (long)Math.Max(1, _openFrameEvery));
    private void FramesChanged() { Raise(nameof(OpenFramesKept)); Raise(nameof(OpenButton)); }
    /// <summary>Sets the frames to read (the Open page's from / to / every; to past the end reads to the end).</summary>
    public void SetOpenFrames(long from, long to, long every)
    {
        OpenFrameFrom = from;
        OpenFrameTo = to < 0 ? OpenFrameMax : to;
        OpenFrameEvery = every;
    }

    /// <summary>Shows what the file holds (read off the UI thread); a dump picks up a data file beside it as its topology.</summary>
    public void PreviewOpen(string path, string? topology = null)
    {
        if (IsReaderFile(path)) { OpenReader(path); return; }
        if (path.EndsWith(".capstable", StringComparison.OrdinalIgnoreCase)) { Status = StudyLoad(path); return; }
        if (topology == null && NeedsImport(path)) { ShowImport(path); return; }
        if (_module != 27) _returnModule = _module;
        _openPath = path;
        _openTopo = topology ?? TopologyFor(path);
        FindGmxSet(path);
        Raise(nameof(OpenPath)); Raise(nameof(OpenTopology)); Raise(nameof(OpenFileName));
        SetModule(27);
        var gen = ++_openGen;
        OpenBusy = true;
        OpenTitle = "Reading " + OpenFileName;
        var (p, t) = (_openPath, _openTopo);
        Task.Run(() =>
        {
            string? json = null, error = null;
            try { json = CapsDocument.InspectFile(p, t.Length > 0 ? t : null); } catch (Exception e) { error = e.Message; }
            Avalonia.Threading.Dispatcher.UIThread.Post(() =>
            {
                if (gen != _openGen) return;
                OpenBusy = false;
                FillOpenPreview(json, error);
            });
        });
    }

    /// <summary>The preview read at once on this thread (the self-test; the page reads it in the background).</summary>
    internal void PreviewOpenNow(string path, string? topology)
    {
        _openPath = path;
        _openTopo = topology ?? "";
        string? json = null, error = null;
        try { json = CapsDocument.InspectFile(path, _openTopo.Length > 0 ? _openTopo : null); } catch (Exception e) { error = e.Message; }
        FillOpenPreview(json, error);
    }

    private void FillOpenPreview(string? json, string? error)
    {
        OpenHead.Clear(); OpenColumns.Clear(); OpenTypes.Clear(); OpenNotes.Clear();
        if (json == null) { OpenTitle = "Cannot read " + OpenFileName; OpenFormat = "unknown"; OpenNotes.Add(error ?? ""); Raise(nameof(OpenButton)); return; }
        var j = JsonNode.Parse(json)!;
        var inv = CultureInfo.InvariantCulture;
        OpenFormat = (string?)j["format_name"] ?? "";
        OpenTitle = "Open · " + OpenFormat;
        OpenIsDump = (string?)j["format"] == "lammps-dump";
        int k = 0;
        foreach (var l in (JsonArray)j["head"]!) OpenHead.Add(new PreviewLine((++k).ToString(inv), (string?)l ?? ""));
        int used = 0, all = 0;
        foreach (var c in (JsonArray)j["columns"]!)
        {
            var u = (bool?)c!["used"] ?? true;
            used += u ? 1 : 0; ++all;
            OpenColumns.Add(new ColumnRow((string?)c["name"] ?? "", (string?)c["maps_to"] ?? "", (string?)c["kind"] ?? "", u));
        }
        OpenMapped = all > 0 ? $"{used} of {all} mapped" : "";
        foreach (var t in (JsonArray)j["types"]!)
            OpenTypes.Add(new TypeRow(((double?)t!["type"] ?? 0) > 0 ? ((double?)t["type"] ?? 0).ToString("0", inv) : "—", (string?)t["label"] ?? "",
                                      ((double?)t["mass"] ?? 0).ToString("F3", inv), (string?)t["element"] ?? ""));
        OpenAtoms = ((double?)j["atoms"] ?? 0).ToString("N0", inv);
        OpenFrames = ((double?)j["frames"] ?? 0).ToString("0", inv);
        _openFrameFrom = 0; _openFrameTo = OpenFrameMax; _openFrameEvery = 1;
        foreach (var n in new[] { nameof(OpenManyFrames), nameof(OpenFrameMax), nameof(OpenFrameFrom), nameof(OpenFrameTo), nameof(OpenFrameEvery) }) Raise(n);
        FramesChanged();
        OpenBonds = (string?)j["bonds_from"] ?? "";
        OpenUnits = (string?)j["units"] ?? "";
        foreach (var n in (JsonArray)j["notes"]!) OpenNotes.Add((string?)n ?? "");
        Raise(nameof(OpenButton));
    }

    // ---- GROMACS set (design/boards/OpenGromacs): coordinates, topology with its includes, and trajectory found together
    // by name next to whichever one was dropped, opened as one document
    public ObservableCollection<GmxSetRow> GmxSet { get; } = new();
    private string _gmxGro = "", _gmxTop = "", _gmxTraj = "";
    public bool HasGmxSet => GmxSet.Count > 0;
    public string GmxSetTitle => $"Files in this set · {GmxSet.Count(r => r.Found && r.Role != "include")} of 3 found";
    public bool CanOpenGmxSet => _gmxGro.Length > 0 && (_gmxTop.Length > 0 || _gmxTraj.Length > 0) || _gmxTraj.Length > 0 && _gmxTop.Length > 0;

    internal void FindGmxSet(string path)
    {
        GmxSet.Clear();
        _gmxGro = _gmxTop = _gmxTraj = "";
        var ext = Path.GetExtension(path).ToLowerInvariant();
        if (ext is not (".gro" or ".top" or ".xtc" or ".trr")) { RaiseGmx(); return; }
        var dir = Path.GetDirectoryName(Path.GetFullPath(path)) ?? ".";
        string Find(params string[] exts)
        {
            foreach (var e in exts) { var same = Path.ChangeExtension(path, e); if (File.Exists(same)) return same; }
            foreach (var e in exts)
            {
                try { var m = Directory.GetFiles(dir, "*" + e); if (m.Length == 1) return m[0]; } catch { }
            }
            foreach (var n in new[] { "topol.top", "conf.gro", "confout.gro", "traj_comp.xtc", "traj.xtc", "traj.trr" })
                if (exts.Contains(Path.GetExtension(n)) && File.Exists(Path.Combine(dir, n))) return Path.Combine(dir, n);
            return "";
        }
        _gmxGro = ext == ".gro" ? path : Find(".gro");
        _gmxTop = ext == ".top" ? path : Find(".top");
        _gmxTraj = ext is ".xtc" or ".trr" ? path : Find(".xtc", ".trr");
        string Size(string f) { try { var b = new FileInfo(f).Length; return b > 1 << 20 ? $"{b / 1048576.0:0.#} MB" : $"{b / 1024.0:0} KB"; } catch { return ""; } }
        GmxSet.Add(new GmxSetRow("coordinates", _gmxGro.Length > 0 ? Path.GetFileName(_gmxGro) : "no .gro", _gmxGro.Length > 0, _gmxGro.Length > 0 ? "coordinates · box · " + Size(_gmxGro) : "the first trajectory frame stands in"));
        GmxSet.Add(new GmxSetRow("topology", _gmxTop.Length > 0 ? Path.GetFileName(_gmxTop) : "no .top", _gmxTop.Length > 0, _gmxTop.Length > 0 ? "bonds, types, charges" : "bonds are perceived from distances"));
        if (_gmxTop.Length > 0)
            try
            {
                // its #include files: found beside it (force-field includes from the GROMACS share are listed as not here)
                foreach (var l in File.ReadLines(_gmxTop).Take(4000))
                {
                    var t = l.Trim();
                    if (!t.StartsWith("#include")) continue;
                    var inc = t[8..].Trim().Trim('"', '<', '>');
                    var at = Path.Combine(Path.GetDirectoryName(_gmxTop) ?? ".", inc);
                    GmxSet.Add(new GmxSetRow("include", "  + " + inc, File.Exists(at), File.Exists(at) ? "included by the topology" : "not beside the topology (a GROMACS share file is not read)"));
                }
            }
            catch { }
        GmxSet.Add(new GmxSetRow("trajectory", _gmxTraj.Length > 0 ? Path.GetFileName(_gmxTraj) : "no .xtc / .trr", _gmxTraj.Length > 0, _gmxTraj.Length > 0 ? "frames · " + Size(_gmxTraj) : "one frame"));
        RaiseGmx();
    }
    private void RaiseGmx() { Raise(nameof(HasGmxSet)); Raise(nameof(GmxSetTitle)); Raise(nameof(CanOpenGmxSet)); }

    /// <summary>Open set: the trajectory with the topology (or the coordinates), else the coordinates with the topology.</summary>
    public void OpenGmxSet()
    {
        if (!CanOpenGmxSet) return;
        var main = _gmxTraj.Length > 0 ? _gmxTraj : _gmxGro;
        var topo = _gmxTop.Length > 0 ? _gmxTop : _gmxGro;
        if (main == topo) topo = "";
        SetModule(_returnModule is 27 ? 8 : _returnModule);
        OpenRequested?.Invoke(main + (topo.Length > 0 ? "\n" + topo : ""));
    }

    public void SetOpenTopology(string? topology)
    {
        PreviewOpen(_openPath, topology ?? "");
    }

    public void CancelOpen() => SetModule(_returnModule is 27 ? 8 : _returnModule);

    /// <summary>Opens the previewed file (with its topology) and goes back to where it was asked for.</summary>
    public void ConfirmOpen()
    {
        var back = _returnModule is 27 ? 8 : _returnModule;
        var sel = OpenSelection();
        if (OpenManyFrames && !sel.All)
        {
            if (sel.Kept(OpenFrameCount) == 0) { Status = "No frames in the range asked for"; return; }
            _pendingFrames = sel;
        }
        SetModule(back);
        OpenRequested?.Invoke(_openPath + (_openTopo.Length > 0 ? "\n" + _openTopo : ""));
    }
}

/// <summary>A file of a GROMACS set (design/boards/OpenGromacs): its role, name, whether it is there, what it carries.</summary>
public sealed record GmxSetRow(string Role, string Name, bool Found, string Detail)
{
    public Avalonia.Media.IBrush Brush => Tokens.Brush(Found ? "OkB" : "DimB");
}
