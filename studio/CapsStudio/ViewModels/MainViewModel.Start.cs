using System.Collections.ObjectModel;
using System.Globalization;
using CapsStudio.Interop;

namespace CapsStudio.ViewModels;

/// <summary>A builder on the Start board: where it leads (a module), or not yet available.</summary>
public sealed record StartBuilder(string Title, string Text, string Icon, int Module, string Note)
{
    public bool Available => Module >= 0;
}

/// <summary>A guided start: opens the sample (when needed) and goes to a module.</summary>
public sealed record StartGuide(string Title, string Where, int Module, bool Sample, string Action = "");

/// <summary>Start (design/boards/Start): the Studio with nothing open. A quick-start box that recognises a file
/// path, a SMILES string or a module; the builders; recent work; guided starts; this machine.</summary>
public sealed partial class MainViewModel
{
    public ObservableCollection<RecentItem> Recent { get; } = new();
    public bool HasRecent => Recent.Count > 0;

    public StartBuilder[] Builders { get; } =
    [
        new("Molecule", "From SMILES or a 2D sketch", "hex", 9, "Draw or type a molecule; 3D with conformers, cleaned with a force field"),
        new("Polymer", "Repeat units, sequence, tacticity", "grow", 13, "Homopolymers and copolymers from a SMILES library of repeat units"),
        new("Crystal", "Space group, lattice, CIF import", "cube", 29, "A space group, a lattice and an asymmetric unit expanded into a cell; symmetry found from a CIF"),
        new("Amorphous cell", "Grow and pack a periodic cell", "pack", 0, "Grow chains into a periodic cell at a target density"),
        new("Surface or interface", "Cleave, stack, add vacuum", "layers", 14, "Cleave a crystal (CIF) along (hkl) and grow a polymer film on it: fibre–rubber interfaces"),
        new("Solvated system", "Box, solvent model, ions", "flask", 31, "Water (TIP3P, SPC/E, TIP4P/2005) or a swelling solvent and salt around the open structure (CAPS Pack)"),
    ];

    /// <summary>Fibre–rubber composites, the builders and protocols in the order they are used.</summary>
    public StartGuide[] RubberGuides { get; } =
    [
        new("Rubber film on a glass-fibre surface", "Surface", 14, false),
        new("Silica, CNT or graphene filler in rubber", "Nano", 15, false),
        new("Sulfur-cure a natural-rubber cell", "React", 6, false, "sulfur"),
        new("Adhesion, orientation and pull-out", "Analyze", 1, false, "interface"),
    ];

    /// <summary>Runs a guided start: opens the module and sets it up for the task.</summary>
    public void StartGuideAction(StartGuide g)
    {
        if (g.Module == 14) { OpenSurface(); return; }
        if (g.Module == 15) { OpenNano(); return; }
        SetModule(g.Module);
        if (g.Action == "sulfur") RxSet = 2;
        if (g.Action == "interface")
            foreach (var c in Analyze.Groups.SelectMany(x => x.Chips)) c.IsOn = c.Id is "zprofile" or "adhesion" or "orientation" or "pull_shear";
    }

    public StartGuide[] Guides { get; } =
    [
        new("Grow an amorphous PS cell", "Grow", 0, false),
        new("Type the sample with GAFF2", "Field", 7, true),
        new("Relax the sample cell", "Relax", 2, true),
        new("Equilibrate with the 21-step protocol", "Equilibrate", 4, true),
        new("Density, g(r), Tg and moduli", "Analyze", 1, true),
    ];

    public string MachineThreads { get; } = Environment.ProcessorCount.ToString(CultureInfo.InvariantCulture);
    public string MachineRender { get; } = "CPU · software renderer";
    public string MachineOs { get; } = System.Runtime.InteropServices.RuntimeInformation.OSDescription;
    public string VersionText { get; } = $"CAPS {typeof(MainViewModel).Assembly.GetName().Version?.ToString(3)} · ABI {Native.AbiVersion()}";

    public void LoadRecent()
    {
        Recent.Clear();
        foreach (var r in RecentFiles.Load().Take(6)) Recent.Add(r);
        Raise(nameof(HasRecent));
        Raise(nameof(RecentFew));
        Raise(nameof(RecentCount));
    }

    /// <summary>Records the open document in the recent list with a thumbnail rendered off the UI thread.</summary>
    private void Remember(string path, string? topology)
    {
        var doc = _doc;
        if (doc == null) return;
        var s = doc.Summary();
        var detail = string.Format(CultureInfo.InvariantCulture, "{0:N0} atoms{1}", s.Atoms, s.Frames > 1 ? $" · {s.Frames:N0} frames" : "");
        var opt = new CapsRenderOpts { Width = 720, Height = 340, Supersample = 2, Background = 2, Style = 0, ColourBy = 0, Outlines = 1, DepthCue = 1, ShowCell = 1,
                                       Highlight0 = -1, Highlight1 = -1, Highlight2 = -1, Highlight3 = -1 };
        var cam = new CapsCamera { Yaw = 0.55, Pitch = 0.40, Zoom = 1.0 };
        Task.Run(() =>
        {
            RecentFiles.Touch(path, topology, detail, thumb => { try { doc.ExportPng(cam, opt, thumb); } catch { /* closed meanwhile */ } });
            Avalonia.Threading.Dispatcher.UIThread.Post(LoadRecent);
        });
    }

    /// <summary>Closes the document and goes back to Start.</summary>
    public void CloseDocument()
    {
        if (_doc == null) return;
        if (Busy) { Status = "Wait for the run to finish (or cancel it) before closing"; return; }
        var d = _doc;
        ClearFocus();
        Document = null;
        d.Dispose();
        Field.Reset();
        Analyze.Load("");
        Title = "";
        _selection.Clear();
        Frames = 1;
        Raise(nameof(FrameMax));
        IsPlaying = false;
        SetModule(8);
        Status = "Closed · Start";
        LoadRecent();
        RenderRequested?.Invoke();
    }

    // ---------------------------------------------------------------- quick start
    private string _quick = "";
    /// <summary>What the quick-start box recognised: 0 nothing, 1 a file, 2 SMILES, 3 a module, 4 not recognised.</summary>
    public int QuickKind { get; private set; }
    public string QuickBadge { get; private set; } = "";
    public string QuickAction { get; private set; } = "Open";
    public bool QuickReady => QuickKind is 1 or 2 or 3;
    public bool QuickHasBadge => QuickBadge.Length > 0;
    private int _quickModule = -1;

    public string QuickText
    {
        get => _quick;
        set
        {
            if (!Set(ref _quick, value)) return;
            var t = value.Trim().Trim('"', '\'');
            _quickModule = -1;
            if (t.Length == 0) (QuickKind, QuickBadge, QuickAction) = (0, "", "Open");
            else if (LooksLikePath(t) && File.Exists(ExpandHome(t))) (QuickKind, QuickBadge, QuickAction) = (1, "File found", "Open");
            else if (FindModule(t) is { } m) { _quickModule = m.Module; (QuickKind, QuickBadge, QuickAction) = (3, m.Name, "Go"); }
            else if (IsSmiles(t)) (QuickKind, QuickBadge, QuickAction) = (2, "SMILES detected", "Build 3D");
            else if (LooksLikePath(t)) (QuickKind, QuickBadge, QuickAction) = (4, "No such file", "Open");
            else (QuickKind, QuickBadge, QuickAction) = (4, "Not recognised", "Open");
            Raise(nameof(QuickKind)); Raise(nameof(QuickBadge)); Raise(nameof(QuickAction)); Raise(nameof(QuickReady)); Raise(nameof(QuickHasBadge));
            Raise(nameof(QuickOk)); Raise(nameof(QuickBad));
        }
    }
    public bool QuickOk => QuickKind is 1 or 2 or 3;
    public bool QuickBad => QuickKind == 4;

    /// <summary>Acts on the quick-start box; returns the file to open (the window opens it), or null.</summary>
    public string? QuickGo()
    {
        switch (QuickKind)
        {
            case 1: return ExpandHome(_quick.Trim().Trim('"', '\''));
            case 3: SetModule(_quickModule); return null;
            case 2: OpenBuilder(_quick.Trim()); return null;
            default: return null;
        }
    }

    private static string ExpandHome(string p) =>
        p.StartsWith('~') ? Environment.GetFolderPath(Environment.SpecialFolder.UserProfile) + p[1..] : p;

    private static bool LooksLikePath(string t) => t.Contains('/') && (t.StartsWith('/') || t.StartsWith('~') || t.StartsWith('.')) || t.Contains('\\') ||
        System.IO.Path.GetExtension(t) is ".data" or ".lmp" or ".lammpstrj" or ".dump" or ".gro" or ".pdb" or ".xyz" or ".extxyz" or ".mol2";

    private static readonly (string Word, int Module, string Name)[] ModuleWords =
    [
        ("studio", 8, "Go to Studio"), ("grow", 0, "Go to Grow"), ("amorphous", 0, "Go to Grow"), ("polymer", 0, "Go to Grow"),
        ("pack", 5, "Go to Pack"), ("packmol", 5, "Go to Pack"), ("solvate", 5, "Go to Pack"),
        ("relax", 2, "Go to Relax"), ("minimise", 2, "Go to Relax"), ("minimize", 2, "Go to Relax"),
        ("dynamics", 3, "Go to Dynamics"), ("md", 3, "Go to Dynamics"), ("lammps", 3, "Go to Dynamics"),
        ("equilibrate", 4, "Go to Equilibrate"), ("protocol", 4, "Go to Equilibrate"),
        ("react", 6, "Go to React"), ("crosslink", 6, "Go to React"),
        ("analyze", 1, "Go to Analyze"), ("analyse", 1, "Go to Analyze"), ("tg", 1, "Go to Analyze"), ("mechanics", 1, "Go to Analyze"), ("rdf", 1, "Go to Analyze"),
        ("field", 7, "Go to Field"), ("force field", 7, "Go to Field"), ("gaff", 7, "Go to Field"), ("typing", 7, "Go to Field"),
    ];

    private static (int Module, string Name)? FindModule(string t)
    {
        var w = t.ToLowerInvariant();
        foreach (var m in ModuleWords)
            if (w == m.Word || (w.Length >= 3 && m.Word.StartsWith(w, StringComparison.Ordinal))) return (m.Module, m.Name);
        return null;
    }

    /// <summary>A SMILES string by its grammar (organic subset, bracket atoms, bonds, branches, ring closures).</summary>
    public static bool IsSmiles(string s)
    {
        if (s.Length == 0 || s.Contains(' ')) return false;
        int depth = 0, atoms = 0;
        for (var i = 0; i < s.Length;)
        {
            var c = s[i];
            if (c == '[')
            {
                var j = s.IndexOf(']', i);
                if (j < 0) return false;
                var inner = s[(i + 1)..j];
                if (inner.Length == 0 || !inner.Any(char.IsLetter)) return false;
                atoms++; i = j + 1; continue;
            }
            if (i + 1 < s.Length && (s.Substring(i, 2) is "Cl" or "Br")) { atoms++; i += 2; continue; }
            if ("BCNOPSFIbcnops*".IndexOf(c) >= 0) { atoms++; i++; continue; }
            if (c == '(') { depth++; i++; continue; }
            if (c == ')') { if (--depth < 0) return false; i++; continue; }
            if (char.IsDigit(c) || "-=#$:/\\.%+@".IndexOf(c) >= 0) { i++; continue; }
            return false;
        }
        return depth == 0 && atoms > 0;
    }
}
