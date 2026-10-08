using System.Text.Json;
using System.Text.Json.Serialization;

namespace CapsStudio.ViewModels;

/// <summary>User settings (design/boards/Settings), saved to ~/.caps/settings.json. A missing or unreadable file gives
/// the defaults; unknown fields are ignored, so older and newer Studios share the file.</summary>
public sealed class AppSettings
{
    public string Theme { get; set; } = "dark";          // dark | light | system
    public double Scale { get; set; } = 1.0;              // interface scale, 0.8 … 1.5
    public int Palette { get; set; }                      // 0 CAPS, 1 Okabe–Ito, 2 monochrome
    public int Threads { get; set; }                      // 0 = automatic
    public int Background { get; set; }                   // older Studios' view background (0 dark, 1 white): kept for them, not read
    public string ViewBackdrop { get; set; } = "theme";   // view background: theme (follows the Studio's theme) | dark | white
    public bool Outlines { get; set; } = true;
    public bool DepthCue { get; set; } = true;
    public bool AutoClean { get; set; }                     // a UFF clean-up of the edited atoms after every builder edit (A)
    public bool HAutopilot { get; set; }                    // after every sketch edit: missing hydrogens added, surplus ones removed
    public bool NotifyRuns { get; set; } = true;            // a system notification when a run of a minute or more ends
    public bool GpuView { get; set; } = true;             // the 3D view on the GPU (OpenGL) where available
    // the main window as it was closed (device pixels); restored only onto a screen that is connected now
    public int WindowX { get; set; }
    public int WindowY { get; set; }
    public double WindowW { get; set; }
    public double WindowH { get; set; }
    public bool WindowMax { get; set; }
    public bool WindowSaved { get; set; }
    public int Style { get; set; }                        // ball & stick …
    public bool AutoStyle { get; set; } = true;           // display style by model size (design/boards/DisplayStyles)
    public long AutoNoH { get; set; } = 20_000;           // above this many atoms: No H
    public long AutoBackbone { get; set; } = 200_000;     // above: Backbone (with a 30 Å all-atom lens)
    public long AutoLodAtoms { get; set; } = 2_000_000;   // above: Backbone + level of detail
    /// <summary>Saved pipeline expressions (design/boards/ExpressionSelect), counted on every frame shown.</summary>
    public List<string> PipelineExpressions { get; set; } =
        ["Type == 2", "Type == 2 && Position.Z > 13", "Charge < -0.05", "MoleculeIdentifier == 3", "Element == \"H\" && Position.X < 10", "Monomer <= 2 || Monomer >= 7"];
    private string _forceField = "gaff-amber25";
    /// <summary>A library id; ids saved before the force fields got CAPS's own names drop their old suffix.</summary>
    public string ForceField { get => _forceField; set => _forceField = value.EndsWith("-dlfield", StringComparison.Ordinal) ? value[..^8] : value; }
    public int Electrostatics { get; set; }               // 0 damped shifted force, 1 particle-mesh Ewald
    public bool LjTail { get; set; } = true;              // LJ long-range tail correction in Dynamics and Equilibrate
    public double EwaldRtol { get; set; } = 1e-5;
    public double PmeSpacing { get; set; } = 1.0;
    public int PmeOrder { get; set; } = 5;
    public int FigureBackground { get; set; } = 1;        // export figure: 0 dark, 1 white, 2 transparent
    public bool FigurePerProject { get; set; } = true;    // export figure: the background kept with the open project
    public int FigureBits { get; set; } = 8;              // export figure: PNG bit depth (8 or 16)
    public int FigurePreset { get; set; }                 // journal single column …
    public int ReaderVerbosity { get; set; } = 1;         // keyboard walk: 0 brief, 1 full
    public bool AnnounceFrames { get; set; }
    public bool SelectionSound { get; set; }
    /// <summary>Visualize: steps slower than SlowStepMs (and Python steps) wait while the frame slider is dragged.</summary>
    public bool PauseSlowSteps { get; set; } = true;
    public bool EngineNotice { get; set; } = true;   // say when a project was last saved with another CAPS version   // a short system sound when atoms are picked or selected
    public bool ShowAnnouncement { get; set; } = true;    // the announcement bar in the view, for sighted keyboard users
    public List<RemoteHost> Hosts { get; set; } = new();   // Compute & remote: SSH hosts (no credentials: the SSH agent holds them)
    public string JobTemplate { get; set; } = RemoteHost.DefaultTemplate;
    public List<MyFragment> MyFragments { get; set; } = new();   // the fragment library's "My fragments"
    public List<MyFragment> PackMolecules { get; set; } = new();
    public string ProjectParent { get; set; } = "";
    public LookSettings Look { get; set; } = new();   // Look: Keep as my look
    public List<BondRuleSet> BondRuleSets { get; set; } = new();   // Bond rules: rule sets kept by name   // where the last new project was made (New project starts there)   // Pack › Add molecule: "Your molecules" (whole molecules by SMILES)
    public List<SavedQueryData> SavedQueries { get; set; } = SavedQueryData.Defaults();   // Select by query (SmartSelect)
    /// <summary>The first-run tour was finished or skipped.</summary>
    public bool TourDone { get; set; }
    /// <summary>Keys the user gave palette commands: command id → gesture ("Meta+Shift+G").</summary>
    public Dictionary<string, string> Shortcuts { get; set; } = new();
    /// <summary>Pane sizes the user dragged or hid (Panes): grid key → each track's size ("300|*|child:230").</summary>
    public Dictionary<string, string> PaneSizes { get; set; } = new();
    /// <summary>Tool shelves (design/boards/Shelves): the workspace in use, each workspace's layout, the shelves the user
    /// made, and whether shelves are locked in place.</summary>
    public string Workspace { get; set; } = "Sketch";
    public Dictionary<string, List<ShelfState>> ShelfLayouts { get; set; } = new();
    public List<CustomShelf> CustomShelves { get; set; } = new();
    public bool ShelvesLocked { get; set; }
    public bool ShelfHintShown { get; set; }                // the one-time hint about shelf tabs was shown
    /// <summary>Names the user gave projects: folder ("" for the session) → name shown in the project tree.</summary>
    public Dictionary<string, string> ProjectNames { get; set; } = new();
    public string ReduceMotion { get; set; } = "system";   // system (follow the OS) | on | off
    public bool HighContrast { get; set; }                  // Accessibility: strong outlines in the view, firmer borders and secondary text
    public bool FocusRingsAlways { get; set; }              // Accessibility: a focus ring on every focused control, not only after Tab
    public bool AnnounceProgress { get; set; }              // Accessibility: job progress read out at each quarter
    public string PythonExe { get; set; } = "";             // Python & scripting: the interpreter ("" CAPS_PYTHON, else python3)
    public int UnitSystem { get; set; }                    // Settings › Units: 0 CAPS, 1 SI-derived, 2 LAMMPS metal
    public string UnitEnergy { get; set; } = "kcal/mol";
    public string UnitLength { get; set; } = "Å";
    public string UnitPressure { get; set; } = "atm";
    public string UnitTime { get; set; } = "ps";

    public static string Folder => Path.Combine(Environment.GetFolderPath(Environment.SpecialFolder.UserProfile), ".caps");
    /// <summary>Tests and screenshots point elsewhere so they never change the user's file.</summary>
    public static string? Override { get; set; }
    public static string FilePath => Override ?? Path.Combine(Folder, "settings.json");
    public static string DisplayPath => Override ?? "~/.caps/settings.json";

    private static readonly JsonSerializerOptions Json = new() { WriteIndented = true, PropertyNamingPolicy = JsonNamingPolicy.SnakeCaseLower };

    public static AppSettings Load(string? path = null)
    {
        try
        {
            var p = path ?? FilePath;
            if (File.Exists(p)) return (JsonSerializer.Deserialize<AppSettings>(File.ReadAllText(p), Json) ?? new AppSettings()).Clamped();
        }
        catch { /* defaults */ }
        return new AppSettings();
    }

    public void Save(string? path = null)
    {
        try
        {
            var p = path ?? FilePath;
            Directory.CreateDirectory(Path.GetDirectoryName(p)!);
            File.WriteAllText(p, JsonSerializer.Serialize(this, Json));
        }
        catch { /* a read-only home keeps the settings for this session only */ }
    }

    private AppSettings Clamped()
    {
        Scale = Math.Clamp(Scale, 0.8, 1.5);
        Palette = Math.Clamp(Palette, 0, 2);
        Threads = Math.Clamp(Threads, 0, 64);
        Background = Math.Clamp(Background, 0, 1);
        if (ViewBackdrop is not ("theme" or "dark" or "white")) ViewBackdrop = "theme";
        Style = Math.Clamp(Style, 0, 4);
        Electrostatics = Math.Clamp(Electrostatics, 0, 1);
        EwaldRtol = Math.Clamp(EwaldRtol, 1e-10, 1e-2);
        PmeSpacing = Math.Clamp(PmeSpacing, 0.3, 3.0);
        PmeOrder = Math.Clamp(PmeOrder, 3, 10);
        if (Theme is not ("dark" or "light" or "system")) Theme = "dark";
        return this;
    }
}

/// <summary>A remote machine jobs can run on (design/boards/RemoteCompute): reached with SSH through the user's agent.</summary>
public sealed class RemoteHost
{
    public string Name { get; set; } = "hpc-login";
    public string Hostname { get; set; } = "";
    public string User { get; set; } = "";
    public int Port { get; set; } = 22;
    public string Scheduler { get; set; } = "SLURM";      // SLURM | PBS | none
    public string Partition { get; set; } = "";
    public string WorkDir { get; set; } = "/scratch/$USER/caps";

    public const string DefaultTemplate = "#!/bin/bash\n#SBATCH --job-name=caps-{job}\n#SBATCH --partition={partition}\n#SBATCH --cpus-per-task=8\n" +
                                          "#SBATCH --time=24:00:00\n#SBATCH --output=caps-%j.log\ncd {workdir}/{job}\ncaps run {recipe}\n";
}

/// <summary>A fragment the user saved: a name and SMILES with * attachment points.</summary>
/// <summary>The view's sizes kept as the user's own look (design/boards/Look).</summary>
public sealed class LookSettings
{
    public double AtomScale { get; set; } = 0.28;
    public double Stick { get; set; } = 0.14;
    public double Space { get; set; } = 1.0;
    public double Line { get; set; } = 1.4;
    public bool Orders { get; set; }
    public bool Arcs { get; set; }
}

public sealed class MyFragment
{
    public string Name { get; set; } = "";
    public string Smiles { get; set; } = "";
}

/// <summary>A named selection query (design/boards/SmartSelect).</summary>
public sealed class SavedQueryData
{
    public string Name { get; set; } = "";
    public string Query { get; set; } = "";
    public static List<SavedQueryData> Defaults() =>
    [
        new() { Name = "Aromatic rings", Query = "smarts \"c1ccccc1\"" },
        new() { Name = "Backbone CH₂", Query = "smarts \"[CH2;!R]\"" },
        new() { Name = "Stereo centres", Query = "stereo *" },
        new() { Name = "Hydrogens", Query = "element H" },
        new() { Name = "Near ring 1", Query = "within 5 of ring 1" },
        new() { Name = "… heavy only", Query = "within 5 of ring 1 and not element H" },
    ];
}

/// <summary>Where a shelf is: dock top-left | top-right | top-2 | left | right | bottom | float | window | hidden; its
/// order in the dock; its place when floating (window coordinates, or screen for a popped-out one); folded to a puck.</summary>
public sealed class ShelfState
{
    public string Id { get; set; } = "";
    public string Dock { get; set; } = "top-left";
    public int Order { get; set; }
    public double X { get; set; }
    public double Y { get; set; }
    public bool Folded { get; set; }
}

/// <summary>A shelf the user made: a name, a glyph and the tools on it (ids of the built-in tools).</summary>
public sealed class CustomShelf
{
    public string Id { get; set; } = "";
    public string Name { get; set; } = "";
    public string Glyph { get; set; } = "pin";
    public List<string> Tools { get; set; } = new();
    /// <summary>Tool size in px: 28 small, 36 medium, 44 large (design/boards/ShelfEditor).</summary>
    public int Size { get; set; } = 36;
    /// <summary>Names in a tooltip on hover (true) or written under each tool (false).</summary>
    public bool LabelsOnHover { get; set; } = true;
}
