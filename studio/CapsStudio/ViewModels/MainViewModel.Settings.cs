using CapsStudio.Interop;

namespace CapsStudio.ViewModels;

/// <summary>A choice in the Settings page's colour palette list, with its swatches.</summary>
public sealed record PaletteChoice(int Index, string Name, string[] Swatches);

/// <summary>Settings (design/boards/Settings): appearance, the 3D view, force fields, compute, files. Every change
/// applies at once and is saved to ~/.caps/settings.json.</summary>
public sealed partial class MainViewModel
{
    public bool IsSettings => _module == 10;
    private AppSettings _settings = new();
    public AppSettings Settings => _settings;

    /// <summary>The window scales the interface when this changes.</summary>
    public event Action<double>? ScaleChanged;

    public static readonly string[] SettingsTabs = ["Appearance", "3D view", "Force fields", "Compute & remote", "Files"];
    public static readonly string[] SettingsIcons = ["eye", "cube", "tag", "server", "folder"];
    private int _settingsTab;
    public int SettingsTab { get => _settingsTab; set => Set(ref _settingsTab, value); }

    public PaletteChoice[] PaletteChoices { get; } =
    [
        new(0, "Element (CAPS default)", ["#8B969E", "#E9ECEF", "#E35049", "#2271DB", "#D6AC5C", "#7DC884"]),
        new(1, "Okabe–Ito (colour-blind safe)", ["#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00"]),
        new(2, "Monochrome (lightness only)", ["#F0F0F0", "#B4B4B4", "#7C7C7C", "#4A4A4A", "#D8D8D8", "#9C9C9C"]),
    ];

    /// <summary>Reads the settings file and applies it (the window calls this once at start).</summary>
    public void LoadSettings()
    {
        _settings = AppSettings.Load();
        ApplyAll();
    }

    private void ApplyAll()
    {
        Tokens.UseTheme(_settings.Theme);
        try { Native.SetPalette(_settings.Palette); Native.SetThreads(_settings.Threads); ApplyElectrostatics(); } catch { /* an older core: defaults */ }
        _viewBackground = _settings.Background; Raise(nameof(ViewBackground)); Raise(nameof(ViewIsLight));
        _outlines = _settings.Outlines; Raise(nameof(Outlines));
        _depthCue = _settings.DepthCue; Raise(nameof(DepthCue));
        _style = _settings.Style; Raise(nameof(StyleIndex)); Raise(nameof(StyleText));
        var ff = Field.Library.ToList().FindIndex(x => x.Id == _settings.ForceField);
        if (ff >= 0) Field.FfIndex = ff;
        var clean = CleanChoices.ToList().FindIndex(c => c.File != null && Path.GetFileNameWithoutExtension(c.File) == _settings.ForceField);
        if (clean >= 0) _molClean = clean;
        ScaleChanged?.Invoke(_settings.Scale);
        LoadHosts();
        Raise(nameof(JobTemplate));
        foreach (var n in new[] { nameof(SetTheme), nameof(SetScale), nameof(ScaleText), nameof(SetPalette), nameof(SetThreads), nameof(ThreadsText),
                                  nameof(SetBackground), nameof(SetOutlines), nameof(SetDepthCue), nameof(SetStyle), nameof(SetForceField),
                                  nameof(SetElectrostatics), nameof(PmeOn), nameof(SetEwaldExponent), nameof(SetPmeSpacing), nameof(SetPmeOrder) })
            Raise(n);
        RenderRequested?.Invoke();
        MolViewChanged?.Invoke();
    }

    private void Changed(string what)
    {
        _settings.Save();
        Status = $"{what} · saved to {AppSettings.DisplayPath}";
    }

    public string SetTheme
    {
        get => _settings.Theme;
        set
        {
            if (_settings.Theme == value) return;
            // the view follows Paper (white) and Graphite (dark) when it still has the old theme's background
            var wasLight = _settings.Theme == "light";
            _settings.Theme = value;
            Tokens.UseTheme(value);
            if (value is "light" or "dark" && _settings.Background == (wasLight ? 1 : 0)) SetBackground = value == "light" ? 1 : 0;
            Raise(); RenderRequested?.Invoke(); MolViewChanged?.Invoke(); Changed("Theme");
        }
    }
    public bool ThemeDark { get => _settings.Theme == "dark"; set { if (value) { SetTheme = "dark"; RaiseTheme(); } } }
    public bool ThemeLight { get => _settings.Theme == "light"; set { if (value) { SetTheme = "light"; RaiseTheme(); } } }
    public bool ThemeSystem { get => _settings.Theme == "system"; set { if (value) { SetTheme = "system"; RaiseTheme(); } } }
    private void RaiseTheme() { Raise(nameof(ThemeDark)); Raise(nameof(ThemeLight)); Raise(nameof(ThemeSystem)); }

    public double SetScale
    {
        get => _settings.Scale * 100;
        set
        {
            var v = Math.Clamp(Math.Round(value / 10) * 10, 80, 150) / 100;
            if (Math.Abs(_settings.Scale - v) < 1e-9) return;
            _settings.Scale = v;
            Raise(); Raise(nameof(ScaleText));
            ScaleChanged?.Invoke(v);
            Changed("Interface scale");
        }
    }
    public string ScaleText => $"{Math.Round(_settings.Scale * 100)} %";
    public void StepScale(int dir) => SetScale = _settings.Scale * 100 + 10 * dir;

    public int SetPalette
    {
        get => _settings.Palette;
        set
        {
            if (_settings.Palette == value) return;
            _settings.Palette = value;
            try { Native.SetPalette(value); } catch { }
            Raise();
            RefreshLegend();
            RenderRequested?.Invoke();
            MolViewChanged?.Invoke();
            Changed("Colour palette");
        }
    }

    public int SetThreads
    {
        get => _settings.Threads;
        set
        {
            var v = Math.Clamp(value, 0, 64);
            if (_settings.Threads == v) return;
            _settings.Threads = v;
            try { Native.SetThreads(v); } catch { }
            Raise(); Raise(nameof(ThreadsText));
            Changed("Threads");
        }
    }
    public decimal SetThreadsD { get => SetThreads; set => SetThreads = (int)value; }
    public string ThreadsText => _settings.Threads == 0 ? $"automatic: {Math.Min(16, Environment.ProcessorCount)} of {Environment.ProcessorCount}" : $"{_settings.Threads} of {Environment.ProcessorCount}";

    public int SetBackground { get => _settings.Background; set { if (_settings.Background == value) return; _settings.Background = value; ViewBackground = value; Raise(); Changed("View background"); } }
    public bool SetOutlines { get => _settings.Outlines; set { if (_settings.Outlines == value) return; _settings.Outlines = value; Outlines = value; Raise(); Changed("Outlines"); } }
    public bool SetDepthCue { get => _settings.DepthCue; set { if (_settings.DepthCue == value) return; _settings.DepthCue = value; DepthCue = value; Raise(); Changed("Depth cue"); } }
    public int SetStyle { get => _settings.Style; set { if (_settings.Style == value) return; _settings.Style = value; StyleIndex = value; Raise(); Changed("Display style"); } }

    /// <summary>The default force field: the Field page's first choice and the molecule builder's clean-up.</summary>
    public int SetForceField
    {
        get => Math.Max(0, Field.Library.ToList().FindIndex(x => x.Id == _settings.ForceField));
        set
        {
            if (value < 0 || value >= Field.Library.Count) return;
            var id = Field.Library[value].Id;
            if (_settings.ForceField == id) return;
            _settings.ForceField = id;
            Field.FfIndex = value;
            var clean = CleanChoices.ToList().FindIndex(c => c.File != null && Path.GetFileNameWithoutExtension(c.File) == id);
            if (clean >= 0) MolClean = clean;
            Raise();
            Changed("Default force field");
        }
    }

    // ---- electrostatics (Force fields tab): damped shifted force or particle-mesh Ewald
    private void ApplyElectrostatics() => Native.SetElectrostatics(_settings.Electrostatics, _settings.EwaldRtol, _settings.PmeSpacing, _settings.PmeOrder);
    public int SetElectrostatics
    {
        get => _settings.Electrostatics;
        set { if (_settings.Electrostatics == value) return; _settings.Electrostatics = value; ApplyElectrostatics(); Raise(); Raise(nameof(PmeOn)); Changed("Electrostatics"); }
    }
    public bool PmeOn => _settings.Electrostatics == 1;
    public int SetEwaldExponent
    {
        get => (int)Math.Round(Math.Log10(_settings.EwaldRtol));
        set { var v = Math.Pow(10, Math.Clamp(value, -10, -3)); if (Math.Abs(v - _settings.EwaldRtol) < 1e-15) return; _settings.EwaldRtol = v; ApplyElectrostatics(); Raise(); Changed("Ewald tolerance"); }
    }
    public decimal SetPmeSpacing
    {
        get => (decimal)_settings.PmeSpacing;
        set { var v = Math.Clamp((double)value, 0.3, 3.0); if (Math.Abs(v - _settings.PmeSpacing) < 1e-9) return; _settings.PmeSpacing = v; ApplyElectrostatics(); Raise(); Changed("PME grid spacing"); }
    }
    public decimal SetPmeOrder
    {
        get => _settings.PmeOrder;
        set { var v = Math.Clamp((int)value, 3, 10); if (v == _settings.PmeOrder) return; _settings.PmeOrder = v; ApplyElectrostatics(); Raise(); Changed("PME order"); }
    }
    public string ElectrostaticsText => _settings.Electrostatics == 1
        ? string.Format(System.Globalization.CultureInfo.InvariantCulture, "PME · β from erfc(β rc) = 1e{0} · grid ≤ {1:F2} Å · order {2}", SetEwaldExponent, _settings.PmeSpacing, _settings.PmeOrder)
        : "DSF · α 0.2 Å⁻¹";

    public string SettingsPath => AppSettings.DisplayPath;
    public string RecentPath => "~/.caps/recent.json";
    public string LibraryPath => Paths.ForceFields is { } p ? RecentFiles.Tilde(p) : "not found";
    public string SamplesPath => Paths.Samples is { } p ? RecentFiles.Tilde(p) : "not found";

    public void ClearRecent()
    {
        RecentFiles.Clear();
        LoadRecent();
        Status = "Recent files cleared";
    }

    public void ResetSettings()
    {
        _settings = new AppSettings();
        _settings.Save();
        ApplyAll();
        RaiseTheme();
        Status = $"Settings reset to the defaults · saved to {AppSettings.DisplayPath}";
    }

    public void ExportSettings(string path) { _settings.Save(path); Status = $"Settings exported to {path}"; }

    public void ImportSettings(string path)
    {
        _settings = AppSettings.Load(path);
        _settings.Save();
        ApplyAll();
        RaiseTheme();
        Status = $"Settings imported from {Path.GetFileName(path)}";
    }
}
