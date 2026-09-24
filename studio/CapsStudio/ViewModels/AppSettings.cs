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
    public int Background { get; set; }                   // view background: 0 dark, 1 white
    public bool Outlines { get; set; } = true;
    public bool DepthCue { get; set; } = true;
    public int Style { get; set; }                        // ball & stick …
    public string ForceField { get; set; } = "gaff-amber25-dlfield";
    public int Electrostatics { get; set; }               // 0 damped shifted force, 1 particle-mesh Ewald
    public double EwaldRtol { get; set; } = 1e-5;
    public double PmeSpacing { get; set; } = 1.0;
    public int PmeOrder { get; set; } = 5;

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
        Style = Math.Clamp(Style, 0, 4);
        Electrostatics = Math.Clamp(Electrostatics, 0, 1);
        EwaldRtol = Math.Clamp(EwaldRtol, 1e-10, 1e-2);
        PmeSpacing = Math.Clamp(PmeSpacing, 0.3, 3.0);
        PmeOrder = Math.Clamp(PmeOrder, 3, 10);
        if (Theme is not ("dark" or "light" or "system")) Theme = "dark";
        return this;
    }
}
