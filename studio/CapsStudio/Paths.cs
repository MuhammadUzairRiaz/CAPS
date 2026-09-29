namespace CapsStudio;

/// <summary>Where CAPS keeps its shipped resources (data/ with force fields, typing rules and reference values;
/// samples/). Looked up in this order: $CAPS_HOME; next to the executable (Windows installer, AppImage, tarball);
/// ../Resources (macOS .app bundle); ../share/caps (Linux /usr prefix: .deb, PKGBUILD, Flatpak /app); and the
/// development tree (studio/CapsStudio/bin/... → the repository).</summary>
public static class Paths
{
    private static IEnumerable<string> Roots()
    {
        var env = Environment.GetEnvironmentVariable("CAPS_HOME");
        if (!string.IsNullOrEmpty(env)) yield return env;
        var baseDir = AppContext.BaseDirectory.TrimEnd(Path.DirectorySeparatorChar, Path.AltDirectorySeparatorChar);
        yield return baseDir;
        var parent = Path.GetDirectoryName(baseDir);
        if (parent != null)
        {
            yield return Path.Combine(parent, "Resources");
            yield return Path.Combine(parent, "share", "caps");
        }
        for (var dir = parent; dir != null; dir = Path.GetDirectoryName(dir)) yield return dir;
    }

    private static string? Find(string relative, string probe)
    {
        foreach (var root in Roots())
        {
            var cand = Path.Combine(root, relative);
            if (File.Exists(Path.Combine(cand, probe))) return cand;
        }
        return null;
    }

    /// <summary>The force-field library directory (data/forcefields), or null.</summary>
    public static string? ForceFields => Find(Path.Combine("data", "forcefields"), "catalogue.json");
    /// <summary>The library of literature many-body potentials (data/potentials: LAMMPS files and their catalogue).</summary>
    public static string? Potentials => Find(Path.Combine("data", "potentials"), "catalogue.json");
    /// <summary>data/reactions/library.json (the reaction schemes), or null.</summary>
    public static string? Reactions => Find(Path.Combine("data", "reactions"), "library.json") is { } d ? Path.Combine(d, "library.json") : null;
    /// <summary>data/reference/polymers.json, or null.</summary>
    public static string? References => Find(Path.Combine("data", "reference"), "polymers.json") is { } d ? Path.Combine(d, "polymers.json") : null;
    public static string? Solvents => Find(Path.Combine("data", "reference"), "solvents.json") is { } d ? Path.Combine(d, "solvents.json") : null;
    /// <summary>data/polymers/library.json (repeat units as SMILES), or null.</summary>
    public static string? Polymers => Find(Path.Combine("data", "polymers"), "library.json") is { } d ? Path.Combine(d, "library.json") : null;
    /// <summary>data/crystals (bulk crystals as CIF, catalogue.json), or null.</summary>
    public static string? Crystals => Find(Path.Combine("data", "crystals"), "catalogue.json");
    /// <summary>The theory manual (data/manual/manual.json).</summary>
    public static string? Manual => Find(Path.Combine("data", "manual"), "manual.json") is { } d ? Path.Combine(d, "manual.json") : null;
    public static string? Fragments => Find(Path.Combine("data", "fragments"), "catalogue.json") is { } d ? Path.Combine(d, "catalogue.json") : null;
    /// <summary>The samples directory, or null.</summary>
    public static string? Samples => Find("samples", "ps_melt.data");
    /// <summary>The caps Python package that Python pipeline steps import (data/python/caps).</summary>
    public static string? Python => Find(Path.Combine("data", "python"), Path.Combine("caps", "runner.py"));
}
