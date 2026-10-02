using System;
using System.Collections.Generic;
using System.Collections.ObjectModel;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Text.Json.Nodes;

namespace CapsStudio.ViewModels;

/// <summary>One unit of the composition calculator's table.</summary>
public sealed record CompRow(string Letter, string Name, string Mol, string Wt, string PerChain);

/// <summary>Polymer builder › composition (a target ratio by moles or weight turned into the unit shares, with what a
/// chain of this DP holds) and your own polymers in the library.</summary>
public sealed partial class MainViewModel
{
    // ---- the composition calculator
    public static readonly string[] CompBases = ["by moles (unit ratio)", "by weight"];
    private int _compBasis;
    public int CompBasis { get => _compBasis; set { if (Set(ref _compBasis, Math.Clamp(value, 0, 1))) CompRefresh(); } }
    public ObservableCollection<CompRow> CompRows { get; } = new();
    private string _compText = "";
    public string CompText { get => _compText; private set => Set(ref _compText, value); }
    public bool CompShown => PolyUnits.Count > 1;

    /// <summary>Mole fractions from the targets: as given by moles; by weight x_i = (w_i / M_i) / Σ (w_j / M_j).</summary>
    private double[]? CompMoleFractions()
    {
        var t = PolyUnits.Select(u => (double)u.Target).ToArray();
        var m = PolyUnits.Select(u => u.Mass).ToArray();
        if (t.Sum() <= 0 || m.Any(x => x <= 0)) return null;
        var x = _compBasis == 1 ? t.Select((w, i) => w / m[i]).ToArray() : t;
        var s = x.Sum();
        return x.Select(v => v / s).ToArray();
    }

    public void CompRefresh()
    {
        CompRows.Clear();
        Raise(nameof(CompShown));
        var x = CompMoleFractions();
        if (x == null) { CompText = PolyUnits.Any(u => u.Mass <= 0) ? "Every unit needs a valid SMILES (its mass) first." : "Give each unit a share."; return; }
        var inv = CultureInfo.InvariantCulture;
        var m = PolyUnits.Select(u => u.Mass).ToArray();
        var mbar = x.Select((v, i) => v * m[i]).Sum();   // mean unit mass
        var dp = _growDp;
        // exact counts by largest remainders (as the exact-composition sequence builds them)
        var exact = x.Select(v => (int)Math.Floor(dp * v)).ToArray();
        foreach (var k in x.Select((v, i) => (r: dp * v - Math.Floor(dp * v), i)).OrderByDescending(p => p.r).Take(dp - exact.Sum()).Select(p => p.i).ToList()) ++exact[k];
        for (var i = 0; i < PolyUnits.Count; ++i)
        {
            var u = PolyUnits[i];
            var sd = Math.Sqrt(dp * x[i] * (1 - x[i]));
            CompRows.Add(new CompRow(u.Letter, u.Name.Length > 0 ? u.Name : u.Smiles,
                (100 * x[i]).ToString("0.0", inv) + " %", (100 * x[i] * m[i] / mbar).ToString("0.0", inv) + " %",
                $"{exact[i]} · {(dp * x[i]).ToString("0.0", inv)} ± {sd.ToString("0.0", inv)}"));
        }
        CompText = string.Format(inv, "DP {0}: a chain of about {1:0} g/mol (mean unit {2:0.0} g/mol). \"Random\" draws each unit on its own, so chains scatter around the ratio (± one standard deviation shown); \"Random, exact composition\" gives every chain the exact counts in random order.",
            dp, dp * mbar, mbar);
    }

    /// <summary>The unit shares from the calculator (mole fractions), and a random sequence when none is chosen.</summary>
    public void ApplyComposition(bool exact)
    {
        var x = CompMoleFractions();
        if (x == null) return;
        for (var i = 0; i < PolyUnits.Count; ++i) PolyUnits[i].Weight = Math.Round((decimal)x[i], 4);
        PolySequence = exact ? 6 : 3;
        Status = $"Unit shares set from the {(_compBasis == 1 ? "weight" : "mole")} ratio: " + string.Join(" : ", PolyUnits.Select(u => u.Weight.ToString("0.###", CultureInfo.InvariantCulture)));
        CompRefresh();
    }

    // ---- your polymers, beside the library's
    // tests and screenshots (a settings override) keep theirs beside it, never in the user's folder
    public static string UserPolymerFile => AppSettings.Override != null
        ? Path.Combine(Path.GetDirectoryName(AppSettings.Override)!, "caps-polymers.json")
        : Path.Combine(AppSettings.Folder, "polymers.json");

    private void LoadUserPolymers()
    {
        try
        {
            if (!File.Exists(UserPolymerFile)) return;
            var j = JsonNode.Parse(File.ReadAllText(UserPolymerFile));
            foreach (var p in (j?["polymers"] as JsonArray ?? []).OfType<JsonObject>())
            {
                var id = (string?)p["id"] ?? "";
                _libById[id] = p;
                PolymerLibrary.Add(new LibraryEntry(id, (string?)p["name"] ?? id, (string?)p["smiles"] ?? "", UserRubber(p, [(string?)p["smiles"] ?? ""]), false, null, true));
            }
            foreach (var c in (j?["copolymers"] as JsonArray ?? []).OfType<JsonObject>())
                PolymerLibrary.Add(new LibraryEntry((string?)c["id"] ?? "", (string?)c["name"] ?? "", "", UserRubber(c, UnitSmiles(c)), true, c, true));
        }
        catch (Exception e) { Status = "Your polymers (" + UserPolymerFile + "): " + e.Message; }
    }

    /// <summary>Your entry is a rubber when saved so (its "rubber" tag), not when saved as "not-rubber"; an entry saved
    /// before the tag existed is judged from its units.</summary>
    private bool UserRubber(JsonObject o, IEnumerable<string> units)
    {
        var tags = (o["tags"] as JsonArray)?.Select(x => (string?)x).ToList() ?? [];
        if (tags.Contains("rubber")) return true;
        if (tags.Contains("not-rubber")) return false;
        return units.Any(IsRubberUnit);
    }

    private IEnumerable<string> UnitSmiles(JsonObject c) =>
        (c["units"] as JsonArray ?? []).Select(x => (string?)x ?? "").Select(u => _libById.TryGetValue(u, out var p) ? (string?)p["smiles"] ?? "" : u);

    /// <summary>The repeat units of the library's rubbers (each rubber homopolymer, and every unit of a rubber copolymer),
    /// written alike: "[*]" as "*", cis/trans marks dropped (a unit drawn without them is still the same backbone).</summary>
    private HashSet<string>? _rubberUnits;
    internal static string UnitKey(string smiles) => smiles.Trim().Replace("[*]", "*").Replace("/", "").Replace("\\", "");
    public bool IsRubberUnit(string smiles)
    {
        if (_rubberUnits == null && PolymerLibrary.Count > 0)
        {
            _rubberUnits = [];
            foreach (var e in PolymerLibrary.Where(e => e.Rubber && !e.User))
                foreach (var u in e.Copolymer && e.Preset != null ? UnitSmiles(e.Preset) : [e.Smiles])
                    if (u.Length > 0) _rubberUnits.Add(UnitKey(u));
        }
        return smiles.Length > 0 && _rubberUnits != null && _rubberUnits.Contains(UnitKey(smiles));
    }

    /// <summary>Whether the polymer on the builder is saved as a rubber: guessed from its units (a unit of one of the
    /// library's rubbers) until you set it.</summary>
    private bool _polyRubber, _polyRubberSet;
    public bool PolySaveRubber { get => _polyRubber; set { _polyRubberSet = true; Set(ref _polyRubber, value); } }
    public string PolyRubberTip => _polyRubberSet ? "Set by you" : _polyRubber
        ? "Guessed: a unit is a unit of one of the library's rubbers (" + string.Join(", ", PolyUnits.Where(u => u.Ok && IsRubberUnit(u.Smiles)).Select(u => u.Name.Length > 0 ? u.Name : u.Letter)) + ")"
        : "Guessed: no unit is a unit of one of the library's rubbers; tick it when this is a rubber";
    private void GuessPolyRubber()
    {
        if (!_polyRubberSet) Set(ref _polyRubber, PolyUnits.Any(u => u.Ok && IsRubberUnit(u.Smiles)), nameof(PolySaveRubber));
        Raise(nameof(PolyRubberTip));
    }

    private JsonObject ReadUserFile()
    {
        try { if (File.Exists(UserPolymerFile) && JsonNode.Parse(File.ReadAllText(UserPolymerFile)) is JsonObject o) return o; } catch { }
        return new JsonObject { ["format"] = "caps-polymers", ["polymers"] = new JsonArray(), ["copolymers"] = new JsonArray() };
    }

    /// <summary>The polymer on the builder saved to your library: one unit as a repeat unit, several as a copolymer with
    /// its sequence and shares (or blocks, or pattern).</summary>
    public string SavePolymerToLibrary()
    {
        var units = PolyUnits.Where(u => u.Ok).ToList();
        if (units.Count == 0) return "Nothing to save: give the repeat unit a valid SMILES.";
        var name = (_polyName.Trim().Length > 0 ? _polyName.Trim() : units.Count == 1 ? units[0].Name : "").Trim();
        if (name.Length == 0) return "Give the polymer a name first.";
        var file = ReadUserFile();
        var polys = file["polymers"] as JsonArray ?? new JsonArray();
        var copos = file["copolymers"] as JsonArray ?? new JsonArray();
        file["polymers"] = polys;
        file["copolymers"] = copos;
        int Next(JsonArray a, string prefix) =>
            a.OfType<JsonObject>().Select(o => int.TryParse(((string?)o["id"] ?? "").Replace(prefix, ""), out var k) ? k : 0).DefaultIfEmpty(0).Max() + 1;
        // a name already there is replaced
        void Drop(JsonArray a) { foreach (var o in a.OfType<JsonObject>().Where(o => string.Equals((string?)o["name"], name, StringComparison.OrdinalIgnoreCase)).ToList()) a.Remove(o); }
        Drop(polys);
        Drop(copos);
        JsonArray Tags() => new("user", _polyRubber ? "rubber" : "not-rubber");
        if (units.Count == 1)
            polys.Add(new JsonObject { ["id"] = $"U{Next(polys, "U"):000}", ["name"] = name, ["smiles"] = units[0].Smiles, ["tags"] = Tags() });
        else
        {
            // each unit saved with its SMILES (a unit also saved alone keeps its own entry)
            var c = new JsonObject
            {
                ["id"] = $"UC{Next(copos, "UC"):000}", ["name"] = name,
                ["units"] = new JsonArray(units.Select(u => (JsonNode)u.Smiles).ToArray()),
                ["unit_names"] = new JsonArray(units.Select(u => (JsonNode)(u.Name.Length > 0 ? u.Name : u.Smiles)).ToArray()),
                ["sequence"] = SequenceIds[Math.Clamp(_polySeq, 0, SequenceIds.Length - 1)], ["tags"] = Tags(),
            };
            if (_polySeq is 3 or 6) c["weights"] = new JsonArray(units.Select(u => (JsonNode)(double)u.Weight).ToArray());
            if (_polySeq == 2) c["blocks"] = new JsonArray(units.Select(u => (JsonNode)(double)u.Block).ToArray());
            if (_polySeq == 5) c["pattern"] = PolyPattern;
            copos.Add(c);
        }
        Directory.CreateDirectory(Path.GetDirectoryName(UserPolymerFile)!);
        File.WriteAllText(UserPolymerFile, file.ToJsonString(new System.Text.Json.JsonSerializerOptions { WriteIndented = true }));
        ReloadPolymerLibrary();
        // shown at once: the search and the Rubbers filter let it through
        if (_libRubber && !_polyRubber) LibraryRubberOnly = false;
        if (_libQuery.Trim().Length > 0 && !name.Contains(_libQuery.Trim(), StringComparison.OrdinalIgnoreCase)) LibraryQuery = "";
        return $"Saved \"{name}\" to your polymers{(_polyRubber ? " as a rubber" : "")} ({UserPolymerFile})";
    }

    /// <summary>One of your polymers removed from your library (the built-in ones stay).</summary>
    public string RemoveUserPolymer(LibraryEntry e)
    {
        if (!e.User) return "Only your own polymers can be removed.";
        var file = ReadUserFile();
        foreach (var key in new[] { "polymers", "copolymers" })
            if (file[key] is JsonArray a)
                foreach (var o in a.OfType<JsonObject>().Where(o => (string?)o["id"] == e.Id).ToList()) a.Remove(o);
        File.WriteAllText(UserPolymerFile, file.ToJsonString(new System.Text.Json.JsonSerializerOptions { WriteIndented = true }));
        ReloadPolymerLibrary();
        return $"Removed \"{e.Name}\" from your polymers";
    }

    private void ReloadPolymerLibrary()
    {
        PolymerLibrary.Clear();
        _rubberUnits = null;
        LoadPolymerLibrary();
    }
}
