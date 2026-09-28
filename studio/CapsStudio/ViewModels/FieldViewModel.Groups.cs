using System.Collections.ObjectModel;
using System.Globalization;
using System.Text.Json.Nodes;

namespace CapsStudio.ViewModels;

/// <summary>A literature potential of the library (data/potentials/catalogue.json).</summary>
public sealed record PotentialEntry(string Id, string Name, string Style, string File, string[] Elements, string For, string Citation)
{
    public string Label => Id.Length == 0 ? Name : $"{Name} · {string.Join(" ", Elements)}";
    public string Tip => Id.Length == 0 ? "" : $"{For} · {Style} · {Citation}";
    public override string ToString() => Label;
}

/// <summary>One group of Field · by group: molecules, and the force field (or the literature many-body potential) they take.</summary>
public sealed class FieldGroupRow : ObservableObject
{
    public FieldGroupRow(FieldViewModel owner) { Owner = owner; }
    public FieldViewModel Owner { get; }
    public ObservableCollection<FfEntry> Library => Owner.Library;
    private string _name = "", _molecules = "", _file = "";
    private int _kind, _ffIndex, _charges, _style, _units;
    public string Name { get => _name; set => Set(ref _name, value); }
    /// <summary>"1", "2-10, 12" or "rest".</summary>
    public string Molecules { get => _molecules; set => Set(ref _molecules, value); }
    /// <summary>0 a force field of the library, 1 a literature many-body potential (LAMMPS file).</summary>
    public int Kind { get => _kind; set { if (Set(ref _kind, value)) { Raise(nameof(IsForceField)); Raise(nameof(IsPotential)); } } }
    public bool IsForceField => _kind == 0;
    public bool IsPotential => _kind == 1;
    public int FfIndex { get => _ffIndex; set { if (Set(ref _ffIndex, value)) Owner.GroupForceFieldChosen(ForceField); } }
    public int Charges { get => _charges; set => Set(ref _charges, value); }
    public int Style { get => _style; set => Set(ref _style, value); }
    public string File { get => _file; set => Set(ref _file, value); }
    /// <summary>0 the file says (its first line's UNITS:), 1 metal (eV), 2 real (kcal/mol).</summary>
    public int Units { get => _units; set => Set(ref _units, value); }
    public FfEntry? ForceField => _ffIndex >= 0 && _ffIndex < Library.Count ? Library[_ffIndex] : null;
    public List<PotentialEntry> Potentials => Owner.PotentialLibrary;
    private int _pick;
    /// <summary>A potential of the library (0: the user's own file); choosing one sets the style and the file.</summary>
    public int Pick
    {
        get => _pick;
        set
        {
            if (!Set(ref _pick, value) || value <= 0 || value >= Potentials.Count) return;
            var e = Potentials[value];
            var si = Array.IndexOf(FieldViewModel.PotentialStyles, e.Style);
            if (si < 0) throw new InvalidOperationException($"potential style {e.Style} is not in the Studio's list");
            Style = si;
            File = e.File;
            Units = 0;
            Raise(nameof(PickTip));
        }
    }
    public string PickTip => _pick > 0 && _pick < Potentials.Count ? Potentials[_pick].Tip : "Your own file: the NIST Interatomic Potentials Repository, a paper's supplement …";
}

public sealed partial class FieldViewModel
{
    public static readonly string[] GroupKinds = ["Force field", "Literature potential"];
    /// <summary>The many-body styles CAPS writes (manybody.hpp): LAMMPS converts the first ten's files between metal and real units;
    /// AIREBO, AIREBO-M and REBO (carbon, hydrogen) are read in metal units only: the LAMMPS files are then in metal units.</summary>
    public static readonly string[] PotentialStyles = ["tersoff", "tersoff/mod", "tersoff/mod/c", "tersoff/zbl", "sw", "vashishta", "gw", "gw/zbl", "eam/alloy", "eam/fs",
                                                       "airebo", "airebo/morse", "rebo"];
    public static readonly string[] PotentialUnits = ["as the file says", "metal (eV)", "real (kcal/mol)"];
    public static readonly string[] EpsRules = ["geometric √(εᵢεⱼ)", "arithmetic (εᵢ+εⱼ)/2"];
    public static readonly string[] SigmaRules = ["arithmetic (σᵢ+σⱼ)/2", "geometric √(σᵢσⱼ)", "sixth power (class II)"];
    public static readonly string[] Scaling14Modes = ["refuse different scalings", "the first group's for all"];
    public static readonly string[] Cross96Modes = ["refuse 9-6 with 12-6", "12-6 with the 9-6 site's ε and r_min"];

    public ObservableCollection<FieldGroupRow> Groups { get; } = new();
    private List<PotentialEntry>? _potentials;
    /// <summary>The potentials of the library, "own file" first.</summary>
    public List<PotentialEntry> PotentialLibrary => _potentials ??= LoadPotentials();

    private static List<PotentialEntry> LoadPotentials()
    {
        var r = new List<PotentialEntry> { new("", "Own file (choose below)", "", "", [], "", "") };
        if (Paths.Potentials is not { } dir) return r;
        try
        {
            using var js = System.Text.Json.JsonDocument.Parse(System.IO.File.ReadAllText(Path.Combine(dir, "catalogue.json")));
            foreach (var p in js.RootElement.GetProperty("potentials").EnumerateArray())
            {
                string S(string k) => p.TryGetProperty(k, out var v) && v.ValueKind == System.Text.Json.JsonValueKind.String ? v.GetString()! : "";
                r.Add(new PotentialEntry(S("id"), S("name"), S("style"), Path.Combine(dir, S("file")),
                                         p.GetProperty("elements").EnumerateArray().Select(x => x.GetString()!).ToArray(), S("for"), S("citation")));
            }
        }
        catch { }
        return r;
    }

    /// <summary>The library potential for a crystal of these elements: the one covering them with the fewest others
    /// (the first listed on a tie), or none.</summary>
    public int PotentialFor(IReadOnlyCollection<string> elements)
    {
        if (elements.Count == 0) return -1;
        var best = -1;
        for (var k = 1; k < PotentialLibrary.Count; ++k)
            if (elements.All(e => PotentialLibrary[k].Elements.Contains(e)) && (best < 0 || PotentialLibrary[k].Elements.Length < PotentialLibrary[best].Elements.Length))
                best = k;
        return best;
    }
    private bool _groupMode;
    /// <summary>The by-group editor is open.</summary>
    public bool GroupMode { get => _groupMode; set { if (Set(ref _groupMode, value) && value && Groups.Count == 0) SuggestGroups(); } }
    private int _epsRule, _sigmaRule, _scaling14, _cross96;
    public int EpsRule { get => _epsRule; set => Set(ref _epsRule, value); }
    public int SigmaRule { get => _sigmaRule; set { if (Set(ref _sigmaRule, value)) Raise(nameof(EpsRuleEnabled)); } }
    public bool EpsRuleEnabled => _sigmaRule != 2;   // the sixth-power rule sets ε too
    public int Scaling14 { get => _scaling14; set => Set(ref _scaling14, value); }
    public int Cross96 { get => _cross96; set => Set(ref _cross96, value); }
    private string _pairsText = "";
    /// <summary>Explicit cross pairs, one per line: type_a type_b ε (kcal/mol) σ (Å).</summary>
    public string PairsText { get => _pairsText; set => Set(ref _pairsText, value); }
    private bool _grouped;
    /// <summary>The assignment shown was made by group (overrides and entered parameters apply to one force field).</summary>
    public bool IsGrouped { get => _grouped; private set { if (Set(ref _grouped, value)) Raise(nameof(NotGrouped)); } }
    public bool NotGrouped => !_grouped;

    /// <summary>What the rest of the Studio knows about the structure's parts: the held filler, a blend's components
    /// (name, molecules). Set by the main view model.</summary>
    public Func<List<(string Name, string Molecules, bool Crystal, string[] Elements)>>? GroupSuggestions { get; set; }

    /// <summary>A class II force field (PCFF, COMPASS: 9-6) in a group: the cross pairs by its own sixth-power rule
    /// unless another rule was chosen.</summary>
    internal void GroupForceFieldChosen(FfEntry? e)
    {
        if (e != null && _sigmaRule == 0 && (e.Id.Contains("pcff", StringComparison.Ordinal) || e.Id.Contains("compass", StringComparison.Ordinal)))
            SigmaRule = 2;
    }

    public FieldGroupRow AddGroup(string name = "", string molecules = "rest", bool potential = false)
    {
        var g = new FieldGroupRow(this) { Name = name.Length > 0 ? name : $"group {Groups.Count + 1}", Molecules = molecules, Kind = potential ? 1 : 0, FfIndex = Math.Max(0, _ffIndex) };
        Groups.Add(g);
        return g;
    }

    public void RemoveGroup(FieldGroupRow g) => Groups.Remove(g);

    /// <summary>Groups from what is known: the held filler and the rest, a blend's components, else one group for all.</summary>
    public void SuggestGroups()
    {
        Groups.Clear();
        var s = GroupSuggestions?.Invoke() ?? new();
        if (s.Count == 0) s.Add(("all", "rest", false, []));
        foreach (var (name, mols, crystal, elements) in s)
        {
            var g = AddGroup(name, mols, false);
            if (!crystal) continue;
            // a crystal filler: the library's literature potential for its elements, else UFF (every element)
            if (PotentialFor(elements) is var pk && pk > 0) { g.Kind = 1; g.Pick = pk; continue; }
            var k = Library.ToList().FindIndex(e => e.Id == "uff");
            if (k >= 0) g.FfIndex = k;
        }
    }

    /// <summary>The by-group specification (caps_field_assign_groups), or an error message.</summary>
    private (JsonObject? Spec, string Error) GroupSpec()
    {
        if (Groups.Count == 0) return (null, "Add a group.");
        var arr = new JsonArray();
        foreach (var g in Groups)
        {
            if (g.Molecules.Trim().Length == 0) return (null, $"{g.Name}: which molecules? (1, 2-10, or rest)");
            var o = new JsonObject { ["name"] = g.Name.Trim().Length > 0 ? g.Name.Trim() : "group", ["molecules"] = g.Molecules.Trim() };
            if (g.IsPotential)
            {
                if (g.File.Trim().Length == 0) return (null, $"{g.Name}: choose the potential file (LAMMPS's potentials folder, the NIST repository, the paper)");
                var p = new JsonObject { ["style"] = PotentialStyles[Math.Clamp(g.Style, 0, PotentialStyles.Length - 1)], ["file"] = g.File.Trim() };
                if (g.Units > 0) p["units"] = g.Units == 1 ? "metal" : "real";
                o["potential"] = p;
            }
            else
            {
                if (g.ForceField is not { } ff) return (null, $"{g.Name}: choose a force field");
                o["forcefield"] = ff.File;
                o["charges"] = CoreCharges(g.Charges);
            }
            arr.Add(o);
        }
        var pairs = new JsonArray();
        var lineNo = 0;
        foreach (var line in _pairsText.Split('\n'))
        {
            ++lineNo;
            var w = line.Split((char[])[' ', '\t', ','], StringSplitOptions.RemoveEmptyEntries);
            if (w.Length == 0 || w[0].StartsWith('#')) continue;
            if (w.Length != 4 || !double.TryParse(w[2], NumberStyles.Float, Inv, out var eps) || !double.TryParse(w[3], NumberStyles.Float, Inv, out var sig) || eps < 0 || sig <= 0)
                return (null, $"Cross pair line {lineNo}: type_a type_b ε (kcal/mol) σ (Å)");
            pairs.Add(new JsonObject { ["a"] = w[0], ["b"] = w[1], ["eps"] = eps, ["sigma"] = sig });
        }
        return (new JsonObject
        {
            ["groups"] = arr,
            ["eps_rule"] = _epsRule == 1 ? "arithmetic" : "geometric",
            ["sigma_rule"] = _sigmaRule switch { 1 => "geometric", 2 => "sixthpower", _ => "arithmetic" },
            ["scaling14"] = _scaling14 == 1 ? "first" : "refuse",
            ["cross96"] = _cross96 == 1 ? "rmin" : "refuse",
            ["pairs"] = pairs,
        }, "");
    }

    public Task AssignGroups()
    {
        var (spec, err) = GroupSpec();
        if (spec == null) { Log = err; _status(err); return Task.CompletedTask; }
        Recorder?.Invoke(GroupPython(spec));
        _assignedId = "";
        var json = spec.ToJsonString();
        return AssignGroupsRun(json);
    }

    private async Task AssignGroupsRun(string json)
    {
        await Do("Assigned by group", d => d.FieldAssignGroups(json));
        if (IsGrouped && Complete) GroupMode = false;   // the report gets the room; By group opens the groups again
    }

    /// <summary>The same assignment in Python (the macro recorder, Copy as Python).</summary>
    private static string GroupPython(JsonObject spec)
    {
        static string Q(string? s) => "\"" + (s ?? "").Replace("\\", "/").Replace("\"", "\\\"") + "\"";
        var parts = new List<string>();
        foreach (var n in spec["groups"]!.AsArray())
        {
            var g = n!.AsObject();
            var head = $"\"name\": {Q((string?)g["name"])}, \"molecules\": {Q((string?)g["molecules"])}";
            if (g["potential"] is JsonObject p)
                parts.Add($"{{{head}, \"potential\": {{\"style\": {Q((string?)p["style"])}, \"file\": {Q((string?)p["file"])}" + (p["units"] is { } u ? $", \"units\": {Q((string?)u)}" : "") + "}}");
            else
            {
                var c = (int?)g["charges"] ?? 4;
                parts.Add($"{{{head}, \"forcefield\": {Q((string?)g["forcefield"])}, \"charges\": \"{(c switch { 0 => "forcefield", 1 => "gasteiger", 2 => "keep", 3 => "qeq", 5 => "increments", _ => "auto" })}\"}}");
            }
        }
        var pairs = spec["pairs"]!.AsArray().Select(n => string.Format(Inv, "{{\"a\": \"{0}\", \"b\": \"{1}\", \"eps\": {2}, \"sigma\": {3}}}", (string?)n!["a"], (string?)n["b"], (double)n["eps"]!, (double)n["sigma"]!)).ToList();
        return $"doc.field.assign_groups([{string.Join(", ", parts)}], eps_rule=\"{spec["eps_rule"]}\", sigma_rule=\"{spec["sigma_rule"]}\", scaling14=\"{spec["scaling14"]}\", cross96=\"{spec["cross96"]}\"" +
               (pairs.Count > 0 ? $", pairs=[{string.Join(", ", pairs)}]" : "") + ")";
    }
}
