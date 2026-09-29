using System.Collections.ObjectModel;
using System.ComponentModel;
using System.Globalization;
using System.Text.Json.Nodes;
using CapsStudio.Interop;

namespace CapsStudio.ViewModels;

/// <summary>A repeat unit or copolymer preset from data/polymers/library.json.</summary>
public sealed record LibraryEntry(string Id, string Name, string Smiles, bool Rubber, bool Copolymer, JsonObject? Preset, bool User = false)
{
    public string Kind => User ? (Copolymer ? "yours · copolymer" : "yours") : Copolymer ? "copolymer" : Rubber ? "rubber" : "";
    public bool HasKind => Kind.Length > 0;
}

/// <summary>One repeat unit of the chain being built (A, B, C …).</summary>
public sealed class PolyUnit : INotifyPropertyChanged
{
    public event PropertyChangedEventHandler? PropertyChanged;
    private void Raise(string n) => PropertyChanged?.Invoke(this, new PropertyChangedEventArgs(n));
    public Action? Changed;

    private string _letter = "A";
    public string Letter { get => _letter; set { _letter = value; Raise(nameof(Letter)); } }
    private string _name = "", _smiles = "", _info = "", _error = "";
    private decimal _weight = 1, _block = 10;
    public string Name { get => _name; set { _name = value; Raise(nameof(Name)); Changed?.Invoke(); } }
    public string Smiles
    {
        get => _smiles;
        set { _smiles = value; Raise(nameof(Smiles)); Describe(); Changed?.Invoke(); }
    }
    public string Info { get => _info; private set { _info = value; Raise(nameof(Info)); } }
    public string Error { get => _error; private set { _error = value; Raise(nameof(Error)); Raise(nameof(HasError)); Raise(nameof(Ok)); } }
    public bool HasError => _error.Length > 0;
    public bool Ok => _error.Length == 0 && _smiles.Length > 0;
    public int Stereocentres { get; private set; }
    /// <summary>The unit's molar mass in the chain (g/mol), for the composition calculator.</summary>
    public double Mass { get; private set; }
    /// <summary>The atoms at the head (first *) and tail (second *), where the chain's end groups bond.</summary>
    public string HeadElement { get; private set; } = "";
    public string TailElement { get; private set; } = "";
    /// <summary>Whether the atom at the first / second * carries =O (written C(=O) next to the *).</summary>
    public bool EndCarbonyl(bool head)
    {
        var t = _smiles.Replace("[*]", "*");
        var stars = t.Select((c, i) => (c, i)).Where(p => p.c == '*').Select(p => p.i).ToList();
        if (stars.Count < 2) return false;
        if (head)
        {
            var after = t[(stars[0] + 1)..];   // *C(=O)… : the atom right after the head *
            return after.StartsWith("C(=O)");
        }
        var before = t[..stars[1]];   // …C(=O)* or …C(*)=O
        return before.EndsWith("C(=O)") || before.EndsWith("C(=O)(") || t[stars[1]..].StartsWith("*)=O") && before.EndsWith("C(");
    }
    private decimal _target = 50;
    /// <summary>The share this unit should have (mole or weight %, as the composition calculator is set).</summary>
    public decimal Target { get => _target; set { _target = Math.Max(0, value); Raise(nameof(Target)); TargetChanged?.Invoke(); } }
    public Action? TargetChanged;
    public decimal Weight { get => _weight; set { _weight = Math.Max(0, value); Raise(nameof(Weight)); Changed?.Invoke(); } }
    public decimal Block { get => _block; set { _block = Math.Max(1, value); Raise(nameof(Block)); Changed?.Invoke(); } }

    private void Describe()
    {
        try
        {
            var r = JsonNode.Parse(CapsDocument.UnitInfo(_smiles))!;
            if ((bool?)r["ok"] == true)
            {
                Stereocentres = (int?)r["stereocentres"] ?? 0;
                Mass = (double?)r["mass"] ?? 0;
                HeadElement = (string?)r["head_element"] ?? "";
                TailElement = (string?)r["tail_element"] ?? "";
                Info = string.Format(CultureInfo.InvariantCulture, "{0} · {1:F2} g/mol · head {2} → tail {3}{4}", Sub((string?)r["formula"] ?? ""),
                    (double?)r["mass"] ?? 0, r["head_element"], r["tail_element"], Stereocentres > 0 ? " · stereocentre" : "");
                if ((string?)r["note"] is { } note) Info += "\n" + char.ToUpper(note[0]) + note[1..];
                Error = "";
            }
            else { Info = ""; Error = (string?)r["error"] ?? "not a repeat unit"; }
        }
        catch (Exception e) { Info = ""; Error = e.Message; }
    }

    internal static string Sub(string f) => string.Concat(f.Select(c => c is >= '0' and <= '9' ? (char)('₀' + (c - '0')) : c));
}

/// <summary>The polymer builder (design/boards/PolymerBuilder): repeat units from the library or SMILES, the sequence
/// (homopolymer, alternating, block, random, gradient, pattern) over any number of units, tacticity, DP; one chain
/// previewed; sent to Grow for a cell.</summary>
public sealed partial class MainViewModel
{
    public bool IsPolymer => _module == 13;
    public static readonly string[] SequenceKinds = ["Homopolymer", "Alternating", "Block", "Random", "Gradient", "Pattern", "Random, exact composition"];
    private static readonly string[] SequenceIds = ["homopolymer", "alternating", "block", "random", "gradient", "pattern", "shuffled"];

    public List<LibraryEntry> PolymerLibrary { get; } = new();
    public ObservableCollection<LibraryEntry> LibraryShown { get; } = new();
    private string _libQuery = "";
    private bool _libRubber;
    public string LibraryQuery { get => _libQuery; set { if (Set(ref _libQuery, value)) FilterLibrary(); } }
    public bool LibraryRubberOnly { get => _libRubber; set { if (Set(ref _libRubber, value)) FilterLibrary(); } }
    public string LibraryCount => $"{PolymerLibrary.Count(e => !e.Copolymer)} repeat units · {PolymerLibrary.Count(e => e.Copolymer)} copolymers";

    public ObservableCollection<PolyUnit> PolyUnits { get; } = new();
    private int _polySeq;
    public int PolySequence { get => _polySeq; set { if (Set(ref _polySeq, value)) { RaisePolyKinds(); PolyChanged(); } } }
    public bool PolyShowWeights => _polySeq is 3 or 6;
    public bool PolyShowBlocks => _polySeq == 2;
    public bool PolyShowPattern => _polySeq == 5;
    // ---- architecture: linear, branched (random side chains), star, comb, dendrimer
    private static readonly string[] ArchIds = ["linear", "branched", "star", "comb", "dendrimer"];
    private int _polyArch;
    private decimal _polyArms = 4, _polyArmDp = 5, _polySpacing = 4, _polyBranchP = 0.1m;
    public int PolyArch
    {
        get => _polyArch;
        set
        {
            if (!Set(ref _polyArch, Math.Clamp(value, 0, 4))) return;
            Raise(nameof(PolyIsLinear)); Raise(nameof(PolyIsStar)); Raise(nameof(PolyHasCore)); Raise(nameof(PolyIsDendrimer)); Raise(nameof(PolySideLabel)); Raise(nameof(PolyHasSideChains)); Raise(nameof(PolyIsComb)); Raise(nameof(PolyIsBranched));
            Raise(nameof(PolyArchNote)); Raise(nameof(PolyPreviewCaption));
            PolyChanged();
        }
    }
    public bool PolyIsLinear => _polyArch == 0;
    public string PolyPreviewCaption => _polyArch == 0 ? "Preview · one chain" : "Preview · one " + ArchIds[_polyArch] + " molecule";
    public bool PolyIsBranched => _polyArch == 1;
    public bool PolyIsStar => _polyArch == 2;
    public bool PolyIsComb => _polyArch == 3;
    public bool PolyIsDendrimer => _polyArch == 4;
    public bool PolyHasCore => _polyArch is 2 or 4;
    public bool PolyHasSideChains => _polyArch is 1 or 3 or 4;
    public string PolySideLabel => _polyArch == 4 ? "Units per branch" : "Units per side chain";
    private decimal _polyGenerations = 2;
    public decimal PolyGenerations { get => _polyGenerations; set { if (Set(ref _polyGenerations, Math.Clamp(value, 1, 6))) { Raise(nameof(PolyArchNote)); PolyChanged(); } } }
    public decimal PolyArms { get => _polyArms; set { if (Set(ref _polyArms, Math.Clamp(value, 3, 4))) { Raise(nameof(PolyArchNote)); PolyChanged(); } } }
    public decimal PolyArmDp { get => _polyArmDp; set { if (Set(ref _polyArmDp, Math.Clamp(value, 1, 500))) PolyChanged(); } }
    public decimal PolySpacing { get => _polySpacing; set { if (Set(ref _polySpacing, Math.Clamp(value, 1, 500))) PolyChanged(); } }
    public decimal PolyBranchP { get => _polyBranchP; set { if (Set(ref _polyBranchP, Math.Clamp(value, 0, 1))) PolyChanged(); } }
    public string PolyArchNote => _polyArch switch
    {
        2 => "Arms of DP units on one core carbon, the first arm's head atom (star SBR and BR are coupled on silicon or tin: here the core is carbon). Four arms need a CH₂ or CH₃ head.",
        3 => "A side chain of the chain's own units on every n-th backbone unit, on a hydrogen of its head or tail atom (or of the unit beside it when that has clearly more room).",
        1 => "Long-chain branches: each backbone unit carries a side chain with this probability.",
        4 => string.Format(System.Globalization.CultureInfo.InvariantCulture,
                 "A star core of {0} arms of DP units; every free end splits in two branches, generation after generation: {1} branches, {2} end groups per molecule. Both branches hang on the end unit (the second on the unit before it when the end is crowded).",
                 (int)_polyArms, (int)_polyArms * ((1 << ((int)_polyGenerations + 1)) - 2), (int)_polyArms * (1 << (int)_polyGenerations)),
        _ => "",
    } + (_polyArch == 0 ? "" : " Branch points are crowded: they may grow at a reduced contact scale; relax with push-off before dynamics.");
    private string _polyPattern = "AB";
    public string PolyPattern { get => _polyPattern; set { if (Set(ref _polyPattern, value)) PolyChanged(); } }
    private string _polyName = "Polystyrene";
    public string PolyName { get => _polyName; set { if (Set(ref _polyName, value)) PolyChanged(); } }
    private string _polyPreview = "", _polyStrip = "", _polyError = "";
    public string PolyPreview { get => _polyPreview; private set => Set(ref _polyPreview, value); }
    public string PolyError { get => _polyError; private set { if (Set(ref _polyError, value)) Raise(nameof(PolyHasError)); } }
    public bool PolyHasError => _polyError.Length > 0;
    public int[] PolyStripUnits { get; private set; } = [];
    public event Action? PolyStripChanged;
    private CapsDocument? _polyDoc;
    public CapsDocument? PolyDoc { get => _polyDoc; private set => Set(ref _polyDoc, value); }
    private bool _polyBuilding;
    public bool PolyBuilding { get => _polyBuilding; private set => Set(ref _polyBuilding, value); }
    public bool PolyCanAdd => PolyUnits.Count < 8;

    private void RaisePolyKinds() { Raise(nameof(PolyShowWeights)); Raise(nameof(PolyShowBlocks)); Raise(nameof(PolyShowPattern)); }

    public void LoadPolymerLibrary()
    {
        if (PolymerLibrary.Count > 0) return;
        try
        {
            var path = Paths.Polymers;
            if (path == null) return;
            var j = JsonNode.Parse(File.ReadAllText(path))!;
            var byId = new Dictionary<string, JsonObject>();
            foreach (var p in (j["polymers"] as JsonArray ?? []).OfType<JsonObject>())
            {
                var id = (string?)p["id"] ?? "";
                byId[id] = p;
                var tags = (p["tags"] as JsonArray)?.Select(x => (string?)x).ToList() ?? [];
                PolymerLibrary.Add(new LibraryEntry(id, (string?)p["name"] ?? id, (string?)p["smiles"] ?? "", tags.Contains("rubber"), false, null));
            }
            foreach (var c in (j["copolymers"] as JsonArray ?? []).OfType<JsonObject>())
            {
                var tags = (c["tags"] as JsonArray)?.Select(x => (string?)x).ToList() ?? [];
                PolymerLibrary.Add(new LibraryEntry((string?)c["id"] ?? "", (string?)c["name"] ?? "", "", tags.Contains("rubber"), true, c));
            }
            _libById = byId;
        }
        catch (Exception e) { Status = "Polymer library: " + e.Message; }
        LoadUserPolymers();
        FilterLibrary();
        Raise(nameof(LibraryCount));
        if (PolyUnits.Count == 0) AddPolyUnit("Styrene", "*CC(*)c1ccccc1");
    }
    private Dictionary<string, JsonObject> _libById = new();

    private void FilterLibrary()
    {
        LibraryShown.Clear();
        var q = _libQuery.Trim().ToLowerInvariant();
        foreach (var e in PolymerLibrary)
        {
            if (_libRubber && !e.Rubber) continue;
            if (q.Length > 0 && !e.Name.ToLowerInvariant().Contains(q) && !e.Smiles.ToLowerInvariant().Contains(q) && !e.Id.ToLowerInvariant().Contains(q)) continue;
            LibraryShown.Add(e);
        }
    }

    public void AddPolyUnit(string name = "", string smiles = "*CC*")
    {
        if (PolyUnits.Count >= 8) return;
        var u = new PolyUnit { Letter = ((char)('A' + PolyUnits.Count)).ToString(), Changed = PolyChanged, TargetChanged = CompRefresh };
        u.Name = name;
        u.Smiles = smiles;
        PolyUnits.Add(u);
        Raise(nameof(PolyCanAdd));
        PolyChanged();
    }

    public void RemovePolyUnit(PolyUnit u)
    {
        if (PolyUnits.Count <= 1) return;
        PolyUnits.Remove(u);
        for (var k = 0; k < PolyUnits.Count; ++k) PolyUnits[k].Letter = ((char)('A' + k)).ToString();
        Raise(nameof(PolyCanAdd));
        PolyChanged();
    }

    /// <summary>A library entry: a repeat unit goes into the selected slot (or A), a copolymer preset replaces all.</summary>
    public void UseLibrary(LibraryEntry e, PolyUnit? slot)
    {
        if (!e.Copolymer)
        {
            var target = slot ?? PolyUnits.FirstOrDefault();
            if (target == null) AddPolyUnit(e.Name, e.Smiles);
            else { target.Name = e.Name; target.Smiles = e.Smiles; }
            if (PolyUnits.Count == 1) _polyName = e.Name;
            Raise(nameof(PolyName));
            PolyChanged();
            return;
        }
        var p = e.Preset!;
        PolyUnits.Clear();
        var units = (p["units"] as JsonArray ?? []).Select(x => (string?)x ?? "").ToList();
        var names = (p["unit_names"] as JsonArray)?.Select(x => (string?)x ?? "").ToList();
        for (var k = 0; k < units.Count; ++k)
        {
            var u = units[k];
            var (nm, smi) = _libById.TryGetValue(u, out var lp) ? ((string?)lp["name"] ?? u, (string?)lp["smiles"] ?? "") : (names != null && k < names.Count ? names[k] : u, u);
            AddPolyUnit(nm, smi);
        }
        var seq = Array.IndexOf(SequenceIds, (string?)p["sequence"] ?? "homopolymer");
        _polySeq = Math.Max(0, seq);
        RaisePolyKinds();
        Raise(nameof(PolySequence));
        if (p["weights"] is JsonArray w)
        {
            // the library's shares are unit (mole) fractions: the calculator starts from them
            for (var k = 0; k < Math.Min(w.Count, PolyUnits.Count); ++k)
            {
                PolyUnits[k].Weight = (decimal)((double?)w[k] ?? 1);
                PolyUnits[k].Target = Math.Round(100 * PolyUnits[k].Weight / Math.Max(1e-9m, w.Sum(x => (decimal)((double?)x ?? 1))), 2);
            }
            _compBasis = 0;
            Raise(nameof(CompBasis));
        }
        if (p["blocks"] is JsonArray b)
            for (var k = 0; k < Math.Min(b.Count, PolyUnits.Count); ++k) PolyUnits[k].Block = (decimal)((double?)b[k] ?? 10);
        if ((string?)p["pattern"] is { Length: > 0 } pat) _polyPattern = pat;
        Raise(nameof(PolyPattern));
        _polyName = e.Name;
        Raise(nameof(PolyName));
        PolyChanged();
    }

    private string _polyEndsText = "";
    /// <summary>What the chain ends become: the atom each end group bonds to, and a warning for a chemically wrong choice.</summary>
    public string PolyEndsText { get => _polyEndsText; private set => Set(ref _polyEndsText, value); }
    private string EndsText()
    {
        if (PolyStripUnits.Length == 0 || PolyUnits.Count == 0) return "";
        var first = PolyUnits[Math.Clamp(PolyStripUnits[0], 0, PolyUnits.Count - 1)];
        var last = PolyUnits[Math.Clamp(PolyStripUnits[^1], 0, PolyUnits.Count - 1)];
        string End(bool head, PolyUnit u, int cap)
        {
            var el = head ? u.HeadElement : u.TailElement;
            var carbonyl = u.EndCarbonyl(head);
            var atom = carbonyl ? "C=O" : el;
            var g = EndGroups[cap];
            var result = (g, el, carbonyl) switch
            {
                ("hydrogen", "O", _) => " → –OH",
                ("hydrogen", "N", _) => " → –NH",
                ("hydrogen", "C", true) => " → aldehyde –CHO",
                ("hydroxyl", "C", true) => " → acid –COOH",
                ("hydroxyl", "C", false) => " → alcohol –C–OH",
                ("methyl", "C", true) => " → methyl ketone",
                ("amine", "C", true) => " → amide –C(=O)NH₂",
                _ => "",
            };
            var warn = (g is "hydroxyl" or "amine" && el is "O" or "N" or "S") ? $"  ⚠ {el}–{(g == "hydroxyl" ? "O" : "N")}: this end is already {el}; hydrogen gives –{el}H"
                : (g == "carboxyl" && carbonyl) ? "  ⚠ C(=O)–COOH: this end is already a carbonyl; hydroxyl gives –COOH" : "";
            return $"{(head ? "Head" : "Tail")}: {g} on {atom}{result}{warn}";
        }
        return End(true, first, _headCap) + "\n" + End(false, last, _tailCap);
    }

    // end groups in place of the chain ends' hydrogens (core chain_end_smiles)
    public static readonly string[] EndGroups = ["hydrogen", "methyl", "ethyl", "tert-butyl", "sec-butyl", "phenyl", "hydroxyl", "carboxyl", "vinyl", "amine"];
    private int _headCap, _tailCap;
    public int HeadCap { get => _headCap; set { if (Set(ref _headCap, Math.Clamp(value, 0, EndGroups.Length - 1))) PolyChanged(); } }
    public int TailCap { get => _tailCap; set { if (Set(ref _tailCap, Math.Clamp(value, 0, EndGroups.Length - 1))) PolyChanged(); } }

    // how units join (core Linkage): head-to-tail, head-to-head (every second unit reversed), random reversals
    public static readonly string[] LinkageIds = ["head-to-tail", "head-to-head", "random"];
    private int _polyLinkage;
    private decimal _polyInversion = 5;
    public int PolyLinkage { get => _polyLinkage; set { if (Set(ref _polyLinkage, Math.Clamp(value, 0, 2))) { Raise(nameof(PolyRandomLinkage)); PolyChanged(); } } }
    public bool PolyRandomLinkage => _polyLinkage == 2;
    /// <summary>Random linkage: the share of reversed units, % (regio-defects: PVDF ≈ 3–6 %).</summary>
    public decimal PolyInversion { get => _polyInversion; set { if (Set(ref _polyInversion, Math.Clamp(value, 0, 100))) PolyChanged(); } }

    /// <summary>The chain as the core reads it (caps_chain_preview / caps_grow_chains).</summary>
    public string PolySpecJson(int? dp = null)
    {
        var o = new JsonObject
        {
            ["units"] = new JsonArray(PolyUnits.Select(u => (JsonNode)new JsonObject { ["name"] = u.Name, ["smiles"] = u.Smiles }).ToArray()),
            ["sequence"] = SequenceIds[Math.Clamp(_polySeq, 0, SequenceIds.Length - 1)],
            ["dp"] = dp ?? _growDp,
            ["weights"] = new JsonArray(PolyUnits.Select(u => (JsonNode)(double)u.Weight).ToArray()),
            ["blocks"] = new JsonArray(PolyUnits.Select(u => (JsonNode)(int)u.Block).ToArray()),
            ["pattern"] = _polyPattern,
        };
        if (_headCap > 0) o["head_cap"] = EndGroups[_headCap];
        if (_tailCap > 0) o["tail_cap"] = EndGroups[_tailCap];
        if (_polyLinkage > 0)
        {
            o["linkage"] = LinkageIds[_polyLinkage];
            if (_polyLinkage == 2) o["inversion"] = (double)_polyInversion / 100;
        }
        if (_polyArch != 0)
        {
            o["architecture"] = ArchIds[_polyArch];
            o["arms"] = (int)_polyArms;
            o["arm_dp"] = (int)_polyArmDp;
            o["spacing"] = (int)_polySpacing;
            o["branch_probability"] = (double)_polyBranchP;
            o["generations"] = (int)_polyGenerations;
        }
        // unit templates cleaned with the default force field when it types them (GAFF2 unless Settings says otherwise)
        if (CleanChoices.FirstOrDefault(c => c.File != null && Path.GetFileNameWithoutExtension(c.File) == _settings.ForceField)?.File is { } clean) o["forcefield"] = clean;
        else if (CleanChoices.FirstOrDefault(c => c.File != null)?.File is { } first) o["forcefield"] = first;
        return o.ToJsonString();
    }

    private void PolyChanged()
    {
        CompRefresh();
        if (PolyUnits.Count == 0) return;
        if (PolyUnits.FirstOrDefault(u => !u.Ok) is { } bad) { PolyError = $"{bad.Letter}: {(bad.HasError ? bad.Error : "no SMILES")}"; return; }
        try
        {
            var r = JsonNode.Parse(CapsDocument.ChainPreview(PolySpecJson(), (ulong)(_polyDraw = PolySeed.Take())))!;
            if ((bool?)r["ok"] != true) { PolyError = (string?)r["error"] ?? "cannot build this chain"; return; }
            PolyError = "";
            PolyStripUnits = (r["sequence"] as JsonArray ?? []).Select(x => (int?)x ?? 0).ToArray();
            SchedulePolyPreview();
            var counts = PolyUnits.Select((u, k) => PolyStripUnits.Count(x => x == k)).ToArray();
            PolyPreview = string.Format(CultureInfo.InvariantCulture, "{0} · {1:N0} g/mol · {2:N0} atoms per chain · {3}",
                PolyUnit.Sub((string?)r["formula"] ?? ""), (double?)r["mass"] ?? 0, (int?)r["atoms"] ?? 0,
                string.Join(" · ", PolyUnits.Select((u, k) => $"{u.Letter} {counts[k]}")));
            var reversed = (r["inverted"] as JsonArray ?? []).Count(x => (int?)x == 1);
            if (_polyLinkage > 0) PolyPreview += $" · {reversed} reversed";
            _polyAtoms = (int?)r["atoms"] ?? 0;
            _polyMass = (double?)r["mass"] ?? 0;
            PolyEndsText = EndsText();
            if (r["molecule"] is JsonObject m)   // a branched molecule: its arms, atoms and mass
            {
                var arms = (double?)m["branches"] ?? (double?)m["arms"] ?? 0;   // a dendrimer: its branches besides the core arms
                _polyAtoms = (int)Math.Round((double?)m["atoms"] ?? _polyAtoms);
                _polyMass = (double?)m["mass"] ?? _polyMass;
                PolyPreview += string.Format(CultureInfo.InvariantCulture, "\nper molecule: {0} {1} · {4}{2:N0} atoms · {4}{3:N0} g/mol",
                    _polyArch == 1 ? "≈ " + arms.ToString("0.#", CultureInfo.InvariantCulture) : arms.ToString("0", CultureInfo.InvariantCulture),
                    _polyArch switch { 2 => "more arms", 4 => $"branches on {(int)_polyArms} core arms", _ => "side chains" }, _polyAtoms, _polyMass, PolyUnits.Count > 1 || _polyArch == 1 ? "≈ " : "");   // each arm draws its own sequence
            }
            PolyStripChanged?.Invoke();
        }
        catch (Exception e) { PolyError = e.Message; }
    }
    private int _polyAtoms;
    private double _polyMass;

    /// <summary>One chain (the preview) grown on its own.</summary>
    // the 3D chain follows the builder: rebuilt 0.6 s after the last change (while the page is open), and again after a
    // build that was running when something changed
    public bool PolyAutoPreview { get; set; }
    private int _polyPreviewGen;
    private bool _polyPreviewAgain;
    private void SchedulePolyPreview()
    {
        if (!PolyAutoPreview || _module != 13) return;   // the Polymer builder on screen
        var gen = ++_polyPreviewGen;
        Avalonia.Threading.DispatcherTimer.RunOnce(() =>
        {
            if (gen != _polyPreviewGen || !PolyAutoPreview || _module != 13 || PolyHasError) return;
            if (_polyBuilding) { _polyPreviewAgain = true; return; }
            _ = BuildPolyPreview();
        }, TimeSpan.FromMilliseconds(600));
    }

    public async Task BuildPolyPreview()
    {
        if (PolyHasError || _polyBuilding) return;
        _polyPreviewAgain = false;
        PolyBuilding = true;
        var spec = PolySpecJson(Math.Min(_growDp, 60));
        var tact = _growTact;
        try
        {
            var scale = _polyArch == 0 ? 0.8 : -0.8;   // branched: step down where a branch point is crowded
            var seed = (ulong)_polyDraw;   // the chain the strip shows
            var (doc, _) = await Task.Run(() => CapsDocument.GrowChains(spec, new CapsGrowOpts { Chains = 1, Dp = 0, Tacticity = tact, Seed = seed, Density = 0.02, ContactScale = scale, Curve = 1 }, null, "chain"));
            var old = _polyDoc;
            PolyDoc = doc;
            old?.Dispose();
        }
        catch (Exception e) { PolyError = "Preview: " + e.Message; }
        finally { PolyBuilding = false; }
        if (_polyPreviewAgain && PolyAutoPreview) { _polyPreviewAgain = false; _ = BuildPolyPreview(); }
    }

    /// <summary>One chain of the spec opened as the Studio document.</summary>
    public async Task BuildPolymerInStudio()
    {
        if (PolyHasError || !Idle) return;
        var spec = PolySpecJson();
        var tact = _growTact;
        PolyBuilding = true;
        try
        {
            var scale = _polyArch == 0 ? 0.8 : -0.8;
            var seed = (ulong)_polyDraw;
            var (doc, _) = await Task.Run(() => CapsDocument.GrowChains(spec, new CapsGrowOpts { Chains = 1, Dp = 0, Tacticity = tact, Seed = seed, Density = 0.02, ContactScale = scale, Curve = 1 }, null, "chain"));
            Show(doc, (_polyName.Length > 0 ? _polyName : "polymer") + (_polyArch == 0 ? $" · 1 chain × {_growDp}" : $" · 1 {ArchIds[_polyArch]} molecule"));
            GrownUnsaved = true;
            SetModule(8);
        }
        catch (Exception e) { PolyError = e.Message; }
        finally { PolyBuilding = false; }
    }

    // ---- Grow uses the chain spec once one is sent
    private string? _growSpec;
    private string _growSpecName = "PS";
    private int _growSpecAtoms;
    private double _growSpecMass;
    public bool GrowHasSpec => _growSpec != null;
    public string GrowComponentName => _growSpec == null ? "PS " + GrowTacticityName : _growSpecName;

    public void SendPolymerToGrow()
    {
        if (PolyHasError) return;
        _growSpec = PolySpecJson();
        _growSpecName = PolyUnits.Count == 1 ? (PolyUnits[0].Name.Length > 0 ? PolyUnits[0].Name : "polymer") : (_polyName.Length > 0 ? _polyName : "copolymer");
        _growSpecAtoms = _polyAtoms;
        _growSpecMass = _polyMass;
        Raise(nameof(GrowHasSpec)); Raise(nameof(GrowComponentName)); Raise(nameof(GrowAtomsText)); Raise(nameof(GrowEstimate)); Raise(nameof(GrowCommand));
        SetModule(0);
        Status = $"Grow builds {_growSpecName} now · set chains, density and contact scale, then Grow";
    }

    public void UsePolystyreneInGrow()
    {
        _growSpec = null;
        _growSpecName = "PS";   // the cell's name follows the polymer grown
        Raise(nameof(GrowHasSpec)); Raise(nameof(GrowComponentName)); Raise(nameof(GrowAtomsText)); Raise(nameof(GrowEstimate)); Raise(nameof(GrowCommand));
    }

    /// <summary>Atoms and mass per chain for the chosen polymer at the current DP (polystyrene when none).</summary>
    private (long Atoms, double Mass) GrowChainSize()
    {
        if (_growSpec == null) return (16L * _growDp + 2, _growDp * (8 * 12.011 + 8 * 1.008) + 2 * 1.008);
        try
        {
            var spec = JsonNode.Parse(_growSpec)!.AsObject();
            spec["dp"] = _growDp;
            var r = JsonNode.Parse(CapsDocument.ChainPreview(spec.ToJsonString(), (ulong)_growSeed))!;
            return ((int?)r["atoms"] ?? 0, (double?)r["mass"] ?? 0);
        }
        catch { return (_growSpecAtoms, _growSpecMass); }
    }
}
