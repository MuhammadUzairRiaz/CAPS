using System;
using System.Collections.Generic;
using System.Collections.ObjectModel;
using System.Globalization;
using System.Linq;
using System.Text.Json;
using System.Text.Json.Nodes;
using System.Threading.Tasks;
using CapsStudio.Interop;

namespace CapsStudio.ViewModels;

/// <summary>One reaction of the React text, with its relative rate when several run together.</summary>
public sealed class RxReaction : ObservableObject
{
    public RxReaction(string name, decimal weight) { Name = name; _weight = weight; }
    public string Name { get; }
    private decimal _weight;
    public decimal WeightD { get => _weight; set => Set(ref _weight, Math.Clamp(value, 0.01m, 1000m)); }
}

/// <summary>React › the network (what links, how many, with which force field), live progress, sulfur donors and the
/// LAMMPS fix bond/react set both ways.</summary>
public sealed partial class MainViewModel
{
    // ---- reactions in the text, and more added from the built-ins
    public static readonly (string Name, string Label)[] RxBuiltins =
    [
        ("cc_crosslink", "C–C crosslink (H₂ leaves)"), ("sulfur_allylic", "Sulfur donor S–H + allylic C–H → C–S"),
        ("peroxide_allylic", "Peroxide: allylic C–C"), ("polysulfide_allylic", "Silane polysulfide → C–S"),
        ("enr_acid_ester", "ENR epoxide + COOH → β-hydroxy ester"), ("anhydride_alcohol", "Anhydride (MAH) + OH → half-ester acid"),
        ("ester_condensation", "COOH + OH → ester + H₂O"), ("epoxy_amine_primary", "Epoxide + primary amine"),
        ("epoxy_amine_secondary", "Epoxide + secondary amine"),
    ];
    public static string[] RxBuiltinLabels => RxBuiltins.Select(b => b.Label).ToArray();
    private int _rxAddIndex;
    public int RxAddIndex { get => _rxAddIndex; set => Set(ref _rxAddIndex, Math.Clamp(value, 0, RxBuiltins.Length - 1)); }
    /// <summary>Another built-in reaction added to the text: several reactions run together (ENR–MAH–ENR with ENR–PBS …).</summary>
    public void AddReaction()
    {
        try
        {
            var add = CapsDocument.ReactionTemplate(RxBuiltins[_rxAddIndex].Name);
            if (RxReactions.Any(r => r.Name == RxBuiltins[_rxAddIndex].Name)) { RxLog = RxBuiltins[_rxAddIndex].Name + " is already in the list"; return; }
            _rxSet = ReactionSets.Length - 1;   // custom from here on
            Raise(nameof(RxSet));
            RxText = RxText.TrimEnd() + "\n\n" + add;
        }
        catch (Exception e) { RxLog = e.Message; }
    }
    public ObservableCollection<RxReaction> RxReactions { get; } = new();
    public bool RxSeveral => RxReactions.Count > 1;
    private bool _rxByWeights;
    /// <summary>Several reactions: false — the closest pairs first whatever the reaction (automatic, by distance); true —
    /// each pick draws a reaction in proportion to its weight, then its closest free pair.</summary>
    public bool RxByWeights { get => _rxByWeights; set => Set(ref _rxByWeights, value); }
    public bool RxByDistance { get => !_rxByWeights; set { if (value) RxByWeights = false; } }
    private void RefreshRxReactions()
    {
        var names = System.Text.RegularExpressions.Regex.Matches(_rxText, @"(?m)^\s*reaction\s+(\S+)").Select(m => m.Groups[1].Value).ToList();
        if (names.SequenceEqual(RxReactions.Select(r => r.Name))) return;
        var old = RxReactions.ToDictionary(r => r.Name, r => r.WeightD);
        RxReactions.Clear();
        foreach (var n in names) RxReactions.Add(new RxReaction(n, old.TryGetValue(n, out var w) ? w : 1m));
        Raise(nameof(RxSeveral));
    }

    // ---- the network
    private bool _rxBetween = true, _rxKeepBy, _rxUseField = true, _rxAutoCapture = true;
    private double _rxCaptureMax = 8, _rxTargetValue = 1;
    private int _rxTargetKind;
    public bool RxBetweenChains { get => _rxBetween; set => Set(ref _rxBetween, value); }
    public bool RxKeepByproducts { get => _rxKeepBy; set => Set(ref _rxKeepBy, value); }
    public bool RxUseField { get => _rxUseField; set { if (Set(ref _rxUseField, value)) Raise(nameof(RxFieldText)); } }
    public bool RxAutoCapture { get => _rxAutoCapture; set => Set(ref _rxAutoCapture, value); }
    public decimal RxCaptureMaxD { get => (decimal)_rxCaptureMax; set => Set(ref _rxCaptureMax, (double)Math.Clamp(value, 1m, 20m), nameof(RxCaptureMaxD)); }
    public static readonly string[] RxTargetKinds =
        ["Conversion of the reactive sites", "Links between chains (number)", "Links per chain", "Crosslink density ν (mol/m³)",
         "Mc — mass between crosslinks (g/mol)", "Degree of crosslinking DC (%)"];
    public int RxTargetKind
    {
        get => _rxTargetKind;
        set
        {
            if (!Set(ref _rxTargetKind, Math.Clamp(value, 0, RxTargetKinds.Length - 1))) return;
            _rxTargetValue = _rxTargetKind switch { 1 => 10, 2 => 1, 3 => 100, 4 => 5000, 5 => 5, _ => 1 };
            Raise(nameof(RxTargetValueD)); Raise(nameof(RxTargetIsConversion)); Raise(nameof(RxTargetHelp));
        }
    }
    public bool RxTargetIsConversion => _rxTargetKind == 0;
    public decimal RxTargetValueD { get => (decimal)_rxTargetValue; set => Set(ref _rxTargetValue, (double)Math.Max(0.0001m, value), nameof(RxTargetValueD)); }
    public string RxTargetHelp => _rxTargetKind switch
    {
        0 => "Stops when this fraction of the counted reactive groups has reacted (Target conversion below).",
        1 => "A link joins two different chains (a curative bridging two chains is one link).",
        2 => "Each link joins two chains: links per chain = 2 × links / chains.",
        3 => "ν = links / (V N_A), the cell's volume.",
        4 => "Mc = chain mass / (2 × links): strands between tetrafunctional junctions, chain ends not subtracted.",
        _ => "DC = 2 × links / monomers × 100 % (Vasilev et al. 2021; Alamfard et al. 2023); monomers are the chains' repeat units.",
    };
    public string RxFieldText => !_rxUseField ? "Built-in default (GAFF for C and H, UFF otherwise) during the run; assign a force field afterwards"
        : Field.Assigned ? $"{Field.ForceFieldName}: types, charges and parameters after every cycle, and assigned to the network at the end"
        : "No force field assigned: the built-in default runs (assign one in the Force field step to crosslink with it)";

    // ---- live progress
    private string _rxLiveText = "";
    public string RxLiveText { get => _rxLiveText; private set => Set(ref _rxLiveText, value); }
    private string _rxNetworkText = "";
    public string RxNetworkText { get => _rxNetworkText; private set => Set(ref _rxNetworkText, value); }
    private Action<CapsDocument, string> ReactLiveReceiver()
    {
        var ticket = ++_runLiveTicket;
        var inv = CultureInfo.InvariantCulture;
        return (live, stats) => Avalonia.Threading.Dispatcher.UIThread.Post(() =>
        {
            if (ticket != _runLiveTicket) { live.Dispose(); return; }
            var old = RunLiveDoc;
            RunLiveDoc = live;
            old?.Dispose();
            try
            {
                var j = JsonNode.Parse(stats);
                int links = (int?)j?["crosslinks"] ?? 0, target = (int?)j?["target"] ?? 0, cycle = (int?)j?["cycle"] ?? 0;
                double nu = (double?)j?["density"] ?? 0, dc = (double?)j?["degree"] ?? 0, conv = (double?)j?["conversion"] ?? 0;
                var t = target > 0 ? string.Format(inv, "{0} / {1} links ({2:F0} %)", links, target, 100.0 * links / target) : $"{links} links";
                RunLiveText = string.Format(inv, "cycle {0} · {1} · ν {2:G4} mol/m³", cycle, t, nu) + (dc > 0 ? string.Format(inv, " · DC {0:F2} %", dc) : "")
                              + string.Format(inv, " · α {0:F3}", conv);
                RxLiveText = RunLiveText;
            }
            catch { RunLiveText = "reacting"; }
        });
    }

    /// <summary>The network of the last run, from the core's summary.</summary>
    private void LoadReactSummary(CapsDocument doc)
    {
        try
        {
            var j = JsonNode.Parse(doc.ReactSummary());
            if (j == null) return;
            var inv = CultureInfo.InvariantCulture;
            int chains = (int?)j["chains"] ?? 0, links = (int?)j["crosslinks"] ?? 0, target = (int?)j["target"] ?? 0, mono = (int?)j["monomers"] ?? 0,
                by = (int?)j["byproducts"] ?? 0, loops = (int?)j["intrachain"] ?? 0;
            var lines = new List<string>
            {
                target > 0 ? string.Format(inv, "Links between chains: {0} of the {1} asked ({2:F0} %)", links, target, 100.0 * links / target)
                           : string.Format(inv, "Links between chains: {0}", links),
                string.Format(inv, "Chains {0} · per chain {1:F2} · ν {2:G4} mol/m³ · Mc ≈ {3:G4} g/mol", chains, (double?)j["per_chain"] ?? 0, (double?)j["density"] ?? 0, (double?)j["mc"] ?? 0),
            };
            if (mono > 0) lines.Add(string.Format(inv, "Degree of crosslinking DC {0:F2} % ({1} repeat units)", (double?)j["degree"] ?? 0, mono));
            if (loops > 0) lines.Add($"{loops} links closed within one chain (loops, not counted)");
            if (by > 0) lines.Add($"{by} byproduct molecules {(_rxKeepBy ? "kept in the cell" : "removed")}");
            lines.Add("Force field during the run: " + ((string?)j["field"] ?? "—"));
            lines.Add("After: " + ((string?)j["field_after"] ?? "—"));
            RxNetworkText = string.Join("\n", lines);
        }
        catch { }
    }

    // ---- sulfur donors: H–Sx–H with x fixed or spread over a range, dosed in phr of the rubber
    private int _rxSxMin = 2, _rxSxMax = 2;
    private double _rxPhr = 2.5;
    public decimal RxSxMinD { get => _rxSxMin; set { _rxSxMin = Math.Clamp((int)value, 1, 8); if (_rxSxMax < _rxSxMin) { _rxSxMax = _rxSxMin; Raise(nameof(RxSxMaxD)); } Raise(); Raise(nameof(RxSulfurText)); } }
    public decimal RxSxMaxD { get => _rxSxMax; set { _rxSxMax = Math.Clamp((int)value, 1, 8); if (_rxSxMin > _rxSxMax) { _rxSxMin = _rxSxMax; Raise(nameof(RxSxMinD)); } Raise(); Raise(nameof(RxSulfurText)); } }
    public decimal RxPhrD { get => (decimal)_rxPhr; set { _rxPhr = (double)Math.Clamp(value, 0.01m, 100m); Raise(); Raise(nameof(RxSulfurText)); } }
    /// <summary>The donors a dose makes: sulfur atoms = phr / 100 × the cell's mass / 32.06, in H–Sx–H of x drawn evenly
    /// from the range (one x when the range is one number).</summary>
    public string RxSulfurText
    {
        get
        {
            if (_doc == null) return "";
            var plan = SulfurPlan();
            var s = plan.Sum(p => p.X * p.N);
            return plan.Count == 0 ? "No donor at this dose: raise the phr or build a larger cell"
                : $"{s} S atoms in {plan.Sum(p => p.N)} donors: " + string.Join(", ", plan.Select(p => $"{p.N} × S{p.X}"));
        }
    }
    private List<(int X, int N)> SulfurPlan()
    {
        if (_doc == null) return [];
        var mass = _doc.Summary().TotalMass;
        var atoms = (int)Math.Round(_rxPhr / 100.0 * mass / 32.06);
        var rng = new Random(_rxSeed);
        var count = new SortedDictionary<int, int>();
        for (var left = atoms; left > 0;)
        {
            var x = Math.Min(left, _rxSxMin + rng.Next(_rxSxMax - _rxSxMin + 1));
            if (x < _rxSxMin && count.Count > 0) break;   // the remainder is shorter than the shortest donor asked for
            count[x] = count.GetValueOrDefault(x) + 1;
            left -= x;
        }
        return count.Select(kv => (kv.Key, kv.Value)).ToList();
    }
    public async Task InsertSulfurDonors()
    {
        if (_doc == null || !Idle) return;
        var doc = _doc;
        var plan = SulfurPlan();
        if (plan.Count == 0) { RxLog = RxSulfurText; return; }
        try
        {
            Status = "Inserting sulfur donors…";
            var seed = (ulong)_rxSeed;
            var reps = new List<string>();
            foreach (var (x, n) in plan)
                reps.Add(await Task.Run(() => doc.InsertMolecules(new string('S', x), n, 2.0, seed++)));
            Field.Reset();
            AfterRun(doc, " · +S");
            RxLog = $"Inserted {RxSulfurText}\n" + string.Join("\n", reps) +
                    "\nThe sulfur cure template bonds each S–H end to an allylic carbon: a donor reaching two chains becomes a C–Sx–C bridge. Assign the force field again before crosslinking with it.";
            Raise(nameof(RxFieldText));
        }
        catch (Exception e) { RxLog = "Could not insert: " + e.Message; Status = "Could not insert the sulfur donors"; }
    }

    // ---- LAMMPS fix bond/react
    private int _rxVariants = 6;
    private string _rxBrStem = "react";
    private double _rxBrSteps = 100000;
    public decimal RxVariantsD { get => _rxVariants; set => Set(ref _rxVariants, Math.Clamp((int)value, 1, 50), nameof(RxVariantsD)); }
    public string RxBrStem { get => _rxBrStem; set => Set(ref _rxBrStem, string.IsNullOrWhiteSpace(value) ? "react" : value.Trim()); }
    public decimal RxBrStepsD { get => (decimal)_rxBrSteps; set => Set(ref _rxBrSteps, (double)Math.Clamp(value, 1m, 1e10m), nameof(RxBrStepsD)); }
    private string _rxBrText = "Writes pre- and post-reaction templates and map files cut from reaction sites of this structure, typed with its force field, " +
                               "a data file holding every type the reactions create, and an input with fix bond/react.";
    public string RxBrText { get => _rxBrText; private set => Set(ref _rxBrText, value); }
    public async Task ExportBondReact(string dir)
    {
        if (_doc == null || !Idle) return;
        var doc = _doc;
        var inv = CultureInfo.InvariantCulture;
        var weights = RxSeveral && _rxByWeights ? "[" + string.Join(",", RxReactions.Select(r => r.WeightD.ToString(inv))) + "]" : "[]";
        var opts = string.Format(inv, "{{\"stem\":{0},\"variants\":{1},\"keep_byproducts\":{2},\"between_chains\":{3},\"weights\":{4},\"nevery\":100,\"temperature\":{5},\"steps\":{6},\"seed\":{7}}}",
            JsonSerializer.Serialize(_rxBrStem), _rxVariants, _rxKeepBy ? "true" : "false", _rxBetween ? "true" : "false", weights, _rxTemp, _rxBrSteps, Math.Max(1, _rxSeed));
        try
        {
            Status = "Writing the fix bond/react set…";
            var text = _rxText;
            var rep = await Task.Run(() => doc.BondReactExport(text, dir, opts));
            var j = JsonNode.Parse(rep);
            var files = j?["files"]?.AsArray().Count ?? 0;
            var notes = string.Join("\n", j?["notes"]?.AsArray().Select(n => (string?)n ?? "") ?? []);
            var variants = string.Join("\n", j?["variants"]?.AsArray().Select(v => $"{(string?)v?["name"]}: {(int?)v?["sites"]} pairs · {(int?)v?["pre_atoms"]} atoms, {(int?)v?["edge"]} edge, {(int?)v?["deleted"]} deleted") ?? []);
            RxBrText = $"{files} files in {dir}\n{notes}\n{variants}\nRun: lmp -in {_rxBrStem}.in";
            Status = $"Wrote the fix bond/react set to {dir}";
        }
        catch (Exception e) { RxBrText = "Could not write: " + e.Message; Status = "Could not write the fix bond/react set"; }
    }
    /// <summary>A pre / post / map set (and the data file whose Masses name the types) read back into the React text.</summary>
    public void ImportBondReact(string pre, string post, string map, string masses)
    {
        try
        {
            var name = System.IO.Path.GetFileNameWithoutExtension(map).Replace("_map", "");
            var j = JsonNode.Parse(CapsDocument.BondReactImport(pre, post, map, masses, name, 0));
            var text = (string?)j?["text"] ?? "";
            _rxSet = ReactionSets.Length - 1;
            Raise(nameof(RxSet));
            RxText = (RxText.Trim().Length == 0 ? "" : RxText.TrimEnd() + "\n\n") + text;
            RxLog = "Read " + System.IO.Path.GetFileName(map) + ":\n" + string.Join("\n", j?["notes"]?.AsArray().Select(n => (string?)n ?? "") ?? []);
        }
        catch (Exception e) { RxLog = "Could not read the templates: " + e.Message; }
    }
}
