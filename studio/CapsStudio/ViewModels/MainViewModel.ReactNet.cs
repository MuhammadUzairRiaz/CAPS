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
            RxTargetChosen(_rxTargetKind);
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

    /// <summary>One link of the last run: each side's chain and repeat unit, the bridging molecule.</summary>
    public sealed record RxLink(int ChainA, int UnitA, int ChainB, int UnitB, int Via, string ViaName, int Cycle);
    public List<RxLink> RxLinks { get; } = new();
    public int RxChainCount { get; private set; }
    public int RxUnitsPerChain { get; private set; }

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
            RxLinks.Clear();
            foreach (var l in j["links"]?.AsArray() ?? [])
                if (l != null)
                    RxLinks.Add(new RxLink((int?)l["chain_a"] ?? 0, (int?)l["unit_a"] ?? 0, (int?)l["chain_b"] ?? 0, (int?)l["unit_b"] ?? 0, (int?)l["via"] ?? 0,
                                           (string?)l["via_name"] ?? "", (int?)l["cycle"] ?? 0));
            RxChainCount = chains;
            RxUnitsPerChain = chains > 0 && mono > 0 ? (mono + chains - 1) / chains : 0;
            foreach (var l in RxLinks.Take(12))
                lines.Add($"  chain {l.ChainA} unit {l.UnitA} — {(l.Via > 0 ? $"{l.ViaName} #{l.Via}" : "direct")} — chain {l.ChainB} unit {l.UnitB}");
            if (RxLinks.Count > 12) lines.Add($"  … {RxLinks.Count - 12} more (the link map shows all)");
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
    // the force field the reactions are typed with (the assigned one, or a library force field with automatic typing),
    // the LAMMPS styles (the force field's own, as the engine export writes them, or CAPS-exact), and the cure options
    public List<string> RxBrFfChoices => new[] { "The assigned force field" + (Field.Assigned ? $" ({Field.ForceFieldName})" : "") }
        .Concat(Field.Library.Where(e => e.AutoTyping).Select(e => e.Label)).ToList();
    private int _rxBrFf, _rxBrStyles, _rxBrKspace, _rxBrMolIds;
    public int RxBrFfIndex { get => _rxBrFf; set => Set(ref _rxBrFf, Math.Max(0, value)); }
    public static string[] RxBrStyleChoices => ["The force field's own (as the engine export)", "CAPS-exact forms"];
    public int RxBrStylesIndex { get => _rxBrStyles; set => Set(ref _rxBrStyles, Math.Clamp(value, 0, 1)); }
    public static string[] RxBrKspaceChoices => ["auto", "pppm", "ewald", "dsf", "cut"];
    public int RxBrKspaceIndex { get => _rxBrKspace; set => Set(ref _rxBrKspace, Math.Clamp(value, 0, 4)); }
    public static string[] RxBrMolIdChoices => ["Reset to the bonded pieces (LAMMPS default)", "Keep the data file's (reset_mol_ids no)", "Chains keep their ids (molmap)"];
    public int RxBrMolIdsIndex { get => _rxBrMolIds; set => Set(ref _rxBrMolIds, Math.Clamp(value, 0, 2)); }
    private string _rxBrSurvey = "", _rxBrGroups = "", _rxBrTargets = "", _rxBrLinks = "";
    public string RxBrSurvey { get => _rxBrSurvey; set => Set(ref _rxBrSurvey, value ?? ""); }
    public string RxBrGroups { get => _rxBrGroups; set => Set(ref _rxBrGroups, value ?? ""); }
    public string RxBrTargets { get => _rxBrTargets; set => Set(ref _rxBrTargets, value ?? ""); }
    public string RxBrLinks { get => _rxBrLinks; set => Set(ref _rxBrLinks, value ?? ""); }
    private double _rxBrLimiting, _rxBrCheck = 1000, _rxBrMax = 2000000;
    public decimal RxBrLimitingD { get => (decimal)_rxBrLimiting; set => Set(ref _rxBrLimiting, (double)Math.Clamp(value, 0m, 1e9m), nameof(RxBrLimitingD)); }
    public decimal RxBrCheckD { get => (decimal)_rxBrCheck; set => Set(ref _rxBrCheck, (double)Math.Clamp(value, 1m, 1e9m), nameof(RxBrCheckD)); }
    public decimal RxBrMaxStepsD { get => (decimal)_rxBrMax; set => Set(ref _rxBrMax, (double)Math.Clamp(value, 1m, 1e12m), nameof(RxBrMaxStepsD)); }
    private double _rxBrStabilize = 200, _rxBrHTransfer = 3.5;
    public decimal RxBrStabilizeD { get => (decimal)_rxBrStabilize; set => Set(ref _rxBrStabilize, (double)Math.Clamp(value, 1m, 100000m), nameof(RxBrStabilizeD)); }
    public decimal RxBrHTransferD { get => (decimal)_rxBrHTransfer; set => Set(ref _rxBrHTransfer, (double)Math.Clamp(value, 0m, 20m), nameof(RxBrHTransferD)); }

    private static JsonArray NumberList(string text)
    {
        var a = new JsonArray();
        foreach (var w in text.Split(new[] { ',', ';', ' ' }, StringSplitOptions.RemoveEmptyEntries))
            a.Add(double.Parse(w.TrimEnd('%'), NumberStyles.Float, CultureInfo.InvariantCulture));
        return a;
    }
    /// <summary>The options JSON of the fix bond/react export from the page.</summary>
    public string BondReactOptions()
    {
        var o = new JsonObject
        {
            ["stem"] = _rxBrStem, ["variants"] = _rxVariants, ["keep_byproducts"] = _rxKeepBy, ["between_chains"] = _rxBetween,
            ["nevery"] = 100, ["temperature"] = _rxTemp, ["steps"] = _rxBrSteps, ["seed"] = Math.Max(1, _rxSeed),
            ["lammps_styles"] = _rxBrStyles == 0 ? "native" : "exact", ["kspace"] = RxBrKspaceChoices[_rxBrKspace],
            ["mol_ids"] = _rxBrMolIds switch { 1 => "keep", 2 => "molmap", _ => "reset" },
            ["stabilize_steps"] = (int)_rxBrStabilize, ["h_transfer_max"] = _rxBrHTransfer,
        };
        if (RxSeveral && _rxByWeights) o["weights"] = new JsonArray(RxReactions.Select(r => (JsonNode)(double)r.WeightD).ToArray());
        var choices = Field.Library.Where(e => e.AutoTyping).ToList();
        if (_rxBrFf > 0 && _rxBrFf - 1 < choices.Count) o["forcefield"] = choices[_rxBrFf - 1].File;
        if (_rxBrSurvey.Trim().Length > 0) o["survey_after"] = NumberList(_rxBrSurvey);
        // "PBS: 7-10; ENR50: 1-6; MAH: 11-13" → type groups by molecule
        if (_rxBrGroups.Trim().Length > 0)
        {
            var g = new JsonArray();
            foreach (var part in _rxBrGroups.Split(';', StringSplitOptions.RemoveEmptyEntries))
            {
                var kv = part.Split(':', 2);
                if (kv.Length != 2) throw new FormatException("type groups: name: molecules; … (PBS: 7-10; ENR50: 1-6)");
                g.Add(new JsonObject { ["name"] = kv[0].Trim(), ["molecules"] = kv[1].Trim() });
            }
            o["type_groups"] = g;
        }
        if (_rxBrTargets.Trim().Length > 0)
        {
            o["targets"] = NumberList(_rxBrTargets);
            o["limiting"] = _rxBrLimiting;
            o["check_every"] = _rxBrCheck;
            o["max_steps"] = _rxBrMax;
            if (_rxBrLinks.Trim().Length > 0)
                o["link_reactions"] = new JsonArray(_rxBrLinks.Split(new[] { ',', ' ' }, StringSplitOptions.RemoveEmptyEntries).Select(x => (JsonNode)x).ToArray());
        }
        return o.ToJsonString();
    }

    public async Task<bool> ExportBondReact(string dir)
    {
        if (_doc == null || !Idle) return false;
        var doc = _doc;
        try
        {
            var opts = BondReactOptions();
            Status = "Writing the fix bond/react set…";
            var text = _rxText;
            var rep = await Task.Run(() => doc.BondReactExport(text, dir, opts));
            var j = JsonNode.Parse(rep);
            var files = j?["files"]?.AsArray().Count ?? 0;
            var notes = string.Join("\n", j?["notes"]?.AsArray().Select(n => (string?)n ?? "") ?? []);
            var variants = string.Join("\n", j?["variants"]?.AsArray().Select(v => $"{(string?)v?["name"]}: {(int?)v?["sites"]} pairs · {(int?)v?["pre_atoms"]} atoms, {(int?)v?["edge"]} edge, {(int?)v?["deleted"]} deleted") ?? []);
            var steps = string.Join("\n", j?["steps"]?.AsArray().Select(v => $"{(string?)v?["reaction"]}{((string?)v?["step"] is { Length: > 0 } st ? " (" + st + " step)" : "")}: " +
                                                                              $"{(int?)v?["covered"]} of {(int?)v?["candidates"]} sites covered by {(int?)v?["variants"]} templates") ?? []);
            var links = j?["link_reactions"]?.AsArray() is { Count: > 0 } la ? "\nLinks counted: " + string.Join(", ", la.Select(x => (string?)x)) : "";
            RxBrText = $"{files} files in {dir}\n{notes}\n{steps}\n{variants}{links}\nRun: lmp -in {_rxBrStem}.in (in that folder)";
            Status = $"Wrote the fix bond/react set to {dir}";
            return true;
        }
        catch (Exception e) { RxBrText = "Could not write: " + e.Message; Status = "Could not write the fix bond/react set"; return false; }
    }

    // ---- the set run in LAMMPS (fix bond/react) from the page: the progress read from its log and crosslink_progress.dat
    private System.Diagnostics.Process? _rxLmp;
    private bool _rxLmpRunning;
    public bool RxLmpRunning { get => _rxLmpRunning; private set { if (Set(ref _rxLmpRunning, value)) Raise(nameof(RxLmpIdle)); } }
    public bool RxLmpIdle => !_rxLmpRunning;
    private string _rxLmpText = "";
    public string RxLmpText { get => _rxLmpText; private set => Set(ref _rxLmpText, value); }
    public async Task RunBondReact(string dir)
    {
        if (_rxLmpRunning) return;
        var lmp = FindLammps();
        if (lmp == null) { RxLmpText = "LAMMPS was not found (lmp on the PATH, or set LAMMPS_EXE): the files are written, run them where LAMMPS is."; }
        if (!await ExportBondReact(dir) || lmp == null) return;
        var stem = _rxBrStem;
        try
        {
            // the packages the deck needs: REACTION always, CLASS2 for class II styles (PCFF, COMPASS), KSPACE for PPPM / Ewald
            var deck = await System.IO.File.ReadAllTextAsync(System.IO.Path.Combine(dir, stem + ".in"));
            var need = new List<string> { "REACTION", "MOLECULE" };
            if (deck.Contains("class2")) need.Add("CLASS2");
            if (deck.Contains("kspace_style")) need.Add("KSPACE");
            var help = new System.Diagnostics.ProcessStartInfo(lmp, "-h") { UseShellExecute = false, RedirectStandardOutput = true, CreateNoWindow = true };
            using (var h = System.Diagnostics.Process.Start(help))
            {
                var text = h == null ? "" : await h.StandardOutput.ReadToEndAsync();
                if (h != null) await h.WaitForExitAsync();
                var at = text.IndexOf("Installed packages", StringComparison.Ordinal);
                var end = at < 0 ? -1 : text.IndexOf("List of individual", at, StringComparison.Ordinal);
                var pk = at < 0 ? "" : text[at..(end < 0 ? text.Length : end)];
                var words = pk.Split((char[]?)null, StringSplitOptions.RemoveEmptyEntries).ToHashSet();
                var missing = need.Where(n => !words.Contains(n)).ToList();
                if (at >= 0 && missing.Count > 0)
                {
                    RxLmpText = $"{lmp} was built without {string.Join(", ", missing)}: point LAMMPS_EXE at a LAMMPS built with them " +
                                "(cmake -D PKG_" + string.Join("=on -D PKG_", missing) + "=on). The files are written; nothing was run.";
                    return;
                }
            }
            var psi = new System.Diagnostics.ProcessStartInfo(lmp) { WorkingDirectory = dir, UseShellExecute = false, RedirectStandardOutput = true, RedirectStandardError = true, CreateNoWindow = true };
            foreach (var a in new[] { "-in", stem + ".in", "-log", "log.lammps" }) psi.ArgumentList.Add(a);
            var p = System.Diagnostics.Process.Start(psi) ?? throw new InvalidOperationException("LAMMPS did not start");
            _rxLmp = p;
            RxLmpRunning = true;
            var lastLine = "";
            p.OutputDataReceived += (_, e) => { if (e.Data is { Length: > 0 } l) lastLine = l; };
            p.ErrorDataReceived += (_, e) => { if (e.Data is { Length: > 0 } l) lastLine = l; };
            p.BeginOutputReadLine();
            p.BeginErrorReadLine();
            Status = $"LAMMPS is running {stem}.in in {dir}";
            var progress = System.IO.Path.Combine(dir, "crosslink_progress.dat");
            while (!p.HasExited)
            {
                await Task.Delay(2000);
                var row = "";
                try
                {
                    if (System.IO.File.Exists(progress))
                        row = System.IO.File.ReadLines(progress).LastOrDefault(l => !l.StartsWith('#')) ?? "";
                }
                catch (System.IO.IOException) { }
                RxLmpText = "LAMMPS running · " + (row.Length > 0 ? "step, time (ps), reactions…, links, crosslink density (%), ν, bonds: " + row : lastLine);
            }
            var ok = p.ExitCode == 0;
            var xl = System.IO.Directory.GetFiles(dir, stem + "_XL*.data").Select(System.IO.Path.GetFileName).ToList();
            RxLmpText = (ok ? "LAMMPS finished" : $"LAMMPS stopped (exit {p.ExitCode}): {lastLine}") +
                        (xl.Count > 0 ? "\nWritten at the targets: " + string.Join(", ", xl) : "") + $"\nLog: {System.IO.Path.Combine(dir, "log.lammps")}";
            Status = ok ? "LAMMPS fix bond/react run finished" : "LAMMPS stopped";
        }
        catch (Exception e) { RxLmpText = "Could not run LAMMPS: " + e.Message; }
        finally { RxLmpRunning = false; _rxLmp = null; }
    }
    public void StopBondReact()
    {
        try { if (_rxLmp is { HasExited: false } p) p.Kill(true); } catch (Exception) { }
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

/// <summary>One reaction of the library: its scheme and roles, and the CAPS template it becomes (or why not).</summary>
public sealed record RxLibEntry(string Id, string Name, string Category, string Description, string Scheme, string Roles, string Template, string Error,
                                string[] Steps)
{
    public bool Runs => Template.Length > 0;
    public string Status => Runs ? "template ready" : Steps.Length > 0 ? "in steps: " + string.Join(" → ", Steps) : Error;
}

public sealed partial class MainViewModel
{
    // ---- the reaction library: CAPS's cure templates and the worked schemes (data/reactions)
    private List<RxLibEntry>? _rxLib;
    private int _rxLibCat;
    private RxLibEntry? _rxLibSel;
    private bool _rxLibOpen;
    public bool RxLibOpen { get => _rxLibOpen; set { if (Set(ref _rxLibOpen, value) && value) LoadRxLib(); } }
    private readonly ObservableCollection<string> _rxLibCats = new();
    /// <summary>The reaction choice's categories (CAPS's cures first, then the library's classes); loaded on first use.</summary>
    public ObservableCollection<string> RxLibCategories { get { LoadRxLib(); return _rxLibCats; } }
    public ObservableCollection<RxLibEntry> RxLibItems { get; } = new();
    public int RxLibCategory { get => _rxLibCat; set { if (Set(ref _rxLibCat, value)) FilterRxLib(); } }
    private bool _rxLibFiltering;
    private string _rxTextBeforePick = "";   // the reactions before the last pick: Add puts the picked one beside them
    /// <summary>The chosen reaction: picking one loads it in place of the text (Add keeps the others).</summary>
    public RxLibEntry? RxLibSelected
    {
        get => _rxLibSel;
        set
        {
            if (!Set(ref _rxLibSel, value)) return;
            Raise(nameof(RxLibHasSelection)); Raise(nameof(RxChoiceNote)); Raise(nameof(RxChoiceHasScheme));
            if (value != null && !_rxLibFiltering) { _rxTextBeforePick = _rxText; UseRxLib(false); }
        }
    }
    public bool RxLibHasSelection => _rxLibSel != null;
    /// <summary>What the chosen reaction does (a cure's note, or a library reaction's description).</summary>
    public string RxChoiceNote => _rxLibSel?.Description is { Length: > 0 } d ? d : RxSetNote;
    public bool RxChoiceHasScheme => _rxLibSel?.Scheme is { Length: > 0 };
    private const string CuresCategory = "Cures & crosslinks";

    private void LoadRxLib()
    {
        if (_rxLib != null) return;
        _rxLib = new List<RxLibEntry>();
        for (var k = 0; k < ReactionSets.Length - 1; ++k)   // the cures (all but Custom), by their short names
            _rxLib.Add(new RxLibEntry($"set:{k}", ReactionSets[k], CuresCategory, RxSetNotes[k], "", "", "set", "", []));
        try
        {
            if (Paths.Reactions is { } path)
            {
                var j = JsonNode.Parse(CapsDocument.ReactionLibrary(path));
                foreach (var e in j?["reactions"]?.AsArray() ?? [])
                {
                    if (e == null) continue;
                    string Side(string k) => string.Join(" + ", e[k]?.AsArray().Select(x => $"{(string?)x?["smiles"]} ({(string?)x?["label"]})") ?? []);
                    var roles = string.Join("\n", e["tags"]?.AsArray().Select(t => $"{(string?)t?["map"]}  {(string?)t?["role"]}") ?? []);
                    var note = (string?)e["note"] ?? "";
                    _rxLib.Add(new RxLibEntry((string?)e["id"] ?? "", (string?)e["name"] ?? "", (string?)e["category"] ?? "", (string?)e["description"] ?? "",
                        Side("reactants") + "\n  →  " + Side("products"), roles + (note.Length > 0 ? "\n" + note : ""), (string?)e["template"] ?? "",
                        (string?)e["error"] ?? "", e["steps"]?.AsArray().Select(x => (string?)x ?? "").ToArray() ?? []));
                }
            }
        }
        catch (Exception e) { RxLog = "Could not read the reaction library: " + e.Message; }
        _rxLibCats.Clear();
        foreach (var c in _rxLib.Select(x => x.Category).Distinct()) _rxLibCats.Add(c);
        _rxLibCat = 0;
        Raise(nameof(RxLibCategory));
        FilterRxLib();
    }
    /// <summary>The pickers show the cure in use (chosen elsewhere: a guide, a script, a session).</summary>
    private void SyncRxChoice()
    {
        if (_rxLib == null) return;
        var id = $"set:{_rxSet}";
        if (_rxLib.FirstOrDefault(x => x.Id == id) is not { } e) { Raise(nameof(RxChoiceNote)); return; }
        var cat = _rxLibCats.IndexOf(e.Category);
        if (cat >= 0 && cat != _rxLibCat) { _rxLibCat = cat; Raise(nameof(RxLibCategory)); FilterRxLib(); }
        _rxLibFiltering = true;
        RxLibSelected = RxLibItems.FirstOrDefault(x => x.Id == id);
        _rxLibFiltering = false;
        Raise(nameof(RxChoiceNote));
    }

    private void FilterRxLib()
    {
        _rxLibFiltering = true;
        RxLibItems.Clear();
        if (_rxLib != null && _rxLibCat >= 0 && _rxLibCat < _rxLibCats.Count)
            foreach (var e in _rxLib.Where(x => x.Category == _rxLibCats[_rxLibCat])) RxLibItems.Add(e);
        // the reaction in use stays shown when it is in this category; switching category loads nothing by itself
        RxLibSelected = RxLibItems.FirstOrDefault(x => x.Id == $"set:{_rxSet}") ;
        _rxLibFiltering = false;
    }
    /// <summary>The selected reaction's template(s): in place of the text, or added to it; a scheme of three molecules adds its steps.</summary>
    public void UseRxLib(bool add)
    {
        if (_rxLibSel is not { } e || _rxLib == null) return;
        if (e.Id.StartsWith("set:", StringComparison.Ordinal) && int.TryParse(e.Id[4..], out var set))
        {
            if (!add) { RxSet = set; RxLog = e.Name + ": " + e.Description; return; }
            var keep = _rxTextBeforePick.Trim().Length > 0 ? _rxTextBeforePick : _rxText;
            var old = _rxSet;
            _rxSet = set;
            LoadReactionSet();   // the set's text …
            var more = _rxText;
            _rxSet = ReactionSets.Length - 1;
            Raise(nameof(RxSet)); Raise(nameof(RxSetNote));
            RxText = (keep.TrimEnd() + "\n\n" + more).Trim() + "\n";   // … added to what was there
            RxLog = e.Name + " added: the reactions run together";
            _ = old;
            return;
        }
        var texts = new List<string>();
        if (e.Runs) texts.Add(e.Template);
        else foreach (var id in e.Steps.Distinct()) if (_rxLib.FirstOrDefault(x => x.Id == id && x.Runs) is { } st) texts.Add(st.Template);
        if (texts.Count == 0) { RxLog = e.Name + ": " + e.Error; return; }
        _rxSet = ReactionSets.Length - 1;
        Raise(nameof(RxSet)); Raise(nameof(RxSetNote));
        if (add && _rxTextBeforePick.Trim().Length > 0) RxText = _rxTextBeforePick;   // beside the reactions there were before the pick
        var have = new HashSet<string>(RxReactions.Select(r => r.Name));
        var fresh = add ? texts.Where(t => !have.Contains(System.Text.RegularExpressions.Regex.Match(t, @"reaction\s+(\S+)").Groups[1].Value)).ToList() : texts;
        RxText = add ? (RxText.TrimEnd() + "\n\n" + string.Join("\n", fresh)).Trim() + "\n" : string.Join("\n", fresh);
        RxLog = $"{e.Name}: {(e.Runs ? "template" : "steps " + string.Join(" → ", e.Steps))} {(add ? "added" : "loaded")} · " + e.Description;
    }
}
