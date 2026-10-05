using System.Collections.ObjectModel;
using System.Globalization;
using System.Text.Json.Nodes;

namespace CapsStudio.ViewModels;

/// <summary>One pair of elements in Bond rules: its distances, its cut-off (or never bonded).</summary>
public sealed class RulePair : ObservableObject
{
    public string A { get; init; } = "";
    public string B { get; init; } = "";
    /// <summary>As chemists write it: C first, H last, the others alphabetically (C – H, B – N, O – Zn).</summary>
    public string Name
    {
        get
        {
            int Rank(string e) => e == "C" ? 0 : e == "H" ? 2 : 1;
            var (x, y) = Rank(A) < Rank(B) || (Rank(A) == Rank(B) && string.CompareOrdinal(A, B) <= 0) ? (A, B) : (B, A);
            return $"{x} – {y}";
        }
    }
    public int[] Counts { get; init; } = [];
    public double Lo { get; init; } = 0.8;
    public double Bin { get; init; } = 0.04;
    public double Hi => Lo + Counts.Length * Bin;
    public int Bonded { get; init; }
    public double Suggested { get; init; }
    public bool Ionic { get; init; }
    public event Action? Changed;
    /// <summary>The cut-off was let go: the counts and the view's preview follow.</summary>
    public Action? Committed;
    private double _cut;
    private bool _never;
    public double Cutoff { get => _cut; set { if (Set(ref _cut, Math.Clamp(value, Lo, Hi))) { Raise(nameof(Label)); Changed?.Invoke(); } } }
    public bool Never { get => _never; set { if (Set(ref _never, value)) { Raise(nameof(Label)); Changed?.Invoke(); Committed?.Invoke(); } } }
    public int Within => Counts.Select((c, k) => Lo + (k + 0.5) * Bin <= _cut ? c : 0).Sum();
    public string Label => _never ? (Ionic ? "ionic: none" : "never bond") : string.Format(CultureInfo.InvariantCulture, "≤ {0:0.00} Å · {1:N0}", _cut, Within);
}

/// <summary>A rule set kept by name (Settings): every pair's cut-off or never.</summary>
public sealed class BondRuleSet
{
    public string Name { get; set; } = "";
    public List<BondRuleItem> Rules { get; set; } = new();
}
public sealed class BondRuleItem
{
    public string A { get; set; } = "";
    public string B { get; set; } = "";
    public double Max { get; set; }
    public bool Never { get; set; }
}

/// <summary>Bond rules (design/boards/BondRules): bonds from what the structure actually holds — per pair of elements a
/// histogram of distances (0.8 – 3.2 Å) with the cut-off on it, dragged where the gap is; ionic pairs never bonded; the
/// bonds it makes counted live and the new ones lit in the view; rule sets kept by name; Bond it is one undoable step.</summary>
public sealed partial class MainViewModel
{
    public ObservableCollection<RulePair> RulePairs { get; } = new();
    private bool _rulesOpen;
    public bool RulesOpen
    {
        get => _rulesOpen;
        set
        {
            if (!Set(ref _rulesOpen, value)) return;
            if (value) { BrushOpen = false; AppearanceOpen = false; SelectionOpen = false; InteractionsOpen = false; HistoryOpen = false; StatesOpen = false; LensOpen = false; LodOpen = false; LoadRulePairs(); }
            else ClearRulesPreview();
            Raise(nameof(ShowRulesPanel));
            Raise(nameof(ShowStudioTabs));
        }
    }
    public bool ShowRulesPanel => IsStudio && _rulesOpen && _doc != null && !_appOpen && !_selOpen && !_ixOpen && !_lodOpen && !_histOpen && !_lensOpen && !_stOpen;
    private string _rulesSummary = "", _rulesNote = "", _ruleSetName = "";
    public string RulesSummary { get => _rulesSummary; private set => Set(ref _rulesSummary, value); }
    public string RulesNote { get => _rulesNote; private set => Set(ref _rulesNote, value); }
    public string RulesTitle => "Bond rules · " + Title.Replace(" (unsaved)", "");
    public List<string> RuleSetNames => ["The structure's gaps", .. _settings.BondRuleSets.Select(r => r.Name)];
    private int _ruleSetIndex;
    public int RuleSetIndex { get => _ruleSetIndex; set { if (Set(ref _ruleSetIndex, value) && value > 0) UseRuleSet(_settings.BondRuleSets[value - 1]); else if (value == 0) LoadRulePairs(); } }
    public string RuleSetName { get => _ruleSetName; set => Set(ref _ruleSetName, value ?? ""); }

    public void OpenRules() { SetModule(8); RulesOpen = true; }

    /// <summary>The structure's pairs with their histograms; each cut-off where its gap is (ionic pairs never).</summary>
    public void LoadRulePairs()
    {
        RulePairs.Clear();
        Raise(nameof(RulesTitle));
        if (_doc == null) return;
        try
        {
            var j = JsonNode.Parse(_doc.PairHistograms("{}"))!;
            foreach (var p in (j["pairs"] as JsonArray ?? []).OfType<JsonObject>())
            {
                var counts = (p["counts"] as JsonArray ?? []).Select(x => (int)((double?)x ?? 0)).ToArray();
                if (counts.Sum() == 0 && (int)((double?)p["bonded"] ?? 0) == 0) continue;   // no two of them within 3.2 Å
                var pair = new RulePair
                {
                    A = (string?)p["a"] ?? "", B = (string?)p["b"] ?? "", Counts = counts, Lo = (double?)p["lo"] ?? 0.8, Bin = (double?)p["bin"] ?? 0.04,
                    Bonded = (int)((double?)p["bonded"] ?? 0), Suggested = (double?)p["suggested"] ?? 1.8, Ionic = (bool?)p["ionic"] ?? false,
                };
                pair.Cutoff = pair.Suggested;
                pair.Never = pair.Ionic && pair.Bonded == 0;
                pair.Committed = PreviewRules;
                RulePairs.Add(pair);
            }
        }
        catch (Exception e) { RulesNote = e.Message; }
        PreviewRules();
    }

    private JsonObject RulesJson(bool preview) => new()
    {
        ["rules"] = new JsonArray(RulePairs.Select(p => (JsonNode)new JsonObject { ["a"] = p.A, ["b"] = p.B, ["max"] = p.Cutoff, ["never"] = p.Never }).ToArray()),
        ["preview"] = preview,
    };

    /// <summary>The bonds the rules make, against the structure's: counted, the new ones lit in the view.</summary>
    public void PreviewRules()
    {
        if (_doc == null) return;
        try
        {
            var r = JsonNode.Parse(_doc.BondRulesPreview(RulesJson(true).ToJsonString()))!;
            var inv = CultureInfo.InvariantCulture;
            int now = (int)((double?)r["bonds"] ?? 0), after = (int)((double?)r["after"] ?? 0), added = (int)((double?)r["added"] ?? 0), removed = (int)((double?)r["removed"] ?? 0);
            RulesSummary = string.Format(inv, "{0:N0} bonds → {1:N0}", now, after);
            var per = (r["pairs"] as JsonObject)?.Select(kv => $"{((double?)kv.Value ?? 0):+0;-0} {kv.Key}").ToList() ?? [];
            RulesNote = after == now && added == 0 ? "the same bonds as now"
                : string.Join(" · ", per) + (added > 0 ? " · the new ones lit in the view" : "") + (removed > 0 ? $" · {removed:N0} taken away" : "");
        }
        catch (Exception e) { RulesNote = e.Message; }
        RenderRequested?.Invoke();
    }

    private void ClearRulesPreview()
    {
        try { _doc?.BondRulesPreview("{}"); } catch { }
        RenderRequested?.Invoke();
    }

    /// <summary>Bond it: the structure's bonds from the rules (one undoable step); its force field must be assigned again.</summary>
    public void ApplyBondRules()
    {
        if (_doc == null) return;
        RunEdit(new { op = "bond_rules", rules = RulesJson(false)["rules"] });
        LoadRulePairs();
    }

    public void SaveRuleSet()
    {
        var name = _ruleSetName.Trim();
        if (name.Length == 0) { Status = "Name the rule set first"; return; }
        _settings.BondRuleSets.RemoveAll(r => r.Name == name);
        _settings.BondRuleSets.Add(new BondRuleSet { Name = name, Rules = RulePairs.Select(p => new BondRuleItem { A = p.A, B = p.B, Max = p.Cutoff, Never = p.Never }).ToList() });
        _settings.Save();
        Raise(nameof(RuleSetNames));
        _ruleSetIndex = _settings.BondRuleSets.FindIndex(r => r.Name == name) + 1;
        Raise(nameof(RuleSetIndex));
        Status = $"Rule set “{name}” kept: choose it for any structure";
    }

    /// <summary>A kept rule set on this structure's pairs (pairs it does not name keep their gap).</summary>
    private void UseRuleSet(BondRuleSet set)
    {
        foreach (var p in RulePairs)
            if (set.Rules.FirstOrDefault(r => (r.A == p.A && r.B == p.B) || (r.A == p.B && r.B == p.A)) is { } r)
            {
                p.Committed = null;
                p.Never = r.Never;
                p.Cutoff = r.Max;
                p.Committed = PreviewRules;
            }
        RuleSetName = set.Name;
        PreviewRules();
    }
}
