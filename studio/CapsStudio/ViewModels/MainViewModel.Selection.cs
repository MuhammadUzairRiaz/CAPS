using System.Collections.ObjectModel;
using System.Globalization;
using System.Text.Json.Nodes;
using CapsStudio.Interop;

namespace CapsStudio.ViewModels;

/// <summary>A named selection: its atoms and a colour for its dot.</summary>
public sealed record NamedSet(string Name, int[] Atoms, string Colour)
{
    public string CountText => $"{Atoms.Length:N0} atoms";
}

/// <summary>A dyad of the tacticity strip: m (meso) or r (racemo) between two consecutive stereocentres.</summary>
public sealed record DyadCell(string Kind, string Label)
{
    public bool IsMeso => Kind == "m";
}

/// <summary>Selection &amp; stereo (design/boards/SelectionStereo): select by SMARTS, element, type, charge range,
/// distance or along bonds; named sets; stereocentres and tacticity (m/r dyads, triads) with make-isotactic and
/// make-syndiotactic.</summary>
public sealed partial class MainViewModel
{
    public static readonly string[] SelectModes = ["SMARTS", "Element", "Atom type", "Charge range", "Within distance", "Grow along bonds"];
    private static readonly string[] SelectModeCodes = ["smarts", "element", "type", "charge", "within", "grow"];
    private static readonly string[] SetColours = ["#6CC4D8", "#F5A524", "#9AA1A8", "#9B7BD6", "#7CC784", "#E07A5F"];

    private bool _selOpen;
    public bool SelectionOpen { get => _selOpen; set { if (Set(ref _selOpen, value)) { if (value) { AppearanceOpen = false; InteractionsOpen = false; LodOpen = false; HistoryOpen = false; RefreshStereo(); } RaiseSelectionPanel(); } } }
    public bool ShowSelectionPanel => IsStudio && _selOpen && _doc != null;
    private void RaiseSelectionPanel() { Raise(nameof(ShowSelectionPanel)); Raise(nameof(ShowStudioTabs)); }

    private int _selMode;
    public int SelectMode { get => _selMode; set { if (Set(ref _selMode, Math.Clamp(value, 0, 5))) { Raise(nameof(SelectPatternLabel)); Raise(nameof(SelectPatternHint)); Raise(nameof(SelectModeIsSmarts)); Raise(nameof(SelectModeIsCharge)); } } }
    public bool SelectModeIsSmarts => _selMode == 0;
    public bool SelectModeIsCharge => _selMode == 3;
    private string _selPattern = "c1ccccc1";
    public string SelectPattern { get => _selPattern; set => Set(ref _selPattern, value ?? ""); }
    private decimal _selLo = -0.2m, _selHi = 0.2m, _selDist = 3.0m, _selSteps = 1;
    public decimal SelectChargeLo { get => _selLo; set => Set(ref _selLo, value); }
    public decimal SelectChargeHi { get => _selHi; set => Set(ref _selHi, value); }
    public decimal SelectDistance { get => _selDist; set => Set(ref _selDist, Math.Clamp(value, 0.5m, 100m)); }
    public decimal SelectSteps { get => _selSteps; set => Set(ref _selSteps, Math.Clamp(value, 1, 50)); }
    public string SelectPatternLabel => _selMode switch { 0 => "SMARTS pattern", 1 => "Elements", 2 => "Atom type", 3 => "Charge range · e", 4 => "Distance from the selection · Å", _ => "Bonds to grow by" };
    public string SelectPatternHint => _selMode switch
    {
        0 => "The pattern's first atom is selected: c1ccccc1 takes every ring atom",
        1 => "Symbols, e.g. C N O",
        2 => "A type label (c3, ca) or number",
        3 => "Atoms whose partial charge lies in the range",
        4 => "Every atom within the distance of an atom already selected",
        _ => "The selection grown along its bonds",
    };

    private int _selCount;
    public int SelectedCount { get => _selCount; private set { if (Set(ref _selCount, value)) { Raise(nameof(SelectedChip)); Raise(nameof(SelectionStatus)); } } }
    public string SelectedChip => $"{_selCount:N0} selected";
    public string SelectionStatus => $"{_doc?.Summary().Atoms ?? 0:N0} atoms · {_selCount:N0} selected";
    private string _selError = "";
    public string SelectError { get => _selError; private set { if (Set(ref _selError, value)) Raise(nameof(SelectHasError)); } }
    public bool SelectHasError => _selError.Length > 0;
    private string _selHud = "";
    public string SelectHud { get => _selHud; private set => Set(ref _selHud, value); }

    /// <summary>Runs the chosen mode: op replace | add | subtract | intersect | invert.</summary>
    public void RunSelect(string op)
    {
        if (_doc == null) return;
        var o = new JsonObject { ["mode"] = SelectModeCodes[_selMode], ["pattern"] = _selPattern, ["op"] = op };
        if (_selMode == 3) { o["lo"] = (double)_selLo; o["hi"] = (double)_selHi; }
        if (_selMode == 4) o["distance"] = (double)_selDist;
        if (_selMode == 5) o["steps"] = (int)_selSteps;
        if (op == "invert") o = new JsonObject { ["mode"] = "all", ["op"] = "invert" };
        var r = JsonNode.Parse(_doc.Select(o.ToJsonString()))!;
        if (r["ok"]?.GetValue<bool>() != true) { SelectError = r["error"]?.GetValue<string>() ?? "cannot select"; return; }
        SelectError = "";
        SelectedCount = (int)r["count"]!.GetValue<double>();
        SelectHud = op == "invert" ? "inverted" : $"{SelectModes[_selMode]} {(_selMode is 0 or 1 or 2 ? _selPattern : "")}".Trim();
        RenderRequested?.Invoke();
    }

    public void ClearDocSelection()
    {
        if (_doc == null) return;
        _doc.Select("{\"mode\":\"none\"}");
        SelectedCount = 0;
        SelectHud = "";
        RenderRequested?.Invoke();
    }

    public ObservableCollection<NamedSet> NamedSets { get; } = new();

    public void SaveSelectionAsSet(string? name = null)
    {
        if (_doc == null) return;
        var j = JsonNode.Parse(_doc.SelectionJson())!;
        var idx = j["indices"]!.AsArray().Select(x => (int)x!.GetValue<double>()).ToArray();
        if (idx.Length == 0) { SelectError = "Nothing is selected"; return; }
        name ??= _selMode is 0 or 1 or 2 && _selPattern.Length > 0 ? _selPattern : $"set-{NamedSets.Count + 1}";
        foreach (var old in NamedSets.Where(s => s.Name == name).ToList()) NamedSets.Remove(old);
        NamedSets.Add(new NamedSet(name, idx, SetColours[NamedSets.Count % SetColours.Length]));
    }

    public void UseNamedSet(NamedSet set)
    {
        if (_doc == null) return;
        var r = JsonNode.Parse(_doc.Select(new JsonObject { ["mode"] = "indices", ["atoms"] = new JsonArray(set.Atoms.Select(a => (JsonNode)a).ToArray()), ["op"] = "replace" }.ToJsonString()))!;
        SelectedCount = (int)(r["count"]?.GetValue<double>() ?? 0);
        SelectHud = set.Name;
        RenderRequested?.Invoke();
    }

    // ---------------------------------------------------------------- stereo

    public ObservableCollection<DyadCell> Dyads { get; } = new();
    private string _tacLabel = "", _dyadCounts = "", _triadCounts = "";
    private int _stereoCentres;
    public int StereoCentres { get => _stereoCentres; private set => Set(ref _stereoCentres, value); }
    public string TacticityLabel { get => _tacLabel; private set => Set(ref _tacLabel, value); }
    public string DyadCounts { get => _dyadCounts; private set => Set(ref _dyadCounts, value); }
    public string TriadCounts { get => _triadCounts; private set => Set(ref _triadCounts, value); }
    public bool HasDyads => Dyads.Count > 0;
    public string StereoChainText { get; private set; } = "";

    public void RefreshStereo()
    {
        Dyads.Clear();
        if (_doc == null) return;
        try
        {
            var t = JsonNode.Parse(_doc.Tacticity())!;
            StereoCentres = (int)t["centres"]!.GetValue<double>();
            var label = t["label"]!.GetValue<string>();
            TacticityLabel = label.Length == 0 ? "no vinyl stereocentres" : label;
            var chains = t["chains"]!.AsArray();
            // the first chain's dyads in the strip (up to 60)
            if (chains.Count > 0)
            {
                var d = chains[0]!["dyads"]!.GetValue<string>();
                for (var k = 0; k < d.Length && k < 60; k++) Dyads.Add(new DyadCell(d[k].ToString(), $"{k + 1}–{k + 2}"));
                StereoChainText = chains.Count > 1 ? $"chain 1 of {chains.Count}" : "the chain";
            }
            DyadCounts = $"m = {t["m"]!.GetValue<double>():0} · r = {t["r"]!.GetValue<double>():0}";
            TriadCounts = $"triads mm {t["mm"]!.GetValue<double>():0} · mr {t["mr"]!.GetValue<double>():0} · rr {t["rr"]!.GetValue<double>():0}";
        }
        catch (Exception e) { SelectError = e.Message; }
        Raise(nameof(HasDyads));
        Raise(nameof(StereoChainText));
    }

    public async Task MakeTactic(bool iso)
    {
        if (_doc == null) return;
        var doc = _doc;
        Status = iso ? "Making the chains isotactic…" : "Making the chains syndiotactic…";
        var text = await Task.Run(() => doc.Edit(iso ? "{\"op\":\"tacticity\",\"to\":\"isotactic\"}" : "{\"op\":\"tacticity\",\"to\":\"syndiotactic\"}"));
        var r = JsonNode.Parse(text)!;
        if (r["ok"]?.GetValue<bool>() != true) { SelectError = r["error"]?.GetValue<string>() ?? "cannot change the tacticity"; Status = SelectError; return; }
        Record($"doc.edit(op=\"tacticity\", to=\"{(iso ? "isotactic" : "syndiotactic")}\")");
        AfterEdit(r["what"]!.GetValue<string>());
        RefreshStereo();
    }
}
