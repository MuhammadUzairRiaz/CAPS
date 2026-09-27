using System;
using System.Collections.ObjectModel;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Text;
using System.Text.Json.Nodes;

namespace CapsStudio.ViewModels;

/// <summary>A state of the structure the compare panel can pick: the start, the structure now, a frame or a snapshot.</summary>
public sealed record StateChoice(string Label, string Kind, int Index)
{
    public override string ToString() => Label;
    public JsonObject Json() => new() { ["kind"] = Kind, ["index"] = Index };
}
public sealed record ShiftRow(string Label, string Where, string Shift);

/// <summary>Compare states (design/boards/Compare): one state of the structure superposed on another (Horn 1987
/// quaternions), RMSD over all, heavy and backbone atoms, each atom's shift in the view and the largest ones.</summary>
public partial class MainViewModel
{
    private bool _stOpen, _stColour = true;
    private int _stFit;   // 0 all, 1 heavy, 2 backbone, 3 selection
    private static readonly string[] StFitIds = ["all", "heavy", "backbone", "selection"];
    private StateChoice? _stRef, _stMov;
    private string _stRmsdAll = "—", _stRmsdHeavy = "—", _stRmsdBack = "—", _stInfo = "", _stError = "";
    private double[] _stShifts = [];

    public bool StatesOpen
    {
        get => _stOpen;
        set
        {
            if (!Set(ref _stOpen, value)) return;
            if (value) { AppearanceOpen = false; SelectionOpen = false; InteractionsOpen = false; LodOpen = false; LensOpen = false; HistoryOpen = false; RefreshStateChoices(); CompareStates(); }
            else ClearStateColours();
            Raise(nameof(ShowStatesPanel)); Raise(nameof(ShowStudioTabs));
        }
    }
    public bool ShowStatesPanel => IsStudio && _stOpen && _doc != null;
    public ObservableCollection<StateChoice> StateChoices { get; } = new();
    public ObservableCollection<ShiftRow> LargestShifts { get; } = new();
    public StateChoice? CompareReference { get => _stRef; set { if (Set(ref _stRef, value)) CompareStates(); } }
    public StateChoice? CompareMoving { get => _stMov; set { if (Set(ref _stMov, value)) CompareStates(); } }
    public int CompareFit { get => _stFit; set { if (Set(ref _stFit, Math.Clamp(value, 0, 3))) CompareStates(); } }
    public bool CompareColour { get => _stColour; set { if (Set(ref _stColour, value)) { if (value) CompareStates(); else ClearStateColours(); } } }
    public string CompareRmsdAll { get => _stRmsdAll; private set => Set(ref _stRmsdAll, value); }
    public string CompareRmsdHeavy { get => _stRmsdHeavy; private set => Set(ref _stRmsdHeavy, value); }
    public string CompareRmsdBackbone { get => _stRmsdBack; private set => Set(ref _stRmsdBack, value); }
    public string CompareInfo { get => _stInfo; private set => Set(ref _stInfo, value); }
    public string CompareError { get => _stError; private set { if (Set(ref _stError, value)) Raise(nameof(CompareHasError)); } }
    public bool CompareHasError => _stError.Length > 0;
    /// <summary>Each atom's shift after the fit (Å), for the plot under the panel.</summary>
    public double[] CompareShifts => _stShifts;
    public event Action? CompareStatesChanged;

    /// <summary>The states to pick from: the start, the structure now, each frame (up to 200) and each snapshot.</summary>
    public void RefreshStateChoices()
    {
        var keepRef = _stRef;
        var keepMov = _stMov;
        StateChoices.Clear();
        if (_doc == null) return;
        StateChoices.Add(new StateChoice("Start (as opened, before edits)", "start", 0));
        StateChoices.Add(new StateChoice("Now (the structure shown)", "current", 0));
        var frames = (int)_doc.Summary().Frames;
        if (frames > 1)
            for (var f = 0; f < Math.Min(frames, 200); f++) StateChoices.Add(new StateChoice($"Frame {f + 1}", "frame", f));
        try
        {
            var h = JsonNode.Parse(_doc.History())!;
            var i = 0;
            foreach (var s in h["snapshots"] as JsonArray ?? new JsonArray())
                StateChoices.Add(new StateChoice($"Snapshot · {(string?)s!["name"] ?? ""}", "snapshot", i++));
        }
        catch { /* no history */ }
        StateChoice? Find(StateChoice? c) => c == null ? null : StateChoices.FirstOrDefault(x => x.Kind == c.Kind && x.Index == c.Index);
        _stRef = Find(keepRef) ?? (frames > 1 ? StateChoices.First(x => x.Kind == "frame") : StateChoices[0]);
        _stMov = Find(keepMov) ?? StateChoices.First(x => x.Kind == "current");
        Raise(nameof(CompareReference)); Raise(nameof(CompareMoving));
    }

    /// <summary>Compare › open the panel with the structure now against this snapshot.</summary>
    public void CompareWithSnapshot(SnapshotRow r)
    {
        if (_doc == null) return;
        if (!IsStudio) SetModule(8);
        _stOpen = false;
        StatesOpen = true;
        var pick = StateChoices.FirstOrDefault(x => x.Kind == "snapshot" && x.Index == r.Index);
        if (pick != null) CompareReference = pick;
    }

    /// <summary>Compare › both states side by side: the reference on the right of the split view, the moving state on the
    /// left (the structure shown, at that frame when it is one).</summary>
    public void OpenStatesInSplit()
    {
        if (_doc == null || _stRef == null || _stMov == null) return;
        try
        {
            var right = _doc.StateDocument(_stRef.Json().ToJsonString(), "reference");
            if (_stMov.Kind == "frame") Frame = _stMov.Index;
            else if (_stMov.Kind != "current")
            {
                right.Dispose();
                CompareError = "the split view shows the structure itself on the left: choose Now or a frame as the moving state";
                return;
            }
            OpenSplit();
            SetSplitBDocument(right, $"{_stRef.Label} (reference)");
            SyncCamera = true;
            Status = $"Split view: {_stMov.Label} (left) against {_stRef.Label} (right), cameras synced";
        }
        catch (Exception e) { CompareError = e.Message; }
    }

    public void CompareStates()
    {
        if (_doc == null || !_stOpen || _stRef == null || _stMov == null) return;
        var o = new JsonObject
        {
            ["reference"] = _stRef.Json(), ["moving"] = _stMov.Json(), ["fit"] = StFitIds[_stFit],
            ["largest"] = 12, ["per_atom"] = 1, ["colour"] = _stColour && _stMov.Kind == "current" ? 1 : 0,
        };
        JsonNode r;
        try { r = JsonNode.Parse(_doc.CompareStates(o.ToJsonString()))!; }
        catch (Exception e) { r = new JsonObject { ["ok"] = false, ["error"] = e.Message }; }
        LargestShifts.Clear();
        var inv = CultureInfo.InvariantCulture;
        if (r["ok"]?.GetValue<bool>() != true)
        {
            CompareError = (string?)r["error"] ?? "could not compare";
            CompareRmsdAll = CompareRmsdHeavy = CompareRmsdBackbone = "—";
            _stShifts = [];
            CompareInfo = "";
            CompareStatesChanged?.Invoke();
            return;
        }
        CompareError = "";
        string A(JsonNode? n) => n == null ? "—" : n.GetValue<double>().ToString("0.000", inv) + " Å";
        CompareRmsdAll = A(r["rmsd"]?["all"]);
        CompareRmsdHeavy = A(r["rmsd"]?["heavy"]);
        CompareRmsdBackbone = r["rmsd"]?["backbone"] is { } b ? A(b) : "no backbone";
        foreach (var e in r["largest"] as JsonArray ?? new JsonArray())
            LargestShifts.Add(new ShiftRow((string?)e!["label"] ?? "",
                (e["backbone"]?.GetValue<bool>() == true ? "backbone · " : "") + $"molecule {(int)(e["molecule"]?.GetValue<double>() ?? 0)}",
                e["shift"]!.GetValue<double>().ToString("0.00", inv)));
        _stShifts = (r["shifts"] as JsonArray ?? new JsonArray()).Select(x => x!.GetValue<double>()).ToArray();
        CompareInfo = string.Format(inv, "{0:N0} atoms · fit on {1} ({2:N0} atoms) · mean shift {3:0.00} Å, largest {4:0.00} Å · Horn 1987 quaternion superposition",
            r["atoms"]!.GetValue<double>(), StFitIds[_stFit], r["fitted"]!.GetValue<double>(), r["mean"]!.GetValue<double>(), r["max"]!.GetValue<double>());
        if (_stColour && _stMov.Kind != "current") CompareInfo += " · the view colours shifts when the moving state is Now";
        CompareStatesChanged?.Invoke();
        RequestRenderAfterColours();
    }

    private void ClearStateColours()
    {
        try { _doc?.CompareStates("{\"op\":\"clear\"}"); } catch { }
        RequestRenderAfterColours();
    }

    private void RequestRenderAfterColours() => RenderRequested?.Invoke();

    /// <summary>Atom, element, molecule and shift of every atom as CSV.</summary>
    public void ExportShiftsCsv(string path)
    {
        if (_doc == null || _stShifts.Length == 0) return;
        var inv = CultureInfo.InvariantCulture;
        var sb = new StringBuilder();
        sb.Append(string.Format(inv, "# shift after superposing {0} on {1} (fit on {2}); RMSD all {3}, heavy {4}, backbone {5}\n",
            _stMov?.Label, _stRef?.Label, StFitIds[_stFit], CompareRmsdAll, CompareRmsdHeavy, CompareRmsdBackbone));
        sb.Append("atom,element,molecule,shift_A\n");
        for (var i = 0; i < _stShifts.Length; i++)
        {
            var a = _doc.Atom(i);
            sb.Append(string.Format(inv, "{0},{1},{2},{3:0.#####}\n", i + 1, a.ElementSymbol, a.Mol, _stShifts[i]));
        }
        File.WriteAllText(path, sb.ToString());
        Status = $"Wrote {Path.GetFileName(path)}";
    }
}
