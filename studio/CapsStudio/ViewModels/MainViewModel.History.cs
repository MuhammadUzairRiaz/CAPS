using System;
using System.Collections.ObjectModel;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Text.Json.Nodes;
using Avalonia.Media.Imaging;
using CapsStudio.Interop;

namespace CapsStudio.ViewModels;

public sealed record HistoryRow(int Number, string What, string Atoms, string State, int Step)
{
    public bool IsCurrent => State == "current";
    public bool IsUndone => State == "undone";
    public string Note => State == "current" ? "you are here" : State == "undone" ? "undone" : "";
}
public sealed record BranchRow(int Index, string Title, string Detail);
public sealed class SnapshotRow
{
    public required int Index { get; init; }
    public required string Name { get; init; }
    public required string Detail { get; init; }
    public Bitmap? Thumb { get; set; }
}

/// <summary>History & snapshots (design/boards/History): every edit with the atoms after it, "you are here", undone steps
/// (kept as a branch when an edit follows an undo — nothing is lost until it is deleted), and named snapshots to go back
/// to or compare with the structure now.</summary>
public partial class MainViewModel
{
    private bool _histOpen;
    public bool HistoryOpen
    {
        get => _histOpen;
        set
        {
            if (!Set(ref _histOpen, value)) return;
            if (value) { AppearanceOpen = false; SelectionOpen = false; InteractionsOpen = false; LodOpen = false; RefreshHistoryPanel(); }
            Raise(nameof(ShowHistoryPanel)); Raise(nameof(ShowStudioTabs));
        }
    }
    public bool ShowHistoryPanel => IsStudio && _histOpen && _doc != null;
    public ObservableCollection<HistoryRow> HistoryRows { get; } = new();
    public ObservableCollection<BranchRow> HistoryBranches { get; } = new();
    public ObservableCollection<SnapshotRow> Snapshots { get; } = new();
    private readonly System.Collections.Generic.Dictionary<string, Bitmap> _snapThumbs = new();
    private string _historySummary = "";
    public string HistorySummary { get => _historySummary; private set => Set(ref _historySummary, value); }
    public bool HasBranches => HistoryBranches.Count > 0;
    public bool HasSnapshots => Snapshots.Count > 0;

    public void RefreshHistoryPanel()
    {
        HistoryRows.Clear(); HistoryBranches.Clear(); Snapshots.Clear();
        if (_doc == null) return;
        JsonNode h;
        try { h = JsonNode.Parse(_doc.History())!; } catch { return; }
        var inv = CultureInfo.InvariantCulture;
        string At(JsonNode? n) => $"{(n?.GetValue<double>() ?? 0).ToString("N0", inv)} atoms";
        var steps = h["steps"] as JsonArray ?? new JsonArray();
        var anyDone = steps.Any(s => (string?)s!["state"] != "undone");
        HistoryRows.Add(new HistoryRow(1, "Opened · " + Path.GetFileName(_doc.Path ?? Title), At(h["start_atoms"]), anyDone ? "done" : "current", 0));
        var k = 1;
        foreach (var s in steps)
        {
            ++k;
            HistoryRows.Add(new HistoryRow(k, (string?)s!["what"] ?? "", At(s["atoms"]), (string?)s["state"] ?? "done", k - 1));
        }
        var bi = 0;
        foreach (var b in h["branches"] as JsonArray ?? new JsonArray())
        {
            var st = ((JsonArray)b!["steps"]!).Select(x => (string?)x ?? "").ToArray();
            HistoryBranches.Add(new BranchRow(bi++, $"Branch · {st.Length} step{(st.Length == 1 ? "" : "s")}", string.Join(" → ", st)));
        }
        var si = 0;
        foreach (var s in h["snapshots"] as JsonArray ?? new JsonArray())
        {
            var name = (string?)s!["name"] ?? "";
            var detail = $"{At(s["atoms"])} · step {(int)(s["step"]?.GetValue<double>() ?? 0)}" + (s["on_branch"]?.GetValue<bool>() == true ? " · on a branch" : "");
            Snapshots.Add(new SnapshotRow { Index = si++, Name = name, Detail = detail, Thumb = _snapThumbs.GetValueOrDefault(name) });
        }
        var undone = HistoryRows.Count(r => r.IsUndone);
        HistorySummary = $"step {HistoryRows.Count - undone} of {HistoryRows.Count}" + (HistoryBranches.Count > 0 ? $" · {HistoryBranches.Count} branch{(HistoryBranches.Count == 1 ? "" : "es")}" : "");
        foreach (var n in new[] { nameof(HasBranches), nameof(HasSnapshots) }) Raise(n);
    }

    /// <summary>Undo or redo until row `step` (0 = as opened) is the current state.</summary>
    public void JumpToStep(int step)
    {
        if (_doc == null) return;
        var current = HistoryRows.LastOrDefault(r => !r.IsUndone)?.Step ?? 0;
        var ok = true;
        while (ok && current > step) { ok = _doc.Undo(false); current--; }
        while (ok && current < step) { ok = _doc.Undo(true); current++; }
        AfterEdit(current == step ? $"History: step {step + 1}" : "History: could not reach that step");
    }

    public void TakeSnapshot(string? name = null)
    {
        if (_doc == null) return;
        name ??= $"snapshot {Snapshots.Count + 1}";
        _doc.Snapshot(new JsonObject { ["op"] = "take", ["name"] = name }.ToJsonString());
        try
        {
            var o = new CapsRenderOpts { Width = 96, Height = 64, Supersample = 2, Background = 2, Style = _style, ColourBy = _colour, Outlines = 1, DepthCue = 1, ShowCell = 0,
                                         Highlight0 = -1, Highlight1 = -1, Highlight2 = -1, Highlight3 = -1 };
            var rgba = new byte[96 * 64 * 4];
            _doc.Render(Camera, o, rgba);
            _snapThumbs[name] = ToBitmap(rgba, 96, 64);
        }
        catch { /* a thumbnail is optional */ }
        RefreshHistoryPanel();
        Status = $"Snapshot “{name}” taken";
    }

    public void RestoreSnapshot(SnapshotRow r)
    {
        if (_doc == null) return;
        _doc.Snapshot(new JsonObject { ["op"] = "restore", ["index"] = r.Index }.ToJsonString());
        AfterEdit($"Restored snapshot “{r.Name}” (undo goes back)");
    }

    public void DeleteSnapshot(SnapshotRow r)
    {
        _doc?.Snapshot(new JsonObject { ["op"] = "delete", ["index"] = r.Index }.ToJsonString());
        _snapThumbs.Remove(r.Name);
        RefreshHistoryPanel();
    }

    /// <summary>Compare › the structure now (A) against the snapshot (B), both written to temporary files.</summary>
    public async void CompareSnapshot(SnapshotRow r)
    {
        if (_doc == null) return;
        var dir = Path.Combine(Path.GetTempPath(), "caps-snapshots");
        Directory.CreateDirectory(dir);
        var a = Path.Combine(dir, "now.data");
        var b = Path.Combine(dir, string.Concat(r.Name.Select(c => char.IsLetterOrDigit(c) ? c : '_')) + ".data");
        try
        {
            _doc.Save(a);
            _doc.Snapshot(new JsonObject { ["op"] = "save", ["index"] = r.Index, ["path"] = b }.ToJsonString());
        }
        catch (Exception e) { Status = "Could not write the snapshot: " + e.Message; return; }
        SetModule(23);
        await SetCompareInput(true, a);
        await SetCompareInput(false, b);
        Status = $"Comparing the structure now with snapshot “{r.Name}”";
    }

    /// <summary>The branch's steps become the redo steps (the current redo steps become a branch).</summary>
    public void RestoreBranch(BranchRow r)
    {
        _doc?.Snapshot(new JsonObject { ["op"] = "branch", ["index"] = r.Index }.ToJsonString());
        RefreshHistory();
        RefreshHistoryPanel();
        Status = "Branch restored: redo steps through it (⇧⌘Z)";
    }

    public void DeleteBranch(BranchRow r)
    {
        _doc?.Snapshot(new JsonObject { ["op"] = "drop_branch", ["index"] = r.Index }.ToJsonString());
        RefreshHistoryPanel();
    }
}
