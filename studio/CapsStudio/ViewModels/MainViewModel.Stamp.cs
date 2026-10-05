using System.Collections.ObjectModel;
using System.Globalization;
using System.Text.Json.Nodes;
using CapsStudio.Interop;

namespace CapsStudio.ViewModels;

/// <summary>A piece in the clipboard tray: atoms copied out (or a file dropped in) with their bonds, types and charges.</summary>
public sealed class TrayPiece : ObservableObject
{
    public string Name { get; init; } = "";
    public string Detail { get; init; } = "";
    public string Json { get; init; } = "";
    private bool _armed;
    public bool Armed { get => _armed; set => Set(ref _armed, value); }
}

/// <summary>Clipboard tray and stamps (design/boards/Stamp): Copy keeps the last pieces in a tray under the view; one
/// picked becomes a stamp under the pointer — click to place it, the wheel turns it, ⌥ click places many — pushed clear
/// of what is there. A file dragged onto the window: into this structure as a stamp, or open as a new structure.</summary>
public sealed partial class MainViewModel
{
    public ObservableCollection<TrayPiece> Tray { get; } = new();
    public bool HasTray => Tray.Count > 0;
    private bool _trayShown = true;
    public bool TrayShown { get => _trayShown; set { if (Set(ref _trayShown, value)) Raise(nameof(ShowTray)); } }
    public bool ShowTray => _trayShown && Tray.Count > 0 && IsStudio && _doc != null;
    private const int TrayMax = 8;
    private static string TrayFile => Path.Combine(AppSettings.Override != null ? Path.GetDirectoryName(AppSettings.Override)! : AppSettings.Folder, "tray.json");

    private TrayPiece? _stamp;
    private double _stampDegrees;
    public bool StampArmed => _stamp != null;
    public string StampHint => _stamp == null ? "" : string.Format(CultureInfo.InvariantCulture,
        "stamp {0} · click to place · the wheel turns it ({1:0}°) · ⌥ click places many · Esc stops · 1.5 Å clear", _stamp.Name, _stampDegrees);

    public void LoadTray()
    {
        Tray.Clear();
        try
        {
            if (File.Exists(TrayFile) && JsonNode.Parse(File.ReadAllText(TrayFile))?["pieces"] is JsonArray a)
                foreach (var p in a.OfType<JsonObject>())
                    Tray.Add(new TrayPiece { Name = (string?)p["name"] ?? "piece", Detail = (string?)p["detail"] ?? "", Json = (string?)p["json"] ?? "" });
        }
        catch { }
        RaiseTray();
    }

    private void SaveTray()
    {
        try
        {
            Directory.CreateDirectory(Path.GetDirectoryName(TrayFile)!);
            var a = new JsonArray(Tray.Select(p => (JsonNode)new JsonObject { ["name"] = p.Name, ["detail"] = p.Detail, ["json"] = p.Json }).ToArray());
            File.WriteAllText(TrayFile, new JsonObject { ["format"] = "caps-tray", ["pieces"] = a }.ToJsonString());
        }
        catch { }
    }

    private void RaiseTray() { Raise(nameof(HasTray)); Raise(nameof(ShowTray)); }

    private void AddToTray(string json, string source)
    {
        var j = JsonNode.Parse(json)!;
        var atoms = (j["atoms"] as JsonArray)?.Count ?? 0;
        var mols = (int)((double?)j["molecules"] ?? 1);
        var name = (string?)j["name"] ?? "piece";
        Tray.Insert(0, new TrayPiece
        {
            Name = name, Json = json,
            Detail = $"{atoms:N0} atoms" + (mols > 1 ? $" · {mols} molecules" : "") + (source.Length > 0 ? " · from " + source : ""),
        });
        while (Tray.Count > TrayMax) Tray.RemoveAt(Tray.Count - 1);
        TrayShown = true;
        SaveTray();
        RaiseTray();
    }

    /// <summary>Copy: the selection (else the whole structure) kept in the tray, named by what it is.</summary>
    public void CopyToTray()
    {
        if (_doc == null) return;
        var a = SelectionAtoms();
        try
        {
            var name = a.Length == 0 ? Title.Replace(" (unsaved)", "") : DescribeAtoms(a);
            var json = _doc.Piece(new JsonObject { ["atoms"] = a.Length == 0 ? "all" : new JsonArray(a.Select(i => (JsonNode)i).ToArray()), ["name"] = name }.ToJsonString());
            AddToTray(json, Title.Replace(" (unsaved)", ""));
            Status = $"{(a.Length == 0 ? "The structure" : $"{a.Length:N0} atoms")} in the clipboard tray · pick it there to place it as a stamp";
        }
        catch (Exception e) { Status = "Could not copy: " + e.Message; }
    }

    public void RemoveFromTray(TrayPiece p)
    {
        if (p == _stamp) DisarmStamp();
        Tray.Remove(p);
        SaveTray();
        RaiseTray();
    }

    /// <summary>A piece picked in the tray: the stamp under the pointer (picked again: put away).</summary>
    public void ArmStamp(TrayPiece p)
    {
        if (_doc == null) { Status = "Open a structure to place the piece in"; return; }
        if (_stamp == p) { DisarmStamp(); return; }
        if (_stamp != null) _stamp.Armed = false;
        _stamp = p;
        p.Armed = true;
        Raise(nameof(StampArmed)); Raise(nameof(StampHint));
        Status = StampHint;
    }

    public void DisarmStamp()
    {
        if (_stamp != null) _stamp.Armed = false;
        _stamp = null;
        _stampDegrees = 0;
        Raise(nameof(StampArmed)); Raise(nameof(StampHint));
    }

    /// <summary>The wheel turns the stamp about the view axis, 15° a notch.</summary>
    public void TurnStamp(double notches)
    {
        _stampDegrees = ((_stampDegrees + 15 * Math.Round(notches)) % 360 + 360) % 360;
        Raise(nameof(StampHint));
        Status = StampHint;
    }

    /// <summary>The stamp placed at a point (Å) of the view plane, turned about the view axis; many: it stays armed.</summary>
    public void PlaceStamp(double[] at, double[] viewAxis, bool many)
    {
        if (_stamp == null || _doc == null || at.Length != 3) return;
        RunEdit(new { op = "stamp", piece = _stamp.Json, at, axis = viewAxis, degrees = _stampDegrees, clear = 1.5, name = _stamp.Name });
        RefreshTags();
        if (!many) DisarmStamp();
    }

    /// <summary>A file dropped on the into-this-structure target: it becomes the stamp (and joins the tray).</summary>
    public void DropIntoStructure(string path)
    {
        if (_doc == null) { Status = "Open a structure first, or drop the file on Open as a new structure"; return; }
        try
        {
            AddToTray(CapsDocument.PieceFromFile(path), Path.GetFileName(path));
            ArmStamp(Tray[0]);
        }
        catch (Exception e) { Status = $"Could not read {Path.GetFileName(path)}: {e.Message}"; }
    }

    // ---- the drop targets shown while a file is dragged over the window
    private bool _dropTargets;
    public bool DropTargets { get => _dropTargets; set => Set(ref _dropTargets, value); }
    public string DropIntoText => _doc == null ? "" : "Drop into " + Title.Replace(" (unsaved)", "");
}
