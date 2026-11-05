using System.Globalization;
using System.Text.Json.Nodes;

namespace CapsStudio.ViewModels;

/// <summary>Unit cell › Space group (core spacegroup.hpp): the group kept with the structure (read from its CIF, built, or
/// found), whether it still holds, the cell fitted to its crystal system on Apply, and moves mirrored to every symmetry
/// image of an atom.</summary>
public sealed partial class MainViewModel
{
    private string _sgText = "", _sgNote = "";
    private bool _sgKept, _sgFitCell = true, _sgMirror = true;
    public string SgText { get => _sgText; private set => Set(ref _sgText, value); }
    public string SgNote { get => _sgNote; private set => Set(ref _sgNote, value); }
    public bool SgKept { get => _sgKept; private set => Set(ref _sgKept, value); }
    /// <summary>Apply fits the six numbers to the kept group's crystal system (cubic a = b = c, hexagonal γ = 120° …).</summary>
    public bool SgFitCell { get => _sgFitCell; set => Set(ref _sgFitCell, value); }
    /// <summary>The move tool takes each atom's symmetry images with it while the group holds.</summary>
    public bool SgMirror { get => _sgMirror; set => Set(ref _sgMirror, value); }

    private void SpaceGroupCall(string op)
    {
        if (_doc == null) { SgText = ""; SgKept = false; return; }
        var inv = CultureInfo.InvariantCulture;
        try
        {
            var j = JsonNode.Parse(_doc.Symmetry(new JsonObject { ["op"] = op, ["tolerance"] = 0.1 }.ToJsonString()))!;
            if (j["ok"]?.GetValue<bool>() != true) { SgNote = (string?)j["error"] ?? "symmetry failed"; return; }
            SgKept = j["kept"]!.GetValue<bool>();
            if (!SgKept) { SgText = "None kept: the structure is P 1"; SgNote = "Find looks for the highest group its atoms have (0.1 Å)."; return; }
            var holds = j["holds"]!.GetValue<bool>();
            var number = (int)j["number"]!.GetValue<double>();
            SgText = $"{(string?)j["name"]}{(number > 0 ? $" (No. {number}, {(string?)j["system"]})" : "")} · {(int)j["operations"]!.GetValue<double>()} operations";
            SgNote = holds
                ? $"Holds · a CIF lists {(int)j["asymmetric"]!.GetValue<double>()} atoms with the operations"
                : $"Broken by an edit (an image {j["worst"]!.GetValue<double>().ToString("0.00", inv)} Å off): a CIF is written as P 1";
        }
        catch (Exception e) { SgNote = e.Message; }
    }

    public void RefreshSpaceGroup() => SpaceGroupCall("info");
    public void FindSpaceGroup() { if (!EditLocked()) { SpaceGroupCall("find"); AfterEdit("Space group kept: " + SgText + " (undo with ⌘Z)"); } }
    public void DropSpaceGroup() { if (!EditLocked()) { SpaceGroupCall("clear"); AfterEdit("Space group dropped (undo with ⌘Z)"); } }
}
