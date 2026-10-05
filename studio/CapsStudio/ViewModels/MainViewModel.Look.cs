using System.Globalization;
using System.Text.Json.Nodes;

namespace CapsStudio.ViewModels;

/// <summary>Look (design/boards/Look): a button on the view opens the sizes — atoms, sticks, space-filling spheres, lines —
/// for everything or only the selection, bond orders drawn (double and triple as parallel sticks, aromatic with a dashed
/// inner one) and arcs on angle and torsion monitors; the view follows live. Keep as my look makes them the default.</summary>
public sealed partial class MainViewModel
{
    private double _lookAtom = 0.28, _lookStick = 0.14, _lookSpace = 1.0, _lookLine = 1.4, _lookSel = 1.0;
    private bool _lookOrders, _lookArcs;
    private int _lookScope;   // 0 all, 1 the selection
    private bool _lookLoaded;

    public int LookScope { get => _lookScope; set { if (Set(ref _lookScope, Math.Clamp(value, 0, 1))) { Raise(nameof(LookAll)); Raise(nameof(LookSelection)); } } }
    public bool LookAll => _lookScope == 0;
    public bool LookSelection => _lookScope == 1;
    public double LookAtom { get => _lookAtom; set { if (Set(ref _lookAtom, Math.Round(value, 2))) { Raise(nameof(LookAtomText)); ApplyLook(); } } }
    public double LookStick { get => _lookStick; set { if (Set(ref _lookStick, Math.Round(value, 2))) { Raise(nameof(LookStickText)); ApplyLook(); } } }
    public double LookSpace { get => _lookSpace; set { if (Set(ref _lookSpace, Math.Round(value, 2))) { Raise(nameof(LookSpaceText)); ApplyLook(); } } }
    public double LookLine { get => _lookLine; set { if (Set(ref _lookLine, Math.Round(value, 1))) { Raise(nameof(LookLineText)); ApplyLook(); } } }
    public double LookSel { get => _lookSel; set { if (Set(ref _lookSel, Math.Round(value, 2))) { Raise(nameof(LookSelText)); ApplyLookSelection(); } } }
    public bool LookOrders { get => _lookOrders; set { if (Set(ref _lookOrders, value)) ApplyLook(); } }
    public bool LookArcs { get => _lookArcs; set { if (Set(ref _lookArcs, value)) RenderRequested?.Invoke(); } }
    public string LookAtomText => _lookAtom.ToString("0.00", Inv) + " × r vdW";
    public string LookStickText => _lookStick.ToString("0.00", Inv) + " Å";
    public string LookSpaceText => _lookSpace.ToString("0.00", Inv) + " × r vdW";
    public string LookLineText => _lookLine.ToString("0.0", Inv) + " px";
    public string LookSelText => _lookSel.ToString("0.00", Inv) + " ×";

    /// <summary>The user's own look as the starting sizes (once, from the settings).</summary>
    private void LoadLook()
    {
        if (_lookLoaded) return;
        _lookLoaded = true;
        var l = _settings.Look;
        _lookAtom = l.AtomScale; _lookStick = l.Stick; _lookSpace = l.Space; _lookLine = l.Line; _lookOrders = l.Orders; _lookArcs = l.Arcs;
        foreach (var n in new[] { nameof(LookAtom), nameof(LookStick), nameof(LookSpace), nameof(LookLine), nameof(LookOrders), nameof(LookArcs),
                                  nameof(LookAtomText), nameof(LookStickText), nameof(LookSpaceText), nameof(LookLineText) }) Raise(n);
    }

    /// <summary>The sizes to the open structure (the view redraws).</summary>
    public void ApplyLook()
    {
        LoadLook();
        if (_doc == null || _doc.IsDisposed) return;
        try
        {
            _doc.SetLook(new JsonObject
            {
                ["atom_scale"] = _lookAtom, ["bond_radius"] = _lookStick, ["space_scale"] = _lookSpace, ["line_px"] = _lookLine, ["bond_orders"] = _lookOrders,
            }.ToJsonString());
        }
        catch (Exception e) { Status = "Look: " + e.Message; }
        RenderRequested?.Invoke();
    }

    /// <summary>Selection scope: the selected atoms (and the sticks between them) at their own size.</summary>
    private void ApplyLookSelection()
    {
        if (_doc == null) return;
        var a = SelectionAtoms();
        if (a.Length == 0) { Status = "Select the atoms to size on their own (Look › Selection)"; return; }
        try { _doc.SetLook(new JsonObject { ["factor"] = new JsonObject { ["atoms"] = new JsonArray(a.Select(i => (JsonNode)i).ToArray()), ["value"] = _lookSel } }.ToJsonString()); }
        catch (Exception e) { Status = "Look: " + e.Message; }
        RenderRequested?.Invoke();
    }

    /// <summary>Keep as my look: these sizes are every new structure's.</summary>
    public void KeepLook()
    {
        _settings.Look = new LookSettings { AtomScale = _lookAtom, Stick = _lookStick, Space = _lookSpace, Line = _lookLine, Orders = _lookOrders, Arcs = _lookArcs };
        _settings.Save();
        Status = "This look is your default now (Settings keep it)";
    }

    /// <summary>Reset: CAPS's sizes, no bond orders or arcs, the selection's own size gone.</summary>
    public void ResetLook()
    {
        _lookAtom = 0.28; _lookStick = 0.14; _lookSpace = 1.0; _lookLine = 1.4; _lookOrders = false; _lookArcs = false; _lookSel = 1.0;
        foreach (var n in new[] { nameof(LookAtom), nameof(LookStick), nameof(LookSpace), nameof(LookLine), nameof(LookOrders), nameof(LookArcs), nameof(LookSel),
                                  nameof(LookAtomText), nameof(LookStickText), nameof(LookSpaceText), nameof(LookLineText), nameof(LookSelText) }) Raise(n);
        try { _doc?.SetLook("{\"clear_factors\":true}"); } catch { }
        ApplyLook();
    }
}
