using System.Globalization;
using System.Text.Json.Nodes;

namespace CapsStudio.ViewModels;

/// <summary>Atom bubble (design/boards/AtomBubble): at the picked atoms, their element ring, type, partial charge, formal
/// charge and a nudge. With several atoms a value they do not share reads "mixed" and stays each atom's own until set.
/// Enter applies; each change is one undo step. The orient puck turns the view along the selection's long axis, onto its
/// plane, or side-on.</summary>
public sealed partial class MainViewModel
{
    public static readonly string[] BubbleElements = ["C", "H", "N", "O", "S", "Si", "B", "Zn"];
    private int[] _bubbleAtoms = [];
    private string _bubbleTitle = "", _bubbleType = "", _bubbleCharge = "", _bubbleFormal = "", _bubbleElement = "";
    public string BubbleTitle { get => _bubbleTitle; private set => Set(ref _bubbleTitle, value); }
    public string BubbleType { get => _bubbleType; private set => Set(ref _bubbleType, value); }
    public string BubbleElement { get => _bubbleElement; private set => Set(ref _bubbleElement, value); }
    public string BubbleCharge { get => _bubbleCharge; set => Set(ref _bubbleCharge, value ?? ""); }
    public string BubbleFormal { get => _bubbleFormal; set => Set(ref _bubbleFormal, value ?? ""); }

    /// <summary>The bubble filled from the selection (the atoms it edits).</summary>
    public bool OpenBubble()
    {
        var a = SelectionAtoms();
        if (_doc == null || a.Length == 0) { Status = "Pick the atoms to edit"; return false; }
        _bubbleAtoms = a;
        var inv = CultureInfo.InvariantCulture;
        var atoms = a.Take(5000).Select(i => _doc.Atom(i)).ToList();
        string Same(IEnumerable<string> v) { var d = v.Distinct().ToList(); return d.Count == 1 ? d[0] : "mixed"; }
        BubbleElement = Same(atoms.Select(x => x.ElementSymbol));
        string[] types = [];
        try { types = (JsonNode.Parse(_doc.AtomLabels("type")) as JsonArray ?? JsonNode.Parse(_doc.AtomLabels("type"))?["labels"] as JsonArray ?? []).Select(x => (string?)x ?? "").ToArray(); } catch { }
        BubbleType = types.Length > 0 ? Same(a.Where(i => i < types.Length).Select(i => types[i])) : "";
        BubbleCharge = Same(atoms.Select(x => x.Charge.ToString("0.000", inv)));
        BubbleFormal = Same(atoms.Select(x => Math.Abs(x.Charge - Math.Round(x.Charge)) < 0.05 ? ((int)Math.Round(x.Charge)).ToString("+0;-0;0", inv) : "—"));
        var where = DescribeAtoms(a);
        BubbleTitle = a.Length == 1 ? $"{atoms[0].ElementSymbol}{a[0] + 1} · {where}" : $"{a.Length:N0} atoms · {where}";
        return true;
    }

    public void BubbleSetElement(string symbol)
    {
        if (_bubbleAtoms.Length == 0) return;
        RunEdit(new { op = "element", atoms = _bubbleAtoms, element = symbol });
        Reselect(_bubbleAtoms);
        OpenBubble();
    }

    /// <summary>Enter in the charge box: every atom of the bubble gets it ("mixed" left: each keeps its own).</summary>
    public void BubbleApplyCharge()
    {
        if (_bubbleAtoms.Length == 0 || !double.TryParse(_bubbleCharge.Replace('−', '-'), NumberStyles.Float, CultureInfo.InvariantCulture, out var q)) return;
        RunEdit(new { op = "charge", atoms = _bubbleAtoms, charge = q });
        Reselect(_bubbleAtoms);
        OpenBubble();
    }

    /// <summary>A formal charge: the atom's charge as that whole number (CAPS counts hydrogens and valence from it).</summary>
    public void BubbleApplyFormal()
    {
        if (_bubbleAtoms.Length == 0 || !int.TryParse(_bubbleFormal.Replace('−', '-').Replace("+", ""), NumberStyles.Integer, CultureInfo.InvariantCulture, out var f)) return;
        RunEdit(new { op = "charge", atoms = _bubbleAtoms, charge = (double)f });
        Reselect(_bubbleAtoms);
        OpenBubble();
    }

    /// <summary>Nudge: the atoms moved 0.1 Å along a world axis (x, y, z; sign ±1).</summary>
    public void BubbleNudge(int axis, int sign)
    {
        if (_bubbleAtoms.Length == 0) return;
        var by = new double[3];
        by[Math.Clamp(axis, 0, 2)] = 0.1 * Math.Sign(sign);
        TranslateAtoms(_bubbleAtoms, by);
        Reselect(_bubbleAtoms);
    }

    /// <summary>The atoms selected again after an edit (an edit clears the selection).</summary>
    private void Reselect(int[] atoms)
    {
        if (_doc == null) return;
        var n = (int)_doc.Summary().Atoms;
        var keep = atoms.Where(i => i < n).ToArray();
        if (keep.Length == 0) return;
        var r = JsonNode.Parse(_doc.Select(new JsonObject { ["mode"] = "indices", ["atoms"] = new JsonArray(keep.Select(i => (JsonNode)i).ToArray()), ["op"] = "replace" }.ToJsonString()))!;
        SelectedCount = (int)((double?)r["count"] ?? 0);
        RefreshSelBar();
        RenderRequested?.Invoke();
    }

    // ---------------------------------------------------------------- the orient puck

    /// <summary>The view along the selection's long axis ("along"), down onto its plane ("onto") or side-on ("side"),
    /// framed on it. The renderer turns by yaw about y then pitch about x: a direction u faces the viewer when
    /// yaw = atan2(−uₓ, u_z) and pitch = atan2(u_y, √(uₓ² + u_z²)).</summary>
    public void OrientView(string how)
    {
        if (_doc == null) return;
        var a = SelectionAtoms();
        if (a.Length < 3) { Status = "Select at least three atoms to turn the view to"; return; }
        try
        {
            var g = JsonNode.Parse(_doc.ProbeGeometry(new JsonObject { ["kind"] = "ellipsoid", ["atoms"] = new JsonArray(a.Select(i => (JsonNode)i).ToArray()) }.ToJsonString()))!;
            var axes = (g["axes"] as JsonArray)!.Select(v => (v as JsonArray)!.Select(x => (double)x!).ToArray()).ToArray();
            var u = how switch { "along" => axes[0], "onto" => axes[2], _ => axes[1] };
            var r = Math.Sqrt(u[0] * u[0] + u[2] * u[2]);
            Camera.Yaw = Math.Atan2(-u[0], u[2]);
            Camera.Pitch = Math.Atan2(u[1], r);
            FrameSelection();
            Status = how switch { "along" => "Looking along the selection's long axis", "onto" => "Looking down onto the selection's plane", _ => "The selection side-on" };
        }
        catch (Exception e) { Status = "Could not turn the view: " + e.Message; }
        RenderRequested?.Invoke();
    }
}
