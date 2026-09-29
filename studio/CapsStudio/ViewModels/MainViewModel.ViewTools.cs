using System.Collections.ObjectModel;
using System.Globalization;
using System.Text.Json.Nodes;
using CapsStudio.Interop;

namespace CapsStudio.ViewModels;

/// <summary>A measurement pinned to the view: its atoms and its value, kept up to date with the frame and edits.</summary>
public sealed class MonitorRow : ObservableObject
{
    public int[] Atoms { get; init; } = [];
    public string Ids { get; init; } = "";
    public string Kind => Atoms.Length switch { 2 => "d", 3 => "∠", _ => "φ" };
    private string _value = "";
    public string Value { get => _value; set => Set(ref _value, value); }
    public string Text => $"{Kind} {Ids}  {Value}";
}

/// <summary>A pinned measurement as the overlay draws it: a dashed line for a distance, the label at the atoms' middle.</summary>
public readonly record struct MonitorMark(double X0, double Y0, double X1, double Y1, bool Line, string Text);

/// <summary>The Studio view's lasso, translate and pin-monitor tools (the Main board's toolbar).</summary>
public sealed partial class MainViewModel
{
    // EditTool 4 lasso select, 5 translate
    public bool IsLassoTool => _editTool == 4;
    public bool IsMoveTool => _editTool == 5;

    /// <summary>The atoms whose centres fall inside the lasso (visible ones) become the selection, or join it.</summary>
    public void LassoSelect(IReadOnlyCollection<int> atoms, bool add)
    {
        if (_doc == null) return;
        if (atoms.Count == 0 && !add) { ClearDocSelection(); Status = "Lasso: nothing inside; the selection is cleared"; return; }
        var r = JsonNode.Parse(_doc.Select(new JsonObject
        {
            ["mode"] = "indices", ["atoms"] = new JsonArray(atoms.Select(a => (JsonNode)a).ToArray()), ["op"] = add ? "add" : "replace",
        }.ToJsonString()))!;
        if (r["ok"]?.GetValue<bool>() != true) { Status = "Lasso: " + (r["error"]?.GetValue<string>() ?? "cannot select"); return; }
        SelectedCount = (int)r["count"]!.GetValue<double>();
        SelectHud = "lasso";
        Status = $"Lasso: {atoms.Count:N0} atoms {(add ? "added" : "selected")} · {SelectedCount:N0} selected";
        RenderRequested?.Invoke();
    }

    /// <summary>What a drag with the move tool carries: the selection, else the molecule of the atom under the cursor.</summary>
    public int[] MoveSet(int hit)
    {
        if (_doc == null) return [];
        var sel = JsonNode.Parse(_doc.SelectionJson())!["indices"]!.AsArray().Select(x => (int)x!.GetValue<double>()).ToArray();
        if (sel.Length > 0) return sel;
        if (hit < 0) return [];
        var mol = _doc.Atom(hit).Mol;
        var n = (int)_doc.Summary().Atoms;
        return Enumerable.Range(0, n).Where(i => _doc.Atom(i).Mol == mol).ToArray();
    }

    /// <summary>Moves the atoms rigidly by the world vector (one undoable step).</summary>
    public void TranslateAtoms(int[] atoms, double[] by)
    {
        if (atoms.Length == 0 || by.Length != 3) return;
        RunEdit(new { op = "translate", atoms, by });
    }

    /// <summary>A screen drag (view points) as a move in the view plane, from the projection of the atoms: the map
    /// screen = A·world + b is fitted by least squares to the visible atoms' projections, and the smallest world move
    /// giving the drag is Aᵀ(AAᵀ)⁻¹ d (exact for the orthographic view).</summary>
    public static double[]? ScreenToWorld(IReadOnlyList<(double X, double Y, double Z, double Sx, double Sy)> pts, double dx, double dy)
    {
        if (pts.Count < 4) return null;
        // normal equations for [a1 a2 a3 b] · [x y z 1] = s, once for sx and once for sy
        var m = new double[4, 4];
        double[] vx = new double[4], vy = new double[4];
        foreach (var p in pts)
        {
            double[] w = [p.X, p.Y, p.Z, 1];
            for (var r = 0; r < 4; ++r)
            {
                for (var c = 0; c < 4; ++c) m[r, c] += w[r] * w[c];
                vx[r] += w[r] * p.Sx;
                vy[r] += w[r] * p.Sy;
            }
        }
        var ax = Solve4(m, vx);
        var ay = Solve4(m, vy);
        if (ax == null || ay == null) return null;
        // A = [ax0 ax1 ax2; ay0 ay1 ay2]; Δw = Aᵀ (A Aᵀ)⁻¹ d
        double g11 = ax[0] * ax[0] + ax[1] * ax[1] + ax[2] * ax[2], g22 = ay[0] * ay[0] + ay[1] * ay[1] + ay[2] * ay[2];
        double g12 = ax[0] * ay[0] + ax[1] * ay[1] + ax[2] * ay[2];
        var det = g11 * g22 - g12 * g12;
        if (Math.Abs(det) < 1e-12) return null;
        var u = (g22 * dx - g12 * dy) / det;
        var v = (-g12 * dx + g11 * dy) / det;
        return [ax[0] * u + ay[0] * v, ax[1] * u + ay[1] * v, ax[2] * u + ay[2] * v];
    }

    private static double[]? Solve4(double[,] a0, double[] b0)
    {
        var a = (double[,])a0.Clone();
        var b = (double[])b0.Clone();
        for (var c = 0; c < 4; ++c)
        {
            var piv = c;
            for (var r = c + 1; r < 4; ++r) if (Math.Abs(a[r, c]) > Math.Abs(a[piv, c])) piv = r;
            if (Math.Abs(a[piv, c]) < 1e-12) return null;
            if (piv != c)
            {
                for (var k = 0; k < 4; ++k) (a[c, k], a[piv, k]) = (a[piv, k], a[c, k]);
                (b[c], b[piv]) = (b[piv], b[c]);
            }
            for (var r = 0; r < 4; ++r)
            {
                if (r == c) continue;
                var f = a[r, c] / a[c, c];
                for (var k = c; k < 4; ++k) a[r, k] -= f * a[c, k];
                b[r] -= f * b[c];
            }
        }
        return [b[0] / a[0, 0], b[1] / a[1, 1], b[2] / a[2, 2], b[3] / a[3, 3]];
    }

    // ---------------------------------------------------------------- pinned monitors
    public ObservableCollection<MonitorRow> Monitors { get; } = new();
    public bool HasMonitors => Monitors.Count > 0;

    /// <summary>The distance, angle or dihedral of the 2–4 picked atoms, pinned: it follows the frames and edits.</summary>
    public void PinMeasurement()
    {
        if (_doc == null || _selection.Count < 2) { Status = "Pick 2–4 atoms (⇧ click) to measure, then pin the measurement"; return; }
        var atoms = _selection.Take(4).ToArray();
        if (Monitors.Any(m => m.Atoms.SequenceEqual(atoms))) { Status = "That measurement is pinned already"; return; }
        var row = new MonitorRow { Atoms = atoms, Ids = string.Join("–", atoms.Select(i => _doc.Atom(i).Id)) };
        Monitors.Add(row);
        Raise(nameof(HasMonitors)); Raise(nameof(MonitorsVisible));
        RefreshMonitors();
        Status = $"Pinned {row.Text.Trim()} · it updates with the frames";
        RenderRequested?.Invoke();
    }

    public void UnpinMonitor(MonitorRow m)
    {
        Monitors.Remove(m);
        Raise(nameof(HasMonitors)); Raise(nameof(MonitorsVisible));
        RenderRequested?.Invoke();
    }

    /// <summary>Recomputes each monitor on the shown frame; monitors whose atoms are gone (after a delete) are dropped.</summary>
    public void RefreshMonitors()
    {
        if (_doc == null || Monitors.Count == 0) return;
        var n = (int)_doc.Summary().Atoms;
        foreach (var m in Monitors.Where(m => m.Atoms.Any(a => a >= n)).ToList()) Monitors.Remove(m);
        Raise(nameof(HasMonitors)); Raise(nameof(MonitorsVisible));
        foreach (var m in Monitors)
        {
            var v = _doc.Measure(m.Atoms);
            m.Value = m.Atoms.Length == 2 ? v.ToString("0.000", CultureInfo.InvariantCulture) + " Å" : v.ToString("0.0", CultureInfo.InvariantCulture) + "°";
        }
    }

    /// <summary>The monitors as the overlay draws them, from the atoms' screen positions in this render.</summary>
    public List<MonitorMark> MonitorMarks(CapsCamera cam, CapsRenderOpts opt, double scaling)
    {
        var list = new List<MonitorMark>();
        if (_doc == null || Monitors.Count == 0) return list;
        RefreshMonitors();
        var p = _doc.ProjectAtoms(cam, opt, (int)_doc.Summary().Atoms);
        foreach (var m in Monitors)
        {
            if (m.Atoms.Any(a => 3 * a + 2 >= p.Length)) continue;
            var xs = m.Atoms.Select(a => p[3 * a] / scaling).ToArray();
            var ys = m.Atoms.Select(a => p[3 * a + 1] / scaling).ToArray();
            if (m.Atoms.Length == 2) list.Add(new MonitorMark(xs[0], ys[0], xs[1], ys[1], true, m.Value));
            else list.Add(new MonitorMark(xs.Average(), ys.Average(), 0, 0, false, $"{m.Kind} {m.Value}"));
        }
        return list;
    }
}
