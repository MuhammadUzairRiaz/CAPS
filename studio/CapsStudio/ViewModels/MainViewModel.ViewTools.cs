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
    // EditTool 4 lasso select, 5 translate, 6 rotate
    public bool IsLassoTool => _editTool == 4;
    public bool IsMoveTool => _editTool == 5;
    public bool IsRotateTool => _editTool == 6;

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

    /// <summary>Rotates the atoms rigidly about the axis through their centre (one undoable step).</summary>
    public void RotateAtoms(int[] atoms, double[] axis, double degrees)
    {
        if (atoms.Length == 0 || axis.Length != 3 || Math.Abs(degrees) < 1e-6) return;
        RunEdit(new { op = "rotate", atoms, axis, degrees });
    }

    /// <summary>A rotate-tool drag as an axis and an angle. Axis locked ('x', 'y', 'z'): that world axis, 0.5° per point
    /// of horizontal drag. Otherwise the view axis through the atoms' centre, by the angle the cursor swept around the
    /// centre's projection; its sense is checked on the fitted projection, so the atoms turn the way the drag goes.
    /// snap: whole multiples of 15°.</summary>
    public static (double[] Axis, double Degrees)? RotationFromDrag(IReadOnlyList<(double X, double Y, double Z, double Sx, double Sy)> pts, double[] centre,
                                                                   double pressX, double pressY, double x, double y, char axisLock, bool snap)
    {
        double deg;
        double[] axis;
        if (axisLock is 'x' or 'y' or 'z')
        {
            axis = [axisLock == 'x' ? 1 : 0, axisLock == 'y' ? 1 : 0, axisLock == 'z' ? 1 : 0];
            deg = (x - pressX) * 0.5;
        }
        else
        {
            if (Fit(pts) is not { } f) return null;
            var (ax, ay) = f;
            // the view axis: normal to the screen's two world gradients
            double[] n = [ax[1] * ay[2] - ax[2] * ay[1], ax[2] * ay[0] - ax[0] * ay[2], ax[0] * ay[1] - ax[1] * ay[0]];
            var nn = Math.Sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
            if (nn < 1e-12) return null;
            for (var k = 0; k < 3; ++k) n[k] /= nn;
            double Sx(double[] w) => ax[0] * w[0] + ax[1] * w[1] + ax[2] * w[2] + ax[3];
            double Sy(double[] w) => ay[0] * w[0] + ay[1] * w[1] + ay[2] * w[2] + ay[3];
            var cx = Sx(centre); var cy = Sy(centre);
            deg = (Math.Atan2(y - cy, x - cx) - Math.Atan2(pressY - cy, pressX - cx)) * 180 / Math.PI;
            if (deg > 180) deg -= 360;
            if (deg < -180) deg += 360;
            // the sense: a +10° turn of a vector across the view about n, seen on the screen
            var gx = Math.Sqrt(ax[0] * ax[0] + ax[1] * ax[1] + ax[2] * ax[2]);
            double[] t = [ax[0] / gx, ax[1] / gx, ax[2] / gx];
            var th = 10 * Math.PI / 180;
            double[] cxn = [n[1] * t[2] - n[2] * t[1], n[2] * t[0] - n[0] * t[2], n[0] * t[1] - n[1] * t[0]];
            double[] rt = [t[0] * Math.Cos(th) + cxn[0] * Math.Sin(th), t[1] * Math.Cos(th) + cxn[1] * Math.Sin(th), t[2] * Math.Cos(th) + cxn[2] * Math.Sin(th)];
            double[] o = [0, 0, 0];
            double a0 = Math.Atan2(Sy(t) - Sy(o), Sx(t) - Sx(o)), a1 = Math.Atan2(Sy(rt) - Sy(o), Sx(rt) - Sx(o));
            var seen = a1 - a0;
            if (seen > Math.PI) seen -= 2 * Math.PI;
            if (seen < -Math.PI) seen += 2 * Math.PI;
            if (seen < 0) for (var k = 0; k < 3; ++k) n[k] = -n[k];
            axis = n;
        }
        if (snap) deg = Math.Round(deg / 15) * 15;
        return (axis, deg);
    }

    // the screen map s = A·w + b fitted by least squares to the atoms' projections: (a_x, b_x), (a_y, b_y) as 4 numbers each
    private static (double[] Ax, double[] Ay)? Fit(IReadOnlyList<(double X, double Y, double Z, double Sx, double Sy)> pts)
    {
        if (pts.Count < 4) return null;
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
        return ax == null || ay == null ? null : (ax, ay);
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
