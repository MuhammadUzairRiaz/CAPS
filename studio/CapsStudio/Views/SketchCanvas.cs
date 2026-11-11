using System.Globalization;
using System.Text.Json;
using System.Text.Json.Nodes;
using Avalonia;
using Avalonia.Controls;
using Avalonia.Input;
using Avalonia.Media;
using CapsStudio.Interop;

namespace CapsStudio.Views;

/// <summary>The molecule sketch (design/boards/Sketch): skeletal drawing on a dot grid. Draw: click places an atom (or
/// a bonded one when clicked on an atom), drag grows a bond, click a bond cycles its order, click an atom with another
/// element changes it. Ring: a benzene ring on empty space, on an atom (spiro) or fused on a bond. Erase removes an
/// atom or a bond. Charge cycles 0, +1, −1. Every edit writes the SMILES (through the core) and raises Edited.</summary>
public sealed class SketchCanvas : Control
{
    public enum Tool { Draw, Ring, Erase, Charge, Bond, Cut, Move }

    private sealed class A
    {
        public int Z = 6, Charge, Isotope, HCount = -1, Chiral, Map;
        public bool Aromatic, Bracket;
        public List<int> Order = new();
        public double X, Y;
    }
    private sealed class B { public int A, B2, Order = 1, Dir; }

    private readonly List<A> _atoms = new();
    private readonly List<B> _bonds = new();
    public Tool Mode { get => _mode; set { _mode = value; _pending = -1; InvalidateVisual(); } }
    private Tool _mode = Tool.Draw;
    public int Element { get; set; } = 6;
    /// <summary>The Bond tool's order: 1, 2 or 3.</summary>
    public int BondOrder { get; set; } = 1;
    /// <summary>The Ring tool's size (3 to 8) and whether a six-ring is benzene (alternating double bonds) or cyclohexane.</summary>
    public int RingSize { get; set; } = 6;
    public bool RingAromatic { get; set; } = true;
    public bool CanUndo => _undo.Count > 0;
    public bool CanRedo => _redo.Count > 0;
    /// <summary>Raised when undo or redo becomes possible or impossible.</summary>
    public event Action? HistoryChanged;
    private readonly List<string> _undo = new(), _redo = new();
    private string? _before;                 // the sketch when the pointer went down, pushed when that gesture edits it
    private Dictionary<A, (int Sum, int Pi)>? _beforeSums; // each atom's bond-order sum (and aromatic pi share) when the pointer went down (bracket hydrogens follow it)
    private int _moveAtom = -1;
    private int _pending = -1;              // Bond tool: the first of two atoms clicked
    private bool _panning;
    /// <summary>Raised with the new SMILES after an edit (empty when the sketch is empty).</summary>
    public event Action<string>? Edited;

    private double _scale = 42;              // pixels per bond
    private Point _origin;                  // canvas position of model (0, 0)
    private int _hoverAtom = -1, _hoverBond = -1;
    private int _dragFrom = -1;
    private Point _dragAt, _pressAt;
    private bool _dragging;

    public SketchCanvas()
    {
        Focusable = true;
        ClipToBounds = true;
    }

    // ---------------------------------------------------------------- model in and out
    /// <summary>Loads a caps_smiles_depict graph (empty string clears) and fits it to the view.</summary>
    public void Load(string json)
    {
        if (_atoms.Count > 0) Remember(Snapshot());
        LoadGraph(json);
        Fit();
        InvalidateVisual();
    }

    private void LoadGraph(string json)
    {
        _atoms.Clear();
        _bonds.Clear();
        if (json.Length > 0)
        {
            using var js = JsonDocument.Parse(json);
            var r = js.RootElement;
            if (r.TryGetProperty("atoms", out var atoms))
            {
                foreach (var a in atoms.EnumerateArray())
                    _atoms.Add(new A
                    {
                        Z = a.GetProperty("z").GetInt32(), Charge = a.GetProperty("charge").GetInt32(), Isotope = a.GetProperty("isotope").GetInt32(),
                        HCount = a.GetProperty("hcount").GetInt32(), Chiral = a.GetProperty("chiral").GetInt32(), Map = a.GetProperty("map").GetInt32(),
                        Aromatic = a.GetProperty("aromatic").GetBoolean(), Bracket = a.GetProperty("bracket").GetBoolean(),
                        Order = a.GetProperty("order").EnumerateArray().Select(x => x.GetInt32()).ToList(),
                        X = a.GetProperty("x").GetDouble(), Y = -a.GetProperty("y").GetDouble(),
                    });
                foreach (var b in r.GetProperty("bonds").EnumerateArray())
                    _bonds.Add(new B { A = b.GetProperty("a").GetInt32(), B2 = b.GetProperty("b").GetInt32(), Order = b.GetProperty("order").GetInt32(), Dir = b.GetProperty("dir").GetInt32() });
            }
        }
    }

    private void Fit()
    {
        if (Bounds.Width < 10) { _pendingFit = true; return; }
        _pendingFit = false;
        if (_atoms.Count == 0) { _scale = 42; _origin = new Point(Bounds.Width / 2, Bounds.Height / 2); return; }
        double x0 = _atoms.Min(a => a.X), x1 = _atoms.Max(a => a.X), y0 = _atoms.Min(a => a.Y), y1 = _atoms.Max(a => a.Y);
        var w = Math.Max(1, x1 - x0) + 2;
        var h = Math.Max(1, y1 - y0) + 2;
        _scale = Math.Clamp(Math.Min(Bounds.Width / w, Bounds.Height / h), 16, 48);
        _origin = new Point(Bounds.Width / 2 - (x0 + x1) / 2 * _scale, Bounds.Height / 2 - (y0 + y1) / 2 * _scale);
    }
    private bool _pendingFit = true;

    protected override void OnSizeChanged(SizeChangedEventArgs e)
    {
        base.OnSizeChanged(e);
        if (_pendingFit) Fit();
        InvalidateVisual();
    }

    private string Graph()
    {
        var atoms = new JsonArray();
        for (var i = 0; i < _atoms.Count; i++)
        {
            var a = _atoms[i];
            // an element outside the organic subset (Si, Se, a metal …) is written in brackets with its hydrogens counted,
            // as [SiH4] or [SeH2]: a bare [Se] would mean no hydrogens at all
            var bracket = a.Bracket || !(a.Z is 5 or 6 or 7 or 8 or 9 or 15 or 16 or 17 or 35 or 53);
            var hcount = a.Bracket ? a.HCount : bracket ? ImplicitH(i) : a.HCount;
            atoms.Add(new JsonObject
            {
                ["z"] = a.Z, ["charge"] = a.Charge, ["isotope"] = a.Isotope, ["hcount"] = hcount, ["aromatic"] = a.Aromatic,
                ["bracket"] = bracket, ["chiral"] = a.Chiral, ["map"] = a.Map, ["order"] = new JsonArray(a.Order.Select(v => (JsonNode)v).ToArray()),
            });
        }
        var bonds = new JsonArray();
        foreach (var b in _bonds) bonds.Add(new JsonObject { ["a"] = b.A, ["b"] = b.B2, ["order"] = b.Order, ["dir"] = b.Dir });
        return new JsonObject { ["atoms"] = atoms, ["bonds"] = bonds }.ToJsonString();
    }

    private static string Sym(int z) => z >= 1 && z <= ViewModels.MainViewModel.Syms.Length ? ViewModels.MainViewModel.Syms[z - 1] : "*";

    private string Snapshot() => JsonSerializer.Serialize(new
    {
        atoms = _atoms.Select(a => new[] { a.Z, a.Charge, a.Isotope, a.HCount, a.Chiral, a.Map, a.Aromatic ? 1 : 0, a.Bracket ? 1 : 0 }).ToArray(),
        xy = _atoms.Select(a => new[] { a.X, a.Y }).ToArray(),
        order = _atoms.Select(a => a.Order.ToArray()).ToArray(),
        bonds = _bonds.Select(b => new[] { b.A, b.B2, b.Order, b.Dir }).ToArray(),
    });

    private void Restore(string snap)
    {
        using var js = JsonDocument.Parse(snap);
        var r = js.RootElement;
        _atoms.Clear();
        _bonds.Clear();
        var xy = r.GetProperty("xy").EnumerateArray().ToList();
        var ord = r.GetProperty("order").EnumerateArray().ToList();
        var k = 0;
        foreach (var a in r.GetProperty("atoms").EnumerateArray())
        {
            var v = a.EnumerateArray().Select(x => x.GetInt32()).ToArray();
            var p = xy[k].EnumerateArray().Select(x => x.GetDouble()).ToArray();
            _atoms.Add(new A { Z = v[0], Charge = v[1], Isotope = v[2], HCount = v[3], Chiral = v[4], Map = v[5], Aromatic = v[6] != 0, Bracket = v[7] != 0,
                               X = p[0], Y = p[1], Order = ord[k].EnumerateArray().Select(x => x.GetInt32()).ToList() });
            k++;
        }
        foreach (var b in r.GetProperty("bonds").EnumerateArray())
        {
            var v = b.EnumerateArray().Select(x => x.GetInt32()).ToArray();
            _bonds.Add(new B { A = v[0], B2 = v[1], Order = v[2], Dir = v[3] });
        }
    }

    private void Remember(string snap)
    {
        _undo.Add(snap);
        if (_undo.Count > 200) _undo.RemoveAt(0);
        _redo.Clear();
        HistoryChanged?.Invoke();
    }

    /// <summary>Undo (or redo) the last sketch edit.</summary>
    public void Undo(bool redo = false)
    {
        var from = redo ? _redo : _undo;
        if (from.Count == 0) return;
        (redo ? _undo : _redo).Add(Snapshot());
        Restore(from[^1]);
        from.RemoveAt(from.Count - 1);
        _hoverAtom = _hoverBond = -1;
        Emit();
        HistoryChanged?.Invoke();
    }

    /// <summary>Clears the sketch (one step back with Undo).</summary>
    public void Clear()
    {
        if (_atoms.Count == 0) return;
        Remember(Snapshot());
        _atoms.Clear();
        _bonds.Clear();
        Fit();
        Emit();
    }

    /// <summary>Lays the drawing out again (the core's 2D depiction of its SMILES), keeping it as one undo step.</summary>
    public void Tidy(string depictJson)
    {
        if (_atoms.Count == 0 || depictJson.Length == 0) return;
        Remember(Snapshot());
        LoadGraph(depictJson);
        Fit();
        InvalidateVisual();
    }

    private void Emit()
    {
        string smiles;
        try { smiles = _atoms.Count == 0 ? "" : CapsDocument.SmilesWrite(Graph()); }
        catch { smiles = ""; }
        InvalidateVisual();
        Edited?.Invoke(smiles);
    }

    /// <summary>Raised with a short message when an edit was refused (the status line shows it).</summary>
    public event Action<string>? Notice;

    /// <summary>The most bonds an uncharged atom takes (its highest common valence); 0: no limit kept (metals, others).</summary>
    private static int MaxValence(int z) => z switch
    {
        1 or 9 or 17 or 35 or 53 => 1, 8 => 2, 5 or 7 => 3, 6 or 14 or 32 or 50 => 4, 15 or 33 => 5, 16 or 34 or 52 => 6, _ => 0,
    };

    /// <summary>An atom of `touched` holding more bonds than its valence allows, or −1.</summary>
    private int OverValence(int[] touched)
    {
        foreach (var i in touched)
        {
            if (i < 0 || i >= _atoms.Count) continue;
            var a = _atoms[i];
            var max = MaxValence(a.Z);
            if (max == 0 || a.Charge != 0 || a.Aromatic) continue;
            var sum = _bonds.Where(b => b.A == i || b.B2 == i).Sum(b => b.Order == 4 ? 1 : b.Order);
            if (sum > max) return i;
        }
        return -1;
    }

    private void Changed(params int[] touched)
    {
        // an edit that gives an atom more bonds than its valence (a double bond on a ring carbon that has one already, a
        // fifth bond on carbon) is taken back, and the status line says why
        if (OverValence(touched) is var bad and >= 0 && _before != null)
        {
            var sym = Sym(_atoms[bad].Z);
            var max = MaxValence(_atoms[bad].Z);
            Restore(_before);
            _before = null;
            _pending = -1;
            InvalidateVisual();
            Notice?.Invoke($"Not drawn: {sym} takes {max} bonds and that would give it more. Remove a bond from it first, or use a lower bond order.");
            return;
        }
        if (_before != null) { Remember(_before); _before = null; }
        // a bracket atom carries its hydrogens as a count: each bond order it gains takes one hydrogen and each it loses
        // gives one back, so [C@H] with a new bond becomes [C] and [O-] cut from its carbon becomes [OH-]
        if (_beforeSums != null)
        {
            for (var i = 0; i < _atoms.Count; i++)
            {
                var a = _atoms[i];
                if (!a.Bracket || a.HCount < 0 || !_beforeSums.TryGetValue(a, out var was)) continue;
                // an atom still aromatic keeps its pi share; one written out as Kekulé counts the double bond it now has
                var now = BondSum(i) + (a.Aromatic && HasAromaticBond(i) ? was.Pi : 0);
                a.HCount = Math.Max(0, a.HCount - (now - was.Sum));
            }
            _beforeSums = null;
        }
        // an edit next to a stereocentre or a marked double bond drops that stereo (it would no longer mean the same)
        foreach (var i in touched)
            if (i >= 0 && i < _atoms.Count)
            {
                _atoms[i].Chiral = 0;
                foreach (var b in _bonds) if (b.A == i || b.B2 == i) b.Dir = 0;
            }
        Emit();
    }

    // ---------------------------------------------------------------- geometry helpers
    private Point ToScreen(double x, double y) => new(_origin.X + x * _scale, _origin.Y + y * _scale);
    private (double X, double Y) ToModel(Point p) => ((p.X - _origin.X) / _scale, (p.Y - _origin.Y) / _scale);

    private int AtomAt(Point p)
    {
        var best = -1;
        var bd = 0.35 * _scale;
        for (var i = 0; i < _atoms.Count; i++)
        {
            var d = Dist(ToScreen(_atoms[i].X, _atoms[i].Y), p);
            if (d < bd) { bd = d; best = i; }
        }
        return best;
    }

    private int BondAt(Point p)
    {
        for (var k = 0; k < _bonds.Count; k++)
        {
            var a = ToScreen(_atoms[_bonds[k].A].X, _atoms[_bonds[k].A].Y);
            var b = ToScreen(_atoms[_bonds[k].B2].X, _atoms[_bonds[k].B2].Y);
            var ab = b - a;
            var len2 = ab.X * ab.X + ab.Y * ab.Y;
            if (len2 < 1) continue;
            var t = Math.Clamp(((p.X - a.X) * ab.X + (p.Y - a.Y) * ab.Y) / len2, 0, 1);
            if (t < 0.15 || t > 0.85) continue;
            var q = new Point(a.X + t * ab.X, a.Y + t * ab.Y);
            if (Dist(q, p) < 0.18 * _scale) return k;
        }
        return -1;
    }

    /// <summary>What a Ring click lands on: an atom (spiro), a bond (fused), or neither (a free ring). Near the drawing a
    /// click always means a bond: the nearest one, by its middle, within 1.8 bond lengths, so a ring is never dropped on
    /// top of the ones drawn (the bond itself is a thin target when the drawing is zoomed out).</summary>
    private (int Atom, int Bond) RingTarget(Point p)
    {
        var at = AtomAt(p);
        if (at >= 0) return (at, -1);
        bool Edge(int k) => Neighbours(_bonds[k].A).Count <= 2 && Neighbours(_bonds[k].B2).Count <= 2;
        var bond = BondAt(p);
        if (bond >= 0 && Edge(bond)) return (-1, bond);
        var best = -1;
        var bd = 1.8 * _scale;
        for (var k = 0; k < _bonds.Count; k++)
        {
            // an edge a ring can fuse onto: both atoms with room (a shared bond inside a ring system has none)
            if (!Edge(k)) continue;
            var d = Dist(BondMid(k), p);
            if (d < bd) { bd = d; best = k; }
        }
        if (best >= 0) return (-1, best);
        // still near an atom (a lone one, or the end of a chain): spiro on it rather than a ring over it
        var near = -1;
        var nd = 1.0 * _scale;
        for (var i = 0; i < _atoms.Count; i++)
        {
            if (Neighbours(i).Count > 2) continue;   // a spiro ring needs two free bonds on the atom
            var d = Dist(ToScreen(_atoms[i].X, _atoms[i].Y), p);
            if (d < nd) { nd = d; near = i; }
        }
        return (near, -1);
    }

    private Point BondMid(int k)
    {
        var a = ToScreen(_atoms[_bonds[k].A].X, _atoms[_bonds[k].A].Y);
        var b = ToScreen(_atoms[_bonds[k].B2].X, _atoms[_bonds[k].B2].Y);
        return new Point((a.X + b.X) / 2, (a.Y + b.Y) / 2);
    }

    private static double Dist(Point a, Point b) => Math.Sqrt((a.X - b.X) * (a.X - b.X) + (a.Y - b.Y) * (a.Y - b.Y));

    private List<int> Neighbours(int i)
    {
        var n = new List<int>();
        foreach (var b in _bonds) { if (b.A == i) n.Add(b.B2); else if (b.B2 == i) n.Add(b.A); }
        return n;
    }

    /// <summary>The free direction for a new bond on atom i (120° from its neighbours, or straight on for one).</summary>
    private (double X, double Y) FreeDirection(int i)
    {
        var a = _atoms[i];
        var n = Neighbours(i);
        if (n.Count == 0) return (Math.Cos(-Math.PI / 6), Math.Sin(-Math.PI / 6));
        var angles = n.Select(k => Math.Atan2(_atoms[k].Y - a.Y, _atoms[k].X - a.X)).OrderBy(t => t).ToList();
        if (n.Count == 1)
        {
            // zig-zag: 120° from the one neighbour, on the side away from its own other neighbours
            var t0 = angles[0];
            var c1 = t0 + 2 * Math.PI / 3;
            var c2 = t0 - 2 * Math.PI / 3;
            var nn = Neighbours(n[0]).Where(k => k != i).ToList();
            if (nn.Count > 0)
            {
                var far = nn.Select(k => Math.Atan2(_atoms[k].Y - a.Y, _atoms[k].X - a.X)).First();
                return Gap(c1, far) > Gap(c2, far) ? (Math.Cos(c1), Math.Sin(c1)) : (Math.Cos(c2), Math.Sin(c2));
            }
            return (Math.Cos(c2), Math.Sin(c2));
        }
        // the middle of the widest gap between neighbours
        double best = 0, dir = 0;
        for (var k = 0; k < angles.Count; k++)
        {
            var t1 = angles[k];
            var t2 = k + 1 < angles.Count ? angles[k + 1] : angles[0] + 2 * Math.PI;
            if (t2 - t1 > best) { best = t2 - t1; dir = (t1 + t2) / 2; }
        }
        return (Math.Cos(dir), Math.Sin(dir));
    }

    private static double Gap(double a, double b)
    {
        var d = Math.Abs(a - b) % (2 * Math.PI);
        return d > Math.PI ? 2 * Math.PI - d : d;
    }

    private int AddAtom(double x, double y)
    {
        _atoms.Add(new A { Z = Element, X = x, Y = y });
        return _atoms.Count - 1;
    }

    private void AddBond(int a, int b, int order = 1)
    {
        if (a == b || _bonds.Any(x => (x.A == a && x.B2 == b) || (x.A == b && x.B2 == a))) return;
        _bonds.Add(new B { A = a, B2 = b, Order = order });
        foreach (var k in new[] { a, b }) _atoms[k].Order.Clear();
    }

    private void RemoveAtom(int i)
    {
        _bonds.RemoveAll(b => b.A == i || b.B2 == i);
        _atoms.RemoveAt(i);
        foreach (var b in _bonds) { if (b.A > i) b.A--; if (b.B2 > i) b.B2--; }
        foreach (var a in _atoms) { a.Order.Clear(); a.Chiral = 0; }
    }

    /// <summary>A ring of RingSize atoms through atom `at` (spiro), fused on bond `bond`, or free at (x, y). A six-ring with
    /// RingAromatic on is benzene: Kekulé double bonds, shared with the ring it fuses onto (naphthalene and beyond).</summary>
    private void AddRing(double x, double y, int at, int bond)
    {
        var n = Math.Clamp(RingSize, 3, 8);
        var step0 = 2 * Math.PI / n;
        var apothem = 0.5 / Math.Tan(Math.PI / n);    // bond midpoint to centre, unit bonds
        var radius = 0.5 / Math.Sin(Math.PI / n);     // atom to centre
        var ids = new List<int>();
        if (bond >= 0 && (_atoms[_bonds[bond].A].Aromatic || _atoms[_bonds[bond].B2].Aromatic)) Dearomatise(_bonds[bond].A, _bonds[bond].B2);
        if (at >= 0 && _atoms[at].Aromatic) Dearomatise(at);
        if (bond >= 0)
        {
            var b = _bonds[bond];
            var pa = _atoms[b.A];
            var pb = _atoms[b.B2];
            // the ring goes on the side with fewer atoms
            double mx = (pa.X + pb.X) / 2, my = (pa.Y + pb.Y) / 2, nx = -(pb.Y - pa.Y), ny = pb.X - pa.X;
            var side = _atoms.Sum(a => Math.Sign((a.X - mx) * nx + (a.Y - my) * ny));
            if (side > 0) { nx = -nx; ny = -ny; }
            var len = Math.Sqrt(nx * nx + ny * ny);
            var bl = Math.Sqrt((pb.X - pa.X) * (pb.X - pa.X) + (pb.Y - pa.Y) * (pb.Y - pa.Y));
            var cx = mx + nx / len * apothem * bl;
            var cy = my + ny / len * apothem * bl;
            var r = radius * bl;
            var t0 = Math.Atan2(pa.Y - cy, pa.X - cx);
            var t1 = Math.Atan2(pb.Y - cy, pb.X - cx);
            var step = Gap(t0 + step0, t1) < 0.1 ? step0 : -step0;
            ids.Add(b.A);
            ids.Add(b.B2);
            for (var k = 2; k < n; k++) ids.Add(AddAtomAt(cx + r * Math.Cos(t0 + k * step), cy + r * Math.Sin(t0 + k * step)));
        }
        else if (at >= 0)
        {
            var (dx, dy) = FreeDirection(at);
            var cx = _atoms[at].X + dx * radius;
            var cy = _atoms[at].Y + dy * radius;
            var t0 = Math.Atan2(_atoms[at].Y - cy, _atoms[at].X - cx);
            ids.Add(at);
            for (var k = 1; k < n; k++) ids.Add(AddAtomAt(cx + radius * Math.Cos(t0 + k * step0), cy + radius * Math.Sin(t0 + k * step0)));
        }
        else
        {
            var tStart = Math.PI / 2 + Math.PI / n;      // a flat bottom edge
            for (var k = 0; k < n; k++) ids.Add(AddAtomAt(x + radius * Math.Cos(tStart + k * step0), y + radius * Math.Sin(tStart + k * step0)));
        }
        for (var k = 0; k < n; k++)
        {
            var a = ids[k];
            var b = ids[(k + 1) % n];
            if (!_bonds.Any(q => (q.A == a && q.B2 == b) || (q.A == b && q.B2 == a))) AddBond(a, b, 1);
        }
        if (n == 6 && RingAromatic) KekuleFill(ids);
        Changed(ids.ToArray());
    }

    /// <summary>Double bonds for a new ring: every ring atom that can still take a double bond gets one. A free atom of the
    /// new ring looks for an alternating path (single, double, single, …) through the ring bonds to another free atom and the
    /// path flips, so a ring fused on benzene shifts the old ring's double bonds where it must (naphthalene, anthracene,
    /// phenanthrene come out fully Kekulé, with one H on each CH). An atom with no partner left (a spiro centre, a ring fused
    /// on a saturated bond) keeps its hydrogens.</summary>
    private void KekuleFill(List<int> ids)
    {
        var ringBonds = new HashSet<B>();
        foreach (var r in SmallRings())
            foreach (var b in _bonds)
                if (r.Contains(b.A) && r.Contains(b.B2)) ringBonds.Add(b);
        bool CanTake(int i) => _atoms[i].Z is 6 or 7 && !HasDouble(i) && ImplicitH(i) > 0 && _atoms[i].Charge == 0;
        B? DoubleOf(int i) => _bonds.FirstOrDefault(b => (b.A == i || b.B2 == i) && b.Order == 2 && ringBonds.Contains(b));
        int Other(B b, int i) => b.A == i ? b.B2 : b.A;
        foreach (var start in ids)
        {
            if (!CanTake(start)) continue;
            var seen = new HashSet<int> { start };
            var path = new List<B>();
            bool Search(int u)
            {
                foreach (var b in _bonds.Where(b => (b.A == u || b.B2 == u) && b.Order == 1 && ringBonds.Contains(b)))
                {
                    var v = Other(b, u);
                    if (!seen.Add(v)) continue;
                    path.Add(b);
                    if (CanTake(v)) return true;
                    if (DoubleOf(v) is { } d && Other(d, v) is var w && seen.Add(w))
                    {
                        path.Add(d);
                        if (Search(w)) return true;
                        path.RemoveAt(path.Count - 1);
                    }
                    path.RemoveAt(path.Count - 1);
                }
                return false;
            }
            if (Search(start))
                foreach (var b in path) b.Order = b.Order == 1 ? 2 : 1;
        }
    }

    private bool HasDouble(int i) => _bonds.Any(b => (b.A == i || b.B2 == i) && (b.Order == 2 || b.Order == 4));

    /// <summary>An atom at (x, y), or the existing one there (rings fuse onto what is drawn).</summary>
    private int AddAtomAt(double x, double y, int z = 6)
    {
        for (var i = 0; i < _atoms.Count; i++)
            if (Math.Abs(_atoms[i].X - x) < 0.2 && Math.Abs(_atoms[i].Y - y) < 0.2) return i;
        _atoms.Add(new A { Z = z, X = x, Y = y });
        return _atoms.Count - 1;
    }

    // ---------------------------------------------------------------- for the self-test (it clicks these points with the pointer)
    internal int AtomCount => _atoms.Count;
    internal int BondCount => _bonds.Count;
    internal Point AtomPoint(int i) => ToScreen(_atoms[i].X, _atoms[i].Y);
    internal Point BondPoint(int k) { var a = AtomPoint(_bonds[k].A); var b = AtomPoint(_bonds[k].B2); return new Point((a.X + b.X) / 2, (a.Y + b.Y) / 2); }
    /// <summary>A bond of the outer edge (a ring bond whose atoms have only two neighbours each), for fusing the next ring.</summary>
    internal int OuterBond(double towardsX) =>
        Enumerable.Range(0, _bonds.Count).Where(k => Neighbours(_bonds[k].A).Count == 2 && Neighbours(_bonds[k].B2).Count == 2)
                  .OrderByDescending(k => (_atoms[_bonds[k].A].X + _atoms[_bonds[k].B2].X) * towardsX).FirstOrDefault(-1);
    internal int Neighbours2(int i) => Neighbours(i).Count;
    internal Point Centre() { if (_atoms.Count == 0) return default; return ToScreen(_atoms.Average(a => a.X), _atoms.Average(a => a.Y)); }
    internal double BondPixels => _scale;
    /// <summary>The closest two atoms (in bond lengths): below about 0.5 two rings were drawn over each other.</summary>
    internal double ClosestPair()
    {
        var best = double.MaxValue;
        for (var i = 0; i < _atoms.Count; i++)
            for (var j = i + 1; j < _atoms.Count; j++)
                best = Math.Min(best, Math.Sqrt((_atoms[i].X - _atoms[j].X) * (_atoms[i].X - _atoms[j].X) + (_atoms[i].Y - _atoms[j].Y) * (_atoms[i].Y - _atoms[j].Y)));
        return best;
    }
    internal int MaxHydrogens => _atoms.Count == 0 ? 0 : Enumerable.Range(0, _atoms.Count).Max(ImplicitH);

    // ---------------------------------------------------------------- input
    protected override void OnKeyDown(KeyEventArgs e)
    {
        base.OnKeyDown(e);
        var cmd = OperatingSystem.IsMacOS() ? KeyModifiers.Meta : KeyModifiers.Control;
        if (e.Key == Key.Z && e.KeyModifiers == cmd) { Undo(); e.Handled = true; }
        else if ((e.Key == Key.Z && e.KeyModifiers == (cmd | KeyModifiers.Shift)) || (e.Key == Key.Y && e.KeyModifiers == cmd)) { Undo(redo: true); e.Handled = true; }
        else if (e.Key == Key.Escape && _pending >= 0) { _pending = -1; InvalidateVisual(); e.Handled = true; }
    }

    protected override void OnPointerPressed(PointerPressedEventArgs e)
    {
        base.OnPointerPressed(e);
        Focus();
        var p = e.GetPosition(this);
        _pressAt = p;
        _dragAt = p;
        _before = Snapshot();
        _beforeSums = new Dictionary<A, (int Sum, int Pi)>();
        for (var i = 0; i < _atoms.Count; i++) { var pi = AromaticPi(i); _beforeSums[_atoms[i]] = (BondSum(i) + pi, pi); }
        _dragFrom = Mode is Tool.Draw or Tool.Bond ? AtomAt(p) : -1;
        _moveAtom = Mode == Tool.Move ? AtomAt(p) : -1;
        _panning = Mode == Tool.Move && _moveAtom < 0;
        _dragging = false;
        e.Pointer.Capture(this);
    }

    protected override void OnPointerMoved(PointerEventArgs e)
    {
        base.OnPointerMoved(e);
        var p = e.GetPosition(this);
        if (Mode == Tool.Move && e.GetCurrentPoint(this).Properties.IsLeftButtonPressed && (_moveAtom >= 0 || _panning))
        {
            var (dx, dy) = ((p.X - _dragAt.X) / _scale, (p.Y - _dragAt.Y) / _scale);
            if (_moveAtom >= 0) { _atoms[_moveAtom].X += dx; _atoms[_moveAtom].Y += dy; }
            else _origin = new Point(_origin.X + p.X - _dragAt.X, _origin.Y + p.Y - _dragAt.Y);
            _dragAt = p;
            if (Dist(p, _pressAt) > 3) _dragging = true;
            InvalidateVisual();
            return;
        }
        if (_dragFrom >= 0 && e.GetCurrentPoint(this).Properties.IsLeftButtonPressed)
        {
            _dragAt = p;
            if (Dist(p, _pressAt) > 6) _dragging = true;
            InvalidateVisual();
            return;
        }
        var (ha, hb) = Mode == Tool.Ring ? RingTarget(p) : (AtomAt(p), -1);
        if (Mode != Tool.Ring && ha < 0) hb = BondAt(p);
        if (ha != _hoverAtom || hb != _hoverBond) { _hoverAtom = ha; _hoverBond = hb; InvalidateVisual(); }
    }

    protected override void OnPointerExited(PointerEventArgs e)
    {
        base.OnPointerExited(e);
        _hoverAtom = _hoverBond = -1;
        InvalidateVisual();
    }

    protected override void OnPointerReleased(PointerReleasedEventArgs e)
    {
        base.OnPointerReleased(e);
        e.Pointer.Capture(null);
        var p = e.GetPosition(this);
        var from = _dragFrom;
        _dragFrom = -1;
        var at = AtomAt(p);
        var bond = at < 0 ? BondAt(p) : -1;
        // the click is where the pointer went down: a small slip while pressing (a trackpad tap) does not move it
        if (!_dragging && Mode is not (Tool.Move) && Dist(p, _pressAt) < 8)
        {
            p = _pressAt;
            at = AtomAt(p);
            bond = at < 0 ? BondAt(p) : -1;
        }
        switch (Mode)
        {
            case Tool.Move:
                // an atom moved: one undo step (the SMILES is the same, so nothing else changes)
                if (_moveAtom >= 0 && _dragging && _before != null) { Remember(_before); _before = null; }
                _moveAtom = -1;
                _panning = false;
                break;
            case Tool.Bond:
                {
                    var order = Math.Clamp(BondOrder, 1, 3);
                    if (from >= 0 && _dragging)
                    {
                        if (at >= 0 && at != from)
                        {
                            var old = _bonds.FirstOrDefault(x => (x.A == from && x.B2 == at) || (x.A == at && x.B2 == from));
                            if (old != null) SetOrder(old, order); else AddBond(from, at, order);
                            Changed(from, at);
                        }
                        else if (at < 0)
                        {
                            var (mx, my) = ToModel(p);
                            var t = Math.Atan2(my - _atoms[from].Y, mx - _atoms[from].X);
                            t = Math.Round(t / (Math.PI / 6)) * (Math.PI / 6);
                            var k = AddAtomAt(_atoms[from].X + Math.Cos(t), _atoms[from].Y + Math.Sin(t), Element);
                            AddBond(from, k, order);
                            Changed(from, k);
                        }
                    }
                    else if (at >= 0 && _pending >= 0 && _pending != at && _pending < _atoms.Count)
                    {
                        // the second of two clicked atoms: bond them (as the Studio's Draw bond does)
                        var first = _pending;
                        _pending = -1;
                        var old = _bonds.FirstOrDefault(x => (x.A == first && x.B2 == at) || (x.A == at && x.B2 == first));
                        if (old != null) SetOrder(old, order); else AddBond(first, at, order);
                        Changed(first, at);
                    }
                    else if (at >= 0 && _pending == at)
                    {
                        // the same atom twice: a new atom bonded to it
                        _pending = -1;
                        var (dx, dy) = FreeDirection(at);
                        var k = AddAtom(_atoms[at].X + dx, _atoms[at].Y + dy);
                        AddBond(at, k, order);
                        Changed(at, k);
                    }
                    else if (at >= 0) _pending = at;
                    else if (bond >= 0) { _pending = -1; SetOrder(_bonds[bond], order); Changed(_bonds[bond].A, _bonds[bond].B2); }
                    else if (_pending >= 0) _pending = -1;
                    else
                    {
                        // empty space: a two-atom fragment with this bond
                        var (mx, my) = ToModel(p);
                        var a0 = AddAtom(mx - 0.5, my);
                        var a1 = AddAtom(mx + 0.5, my);
                        AddBond(a0, a1, order);
                        Changed(a0, a1);
                    }
                }
                break;
            case Tool.Cut:
                // cut a bond and keep both atoms (Erase removes atoms too)
                if (bond >= 0) { var b = _bonds[bond]; _bonds.RemoveAt(bond); _atoms[b.A].Order.Clear(); _atoms[b.B2].Order.Clear(); Changed(b.A, b.B2); }
                break;
            case Tool.Draw:
                if (from >= 0 && _dragging)
                {
                    if (at >= 0 && at != from) { AddBond(from, at); Changed(from, at); }
                    else if (at < 0)
                    {
                        // snap to 30° steps, one bond long
                        var (mx, my) = ToModel(p);
                        var t = Math.Atan2(my - _atoms[from].Y, mx - _atoms[from].X);
                        t = Math.Round(t / (Math.PI / 6)) * (Math.PI / 6);
                        var k = AddAtom(_atoms[from].X + Math.Cos(t), _atoms[from].Y + Math.Sin(t));
                        AddBond(from, k);
                        Changed(from, k);
                    }
                }
                else if (at >= 0)
                {
                    if (_atoms[at].Z != Element)
                    {
                        _atoms[at].Z = Element;
                        _atoms[at].Aromatic = _atoms[at].Aromatic && Element is 6 or 7 or 8 or 16;
                        _atoms[at].Bracket = false;
                        _atoms[at].HCount = -1;
                        Changed(at);
                    }
                    else
                    {
                        var (dx, dy) = FreeDirection(at);
                        var k = AddAtom(_atoms[at].X + dx, _atoms[at].Y + dy);
                        AddBond(at, k);
                        Changed(at, k);
                    }
                }
                else if (bond >= 0)
                {
                    var b = _bonds[bond];
                    static int Next(int o) => o switch { 1 => 2, 2 => 3, 4 => 2, _ => 1 };
                    var was = b.Order;
                    var first = Next(was);
                    // a step the atoms have no room for (a triple bond on a carbon that has two other bonds) goes on round
                    // the cycle, so the bond falls back to single instead of sticking at double
                    b.Order = first;
                    while (b.Order != 1 && OverValence([b.A, b.B2]) >= 0) b.Order = Next(b.Order);
                    // no order fits but the one it had: try the first step, which Changed takes back and says why
                    if (b.Order == was) b.Order = first;
                    if (b.Order != 4 && (_atoms[b.A].Aromatic || _atoms[b.B2].Aromatic)) Dearomatise(b.A, b.B2);
                    Changed(b.A, b.B2);
                }
                else
                {
                    var (mx, my) = ToModel(p);
                    var k = AddAtom(mx, my);
                    Changed(k);
                }
                break;
            case Tool.Ring:
                {
                    var (mx, my) = ToModel(p);
                    var (ra, rb) = RingTarget(p);
                    AddRing(mx, my, ra, rb);
                }
                break;
            case Tool.Erase:
                if (at >= 0) { var n = Neighbours(at); RemoveAtom(at); Changed(); }
                else if (bond >= 0) { var b = _bonds[bond]; _bonds.RemoveAt(bond); _atoms[b.A].Order.Clear(); _atoms[b.B2].Order.Clear(); Changed(b.A, b.B2); }
                break;
            case Tool.Charge:
                if (at >= 0)
                {
                    var a = _atoms[at];
                    a.Charge = a.Charge switch { 0 => 1, 1 => -1, _ => 0 };
                    a.Bracket = a.Charge != 0 || a.Isotope != 0;
                    // the hydrogens follow the charge: a charged atom takes the valence of the element it is isoelectronic
                    // with (N+ as C: [NH3+]; O- as F: [O-]; C+ as B: [CH3+]; C- as N: [CH3-])
                    a.HCount = a.Bracket ? ValenceH(at, a.Z - a.Charge) : -1;
                    Changed(at);
                }
                break;
        }
        _dragging = false;
        _before = null;
        _beforeSums = null;
        _hoverAtom = AtomAt(p);
        InvalidateVisual();
    }

    /// <summary>A bond set to an order; an aromatic ring it belongs to is written out as Kekulé first.</summary>
    private void SetOrder(B b, int order)
    {
        if (_atoms[b.A].Aromatic || _atoms[b.B2].Aromatic) Dearomatise(b.A, b.B2);
        b.Order = order;
        _atoms[b.A].Order.Clear();
        _atoms[b.B2].Order.Clear();
    }

    /// <summary>An aromatic ring touched by a bond edit is written out as a Kekulé structure.</summary>
    private void Dearomatise(params int[] seeds)
    {
        var ring = new HashSet<int>();
        var stack = new Stack<int>(seeds);
        while (stack.Count > 0)
        {
            var u = stack.Pop();
            if (!_atoms[u].Aromatic || !ring.Add(u)) continue;
            foreach (var v in Neighbours(u)) stack.Push(v);
        }
        foreach (var i in ring) { _atoms[i].Aromatic = false; _atoms[i].Order.Clear(); _atoms[i].Chiral = 0; }
        foreach (var b in _bonds.Where(b => b.Order == 4 && ring.Contains(b.A) && ring.Contains(b.B2))) b.Order = 1;
        // alternate doubles greedily
        foreach (var b in _bonds.Where(b => ring.Contains(b.A) && ring.Contains(b.B2) && b.Order == 1))
            if (!HasDouble(b.A) && !HasDouble(b.B2) && ImplicitH(b.A) > 0 && ImplicitH(b.B2) > 0) b.Order = 2;
    }

    private static int[] Valences(int z) => z switch
    {
        5 => [3], 6 => [4], 7 => [3, 5], 8 => [2], 15 => [3, 5], 16 => [2, 4, 6], 9 or 17 or 35 or 53 => [1],
        // outside the SMILES organic subset: their usual valences, so a drawn atom carries its hydrogens (SiH4, GeH4, SnH4,
        // AsH3, SbH3, SeH2, TeH2, BiH3)
        14 or 32 or 50 => [4], 33 or 51 or 83 => [3, 5], 34 or 52 => [2, 4, 6], 1 => [1], _ => [],
    };

    /// <summary>Hydrogens on atom i as the SMILES reader will give them (for labels).</summary>
    private int ImplicitH(int i)
    {
        var a = _atoms[i];
        if (a.HCount >= 0 && a.Bracket) return a.HCount;
        return ValenceH(i, a.Z);
    }

    /// <summary>Hydrogens atom i takes from the valence of element z (its own, or the one a charged atom is isoelectronic
    /// with) and the bonds it has.</summary>
    private int ValenceH(int i, int z)
    {
        var a = _atoms[i];
        var sum = 0;
        var arom = 0;
        foreach (var b in _bonds) if (b.A == i || b.B2 == i) { sum += b.Order == 4 ? 1 : b.Order; arom += b.Order == 4 ? 1 : 0; }
        if (a.Aromatic && arom > 0)
        {
            if (z is 8 or 16) return 0;
            return Math.Max(0, (z == 6 ? 4 : 3) - (sum + 1));
        }
        foreach (var v in Valences(z)) if (v >= sum) return v - sum;
        return 0;
    }

    /// <summary>The bond orders on atom i, an aromatic bond counted as 1.</summary>
    private int BondSum(int i)
    {
        var sum = 0;
        foreach (var b in _bonds) if (b.A == i || b.B2 == i) sum += b.Order == 4 ? 1 : b.Order;
        return sum;
    }

    private bool HasAromaticBond(int i) => _bonds.Any(b => (b.A == i || b.B2 == i) && b.Order == 4);

    /// <summary>The double bond an aromatic atom has in its Kekulé form (1 or 0), from its valence, bonds and hydrogens:
    /// 1 for benzene's c, pyridine's n and [nH+]; 0 for pyrrole's [nH], [se] and [cH-], whose ring bonds are all single.</summary>
    private int AromaticPi(int i)
    {
        var a = _atoms[i];
        if (!a.Aromatic || !HasAromaticBond(i)) return 0;
        var used = BondSum(i) + ImplicitH(i);
        foreach (var v in Valences(a.Z - a.Charge)) if (v >= used) return Math.Min(1, v - used);
        return 0;
    }

    // ---------------------------------------------------------------- drawing
    private static readonly Dictionary<int, string> Symbols = new()
    {
        [1] = "H", [5] = "B", [6] = "C", [7] = "N", [8] = "O", [9] = "F", [14] = "Si", [15] = "P", [16] = "S", [17] = "Cl", [35] = "Br", [53] = "I",
        [11] = "Na", [19] = "K", [3] = "Li", [12] = "Mg", [20] = "Ca", [26] = "Fe", [29] = "Cu", [30] = "Zn",
    };

    private static IBrush ElementBrush(int z) => z switch
    {
        7 => new SolidColorBrush(Color.Parse("#5B93E6")),
        8 => new SolidColorBrush(Color.Parse("#EE6A63")),
        9 or 17 => new SolidColorBrush(Color.Parse("#7DC884")),
        16 => new SolidColorBrush(Color.Parse("#D6AC5C")),
        35 => new SolidColorBrush(Color.Parse("#DE775D")),
        15 => new SolidColorBrush(Color.Parse("#F0A83C")),
        _ => Tokens.Brush("TextB"),
    };

    private bool Labelled(int i)
    {
        var a = _atoms[i];
        return a.Z != 6 || a.Charge != 0 || a.Isotope != 0 || Neighbours(i).Count == 0;
    }

    public override void Render(DrawingContext ctx)
    {
        var bg = Tokens.Brush("Bg0B");
        ctx.FillRectangle(bg, new Rect(Bounds.Size));
        // dot grid
        var dot = Tokens.Brush("Bg3B");
        for (var x = 12.0; x < Bounds.Width; x += 24)
            for (var y = 12.0; y < Bounds.Height; y += 24)
                ctx.FillRectangle(dot, new Rect(x, y, 1.5, 1.5));
        var ink = Tokens.Brush("TextB");
        var pen = new Pen(ink, 1.6, lineCap: PenLineCap.Round);
        var thin = new Pen(ink, 1.4, lineCap: PenLineCap.Round);
        var dash = new Pen(ink, 1.2, new DashStyle([2.5, 2.5], 0), PenLineCap.Flat);
        var rings = SmallRings();

        for (var k = 0; k < _bonds.Count; k++)
        {
            var b = _bonds[k];
            var pa = ToScreen(_atoms[b.A].X, _atoms[b.A].Y);
            var pb = ToScreen(_atoms[b.B2].X, _atoms[b.B2].Y);
            // shorten at labelled atoms
            var d = pb - pa;
            var len = Math.Sqrt(d.X * d.X + d.Y * d.Y);
            if (len < 1) continue;
            var u = new Point(d.X / len, d.Y / len);
            var sa = Labelled(b.A) ? 0.3 * _scale : 0;
            var sb = Labelled(b.B2) ? 0.3 * _scale : 0;
            var a1 = pa + u * sa;
            var b1 = pb - u * sb;
            var nrm = new Point(-u.Y, u.X);
            var hover = k == _hoverBond;
            var p = hover ? new Pen(Tokens.Brush("SelB"), 2.2, lineCap: PenLineCap.Round) : pen;
            // the inner side: toward the centre of a ring the bond is in
            var centre = rings.FirstOrDefault(r => r.Contains(b.A) && r.Contains(b.B2));
            if (b.Order == 1) ctx.DrawLine(p, a1, b1);
            else if (b.Order == 3)
            {
                ctx.DrawLine(p, a1, b1);
                var o = nrm * (0.14 * _scale);
                ctx.DrawLine(thin, a1 + o, b1 + o);
                ctx.DrawLine(thin, a1 - o, b1 - o);
            }
            else if (centre != null || (Neighbours(b.A).Count > 1 && Neighbours(b.B2).Count > 1 && b.Order == 4))
            {
                ctx.DrawLine(p, a1, b1);
                var side = 1.0;
                if (centre != null)
                {
                    var cx = centre.Average(i => _atoms[i].X);
                    var cy = centre.Average(i => _atoms[i].Y);
                    var c = ToScreen(cx, cy);
                    side = Math.Sign((c.X - pa.X) * nrm.X + (c.Y - pa.Y) * nrm.Y);
                    if (side == 0) side = 1;
                }
                var o = nrm * (0.17 * _scale * side);
                var shrink = u * (0.14 * _scale);
                ctx.DrawLine(b.Order == 4 ? dash : thin, a1 + o + shrink, b1 + o - shrink);
            }
            else
            {
                // a double bond outside rings: centred pair (a terminal =O or =CH2 hangs straight)
                var o = nrm * (0.085 * _scale);
                ctx.DrawLine(p, a1 + o, b1 + o);
                ctx.DrawLine(p, a1 - o, b1 - o);
            }
        }

        // drag preview
        if (_dragFrom >= 0 && _dragging)
        {
            var pa = ToScreen(_atoms[_dragFrom].X, _atoms[_dragFrom].Y);
            ctx.DrawLine(new Pen(Tokens.Brush("SelB"), 1.6, new DashStyle([4, 3], 0)), pa, _dragAt);
        }

        // atom labels
        var font = new Typeface(Tokens.Mono);
        for (var i = 0; i < _atoms.Count; i++)
        {
            var a = _atoms[i];
            var c = ToScreen(a.X, a.Y);
            if (Labelled(i))
            {
                var h = ImplicitH(i);
                var sym = Sym(a.Z);
                var text = sym + (h > 0 ? "H" + (h > 1 ? Sub(h) : "") : "");
                if (a.Charge != 0) text += a.Charge > 0 ? (a.Charge > 1 ? $"{a.Charge}+" : "+") : (a.Charge < -1 ? $"{-a.Charge}−" : "−");
                var ft = new FormattedText(text, CultureInfo.InvariantCulture, FlowDirection.LeftToRight, font, 0.36 * _scale + 3, ElementBrush(a.Z));
                ctx.FillRectangle(bg, new Rect(c.X - ft.Width / 2 - 2, c.Y - ft.Height / 2, ft.Width + 4, ft.Height));
                ctx.DrawText(ft, new Point(c.X - ft.Width / 2, c.Y - ft.Height / 2));
            }
            if (a.Chiral != 0)
            {
                var ft = new FormattedText(a.Chiral == 1 ? "@" : "@@", CultureInfo.InvariantCulture, FlowDirection.LeftToRight, font, 9, Tokens.Brush("DimB"));
                ctx.DrawText(ft, new Point(c.X + 5, c.Y + 3));
            }
        }

        // Bond tool: the first atom of a pair, ringed until the second is clicked
        if (Mode == Tool.Bond && _pending >= 0 && _pending < _atoms.Count)
        {
            var c = ToScreen(_atoms[_pending].X, _atoms[_pending].Y);
            ctx.DrawEllipse(null, new Pen(Tokens.Brush("AccB"), 2), c, 0.34 * _scale, 0.34 * _scale);
        }

        // hover: the atom and what a click will do
        if (_hoverAtom >= 0 && _hoverAtom < _atoms.Count && _dragFrom < 0)
        {
            var a = _atoms[_hoverAtom];
            var c = ToScreen(a.X, a.Y);
            var sel = Tokens.Brush("SelB");
            ctx.DrawEllipse(null, new Pen(sel, 1.5), c, 0.3 * _scale, 0.3 * _scale);
            var h = ImplicitH(_hoverAtom);
            var sym = Sym(a.Z);
            var what = Mode switch
            {
                Tool.Draw => a.Z != Element ? $"click: make it {Sym(Element)}" : "drag to grow",
                Tool.Ring => "click: spiro ring",
                Tool.Erase => "click: remove",
                Tool.Bond => _pending >= 0 && _pending != _hoverAtom ? "click: bond to the first atom" : _pending == _hoverAtom ? "click again: a new atom on it" : "click, then another atom: bond them",
                Tool.Cut => "click a bond to cut it",
                Tool.Move => "drag to move",
                _ => "click: charge",
            };
            var hint = new FormattedText($"{sym}{(h > 0 ? "H" + (h > 1 ? Sub(h) : "") : "")} · {what}", CultureInfo.InvariantCulture, FlowDirection.LeftToRight, font, 11, sel);
            ctx.DrawText(hint, new Point(c.X + 0.4 * _scale, c.Y - 0.55 * _scale - hint.Height / 2));
        }

        if (_atoms.Count == 0)
        {
            var ft = new FormattedText("Click to place an atom, or type a SMILES on the right", CultureInfo.InvariantCulture, FlowDirection.LeftToRight,
                new Typeface(Tokens.Sans), 13, Tokens.Brush("DimB"));
            ctx.DrawText(ft, new Point((Bounds.Width - ft.Width) / 2, (Bounds.Height - ft.Height) / 2));
        }
    }

    private static string Sub(int n) => string.Concat(n.ToString(CultureInfo.InvariantCulture).Select(c => (char)('₀' + (c - '0'))));

    /// <summary>Rings up to eight atoms (the smallest through each bond), for the inner lines of ring double bonds.</summary>
    private List<List<int>> SmallRings()
    {
        var rings = new List<List<int>>();
        foreach (var b in _bonds)
        {
            // shortest path b.A → b.B2 without the bond itself
            var prev = new Dictionary<int, int> { [b.A] = -1 };
            var q = new Queue<int>();
            q.Enqueue(b.A);
            var found = false;
            while (q.Count > 0 && !found)
            {
                var u = q.Dequeue();
                foreach (var v in Neighbours(u))
                {
                    if (u == b.A && v == b.B2) continue;
                    if (prev.ContainsKey(v)) continue;
                    prev[v] = u;
                    if (v == b.B2) { found = true; break; }
                    q.Enqueue(v);
                }
            }
            if (!found) continue;
            var path = new List<int>();
            for (var v = b.B2; v != -1; v = prev[v]) path.Add(v);
            if (path.Count > 8) continue;
            path.Sort();
            if (!rings.Any(r => r.SequenceEqual(path))) rings.Add(path);
        }
        return rings;
    }
}
