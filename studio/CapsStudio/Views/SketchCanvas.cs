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
    public enum Tool { Draw, Ring, Erase, Charge }

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
    public Tool Mode { get; set; } = Tool.Draw;
    public int Element { get; set; } = 6;
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
        Fit();
        InvalidateVisual();
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
        foreach (var a in _atoms)
            atoms.Add(new JsonObject
            {
                ["z"] = a.Z, ["charge"] = a.Charge, ["isotope"] = a.Isotope, ["hcount"] = a.HCount, ["aromatic"] = a.Aromatic,
                ["bracket"] = a.Bracket, ["chiral"] = a.Chiral, ["map"] = a.Map, ["order"] = new JsonArray(a.Order.Select(v => (JsonNode)v).ToArray()),
            });
        var bonds = new JsonArray();
        foreach (var b in _bonds) bonds.Add(new JsonObject { ["a"] = b.A, ["b"] = b.B2, ["order"] = b.Order, ["dir"] = b.Dir });
        return new JsonObject { ["atoms"] = atoms, ["bonds"] = bonds }.ToJsonString();
    }

    private void Changed(params int[] touched)
    {
        // an edit next to a stereocentre or a marked double bond drops that stereo (it would no longer mean the same)
        foreach (var i in touched)
            if (i >= 0 && i < _atoms.Count)
            {
                _atoms[i].Chiral = 0;
                foreach (var b in _bonds) if (b.A == i || b.B2 == i) b.Dir = 0;
            }
        string smiles;
        try { smiles = _atoms.Count == 0 ? "" : CapsDocument.SmilesWrite(Graph()); }
        catch { smiles = ""; }
        InvalidateVisual();
        Edited?.Invoke(smiles);
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

    /// <summary>A benzene ring (Kekulé) through atom `at` or on bond `bond`, or free at (x, y).</summary>
    private void AddRing(double x, double y, int at, int bond)
    {
        var ids = new List<int>();
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
            var cx = mx + nx / len * Math.Sqrt(3) / 2;
            var cy = my + ny / len * Math.Sqrt(3) / 2;
            var t0 = Math.Atan2(pa.Y - cy, pa.X - cx);
            var t1 = Math.Atan2(pb.Y - cy, pb.X - cx);
            var step = Gap(t0 + Math.PI / 3, t1) < 0.1 ? Math.PI / 3 : -Math.PI / 3;
            ids.Add(b.A);
            ids.Add(b.B2);
            for (var k = 2; k < 6; k++) ids.Add(AddAtomAt(cx + Math.Cos(t0 + k * step), cy + Math.Sin(t0 + k * step)));
        }
        else if (at >= 0)
        {
            var (dx, dy) = FreeDirection(at);
            var cx = _atoms[at].X + dx;
            var cy = _atoms[at].Y + dy;
            var t0 = Math.Atan2(_atoms[at].Y - cy, _atoms[at].X - cx);
            ids.Add(at);
            for (var k = 1; k < 6; k++) ids.Add(AddAtomAt(cx + Math.Cos(t0 + k * Math.PI / 3), cy + Math.Sin(t0 + k * Math.PI / 3)));
        }
        else
        {
            for (var k = 0; k < 6; k++) ids.Add(AddAtomAt(x + Math.Cos(Math.PI / 6 + k * Math.PI / 3), y + Math.Sin(Math.PI / 6 + k * Math.PI / 3)));
        }
        // alternate double bonds, keeping a double bond the atoms already have
        for (var k = 0; k < 6; k++)
        {
            var a = ids[k];
            var b = ids[(k + 1) % 6];
            if (!_bonds.Any(x => (x.A == a && x.B2 == b) || (x.A == b && x.B2 == a))) AddBond(a, b, 1);
        }
        var startOdd = HasDouble(ids[0]) || HasDouble(ids[1]) ? 1 : 0;
        for (var k = startOdd; k < 6; k += 2)
        {
            var a = ids[k];
            var b = ids[(k + 1) % 6];
            var bb = _bonds.First(x => (x.A == a && x.B2 == b) || (x.A == b && x.B2 == a));
            if (!HasDouble(a) && !HasDouble(b)) bb.Order = 2;
        }
        Changed(ids.ToArray());
    }

    private bool HasDouble(int i) => _bonds.Any(b => (b.A == i || b.B2 == i) && (b.Order == 2 || b.Order == 4));

    /// <summary>An atom at (x, y), or the existing one there (rings fuse onto what is drawn).</summary>
    private int AddAtomAt(double x, double y)
    {
        for (var i = 0; i < _atoms.Count; i++)
            if (Math.Abs(_atoms[i].X - x) < 0.2 && Math.Abs(_atoms[i].Y - y) < 0.2) return i;
        _atoms.Add(new A { Z = 6, X = x, Y = y });
        return _atoms.Count - 1;
    }

    // ---------------------------------------------------------------- input
    protected override void OnPointerPressed(PointerPressedEventArgs e)
    {
        base.OnPointerPressed(e);
        Focus();
        var p = e.GetPosition(this);
        _pressAt = p;
        _dragAt = p;
        _dragFrom = Mode == Tool.Draw ? AtomAt(p) : -1;
        _dragging = false;
        e.Pointer.Capture(this);
    }

    protected override void OnPointerMoved(PointerEventArgs e)
    {
        base.OnPointerMoved(e);
        var p = e.GetPosition(this);
        if (_dragFrom >= 0 && e.GetCurrentPoint(this).Properties.IsLeftButtonPressed)
        {
            _dragAt = p;
            if (Dist(p, _pressAt) > 6) _dragging = true;
            InvalidateVisual();
            return;
        }
        var ha = AtomAt(p);
        var hb = ha < 0 ? BondAt(p) : -1;
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
        switch (Mode)
        {
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
                    b.Order = b.Order switch { 1 => 2, 2 => 3, 4 => 2, _ => 1 };
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
                    AddRing(mx, my, at, bond);
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
                    a.HCount = a.Bracket ? ImplicitH(at) : -1;
                    Changed(at);
                }
                break;
        }
        _dragging = false;
        _hoverAtom = AtomAt(p);
        InvalidateVisual();
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
        5 => [3], 6 => [4], 7 => [3, 5], 8 => [2], 15 => [3, 5], 16 => [2, 4, 6], 9 or 17 or 35 or 53 => [1], _ => [],
    };

    /// <summary>Hydrogens on atom i as the SMILES reader will give them (for labels).</summary>
    private int ImplicitH(int i)
    {
        var a = _atoms[i];
        if (a.HCount >= 0 && a.Bracket) return a.HCount;
        var sum = 0;
        var arom = 0;
        foreach (var b in _bonds) if (b.A == i || b.B2 == i) { sum += b.Order == 4 ? 1 : b.Order; arom += b.Order == 4 ? 1 : 0; }
        if (a.Aromatic && arom > 0)
        {
            if (a.Z is 8 or 16) return 0;
            return Math.Max(0, (a.Z == 6 ? 4 : 3) - (sum + 1));
        }
        foreach (var v in Valences(a.Z)) if (v >= sum) return v - sum;
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
                var sym = Symbols.TryGetValue(a.Z, out var s) ? s : "*";
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

        // hover: the atom and what a click will do
        if (_hoverAtom >= 0 && _hoverAtom < _atoms.Count && _dragFrom < 0)
        {
            var a = _atoms[_hoverAtom];
            var c = ToScreen(a.X, a.Y);
            var sel = Tokens.Brush("SelB");
            ctx.DrawEllipse(null, new Pen(sel, 1.5), c, 0.3 * _scale, 0.3 * _scale);
            var h = ImplicitH(_hoverAtom);
            var sym = Symbols.TryGetValue(a.Z, out var s) ? s : "*";
            var what = Mode switch
            {
                Tool.Draw => a.Z != Element ? $"click: make it {(Symbols.TryGetValue(Element, out var e) ? e : "?")}" : "drag to grow",
                Tool.Ring => "click: spiro ring",
                Tool.Erase => "click: remove",
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
