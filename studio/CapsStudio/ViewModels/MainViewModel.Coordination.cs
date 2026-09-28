using System.Globalization;

namespace CapsStudio.ViewModels;

/// <summary>Coordination geometry of an atom (design spec §11: linear through octahedral), coordinate (dative) bonds for
/// the bond tool, and the out-of-plane angle in the live monitor.</summary>
public sealed partial class MainViewModel
{
    public static readonly string[] CoordinationNames = ["linear", "trigonal planar", "tetrahedral", "square planar", "trigonal bipyramidal", "square pyramidal", "octahedral"];
    private static readonly string[] CoordinationIds = ["linear", "trigonal", "tetrahedral", "square_planar", "trigonal_bipyramidal", "square_pyramidal", "octahedral"];
    private int _coordination = 2;
    public int CoordinationIndex { get => _coordination; set => Set(ref _coordination, Math.Clamp(value, 0, CoordinationNames.Length - 1)); }

    /// <summary>The picked atom's neighbours placed at the ideal directions of the chosen geometry (undo with ⌘Z).</summary>
    public void ApplyCoordination()
    {
        if (_selection.Count != 1 || _doc == null) { Status = "Pick the centre atom (one atom)"; return; }
        var c = _selection[0];
        if (RunEdit(new { op = "set_coordination", atom = c, geometry = CoordinationIds[_coordination] }) != null) Pick(c);
    }

    /// <summary>Bond tool: the kind of bond drawn — single, double, triple, aromatic or coordinate (dative, a ligand's lone
    /// pair to a metal: in no atom's valence).</summary>
    public static readonly string[] BondKinds = ["Single", "Double", "Triple", "Aromatic", "Coordinate (dative)"];
    private int _bondKind;
    public int BondKind { get => _bondKind; set => Set(ref _bondKind, Math.Clamp(value, 0, BondKinds.Length - 1)); }
    private int BondOrderCode => _bondKind switch { 1 => 2, 2 => 3, 3 => 4, 4 => 9, _ => 1 };

    // ideal directions (as the core's coordination_directions)
    private static readonly Dictionary<string, double[][]> IdealSites = BuildSites();
    private static Dictionary<string, double[][]> BuildSites()
    {
        double r3 = 1 / Math.Sqrt(3), c = -0.5, s = Math.Sqrt(3) / 2;
        return new()
        {
            ["linear"] = [[0, 0, 1], [0, 0, -1]],
            ["trigonal planar"] = [[1, 0, 0], [c, s, 0], [c, -s, 0]],
            ["tetrahedral"] = [[r3, r3, r3], [r3, -r3, -r3], [-r3, r3, -r3], [-r3, -r3, r3]],
            ["square planar"] = [[1, 0, 0], [0, 1, 0], [-1, 0, 0], [0, -1, 0]],
            ["trigonal bipyramidal"] = [[0, 0, 1], [0, 0, -1], [1, 0, 0], [c, s, 0], [c, -s, 0]],
            ["square pyramidal"] = [[0, 0, 1], [1, 0, 0], [0, 1, 0], [-1, 0, 0], [0, -1, 0]],
            ["octahedral"] = [[1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0], [0, 0, 1], [0, 0, -1]],
        };
    }
    private static double AngleDeg(double[] a, double[] b)
    {
        double d = 0, na = 0, nb = 0;
        for (int k = 0; k < 3; ++k) { d += a[k] * b[k]; na += a[k] * a[k]; nb += b[k] * b[k]; }
        return Math.Acos(Math.Clamp(d / Math.Sqrt(na * nb), -1, 1)) * 180 / Math.PI;
    }

    /// <summary>The picked atom's coordination: its bonded neighbours and the ideal shape (with as many or more sites) whose
    /// sorted ligand–centre–ligand angles are closest (RMS, degrees); "" for fewer than two neighbours.</summary>
    private string CoordinationOf(int centre)
    {
        if (_doc == null) return "";
        var nb = _doc.Bonded(centre);
        if (nb.Length < 2) return nb.Length == 1 ? "1 neighbour" : "";
        var a = _doc.Atom(centre);
        var u = nb.Select(i => { var b = _doc.Atom(i); return new[] { b.X - a.X, b.Y - a.Y, b.Z - a.Z }; }).ToArray();
        var have = new List<double>();
        for (int i = 0; i < u.Length; ++i) for (int j = i + 1; j < u.Length; ++j) have.Add(AngleDeg(u[i], u[j]));
        have.Sort();
        var (best, rms) = ("", double.MaxValue);
        foreach (var (name, sites) in IdealSites)
        {
            if (sites.Length < u.Length) continue;
            foreach (var subset in Subsets(sites.Length, u.Length))
            {
                var ideal = new List<double>();
                for (int i = 0; i < subset.Length; ++i) for (int j = i + 1; j < subset.Length; ++j) ideal.Add(AngleDeg(sites[subset[i]], sites[subset[j]]));
                ideal.Sort();
                var r = Math.Sqrt(have.Zip(ideal, (x, y) => (x - y) * (x - y)).Average());
                // a shape with its own count wins a tie over a larger one with sites left empty
                if (r < rms - 1e-6 || (Math.Abs(r - rms) <= 1e-6 && sites.Length == u.Length)) (best, rms) = (sites.Length == u.Length ? name : $"{name}, {sites.Length - u.Length} site{(sites.Length - u.Length > 1 ? "s" : "")} empty", r);
            }
        }
        return string.Format(CultureInfo.InvariantCulture, "{0} neighbours · {1} (RMS {2:0.0}°)", u.Length, best, rms);
    }
    private static IEnumerable<int[]> Subsets(int n, int k)
    {
        var idx = Enumerable.Range(0, k).ToArray();
        while (true)
        {
            yield return (int[])idx.Clone();
            int i = k - 1;
            while (i >= 0 && idx[i] == n - k + i) --i;
            if (i < 0) yield break;
            ++idx[i];
            for (int j = i + 1; j < k; ++j) idx[j] = idx[j - 1] + 1;
        }
    }

    /// <summary>Four picks with the first bonded to the other three: the out-of-plane (Wilson) angle of the last bond from
    /// the plane of the centre and the two others, and the mean of the three such angles (as class II impropers use).</summary>
    private string OutOfPlane(int[] p)
    {
        if (_doc == null || p.Length != 4) return "";
        var nb = _doc.Bonded(p[0]);
        if (!nb.Contains(p[1]) || !nb.Contains(p[2]) || !nb.Contains(p[3])) return "";
        var c = _doc.Atom(p[0]);
        double[] V(int i) { var b = _doc.Atom(i); return [b.X - c.X, b.Y - c.Y, b.Z - c.Z]; }
        static double[] Cross(double[] x, double[] y) => [x[1] * y[2] - x[2] * y[1], x[2] * y[0] - x[0] * y[2], x[0] * y[1] - x[1] * y[0]];
        static double Wilson(double[] a, double[] b, double[] d)   // the angle of d from the plane of a, b
        {
            var n = Cross(a, b);
            double nn = Math.Sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]), dd = Math.Sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);
            return nn < 1e-12 || dd < 1e-12 ? 0 : Math.Asin(Math.Clamp((n[0] * d[0] + n[1] * d[1] + n[2] * d[2]) / (nn * dd), -1, 1)) * 180 / Math.PI;
        }
        double[] a = V(p[1]), b = V(p[2]), d = V(p[3]);
        var chi = Wilson(a, b, d);
        var mean = (Math.Abs(Wilson(a, b, d)) + Math.Abs(Wilson(b, d, a)) + Math.Abs(Wilson(d, a, b))) / 3;
        return string.Format(CultureInfo.InvariantCulture, "Out of plane {0} from {1}–{2}–{3}: {4:F2}° (mean of three {5:F2}°)",
                             _doc.Atom(p[3]).Id, _doc.Atom(p[1]).Id, _doc.Atom(p[0]).Id, _doc.Atom(p[2]).Id, chi, mean);
    }
}
