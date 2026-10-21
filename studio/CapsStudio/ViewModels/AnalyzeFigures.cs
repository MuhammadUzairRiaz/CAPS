using System.Globalization;
using CapsStudio.Views;

namespace CapsStudio.ViewModels;

/// <summary>Turns the curves of caps analyze into figures: each property's own plot (g(r) with n(r), the MSD on log axes
/// with the Einstein fit window, V(T) with the two lines and Tg, a time series with its block means and mean ± standard
/// error), with axis titles written as "quantity / unit". Everything drawn comes from the report: the curves, the values
/// and errors, and fits recomputed from the same points (the Einstein line, the KWW decays from their reported τ and β).</summary>
public static class AnalyzeFigures
{
    private static readonly CultureInfo Inv = CultureInfo.InvariantCulture;

    // ------------------------------------------------------------------ axis titles
    private static readonly Dictionary<string, string> Symbols = new()
    {
        ["time"] = "*t*", ["t"] = "*t*", ["τ"] = "*τ*", ["r"] = "*r*", ["R"] = "*R*", ["Rg"] = "*R*_{g}", ["T"] = "*T*", ["q"] = "*q*",
        ["k"] = "*k*", ["n"] = "*n*", ["density"] = "*ρ*", ["g(r)"] = "*g*(*r*)", ["specific volume"] = "*v*", ["coordination number"] = "*n*(*r*)",
        ["strain"] = "*ε*", ["axial strain"] = "*ε*_{axial}", ["stress"] = "*σ*", ["η"] = "*η*", ["γ̇"] = "*γ̇*", ["ε"] = "*ε*", ["J"] = "*J*",
        ["φ"] = "*φ*", ["θ"] = "*θ*", ["P₂"] = "*P*₂", ["S"] = "*S*", ["wavenumber"] = "*ν̃*", ["ν̃"] = "*ν̃*", ["pore diameter"] = "*d*",
        ["probe radius"] = "*r*_{probe}", ["displacement"] = "*d*", ["force"] = "*F*", ["count"] = "*N*_{HB}", ["C(t)"] = "*C*(*t*)",
        ["C(τ)/C(0)"] = "*C*(*τ*) / *C*(0)", ["C_M(t)"] = "*C*_{M}(*t*)", ["α₂"] = "*α*₂", ["CED"] = "CED", ["MSD"] = "MSD",
        ["d ln MSD / d ln t"] = "d ln MSD / d ln *t*", ["⟨u_i · u_{i+k}⟩"] = "⟨**u**_{*i*} · **u**_{*i*+*k*}⟩",
        ["C_n = ⟨R²(n)⟩ / (n ⟨b²⟩)"] = "*C*_{*n*}", ["autocorrelation"] = "autocorrelation", ["probability density"] = "probability density",
        ["temperature"] = "*T*", ["pressure"] = "*P*", ["volume"] = "*V*", ["energy"] = "*E*", ["E"] = "*E*", ["potential energy"] = "*E*_{pot}", ["potential"] = "*E*_{pot}",
        ["kinetic energy"] = "*E*_{kin}", ["kinetic"] = "*E*_{kin}", ["total energy"] = "*E*_{tot}", ["total"] = "*E*_{tot}", ["conserved"] = "*E*_{cons}", ["step"] = "step",
        ["E_pot"] = "*E*_{pot}", ["P"] = "*P*", ["V"] = "*V*", ["radius of gyration"] = "*R*_{g}", ["end-to-end distance"] = "*R*_{ee}", ["largest force"] = "|*F*|_{max}",
    };

    /// <summary>"name (unit)" from the core as "symbol / unit" with the unit in powers (g/cm³ → g cm⁻³).</summary>
    public static string Title(string core)
    {
        core = core.Trim();
        string name = core, unit = "";
        if (core.EndsWith(')'))
        {
            var depth = 0;
            for (var i = core.Length - 1; i >= 0; --i)
            {
                if (core[i] == ')') depth++;
                else if (core[i] == '(' && --depth == 0)
                {
                    if (i > 0 && core[i - 1] == ' ') { name = core[..(i - 1)].Trim(); unit = core[(i + 1)..^1]; }
                    break;
                }
            }
        }
        if (unit.StartsWith("per ", StringComparison.Ordinal)) unit = "1/" + unit[4..];
        var sym = Symbols.TryGetValue(name, out var s) ? s
            : name.Length > 2 && Symbols.TryGetValue(name.ToLowerInvariant(), out var s2) ? s2   // "Time (ps)", "Temperature (K)" from a table's header
            : name.Length <= 2 && name.All(char.IsLetter) ? "*" + name + "*" : name;
        if (unit is "" or "arb.") return unit == "arb." ? sym + " / arb. units" : sym;
        if (unit == "relative") return sym + " (relative)";
        return sym + " / " + Unit(unit);
    }

    private const string Sup = "⁰¹²³⁴⁵⁶⁷⁸⁹";
    /// <summary>a/b/c → a b⁻¹ c⁻¹; 1/x → x⁻¹; a·b → a b; ° → deg in a denominator.</summary>
    public static string Unit(string u)
    {
        u = u.Replace("^½", "^{1/2}");
        var parts = u.Split('/');
        var num = parts[0].Trim().Replace("·", " ");
        var outp = new List<string>();
        if (num != "1" && num.Length > 0) outp.Add(num);
        foreach (var raw in parts.Skip(1))
            foreach (var f in raw.Trim().Trim('(', ')').Split('·', ' ', StringSplitOptions.RemoveEmptyEntries))
            {
                var g = f == "°" ? "deg" : f;
                var k = g.Length;
                while (k > 0 && Sup.Contains(g[k - 1])) k--;
                var basePart = g[..k];
                var exp = g[k..];
                outp.Add(basePart + "⁻" + (exp.Length > 0 ? exp : "¹"));
            }
        return string.Join(" ", outp);
    }

    private static string Num(double v) => ResultCard.Num(v);
    private static string ValErr(ResultCard c)
    {
        if (double.IsNaN(c.Error) || c.Error <= 0) return Num(c.Value);
        var dec = Math.Clamp(1 - (int)Math.Floor(Math.Log10(c.Error)), 0, 8);
        return c.Value.ToString("F" + dec, Inv) + " ± " + c.Error.ToString("F" + dec, Inv);
    }
    private static double Extra(ResultCard c, string prefix) =>
        c.Extra.Where(e => e.Key.StartsWith(prefix, StringComparison.Ordinal)).Select(e => e.Value).DefaultIfEmpty(double.NaN).First();

    private static SciSeries S(SeriesItem it, string name, SciKind kind = SciKind.Line, int colour = 0) => new(name, it.X, it.Y) { Kind = kind, Colour = colour };

    private static SciFigure Generic(SeriesItem it) => new()
    {
        XLabel = Title(it.XLabel), YLabel = Title(it.YLabel), XLog = it.LogLog, YLog = it.LogLog, RefY = it.RefY,
        Series = [S(it, "", it.X.Length <= 40 ? SciKind.LineMarkers : SciKind.Line)],
    };

    /// <summary>Block means of a time series the way the core takes them (equal consecutive blocks, the tail left out).</summary>
    private static (double[] T, double[] M) Blocks(double[] x, double[] y, int nb)
    {
        var n = Math.Min(x.Length, y.Length);
        nb = Math.Min(nb, n);
        if (nb < 2) return ([], []);
        var per = n / nb;
        var t = new double[nb];
        var m = new double[nb];
        for (var b = 0; b < nb; ++b)
        {
            double sx = 0, sy = 0;
            for (var k = b * per; k < (b + 1) * per; ++k) { sx += x[k]; sy += y[k]; }
            t[b] = sx / per;
            m[b] = sy / per;
        }
        return (t, m);
    }

    private static SciFigure TimeSeries(ResultCard c, SeriesItem it, string unitText, int blocks)
    {
        var (bt, bm) = Blocks(it.X, it.Y, blocks);
        var series = new List<SciSeries> { S(it, "per frame", SciKind.Line, -1) with { Width = 0.9 } };
        if (bt.Length >= 2) series.Add(new SciSeries($"block means ({bt.Length} blocks)", bt, bm) { Kind = SciKind.Markers, Colour = 0, Square = true });
        var hasErr = double.IsFinite(c.Error) && c.Error > 0;
        return new SciFigure
        {
            XLabel = Title(it.XLabel), YLabel = Title(it.YLabel), Series = [.. series],
            YBand = double.IsFinite(c.Value) ? (c.Value - (hasErr ? c.Error : 0), c.Value + (hasErr ? c.Error : 0), "mean" + (hasErr ? " ± s.e." : "") + ": " + ValErr(c) + " " + unitText) : null,
        };
    }

    private static (double A, double B) LinFit(IReadOnlyList<double> x, IReadOnlyList<double> y)
    {
        var n = x.Count;
        if (n < 2) return (double.NaN, double.NaN);
        double mx = x.Average(), my = y.Average(), sxx = 0, sxy = 0;
        for (var i = 0; i < n; ++i) { sxx += (x[i] - mx) * (x[i] - mx); sxy += (x[i] - mx) * (y[i] - my); }
        var b = sxy / sxx;
        return (my - b * mx, b);
    }

    // Lanczos Γ (g = 7, n = 9), for the KWW integral τ = τ_KWW Γ(1/β) / β
    private static double Gamma(double z)
    {
        double[] p = [0.99999999999980993, 676.5203681218851, -1259.1392167224028, 771.32342877765313, -176.61502916214059,
            12.507343278686905, -0.13857109526572012, 9.9843695780195716e-6, 1.5056327351493116e-7];
        if (z < 0.5) return Math.PI / (Math.Sin(Math.PI * z) * Gamma(1 - z));
        z -= 1;
        var a = p[0];
        var t = z + 7.5;
        for (var i = 1; i < 9; ++i) a += p[i] / (z + i);
        return Math.Sqrt(2 * Math.PI) * Math.Pow(t, z + 0.5) * Math.Exp(-t) * a;
    }

    private static SciSeries Kww(double[] t, double tau, double beta, string name, int colour)
    {
        var tmax = t.Where(double.IsFinite).DefaultIfEmpty(1).Max();
        var xs = Enumerable.Range(0, 121).Select(k => tmax * k / 120.0).ToArray();
        return new SciSeries(name, xs, xs.Select(x => Math.Exp(-Math.Pow(x / tau, beta))).ToArray()) { Kind = SciKind.Dashed, Colour = colour, Width = 1.3 };
    }

    /// <summary>The figures of one property: the entries for the curve list (the first is also the card's) and the card
    /// figure. all: every card and its raw curves, for properties drawn over another's curve (D over the MSD).</summary>
    public static List<SeriesItem> Build(ResultCard card, List<SeriesItem> raw, IReadOnlyDictionary<string, (ResultCard Card, List<SeriesItem> Raw)> all, int blocks)
    {
        var id = card.Id;
        var unit = card.UnitText.Length > 0 ? Unit(card.Unit) : "";
        var outp = new List<SeriesItem>();
        void Add(SeriesItem basis, string caption, SciFigure f) => outp.Add(basis with { Caption = caption, Figure = f });

        switch (id)
        {
            case "density" or "ced" when raw.Count >= 1 && raw[0].X.Length >= 2:
                Add(raw[0], raw[0].Label + " against time, block means and mean ± s.e.", TimeSeries(card, raw[0], unit, blocks));
                foreach (var r in raw.Skip(1)) Add(r, r.Label, Generic(r));
                break;

            case "rdf" when raw.Count >= 1:
            {
                var g = raw[0];
                Add(g, "g(r)", new SciFigure
                {
                    XLabel = "*r* / Å", YLabel = "*g*(*r*)", RefY = 1, YZero = true, Series = [S(g, "*g*(*r*)")],
                    Notes = double.IsFinite(card.Value) ? [$"first peak at {Num(card.Value)} Å"] : [],
                });
                if (raw.Count > 1)
                {
                    var n = raw[1];
                    var rmin = Extra(card, "first minimum");
                    var cnum = Extra(card, "coordination number at the first minimum");
                    Add(n, "running coordination number n(r)", new SciFigure
                    {
                        XLabel = "*r* / Å", YLabel = "*n*(*r*)", YZero = true, Series = [S(n, "*n*(*r*)", SciKind.Line, 1)],
                        Marks = double.IsFinite(rmin) ? [new SciMark(rmin, $"first minimum, *n* = {Num(cnum)}")] : [],
                    });
                }
                break;
            }

            case "sq" or "xray" or "electron" or "neutron" when raw.Count >= 1:
                Add(raw[0], raw[0].Label, Generic(raw[0]) with { Series = [S(raw[0], "")], YZero = true });
                break;

            case "rg" or "ree" when raw.Count >= 1:
            {
                var h = raw[0];
                var sym = id == "rg" ? "*R*_{g}" : "*R*_{ee}";
                Add(h, raw[0].Label + " over all chains and frames", new SciFigure
                {
                    XLabel = sym + " / Å", YLabel = "*P*(" + sym + ") / Å⁻¹", YZero = true,
                    Series = [S(h, "", SciKind.Steps) with { Fill = true, Width = 1.2 }],
                    Marks = double.IsFinite(card.Value) ? [new SciMark(card.Value, $"⟨{sym}²⟩^{{1/2}} = {ValErr(card)} Å")] : [],
                });
                break;
            }

            case "cn" when raw.Count >= 1:
            {
                var s = raw[0];
                var series = new List<SciSeries> { S(s, "*C*_{*n*}", SciKind.LineMarkers) with { Open = true, Width = 1.0 } };
                var xs = s.X.Where(double.IsFinite).ToArray();
                if (double.IsFinite(card.Value) && xs.Length > 1)
                    series.Add(new SciSeries($"*C*_{{∞}} = {Num(card.Value)} (fit in 1/*n*)", [xs.Min(), xs.Max()], [card.Value, card.Value]) { Kind = SciKind.Dashed, Colour = 1, Width = 1.3 });
                Add(s, "C_n against n and C∞", new SciFigure { XLabel = "*n* (backbone bonds)", YLabel = "*C*_{*n*}", Series = [.. series], YZero = true });
                break;
            }

            case "msd" when raw.Count >= 2:
            {
                var atoms = raw[0];
                var com = raw[1];
                var series = new List<SciSeries>
                {
                    S(atoms, "all atoms", SciKind.Line, 1),
                    S(com, "chain centres of mass", SciKind.LineMarkers, 0) with { Width = 1.4 },
                };
                // a guide ∝ t through the last decade of the centres' curve (slope 1 on log axes: diffusive)
                var pts = com.X.Zip(com.Y).Where(p => p.First > 0 && p.Second > 0).ToArray();
                if (pts.Length > 3)
                {
                    var (t1, m1) = pts[^1];
                    var t0 = Math.Max(pts[0].First, t1 / 10);
                    series.Add(new SciSeries("∝ *t* (diffusive)", [t0, t1], [m1 * 0.55 * t0 / t1, m1 * 0.55]) { Kind = SciKind.Dotted, Colour = -1, Width = 1.3 });
                    var ta = pts[0].First;
                    var tb = Math.Min(t1, ta * 10);
                    var ma = pts[0].Second;
                    series.Add(new SciSeries("∝ *t*² (ballistic)", [ta, tb], [ma * 1.8, ma * 1.8 * (tb / ta) * (tb / ta)]) { Kind = SciKind.Dashed, Colour = -1, Width = 1.0 });
                }
                (double, double, string)? window = null;
                if (all.TryGetValue("diffusion", out var d) && Extra(d.Card, "fit from") is var a && Extra(d.Card, "fit to") is var bto && double.IsFinite(a) && double.IsFinite(bto))
                    window = (a, bto, "Einstein fit window");
                Add(com, "MSD of atoms and chain centres (log–log)", new SciFigure { XLabel = "*t* / ps", YLabel = "MSD / Å²", XLog = true, YLog = true, Series = [.. series], XBand = window });
                if (raw.Count > 2)
                {
                    var sl = raw[2];
                    Add(sl, "local slope d ln MSD / d ln t", new SciFigure
                    {
                        XLabel = "*t* / ps", YLabel = "d ln MSD / d ln *t*", XLog = true, RefY = 1,
                        Series = [S(sl, "chain centres", SciKind.LineMarkers) with { Width = 1.2 }], XBand = window,
                    });
                }
                break;
            }

            case "relaxation" when raw.Count >= 1:
            {
                var series = new List<SciSeries>();
                for (var k = 0; k < raw.Count; ++k) series.Add(S(raw[k], raw[k].Label.Replace("P₂(t)", "*P*₂(*t*)").Replace("⟨R(t)·R(0)⟩/⟨R²⟩", "⟨**R**(*t*)·**R**(0)⟩/⟨*R*²⟩"), SciKind.Line, k));
                var tau = Extra(card, "segmental τ_KWW");
                var beta = Extra(card, "segmental β");
                if (double.IsFinite(tau) && double.IsFinite(beta) && beta > 0)
                    series.Add(Kww(raw[0].X, tau, beta, $"KWW, *τ* = {Num(tau)} ps, *β* = {Num(beta)}", 0));
                var te = Extra(card, "end-to-end τ (ps)");
                var be = Extra(card, "end-to-end β");
                if (raw.Count > 1 && double.IsFinite(te) && double.IsFinite(be) && be > 0)
                    series.Add(Kww(raw[1].X, te * be / Gamma(1 / be), be, $"KWW, *β* = {Num(be)}", 1));
                Add(raw[0], "orientational autocorrelations with KWW fits", new SciFigure { XLabel = "*t* / ps", YLabel = "*C*(*t*)", RefY = 0, Series = [.. series] });
                break;
            }

            case "tg" when raw.Count >= 1:
            {
                var v = raw[0];
                var energy = v.YLabel.StartsWith("potential", StringComparison.Ordinal);
                var tg = card.Value;
                var lines = (Tg: tg, Err: card.Error, AtTg: Extra(card, energy ? "potential energy per atom at Tg" : "specific volume at Tg"),
                             SLow: Extra(card, "slope below Tg"), SHigh: Extra(card, "slope above Tg"), Replicas: 1);
                var ok = double.IsFinite(lines.Tg) && double.IsFinite(lines.AtTg) && double.IsFinite(lines.SLow) && double.IsFinite(lines.SHigh);
                var fig = Cooling(v.X.Zip(v.Y).ToArray(), null, ok ? lines : null, energy);
                if (!ok && raw.Count > 1) fig = fig with { Series = [.. fig.Series, S(raw[1], "two-line fit", SciKind.Line, 1)] };
                // the entry keeps the points with the fit through them (the Glass page pools replicas from it)
                if (raw.Count > 1) v = v with { Markers = true, OverlayX = raw[1].X, OverlayY = raw[1].Y };
                Add(v, (energy ? "potential energy" : "specific volume") + " against T, the two fits and Tg", fig);
                if (raw.Count > 2) Add(raw[2], "density against T", new SciFigure
                {
                    XLabel = "*T* / K", YLabel = "*ρ* / g cm⁻³", Series = [S(raw[2], "", SciKind.LineMarkers) with { Width = 1.0 }],
                    Marks = double.IsFinite(tg) ? [new SciMark(tg, "*T*_{g}")] : [],
                });
                break;
            }

            case "tensile_modulus" when raw.Count >= 2:
                Add(raw[1] with { Label = "stress–strain (points) and smoothed", Markers = true, OverlayX = raw[0].X, OverlayY = raw[0].Y }, "stress–strain (points) and smoothed", new SciFigure
                {
                    XLabel = "*ε*", YLabel = "*σ* / MPa",
                    Series = [S(raw[1], "stress per frame", SciKind.Markers, 1) with { Open = true }, S(raw[0], "smoothed ± 0.5 %", SciKind.Line, 0) with { Width = 1.8 }],
                });
                foreach (var r in raw.Skip(2)) Add(r, r.Label, Generic(r));
                break;

            case "vacf" when raw.Count >= 1:
                Add(raw[0], raw[0].Label, Generic(raw[0]) with { RefY = 0, Series = [S(raw[0], "")] });
                foreach (var r in raw.Skip(1)) Add(r, r.Label, Generic(r) with { Series = [S(r, "", SciKind.Line) with { Fill = false }], YZero = true });
                break;

            default:
                foreach (var r in raw) Add(r, r.Label, Generic(r));
                break;
        }
        return outp;
    }

    /// <summary>A cooling scan as a figure: the hold means (± the replicas' SD), the glass and melt lines (solid where they
    /// were fitted, dashed past Tg) and Tg marked.</summary>
    public static SciFigure Cooling((double X, double Y)[] points, double[]? sd, (double Tg, double Err, double AtTg, double SLow, double SHigh, int Replicas)? lines, bool energy)
    {
        var xs = points.Select(p => p.X).ToArray();
        var series = new List<SciSeries>
        {
            new(lines is { Replicas: > 1 } l0 ? $"hold means ± SD ({l0.Replicas} replicas)" : "hold means (NPT)", xs, points.Select(p => p.Y).ToArray())
                { Kind = SciKind.Markers, Colour = 0, Err = sd is { Length: > 0 } && sd.Length == points.Length ? sd : null },
        };
        var marks = new List<SciMark>();
        if (lines is { } l && double.IsFinite(l.Tg) && xs.Length > 1)
        {
            double tmin = xs.Min(), tmax = xs.Max(), ext = 0.3 * (tmax - tmin);
            double Glass(double t) => l.AtTg + l.SLow * (t - l.Tg);
            double Melt(double t) => l.AtTg + l.SHigh * (t - l.Tg);
            var g1 = Math.Min(tmax, l.Tg + ext);
            var m0 = Math.Max(tmin, l.Tg - ext);
            series.Add(new SciSeries("glass fit", [tmin, l.Tg], [Glass(tmin), Glass(l.Tg)]) { Colour = 1, Width = 1.5 });
            series.Add(new SciSeries("", [l.Tg, g1], [Glass(l.Tg), Glass(g1)]) { Kind = SciKind.Dashed, Colour = 1, Width = 1.2, Legend = false });
            series.Add(new SciSeries("melt fit", [l.Tg, tmax], [Melt(l.Tg), Melt(tmax)]) { Colour = 2, Width = 1.5 });
            series.Add(new SciSeries("", [m0, l.Tg], [Melt(m0), Melt(l.Tg)]) { Kind = SciKind.Dashed, Colour = 2, Width = 1.2, Legend = false });
            var dec = double.IsFinite(l.Err) && l.Err > 0 ? Math.Clamp(1 - (int)Math.Floor(Math.Log10(l.Err)), 0, 3) : 0;
            marks.Add(new SciMark(l.Tg, "*T*_{g} = " + l.Tg.ToString("F" + dec, Inv) + (double.IsFinite(l.Err) && l.Err > 0 ? " ± " + l.Err.ToString("F" + dec, Inv) : "") + " K"));
        }
        return new SciFigure { XLabel = "*T* / K", YLabel = energy ? "*E*_{pot} / kcal mol⁻¹" : "*v* / cm³ g⁻¹", Series = [.. series], Marks = [.. marks], Empty = "" };
    }

    /// <summary>The Einstein figure of the diffusion card: the chain-centre MSD on linear axes, the fit window shaded and
    /// the least-squares line through the points in it (its slope / 6 is D).</summary>
    public static SeriesItem? Diffusion(ResultCard d, IReadOnlyDictionary<string, (ResultCard Card, List<SeriesItem> Raw)> all)
    {
        if (!all.TryGetValue("msd", out var m) || m.Raw.Count < 2) return null;
        var com = m.Raw[1];
        double a = Extra(d, "fit from"), b = Extra(d, "fit to");
        var series = new List<SciSeries> { S(com, "chain centres of mass", SciKind.LineMarkers, 0) with { Width = 1.2 } };
        (double, double, string)? window = null;
        if (double.IsFinite(a) && double.IsFinite(b))
        {
            window = (a, b, "fit window");
            var xi = new List<double>();
            var yi = new List<double>();
            for (var k = 0; k < Math.Min(com.X.Length, com.Y.Length); ++k)
                if (com.X[k] >= a - 1e-9 && com.X[k] <= b + 1e-9) { xi.Add(com.X[k]); yi.Add(com.Y[k]); }
            var (c0, slope) = LinFit(xi, yi);
            if (double.IsFinite(slope))
            {
                var tmax = com.X.Where(double.IsFinite).Max();
                series.Add(new SciSeries($"6*Dt* + *c*, *D* = {ValErr(d)} × 10⁻⁵ cm² s⁻¹", [0, tmax], [c0, c0 + slope * tmax]) { Kind = SciKind.Dashed, Colour = 1, Width = 1.5 });
            }
        }
        return com with { Caption = "Einstein fit to the chain-centre MSD", Property = d.Name, LogLog = false, Figure = new SciFigure { XLabel = "*t* / ps", YLabel = "MSD / Å²", YZero = true, Series = [.. series], XBand = window } };
    }
}
