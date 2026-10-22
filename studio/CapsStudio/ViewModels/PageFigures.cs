using System.Globalization;
using CapsStudio.Interop;
using CapsStudio.Views;

namespace CapsStudio.ViewModels;

/// <summary>The figures of the module pages (React, Equilibrate, Minimise, Dynamics, Mechanics, Scattering, the DFT and
/// polymer pages, the run monitors …), built from the view model's results so the pages only show them and the self-test
/// can check them: every axis titled "quantity / unit" (or a plain name for a count or a ratio), counts and other
/// nonnegative quantities from zero, times in ps, forces and penalties on log axes.</summary>
public static class PageFigures
{
    private static readonly CultureInfo Inv = CultureInfo.InvariantCulture;
    private static string F(double v, string fmt) => v.ToString(fmt, Inv);

    // ------------------------------------------------------------------ titles used on several pages
    public const string Time = "*t* / ps";
    public const string Temperature = "*T* / K";
    public const string Pressure = "*P* / atm";
    public const string Density = "*ρ* / g cm⁻³";
    public const string Energy = "*E* / kcal mol⁻¹";
    public const string PotEnergy = "*E*_{pot} / kcal mol⁻¹";
    public const string MaxForce = "|*F*|_{max} / kcal mol⁻¹ Å⁻¹";
    public const string Iteration = "iteration";
    public const string Cn = "*C*_{*n*}";

    private static SciSeries S(string name, IReadOnlyList<(double X, double Y)> pts, SciKind kind = SciKind.Line, int colour = 0) =>
        SciSeries.Of(name, pts) with { Kind = kind, Colour = colour };
    /// <summary>A curve with points while it is short (a run's first rows), a plain line once it is long.</summary>
    private static SciKind Auto(int n) => n <= 30 ? SciKind.LineMarkers : SciKind.Line;
    private static double[] Xs(IReadOnlyList<(double X, double Y)> p) => p.Select(q => q.X).ToArray();

    // ------------------------------------------------------------------ React
    /// <summary>Links between chains against the reaction cycle, from (0, 0), with the target as a dashed line.</summary>
    public static SciFigure ReactLinks(IReadOnlyList<CapsReactCycle> rows)
    {
        var pts = rows.Count == 0 ? [] : new[] { (0.0, 0.0) }.Concat(rows.Select(r => ((double)r.Cycle, (double)r.Crosslinks))).ToArray();
        var target = rows.Count > 0 && rows[^1].Target > 0 ? rows[^1].Target : 0;
        return new SciFigure
        {
            XLabel = "reaction cycle", YLabel = "links", XMin = 0, YMin = 0, XWhole = true, YWhole = true, Empty = "",
            Series = pts.Length > 1 ? [S("", pts, Auto(pts.Length)) with { Legend = false }] : [],
            Guides = target > 0 ? [new SciGuide(target, $"target, {target} links")] : [],
        };
    }

    /// <summary>The largest cluster's mass fraction against conversion (both 0 … 1), the simulated gel point (where the
    /// reduced weight-average mass peaks) marked.</summary>
    public static SciFigure ReactGel(IReadOnlyList<CapsReactCycle> rows)
    {
        var pts = rows.Select(r => (r.Conversion, r.LargestFraction)).ToArray();
        SciMark[] marks = [];
        if (rows.Count >= 3)
        {
            var best = rows.MaxBy(r => r.ReducedMw);
            if (best.Cycle != rows[^1].Cycle) marks = [new SciMark(best.Conversion, $"gel point, *α* = {F(best.Conversion, "0.000")}")];
        }
        return new SciFigure
        {
            XLabel = "conversion *α*", YLabel = "largest cluster", XMin = 0, XMax = 1, YMin = 0, YMax = 1, Empty = "",
            Series = pts.Length > 1 ? [S("", pts, Auto(pts.Length)) with { Legend = false }] : [], Marks = marks,
        };
    }

    // ------------------------------------------------------------------ thermo time series (Equilibrate, Dynamics, the analysis dock)
    private static SciStage[] Stages((double From, double To, string Label, bool Shade)[] bands) => bands.Select(b => new SciStage(b.From, b.To, b.Label, b.Shade)).ToArray();

    /// <summary>One thermo quantity against time in ps; the protocol's stages shaded (NPT tinted); block means as points;
    /// a target as a dashed line.</summary>
    public static SciFigure Thermo(IReadOnlyList<CapsThermo> rows, Func<CapsThermo, double> y, string yLabel, (double From, double To, string Label, bool Shade)[]? stages = null,
                                   (double X, double Y)[]? blocks = null, double? target = null, string targetLabel = "target", bool fromZero = false)
    {
        var pts = rows.Where(r => double.IsFinite(y(r))).Select(r => (r.TimePs, y(r))).ToArray();
        var series = new List<SciSeries>();
        if (pts.Length > 1) series.Add(S("", pts, SciKind.Line) with { Width = 1.2, Legend = false });
        if (blocks is { Length: > 0 } && pts.Length > 1) series.Add(S("block means", blocks, SciKind.Markers, 1) with { Square = true });
        // labelled stages: headroom above the curve, so the ensemble names sit over empty plot, not over the data
        double? top = null;
        if (stages is { Length: > 0 } && stages.Any(st => st.Label.Length > 0) && pts.Length > 1 && !fromZero)
        {
            var (lo, hi) = (pts.Min(p => p.Item2), pts.Max(p => p.Item2));
            if (hi > lo) top = hi + 0.22 * (hi - lo);
        }
        return new SciFigure
        {
            XLabel = Time, YLabel = yLabel, Series = [.. series], Stages = stages == null ? [] : Stages(stages), YZero = fromZero, YMax = top,
            Guides = target is double t && double.IsFinite(t) ? [new SciGuide(t, targetLabel)] : [], Empty = "",
        };
    }

    public static SciFigure ThermoDensity(IReadOnlyList<CapsThermo> rows, (double From, double To, string Label, bool Shade)[]? stages = null, (double X, double Y)[]? blocks = null) =>
        Thermo(rows.Where(r => r.Density > 0).ToList(), r => r.Density, Density, stages, blocks);
    public static SciFigure ThermoEnergy(IReadOnlyList<CapsThermo> rows, (double From, double To, string Label, bool Shade)[]? stages = null, (double X, double Y)[]? blocks = null) =>
        Thermo(rows, r => r.Potential, PotEnergy, stages, blocks);
    public static SciFigure ThermoTemperature(IReadOnlyList<CapsThermo> rows, double? target) =>
        Thermo(rows, r => r.Temperature, Temperature, target: target, targetLabel: target is double t ? $"target, {F(t, "0")} K" : "", fromZero: true);
    public static SciFigure ThermoPressure(IReadOnlyList<CapsThermo> rows, double? target) =>
        Thermo(rows, r => r.Pressure, Pressure, target: target, targetLabel: target is double t ? $"target, {F(t, "0.##")} atm" : "");

    /// <summary>⟨Rg⟩ of each production block.</summary>
    public static SciFigure RgBlocks((double X, double Y)[] blocks) => new()
    {
        XLabel = "production block", YLabel = "⟨*R*_{g}⟩ / Å", XMin = 0, XWhole = true, Empty = "",
        Series = blocks.Length > 0 ? [S("", blocks, SciKind.LineMarkers, 0) with { Legend = false }] : [],
    };

    /// <summary>C_n = ⟨R²(n)⟩ / (n⟨b²⟩) against n, with the reference (RIS or a target file) dashed.</summary>
    public static SciFigure ChainCn((double X, double Y)[] curve, (double X, double Y)[] reference, string refName = "reference") => new()
    {
        XLabel = "*n* (backbone bonds apart)", YLabel = Cn, YZero = true, XMin = 0, Empty = "",
        Series = curve.Length > 1 ? [S(reference.Length > 1 ? "this cell" : "", curve), .. reference.Length > 1 ? new[] { S(refName, reference, SciKind.Dashed, 1) } : []] : [],
    };

    // ------------------------------------------------------------------ Minimise
    public static SciFigure RelaxEnergy((double X, double Y)[] e) => new()
    {
        XLabel = Iteration, YLabel = Energy, XMin = 0, Empty = "", Series = e.Length > 1 ? [S("", e, Auto(e.Length)) with { Legend = false }] : [],
    };

    /// <summary>The largest force on a log axis with the tolerance dashed.</summary>
    public static SciFigure RelaxForce((double X, double Y)[] f, double tol) => new()
    {
        XLabel = Iteration, YLabel = MaxForce, YLog = true, XMin = 0, Empty = "",
        Series = f.Length > 1 ? [S("", f.Where(p => p.Y > 0).ToArray(), Auto(f.Length)) with { Legend = false }] : [],
        Guides = tol > 0 ? [new SciGuide(tol, $"tolerance, {F(tol, "0.###")}")] : [],
    };

    // ------------------------------------------------------------------ Pack
    public static SciFigure PackPenalty(IReadOnlyList<(double X, double Y)> c) => new()
    {
        XLabel = "packing round", YLabel = "overlap penalty / arb. units", YLog = true, XMin = 0, XWhole = true, Empty = "",
        Series = c.Count > 1 ? [S("", c.Where(p => p.Y > 0).ToArray(), Auto(c.Count)) with { Legend = false }] : [],
    };

    // ------------------------------------------------------------------ Mechanics
    public static SciFigure Mechanics((double X, double Y)[] curve, (double X, double Y)[] fit, double fitTo) => new()
    {
        XLabel = "engineering strain *ε* / %", YLabel = "true stress *σ* / GPa", XMin = 0, Empty = "",
        Series = curve.Length > 1 ? [S("simulation", curve), .. fit.Length > 1 ? new[] { S("linear fit", fit, SciKind.Dashed, 1) } : []] : [],
        XBand = curve.Length > 1 ? (0, fitTo, "fit window") : null, RefY = 0,
    };

    // ------------------------------------------------------------------ Scattering
    public static SciFigure Scattering((double X, double Y)[] xray, (double X, double Y)[] neutron, (double X, double Y)[] experiment)
    {
        var s = new List<SciSeries>();
        if (xray.Length > 1) s.Add(S("X-ray", xray, SciKind.Line, 0));
        if (neutron.Length > 1) s.Add(S("neutron", neutron, SciKind.Dashed, 1));
        if (experiment.Length > 1 && s.Count > 0)
        {
            var lo = s.Min(x => x.X.Min());
            var hi = s.Max(x => x.X.Max());
            var ex = experiment.Where(p => p.X >= lo && p.X <= hi).ToArray();
            if (ex.Length > 1) s.Add(S("experiment", ex, ex.Length > 120 ? SciKind.Dotted : SciKind.Markers, -1) with { Open = true, Width = 1.4 });
        }
        return new SciFigure { XLabel = "*q* / Å⁻¹", YLabel = "*I*(*q*) / arb. units", YZero = true, Series = [.. s], Empty = "" };
    }

    // ------------------------------------------------------------------ Free volume, chain statistics, charges, histograms
    /// <summary>A histogram (bin centres) drawn as filled steps from zero.</summary>
    public static SciFigure Histogram((double X, double Y)[] h, string x, string y) => new()
    {
        XLabel = x, YLabel = y, YMin = 0, Empty = "",
        Series = h.Length > 1 ? [S("", h, SciKind.Steps) with { Fill = true, Width = 1.2, Legend = false }] : [],
    };

    public static SciFigure PoreSizes((double X, double Y)[] psd) => Histogram(psd, "pore diameter *d* / Å", "*P*(*d*) / Å⁻¹");
    public static SciFigure ChainCnStats((double X, double Y)[] cn) => new()
    {
        XLabel = "*n* (bonds in the sub-chain)", YLabel = Cn, YZero = true, XMin = 0, Empty = "",
        Series = cn.Length > 1 ? [S("", cn, Auto(cn.Length)) with { Legend = false }] : [],
    };
    public static SciFigure EndToEnd((double X, double Y)[] h) => Histogram(h, "*R*_{ee} / Å", "*P*(*R*_{ee}) / Å⁻¹");
    public static SciFigure Charges((double X, double Y)[] h) => Histogram(h, "partial charge *q* / *e*", "atoms");

    // ------------------------------------------------------------------ Interface, orientation, semicrystalline
    public static SciFigure Interface(string axis, (string Label, (double X, double Y)[] Data)[] series, (double From, double To)? gap)
    {
        var film = series.FirstOrDefault(s => s.Label == "film").Data ?? [];
        var surf = series.Where(s => s.Label != "film").ToArray();
        var s = new List<SciSeries>();
        if (film.Length > 1) s.Add(S("film", film));
        for (var k = 0; k < surf.Length && k < 2; ++k) if (surf[k].Data.Length > 1) s.Add(S(surf[k].Label, surf[k].Data, k == 0 ? SciKind.Dashed : SciKind.Dotted, k + 1));
        return new SciFigure
        {
            XLabel = $"*{axis}* / Å (normal to the surface)", YLabel = Density, YMin = 0, Series = [.. s], Empty = "",
            XBand = gap is { } g ? (g.From, g.To, "gap") : null,
        };
    }

    public static SciFigure OrderPerFrame((double X, double Y)[] p) => new()
    {
        XLabel = Time, YLabel = "order parameter *S*", Empty = "",
        Series = p.Length > 1 ? [S("", p, Auto(p.Length)) with { Legend = false }] : [],
    };
    public static SciFigure OrderAlong(string axis, (double X, double Y)[] p) => new()
    {
        XLabel = $"*{axis}* / Å", YLabel = $"⟨*P*₂(cos *θ*_{{*{axis}*}})⟩", RefY = 0, Empty = "",
        Series = p.Length > 1 ? [S("", p, Auto(p.Length)) with { Legend = false }] : [],
    };

    public static SciFigure SemiDensity((double X, double Y)[] p) => new()
    {
        XLabel = "*z* / Å", YLabel = Density, YMin = 0, Empty = "", Series = p.Length > 1 ? [S("", p) with { Legend = false }] : [],
    };
    public static SciFigure SemiOrder((double X, double Y)[] p) => new()
    {
        XLabel = "*z* / Å", YLabel = "⟨*P*₂⟩", RefY = 0, Empty = "", Series = p.Length > 1 ? [S("", p) with { Legend = false }] : [],
    };
    public static SciFigure SemiBridges((double X, double Y)[] h) => new()
    {
        XLabel = "bridge length / CH₂ units", YLabel = "bridges", YMin = 0, YWhole = true, Empty = "",
        Series = h.Length > 0 ? [S("", h, SciKind.Bars) with { Legend = false }] : [],
    };

    // ------------------------------------------------------------------ Thermal conductivity
    /// <summary>Green–Kubo: the running κ of each direction against the correlation time; NEMD: the slab temperatures.</summary>
    public static SciFigure Thermal((double X, double Y)[][] curves, bool gk)
    {
        if (gk)
        {
            string[] names = curves.Length == 3 ? ["*x*", "*y*", "*z*"] : curves.Select((_, k) => $"run {k + 1}").ToArray();
            return new SciFigure
            {
                XLabel = "correlation time *t* / ps", YLabel = "running *κ* / W m⁻¹ K⁻¹", XMin = 0, RefY = 0, Empty = "",
                Series = curves.Select((c, k) => S(curves.Length > 1 ? names[k] : "", c, k == 0 ? SciKind.Line : k == 1 ? SciKind.Dashed : SciKind.Dotted, k)).Where(x => x.X.Length > 1).ToArray(),
            };
        }
        var slab = curves.Length > 0 ? curves[0] : [];
        return new SciFigure
        {
            XLabel = "position along the flux / Å", YLabel = Temperature, Empty = "",
            Series = slab.Length > 1 ? [S("", slab, SciKind.LineMarkers) with { Legend = false }] : [],
        };
    }

    // ------------------------------------------------------------------ DFT
    public static SciFigure DftForce(double[] fmax) => new()
    {
        XLabel = "ionic step", YLabel = "*F*_{max} / eV Å⁻¹", YLog = true, XMin = 0, XWhole = true, Empty = "no forces yet",
        Series = fmax.Length > 0 ? [new SciSeries("", fmax.Select((_, i) => (double)(i + 1)).ToArray(), fmax) { Kind = Auto(fmax.Length), Legend = false }] : [],
        Guides = fmax.Length > 0 ? [new SciGuide(0.01, "0.01 eV Å⁻¹")] : [],
    };
    public static SciFigure DftBinding((double X, double Y)[] bars) => new()
    {
        XLabel = "complex (ranked)", YLabel = "*E*_{bind} / eV", RefY = 0, Empty = "no finished complexes",
        Series = bars.Length > 0 ? [S("", bars, SciKind.Bars) with { Legend = false }] : [],
    };
    public static SciFigure DftDos((double X, double Y)[] dos, string group) => new()
    {
        XLabel = "*E* − *E*_{F} / eV", YLabel = "DOS / states eV⁻¹", YMin = 0, Empty = "no DOSCAR yet",
        Series = dos.Length > 1 ? [S(group, dos) with { Legend = false }] : [], Marks = dos.Length > 1 ? [new SciMark(0, "*E*_{F}")] : [],
    };
    public static SciFigure DftWorkFunction((double X, double Y)[] p) => new()
    {
        XLabel = "*z* / Å", YLabel = "*V* − *E*_{F} / eV", RefY = 0, Empty = "no LOCPOT yet",
        Series = p.Length > 1 ? [S("", p) with { Legend = false }] : [],
    };
    /// <summary>The plane-averaged charge density difference (left axis) and its running integral (right axis).</summary>
    public static SciFigure DftCdd((double X, double Y)[] rho, (double X, double Y)[] q) => new()
    {
        XLabel = "*z* / Å", YLabel = "Δ*ρ* / *e* Å⁻¹", Y2Label = "Δ*Q* / *e*", RefY = 0, Empty = "no CDD set yet",
        Series = rho.Length > 1 ? [S("Δ*ρ*", rho), .. q.Length > 1 ? new[] { S("Δ*Q* (right axis)", q, SciKind.Dashed, 1) with { Right = true } } : []] : [],
    };
    public static SciFigure DftMdHeight((double X, double Y)[] p) => new()
    {
        XLabel = Time, YLabel = "COM height / Å", Empty = "no AIMD segments yet",
        Series = p.Length > 1 ? [S("", p, Auto(p.Length)) with { Legend = false }] : [],
    };

    // ------------------------------------------------------------------ Glass schedule, surface area, recipes, torsion
    /// <summary>The cooling schedule of one replica (time in ps).</summary>
    public static SciFigure GlassSchedule((double X, double Y)[] p) => new()
    {
        XLabel = Time, YLabel = Temperature, XMin = 0, Empty = "", Series = p.Length > 1 ? [S("", p) with { Legend = false }] : [],
    };
    public static SciFigure SasaSeries((double X, double Y)[] p) => new()
    {
        XLabel = "frame", YLabel = "SASA / Å²", Empty = "", Series = p.Length > 1 ? [S("", p, Auto(p.Length)) with { Legend = false }] : [],
    };
    public static SciFigure SasaConvergence((double X, double Y)[] p) => new()
    {
        XLabel = "points per atom", YLabel = "SASA / Å²", Empty = "", Series = p.Length > 1 ? [S("", p, SciKind.LineMarkers) with { Legend = false }] : [],
    };
    public static SciFigure RecipeTemperature((double X, double Y)[] p) => new()
    {
        XLabel = Time, YLabel = Temperature, XMin = 0, Empty = "", Series = p.Length > 1 ? [S("", p) with { Legend = false }] : [],
    };
    /// <summary>The pressure of the NPT steps on a log axis; a step at constant volume leaves a gap.</summary>
    public static SciFigure RecipePressure((double X, double Y)[] p, double tEnd) => new()
    {
        XLabel = Time, YLabel = "*P* / bar", YLog = true, XMin = 0, XMax = tEnd > 0 ? tEnd : null, Empty = "no constant-pressure step",
        Series = p.Count(q => double.IsFinite(q.Y)) > 1 ? [S("", p) with { Legend = false }] : [],
    };
    public static SciFigure Torsion((double X, double Y)[] p) => new()
    {
        XLabel = "dihedral *φ* / deg", YLabel = Energy, Empty = "",
        Series = p.Length > 1 ? [S("", p, SciKind.LineMarkers) with { Legend = false }] : [],
    };

    // ------------------------------------------------------------------ tables with columns named "name (unit)"
    /// <summary>A column name from a table ("Time (ps)", "density (g/cm³)") as an axis title.</summary>
    public static string Title(string column) => AnalyzeFigures.Title(column);

    public static SciFigure Columns(string x, string y, (double X, double Y)[] data, (double X, double Y)[]? fit = null, bool markers = false, string empty = "") => new()
    {
        XLabel = Title(x), YLabel = Title(y), Empty = empty,
        Series = data.Length > 1 || markers && data.Length > 0
            ? [S(fit is { Length: > 1 } ? "data" : "", data, markers ? SciKind.Markers : Auto(data.Length)), .. fit is { Length: > 1 } ? new[] { S("fit", fit, SciKind.Dashed, 1) } : []]
            : [],
    };

    /// <summary>A run monitor's curve (Jobs): a force on a log axis, an iteration or cycle count from zero.</summary>
    public static SciFigure Monitor(string x, string y, (double X, double Y)[] data)
    {
        var log = y.Contains("force", StringComparison.Ordinal);
        var pts = log ? data.Where(p => p.Y > 0).ToArray() : data;
        return new SciFigure
        {
            XLabel = Title(x), YLabel = Title(y), YLog = log, XMin = x.StartsWith("time", StringComparison.Ordinal) || x is "iteration" or "reaction cycle" ? 0 : null, Empty = "",
            Series = pts.Length > 1 ? [S("", pts, Auto(pts.Length)) with { Legend = false }] : [],
        };
    }

    /// <summary>One metric over the inputs of a batch, in input order.</summary>
    public static SciFigure Batch(string metric, (double X, double Y)[] pts) => new()
    {
        XLabel = "input", YLabel = metric.Length > 0 ? Title(metric) : "value", Empty = "",
        Series = pts.Length > 0 ? [S("", pts, SciKind.LineMarkers) with { Legend = false }] : [],
    };

    // the trajectory player's three plots keep the colours of their chips (amber, cyan, red)
    private static readonly int[] TrajColours = [1, 0, 5];
    /// <summary>A quantity of the trajectory against time, the current frame's time marked.</summary>
    public static SciFigure TrajectorySeries(string series, (double X, double Y)[] pts, double? cursor, int k) => new()
    {
        XLabel = Time, YLabel = series.Length > 0 ? Title(series) : "value", Empty = "",
        Series = pts.Length > 1 ? [S("", pts, Auto(pts.Length), TrajColours[k % TrajColours.Length]) with { Legend = false }] : [],
        Marks = cursor is double c && pts.Length > 1 ? [new SciMark(c, "")] : [],
    };

    /// <summary>The conversion against the reaction cycle (the analysis dock's network tab).</summary>
    public static SciFigure ReactConversion(IReadOnlyList<CapsReactCycle> rows)
    {
        var pts = rows.Count == 0 ? [] : new[] { (0.0, 0.0) }.Concat(rows.Select(r => ((double)r.Cycle, r.Conversion))).ToArray();
        return new SciFigure
        {
            XLabel = "reaction cycle", YLabel = "conversion *α*", XMin = 0, YMin = 0, XWhole = true, Empty = "",
            Series = pts.Length > 1 ? [S("", pts, Auto(pts.Length)) with { Legend = false }] : [],
        };
    }

    /// <summary>g(r) of the open structure, from r = 0, with g = 1 dashed.</summary>
    public static SciFigure Rdf((double X, double Y)[] g) => new()
    {
        XLabel = "*r* / Å", YLabel = "*g*(*r*)", XMin = 0, YMin = 0, RefY = 1, Empty = "",
        Series = g.Length > 1 ? [S("", g) with { Legend = false }] : [],
    };

    /// <summary>How far each atom moved between two states (atom order).</summary>
    public static SciFigure Shifts(double[] sh) => new()
    {
        XLabel = "atom", YLabel = "shift / Å", XMin = 0, YMin = 0, Empty = "",
        Series = sh.Length > 1 ? [new SciSeries("", sh.Select((_, i) => (double)(i + 1)).ToArray(), sh) { Kind = SciKind.Line, Width = 1.1, Legend = false }] : [],
    };

    /// <summary>A pipeline table: a scatter of points, a series with its running mean over w rows, or a plain curve.</summary>
    public static SciFigure Table(string x, string y, (double X, double Y)[] pts, bool scatter, int w, int colour)
    {
        var series = new List<SciSeries>();
        if (scatter) series.Add(S("", pts, SciKind.Markers, colour) with { Legend = false });
        else if (w >= 2 && pts.Length > 1)
        {
            series.Add(S("series", pts, SciKind.Line, colour) with { Width = 1 });
            series.Add(S($"running mean over {w}", MainViewModel.RunningMean(pts, w), SciKind.Line, colour == 0 ? 1 : 0) with { Width = 2 });
        }
        else if (pts.Length > 1) series.Add(S("", pts, Auto(pts.Length), colour) with { Legend = false });
        return new SciFigure { XLabel = Title(x), YLabel = Title(y), Series = [.. series], Empty = "" };
    }

    /// <summary>A map over two table columns, coloured by a third.</summary>
    public static SciFigure Heat(string x, string y, string z, double[] xs, double[] ys, double[] zs) => new()
    {
        XLabel = Title(x), YLabel = Title(y), Heat = new SciHeat(xs, ys, zs, Title(z)), Empty = "",
    };

    /// <summary>Two inputs of one table side by side: A solid, B dashed.</summary>
    public static SciFigure Compare(string x, string y, (double X, double Y)[] a, (double X, double Y)[] b) => new()
    {
        XLabel = Title(x), YLabel = Title(y), Empty = "",
        Series = [.. a.Length > 1 ? new[] { S("A", a) } : [], .. b.Length > 1 ? new[] { S("B", b, SciKind.Dashed, 1) } : []],
    };

    // ------------------------------------------------------------------ Coarse-grained (the IBI, mapping and analysis pages)
    public static SciFigure CgKinds(string[] categories, (string Name, (double X, double Y)[] Pts)[] bars) => new()
    {
        Categories = categories, XLabel = "bead type", YLabel = "beads", YMin = 0, YWhole = true, Empty = "map a structure to see its beads",
        Series = bars.Select((b, k) => S(b.Name, b.Pts, SciKind.Bars, k)).ToArray(),
    };

    /// <summary>"r (Å)", "θ (°)", "φ (°)" (the bonded table's variable) as an axis title.</summary>
    public static string CgVariable(string unit) => unit.StartsWith('r') ? "*r* / Å" : unit.StartsWith('θ') ? "*θ* / deg" : unit.StartsWith('φ') ? "*φ* / deg" : Title(unit);

    public static SciFigure CgBondedP(string unit, (double X, double Y)[] p, (double Lo, double Hi, double UMax)? range) => new()
    {
        XLabel = CgVariable(unit), YLabel = "*P* (normalised)", YMin = 0, XMin = range?.Lo, XMax = range?.Hi, Empty = "select a table",
        Series = p.Length > 1 ? [S("", p) with { Legend = false }] : [],
    };
    public static SciFigure CgBondedU(string unit, (double X, double Y)[] u, (double Lo, double Hi, double UMax)? range) => new()
    {
        XLabel = CgVariable(unit), YLabel = "*U* / kcal mol⁻¹", XMin = range?.Lo, XMax = range?.Hi, YMax = range?.UMax, Empty = "select a table",
        Series = u.Length > 1 ? [S("", u, SciKind.Line, 1) with { Legend = false }] : [],
    };
    public static SciFigure CgRdf(IEnumerable<(string Name, (double X, double Y)[] Pts)> curves) => new()
    {
        XLabel = "*r* / Å", YLabel = "*g*(*r*)", RefY = 1, YZero = true, Empty = "make the targets to see g(r)",
        Series = curves.Select((c, k) => S(c.Name, c.Pts, k == 0 ? SciKind.Line : SciKind.Dashed, k)).Where(s => s.X.Length > 1).ToArray(),
    };
    public static SciFigure CgPair((double X, double Y)[] table, string tableName, (double X, double Y)[] fit, string fitName)
    {
        var ys = table.Concat(fit).Select(q => q.Y).Where(double.IsFinite).ToArray();
        return new SciFigure
        {
            XLabel = "*r* / Å", YLabel = "*U* / kcal mol⁻¹", RefY = 0, Empty = "start IBI (or pick a pairs.json) to see the tables",
            // the repulsive core reaches thousands of kcal/mol: the figure shows the well
            YMin = ys.Length > 0 ? Math.Min(-0.1, ys.Min() * 1.3) : null, YMax = ys.Length > 0 ? Math.Max(0.5, -3 * ys.Min()) : null,
            Series = [.. table.Length > 1 ? new[] { S(tableName, table) } : [], .. fit.Length > 1 ? new[] { S(fitName, fit, SciKind.Dashed, 1) } : []],
        };
    }
    public static SciFigure CgR2((double X, double Y)[] built, (double X, double Y)[] aa) => new()
    {
        XLabel = "*n* (bonds along the chain)", YLabel = "⟨*R*²(*n*)⟩ *n*⁻¹ / Å²", XMin = 0, YZero = true, Empty = "build a melt to see its chains",
        Series = [.. built.Length > 1 ? new[] { S("built CG chains", built) } : [], .. aa.Length > 1 ? new[] { S("all-atom (mapped)", aa, SciKind.Dashed, -1) } : []],
    };
    public static SciFigure CgStress(IEnumerable<(string Name, (double X, double Y)[] Pts)> curves) => new()
    {
        XLabel = "engineering strain *ε*", YLabel = "*σ* / MPa", RefY = 0, XMin = 0, Empty = "analyse a stress_strain.dat (written by the tension decks' runs)",
        Series = curves.Select((c, k) => S(c.Name, c.Pts, k == 0 ? SciKind.Line : SciKind.Dashed, k) with { Width = k == 0 ? 2 : 1.3 }).Where(s => s.X.Length > 1).ToArray(),
    };
    /// <summary>The bead mean-square displacements g₁, g₂, g₃ on log–log axes (the table holds log₁₀ t).</summary>
    public static SciFigure CgDynamics(IEnumerable<(string Name, (double X, double Y)[] Pts)> curves) => new()
    {
        XLabel = "*t* / fs", YLabel = "*g*_{*i*}(*t*) / Å²", XLog = true, YLog = true, Empty = "run the dynamics of a CG trajectory",
        Series = curves.Select((c, k) => new SciSeries(c.Name.Replace("g₁", "*g*₁").Replace("g₂", "*g*₂").Replace("g₃", "*g*₃"), c.Pts.Select(p => Math.Pow(10, p.X)).ToArray(), c.Pts.Select(p => p.Y).ToArray()) { Colour = k })
            .Where(s => s.X.Length > 1).ToArray(),
    };

    // ------------------------------------------------------------------ polymer-statistics pages
    public static SciFigure Electrostatics((double X, double Y)[] curve, double rc, double tol, bool pme) => new()
    {
        XLabel = "*r* / Å", YLabel = "erfc(*βr*)", YMin = 0, Empty = "",
        Series = curve.Length > 1 ? [S("", curve) with { Legend = false }] : [],
        Marks = [new SciMark(rc, "*r*_{c}")],
        Points = pme ? [new SciPoint(rc, tol, "", 5)] : [],
    };
    public static SciFigure BlendPhase((double X, double Y)[] binodal, (double X, double Y)[] spinodal, (double X, double Y)? critical, double tNow) => new()
    {
        XLabel = "*φ*_{A} (volume fraction of polymer A)", YLabel = Temperature, XMin = 0, XMax = 1, Empty = "",
        Series = [.. binodal.Length > 1 ? new[] { S("binodal (coexistence)", binodal) } : [], .. spinodal.Length > 1 ? new[] { S("spinodal", spinodal, SciKind.Dashed, 1) } : []],
        Points = critical is { } c ? [new SciPoint(c.X, c.Y, "critical point")] : [],
        Guides = double.IsFinite(tNow) ? [new SciGuide(tNow, $"{F(tNow, "0")} K")] : [],
    };
    public static SciFigure Copolymer((double X, double Y)[] curve, SciPoint[] points) => new()
    {
        XLabel = "*f*₁ (M₁ in the feed)", YLabel = "*F*₁ (M₁ in the chain)", XMin = 0, XMax = 1, YMin = 0, YMax = 1, Empty = "",
        Series = [new SciSeries("*f*₁ = *F*₁", [0, 1], [0, 1]) { Kind = SciKind.Dotted, Colour = -1, Width = 1 }, .. curve.Length > 1 ? new[] { S("Mayo–Lewis", curve) } : []],
        Points = points,
    };
    public static SciFigure Pentads(string[] names, double[] model, string modelName, double[] chain) => new()
    {
        Categories = names, XLabel = "pentad", YLabel = "fraction", YMin = 0, Empty = "",
        Series = [new SciSeries(modelName, model.Select((_, i) => (double)i).ToArray(), model) { Kind = SciKind.Bars, Colour = 1 },
                  new SciSeries("this chain", chain.Select((_, i) => (double)i).ToArray(), chain) { Kind = SciKind.Bars, Colour = 0 }],
    };
    /// <summary>The number and weight distributions of N, with the drawn chains as ticks along the bottom.</summary>
    public static SciFigure Polydispersity((double X, double Y)[] number, (double X, double Y)[] weight, double[] drawn) => new()
    {
        XLabel = "degree of polymerisation *N*", YLabel = "fraction", YMin = 0, Empty = "",
        Series = [S("number fraction", number, SciKind.Line, 0), S("weight fraction", weight, SciKind.Dashed, 1)],
        Marks = drawn.Length <= 40 ? drawn.Select(n => new SciMark(n, "")).ToArray() : [],
        Notes = drawn.Length > 0 ? [$"{drawn.Length} chains drawn{(drawn.Length <= 40 ? " (dashed lines)" : "")}"] : [],
    };

    // ------------------------------------------------------------------ checks (the self-test)
    private static readonly string[] Plain = ["iteration", "reaction cycle", "links", "largest cluster", "packing round", "production block",
                                              "frame", "atom", "points per atom", "ionic step", "complex (ranked)", "bead type", "beads", "pentad", "fraction", "atoms", "bridges", "input",
                                              "*C*_{*n*}", "*g*(*r*)", "erfc(*βr*)", "*P* (normalised)", "order parameter *S*", "⟨*P*₂⟩"];

    /// <summary>What is wrong with a figure (empty: nothing): every axis titled as quantity / unit, a plain name for a count
    /// or a ratio, or a symbol of a dimensionless quantity; round major ticks on linear axes; no monospace "x:" titles.</summary>
    public static string Audit(SciFigure fig, double width = 520, double height = 260)
    {
        var bad = new List<string>();
        bool TitleOk(string t)
        {
            if (t.Length == 0) return false;
            if (t.StartsWith("x:", StringComparison.Ordinal) || t.Contains("log₁₀") || t.Contains(" · ")) return false;
            if (t.Contains(" / ")) return t.Split(" / ")[^1].Trim().Length > 0;
            if (System.Text.RegularExpressions.Regex.IsMatch(t, @"^\*[^*]+\*(_\{[^}]*\})? \(")) return true;   // a dimensionless symbol with what it is: "*φ*_{A} (volume fraction of A)"
            if (Plain.Contains(t) || t.StartsWith("conversion", StringComparison.Ordinal) || t.StartsWith("⟨*P*₂", StringComparison.Ordinal)
                || t.StartsWith("*f*₁", StringComparison.Ordinal) || t.StartsWith("*F*₁", StringComparison.Ordinal) || t.StartsWith("engineering strain", StringComparison.Ordinal)
                || t.StartsWith("*N*", StringComparison.Ordinal) || t.StartsWith("degree of polymerisation", StringComparison.Ordinal)) return true;
            return false;
        }
        if (!TitleOk(fig.XLabel)) bad.Add($"x title '{fig.XLabel}'");
        if (!TitleOk(fig.YLabel)) bad.Add($"y title '{fig.YLabel}'");
        if (fig.Series.Length > 0 || fig.Heat != null)
        {
            var t = SciPlot.Ticks(fig, width, height);
            if (t is { } tk)
            {
                if (fig.Categories.Length == 0 && !fig.XLog && !SciPlot.RoundTicks(tk.X)) bad.Add("x ticks " + string.Join(" ", tk.XText));
                if (!fig.YLog && !SciPlot.RoundTicks(tk.Y)) bad.Add("y ticks " + string.Join(" ", tk.YText));
                if (tk.Y.Length < 2 && !fig.YLog) bad.Add("fewer than two y ticks");
            }
        }
        return string.Join("; ", bad);
    }
}
