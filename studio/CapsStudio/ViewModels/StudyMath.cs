using System.Globalization;
using System.Text;

namespace CapsStudio.ViewModels;

/// <summary>Formula columns of a study table: arithmetic over other columns of the same row. Names are bare words
/// (density, Tg_K) or anything in square brackets ([E (GPa)]); + − × ÷ ^, unary minus, parentheses; functions sqrt,
/// ln, log10, log, exp, abs, sin, cos, tan, asin, acos, atan, min, max, pow, round; constants pi, e, R (8.314462618
/// J/(mol K)), kB (1.380649e-23 J/K), NA (6.02214076e23 /mol). A missing or empty cell gives NaN for the row.</summary>
public sealed class Formula
{
    private readonly Func<Func<string, double>, double> _f;
    public IReadOnlyList<string> Names { get; }
    public string Text { get; }
    private Formula(string text, Func<Func<string, double>, double> f, List<string> names) { Text = text; _f = f; Names = names; }
    public double Eval(Func<string, double> value) => _f(value);

    private static readonly Dictionary<string, double> Constants = new()
    {
        ["pi"] = Math.PI, ["e"] = Math.E, ["R"] = 8.314462618, ["kB"] = 1.380649e-23, ["NA"] = 6.02214076e23,
    };

    /// <summary>The formula, or an exception with what is wrong and where.</summary>
    public static Formula Parse(string text)
    {
        var p = new Parser(text);
        var f = p.Sum();
        p.Space();
        if (p.I < text.Length) throw new FormatException($"unexpected '{text[p.I]}' at {p.I + 1}");
        return new Formula(text, f, p.Names.Distinct().ToList());
    }

    private sealed class Parser(string s)
    {
        public int I;
        public readonly List<string> Names = new();
        public void Space() { while (I < s.Length && char.IsWhiteSpace(s[I])) ++I; }
        private bool Eat(char c) { Space(); if (I < s.Length && s[I] == c) { ++I; return true; } return false; }

        public Func<Func<string, double>, double> Sum()
        {
            var a = Product();
            while (true)
            {
                if (Eat('+')) { var l = a; var r = Product(); a = v => l(v) + r(v); }
                else if (Eat('-') || Eat('−')) { var l = a; var r = Product(); a = v => l(v) - r(v); }
                else return a;
            }
        }
        private Func<Func<string, double>, double> Product()
        {
            var a = Power();
            while (true)
            {
                if (Eat('*') || Eat('×') || Eat('·')) { var l = a; var r = Power(); a = v => l(v) * r(v); }
                else if (Eat('/') || Eat('÷')) { var l = a; var r = Power(); a = v => l(v) / r(v); }
                else return a;
            }
        }
        private Func<Func<string, double>, double> Power()
        {
            var a = Unary();
            if (Eat('^')) { var r = Power(); var l = a; return v => Math.Pow(l(v), r(v)); }   // right-associative
            return a;
        }
        private Func<Func<string, double>, double> Unary()
        {
            if (Eat('-') || Eat('−')) { var a = Unary(); return v => -a(v); }
            if (Eat('+')) return Unary();
            return Atom();
        }
        private Func<Func<string, double>, double> Atom()
        {
            Space();
            if (I >= s.Length) throw new FormatException("the formula ends too early");
            var c = s[I];
            if (c == '(')
            {
                ++I;
                var a = Sum();
                if (!Eat(')')) throw new FormatException($"a ')' is missing at {I + 1}");
                return a;
            }
            if (c == '[')
            {
                var end = s.IndexOf(']', I + 1);
                if (end < 0) throw new FormatException("a ']' is missing");
                var name = s[(I + 1)..end].Trim();
                I = end + 1;
                Names.Add(name);
                return v => v(name);
            }
            if (char.IsDigit(c) || c == '.')
            {
                var st = I;
                while (I < s.Length && (char.IsDigit(s[I]) || s[I] == '.')) ++I;
                if (I < s.Length && (s[I] == 'e' || s[I] == 'E') && I + 1 < s.Length && (char.IsDigit(s[I + 1]) || ((s[I + 1] == '-' || s[I + 1] == '+') && I + 2 < s.Length && char.IsDigit(s[I + 2]))))
                {
                    I += 2;
                    while (I < s.Length && char.IsDigit(s[I])) ++I;
                }
                var x = double.Parse(s[st..I], NumberStyles.Float, CultureInfo.InvariantCulture);
                return _ => x;
            }
            if (char.IsLetter(c) || c == '_')
            {
                var st = I;
                while (I < s.Length && (char.IsLetterOrDigit(s[I]) || s[I] == '_')) ++I;
                var w = s[st..I];
                if (Eat('('))
                {
                    var args = new List<Func<Func<string, double>, double>>();
                    if (!Eat(')'))
                    {
                        do args.Add(Sum()); while (Eat(','));
                        if (!Eat(')')) throw new FormatException($"a ')' is missing after the arguments of {w}");
                    }
                    return Function(w, args);
                }
                if (Constants.TryGetValue(w, out var k)) return _ => k;
                Names.Add(w);
                return v => v(w);
            }
            throw new FormatException($"unexpected '{c}' at {I + 1}");
        }
        private static Func<Func<string, double>, double> Function(string w, List<Func<Func<string, double>, double>> a)
        {
            void Need(int n) { if (a.Count != n) throw new FormatException($"{w} takes {n} argument{(n == 1 ? "" : "s")}"); }
            Func<double, double>? one = w switch
            {
                "sqrt" => Math.Sqrt, "ln" => Math.Log, "log10" => Math.Log10, "log" => Math.Log10, "exp" => Math.Exp, "abs" => Math.Abs,
                "sin" => Math.Sin, "cos" => Math.Cos, "tan" => Math.Tan, "asin" => Math.Asin, "acos" => Math.Acos, "atan" => Math.Atan,
                "round" => x => Math.Round(x, MidpointRounding.AwayFromZero), _ => null,
            };
            if (one != null) { Need(1); var f = a[0]; return v => one(f(v)); }
            switch (w)
            {
                case "pow": { Need(2); var (x, y) = (a[0], a[1]); return v => Math.Pow(x(v), y(v)); }
                case "min": if (a.Count == 0) throw new FormatException("min takes arguments"); return v => a.Min(f => f(v));
                case "max": if (a.Count == 0) throw new FormatException("max takes arguments"); return v => a.Max(f => f(v));
            }
            throw new FormatException($"no function '{w}'");
        }
    }
}

/// <summary>Statistics over the columns of a study table (rows with a missing value in any chosen column are left out).</summary>
public static class StudyStats
{
    private static readonly CultureInfo Inv = CultureInfo.InvariantCulture;
    public static string N(double v) => !double.IsFinite(v) ? "—" : v == 0 ? "0" : Math.Abs(v) >= 1e5 || Math.Abs(v) < 1e-3 ? v.ToString("0.###E+0", Inv) : v.ToString("0.####", Inv);

    /// <summary>The rows (as indices) where every column has a finite value.</summary>
    public static int[] Complete(IReadOnlyList<double[]> cols)
    {
        var n = cols.Count == 0 ? 0 : cols.Min(c => c.Length);
        return Enumerable.Range(0, n).Where(i => cols.All(c => double.IsFinite(c[i]))).ToArray();
    }

    public sealed record Describe(int Count, double Mean, double Sd, double Sem, double Min, double Median, double Max);
    public static Describe Of(double[] x)
    {
        var v = x.Where(double.IsFinite).OrderBy(t => t).ToArray();
        if (v.Length == 0) return new Describe(0, double.NaN, double.NaN, double.NaN, double.NaN, double.NaN, double.NaN);
        var m = v.Average();
        var sd = v.Length > 1 ? Math.Sqrt(v.Sum(t => (t - m) * (t - m)) / (v.Length - 1)) : double.NaN;
        var med = v.Length % 2 == 1 ? v[v.Length / 2] : 0.5 * (v[v.Length / 2 - 1] + v[v.Length / 2]);
        return new Describe(v.Length, m, sd, sd / Math.Sqrt(v.Length), v[0], med, v[^1]);
    }

    /// <summary>Pearson r over the rows both have.</summary>
    public static double Pearson(double[] a, double[] b)
    {
        var rows = Complete([a, b]);
        if (rows.Length < 3) return double.NaN;
        double ma = rows.Average(i => a[i]), mb = rows.Average(i => b[i]), sab = 0, saa = 0, sbb = 0;
        foreach (var i in rows) { sab += (a[i] - ma) * (b[i] - mb); saa += (a[i] - ma) * (a[i] - ma); sbb += (b[i] - mb) * (b[i] - mb); }
        return saa > 0 && sbb > 0 ? sab / Math.Sqrt(saa * sbb) : double.NaN;
    }

    public sealed record Regression(string[] Names, double[] Coef, double[] Se, double[] T, double R2, double AdjR2, double Rmse, int Rows, double[] Fitted, int[] Used);

    /// <summary>Least squares y = b0 + Σ b_k x_k (Householder QR); standard errors from σ² (XᵀX)⁻¹ with σ² = RSS/(n − p).</summary>
    public static Regression Fit(double[] y, IReadOnlyList<double[]> xs, IReadOnlyList<string> names)
    {
        var rows = Complete([y, .. xs]);
        var p = xs.Count + 1;
        var n = rows.Length;
        if (n <= p) throw new InvalidOperationException($"{n} complete rows for {p} coefficients: the fit needs more rows than coefficients");
        var X = new double[n, p];
        var Y = new double[n];
        for (var r = 0; r < n; ++r)
        {
            X[r, 0] = 1;
            for (var k = 0; k < xs.Count; ++k) X[r, k + 1] = xs[k][rows[r]];
            Y[r] = y[rows[r]];
        }
        var A = (double[,])X.Clone();
        var b = (double[])Y.Clone();
        var Rd = new double[p];
        for (var k = 0; k < p; ++k)   // Householder: A = QR, b ← Qᵀb
        {
            double norm = 0;
            for (var i = k; i < n; ++i) norm += A[i, k] * A[i, k];
            norm = Math.Sqrt(norm);
            if (norm < 1e-300) throw new InvalidOperationException("the columns are not independent (one is constant or a combination of others)");
            var alpha = A[k, k] > 0 ? -norm : norm;
            var v = new double[n];
            for (var i = k; i < n; ++i) v[i] = A[i, k];
            v[k] -= alpha;
            double vv = 0;
            for (var i = k; i < n; ++i) vv += v[i] * v[i];
            if (vv > 0)
            {
                for (var j = k; j < p; ++j)
                {
                    double s = 0;
                    for (var i = k; i < n; ++i) s += v[i] * A[i, j];
                    s = 2 * s / vv;
                    for (var i = k; i < n; ++i) A[i, j] -= s * v[i];
                }
                double sb = 0;
                for (var i = k; i < n; ++i) sb += v[i] * b[i];
                sb = 2 * sb / vv;
                for (var i = k; i < n; ++i) b[i] -= sb * v[i];
            }
            Rd[k] = A[k, k];
        }
        var scale = Enumerable.Range(0, p).Max(k => Math.Abs(Rd[k]));
        for (var k = 0; k < p; ++k)
            if (Math.Abs(Rd[k]) < 1e-10 * scale) throw new InvalidOperationException("the columns are not independent (one is constant or a combination of others)");
        var coef = new double[p];
        for (var k = p - 1; k >= 0; --k)
        {
            var s = b[k];
            for (var j = k + 1; j < p; ++j) s -= A[k, j] * coef[j];
            coef[k] = s / A[k, k];
        }
        // R⁻¹ for (XᵀX)⁻¹ = R⁻¹ R⁻ᵀ
        var Ri = new double[p, p];
        for (var c = 0; c < p; ++c)
        {
            Ri[c, c] = 1 / A[c, c];
            for (var r = c - 1; r >= 0; --r)
            {
                double s = 0;
                for (var j = r + 1; j <= c; ++j) s += A[r, j] * Ri[j, c];
                Ri[r, c] = -s / A[r, r];
            }
        }
        var fitted = new double[n];
        double rss = 0, my = Y.Average(), tss = 0;
        for (var r = 0; r < n; ++r)
        {
            double f = 0;
            for (var k = 0; k < p; ++k) f += X[r, k] * coef[k];
            fitted[r] = f;
            rss += (Y[r] - f) * (Y[r] - f);
            tss += (Y[r] - my) * (Y[r] - my);
        }
        var s2 = rss / (n - p);
        var se = new double[p];
        for (var k = 0; k < p; ++k)
        {
            double s = 0;
            for (var j = k; j < p; ++j) s += Ri[k, j] * Ri[k, j];
            se[k] = Math.Sqrt(s2 * s);
        }
        var r2 = tss > 0 ? 1 - rss / tss : double.NaN;
        var adj = tss > 0 ? 1 - (rss / (n - p)) / (tss / (n - 1)) : double.NaN;
        return new Regression(["intercept", .. names], coef, se, coef.Zip(se, (c, e) => e > 0 ? c / e : double.NaN).ToArray(), r2, adj, Math.Sqrt(rss / n), n, fitted, rows);
    }

    public sealed record Pca(string[] Names, double[] Eigen, double[] Explained, double[,] Loadings, double[,] Scores, int[] Used);

    /// <summary>PCA of the standardised columns: eigenvectors of the correlation matrix (cyclic Jacobi), largest first;
    /// scores are the standardised rows projected onto them.</summary>
    public static Pca Components(IReadOnlyList<double[]> cols, IReadOnlyList<string> names)
    {
        var rows = Complete(cols);
        var m = cols.Count;
        var n = rows.Length;
        if (m < 2 || n < 3) throw new InvalidOperationException("PCA needs two columns or more and three complete rows or more");
        var Z = new double[n, m];
        for (var k = 0; k < m; ++k)
        {
            var v = rows.Select(i => cols[k][i]).ToArray();
            var mean = v.Average();
            var sd = Math.Sqrt(v.Sum(t => (t - mean) * (t - mean)) / (n - 1));
            if (sd <= 0) throw new InvalidOperationException($"'{names[k]}' is constant");
            for (var r = 0; r < n; ++r) Z[r, k] = (v[r] - mean) / sd;
        }
        var C = new double[m, m];
        for (var a = 0; a < m; ++a)
            for (var b = 0; b < m; ++b)
            {
                double s = 0;
                for (var r = 0; r < n; ++r) s += Z[r, a] * Z[r, b];
                C[a, b] = s / (n - 1);
            }
        var (w, V) = Jacobi(C);
        var order = Enumerable.Range(0, m).OrderByDescending(k => w[k]).ToArray();
        var eig = order.Select(k => Math.Max(0, w[k])).ToArray();
        var tot = eig.Sum();
        var L = new double[m, m];
        for (var c = 0; c < m; ++c)
        {
            var sign = 0.0;   // a fixed sign: the largest loading positive
            var big = 0.0;
            for (var k = 0; k < m; ++k) if (Math.Abs(V[k, order[c]]) > big) { big = Math.Abs(V[k, order[c]]); sign = Math.Sign(V[k, order[c]]); }
            for (var k = 0; k < m; ++k) L[k, c] = sign * V[k, order[c]];
        }
        var S = new double[n, m];
        for (var r = 0; r < n; ++r)
            for (var c = 0; c < m; ++c)
            {
                double s = 0;
                for (var k = 0; k < m; ++k) s += Z[r, k] * L[k, c];
                S[r, c] = s;
            }
        return new Pca(names.ToArray(), eig, eig.Select(e => tot > 0 ? e / tot : 0).ToArray(), L, S, rows);
    }

    /// <summary>Eigenvalues and eigenvectors (columns) of a symmetric matrix by cyclic Jacobi rotations.</summary>
    public static (double[] W, double[,] V) Jacobi(double[,] A0)
    {
        var m = A0.GetLength(0);
        var A = (double[,])A0.Clone();
        var V = new double[m, m];
        for (var i = 0; i < m; ++i) V[i, i] = 1;
        for (var sweep = 0; sweep < 100; ++sweep)
        {
            double off = 0;
            for (var p = 0; p < m; ++p) for (var q = p + 1; q < m; ++q) off += A[p, q] * A[p, q];
            if (off < 1e-22) break;
            for (var p = 0; p < m; ++p)
                for (var q = p + 1; q < m; ++q)
                {
                    if (Math.Abs(A[p, q]) < 1e-300) continue;
                    var theta = (A[q, q] - A[p, p]) / (2 * A[p, q]);
                    var t = Math.Sign(theta == 0 ? 1 : theta) / (Math.Abs(theta) + Math.Sqrt(theta * theta + 1));
                    var c = 1 / Math.Sqrt(t * t + 1);
                    var s = t * c;
                    for (var k = 0; k < m; ++k)
                    {
                        var akp = A[k, p]; var akq = A[k, q];
                        A[k, p] = c * akp - s * akq; A[k, q] = s * akp + c * akq;
                    }
                    for (var k = 0; k < m; ++k)
                    {
                        var apk = A[p, k]; var aqk = A[q, k];
                        A[p, k] = c * apk - s * aqk; A[q, k] = s * apk + c * aqk;
                    }
                    for (var k = 0; k < m; ++k)
                    {
                        var vkp = V[k, p]; var vkq = V[k, q];
                        V[k, p] = c * vkp - s * vkq; V[k, q] = s * vkp + c * vkq;
                    }
                }
        }
        var w = new double[m];
        for (var i = 0; i < m; ++i) w[i] = A[i, i];
        return (w, V);
    }

    /// <summary>A text table with aligned columns (the statistics panel).</summary>
    public static string Grid(IReadOnlyList<string[]> rows)
    {
        if (rows.Count == 0) return "";
        var w = Enumerable.Range(0, rows.Max(r => r.Length)).Select(c => rows.Max(r => c < r.Length ? r[c].Length : 0)).ToArray();
        var sb = new StringBuilder();
        foreach (var r in rows) sb.AppendLine(string.Join("  ", r.Select((t, c) => c == 0 ? t.PadRight(w[c]) : t.PadLeft(w[c]))).TrimEnd());
        return sb.ToString().TrimEnd();
    }
}
