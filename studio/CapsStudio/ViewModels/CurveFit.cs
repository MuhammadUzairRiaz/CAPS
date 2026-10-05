using System.Globalization;

namespace CapsStudio.ViewModels;

/// <summary>Least-squares fits of a curve (D6): linear models solved directly, the others by Levenberg–Marquardt with a
/// numerical Jacobian. Parameter errors from the covariance s² (JᵀJ)⁻¹ at the minimum (s² = SSR / (n − p)), and R².</summary>
public static class CurveFit
{
    public static readonly string[] Models =
    [
        "No fit", "Line  a + b x", "Quadratic  a + b x + c x²", "Power law  a xᵇ", "Exponential decay  a e^(−x/τ) + c",
        "Stretched exponential  a e^(−(x/τ)^β)", "Arrhenius  A e^(−Ea/RT)  (x = T in K)", "Gaussian peaks on a baseline",
    ];

    public sealed record Result(string Model, string[] Names, double[] P, double[] Err, double R2, (double X, double Y)[] Line, string Text);

    /// <summary>Fits model k (index into Models) to the points with lo ≤ x ≤ hi; peaks: the number of Gaussians.</summary>
    public static Result? Fit(int k, double[] x, double[] y, double lo, double hi, int peaks = 1)
    {
        var pts = Enumerable.Range(0, Math.Min(x.Length, y.Length)).Where(i => double.IsFinite(x[i]) && double.IsFinite(y[i]) && x[i] >= lo && x[i] <= hi)
                            .Select(i => (X: x[i], Y: y[i])).ToArray();
        if (k <= 0 || pts.Length < 3) return null;
        var X = pts.Select(p => p.X).ToArray();
        var Y = pts.Select(p => p.Y).ToArray();
        Func<double[], double, double> f;
        string[] names;
        double[] p0;
        switch (k)
        {
            case 1: { var c = Poly(X, Y, 1); f = (p, t) => p[0] + p[1] * t; names = ["a", "b"]; p0 = c; break; }
            case 2: { var c = Poly(X, Y, 2); f = (p, t) => p[0] + p[1] * t + p[2] * t * t; names = ["a", "b", "c"]; p0 = c; break; }
            case 3:
            {
                var ok = pts.Where(q => q.X > 0 && q.Y > 0).ToArray();
                if (ok.Length < 3) return null;
                var c = Poly(ok.Select(q => Math.Log(q.X)).ToArray(), ok.Select(q => Math.Log(q.Y)).ToArray(), 1);
                f = (p, t) => p[0] * Math.Pow(t, p[1]); names = ["a", "b"]; p0 = [Math.Exp(c[0]), c[1]];
                X = ok.Select(q => q.X).ToArray(); Y = ok.Select(q => q.Y).ToArray();
                break;
            }
            case 4:
            {
                var c = Y[^1];
                var a = Y[0] - c;
                var half = Array.FindIndex(Y, v => Math.Abs(v - c) < Math.Abs(a) / Math.E);
                var tau = half > 0 ? X[half] - X[0] : (X[^1] - X[0]) / 3;
                f = (p, t) => p[0] * Math.Exp(-t / p[1]) + p[2]; names = ["a", "τ", "c"]; p0 = [a * Math.Exp(X[0] / Math.Max(1e-12, tau)), Math.Max(1e-9, tau), c];
                break;
            }
            case 5:
            {
                var a = Y[0];
                var e = Array.FindIndex(Y, v => v < a / Math.E);
                var tau = e > 0 ? X[e] : (X[^1] - X[0]) / 2;
                f = (p, t) => p[0] * Math.Exp(-Math.Pow(Math.Max(0, t) / p[1], p[2])); names = ["a", "τ", "β"]; p0 = [a, Math.Max(1e-9, tau), 0.7];
                break;
            }
            case 6:
            {
                var ok = pts.Where(q => q.X > 0 && q.Y > 0).ToArray();
                if (ok.Length < 3) return null;
                var c = Poly(ok.Select(q => 1 / q.X).ToArray(), ok.Select(q => Math.Log(q.Y)).ToArray(), 1);
                f = (p, t) => p[0] * Math.Exp(-p[1] * 1000 / (8.314462618 * t)); names = ["A", "Ea (kJ/mol)"]; p0 = [Math.Exp(c[0]), -c[1] * 8.314462618 / 1000];
                X = ok.Select(q => q.X).ToArray(); Y = ok.Select(q => q.Y).ToArray();
                break;
            }
            case 7:
            {
                var n = Math.Clamp(peaks, 1, 8);
                // starting peaks, greedily: each at the highest point of what the earlier ones leave, its width from where
                // that residual falls to half its height (a narrow peak on a broad halo: the peak first, then the halo)
                var baseLine = Y.Min();
                p0 = new double[1 + 3 * n];
                p0[0] = baseLine;
                names = new string[1 + 3 * n];
                names[0] = "baseline";
                var res = Y.Select(v => v - baseLine).ToArray();
                for (var q = 0; q < n; ++q)
                {
                    var m = 0;
                    for (var i = 1; i < res.Length; ++i) if (res[i] > res[m]) m = i;
                    var h = Math.Max(1e-12, res[m]);
                    int left = m, right = m;
                    while (left > 0 && res[left] > h / 2) --left;
                    while (right < res.Length - 1 && res[right] > h / 2) ++right;
                    var sigma = Math.Max((X[^1] - X[0]) / (4 * X.Length), (X[right] - X[left]) / 2.3548);   // FWHM → σ
                    p0[1 + 3 * q] = h; p0[2 + 3 * q] = X[m]; p0[3 + 3 * q] = sigma;
                    for (var i = 0; i < res.Length; ++i) { var d = (X[i] - X[m]) / sigma; res[i] -= h * Math.Exp(-0.5 * d * d); }
                    names[1 + 3 * q] = $"height {q + 1}"; names[2 + 3 * q] = $"centre {q + 1}"; names[3 + 3 * q] = $"σ {q + 1}";
                }
                f = (p, t) =>
                {
                    var s = p[0];
                    for (var q = 0; q < (p.Length - 1) / 3; ++q) { var d = (t - p[2 + 3 * q]) / Math.Max(1e-12, Math.Abs(p[3 + 3 * q])); s += p[1 + 3 * q] * Math.Exp(-0.5 * d * d); }
                    return s;
                };
                break;
            }
            default: return null;
        }
        var P = k <= 2 ? p0 : Marquardt(f, X, Y, p0);
        var (err, r2) = Errors(f, X, Y, P);
        var x0 = X.Min();
        var x1 = X.Max();
        var line = Enumerable.Range(0, 200).Select(i => x0 + (x1 - x0) * i / 199.0).Select(t => (t, f(P, t))).Where(q => double.IsFinite(q.Item2)).ToArray();
        var inv = CultureInfo.InvariantCulture;
        var text = string.Join(" · ", names.Select((nm, i) => $"{nm} = {Num(P[i])} ± {Num(err[i])}")) + $" · R² = {r2.ToString("0.0000", inv)} · {X.Length} points";
        if (k == 7)
        {
            // areas: √(2π) h σ each, and each peak's share of the peaks' total
            var areas = Enumerable.Range(0, (P.Length - 1) / 3).Select(q => Math.Sqrt(2 * Math.PI) * P[1 + 3 * q] * Math.Abs(P[3 + 3 * q])).ToArray();
            var tot = areas.Sum();
            if (tot > 0) text += " · areas " + string.Join(", ", areas.Select((a, q) => $"{q + 1}: {Num(a)} ({100 * a / tot:0.#} %)"));
        }
        return new Result(Models[k], names, P, err, r2, line, text);
    }

    private static string Num(double v) => !double.IsFinite(v) ? "—" : Math.Abs(v) >= 1e4 || (Math.Abs(v) < 1e-3 && v != 0) ? v.ToString("0.###E+0", CultureInfo.InvariantCulture) : v.ToString("0.####", CultureInfo.InvariantCulture);

    /// <summary>Polynomial of degree d by least squares (normal equations), coefficients from the constant up.</summary>
    public static double[] Poly(double[] x, double[] y, int d)
    {
        var m = d + 1;
        var A = new double[m, m + 1];
        for (var i = 0; i < x.Length; ++i)
        {
            var pw = new double[2 * m];
            pw[0] = 1;
            for (var k = 1; k < 2 * m; ++k) pw[k] = pw[k - 1] * x[i];
            for (var r = 0; r < m; ++r)
            {
                for (var c = 0; c < m; ++c) A[r, c] += pw[r + c];
                A[r, m] += pw[r] * y[i];
            }
        }
        return Solve(A, m);
    }

    private static double[] Solve(double[,] A, int m)
    {
        for (var c = 0; c < m; ++c)
        {
            var p = c;
            for (var r = c + 1; r < m; ++r) if (Math.Abs(A[r, c]) > Math.Abs(A[p, c])) p = r;
            if (Math.Abs(A[p, c]) < 1e-300) return Enumerable.Repeat(double.NaN, m).ToArray();
            for (var j = 0; j <= m; ++j) (A[c, j], A[p, j]) = (A[p, j], A[c, j]);
            for (var r = 0; r < m; ++r)
            {
                if (r == c) continue;
                var f = A[r, c] / A[c, c];
                for (var j = c; j <= m; ++j) A[r, j] -= f * A[c, j];
            }
        }
        return Enumerable.Range(0, m).Select(i => A[i, m] / A[i, i]).ToArray();
    }

    private static double[,] Jacobian(Func<double[], double, double> f, double[] x, double[] p)
    {
        var J = new double[x.Length, p.Length];
        for (var j = 0; j < p.Length; ++j)
        {
            var h = 1e-6 * Math.Max(1e-6, Math.Abs(p[j]));
            var pp = (double[])p.Clone();
            var pm = (double[])p.Clone();
            pp[j] += h; pm[j] -= h;
            for (var i = 0; i < x.Length; ++i) J[i, j] = (f(pp, x[i]) - f(pm, x[i])) / (2 * h);
        }
        return J;
    }

    private static double Ssr(Func<double[], double, double> f, double[] x, double[] y, double[] p)
    {
        var s = 0.0;
        for (var i = 0; i < x.Length; ++i) { var r = y[i] - f(p, x[i]); s += r * r; }
        return double.IsFinite(s) ? s : double.MaxValue;
    }

    /// <summary>Levenberg–Marquardt: (JᵀJ + λ diag JᵀJ) δ = Jᵀr, λ up on a worse step, down on a better one.</summary>
    public static double[] Marquardt(Func<double[], double, double> f, double[] x, double[] y, double[] p0, int iterations = 300)
    {
        var p = (double[])p0.Clone();
        var n = p.Length;
        var lambda = 1e-3;
        var s = Ssr(f, x, y, p);
        for (var it = 0; it < iterations; ++it)
        {
            var J = Jacobian(f, x, p);
            var A = new double[n, n + 1];
            for (var i = 0; i < x.Length; ++i)
            {
                var r = y[i] - f(p, x[i]);
                for (var a = 0; a < n; ++a)
                {
                    for (var b = 0; b < n; ++b) A[a, b] += J[i, a] * J[i, b];
                    A[a, n] += J[i, a] * r;
                }
            }
            var improved = false;
            for (var tries = 0; tries < 12 && !improved; ++tries)
            {
                var B = (double[,])A.Clone();
                for (var a = 0; a < n; ++a) B[a, a] += lambda * Math.Max(1e-12, A[a, a]);
                var d = Solve(B, n);
                if (d.Any(v => !double.IsFinite(v))) { lambda *= 10; continue; }
                var q = p.Zip(d, (u, v) => u + v).ToArray();
                var sq = Ssr(f, x, y, q);
                if (sq < s) { improved = true; var rel = (s - sq) / Math.Max(1e-300, s); p = q; s = sq; lambda = Math.Max(1e-12, lambda / 10); if (rel < 1e-12) return p; }
                else lambda *= 10;
            }
            if (!improved) break;
        }
        return p;
    }

    private static (double[] Err, double R2) Errors(Func<double[], double, double> f, double[] x, double[] y, double[] p)
    {
        var n = p.Length;
        var s = Ssr(f, x, y, p);
        var mean = y.Average();
        var tot = y.Sum(v => (v - mean) * (v - mean));
        var r2 = tot > 0 ? 1 - s / tot : 1;
        var err = Enumerable.Repeat(double.NaN, n).ToArray();
        if (x.Length <= n) return (err, r2);
        var J = Jacobian(f, x, p);
        var s2 = s / (x.Length - n);
        // the inverse of JᵀJ, column by column
        for (var c = 0; c < n; ++c)
        {
            var A = new double[n, n + 1];
            for (var i = 0; i < x.Length; ++i)
                for (var a = 0; a < n; ++a)
                    for (var b = 0; b < n; ++b) A[a, b] += J[i, a] * J[i, b];
            A[c, n] = 1;
            var col = Solve(A, n);
            err[c] = col[c] > 0 ? Math.Sqrt(s2 * col[c]) : double.NaN;
        }
        return (err, r2);
    }
}
