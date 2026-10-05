using System.Globalization;
using System.Text.RegularExpressions;

namespace CapsStudio.ViewModels;

/// <summary>A block of numbers found in a text file: its column names and rows (every row as long as the names).</summary>
public sealed class TextTable
{
    public string Title { get; init; } = "";
    public string[] Columns { get; init; } = [];
    public List<double[]> Rows { get; } = new();
    /// <summary>The line the block starts on (1-based), for jumping to it in the text.</summary>
    public int Line { get; init; }
    public override string ToString() => $"{Title} · {Rows.Count} rows";
    public double[] Column(int k) => Rows.Select(r => k < r.Length ? r[k] : double.NaN).ToArray();
}

/// <summary>Numbers in text outputs: LAMMPS thermo blocks in a log (the header row starting with Step, rows until
/// "Loop time" or the next non-number), GROMACS .xvg (legends from the @ lines), CSV/TSV with an optional header row,
/// and whitespace columns under a commented header (# x y …) as CAPS's own .dat files write them.</summary>
public static class TextTables
{
    private static readonly char[] Ws = [' ', '\t'];

    private static bool Num(string s, out double v)
    {
        if (double.TryParse(s, NumberStyles.Float, CultureInfo.InvariantCulture, out v)) return true;
        if (!s.Equals("-nan", StringComparison.OrdinalIgnoreCase)) return false;
        v = double.NaN;
        return true;
    }

    private static bool AllNum(string[] t, out double[] v)
    {
        v = new double[t.Length];
        if (t.Length == 0) return false;
        for (var i = 0; i < t.Length; ++i) if (!Num(t[i], out v[i])) return false;
        return true;
    }

    public static string Kind(string path, string text)
    {
        var name = Path.GetFileName(path).ToLowerInvariant();
        var ext = Path.GetExtension(name);
        if (ext == ".xvg") return "xvg";
        if (ext is ".csv") return "csv";
        if (ext is ".tsv") return "tsv";
        if (name.StartsWith("log.") || ext == ".log" || text.Contains("LAMMPS (", StringComparison.Ordinal)) return "lammps-log";
        return "columns";
    }

    public static List<TextTable> Parse(string path, string text)
    {
        var lines = text.Replace("\r\n", "\n").Split('\n');
        return Kind(path, text) switch
        {
            "xvg" => Xvg(lines),
            "csv" => Delimited(lines, ','),
            "tsv" => Delimited(lines, '\t'),
            "lammps-log" => LammpsLog(lines) is { Count: > 0 } t ? t : Columns(lines),
            _ => Columns(lines),
        };
    }

    // LAMMPS: "Step Temp E_pair …" (no number in it) then rows of as many numbers; a WARNING line inside a run is skipped
    private static List<TextTable> LammpsLog(string[] lines)
    {
        var res = new List<TextTable>();
        for (var i = 0; i < lines.Length; ++i)
        {
            var h = lines[i].Split(Ws, StringSplitOptions.RemoveEmptyEntries);
            if (h.Length < 2 || h[0] != "Step" || h.Any(x => Num(x, out _))) continue;
            var t = new TextTable { Title = $"Run {res.Count + 1} (line {i + 1})", Columns = h, Line = i + 1 };
            var j = i + 1;
            for (; j < lines.Length; ++j)
            {
                var l = lines[j].Trim();
                if (l.StartsWith("WARNING", StringComparison.Ordinal)) continue;
                var tok = l.Split(Ws, StringSplitOptions.RemoveEmptyEntries);
                if (tok.Length != h.Length || !AllNum(tok, out var v)) break;
                t.Rows.Add(v);
            }
            if (t.Rows.Count > 0) res.Add(t);
            i = j - 1;
        }
        return res;
    }

    // GROMACS xvg: '@ xaxis label "Time (ps)"', '@ s0 legend "Potential"'; comments '#'
    private static List<TextTable> Xvg(string[] lines)
    {
        var legends = new SortedDictionary<int, string>();
        string x = "x", title = "";
        var rows = new List<double[]>();
        var first = 0;
        var rx = new Regex("^@\\s*s(\\d+)\\s+legend\\s+\"(.*)\"");
        for (var i = 0; i < lines.Length; ++i)
        {
            var l = lines[i].Trim();
            if (l.Length == 0 || l[0] == '#') continue;
            if (l[0] == '@')
            {
                var m = rx.Match(l);
                if (m.Success) legends[int.Parse(m.Groups[1].Value, CultureInfo.InvariantCulture)] = m.Groups[2].Value;
                else if (Regex.Match(l, "^@\\s*xaxis\\s+label\\s+\"(.*)\"") is { Success: true } xm) x = xm.Groups[1].Value;
                else if (Regex.Match(l, "^@\\s*title\\s+\"(.*)\"") is { Success: true } tm) title = tm.Groups[1].Value;
                continue;
            }
            if (l[0] == '&') break;   // the next data set: the first is read
            if (AllNum(l.Split(Ws, StringSplitOptions.RemoveEmptyEntries), out var v)) { if (rows.Count == 0) first = i + 1; rows.Add(v); }
        }
        if (rows.Count == 0) return [];
        var n = rows.Max(r => r.Length);
        var cols = new string[n];
        cols[0] = x;
        for (var k = 1; k < n; ++k) cols[k] = legends.TryGetValue(k - 1, out var s) ? s : $"column {k + 1}";
        var t = new TextTable { Title = title.Length > 0 ? title : "Data", Columns = cols, Line = first };
        t.Rows.AddRange(rows.Where(r => r.Length == n));
        return [t];
    }

    private static List<TextTable> Delimited(string[] lines, char sep)
    {
        string[]? head = null;
        var t0 = 0;
        var rows = new List<double[]>();
        for (var i = 0; i < lines.Length; ++i)
        {
            var l = lines[i].Trim();
            if (l.Length == 0 || l[0] == '#') continue;
            var tok = l.Split(sep).Select(s => s.Trim().Trim('"')).ToArray();
            if (AllNum(tok, out var v)) { if (rows.Count == 0) t0 = i + 1; if (head == null || v.Length == head.Length) rows.Add(v); }
            else if (rows.Count == 0) head = tok;
        }
        if (rows.Count == 0) return [];
        var n = rows[0].Length;
        var t = new TextTable { Title = "Table", Columns = head is { } h && h.Length == n ? h : Enumerable.Range(1, n).Select(k => $"column {k}").ToArray(), Line = t0 };
        t.Rows.AddRange(rows.Where(r => r.Length == n));
        return [t];
    }

    // whitespace blocks: runs of rows with the same count of numbers (at least 3 rows, 2 columns); the names from the
    // comment or word line just above, when it has as many words
    private static List<TextTable> Columns(string[] lines)
    {
        var res = new List<TextTable>();
        for (var i = 0; i < lines.Length; ++i)
        {
            var tok = lines[i].Split(Ws, StringSplitOptions.RemoveEmptyEntries);
            if (tok.Length < 2 || !AllNum(tok, out var v0)) continue;
            var t = new List<double[]> { v0 };
            var j = i + 1;
            for (; j < lines.Length; ++j)
            {
                var tk = lines[j].Split(Ws, StringSplitOptions.RemoveEmptyEntries);
                if (tk.Length != tok.Length || !AllNum(tk, out var v)) break;
                t.Add(v);
            }
            if (t.Count >= 3)
            {
                string[] cols = Enumerable.Range(1, tok.Length).Select(k => $"column {k}").ToArray();
                for (var k = i - 1; k >= 0 && k >= i - 2; --k)
                {
                    var w = lines[k].TrimStart('#', ' ', '\t').Split(Ws, StringSplitOptions.RemoveEmptyEntries);
                    if (w.Length == tok.Length && !w.Any(x => Num(x, out _))) { cols = w; break; }
                }
                var tt = new TextTable { Title = $"Block {res.Count + 1} (line {i + 1})", Columns = cols, Line = i + 1 };
                tt.Rows.AddRange(t);
                res.Add(tt);
            }
            i = j - 1;
        }
        return res;
    }

    /// <summary>Mean, its error from five blocks, min, max and last of y over the rows from a fraction of the way in.</summary>
    public static (double Mean, double Err, double Min, double Max, double Last, int N) Stats(double[] y, double from)
    {
        var s = y.Skip((int)Math.Floor(Math.Clamp(from, 0, 0.99) * y.Length)).Where(double.IsFinite).ToArray();
        if (s.Length == 0) return (double.NaN, double.NaN, double.NaN, double.NaN, double.NaN, 0);
        var mean = s.Average();
        var err = double.NaN;
        if (s.Length >= 10)
        {
            var nb = 5;
            var per = s.Length / nb;
            var bm = Enumerable.Range(0, nb).Select(b => s.Skip(b * per).Take(per).Average()).ToArray();
            var m = bm.Average();
            err = Math.Sqrt(bm.Sum(x => (x - m) * (x - m)) / (nb - 1) / nb);
        }
        return (mean, err, s.Min(), s.Max(), s[^1], s.Length);
    }
}
