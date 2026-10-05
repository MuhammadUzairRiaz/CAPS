using System.Collections.ObjectModel;
using System.ComponentModel;
using System.Globalization;
using System.Text;
using System.Text.Json.Nodes;

namespace CapsStudio.ViewModels;

/// <summary>A column of a study table: typed values, or a formula over the other columns of its row.</summary>
public sealed class StudyColumn : INotifyPropertyChanged
{
    public event PropertyChangedEventHandler? PropertyChanged;
    private void Raise(string n) => PropertyChanged?.Invoke(this, new PropertyChangedEventArgs(n));
    private string _name = "", _formula = "", _error = "";
    public string Name { get => _name; set { _name = value ?? ""; Raise(nameof(Name)); Raise(nameof(Header)); } }
    public string Formula { get => _formula; set { _formula = value ?? ""; Raise(nameof(Formula)); Raise(nameof(IsFormula)); Raise(nameof(Header)); } }
    public string Error { get => _error; set { _error = value ?? ""; Raise(nameof(Error)); Raise(nameof(HasError)); } }
    public bool HasError => _error.Length > 0;
    public bool IsFormula => _formula.Length > 0;
    public string Header => IsFormula ? $"{_name} = {_formula}" : _name;
}

/// <summary>One cell: the text typed (or the formula's value), and the number it reads as (NaN when it is not one).</summary>
public sealed class StudyCell : INotifyPropertyChanged
{
    public event PropertyChangedEventHandler? PropertyChanged;
    private string _text = "";
    public StudyRow Row { get; init; } = null!;
    public StudyColumn Column { get; init; } = null!;
    public bool ReadOnly => Column.IsFormula;
    public double Number { get; private set; } = double.NaN;
    public string Text
    {
        get => _text;
        set
        {
            if (ReadOnly || (value ?? "") == _text) return;
            Set(value ?? "");
            Row.Table?.Recompute();
        }
    }
    internal void Set(string text)
    {
        _text = text;
        Number = double.TryParse(text.Trim(), NumberStyles.Float, CultureInfo.InvariantCulture, out var v) ? v : double.NaN;
        PropertyChanged?.Invoke(this, new PropertyChangedEventArgs(nameof(Text)));
    }
    internal void SetComputed(double v)
    {
        Number = v;
        _text = double.IsFinite(v) ? StudyStats.N(v) : "";
        PropertyChanged?.Invoke(this, new PropertyChangedEventArgs(nameof(Text)));
    }
}

/// <summary>A row: its label, the structure file it stands for (opened from the table), and one cell per column.</summary>
public sealed class StudyRow
{
    public StudyTableDoc? Table { get; init; }
    public string Label { get; set; } = "";
    public string Path { get; set; } = "";
    public bool HasPath => Path.Length > 0;
    public ObservableCollection<StudyCell> Cells { get; } = new();
}

/// <summary>The table: columns, rows and the recomputation of formula columns (in order, so a formula may use the
/// formula columns to its left).</summary>
public sealed class StudyTableDoc
{
    public ObservableCollection<StudyColumn> Columns { get; } = new();
    public ObservableCollection<StudyRow> Rows { get; } = new();
    private readonly Dictionary<StudyColumn, Formula?> _parsed = new();
    public event Action? Changed;

    public StudyColumn AddColumn(string name, string formula = "")
    {
        var c = new StudyColumn { Name = Unique(name.Trim().Length > 0 ? name.Trim() : "column"), Formula = formula.Trim() };
        Columns.Add(c);
        foreach (var r in Rows) r.Cells.Add(new StudyCell { Row = r, Column = c });
        Recompute();
        return c;
    }
    public string Unique(string name)
    {
        var n = name;
        for (var k = 2; Columns.Any(c => c.Name == n); ++k) n = $"{name} {k}";
        return n;
    }
    public void RemoveColumn(StudyColumn c)
    {
        var i = Columns.IndexOf(c);
        if (i < 0) return;
        Columns.RemoveAt(i);
        foreach (var r in Rows) r.Cells.RemoveAt(i);
        _parsed.Remove(c);
        Recompute();
    }
    public StudyRow AddRow(string label, string path = "", IDictionary<string, double>? values = null, IDictionary<string, string>? texts = null)
    {
        var r = new StudyRow { Table = this, Label = label, Path = path };
        foreach (var c in Columns)
        {
            var cell = new StudyCell { Row = r, Column = c };
            if (values != null && values.TryGetValue(c.Name, out var v)) cell.Set(double.IsFinite(v) ? v.ToString("R", CultureInfo.InvariantCulture) : "");
            else if (texts != null && texts.TryGetValue(c.Name, out var t)) cell.Set(t);
            r.Cells.Add(cell);
        }
        Rows.Add(r);
        return r;
    }
    public double[] Values(StudyColumn c)
    {
        var i = Columns.IndexOf(c);
        return Rows.Select(r => i >= 0 && i < r.Cells.Count ? r.Cells[i].Number : double.NaN).ToArray();
    }
    public void Recompute()
    {
        for (var k = 0; k < Columns.Count; ++k)
        {
            var c = Columns[k];
            if (!c.IsFormula) { c.Error = ""; continue; }
            Formula? f;
            try
            {
                f = Formula.Parse(c.Formula);
                var unknown = f.Names.Where(n => Columns.Take(k).All(x => x.Name != n)).ToList();
                c.Error = unknown.Count == 0 ? "" : $"no column {string.Join(", ", unknown.Select(u => $"'{u}'"))} to its left";
                if (unknown.Count > 0) f = null;
            }
            catch (Exception e) { c.Error = e.Message; f = null; }
            foreach (var r in Rows)
            {
                if (f == null) { r.Cells[k].SetComputed(double.NaN); continue; }
                double V(string name)
                {
                    for (var j = 0; j < k; ++j) if (Columns[j].Name == name) return r.Cells[j].Number;
                    return double.NaN;
                }
                double v;
                try { v = f.Eval(V); } catch { v = double.NaN; }
                r.Cells[k].SetComputed(v);
            }
        }
        Changed?.Invoke();
    }

    public string Csv()
    {
        var sb = new StringBuilder();
        static string Q(string s) => s.Contains(',') || s.Contains('"') ? "\"" + s.Replace("\"", "\"\"") + "\"" : s;
        sb.AppendLine(string.Join(",", new[] { "row" }.Concat(Columns.Select(c => Q(c.Name)))));
        foreach (var r in Rows) sb.AppendLine(string.Join(",", new[] { Q(r.Label) }.Concat(r.Cells.Select(c => Q(c.Text)))));
        return sb.ToString();
    }

    public JsonObject Json() => new()
    {
        ["caps_study_table"] = 1,
        ["columns"] = new JsonArray(Columns.Select(c => (JsonNode)new JsonObject { ["name"] = c.Name, ["formula"] = c.Formula }).ToArray()),
        ["rows"] = new JsonArray(Rows.Select(r => (JsonNode)new JsonObject
        {
            ["label"] = r.Label, ["path"] = r.Path,
            ["cells"] = new JsonArray(r.Cells.Select(c => (JsonNode?)(c.ReadOnly ? null : JsonValue.Create(c.Text))).ToArray()),
        }).ToArray()),
    };

    public void Load(JsonObject o)
    {
        Rows.Clear();
        Columns.Clear();
        _parsed.Clear();
        foreach (var c in (o["columns"] as JsonArray ?? new JsonArray()).OfType<JsonObject>())
            Columns.Add(new StudyColumn { Name = (string?)c["name"] ?? "column", Formula = (string?)c["formula"] ?? "" });
        foreach (var r in (o["rows"] as JsonArray ?? new JsonArray()).OfType<JsonObject>())
        {
            var row = AddRow((string?)r["label"] ?? "", (string?)r["path"] ?? "");
            if (r["cells"] is JsonArray cells)
                for (var k = 0; k < Math.Min(cells.Count, row.Cells.Count); ++k)
                    if (!row.Cells[k].ReadOnly && cells[k] is JsonValue v) row.Cells[k].Set((string?)v ?? "");
        }
        Recompute();
    }

    /// <summary>Rows from CSV text (the first row names the columns when any of it is not a number; a first column of
    /// words labels the rows). Columns are matched by name, new ones added.</summary>
    public int AddCsv(string text, string labelPrefix = "row")
    {
        var lines = text.Replace("\r\n", "\n").Split('\n').Where(l => l.Trim().Length > 0 && !l.TrimStart().StartsWith('#')).ToList();
        if (lines.Count == 0) return 0;
        var sep = lines[0].Contains('\t') && !lines[0].Contains(',') ? '\t' : lines[0].Contains(',') ? ',' : ';';
        List<string> Split(string l)
        {
            var res = new List<string>();
            var sb = new StringBuilder();
            var q = false;
            for (var i = 0; i < l.Length; ++i)
            {
                var ch = l[i];
                if (ch == '"') { if (q && i + 1 < l.Length && l[i + 1] == '"') { sb.Append('"'); ++i; } else q = !q; }
                else if (ch == sep && !q) { res.Add(sb.ToString().Trim()); sb.Clear(); }
                else sb.Append(ch);
            }
            res.Add(sb.ToString().Trim());
            return res;
        }
        static bool IsNum(string s) => double.TryParse(s, NumberStyles.Float, CultureInfo.InvariantCulture, out _);
        var first = Split(lines[0]);
        var header = first.Any(s => s.Length > 0 && !IsNum(s)) && lines.Count > 1;
        var body = lines.Skip(header ? 1 : 0).Select(Split).ToList();
        var width = Math.Max(first.Count, body.Count > 0 ? body.Max(b => b.Count) : 0);
        var labelCol = body.Count > 0 && body.All(b => b.Count > 0 && b[0].Length > 0 && !IsNum(b[0]));
        var names = Enumerable.Range(0, width).Select(k => header && k < first.Count && first[k].Length > 0 ? first[k] : $"column {k + 1}").ToList();
        var start = labelCol ? 1 : 0;
        foreach (var nm in names.Skip(start)) if (Columns.All(c => c.Name != nm)) AddColumn(nm);
        var n0 = Rows.Count;
        foreach (var b in body)
        {
            var texts = new Dictionary<string, string>();
            for (var k = start; k < b.Count && k < names.Count; ++k) texts[names[k]] = b[k];
            AddRow(labelCol ? b[0] : $"{labelPrefix} {Rows.Count + 1}", "", null, texts);
        }
        Recompute();
        return Rows.Count - n0;
    }
}

/// <summary>Analyze › Study table: results gathered from several structures and runs — the open structure and its
/// Analyze results, the Batch table, a table found by the Reader, a CSV — in one sheet with formula columns, a chart of
/// any column against another (with a fit), and statistics: a description of each column, correlations, multiple linear
/// regression and principal components. Saved as .capstable (with the formulas) or CSV.</summary>
public sealed partial class MainViewModel
{
    public bool IsStudy => _module == 74;
    public StudyTableDoc Study { get; } = new();
    private string _studyFile = "", _studyNewName = "", _studyNewFormula = "", _studyStats = "", _studyError = "";
    private int _studyX, _studyY = 1, _studyFit, _studyMode;
    private bool _studyHooked;
    public ObservableCollection<string> StudyColumnNames { get; } = new();
    public string StudyFile => _studyFile.Length > 0 ? RecentFiles.Tilde(_studyFile) : "not saved";
    public string StudyNewName { get => _studyNewName; set => Set(ref _studyNewName, value ?? ""); }
    public string StudyNewFormula { get => _studyNewFormula; set => Set(ref _studyNewFormula, value ?? ""); }
    public int StudyX { get => _studyX; set { if (Set(ref _studyX, Math.Max(0, value))) StudyPlot(); } }
    public int StudyY { get => _studyY; set { if (Set(ref _studyY, Math.Max(0, value))) StudyPlot(); } }
    public int StudyFit { get => _studyFit; set { if (Set(ref _studyFit, Math.Clamp(value, 0, CurveFit.Models.Length - 1))) StudyPlot(); } }
    public static readonly string[] StudyModes = ["Describe each column", "Correlations (Pearson r)", "Regression: y on the ticked columns", "Principal components of the ticked columns"];
    public int StudyMode { get => _studyMode; set { if (Set(ref _studyMode, Math.Clamp(value, 0, 3))) { Raise(nameof(StudyIsPca)); StudyStatsRun(); } } }
    public bool StudyIsPca => _studyMode == 3;
    public string StudyStatsText { get => _studyStats; private set => Set(ref _studyStats, value); }
    public string StudyError { get => _studyError; private set { if (Set(ref _studyError, value ?? "")) Raise(nameof(StudyHasError)); } }
    public bool StudyHasError => _studyError.Length > 0;
    /// <summary>The columns ticked for regression and PCA (by name).</summary>
    public ObservableCollection<StudyPick> StudyPicks { get; } = new();
    public (double X, double Y)[] StudyCurve { get; private set; } = [];
    public (double X, double Y)[] StudyFitLine { get; private set; } = [];
    public string StudyFitText { get; private set; } = "";
    public string StudyXLabel => _studyX < StudyColumnNames.Count ? StudyColumnNames[_studyX] : "";
    public string StudyYLabel => _studyY < StudyColumnNames.Count ? StudyColumnNames[_studyY] : "";
    public string StudySummary => $"{Study.Rows.Count} rows · {Study.Columns.Count} columns ({Study.Columns.Count(c => c.IsFormula)} formulas) · {StudyFile}";
    public event Action? StudyChanged;

    public void OpenStudy()
    {
        if (!_studyHooked)
        {
            _studyHooked = true;
            Study.Changed += StudyRefresh;
        }
        SetModule(74);
        StudyRefresh();
    }

    private void StudyRefresh()
    {
        var (x, y) = (StudyXLabel, StudyYLabel);
        StudyColumnNames.Clear();
        foreach (var c in Study.Columns) StudyColumnNames.Add(c.Name);
        var picked = StudyPicks.Where(p => p.On).Select(p => p.Name).ToHashSet();
        foreach (var p in StudyPicks) p.PropertyChanged -= OnStudyPick;
        StudyPicks.Clear();
        foreach (var c in Study.Columns)
        {
            var p = new StudyPick { Name = c.Name, On = picked.Contains(c.Name) };
            p.PropertyChanged += OnStudyPick;
            StudyPicks.Add(p);
        }
        _studyX = Math.Max(0, StudyColumnNames.IndexOf(x));
        _studyY = StudyColumnNames.IndexOf(y) is var iy && iy >= 0 ? iy : Math.Min(1, Math.Max(0, StudyColumnNames.Count - 1));
        Raise(nameof(StudyX)); Raise(nameof(StudyY)); Raise(nameof(StudySummary));
        StudyPlot();
        StudyStatsRun();
    }
    private void OnStudyPick(object? s, PropertyChangedEventArgs e) => StudyStatsRun();

    public string StudyAddColumn()
    {
        var name = _studyNewName.Trim().Length > 0 ? _studyNewName.Trim() : _studyNewFormula.Trim().Length > 0 ? "f" : "column";
        if (_studyNewFormula.Trim().Length > 0)
            try { Formula.Parse(_studyNewFormula); }
            catch (Exception e) { StudyError = "Formula: " + e.Message; return StudyError; }
        var c = Study.AddColumn(name, _studyNewFormula);
        StudyError = c.Error.Length > 0 ? $"{c.Name}: {c.Error}" : "";
        StudyNewName = "";
        StudyNewFormula = "";
        return $"Column {c.Name} added";
    }
    public void StudyRemoveColumn(StudyColumn c) { Study.RemoveColumn(c); }
    public void StudyAddEmptyRow() { Study.AddRow($"row {Study.Rows.Count + 1}"); StudyRefresh(); }
    public void StudyRemoveRow(StudyRow r) { Study.Rows.Remove(r); Study.Recompute(); }

    /// <summary>The open structure as a row: its size, mass, density and cell, and the Analyze results shown now.</summary>
    public string StudyAddOpen()
    {
        if (_doc == null) return StudyError = "Open a structure first";
        var s = _doc.Summary();
        var v = new Dictionary<string, double>
        {
            ["atoms"] = s.Atoms, ["molecules"] = s.Molecules, ["mass (g/mol)"] = s.TotalMass, ["charge (e)"] = s.TotalCharge,
        };
        if (s.CellValid != 0) { v["density (g/cm³)"] = s.Density; v["volume (Å³)"] = s.Volume; }
        foreach (var r in Analyze.Results.Where(r => r.HasValue))
        {
            var key = r.Unit.Length > 0 ? $"{r.Name} ({r.Unit})" : r.Name;
            v[key] = r.Value;
            if (double.IsFinite(r.Error) && r.Error > 0) v[key + " ±"] = r.Error;
        }
        foreach (var k in v.Keys) if (Study.Columns.All(c => c.Name != k)) Study.AddColumn(k);
        Study.AddRow(Title, _doc.Path ?? "", v);
        Study.Recompute();
        StudyError = "";
        return $"Added {Title}: {v.Count} values";
    }

    /// <summary>The Batch page's inputs (done ones) with their attributes.</summary>
    public string StudyAddBatch()
    {
        var done = BatchInputs.Where(b => b.Attributes.Count > 0).ToList();
        if (done.Count == 0) return StudyError = "No batch results yet: run Analyze › Batch first";
        foreach (var k in done.SelectMany(b => b.Attributes.Keys).Distinct()) if (Study.Columns.All(c => c.Name != k)) Study.AddColumn(k);
        foreach (var b in done) Study.AddRow(Path.GetFileName(Path.GetDirectoryName(b.Path) ?? b.Path) is { Length: > 0 } d ? d : b.Shown, b.Path, b.Attributes);
        Study.Recompute();
        StudyError = "";
        return $"Added {done.Count} batch rows";
    }

    /// <summary>The Reader's chosen table, row by row.</summary>
    public string StudyAddReader()
    {
        if (_rdTable >= ReaderTables.Count) return StudyError = "No table in the Reader: read a log, xvg or CSV first";
        var t = ReaderTables[_rdTable];
        foreach (var c in t.Columns) if (Study.Columns.All(x => x.Name != c)) Study.AddColumn(c);
        var stem = Path.GetFileName(_rdPath);
        var k = 0;
        foreach (var r in t.Rows)
        {
            var v = new Dictionary<string, double>();
            for (var j = 0; j < t.Columns.Length && j < r.Length; ++j) v[t.Columns[j]] = r[j];
            Study.AddRow($"{stem} {++k}", "", v);
        }
        Study.Recompute();
        StudyError = "";
        return $"Added {t.Rows.Count} rows from {stem}";
    }

    public string StudyImportCsv(string path)
    {
        try { var n = Study.AddCsv(File.ReadAllText(path), Path.GetFileNameWithoutExtension(path)); StudyError = ""; return $"Added {n} rows from {Path.GetFileName(path)}"; }
        catch (Exception e) { return StudyError = e.Message; }
    }
    public string StudyPasteCsv(string text)
    {
        var n = Study.AddCsv(text, "pasted");
        return n > 0 ? $"Pasted {n} rows" : StudyError = "Nothing to paste: copy rows from a spreadsheet or CSV";
    }

    public string StudySave(string path)
    {
        try
        {
            if (path.EndsWith(".csv", StringComparison.OrdinalIgnoreCase)) File.WriteAllText(path, Study.Csv());
            else { File.WriteAllText(path, Study.Json().ToJsonString(new System.Text.Json.JsonSerializerOptions { WriteIndented = true })); _studyFile = path; }
            Raise(nameof(StudyFile)); Raise(nameof(StudySummary));
            return $"Saved {Path.GetFileName(path)}";
        }
        catch (Exception e) { return StudyError = e.Message; }
    }
    public string StudyLoad(string path)
    {
        try
        {
            if (path.EndsWith(".csv", StringComparison.OrdinalIgnoreCase) || path.EndsWith(".tsv", StringComparison.OrdinalIgnoreCase)) return StudyImportCsv(path);
            if (System.Text.Json.Nodes.JsonNode.Parse(File.ReadAllText(path)) is not JsonObject o || o["caps_study_table"] == null) return StudyError = "Not a CAPS study table";
            if (!_studyHooked) { _studyHooked = true; Study.Changed += StudyRefresh; }
            Study.Load(o);
            _studyFile = path;
            Raise(nameof(StudyFile));
            if (_module != 74) SetModule(74);
            StudyRefresh();
            return $"Opened {Path.GetFileName(path)}";
        }
        catch (Exception e) { return StudyError = e.Message; }
    }
    public void StudyOpenRow(StudyRow r) { if (r.HasPath && File.Exists(r.Path)) { Open(r.Path); SetModule(8); } }

    private void StudyPlot()
    {
        StudyCurve = [];
        StudyFitLine = [];
        StudyFitText = "";
        if (_studyX < Study.Columns.Count && _studyY < Study.Columns.Count)
        {
            var x = Study.Values(Study.Columns[_studyX]);
            var y = Study.Values(Study.Columns[_studyY]);
            StudyCurve = x.Zip(y).Where(p => double.IsFinite(p.First) && double.IsFinite(p.Second)).OrderBy(p => p.First).ToArray();
            if (_studyFit > 0 && StudyCurve.Length >= 3)
            {
                var fx = StudyCurve.Select(p => p.X).ToArray();
                var f = CurveFit.Fit(_studyFit, fx, StudyCurve.Select(p => p.Y).ToArray(), fx.Min(), fx.Max(), 1);
                if (f != null) { StudyFitLine = f.Line; StudyFitText = f.Text; } else StudyFitText = "the fit did not converge on these points";
            }
        }
        Raise(nameof(StudyFitText)); Raise(nameof(StudyXLabel)); Raise(nameof(StudyYLabel));
        StudyChanged?.Invoke();
    }

    public void StudyStatsRun()
    {
        try
        {
            var cols = Study.Columns.ToList();
            var ticked = cols.Where(c => StudyPicks.FirstOrDefault(p => p.Name == c.Name)?.On == true).ToList();
            var rows = new List<string[]>();
            switch (_studyMode)
            {
                case 0:
                    rows.Add(["column", "n", "mean", "sd", "sem", "min", "median", "max"]);
                    foreach (var c in cols)
                    {
                        var d = StudyStats.Of(Study.Values(c));
                        if (d.Count == 0) continue;
                        rows.Add([c.Name, d.Count.ToString(CultureInfo.InvariantCulture), StudyStats.N(d.Mean), StudyStats.N(d.Sd), StudyStats.N(d.Sem), StudyStats.N(d.Min), StudyStats.N(d.Median), StudyStats.N(d.Max)]);
                    }
                    StudyStatsText = StudyStats.Grid(rows);
                    break;
                case 1:
                {
                    var use = (ticked.Count >= 2 ? ticked : cols).Where(c => StudyStats.Of(Study.Values(c)).Count >= 3).ToList();
                    var tags = use.Select((c, i) => $"[{i + 1}]").ToArray();
                    rows.Add(["", .. tags]);
                    for (var a = 0; a < use.Count; ++a)
                        rows.Add([$"{tags[a]} {use[a].Name}", .. use.Select(b => StudyStats.N(Math.Round(StudyStats.Pearson(Study.Values(use[a]), Study.Values(b)), 4)))]);
                    StudyStatsText = use.Count < 2 ? "Correlations need two numeric columns with three rows or more" : StudyStats.Grid(rows);
                    break;
                }
                case 2:
                {
                    if (_studyY >= cols.Count) { StudyStatsText = "Choose y (the chart's Up column)"; break; }
                    var yc = cols[_studyY];
                    var xs = ticked.Where(c => c != yc).ToList();
                    if (xs.Count == 0) { StudyStatsText = $"Tick the columns that predict {yc.Name}"; break; }
                    var r = global::CapsStudio.ViewModels.StudyStats.Fit(Study.Values(yc), xs.Select(Study.Values).ToList(), xs.Select(c => c.Name).ToList());
                    rows.Add(["term", "coefficient", "std error", "t"]);
                    for (var k = 0; k < r.Coef.Length; ++k) rows.Add([r.Names[k], StudyStats.N(r.Coef[k]), StudyStats.N(r.Se[k]), StudyStats.N(r.T[k])]);
                    StudyStatsText = $"{yc.Name} = b0 + Σ b·x over {r.Rows} rows\n" + StudyStats.Grid(rows) +
                                 $"\nR² {StudyStats.N(r.R2)} · adjusted R² {StudyStats.N(r.AdjR2)} · RMSE {StudyStats.N(r.Rmse)}";
                    break;
                }
                case 3:
                {
                    if (ticked.Count < 2) { StudyStatsText = "Tick two columns or more for the components"; break; }
                    var p = global::CapsStudio.ViewModels.StudyStats.Components(ticked.Select(Study.Values).ToList(), ticked.Select(c => c.Name).ToList());
                    rows.Add(["", .. Enumerable.Range(1, p.Eigen.Length).Select(k => $"PC{k}")]);
                    rows.Add(["eigenvalue", .. p.Eigen.Select(StudyStats.N)]);
                    rows.Add(["explained %", .. p.Explained.Select(e => (100 * e).ToString("0.0", CultureInfo.InvariantCulture))]);
                    for (var k = 0; k < p.Names.Length; ++k) rows.Add([p.Names[k], .. Enumerable.Range(0, p.Eigen.Length).Select(c => StudyStats.N(Math.Round(p.Loadings[k, c], 4)))]);
                    StudyStatsText = $"Standardised columns over {p.Used.Length} rows; loadings below\n" + StudyStats.Grid(rows);
                    break;
                }
            }
            if (_studyError.StartsWith("Statistics", StringComparison.Ordinal)) StudyError = "";
        }
        catch (Exception e) { StudyStatsText = ""; StudyError = "Statistics: " + e.Message; }
    }

    /// <summary>PCA scores (PC1, PC2) as columns of the table, for charting the rows on the components.</summary>
    public string StudyAddScores()
    {
        var ticked = Study.Columns.Where(c => StudyPicks.FirstOrDefault(p => p.Name == c.Name)?.On == true).ToList();
        try
        {
            var p = global::CapsStudio.ViewModels.StudyStats.Components(ticked.Select(Study.Values).ToList(), ticked.Select(c => c.Name).ToList());
            var c1 = Study.AddColumn("PC1");
            var c2 = Study.AddColumn("PC2");
            var (i1, i2) = (Study.Columns.IndexOf(c1), Study.Columns.IndexOf(c2));
            for (var k = 0; k < p.Used.Length; ++k)
            {
                Study.Rows[p.Used[k]].Cells[i1].Set(p.Scores[k, 0].ToString("R", CultureInfo.InvariantCulture));
                Study.Rows[p.Used[k]].Cells[i2].Set(p.Scores[k, 1].ToString("R", CultureInfo.InvariantCulture));
            }
            Study.Recompute();
            return "Added the PC1 and PC2 scores as columns";
        }
        catch (Exception e) { return StudyError = "Statistics: " + e.Message; }
    }
}

/// <summary>A column ticked (or not) for regression and PCA.</summary>
public sealed class StudyPick : INotifyPropertyChanged
{
    public event PropertyChangedEventHandler? PropertyChanged;
    private bool _on;
    public string Name { get; init; } = "";
    public bool On { get => _on; set { _on = value; PropertyChanged?.Invoke(this, new PropertyChangedEventArgs(nameof(On))); } }
}
