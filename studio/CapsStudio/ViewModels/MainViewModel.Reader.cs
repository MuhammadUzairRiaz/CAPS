using System.Collections.ObjectModel;
using System.Globalization;

namespace CapsStudio.ViewModels;

/// <summary>Studio › Reader: a text output read in CAPS — LAMMPS logs and inputs, GROMACS .mdp and .xvg, CSV, notes.
/// The text with a find bar and a follow mode for a log still being written; the numbers in it (thermo blocks, xvg
/// columns, CSV) as tables, any column drawn against another, with its mean over the tail.</summary>
public sealed partial class MainViewModel
{
    public bool IsReader => _module == 73;
    private static readonly string[] ReaderExts = [".log", ".mdp", ".in", ".txt", ".md", ".csv", ".tsv", ".xvg", ".out"];
    /// <summary>A file the Reader shows instead of opening it as a structure.</summary>
    public static bool IsReaderFile(string path)
    {
        var name = Path.GetFileName(path).ToLowerInvariant();
        return name.StartsWith("log.") || name.StartsWith("in.") || ReaderExts.Contains(Path.GetExtension(name));
    }

    private string _rdPath = "", _rdText = "", _rdInfo = "", _rdFind = "", _rdFindText = "";
    private bool _rdWrap, _rdFollow;
    private long _rdSize;
    private DateTime _rdStamp;
    private Avalonia.Threading.DispatcherTimer? _rdTimer;
    private const int ReaderMaxChars = 4_000_000;   // a longer file shows its last 4 MB (the tables come from the whole file)
    public string ReaderPath => _rdPath;
    public string ReaderName => Path.GetFileName(_rdPath);
    public string ReaderText { get => _rdText; private set => Set(ref _rdText, value); }
    public string ReaderInfo { get => _rdInfo; private set => Set(ref _rdInfo, value); }
    public bool ReaderWrap { get => _rdWrap; set => Set(ref _rdWrap, value); }
    public bool ReaderFollow
    {
        get => _rdFollow;
        set
        {
            if (!Set(ref _rdFollow, value)) return;
            if (value)
            {
                _rdTimer ??= new Avalonia.Threading.DispatcherTimer { Interval = TimeSpan.FromSeconds(2) };
                _rdTimer.Tick -= ReaderTick;
                _rdTimer.Tick += ReaderTick;
                _rdTimer.Start();
            }
            else _rdTimer?.Stop();
        }
    }
    public ObservableCollection<TextTable> ReaderTables { get; } = new();
    public bool ReaderHasTables => ReaderTables.Count > 0;
    private int _rdTable, _rdX, _rdY = 1;
    private double _rdFrom = 0.5;
    public ObservableCollection<string> ReaderColumns { get; } = new();
    public int ReaderTable { get => _rdTable; set { if (Set(ref _rdTable, Math.Clamp(value, 0, Math.Max(0, ReaderTables.Count - 1)))) ReaderPickTable(); } }
    public int ReaderX { get => _rdX; set { if (Set(ref _rdX, Math.Max(0, value))) ReaderPlot(); } }
    public int ReaderY { get => _rdY; set { if (Set(ref _rdY, Math.Max(0, value))) ReaderPlot(); } }
    public decimal ReaderFromD { get => (decimal)(_rdFrom * 100); set { if (Set(ref _rdFrom, Math.Clamp((double)value / 100, 0, 0.99), nameof(ReaderFromD))) ReaderPlot(); } }
    public (double X, double Y)[] ReaderCurve { get; private set; } = [];
    // a fit over the plotted points (the Analyze fitter: line, power law, exponential, KWW, Arrhenius, Gaussian peaks)
    private int _rdFit;
    public int ReaderFit { get => _rdFit; set { if (Set(ref _rdFit, Math.Clamp(value, 0, CurveFit.Models.Length - 1))) ReaderPlot(); } }
    public (double X, double Y)[] ReaderFitLine { get; private set; } = [];
    private string _rdFitText = "";
    public string ReaderFitText { get => _rdFitText; private set => Set(ref _rdFitText, value); }
    public string ReaderXLabel => _rdX < ReaderColumns.Count ? ReaderColumns[_rdX] : "";
    public string ReaderYLabel => _rdY < ReaderColumns.Count ? ReaderColumns[_rdY] : "";
    private string _rdStats = "";
    public string ReaderStats { get => _rdStats; private set => Set(ref _rdStats, value); }
    public event Action? ReaderChanged;
    /// <summary>The view scrolls to the end (follow) or selects a find match: (start, length); (-1, 0) the end.</summary>
    public event Action<int, int>? ReaderSelect;

    public string ReaderFind { get => _rdFind; set { if (Set(ref _rdFind, value ?? "")) { _rdMatch = -1; ReaderFindNext(); } } }
    public string ReaderFindText { get => _rdFindText; private set => Set(ref _rdFindText, value); }
    private int _rdMatch = -1;

    public void OpenReader(string path)
    {
        if (_module != 73) _returnModule = _module;
        _rdPath = path;
        _rdTable = 0;
        ReaderLoad(true);
        SetModule(73);
        Raise(nameof(ReaderPath)); Raise(nameof(ReaderName));
        Remember(path, null);
        Status = $"{ReaderName}: " + (ReaderTables.Count > 0 ? $"{ReaderTables.Count} table(s) of numbers" : "text");
    }

    public void ReaderReload() { if (_rdPath.Length > 0) ReaderLoad(false); }

    private void ReaderTick(object? s, EventArgs e)
    {
        if (!_rdFollow || _module != 73 || _rdPath.Length == 0) return;
        try
        {
            var fi = new FileInfo(_rdPath);
            if (fi.Exists && (fi.Length != _rdSize || fi.LastWriteTimeUtc != _rdStamp)) ReaderLoad(false);
        }
        catch { }
    }

    private void ReaderLoad(bool fresh)
    {
        string text;
        try
        {
            var fi = new FileInfo(_rdPath);
            _rdSize = fi.Length;
            _rdStamp = fi.LastWriteTimeUtc;
            using var fs = new FileStream(_rdPath, FileMode.Open, FileAccess.Read, FileShare.ReadWrite | FileShare.Delete);   // a log still being written
            using var sr = new StreamReader(fs);
            text = sr.ReadToEnd();
        }
        catch (Exception e) { ReaderText = ""; ReaderInfo = "Could not read the file: " + e.Message; return; }
        var inv = CultureInfo.InvariantCulture;
        var lines = 1 + text.Count(c => c == '\n');
        var shown = text.Length > ReaderMaxChars ? text[^ReaderMaxChars..] : text;
        ReaderText = shown;
        var kind = TextTables.Kind(_rdPath, text) switch
        {
            "lammps-log" => "LAMMPS log", "xvg" => "GROMACS xvg", "csv" => "CSV", "tsv" => "TSV",
            _ => Path.GetExtension(_rdPath).ToLowerInvariant() switch { ".mdp" => "GROMACS parameters", ".in" => "LAMMPS input", ".md" => "notes", _ => "text" },
        };
        ReaderInfo = $"{kind} · {lines.ToString("N0", inv)} lines · {(_rdSize / 1024.0).ToString("0.#", inv)} kB · changed {File.GetLastWriteTime(_rdPath):yyyy-MM-dd HH:mm:ss}"
                     + (shown.Length < text.Length ? " · showing the last 4 MB" : "");
        var keepT = fresh ? 0 : _rdTable;
        var (keepX, keepY) = (_rdX, _rdY);
        ReaderTables.Clear();
        foreach (var t in TextTables.Parse(_rdPath, text)) ReaderTables.Add(t);
        Raise(nameof(ReaderHasTables));
        _rdTable = Math.Clamp(keepT, 0, Math.Max(0, ReaderTables.Count - 1));
        Raise(nameof(ReaderTable));
        ReaderPickTable(fresh ? null : (keepX, keepY));
        if (_rdFollow) ReaderSelect?.Invoke(-1, 0);
    }

    private void ReaderPickTable((int X, int Y)? keep = null)
    {
        ReaderColumns.Clear();
        if (_rdTable < ReaderTables.Count) foreach (var c in ReaderTables[_rdTable].Columns) ReaderColumns.Add(c);
        var n = ReaderColumns.Count;
        if (keep is { } k && k.X < n && k.Y < n) (_rdX, _rdY) = k;
        else
        {
            _rdX = 0;
            // the first column that is not the counter: for a LAMMPS log the temperature or energy after Step
            _rdY = n > 1 ? 1 : 0;
        }
        Raise(nameof(ReaderX)); Raise(nameof(ReaderY));
        ReaderPlot();
    }

    private void ReaderPlot()
    {
        var inv = CultureInfo.InvariantCulture;
        if (_rdTable >= ReaderTables.Count || _rdX >= ReaderColumns.Count || _rdY >= ReaderColumns.Count)
        {
            ReaderCurve = [];
            ReaderStats = "";
        }
        else
        {
            var t = ReaderTables[_rdTable];
            var x = t.Column(_rdX);
            var y = t.Column(_rdY);
            ReaderCurve = x.Zip(y).Where(p => double.IsFinite(p.First) && double.IsFinite(p.Second)).ToArray();
            ReaderFitLine = [];
            ReaderFitText = "";
            if (_rdFit > 0 && ReaderCurve.Length >= 3)
            {
                var fx = ReaderCurve.Select(p => p.X).ToArray();
                var fit = CurveFit.Fit(_rdFit, fx, ReaderCurve.Select(p => p.Y).ToArray(), fx.Min(), fx.Max(), 1);
                if (fit != null) { ReaderFitLine = fit.Line; ReaderFitText = fit.Text; }
                else ReaderFitText = "the fit did not converge on these points";
            }
            var st = TextTables.Stats(y, _rdFrom);
            string F(double v) => double.IsFinite(v) ? (Math.Abs(v) >= 1e5 || (Math.Abs(v) < 1e-3 && v != 0) ? v.ToString("0.###E+0", inv) : v.ToString("0.####", inv)) : "—";
            ReaderStats = st.N == 0 ? "no numbers" :
                $"{ReaderYLabel} over the last {(100 - _rdFrom * 100).ToString("0", inv)} % ({st.N} rows): mean {F(st.Mean)}" + (double.IsFinite(st.Err) ? $" ± {F(st.Err)}" : "") +
                $" · min {F(st.Min)} · max {F(st.Max)} · last {F(st.Last)}";
        }
        Raise(nameof(ReaderXLabel)); Raise(nameof(ReaderYLabel));
        ReaderChanged?.Invoke();
    }

    public void ReaderFindNext() => ReaderFindStep(+1);
    public void ReaderFindPrev() => ReaderFindStep(-1);
    private void ReaderFindStep(int dir)
    {
        if (_rdFind.Length == 0 || _rdText.Length == 0) { ReaderFindText = ""; return; }
        var all = new List<int>();
        for (var i = _rdText.IndexOf(_rdFind, StringComparison.OrdinalIgnoreCase); i >= 0 && all.Count < 100000; i = _rdText.IndexOf(_rdFind, i + 1, StringComparison.OrdinalIgnoreCase)) all.Add(i);
        if (all.Count == 0) { ReaderFindText = "not found"; return; }
        _rdMatch = _rdMatch < 0 ? (dir > 0 ? 0 : all.Count - 1) : ((_rdMatch + dir) % all.Count + all.Count) % all.Count;
        ReaderFindText = $"{_rdMatch + 1} of {all.Count}";
        ReaderSelect?.Invoke(all[_rdMatch], _rdFind.Length);
    }

    /// <summary>The chosen table as CSV text (for the clipboard).</summary>
    public string ReaderTableCsv()
    {
        if (_rdTable >= ReaderTables.Count) return "";
        var t = ReaderTables[_rdTable];
        var inv = CultureInfo.InvariantCulture;
        var sb = new System.Text.StringBuilder();
        sb.AppendLine(string.Join(",", t.Columns.Select(c => c.Contains(',') ? $"\"{c}\"" : c)));
        foreach (var r in t.Rows) sb.AppendLine(string.Join(",", r.Select(v => v.ToString("R", inv))));
        return sb.ToString();
    }

    /// <summary>Jumps the text to the chosen table's first line.</summary>
    public void ReaderGoToTable()
    {
        if (_rdTable >= ReaderTables.Count) return;
        var line = ReaderTables[_rdTable].Line;
        var pos = 0;
        for (var k = 1; k < line && pos >= 0; ++k) { pos = _rdText.IndexOf('\n', pos); if (pos >= 0) ++pos; }
        if (pos >= 0) ReaderSelect?.Invoke(pos, 0);
    }
}
