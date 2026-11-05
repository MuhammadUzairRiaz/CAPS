using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using System.Text.Json.Nodes;
using Avalonia.Media;
using CapsStudio.Interop;

namespace CapsStudio.ViewModels;

// Study table extras (core table_tools.hpp): rows ordered by several columns, a filter that shows the rows meeting a
// condition, cells coloured by value, a histogram of a column in the chart, and an SD file of the rows shown.
public partial class MainViewModel
{
    // ---------------------------------------------------------------- order
    private string _studySort = "";
    /// <summary>Columns to order by, first key first, "-" for descending: "Tg, -density".</summary>
    public string StudySort { get => _studySort; set => Set(ref _studySort, value ?? ""); }
    public string StudyApplySort()
    {
        var keys = new JsonArray();
        var desc = new JsonArray();
        foreach (var part in _studySort.Split(',', StringSplitOptions.RemoveEmptyEntries | StringSplitOptions.TrimEntries))
        {
            var d = part.StartsWith('-') || part.StartsWith('−');
            var name = (d ? part[1..] : part).Trim().Trim('[', ']');
            var col = Study.Columns.FirstOrDefault(c => c.Name == name);
            if (col == null) return StudyError = $"no column '{name}'";
            keys.Add(new JsonArray(Study.Values(col).Select(v => (JsonNode?)(double.IsFinite(v) ? JsonValue.Create(v) : null)).ToArray()));
            desc.Add(d);
        }
        if (keys.Count == 0) return StudyError = "name the columns to order by, e.g. Tg, -density";
        var r = JsonNode.Parse(CapsDocument.SortRows(new JsonObject { ["keys"] = keys, ["descending"] = desc }.ToJsonString()))!;
        if (r["ok"]?.GetValue<bool>() != true) return StudyError = (string?)r["error"] ?? "the order failed";
        var rows = Study.Rows.ToList();
        var order = ((JsonArray)r["order"]!).Select(x => (int)x!.GetValue<double>()).ToList();
        for (var k = 0; k < order.Count; ++k)
        {
            var cur = Study.Rows.IndexOf(rows[order[k]]);
            if (cur != k) Study.Rows.Move(cur, k);
        }
        StudyError = "";
        return $"Rows ordered by {_studySort}";
    }

    // ---------------------------------------------------------------- filter
    private string _studyFilter = "";
    /// <summary>A condition on the columns: Tg > 300 and density &lt; 1.2. Empty shows every row.</summary>
    public string StudyFilter { get => _studyFilter; set => Set(ref _studyFilter, value ?? ""); }
    private string _studyFilterText = "";
    public string StudyFilterText { get => _studyFilterText; private set => Set(ref _studyFilterText, value); }
    public string StudyApplyFilter()
    {
        if (_studyFilter.Trim().Length == 0)
        {
            foreach (var r in Study.Rows) r.Visible = true;
            StudyFilterText = "";
            return "Every row shown";
        }
        Formula f;
        try { f = Formula.Parse(_studyFilter); }
        catch (Exception e) { return StudyError = "filter: " + e.Message; }
        var index = Study.Columns.Select((c, i) => (c.Name, i)).ToDictionary(p => p.Name, p => p.i);
        var shown = 0;
        foreach (var r in Study.Rows)
        {
            double V(string name) => index.TryGetValue(name, out var i) && i < r.Cells.Count ? r.Cells[i].Number : double.NaN;
            double v;
            try { v = f.Eval(V); } catch { v = double.NaN; }
            r.Visible = double.IsFinite(v) && v != 0;
            if (r.Visible) ++shown;
        }
        StudyError = "";
        StudyFilterText = $"{shown} of {Study.Rows.Count} rows meet {_studyFilter}";
        return StudyFilterText;
    }

    // ---------------------------------------------------------------- colour by value
    public string[] StudyColourNames => new[] { "No colours" }.Concat(Study.Columns.Select(c => c.Name)).ToArray();
    private int _studyColour;
    public int StudyColour { get => _studyColour; set { if (Set(ref _studyColour, Math.Max(0, value))) StudyApplyColour(); } }
    // viridis at five stops (as the plots), translucent behind the numbers
    private static readonly Color[] Ramp = [Color.FromRgb(68, 1, 84), Color.FromRgb(59, 82, 139), Color.FromRgb(33, 145, 140), Color.FromRgb(94, 201, 98), Color.FromRgb(253, 231, 37)];
    private static IBrush RampBrush(double t)
    {
        t = Math.Clamp(t, 0, 1) * (Ramp.Length - 1);
        var k = Math.Min((int)t, Ramp.Length - 2);
        var u = t - k;
        byte L(byte a, byte b) => (byte)Math.Round(a + (b - a) * u);
        return new SolidColorBrush(Color.FromArgb(0x88, L(Ramp[k].R, Ramp[k + 1].R), L(Ramp[k].G, Ramp[k + 1].G), L(Ramp[k].B, Ramp[k + 1].B)));
    }
    public void StudyApplyColour()
    {
        foreach (var r in Study.Rows) foreach (var c in r.Cells) c.Fill = null;
        if (_studyColour <= 0 || _studyColour > Study.Columns.Count) return;
        var col = Study.Columns[_studyColour - 1];
        var v = Study.Values(col).Where(double.IsFinite).ToArray();
        if (v.Length == 0) return;
        double lo = v.Min(), hi = v.Max();
        var i = Study.Columns.IndexOf(col);
        foreach (var r in Study.Rows)
            if (i < r.Cells.Count && double.IsFinite(r.Cells[i].Number))
                r.Cells[i].Fill = RampBrush(hi > lo ? (r.Cells[i].Number - lo) / (hi - lo) : 0.5);
    }

    // ---------------------------------------------------------------- histogram of the x column
    public static readonly string[] StudyBinRules = ["Freedman–Diaconis", "Scott", "Sturges"];
    private static readonly string[] BinRuleIds = ["fd", "scott", "sturges"];
    private bool _studyHist;
    public bool StudyHist { get => _studyHist; set { if (Set(ref _studyHist, value)) StudyPlot(); } }
    private int _studyBinRule;
    public int StudyBinRule { get => _studyBinRule; set { if (Set(ref _studyBinRule, Math.Clamp(value, 0, 2))) StudyPlot(); } }
    public (double X, double Y)[] StudyHistBars { get; private set; } = [];
    private void StudyHistCompute()
    {
        StudyHistBars = [];
        if (!_studyHist || _studyX >= Study.Columns.Count) return;
        var vals = new JsonArray(Study.Values(Study.Columns[_studyX]).Select(v => (JsonNode?)(double.IsFinite(v) ? JsonValue.Create(v) : null)).ToArray());
        try
        {
            var r = JsonNode.Parse(CapsDocument.Histogram(new JsonObject { ["values"] = vals, ["rule"] = BinRuleIds[_studyBinRule] }.ToJsonString()))!;
            if (r["ok"]?.GetValue<bool>() != true) { StudyFitText = (string?)r["error"] ?? ""; return; }
            var e = ((JsonArray)r["edges"]!).Select(x => x!.GetValue<double>()).ToArray();
            var c = ((JsonArray)r["counts"]!).Select(x => x!.GetValue<double>()).ToArray();
            StudyHistBars = c.Select((n, k) => ((e[k] + e[k + 1]) / 2, n)).ToArray();
            var inv = CultureInfo.InvariantCulture;
            StudyFitText = $"{(int)r["n"]!.GetValue<double>()} values in {c.Length} bins of {r["width"]!.GetValue<double>().ToString("0.####", inv)} ({(string?)r["rule"]})"
                           + ((int)r["left_out"]!.GetValue<double>() > 0 ? $", {(int)r["left_out"]!.GetValue<double>()} without a number left out" : "");
        }
        catch (Exception ex) { StudyFitText = ex.Message; }
    }

    // ---------------------------------------------------------------- SD file of the rows shown
    public string StudyExportSdf(string path)
    {
        var rows = new JsonArray();
        foreach (var r in Study.Rows.Where(r => r.Visible))
        {
            var data = new JsonArray();
            for (var i = 0; i < Study.Columns.Count && i < r.Cells.Count; ++i)
                if (r.Cells[i].Text.Length > 0) data.Add(new JsonArray(Study.Columns[i].Name.Replace('<', '(').Replace('>', ')'), r.Cells[i].Text));
            rows.Add(new JsonObject { ["title"] = r.Label, ["structure"] = r.Path, ["data"] = data });
        }
        var j = JsonNode.Parse(CapsDocument.SdfTable(new JsonObject { ["path"] = path, ["rows"] = rows }.ToJsonString()))!;
        if (j["ok"]?.GetValue<bool>() != true) return StudyError = (string?)j["error"] ?? "the SD file failed";
        var notes = ((JsonArray)j["notes"]!).Select(x => (string?)x ?? "");
        return $"Wrote {(int)j["records"]!.GetValue<double>()} records to {System.IO.Path.GetFileName(path)}" + (notes.Any() ? " · " + string.Join(" · ", notes) : "");
    }
}
