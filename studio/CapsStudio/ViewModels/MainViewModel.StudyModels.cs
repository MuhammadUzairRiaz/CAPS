using System.Globalization;
using System.Text;
using System.Text.Json.Nodes;
using CapsStudio.Interop;

namespace CapsStudio.ViewModels;

/// <summary>Study › statistics by the core (qspr.hpp): stepwise, PLS, GFA and neural-network models of the Up column with
/// cross-validation and y-randomisation, clusters of the rows, discriminant analysis, and descriptors of the rows' structures.</summary>
public sealed partial class MainViewModel
{
    private int _studyK = 3, _studyRun;
    public bool StudyIsCluster => _studyMode == 8;
    public decimal? StudyK { get => _studyK; set { if (value != null && Set(ref _studyK, Math.Clamp((int)value.Value, 2, 20)) && _studyMode == 8) StudyStatsRun(); } }

    // the rows with every chosen column a number
    private static int[] CompleteRows(IReadOnlyList<double[]> cols)
    {
        var n = cols.Count == 0 ? 0 : cols[0].Length;
        return Enumerable.Range(0, n).Where(r => cols.All(c => double.IsFinite(c[r]))).ToArray();
    }

    private void StudyModelRun(List<StudyColumn> cols, List<StudyColumn> ticked)
    {
        var inv = CultureInfo.InvariantCulture;
        string N(JsonNode? v) => v is JsonNode x ? StudyStats.N((double?)x ?? double.NaN) : "–";
        var mode = _studyMode;
        JsonObject json;
        string head;
        int[] used;
        if (mode == 8)
        {
            if (ticked.Count < 1) { StudyStatsText = "Tick the columns to cluster the rows on"; return; }
            var vals = ticked.Select(Study.Values).ToList();
            used = CompleteRows(vals);
            if (used.Length < _studyK) { StudyStatsText = "Fewer complete rows than clusters"; return; }
            json = new JsonObject { ["op"] = "cluster", ["method"] = "ward", ["k"] = _studyK, ["rows"] = new JsonArray(used.Select(r => (JsonNode?)new JsonArray(vals.Select(c => (JsonNode?)c[r]).ToArray())).ToArray()) };
            head = $"Ward clusters of {used.Length} rows on {string.Join(", ", ticked.Select(c => c.Name))} (autoscaled)";
        }
        else
        {
            if (_studyY >= cols.Count) { StudyStatsText = "Choose y (the chart's Up column)"; return; }
            var yc = cols[_studyY];
            var xs = ticked.Where(c => c != yc).ToList();
            if (xs.Count == 0) { StudyStatsText = $"Tick the columns that predict {yc.Name}"; return; }
            var all = xs.Select(Study.Values).Prepend(Study.Values(yc)).ToList();
            used = mode == 9 ? CompleteRows(xs.Select(Study.Values).ToList()).Where(r => Study.Rows[r].Cells[cols.IndexOf(yc)].Text.Trim().Length > 0).ToArray() : CompleteRows(all);
            if (used.Length < 3) { StudyStatsText = "Three complete rows or more"; return; }
            if (mode == 9)
            {
                var yi = cols.IndexOf(yc);
                json = new JsonObject
                {
                    ["op"] = "lda", ["classes"] = new JsonArray(used.Select(r => (JsonNode?)Study.Rows[r].Cells[yi].Text.Trim()).ToArray()),
                    ["rows"] = new JsonArray(used.Select(r => (JsonNode?)new JsonArray(xs.Select(c => (JsonNode?)Study.Values(c)[r]).ToArray())).ToArray()),
                };
                head = $"Discriminant analysis: {yc.Name} as classes, {used.Length} rows, on {string.Join(", ", xs.Select(c => c.Name))}";
            }
            else
            {
                var method = mode switch { 5 => "pls", 6 => "gfa", 7 => "nn", _ => "mlr" };
                json = new JsonObject
                {
                    ["op"] = "model", ["method"] = method, ["stepwise"] = mode == 4, ["squares"] = mode == 4, ["products"] = mode == 4,
                    ["folds"] = used.Length <= 30 ? used.Length : 10, ["randomise"] = 10,
                    ["y"] = new JsonArray(used.Select(r => (JsonNode?)all[0][r]).ToArray()),
                    ["x"] = new JsonArray(xs.Select(c => (JsonNode?)new JsonArray(used.Select(r => (JsonNode?)Study.Values(c)[r]).ToArray())).ToArray()),
                    ["names"] = new JsonArray(xs.Select(c => (JsonNode?)c.Name).ToArray()),
                };
                head = $"{StudyModes[mode]} of {yc.Name} over {used.Length} rows";
            }
        }
        StudyStatsText = head + "\ncomputing…";
        var ticket = ++_studyRun;
        var text = json.ToJsonString();
        Task.Run(() => CapsDocument.Qspr(text)).ContinueWith(t =>
        {
            Avalonia.Threading.Dispatcher.UIThread.Post(() =>
            {
                if (ticket != _studyRun) return;
                try
                {
                    var r = JsonNode.Parse(t.Result)!;
                    if (r["ok"]?.GetValue<bool>() != true) { StudyStatsText = head + "\n" + ((string?)r["error"] ?? "failed"); return; }
                    var sb = new StringBuilder(head + "\n");
                    if (mode == 8)
                    {
                        var lab = ((JsonArray)r["labels"]!).Select(v => (int?)v ?? 0).ToArray();
                        for (var c = 0; c < _studyK; ++c)
                            sb.AppendLine($"cluster {c + 1}: " + string.Join(", ", used.Where((_, i) => lab[i] == c).Select(i => Study.Rows[i].Label.Length > 0 ? Study.Rows[i].Label : $"row {i + 1}")));
                        sb.Append("within sum of squares " + N(r["within_ss"]));
                    }
                    else if (mode == 9)
                    {
                        sb.AppendLine($"right on the rows {100 * ((double?)r["accuracy"] ?? 0):0.#} % · left out one at a time {100 * ((double?)r["accuracy_loo"] ?? 0):0.#} %");
                        var dirs = (JsonArray)r["directions"]!;
                        for (var f = 0; f < dirs.Count; ++f) sb.AppendLine($"function {f + 1}: " + string.Join("  ", ((JsonArray)dirs[f]!).Select(v => N(v))));
                    }
                    else
                    {
                        sb.AppendLine((string?)r["equation"] ?? "");
                        sb.AppendLine($"R² {N(r["r2"])} · adjusted {N(r["r2_adj"])} · RMSE {N(r["rmse"])} · q² {N(r["q2"])} · RMSE(CV) {N(r["rmse_cv"])}");
                        var rr = ((JsonArray)r["random_r2"]!).Select(v => (double?)v ?? 0).ToArray();
                        if (rr.Length > 0) sb.AppendLine($"y-randomisation: shuffled R² mean {StudyStats.N(rr.Average())}, largest {StudyStats.N(rr.Max())}");
                        var cook = ((JsonArray)r["cooks"]!).Select(v => (double?)v ?? 0).ToArray();
                        if (cook.Length == used.Length && cook.Length > 0)
                        {
                            var top = Enumerable.Range(0, cook.Length).OrderByDescending(i => cook[i]).Take(3);
                            sb.AppendLine("most influential (Cook's distance): " + string.Join(", ", top.Select(i => $"{(Study.Rows[used[i]].Label.Length > 0 ? Study.Rows[used[i]].Label : $"row {used[i] + 1}")} {StudyStats.N(cook[i])}")));
                        }
                        foreach (var n in (JsonArray)r["notes"]!) sb.AppendLine((string?)n);
                    }
                    StudyStatsText = sb.ToString().TrimEnd();
                }
                catch (Exception e) { StudyStatsText = head + "\n" + e.Message; }
            });
        });
    }

    /// <summary>Descriptors of the rows' structures (rows that carry a file) as new columns.</summary>
    public async Task StudyAddDescriptors()
    {
        var rows = Study.Rows.Where(r => r.HasPath).ToList();
        if (rows.Count == 0) { StudyError = "Descriptors: no row carries a structure file"; return; }
        string[] keep = ["W", "1chi", "1chiv", "kappa1", "kappa2", "J", "rings", "hbond_donors", "hbond_acceptors", "rotatable_bonds", "vdw_volume", "vdw_area", "SASA", "DPSA1", "FPSA3", "Rg"];
        var json = new JsonObject { ["op"] = "descriptors", ["files"] = new JsonArray(rows.Select(r => (JsonNode?)r.Path).ToArray()) }.ToJsonString();
        Status = "Descriptors…";
        try
        {
            var r = JsonNode.Parse(await Task.Run(() => CapsDocument.Qspr(json)))!;
            if (r["ok"]?.GetValue<bool>() != true) { StudyError = "Descriptors: " + ((string?)r["error"] ?? "failed"); return; }
            var list = (JsonArray)r["rows"]!;
            var cols = keep.ToDictionary(k => k, k => Study.Columns.FirstOrDefault(c => c.Name == k) ?? Study.AddColumn(k));
            var failed = 0;
            for (var i = 0; i < rows.Count; ++i)
            {
                if (list[i]?["descriptors"] is not JsonArray ds) { ++failed; continue; }
                foreach (var d in ds)
                {
                    var name = (string?)d!["name"] ?? "";
                    if (!cols.TryGetValue(name, out var col)) continue;
                    var v = d["value"] is JsonNode x ? (double?)x ?? double.NaN : double.NaN;
                    rows[i].Cells[Study.Columns.IndexOf(col)].Set(double.IsFinite(v) ? v.ToString("R", CultureInfo.InvariantCulture) : "");
                }
            }
            Study.Recompute();
            Status = $"Descriptors of {rows.Count - failed} structures added as {keep.Length} columns" + (failed > 0 ? $" ({failed} could not be read)" : "");
        }
        catch (Exception e) { StudyError = "Descriptors: " + e.Message; }
    }
}
