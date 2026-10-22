using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Text.Json.Nodes;
using System.Text.RegularExpressions;
using CapsStudio.Interop;

namespace CapsStudio.ViewModels;

/// <summary>Runs that finished elsewhere (a cluster job of caps equilibrate or caps tg), opened on the page that would have
/// run them: the same curves, stages, block means, criteria and results as a run started here.</summary>
public sealed partial class MainViewModel
{
    private string _eqOpenedText = "", _gtOpenedText = "";
    /// <summary>What the Equilibrate page shows when it is a finished run read back, not one run here (empty otherwise).</summary>
    public string EqOpenedText { get => _eqOpenedText; private set { if (Set(ref _eqOpenedText, value)) Raise(nameof(EqShowsFinished)); } }
    public bool EqShowsFinished => _eqOpenedText.Length > 0;
    private string _eqOpenedFiles = "";
    /// <summary>The files a finished equilibration was read from (thermo log + output + structure).</summary>
    public string EqOpenedFiles { get => _eqOpenedFiles; private set => Set(ref _eqOpenedFiles, value); }
    public string GtOpenedText { get => _gtOpenedText; private set { if (Set(ref _gtOpenedText, value)) Raise(nameof(GtShowsFinished)); } }
    public bool GtShowsFinished => _gtOpenedText.Length > 0;
    private IReadOnlyList<string> _gtOpenedPaths = [];

    private static readonly Regex ProtocolLine = new(@"^(nvt|npt|nph|nve)\b.*#\s*(\d+)\s*·", RegexOptions.IgnoreCase | RegexOptions.Compiled);
    private static readonly Regex CheckLine = new(@"^check (.+?)\s{2,}(ok|open)\s+change ([-0-9.eE+]+) \(tolerance ([-0-9.eE+]+)\)\s+blocks:(.*)$", RegexOptions.Compiled);
    private static readonly Regex ProgressLine = new(@"^\[\s*\d+/\d+\]\s.*?([-0-9.]+) ps\s+T", RegexOptions.Compiled);
    private static readonly Regex BlocksLine = new(@"after (\d+) production blocks of ([-0-9.eE+]+) ps", RegexOptions.Compiled);

    /// <summary>The files of a finished run: those picked, and when a folder is picked (or one file) the others beside them
    /// (the folder and its out/).</summary>
    private static List<string> RunFiles(IReadOnlyList<string> picked)
    {
        var files = new List<string>();
        void Folder(string d)
        {
            if (!Directory.Exists(d)) return;
            files.AddRange(Directory.GetFiles(d).OrderBy(f => f, StringComparer.Ordinal));
            if (Directory.Exists(Path.Combine(d, "out"))) files.AddRange(Directory.GetFiles(Path.Combine(d, "out")).OrderBy(f => f, StringComparer.Ordinal));
        }
        foreach (var p in picked)
            if (Directory.Exists(p)) Folder(p);
            else if (File.Exists(p)) files.Add(p);
        if (picked.Count == 1 && File.Exists(picked[0])) Folder(Path.GetDirectoryName(Path.GetFullPath(picked[0])) ?? ".");
        return files.Distinct().ToList();
    }

    private static string Head(string path, int bytes = 4096)
    {
        try
        {
            using var f = File.OpenRead(path);
            var b = new byte[Math.Min(bytes, (int)Math.Min(int.MaxValue, f.Length))];
            var n = f.Read(b, 0, b.Length);
            return System.Text.Encoding.UTF8.GetString(b, 0, n);
        }
        catch (Exception) { return ""; }
    }

    /// <summary>Opens an equilibration that ran elsewhere: its thermo log (caps equilibrate --log, step,time_ps,stage,…), the
    /// run's printed output when it lies beside it (the protocol, the convergence checks and the report) and the structure
    /// it wrote. The page shows them as it shows a run started here.</summary>
    public string OpenFinishedEquilibrate(IReadOnlyList<string> picked)
    {
        if (!Idle) return Status = "Wait for the run to finish (or stop it) before opening a finished one";
        var inv = CultureInfo.InvariantCulture;
        var files = RunFiles(picked);
        // the thermo log: the CSV whose rows carry the protocol's stages (an MD log has "md" throughout)
        string? csv = null;
        foreach (var f in files.Where(f => f.EndsWith(".csv", StringComparison.OrdinalIgnoreCase)))
        {
            var head = Head(f);
            if (!head.StartsWith("step,time_ps,stage", StringComparison.Ordinal)) continue;
            var second = head.Split('\n').Skip(1).FirstOrDefault() ?? "";
            var stage = second.Split(',').ElementAtOrDefault(2) ?? "";
            if (stage == "md" && csv != null) continue;
            if (csv == null || stage != "md") csv = f;
            if (stage != "md" && picked.Contains(f)) break;
        }
        if (csv == null) return Status = "No thermo log of an equilibration here: caps equilibrate --log FILE.csv writes it (step, time_ps, stage, …)";
        var rows = new List<CapsThermo>();
        var labels = new List<string>();
        var cols = File.ReadLines(csv).First().Split(',');
        int Col(string n) => Array.IndexOf(cols, n);
        int cStep = Col("step"), cT = Col("time_ps"), cStage = Col("stage"), cTemp = Col("temperature_K"), cPot = Col("potential"), cKin = Col("kinetic"),
            cTot = Col("total"), cCons = Col("conserved"), cP = Col("pressure_atm"), cV = Col("volume_A3"), cRho = Col("density_g_cm3");
        double D(string[] f, int c) => c >= 0 && c < f.Length && double.TryParse(f[c], NumberStyles.Float, inv, out var v) ? v : 0;
        foreach (var line in File.ReadLines(csv).Skip(1))
        {
            var f = line.Split(',');
            if (f.Length < cols.Length) continue;
            rows.Add(new CapsThermo
            {
                Step = (long)D(f, cStep), TimePs = D(f, cT), Temperature = D(f, cTemp), Potential = D(f, cPot), Kinetic = D(f, cKin), Total = D(f, cTot),
                Conserved = D(f, cCons), Pressure = D(f, cP), Volume = D(f, cV), Density = D(f, cRho),
            });
            labels.Add(cStage >= 0 ? f[cStage] : "");
        }
        if (rows.Count < 2) return Status = Path.GetFileName(csv) + ": fewer than two rows";

        // the run's printed output: the one whose progress reaches the log's last time
        string? log = null;
        string[] logLines = [];
        foreach (var f in files.Where(f => f.EndsWith(".log", StringComparison.OrdinalIgnoreCase) || f.EndsWith(".out", StringComparison.OrdinalIgnoreCase) || f.EndsWith(".txt", StringComparison.OrdinalIgnoreCase)))
        {
            if (new FileInfo(f).Length > 20 << 20) continue;
            var lines = File.ReadAllLines(f);
            if (!lines.Any(l => l.StartsWith("# total", StringComparison.Ordinal)) || !lines.Any(l => l.StartsWith("[", StringComparison.Ordinal))) continue;
            var last = lines.Select(l => ProgressLine.Match(l)).Where(m => m.Success).Select(m => double.Parse(m.Groups[1].Value, inv)).DefaultIfEmpty(double.NaN).Last();
            if (log == null || Math.Abs(last - rows[^1].TimePs) < 5) { log = f; logLines = lines; }
            if (Math.Abs(last - rows[^1].TimePs) < 5) break;
        }

        // the structure the run wrote (its "wrote NAME" line), or one picked
        string? structure = picked.FirstOrDefault(p => File.Exists(p) && !p.EndsWith(".csv", StringComparison.OrdinalIgnoreCase) && !p.EndsWith(".log", StringComparison.OrdinalIgnoreCase)
                                                        && !p.EndsWith(".out", StringComparison.OrdinalIgnoreCase) && !p.EndsWith(".txt", StringComparison.OrdinalIgnoreCase));
        if (structure == null && log != null && logLines.LastOrDefault(l => l.StartsWith("wrote ", StringComparison.Ordinal)) is { } wrote)
        {
            var name = wrote[6..].Trim().Split(' ')[0];
            foreach (var dir in new[] { Path.GetDirectoryName(log), Path.GetDirectoryName(csv) }.Distinct())
                if (dir != null && File.Exists(Path.Combine(dir, name))) { structure = Path.Combine(dir, name); break; }
        }
        if (structure != null)
        {
            Open(structure);
            if (_doc?.Path != structure) return Status;   // could not open it: Open said why
        }

        // stages: where each begins, and its ensemble (the protocol's lines, else whether the volume moved in it)
        var ensemble = new Dictionary<int, string>();
        foreach (var l in logLines)
            if (ProtocolLine.Match(l) is { Success: true } m) ensemble[int.Parse(m.Groups[2].Value, inv)] = m.Groups[1].Value.ToUpperInvariant();
        var starts = new List<(int Row, string Label)>();
        for (var i = 0; i < labels.Count; ++i) if (i == 0 || labels[i] != labels[i - 1]) starts.Add((i, labels[i]));
        var bands = starts.Select((s, k) =>
        {
            var end = k + 1 < starts.Count ? starts[k + 1].Row : rows.Count;
            var from = k == 0 ? 0 : rows[s.Row].TimePs;
            var to = k + 1 < starts.Count ? rows[starts[k + 1].Row].TimePs : rows[^1].TimePs;
            string e;
            if (s.Label.StartsWith("production", StringComparison.Ordinal)) e = "NPT";
            else if (int.TryParse(s.Label.Split(' ')[0], NumberStyles.Integer, inv, out var n) && ensemble.TryGetValue(n, out var en)) e = en;
            else
            {
                var vols = rows.Skip(s.Row).Take(end - s.Row).Select(r => r.Volume).ToArray();
                e = vols.Length > 1 && vols.Max() - vols.Min() > 1e-6 * Math.Max(1, vols.Max()) ? "NPT" : "NVT";
            }
            return (from, to, e, e is "NPT" or "NPH");
        }).ToArray();
        // neighbouring stages of one ensemble drawn as one band (heat then cool at constant volume: one NVT stretch)
        var merged = new List<(double From, double To, string Label, bool Shade)>();
        foreach (var band in bands)
            if (merged.Count > 0 && merged[^1].Label == band.Item3) merged[^1] = (merged[^1].From, band.to, band.Item3, band.Item4);
            else merged.Add((band.from, band.to, band.Item3, band.Item4));
        _thermo.Clear();
        _thermo.AddRange(rows);
        EqStageBands = merged.ToArray();
        _eqAccepted = false;

        // convergence: the run's own checks, else the density and energy over the production blocks of the log
        var checks = new JsonArray();
        var blockPs = 0.0;
        var blockCount = 0;
        foreach (var l in logLines)
        {
            if (CheckLine.Match(l) is { Success: true } m)
                checks.Add(new JsonObject
                {
                    ["quantity"] = m.Groups[1].Value.Trim(), ["ok"] = m.Groups[2].Value == "ok",
                    ["change"] = double.Parse(m.Groups[3].Value, inv), ["tolerance"] = double.Parse(m.Groups[4].Value, inv),
                    ["blocks"] = new JsonArray(m.Groups[5].Value.Split(' ', StringSplitOptions.RemoveEmptyEntries).Select(v => (JsonNode)double.Parse(v, inv)).ToArray()),
                });
            if (BlocksLine.Match(l) is { Success: true } b) { blockCount = int.Parse(b.Groups[1].Value, inv); blockPs = double.Parse(b.Groups[2].Value, inv); }
        }
        var prod = starts.Where(s => s.Label.StartsWith("production", StringComparison.Ordinal)).ToList();
        if (blockPs <= 0 && prod.Count > 0)
        {
            var k = starts.IndexOf(prod[0]);
            var end = k + 1 < starts.Count ? rows[starts[k + 1].Row].TimePs : rows[^1].TimePs;
            blockPs = Math.Round(end - rows[prod[0].Row].TimePs + (rows[1].TimePs - rows[0].TimePs), 6);
            blockCount = prod.Count;
        }
        var fromLog = checks.Count > 0;
        if (!fromLog && prod.Count >= 2)
        {
            // as caps equilibrate checks them: block means, the last two changes within the tolerance
            var natoms = structure != null ? _doc?.Summary().Atoms ?? 0 : 0;   // the run's own structure, never another one open here
            (double[] D, double[] E) Means()
            {
                var d = new List<double>();
                var e = new List<double>();
                for (var k = 0; k < starts.Count; ++k)
                {
                    if (!starts[k].Label.StartsWith("production", StringComparison.Ordinal)) continue;
                    var end = k + 1 < starts.Count ? starts[k + 1].Row : rows.Count;
                    var seg = rows.Skip(starts[k].Row).Take(end - starts[k].Row).ToArray();
                    d.Add(seg.Average(r => r.Density));
                    e.Add(natoms > 0 ? seg.Average(r => r.Potential) / natoms : double.NaN);
                }
                return (d.ToArray(), e.ToArray());
            }
            var (dm, em) = Means();
            JsonObject Check(string q, double[] b, double tol, bool relative)
            {
                double Change(int i) => relative ? Math.Abs(b[i] - b[i - 1]) / Math.Max(1e-12, Math.Abs(b[i - 1])) : Math.Abs(b[i] - b[i - 1]);
                var ch = Change(b.Length - 1);
                var ok = b.Length >= 3 && ch <= tol && Change(b.Length - 2) <= tol;
                return new JsonObject { ["quantity"] = q, ["ok"] = ok, ["change"] = ch, ["tolerance"] = tol, ["blocks"] = new JsonArray(b.Select(v => (JsonNode)v).ToArray()) };
            }
            checks.Add(Check("density", dm, 0.005, true));
            if (natoms > 0) checks.Add(Check("potential energy per atom", em, 0.005, false));
        }
        // the run's protocol and blocks become the page's, so Extend continues it as it was run
        var protocol = logLines.Where(l => ProtocolLine.IsMatch(l)).ToArray();
        if (protocol.Length > 0) EqText = string.Join("\n", protocol) + "\n";
        if (prod.Count > 0 && blockPs > 0) { EqUntilConverged = true; EqBlockD = (decimal)blockPs; }
        var json = checks.Count > 0 ? new JsonObject { ["converged"] = checks.All(c => c!["ok"]!.GetValue<bool>()), ["blocks"] = blockCount, ["block_ps"] = blockPs, ["checks"] = checks }.ToJsonString() : "";
        LoadEqChecks(json);
        if (structure == null) EqEnergyBlocks = [];   // per atom in the checks: without the run's structure the count is not known

        // the report: the run's own (its stage table, checks and notes), else what was read
        var names = string.Join(" + ", new[] { csv, log, structure }.Where(f => f != null).Select(f => Path.GetFileName(f!)));
        var tableAt = Array.FindIndex(logLines, l => l.StartsWith("stage ", StringComparison.Ordinal));
        var report = tableAt >= 0 ? string.Join("\n", logLines.Skip(tableAt).Where(l => l.Length > 0)) : "";
        EqLog = $"Finished run read back: {names}\n" + (report.Length > 0 ? report
            : string.Format(inv, "{0} stages · {1:F1} ps · final density {2:F4} g/cm³{3}", starts.Count, rows[^1].TimePs, rows[^1].Density,
                            fromLog || checks.Count == 0 ? "" : "\nno log of the run beside the thermo: Rg and the internal distances are not checked; density and energy as caps equilibrate checks them"));
        EqOpenedFiles = names;
        EqOpenedText = string.Format(inv, "Finished run read back · {0} stages · {1:F0} ps{2}", starts.Count, rows[^1].TimePs, log == null ? " · no log beside it" : "");
        SetModule(4);
        ThermoChanged?.Invoke();
        EqChecksChanged?.Invoke();
        Raise(nameof(EqCanDecide));
        return Status = string.Format(inv, "Finished equilibration opened: {0} stages, {1:F0} ps, final density {2:F4} g/cm³{3}", starts.Count, rows[^1].TimePs, rows[^1].Density,
                                      EqCriteriaText.Length > 0 ? " · " + EqCriteriaText : "");
    }

    /// <summary>Opens cooling scans that ran elsewhere (caps tg): each a run folder, its OUT.tg-points.json, its csv folder or a
    /// table. One is a scan; several are replicas, pooled as the page pools its own (mean ± SD, error bars).</summary>
    public string OpenFinishedGlass(IReadOnlyList<string> paths)
    {
        if (paths.Count == 0) return Status;
        var inv = CultureInfo.InvariantCulture;
        var energy = Analyze.TgProperty == 1;
        var req = new JsonObject
        {
            ["paths"] = new JsonArray(paths.Select(p => (JsonNode)p).ToArray()), ["property"] = Analyze.TgProperty, ["fit"] = Analyze.TgFit,
            ["glassy_max"] = (double)Analyze.TgGlassyMaxD, ["rubbery_min"] = (double)Analyze.TgRubberyMinD,
        };
        JsonObject r;
        try { r = JsonNode.Parse(CapsDocument.TgFit(req.ToJsonString())) as JsonObject ?? new JsonObject(); }
        catch (Exception e) { return GtStatus = Status = "Could not read the cooling scans: " + e.Message; }
        if (r["ok"]?.GetValue<bool>() != true) { OpenGlass(); return GtStatus = Status = "Could not read the cooling scans: " + ((string?)r["error"] ?? "unknown error"); }
        var scans = (r["scans"] as JsonArray ?? []).OfType<JsonObject>().ToList();
        var props = (r["properties"] as JsonArray ?? []).OfType<JsonObject>().ToList();
        var each = props.Where(p => ((string?)p["id"] ?? "").StartsWith("tg_r", StringComparison.Ordinal)).ToList();
        if (each.Count == 0) each = props.Where(p => (string?)p["id"] == "tg").ToList();
        var runs = new List<(double Tg, double Err, double AlphaLow, double AlphaHigh, double Vtg, double SLow, double SHigh, double Rms, (double T, double V)[] Points)>();
        foreach (var p in each)
        {
            var curve = (p["series"] as JsonArray)?.OfType<JsonObject>().FirstOrDefault();
            var extra = p["extra"] as JsonObject ?? new JsonObject();
            double X(string key) => extra.Where(e => e.Key.StartsWith(key, StringComparison.Ordinal)).Select(e => Num(e.Value)).DefaultIfEmpty(double.NaN).First();
            runs.Add((Num(p["value"]), Num(p["error"]), X("expansion below"), X("expansion above"), X(energy ? "potential energy per atom at Tg" : "specific volume at Tg"),
                      X("slope below"), X("slope above"), X("fit residual"), Nums(curve?["x"]).Zip(Nums(curve?["y"])).ToArray()));
        }
        if (runs.Count == 0) { OpenGlass(); return GtStatus = Status = "No cooling points in what was opened"; }
        // the page's schedule as the runs had it (from their logs), so the left column describes them
        var first = scans[0];
        var pts = (first["points"] as JsonArray ?? []).OfType<JsonObject>().Select(q => Num(q["T"])).ToArray();
        var step = Num(first["t_step"]);
        var hold = Num(first["hold_ps"]);
        if (pts.Length > 1)
        {
            Analyze.TgFromD = (decimal)pts[0];
            Analyze.TgToD = (decimal)pts[^1];
            Analyze.TgStepD = (decimal)(step > 0 ? step : Math.Abs(pts[1] - pts[0]));
        }
        if (hold > 0) Analyze.TgPsD = (decimal)hold;
        if (Regex.Match((string?)first["schedule"] ?? "", @"after ([-0-9.eE+]+) ps") is { Success: true } am && double.TryParse(am.Groups[1].Value, NumberStyles.Float, inv, out var anneal))
            Analyze.EqPsD = (decimal)anneal;
        _gtReplicas = runs.Count;
        Raise(nameof(GtReplicas));
        string RunName(JsonObject s)
        {
            var src = (string?)s["source"] ?? "";
            var dir = Path.GetDirectoryName(src) ?? "";
            if (Path.GetFileName(dir) is "csv" or "out") dir = Path.GetDirectoryName(dir) ?? dir;
            return Path.GetFileName(dir);
        }
        var named = string.Join(", ", scans.Select(RunName));
        var withLog = scans.Count(s => ((string?)s["log"] ?? "").Length > 0);
        var rate = step > 0 && hold > 0 ? Sci(step / (hold * 1e-12)) + " K/s" : "a rate not known (no log)";
        OpenGlass();
        PoolGlass(runs, energy, Analyze.TgFitRanges, rate,
                  $"finished runs read back: {named}" + (withLog == scans.Count ? " · schedules from their logs" : withLog > 0 ? $" · {withLog} with a log" : " · no logs beside them"));
        _gtOpenedPaths = paths.ToArray();
        GtOpenedText = $"Finished run{(scans.Count == 1 ? "" : "s")} read back · {named}";
        return Status = $"{(runs.Count == 1 ? "A finished cooling scan" : $"{runs.Count} finished replicas")} opened: Tg {GtTg.Value}";
    }

    /// <summary>The finished scans again, with the page's Tg settings as they are now (property, fit, ranges).</summary>
    public string RefitFinishedGlass() => _gtOpenedPaths.Count > 0 ? OpenFinishedGlass(_gtOpenedPaths) : Status;

    /// <summary>A finished cluster job of Equilibrate or Glass shown on its page from the files it brought back.</summary>
    public bool CanShowJobOnPage(Job? j) => j?.Remote is { } r && j.IsDone && j.Kind is "Equilibrate" or "Glass" && Directory.Exists(r.Local);

    public void ShowJobOnPage(Job? j)
    {
        if (j?.Remote is not { } r || !CanShowJobOnPage(j)) return;
        if (j.Kind == "Glass")
        {
            if (r.Batch.Length > 0) CollectGlassReplicas(r.Batch);   // the batch's replicas pooled, those that never came back named
            else OpenFinishedGlass([r.Local]);
        }
        else OpenFinishedEquilibrate([r.Local]);
    }
}
