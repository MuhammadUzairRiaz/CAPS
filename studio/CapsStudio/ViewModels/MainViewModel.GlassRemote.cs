using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Text.Json.Nodes;
using System.Threading.Tasks;

namespace CapsStudio.ViewModels;

/// <summary>Glass transition on a remote host: one recipe job per replica (seed 1, 2, …), each with the structure and the
/// force field as assigned here (structure.ff.json, nothing typed again there), and the replicas pooled when they are
/// back — the same mean ± SD, fit and slopes as a run on this machine.</summary>
public sealed partial class MainViewModel
{
    public async Task SubmitGlassRemote()
    {
        if (_doc == null || RunWhereIndex == 0) return;
        var h = _settings.Hosts[RunWhereIndex - 1];
        if (h.Hostname.Length == 0) { Status = "The host has no hostname: set it in Settings › Compute & remote"; return; }
        var baseStem = string.Concat(Title.Replace(" (unsaved)", "").Where(char.IsLetterOrDigit).Take(24)) is { Length: > 0 } t ? t : "structure";
        var batch = $"glass-{DateTime.Now:yyyyMMdd-HHmmss}";
        var sent = 0;
        for (var rep = 1; rep <= _gtReplicas; ++rep)
        {
            var k = _jobCounters["Glass"] = _jobCounters.GetValueOrDefault("Glass") + 1;
            var id = $"glass-{k}";
            var stem = $"{baseStem}_tg_r{rep}";
            var job = new Job
            {
                Id = id, Kind = "Glass", Module = 47, Title = $"Glass transition · replica {rep} of {_gtReplicas} · on {h.Name}",
                Document = Title, Atoms = _doc.Summary().Atoms, Provenance = Manifest("Glass"),
            };
            job.Status = "queued";
            job.Remote = new RemoteRun { Host = h.Name, Scheduler = h.Scheduler, Local = Path.Combine(RemoteFolder, $"{id}-{DateTime.Now:yyyyMMdd-HHmmss}"), Stem = stem, Batch = batch };
            Jobs.Insert(0, job);
            Raise(nameof(HasJobs)); Raise(nameof(JobsSummary));
            try
            {
                var local = job.Remote.Local;
                Directory.CreateDirectory(local);
                _doc.Save(Path.Combine(local, "structure.caps.data"));
                var files = new List<string> { "structure.caps.data", "recipe.yaml", "job.sh" };
                string? ff = null;
                if (Field.Assigned) { files.AddRange(SaveForceFieldForHost(local)); ff = "structure.ff.json"; }
                var recipe = GlassRecipe("structure.caps.data", ff, (ulong)rep, stem);
                File.WriteAllText(Path.Combine(local, "recipe.yaml"), recipe);
                await SendPrepared(job, h, id, files.ToArray(), "recipe.yaml");
                ++sent;
            }
            catch (Exception e)
            {
                job.Status = "failed";
                job.Error = e.Message;
                job.Add("Could not send the job: " + e.Message);
                job.Suggestion = "Test the host in Settings › Compute & remote (keys come from your SSH agent; caps must be on the host's PATH).";
                job.SuggestModule = 10;
                job.Ended = DateTime.Now;
                Status = "Could not send the job: " + e.Message;
                break;
            }
        }
        if (sent > 0)
        {
            GtStatus = $"{sent} replica{(sent == 1 ? "" : "s")} sent to {h.Name}{(Field.Assigned ? " with the force field as assigned here" : " (typed there with the default force field)")}; " +
                       "pooled here when all are back (Jobs)";
            Status = GtStatus;
            StartRemotePoll();
        }
        SaveJobs();
        Raise(nameof(JobsSummary)); Raise(nameof(ComputeText));
    }

    /// <summary>The batch's finished replicas read from their properties files (out/STEM.properties.json) and pooled onto
    /// the Glass page; a replica that failed or has not come back is left out and named.</summary>
    public void CollectGlassReplicas(string batch)
    {
        var jobs = Jobs.Where(j => j.Remote?.Batch == batch).OrderBy(j => j.Remote!.Stem, StringComparer.Ordinal).ToList();
        var runs = new List<(double Tg, double Err, double AlphaLow, double AlphaHigh, double Vtg, double SLow, double SHigh, double Rms, (double T, double V)[] Points)>();
        var missing = new List<string>();
        bool? energy = null;
        string rate = "", method = "";
        foreach (var j in jobs)
        {
            var file = Path.Combine(j.Remote!.Local, "out", j.Remote.Stem + ".properties.json");
            var tg = File.Exists(file) ? ReadTgProperty(file) : null;
            if (tg == null) { missing.Add(j.Remote.Stem); continue; }
            var series = tg["series"] as JsonArray;
            var curve = series?.OfType<JsonObject>().FirstOrDefault();
            if (curve == null) { missing.Add(j.Remote.Stem); continue; }
            var label = (string?)curve["label"] ?? "";
            energy ??= label.StartsWith("potential energy", StringComparison.Ordinal);
            var xs = Nums(curve["x"]);
            var ys = Nums(curve["y"]);
            var extra = tg["extra"] as JsonObject ?? new JsonObject();
            double X(string key) => extra.Where(e => e.Key.StartsWith(key, StringComparison.Ordinal)).Select(e => Num(e.Value)).DefaultIfEmpty(double.NaN).First();
            runs.Add((Num(tg["value"]), Num(tg["error"]), X("expansion below"), X("expansion above"), X(energy == true ? "potential energy per atom at Tg" : "specific volume at Tg"),
                      X("slope below"), X("slope above"), X("fit residual"), xs.Zip(ys).ToArray()));
            method = (string?)tg["method"] ?? method;
            if (rate.Length == 0 && tg["notes"] is JsonArray notes)
                foreach (var n in notes.Select(n => (string?)n ?? ""))
                    if (n.StartsWith("effective cooling rate ", StringComparison.Ordinal) &&
                        double.TryParse(n["effective cooling rate ".Length..].Split(' ')[0], NumberStyles.Float, CultureInfo.InvariantCulture, out var v))
                        rate = Sci(v) + " K/s";
        }
        OpenGlass();
        if (runs.Count == 0) { GtStatus = "No replica of this batch has come back with a result" + (missing.Count > 0 ? ": " + string.Join(", ", missing) : ""); return; }
        var fitRanges = method.Contains("glassy and the rubbery range", StringComparison.Ordinal);
        var log = $"on {jobs[0].Remote!.Host}" + (missing.Count > 0 ? $" · left out (no result): {string.Join(", ", missing)}" : "");
        PoolGlass(runs, energy == true, fitRanges, rate.Length > 0 ? rate : GtRateText, log);
        Status = $"Glass transition pooled from {runs.Count} remote replica{(runs.Count == 1 ? "" : "s")}";
    }

    private static JsonObject? ReadTgProperty(string file)
    {
        try { return (JsonNode.Parse(File.ReadAllText(file)) as JsonArray)?.OfType<JsonObject>().FirstOrDefault(p => (string?)p["id"] == "tg"); }
        catch (Exception) { return null; }
    }
    private static double Num(JsonNode? n) => n is JsonValue v && v.TryGetValue<double>(out var d) ? d : double.NaN;
    private static double[] Nums(JsonNode? n) => n is JsonArray a ? a.Select(Num).ToArray() : [];
}
