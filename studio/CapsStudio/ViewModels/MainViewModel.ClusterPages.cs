using System.Globalization;
using System.Text.Json.Nodes;

namespace CapsStudio.ViewModels;

/// <summary>Runs on a cluster from the pages whose engine is a caps command (Mechanics, pull-out, CBMC, Pack, the
/// Coarse-grain stages) and sweeps as one array job: the command the page stands for goes up with every local file and
/// folder it names (renamed to the uploaded copies), and the outputs come back — folders into the places they were meant
/// for on this machine.</summary>
public sealed partial class MainViewModel
{
    private static string N(double v) => v.ToString("0.######", Inv);
    private static string N(decimal v) => v.ToString("0.######", Inv);

    /// <summary>Sends a caps command to the host chosen in Run where. build gets a function that registers a local file
    /// or folder for upload and returns the name the command uses on the host; fetch: output folders (names in the job's
    /// results) and the local folders they come back into.</summary>
    internal async Task<Job?> SubmitCommandJob(string kind, int module, Func<Func<string, string>, string> build, bool withStructure,
                                              List<(string Name, string Local)>? fetch = null, string? stem = null)
    {
        if (RemoteHostNow is not { } h) return null;
        if (withStructure && _doc == null) { Status = "Open a structure first"; return null; }
        if (h.Hostname.Length == 0) { Status = "The host has no hostname: set it in Settings › Compute & remote"; return null; }
        var k = _jobCounters[kind] = _jobCounters.GetValueOrDefault(kind) + 1;
        var id = $"{kind.ToLowerInvariant()}-{k}";
        stem ??= string.Concat(Title.Replace(" (unsaved)", "").Where(char.IsLetterOrDigit).Take(24)) is { Length: > 0 } t ? t : "structure";
        var job = new Job { Id = id, Kind = kind, Module = module >= 0 ? module : _module, Title = $"{kind} · on {h.Name}", Document = withStructure ? Title : "", Atoms = withStructure ? _doc!.Summary().Atoms : 0, Provenance = Manifest(kind) };
        job.Status = "queued";
        job.Remote = new RemoteRun { Host = h.Name, Scheduler = h.Scheduler, Local = Path.Combine(RemoteFolder, $"{id}-{DateTime.Now:yyyyMMdd-HHmmss}"), Stem = stem };
        if (fetch != null) job.Remote.Fetch = fetch.Select(f => f.Name + "|" + f.Local).ToList();
        Jobs.Insert(0, job);
        SelectedJob = job;
        Raise(nameof(HasJobs)); Raise(nameof(JobsSummary));
        try
        {
            var local = job.Remote.Local;
            Directory.CreateDirectory(local);
            var files = new List<string>();
            if (withStructure) { _doc!.Save(Path.Combine(local, "structure.caps.data")); files.Add("structure.caps.data"); }
            var names = new Dictionary<string, string>(StringComparer.Ordinal);
            string Upload(string path)
            {
                var full = Path.GetFullPath(path);
                if (names.TryGetValue(full, out var known)) return known;
                var name = Path.GetFileName(full.TrimEnd(Path.DirectorySeparatorChar));
                for (var n = 2; files.Contains(name); ++n) name = Path.GetFileNameWithoutExtension(full) + "_" + n + Path.GetExtension(full);
                if (Directory.Exists(full)) CopyFolder(full, Path.Combine(local, name));
                else File.Copy(full, Path.Combine(local, name), true);
                files.Add(name);
                names[full] = name;
                return name;
            }
            var command = build(Upload);
            await SendJob(job, h, kind, files.ToArray(), command);
            Status = $"{kind} sent to {h.Name} · {(h.Scheduler == "none" ? "process" : h.Scheduler)} {job.Remote.JobId} · Jobs follows it";
            StartRemotePoll();
        }
        catch (Exception e)
        {
            job.Status = "failed";
            job.Error = e.Message;
            job.Add("Could not send the job: " + e.Message);
            job.Suggestion = "Test the host in Settings › Compute & remote (keys come from your SSH agent; Install / update caps puts caps on the host).";
            job.SuggestModule = 10;
            job.Ended = DateTime.Now;
            Status = "Could not send the job: " + e.Message;
        }
        SaveJobs();
        Raise(nameof(JobsSummary)); Raise(nameof(ComputeText));
        return job;
    }

    private static void CopyFolder(string from, string to)
    {
        Directory.CreateDirectory(to);
        foreach (var f in Directory.GetFiles(from)) File.Copy(f, Path.Combine(to, Path.GetFileName(f)), true);
        foreach (var d in Directory.GetDirectories(from)) CopyFolder(d, Path.Combine(to, Path.GetFileName(d)));
    }

    // ---------------------------------------------------------------- the pages' commands
    /// <summary>Mechanics: the stress–strain pull and the elastic constants chosen on the page, as caps tensile / elastic.</summary>
    public Task<Job?> SubmitMechanicsRemote()
    {
        var a = Analyze;
        var cmds = new List<string>();
        if (_mechTensile)
            cmds.Add($"caps tensile structure.caps.data -o tensile.data --axis {(a.TensAxis >= 3 ? "xyz" : AxisName(a.TensAxis))} --rate {N(a.TensRateD)} --strain {N(a.TensMaxD)} " +
                     $"--temp {N(a.TensTD)}{(a.TensFixed ? " --fixed-lateral" : "")} --json tensile.json --csv curves");
        cmds.Add($"caps elastic structure.caps.data --method {(_mechMethod == 0 ? "strain" : "fluct-run")} --strain {N(a.CijStrainD)} --configs {N(a.CijConfigsD)} --temp {N(a.TensTD)} --json elastic.json");
        return SubmitCommandJob("Mechanics", -1, _ => string.Join(" && ", cmds), true, stem: "tensile");
    }

    /// <summary>The pull-out of a film from its surface (Interface page), as caps pull.</summary>
    public Task<Job?> SubmitPullRemote()
    {
        var a = Analyze;
        var how = a.PullNormalChip.IsOn ? "--normal" : $"--axis {"xyz"[Math.Clamp(a.PullAxis, 0, 2)]}";
        var surface = a.SurfaceMolecules.Trim().Split(',', '-')[0] is { Length: > 0 } m ? m : "1";
        return SubmitCommandJob("Pull", -1, _ => $"caps pull structure.caps.data {how} --distance {N(a.PullDistD)} --rate {N(a.PullRateD)} --temp {N(a.PullTD)} " +
                                                 $"--eq-ps {N(a.PullEqD)} --surface {surface} --csv pull.csv", true);
    }

    /// <summary>Configurational-bias regrowth of chain ends (Equilibrate page), as caps cbmc.</summary>
    public Task<Job?> SubmitCbmcRemote() =>
        SubmitCommandJob("CBMC", 4, _ => $"caps cbmc structure.caps.data -o regrown.data --moves {_cbMoves} --trials {_cbTrials} --max-torsions {_cbTorsions} " +
                                         $"--temp {N(_cbTemp)} --cutoff {N(Math.Min(9.0, _relaxCutoff))} --seed {MdSeedChoice.Take()}{(_relaxCoulomb ? "" : " --no-coulomb")}", true, stem: "regrown");

    /// <summary>Pack: the input as written on the page, its molecule files uploaded and named by their copies.</summary>
    public Task<Job?> SubmitPackRemote()
    {
        if (_packText.Trim().Length == 0) { Status = "Pack: add molecules first"; return Task.FromResult<Job?>(null); }
        return SubmitCommandJob("Pack", 5, upload =>
        {
            var lines = _packText.Replace("\r\n", "\n").Split('\n').Select(line =>
            {
                var words = line.Split(' ', StringSplitOptions.None);
                for (var i = 0; i < words.Length; ++i)
                {
                    var w = words[i].Trim('"');
                    if (w.Length > 2 && (w.Contains('/') || w.Contains('\\') || w.Contains('.')) && File.Exists(w)) words[i] = upload(w);
                }
                return string.Join(' ', words);
            });
            var local = Path.Combine(RemoteFolder, "pack-" + DateTime.Now.ToString("yyyyMMddHHmmssfff", Inv) + ".inp");
            Directory.CreateDirectory(RemoteFolder);
            File.WriteAllText(local, string.Join("\n", lines) + "\n");
            return $"caps pack {upload(local)} -o packed.data";
        }, false, stem: "packed");
    }

    /// <summary>A Coarse-grain stage on the cluster: its command with the same arguments (local files and folders uploaded),
    /// its output folder brought back to where the stage writes it here.</summary>
    internal Task<Job?> SubmitCgRemote(string command, JsonObject args)
    {
        var fetch = new List<(string, string)>();
        return SubmitCommandJob("CG " + command, 80, upload =>
        {
            string Word(string v) => v.Length > 0 && v.IndexOfAny([' ', '\'', '"', '$', '`', '\\', '*', '?', ';', '&', '|', '<', '>', '(', ')']) < 0 ? v : Q(v);
            string Value(string v) => (File.Exists(v) || Directory.Exists(v)) ? upload(v)
                                      : string.Join(",", v.Split(',').Select(p => File.Exists(p) || Directory.Exists(p) ? upload(p) : p));
            var parts = new List<string> { "caps", command };
            if (args["inputs"] is JsonArray ins) parts.AddRange(ins.Select(x => Word(Value((string?)x ?? ""))));
            foreach (var (key, node) in args)
            {
                if (key == "inputs" || node == null) continue;
                var flag = key == "o" ? "-o" : key == "T" ? "-T" : "--" + key.Replace('_', '-');
                if (node is JsonValue jv && jv.TryGetValue<bool>(out var b)) { if (b) parts.Add(flag); continue; }
                var v = node is JsonValue sv && sv.TryGetValue<string>(out var s) ? s : node.ToJsonString();
                if (key == "o")
                {   // the output folder: made on the host under its own name, brought back here
                    var name = Path.GetFileName(v.TrimEnd('/', '\\')) is { Length: > 0 } n ? n : "out";
                    fetch.Add((name, v));
                    parts.Add(flag); parts.Add(Word(name));
                    continue;
                }
                parts.Add(flag); parts.Add(Word(Value(v)));
            }
            return string.Join(' ', parts);
        }, false, fetch, stem: command);
    }

    /// <summary>A sweep's cells as one array job on the host: one task per cell (its recipe), each followed in Jobs.</summary>
    internal async Task<List<Job>> SendRecipeArray(RemoteHost h, string kind, List<(string Title, string Stem, string Recipe)> cells)
    {
        var jobs = new List<Job>();
        var local = Path.Combine(RemoteFolder, $"array-{DateTime.Now:yyyyMMdd-HHmmss}");
        Directory.CreateDirectory(local);
        var files = new List<string>();
        var list = new System.Text.StringBuilder();
        for (var i = 0; i < cells.Count; ++i)
        {
            var f = $"recipe_{i + 1}.json";
            File.WriteAllText(Path.Combine(local, f), cells[i].Recipe);
            files.Add(f);
            var title = string.Concat(cells[i].Stem.Select(c => char.IsLetterOrDigit(c) || c is '-' or '_' or '.' ? c : '_'));
            list.Append(title).Append(" caps run ").Append(f).Append(" --out .\n");
        }
        File.WriteAllText(Path.Combine(local, "tasks.list"), list.ToString());
        await WriteHostProfile(h);
        var stage = $".uploads/array-{DateTime.Now:yyyyMMddHHmmss}";
        var mk = await Tool("ssh", SshArgs(h, $"mkdir -p {RootPath(h, stage)} && cd {RootPath(h, stage)} && pwd"), 30000);
        if (mk.Code != 0) throw new InvalidOperationException("ssh: " + FirstLine(mk.Err, mk.Code));
        var dir = mk.Out.Split('\n').Last().Trim();
        var up = await Tool("scp", ScpArgs(h, files.Append("tasks.list").Select(f => Path.Combine(local, f)), $"{Target(h)}:{dir}/"), 300000);
        if (up.Code != 0) throw new InvalidOperationException("scp: " + FirstLine(up.Err, up.Code));
        var opts = $"--kind {Q(kind.ToLowerInvariant())} --array {Q(dir + "/tasks.list")} --input {Q(string.Join(",", files.Select(f => dir + "/" + f)))} --cpus {RemoteCpus.ToString("0", Inv)}" +
                   (RemoteMem.Length > 0 ? $" --mem {Q(RemoteMem)}" : "") + (RemoteTime.Length > 0 ? $" --time {Q(RemoteTime)}" : "") + (h.Partition.Length > 0 ? $" --partition {Q(h.Partition)}" : "");
        var sub = await Tool("ssh", SshArgs(h, $"{CapsOnHost(h)} job new --host-profile {RootPath(h, "host.json")} {opts} --json --submit; RC=$?; rm -rf {Q(dir)}; exit $RC"), 120000);
        var line = sub.Out.Split('\n').LastOrDefault(l => l.TrimStart().StartsWith('{')) ?? "";
        if (JsonNode.Parse(line) is not JsonObject res || res["submitted"]?.GetValue<bool>() != true)
            throw new InvalidOperationException("caps job new --array: " + FirstLine(sub.Err.Length > 0 ? sub.Err : sub.Out, sub.Code));
        var id = (string?)res["id"] ?? "";
        var tasks = (res["tasks"] as JsonArray ?? []).Select(x => (string?)x ?? "").ToList();
        for (var i = 0; i < cells.Count && i < tasks.Count; ++i)
        {
            var k = _jobCounters[kind] = _jobCounters.GetValueOrDefault(kind) + 1;
            var job = new Job { Id = $"{kind.ToLowerInvariant()}-{k}", Kind = kind, Module = 0, Title = $"{cells[i].Title} · on {h.Name}", Document = "new cell", Provenance = Manifest(kind) };
            job.Remote = new RemoteRun
            {
                Host = h.Name, Scheduler = h.Scheduler, Mode = "job", Dir = tasks[i], Stem = cells[i].Stem,
                JobId = h.Scheduler == "SLURM" ? $"{id}_{i}" : h.Scheduler == "PBS" ? $"{id.Split('[')[0]}[{i}]" : id,
                Local = Path.Combine(RemoteFolder, $"{job.Id}-{DateTime.Now:yyyyMMdd-HHmmss}"),
            };
            Directory.CreateDirectory(job.Remote.Local);
            job.Status = "queued";
            job.Add($"Array task {i} of {id} · {tasks[i]}");
            Jobs.Insert(0, job);
            jobs.Add(job);
        }
        Raise(nameof(HasJobs)); Raise(nameof(JobsSummary));
        SaveJobs();
        StartRemotePoll();
        return jobs;
    }
}
