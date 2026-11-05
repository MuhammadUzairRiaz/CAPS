using System.Globalization;
using System.Text.Json.Nodes;

namespace CapsStudio.ViewModels;

/// <summary>Relax › Reaction path (core reaction_path.hpp): from the open structure (A) to a second file (B), the atoms
/// matched, the images interpolated (IDPP, LST or straight lines) and, with the climbing-image NEB, the barrier with the
/// assigned force field (UFF when none). The images open as a new structure, one frame each.</summary>
public sealed partial class MainViewModel
{
    public static readonly string[] PathMethods = ["IDPP (pair distances, as a band)", "LST (pair distances, image by image)", "Straight lines"];
    private static readonly string[] PathMethodIds = ["idpp", "lst", "linear"];
    private string _pathB = "", _pathResult = "", _pathError = "";
    private int _pathMethod;
    private decimal _pathImages = 7, _pathFmax = 0.05m;
    private bool _pathNeb = true, _pathRelaxEnds = true, _pathMatch = true, _pathRunning;
    public string PathB { get => _pathB; set { if (Set(ref _pathB, value ?? "")) Raise(nameof(PathHasB)); } }
    public bool PathHasB => _pathB.Length > 0;
    public string PathBName => _pathB.Length == 0 ? "no end B chosen" : Path.GetFileName(_pathB);
    public int PathMethod { get => _pathMethod; set => Set(ref _pathMethod, Math.Clamp(value, 0, 2)); }
    public decimal PathImages { get => _pathImages; set => Set(ref _pathImages, Math.Clamp(Math.Round(value), 1, 40)); }
    public decimal PathFmax { get => _pathFmax; set => Set(ref _pathFmax, Math.Clamp(value, 0.001m, 5)); }
    public bool PathNeb { get => _pathNeb; set => Set(ref _pathNeb, value); }
    public bool PathRelaxEnds { get => _pathRelaxEnds; set => Set(ref _pathRelaxEnds, value); }
    public bool PathMatch { get => _pathMatch; set => Set(ref _pathMatch, value); }
    public string PathResult { get => _pathResult; private set => Set(ref _pathResult, value); }
    public string PathError { get => _pathError; private set { if (Set(ref _pathError, value)) Raise(nameof(PathHasError)); } }
    public bool PathHasError => _pathError.Length > 0;
    public bool PathRunning { get => _pathRunning; private set { if (Set(ref _pathRunning, value)) RaiseBusy(); } }

    public void ChoosePathB(string path) { PathB = path; Raise(nameof(PathBName)); }

    public async Task RunReactionPath()
    {
        if (_doc == null || !Idle || _pathB.Length == 0) return;
        var doc = _doc;
        var inv = CultureInfo.InvariantCulture;
        var json = new JsonObject
        {
            ["other"] = _pathB, ["match"] = _pathMatch, ["method"] = PathMethodIds[_pathMethod], ["images"] = (int)_pathImages,
            ["neb"] = _pathNeb, ["relax_ends"] = _pathRelaxEnds, ["fmax"] = (double)_pathFmax,
        }.ToJsonString();
        PathRunning = true;
        PathError = "";
        Status = _pathNeb ? "Running the nudged elastic band…" : "Interpolating the path…";
        try
        {
            var title = $"{Path.GetFileNameWithoutExtension(_pathB)} path ({PathMethodIds[_pathMethod]})";
            var (nd, rep) = await Task.Run(() => doc.ReactionPath(json, title));
            var r = JsonNode.Parse(rep)!;
            if (nd == null || r["ok"]?.GetValue<bool>() != true) { PathError = (string?)r["error"] ?? "the path failed"; Status = "Reaction path: " + PathError; return; }
            double D(string k) => r[k]?.GetValue<double>() ?? 0;
            var head = $"{(int)D("reordered")} atoms of B in another order (RMSD {D("match_rmsd").ToString("0.000", inv)} Å)";
            if (_pathNeb)
            {
                var conv = r["converged"]?.GetValue<bool>() == true ? "converged" : "not converged";
                PathResult = $"{head} · barrier {D("forward").ToString("0.00", inv)} kcal/mol from A, {D("reverse").ToString("0.00", inv)} from B · climbing image {(int)D("climbing")} · {conv} in {(int)D("steps")} steps"
                             + (r["forcefield"] is JsonNode f ? " · " + (string?)f : "");
            }
            else
            {
                var md = ((JsonArray)r["min_distance"]!).Select(x => x!.GetValue<double>()).DefaultIfEmpty(0).Min();
                PathResult = $"{head} · shortest interatomic distance on the path {md.ToString("0.00", inv)} Å";
            }
            var notes = ((JsonArray)r["notes"]!).Select(x => (string?)x ?? "").Where(x => x.StartsWith("no image", StringComparison.Ordinal));
            if (notes.Any()) PathResult += " · " + string.Join(" · ", notes);
            PathRunning = false;
            Show(nd, title);
            Status = "Reaction path: " + PathResult;
        }
        catch (Exception e) { PathError = e.Message; }
        finally { PathRunning = false; }
    }
}
