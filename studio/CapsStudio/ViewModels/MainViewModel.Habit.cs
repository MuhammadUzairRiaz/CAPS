using System.Globalization;
using System.Text;
using System.Text.Json.Nodes;

namespace CapsStudio.ViewModels;

/// <summary>Builders › Crystal › crystal shape (core habit.hpp): the open crystal's faces by BFDH, attachment and surface
/// energies with the assigned force field, the Wulff shapes, and a particle carved in one of them.</summary>
public sealed partial class MainViewModel
{
    public static string[] HabitModels { get; } = ["Spacing (BFDH)", "Growth (attachment energy)", "Equilibrium (surface energy)"];
    private static readonly string[] HabitModelIds = ["bfdh", "growth", "equilibrium"];
    private bool _habRunning, _habHas, _habEnergies, _habSurface;
    private string _habError = "", _habTable = "", _habNote = "", _habShapeText = "";
    private decimal _habIndex = 3, _habRadius = 20;
    private int _habModel;
    private JsonNode? _habLast;
    public bool HabitRunning { get => _habRunning; private set { if (Set(ref _habRunning, value)) RaiseBusy(); } }
    public bool HabitHas { get => _habHas; private set => Set(ref _habHas, value); }
    public bool HabitEnergies { get => _habEnergies; set => Set(ref _habEnergies, value); }
    public bool HabitSurface { get => _habSurface; set => Set(ref _habSurface, value); }
    public decimal? HabitMaxIndex { get => _habIndex; set { if (value != null) Set(ref _habIndex, Math.Clamp(value.Value, 1, 6)); } }
    public decimal? HabitRadius { get => _habRadius; set { if (value != null) Set(ref _habRadius, Math.Clamp(value.Value, 5, 100)); } }
    public int HabitModel { get => _habModel; set { if (Set(ref _habModel, value)) ShowHabitShape(); } }
    public string HabitError { get => _habError; private set { if (Set(ref _habError, value)) Raise(nameof(HabitHasError)); } }
    public bool HabitHasError => _habError.Length > 0;
    public string HabitTable { get => _habTable; private set => Set(ref _habTable, value); }
    public string HabitNote { get => _habNote; private set => Set(ref _habNote, value); }
    public string HabitShapeText { get => _habShapeText; private set => Set(ref _habShapeText, value); }
    /// <summary>The shape on show: vertices, facets, each facet's form, the forms' labels.</summary>
    public (double[][] V, int[][] F, int[] Form, string[] Labels) HabitShape { get; private set; } = ([], [], [], []);
    public event Action? HabitChanged;

    public async Task RunHabit()
    {
        if (_doc == null || !Idle) return;
        var energies = _habEnergies || _habSurface;
        if (energies && BlockedByField("Crystal shape")) return;
        var json = new JsonObject
        {
            ["max_index"] = (int)_habIndex, ["attachment"] = _habEnergies, ["surface"] = _habSurface, ["cutoff"] = Math.Min(12.0, _relaxCutoff),
        }.ToJsonString();
        var doc = _doc!;
        HabitRunning = true;
        HabitError = "";
        Status = "Crystal shape…";
        try
        {
            var res = await Task.Run(() => doc.Habit(json, (st, f) =>
            {
                Avalonia.Threading.Dispatcher.UIThread.Post(() => { if (HabitRunning) Status = $"Crystal shape: {st} · {f * 100:0} %"; });
                return true;
            }));
            var r = JsonNode.Parse(res)!;
            if (r["ok"]?.GetValue<bool>() != true) { HabitError = (string?)r["error"] ?? "habit failed"; Status = "Crystal shape: " + HabitError; return; }
            _habLast = r;
            var inv = CultureInfo.InvariantCulture;
            string N(JsonNode? x, string f) => x is JsonNode v ? ((double?)v ?? 0).ToString(f, inv) : "–";
            var sb = new StringBuilder();
            sb.AppendLine($"{"form",-9}{"d_eff",7}{"E_att",8}{"γ J/m²",8}{"BFDH",6}{"grow",5}{"eq.",5}");
            foreach (var f in (JsonArray)r["faces"]!)
            {
                var shares = new[] { "area_bfdh", "area_growth", "area_equilibrium" }.Select(k => (double?)f![k] ?? 0).ToArray();
                if ((int)((double?)f!["rank"] ?? 0) > 8 && shares.All(x => x == 0)) continue;
                var hkl = string.Join(" ", ((JsonArray)f["hkl"]!).Select(x => (int?)x ?? 0));
                sb.AppendLine(string.Format(inv, "{0,-9}{1,7:0.00}{2,8}{3,8}{4,5:0}%{5,4:0}%{6,4:0}%", "{" + hkl + "}", (double?)f["d_eff"] ?? 0, N(f["e_att"], "0.00"),
                    N(f["gamma"], "0.000"), 100 * shares[0], 100 * shares[1], 100 * shares[2]));
            }
            HabitTable = sb.ToString().TrimEnd();
            var notes = ((JsonArray)r["faces"]!).Where(f => ((string?)f!["note"] ?? "").Length > 0)
                .Select(f => "{" + string.Join(" ", ((JsonArray)f!["hkl"]!).Select(x => (int?)x ?? 0)) + "}: " + (string?)f["note"]).Take(4)
                .Concat(((JsonArray)r["notes"]!).Select(v => (string?)v ?? ""));
            HabitNote = $"{(string?)r["space_group"]} · {(int)((double?)r["units"] ?? 0)} growth units per cell" + (r["e_latt"] is JsonNode el ? $" · bulk {((double?)el ?? 0).ToString("0.00", inv)} kcal/mol per unit" : "") +
                        "\n" + string.Join("\n", notes);
            if (_habModel > 0 && r["shapes"]![HabitModelIds[_habModel]] == null) _habModel = 0;
            Raise(nameof(HabitModel));
            HabitHas = true;
            ShowHabitShape();
            Status = $"Crystal shape: {((JsonArray)r["faces"]!).Count} forms";
        }
        catch (Exception e) { HabitError = e.Message; }
        finally { HabitRunning = false; }
    }

    private void ShowHabitShape()
    {
        if (_habLast is not JsonNode r) return;
        var sh = r["shapes"]![HabitModelIds[_habModel]];
        if (sh == null) { HabitShape = ([], [], [], []); HabitShapeText = "no " + HabitModels[_habModel].ToLowerInvariant() + " shape: compute its energies first"; HabitChanged?.Invoke(); return; }
        var faces = (JsonArray)r["faces"]!;
        var v = ((JsonArray)sh["vertices"]!).Select(p => ((JsonArray)p!).Select(x => (double?)x ?? 0).ToArray()).ToArray();
        var f = ((JsonArray)sh["facets"]!).Select(p => ((JsonArray)p!).Select(x => (int?)x ?? 0).ToArray()).ToArray();
        var form = ((JsonArray)sh["facet_face"]!).Select(x => (int?)x ?? -1).ToArray();
        var labels = faces.Select(x => "{" + string.Join("", ((JsonArray)x!["hkl"]!).Select(q => (int?)q ?? 0)) + "}").ToArray();
        HabitShape = (v, f, form, labels);
        HabitShapeText = string.Format(CultureInfo.InvariantCulture, "{0} facets · aspect ratio {1:0.00}", f.Length, (double?)sh["aspect"] ?? 0);
        HabitChanged?.Invoke();
    }

    /// <summary>A particle carved in the shape on show, opened as a new structure.</summary>
    public async Task HabitParticle()
    {
        if (_doc == null || !Idle || !_habHas) return;
        var path = Path.Combine(Path.GetTempPath(), $"caps-particle-{Environment.ProcessId}-{DateTime.Now:HHmmss}.pdb");
        var json = new JsonObject
        {
            ["max_index"] = (int)_habIndex, ["attachment"] = _habModel == 1, ["surface"] = _habModel == 2, ["cutoff"] = Math.Min(12.0, _relaxCutoff),
            ["particle"] = new JsonObject { ["model"] = HabitModelIds[_habModel], ["radius"] = (double)_habRadius, ["path"] = path },
        }.ToJsonString();
        var doc = _doc!;
        HabitRunning = true;
        HabitError = "";
        try
        {
            var r = JsonNode.Parse(await Task.Run(() => doc.Habit(json)))!;
            if (r["ok"]?.GetValue<bool>() != true) { HabitError = (string?)r["error"] ?? "particle failed"; return; }
            Status = $"Particle: {(int)((double?)r["particle"]!["atoms"] ?? 0)} atoms";
        }
        catch (Exception e) { HabitError = e.Message; return; }
        finally { HabitRunning = false; }
        await OpenProgressive(path, null);
    }
}
