using System.Globalization;
using System.Text;
using System.Text.Json.Nodes;

namespace CapsStudio.ViewModels;

/// <summary>Builders › Crystal › substitutional disorder (core sqs.hpp): a quasirandom alloy on one element's sites, or the
/// symmetry-distinct arrangements of a few substituents with how often each occurs.</summary>
public sealed partial class MainViewModel
{
    private string _disFrom = "", _disTo = "", _disError = "", _disText = "";
    private decimal _disPercent = 50, _disCount = 1, _disPick = 1;
    private bool _disRunning;
    private string[] _disFiles = [];
    public string DisorderFrom { get => _disFrom; set => Set(ref _disFrom, value.Trim()); }
    public string DisorderTo { get => _disTo; set => Set(ref _disTo, value.Trim()); }
    public decimal? DisorderPercent { get => _disPercent; set { if (value != null) Set(ref _disPercent, Math.Clamp(value.Value, 1, 99)); } }
    public decimal? DisorderCount { get => _disCount; set { if (value != null) Set(ref _disCount, Math.Clamp(value.Value, 1, 12)); } }
    public decimal? DisorderPick { get => _disPick; set { if (value != null) Set(ref _disPick, Math.Max(1, value.Value)); } }
    public bool DisorderRunning { get => _disRunning; private set { if (Set(ref _disRunning, value)) RaiseBusy(); } }
    public string DisorderError { get => _disError; private set { if (Set(ref _disError, value)) Raise(nameof(DisorderHasError)); } }
    public bool DisorderHasError => _disError.Length > 0;
    public string DisorderText { get => _disText; private set { if (Set(ref _disText, value)) Raise(nameof(DisorderHasText)); } }
    public bool DisorderHasText => _disText.Length > 0;
    public bool DisorderHasFiles => _disFiles.Length > 0;

    private async Task<JsonNode?> CallDisorder(JsonObject o)
    {
        if (_doc == null || !Idle) return null;
        if (_disFrom.Length == 0 || _disTo.Length == 0) { DisorderError = "Name the element whose sites change and the one that takes them (Cu, Au)."; return null; }
        o["from"] = _disFrom;
        o["to"] = _disTo;
        var doc = _doc!;
        var json = o.ToJsonString();
        DisorderRunning = true;
        DisorderError = "";
        try
        {
            var r = JsonNode.Parse(await Task.Run(() => doc.Disorder(json)))!;
            if (r["ok"]?.GetValue<bool>() != true) { DisorderError = (string?)r["error"] ?? "disorder failed"; return null; }
            return r;
        }
        catch (Exception e) { DisorderError = e.Message; return null; }
        finally { DisorderRunning = false; }
    }

    /// <summary>A special quasirandom structure on the open cell, opened as a new structure.</summary>
    public async Task RunSqs()
    {
        var path = Path.Combine(Path.GetTempPath(), $"caps-sqs-{Environment.ProcessId}-{DateTime.Now:HHmmss}.cif");
        var r = await CallDisorder(new JsonObject { ["op"] = "sqs", ["fraction"] = (double)_disPercent / 100, ["path"] = path });
        if (r == null) return;
        var inv = CultureInfo.InvariantCulture;
        var sb = new StringBuilder();
        sb.AppendLine(string.Format(inv, "{0} of {1} sites · random alloy (1 − 2x)² = {2:0.0000}", (int)((double?)r["substituted"] ?? 0), (int)((double?)r["sites"] ?? 0), (double?)r["target"] ?? 0));
        sb.AppendLine($"{"shell",-7}{"r (Å)",9}{"Π",10}");
        var d = (JsonArray)r["shell_distance"]!;
        var c = (JsonArray)r["correlation"]!;
        for (var m = 0; m < c.Count; ++m) sb.AppendLine(string.Format(inv, "{0,-7}{1,9:0.000}{2,10:0.0000}", m + 1, (double?)d[m] ?? 0, (double?)c[m] ?? 0));
        DisorderText = sb.ToString().TrimEnd();
        _disFiles = [];
        Raise(nameof(DisorderHasFiles));
        Status = "Quasirandom structure opened";
        await OpenProgressive(path, null);
    }

    /// <summary>The symmetry-distinct arrangements of the substituents, listed with their counts; each one can be opened.</summary>
    public async Task RunEnumerate()
    {
        var dir = Path.Combine(Path.GetTempPath(), $"caps-arrangements-{Environment.ProcessId}-{DateTime.Now:HHmmss}");
        var r = await CallDisorder(new JsonObject { ["op"] = "enumerate", ["count"] = (int)_disCount, ["directory"] = dir, ["write"] = 500 });
        if (r == null) return;
        var conf = (JsonArray)r["configurations"]!;
        var mult = (JsonArray)r["multiplicity"]!;
        var sites = (JsonArray)r["sites"]!;
        var sb = new StringBuilder();
        sb.AppendLine($"{(long)((double?)r["total"] ?? 0)} arrangements, {conf.Count} distinct ({(int)((double?)r["operations"] ?? 0)} symmetry operations)");
        sb.AppendLine($"{"#",-5}{"count",8}  atoms");
        for (var k = 0; k < conf.Count && k < 40; ++k)
            sb.AppendLine($"{k + 1,-5}{(long)((double?)mult[k] ?? 0),8}  " + string.Join(" ", ((JsonArray)conf[k]!).Select(x => ((int?)sites[(int?)x ?? 0] ?? 0) + 1)));
        if (conf.Count > 40) sb.AppendLine($"… {conf.Count - 40} more");
        DisorderText = sb.ToString().TrimEnd();
        _disFiles = r["files"] is JsonArray f ? f.Select(x => (string?)x ?? "").ToArray() : [];
        _disPick = 1;
        Raise(nameof(DisorderPick));
        Raise(nameof(DisorderHasFiles));
        Status = $"{conf.Count} distinct arrangements";
    }

    /// <summary>Opens the arrangement picked from the list as a new structure.</summary>
    public async Task OpenArrangement()
    {
        var k = (int)_disPick - 1;
        if (k < 0 || k >= _disFiles.Length) { DisorderError = $"Pick 1 to {_disFiles.Length}."; return; }
        await OpenProgressive(_disFiles[k], null);
    }
}
