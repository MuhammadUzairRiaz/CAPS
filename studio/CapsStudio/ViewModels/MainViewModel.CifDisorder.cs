using System.Globalization;
using System.Runtime.CompilerServices;
using System.Text.Json.Nodes;
using CapsStudio.Interop;

namespace CapsStudio.ViewModels;

/// <summary>Builders › Crystal › CIF disorder (core crystal.hpp CifDisorder): the open structure's CIF read again with its
/// disorder resolved: each assembly's major group, a named group, or a supercell filled at random to the occupancies.</summary>
public sealed partial class MainViewModel
{
    // the file each opened document came from
    private readonly ConditionalWeakTable<CapsDocument, string> _docPaths = new();
    public static string[] CifModes { get; } = ["Major group and occupant", "A named group", "Random fill of a supercell"];
    private int _cifMode;
    private string _cifGroup = "2", _cifInfo = "", _cifError = "", _cifPath = "";
    private decimal _cifNa = 3, _cifNb = 3, _cifNc = 3, _cifSeed = 1;
    public int CifMode { get => _cifMode; set { if (Set(ref _cifMode, Math.Clamp(value, 0, 2))) { Raise(nameof(CifNamed)); Raise(nameof(CifRandom)); } } }
    public bool CifNamed => _cifMode == 1;
    public bool CifRandom => _cifMode == 2;
    public string CifGroup { get => _cifGroup; set => Set(ref _cifGroup, value.Trim()); }
    public decimal CifNa { get => _cifNa; set => Set(ref _cifNa, Math.Clamp(value, 1, 20)); }
    public decimal CifNb { get => _cifNb; set => Set(ref _cifNb, Math.Clamp(value, 1, 20)); }
    public decimal CifNc { get => _cifNc; set => Set(ref _cifNc, Math.Clamp(value, 1, 20)); }
    public decimal CifSeed { get => _cifSeed; set => Set(ref _cifSeed, Math.Clamp(decimal.Truncate(value), 1, int.MaxValue)); }
    public string CifInfo { get => _cifInfo; private set { if (Set(ref _cifInfo, value)) Raise(nameof(CifHasDisorder)); } }
    public bool CifHasDisorder => _cifInfo.Length > 0;
    public string CifError { get => _cifError; private set { if (Set(ref _cifError, value)) Raise(nameof(CifHasError)); } }
    public bool CifHasError => _cifError.Length > 0;

    /// <summary>What the open structure's CIF holds of disorder (nothing shown when it has none).</summary>
    public void RefreshCifDisorder()
    {
        static bool IsCif(string? f) => f != null && f.EndsWith(".cif", StringComparison.OrdinalIgnoreCase) && File.Exists(f);
        _cifPath = _doc == null ? "" : IsCif(_doc.Path) ? _doc.Path : _docPaths.TryGetValue(_doc, out var p) && IsCif(p) ? p : "";
        if (_cifPath.Length == 0) { CifInfo = ""; return; }
        try
        {
            var r = JsonNode.Parse(CapsDocument.CifDisorder(_cifPath))!;
            if ((bool?)r["ok"] != true) { CifInfo = ""; return; }
            int a = (int)((double?)r["assemblies"] ?? 0), g = (int)((double?)r["groups"] ?? 0), pl = (int)((double?)r["places"] ?? 0);
            if (a == 0 && pl == 0) { CifInfo = ""; return; }
            var inv = CultureInfo.InvariantCulture;
            var ex = r["expected"] as JsonObject;
            var comp = ex == null ? "" : " · the occupancies give " + string.Join(" ", ex.Select(kv => kv.Key + ((double?)kv.Value ?? 0).ToString("0.##", inv)));
            CifInfo = $"{Path.GetFileName(_cifPath)}: {a} disorder assembl{(a == 1 ? "y" : "ies")} ({g} groups), {pl} sites with alternatives or vacancies per cell{comp}.";
        }
        catch { CifInfo = ""; }
    }

    /// <summary>The CIF read again with the chosen treatment, as a new structure.</summary>
    public void ReadCifAgain()
    {
        if (_cifPath.Length == 0) return;
        var mode = _cifMode == 0 ? "major" : _cifMode == 1 ? "group" : "random";
        var o = new JsonObject { ["mode"] = mode, ["group"] = _cifGroup, ["supercell"] = new JsonArray((int)_cifNa, (int)_cifNb, (int)_cifNc), ["seed"] = (int)_cifSeed };
        try
        {
            var label = Path.GetFileNameWithoutExtension(_cifPath) + (mode == "random" ? $" (random {(int)_cifNa}×{(int)_cifNb}×{(int)_cifNc})" : mode == "group" ? $" (group {_cifGroup})" : " (major)");
            var (doc, rep) = CapsDocument.OpenCif(_cifPath, o.ToJsonString(), label);
            var path = _cifPath;
            CifError = "";
            Show(doc, label);
            _docPaths.AddOrUpdate(doc, path);
            var r = JsonNode.Parse(rep);
            var inv = CultureInfo.InvariantCulture;
            if (r?["obtained"] is JsonObject ob && r["expected"] is JsonObject ex)
                Status = label + ": " + string.Join(" ", ob.Select(kv => $"{kv.Key} {((double?)kv.Value ?? 0):0} (expected {((double?)ex[kv.Key] ?? 0).ToString("0.#", inv)})"));
        }
        catch (Exception e) { CifError = e.Message; }
    }
}
