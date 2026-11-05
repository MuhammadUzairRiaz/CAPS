using System;
using System.Globalization;
using System.Linq;
using System.Text.Json.Nodes;
using System.Threading.Tasks;
using CapsStudio.Interop;

namespace CapsStudio.ViewModels;

// Polymer page › Dendrimers from units: AB_n units with a core, end and focal group (core dendrimer.hpp through
// caps_dendrimer): perfect shells, slow monomer addition or Flory's statistics; the molecule opens as a new document.
public partial class MainViewModel
{
    public static readonly (string Id, string Name)[] DendPresets =
        [("pamam", "PAMAM (EDA core)"), ("ppi", "PPI (DAB core)"), ("bismpa", "bis-MPA (TMP core)"), ("", "Own units")];
    public string[] DendPresetNames => DendPresets.Select(p => p.Name).ToArray();
    private int _dendPreset;
    public int DendPreset { get => _dendPreset; set { if (Set(ref _dendPreset, value)) { Raise(nameof(DendOwnUnits)); Raise(nameof(DendGenerationLabel)); } } }
    public bool DendOwnUnits => DendPresets[_dendPreset].Id.Length == 0;
    public string DendGenerationLabel => DendOwnUnits ? "Shells" : "Generation";

    public static readonly (string Id, string Name)[] DendGrowths =
        [("perfect", "Perfect (every site)"), ("sma", "Slow monomer addition"), ("flory", "Flory statistics (no core)")];
    public string[] DendGrowthNames => DendGrowths.Select(p => p.Name).ToArray();
    private int _dendGrowth;
    public int DendGrowth { get => _dendGrowth; set { if (Set(ref _dendGrowth, value)) Raise(nameof(DendHyper)); } }
    public bool DendHyper => _dendGrowth > 0;

    private decimal _dendGeneration = 2;
    public decimal DendGeneration { get => _dendGeneration; set => Set(ref _dendGeneration, value); }
    private decimal _dendUnits = 30;
    public decimal DendUnits { get => _dendUnits; set => Set(ref _dendUnits, value); }
    private decimal _dendSeed = 1;
    public decimal DendSeed { get => _dendSeed; set => Set(ref _dendSeed, value); }
    private string _dendCore = "CCC(CO[*:2])(CO[*:2])CO[*:2]";
    public string DendCore { get => _dendCore; set => Set(ref _dendCore, value ?? ""); }
    private string _dendUnit = "[*:1]C(=O)C(C)(CO[*:2])CO[*:2]";
    public string DendUnit { get => _dendUnit; set => Set(ref _dendUnit, value ?? ""); }
    private string _dendEnd = "";
    public string DendEnd { get => _dendEnd; set => Set(ref _dendEnd, value ?? ""); }
    private string _dendFocal = "";
    public string DendFocal { get => _dendFocal; set => Set(ref _dendFocal, value ?? ""); }
    private bool _dendRelax = true;
    public bool DendRelax { get => _dendRelax; set => Set(ref _dendRelax, value); }

    private string _dendResult = "";
    public string DendResult { get => _dendResult; private set { if (Set(ref _dendResult, value)) Raise(nameof(DendHasResult)); } }
    public bool DendHasResult => _dendResult.Length > 0;
    private string _dendError = "";
    public string DendError { get => _dendError; private set { if (Set(ref _dendError, value)) Raise(nameof(DendHasError)); } }
    public bool DendHasError => _dendError.Length > 0;
    private bool _dendRunning;
    public bool DendRunning { get => _dendRunning; private set => Set(ref _dendRunning, value); }

    public async Task BuildDendrimer()
    {
        if (!Idle || _dendRunning) return;
        var preset = DendPresets[_dendPreset].Id;
        var growth = DendGrowths[_dendGrowth].Id;
        var j = new JsonObject
        {
            ["preset"] = preset, ["growth"] = growth, ["units"] = (int)_dendUnits, ["seed"] = (int)_dendSeed, ["relax"] = _dendRelax,
            ["end"] = _dendEnd.Trim(),
        };
        if (preset.Length > 0) j["generation"] = (int)_dendGeneration;
        else { j["core"] = _dendCore.Trim(); j["unit"] = _dendUnit.Trim(); j["shells"] = (int)_dendGeneration; j["focal"] = _dendFocal.Trim(); }
        var json = j.ToJsonString();
        DendRunning = true;
        DendError = "";
        Status = "Building the dendrimer…";
        try
        {
            var title = preset.Length > 0 ? $"{preset.ToUpperInvariant()} G{(int)_dendGeneration}" : $"dendrimer ({growth})";
            if (growth != "perfect") title = (preset.Length > 0 ? preset : "AB") + $" hyperbranched {(int)_dendUnits} ({growth})";
            var (nd, rep) = await Task.Run(() => CapsDocument.Dendrimer(json, title));
            var r = JsonNode.Parse(rep)!;
            if (nd == null || r["ok"]?.GetValue<bool>() != true) { DendError = (string?)r["error"] ?? "the build failed"; Status = "Dendrimer: " + DendError; return; }
            var inv = CultureInfo.InvariantCulture;
            double D(string k) => r[k]?.GetValue<double>() ?? 0;
            var shells = string.Join(" · ", ((JsonArray)r["per_shell"]!).Select(x => ((int)x!.GetValue<double>()).ToString(inv)));
            var text = $"{(int)D("units")} AB{(int)D("branches")} units in {(int)D("shells")} shells ({shells}) · {(string?)r["formula"]} · {D("mass").ToString("0.0", inv)} g/mol · {(int)D("atoms")} atoms";
            if (r["db_frey"] is JsonNode)
                text += $" · D {(int)D("dendritic")}, L {(int)D("linear")}, T {(int)D("terminal")} · DB {D("db_frey").ToString("0.000", inv)} (Frey), {D("db_frechet").ToString("0.000", inv)} (Fréchet); large-N limit {D("db_limit").ToString("0.000", inv)}";
            if ((int)D("ends") > 0) text += $" · {(int)D("ends")} end groups";
            DendResult = text;
            DendRunning = false;
            Show(nd, title);
            Status = "Dendrimer: " + text;
        }
        catch (Exception e) { DendError = e.Message; }
        finally { DendRunning = false; }
    }
}
