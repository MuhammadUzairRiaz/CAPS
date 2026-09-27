using System;
using System.Collections.ObjectModel;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Text.Json.Nodes;
using System.Threading.Tasks;
using CapsStudio.Interop;

namespace CapsStudio.ViewModels;

public sealed record RecipeItem(string Name, string Detail, string? Path, string Text)
{
    public bool IsTemplate => Path == null;
}
public sealed record RecipeStageChip(string Name, string Summary, bool Ok);
public sealed record ScheduleRow(string Number, string Ensemble, string T, string P, string Ps);

/// <summary>Jobs › Recipes (design/boards/RecipeEditor): reusable chains of steps as plain YAML — the same files `caps run`
/// reads — checked as they are edited, the equilibration schedule drawn from the protocol, saved with their SHA-256 (which
/// every run of them records), duplicated, and run on their own or on the open structure.</summary>
public partial class MainViewModel
{
    public bool IsRecipes => _module == 53;
    public static string RecipesFolder => Environment.GetEnvironmentVariable("CAPS_RECIPES") is { Length: > 0 } d ? d
        : Path.Combine(Environment.GetFolderPath(Environment.SpecialFolder.UserProfile), "CAPS", "recipes");
    public ObservableCollection<RecipeItem> Recipes { get; } = new();
    public ObservableCollection<RecipeStageChip> RecipeStages { get; } = new();
    public ObservableCollection<ScheduleRow> RecipeSchedule { get; } = new();
    private RecipeItem? _recipe;
    private string _recipeText = "", _recipeName = "", _recipeValid = "", _recipeSha = "", _recipeProtocol = "";
    private bool _recipeOk, _recipeDirty;
    public RecipeItem? SelectedRecipe
    {
        get => _recipe;
        set { if (Set(ref _recipe, value) && value != null) { _recipeText = value.Text; Raise(nameof(RecipeText)); _recipeDirty = false; CheckRecipeText(); } }
    }
    public string RecipeText { get => _recipeText; set { if (Set(ref _recipeText, value ?? "")) { _recipeDirty = true; CheckRecipeText(); } } }
    public string RecipeName { get => _recipeName; private set => Set(ref _recipeName, value); }
    public string RecipeValid { get => _recipeValid; private set => Set(ref _recipeValid, value); }
    public bool RecipeOk { get => _recipeOk; private set => Set(ref _recipeOk, value); }
    public string RecipeSha { get => _recipeSha; private set => Set(ref _recipeSha, value); }
    public string RecipeProtocol { get => _recipeProtocol; private set { Set(ref _recipeProtocol, value); Raise(nameof(RecipeHasSchedule)); } }
    public bool RecipeHasSchedule => _recipeProtocol.Length > 0;
    public string RecipePathText => _recipe?.Path is { } p ? RecentFiles.Tilde(p) + (_recipeDirty ? " · edited" : "") : "template · Duplicate to keep your own";
    public (double X, double Y)[] RecipeT { get; private set; } = [];
    public (double X, double Y)[] RecipeP { get; private set; } = [];
    public event Action? RecipeChanged;

    private static readonly (string Name, string Text)[] RecipeTemplates =
    [
        ("Amorphous polymer · standard", "# A polystyrene melt: grown, relaxed, equilibrated with the 21-step protocol, then run NPT\nrecipe: 1\nname: amorphous_standard\nbuild:\n  polymer: { smiles: \"*CC(*)c1ccccc1\", dp: 20, chains: 10, tacticity: atactic }\ntype: { forcefield: gaff2 }\ngrow: { density: 0.5, seed: 1 }\nrelax: { method: lbfgs, fmax: 0.5 }\nequilibrate: { protocol: larsen21, t_max: 600, t_final: 300, time_scale: 0.1 }\nmd: { ps: 100, temperature: 300, ensemble: npt }\nanalyze: { properties: [density, rg] }\nexport: [lammps, gromacs]\n"),
        ("Natural rubber melt", "# cis-1,4-polyisoprene chains for fibre–rubber work\nrecipe: 1\nname: nr_melt\nbuild:\n  polymer: { smiles: \"*C/C=C(/C)C*\", dp: 30, chains: 12 }\ntype: { forcefield: default }\ngrow: { density: 0.6, seed: 1 }\nrelax: { method: lbfgs, fmax: 1.0 }\nmd: { ps: 50, temperature: 300, ensemble: npt }\nanalyze: { properties: [density, rg, cn] }\nexport: [lammps]\n"),
        ("Sulfur-cured natural rubber · PCFF", "# NR chains, H–S–S–H donors packed in, cured (C–S–S–C crosslinks), the network typed with PCFF, relaxed, exported to LAMMPS class2\nrecipe: 1\nname: nr_sulfur_cure\nbuild:\n  polymer: { smiles: \"[*]C/C=C(C)\\\\C[*]\", dp: 30, chains: 8 }\ntype: { forcefield: pcff-frc }\ngrow: { density: 0.6, seed: 1 }\nreact:\n  insert: { smiles: SS, count: 24 }\n  templates: [sulfur_allylic]\n  cycles: 60\n  per_cycle: 5\n  md_ps: 1\n  temperature: 450\n  seed: 2\nrelax: { method: lbfgs, fmax: 1.0 }\nexport: [lammps]\n"),
        ("Glass transition", "# A stepwise NPT cooling scan of a saved cell (Run on the open structure uses it instead of the file)\nrecipe: 1\nname: glass_transition\nbuild: { file: cell.data }\ntype: { forcefield: default }\nanalyze:\n  properties: [tg]\n  tg: { t_start: 500, t_end: 250, t_step: 25, ps_per_step: 100, equilibrate_ps: 100, seed: 1 }\n"),
        ("Quick molecule", "recipe: 1\nname: ethanol\nbuild: { molecule: CCO }\ntype: { forcefield: uff, charges: qeq }\nrelax: { fmax: 0.5 }\nexport: [xyz, mol2]\n"),
    ];

    public void OpenRecipes()
    {
        LoadRecipes();
        SetModule(53);
        SelectedRecipe ??= Recipes.FirstOrDefault();
        RecipeChanged?.Invoke();
    }

    public void LoadRecipes(string? select = null)
    {
        Recipes.Clear();
        try
        {
            if (Directory.Exists(RecipesFolder))
                foreach (var f in Directory.GetFiles(RecipesFolder, "*.y*ml").OrderBy(f => f))
                {
                    var text = File.ReadAllText(f);
                    Recipes.Add(new RecipeItem(Path.GetFileNameWithoutExtension(f), StepsText(text) + " · yours", f, text));
                }
        }
        catch { /* an unreadable folder: templates only */ }
        foreach (var (name, text) in RecipeTemplates) Recipes.Add(new RecipeItem(name, StepsText(text) + " · template", null, text));
        if (select != null) SelectedRecipe = Recipes.FirstOrDefault(r => r.Path == select) ?? SelectedRecipe;
    }

    private static string StepsText(string text)
    {
        try
        {
            var n = ((JsonArray?)JsonNode.Parse(CapsDocument.RecipeCheck(text))?["stages"])?.Count ?? 0;
            return $"{n} step{(n == 1 ? "" : "s")}";
        }
        catch { return "—"; }
    }

    private void CheckRecipeText()
    {
        var inv = CultureInfo.InvariantCulture;
        RecipeStages.Clear();
        RecipeSchedule.Clear();
        JsonNode r;
        try { r = JsonNode.Parse(CapsDocument.RecipeCheck(_recipeText))!; }
        catch (Exception e) { RecipeOk = false; RecipeValid = e.Message; return; }
        RecipeOk = r["ok"]?.GetValue<bool>() == true;
        RecipeValid = RecipeOk ? "Pipeline valid" : (string?)r["error"] ?? "invalid";
        RecipeName = (string?)r["name"] ?? "recipe";
        RecipeSha = (string?)r["sha256"] ?? "";
        foreach (var st in (JsonArray)r["stages"]!)
            RecipeStages.Add(new RecipeStageChip(Cap((string)st!["name"]!), (string?)st["summary"] ?? "", st["ok"]?.GetValue<bool>() == true));
        RecipeProtocol = (string?)r["protocol"] ?? "";
        // the schedule: T (and P for NPT steps) against time, and the table
        var t = new System.Collections.Generic.List<(double, double)>();
        var p = new System.Collections.Generic.List<(double, double)>();
        double time = 0;
        var k = 0;
        foreach (var g in (JsonArray?)r["schedule"] ?? new JsonArray())
        {
            var ps = g!["ps"]!.GetValue<double>();
            var t0 = g["t_start"]!.GetValue<double>();
            var t1 = g["t_end"]!.GetValue<double>();
            var bar = g["pressure_bar"]!.GetValue<double>();
            var ens = (string)g["ensemble"]!;
            t.Add((time, t0));
            t.Add((time + ps, t1));
            if (ens == "NPT" && bar > 0) { p.Add((time, Math.Log10(bar))); p.Add((time + ps, Math.Log10(bar))); }
            time += ps;
            RecipeSchedule.Add(new ScheduleRow((++k).ToString(inv), ens, Math.Abs(t0 - t1) < 1e-9 ? t0.ToString("0", inv) : $"{t0:0}→{t1:0}",
                                               ens == "NPT" ? bar.ToString("N0", inv) : "—", ps.ToString("0.#", inv)));
        }
        RecipeT = t.ToArray();
        RecipeP = p.ToArray();
        Raise(nameof(RecipePathText));
        RecipeChanged?.Invoke();
    }

    private static string Cap(string s) => s.Length == 0 ? s : char.ToUpperInvariant(s[0]) + s[1..];

    /// <summary>Writes the recipe (a template is saved as a new file in the recipes folder).</summary>
    public void SaveRecipe()
    {
        try
        {
            Directory.CreateDirectory(RecipesFolder);
            var path = _recipe?.Path ?? UniquePath(RecipeName);
            File.WriteAllText(path, _recipeText);
            _recipeDirty = false;
            LoadRecipes(path);
            Status = $"Saved {Path.GetFileName(path)} · sha256 {RecipeSha[..12]}… — every run of it records this hash";
        }
        catch (Exception e) { Status = "Could not save the recipe: " + e.Message; }
    }

    public void DuplicateRecipe()
    {
        try
        {
            Directory.CreateDirectory(RecipesFolder);
            var path = UniquePath(RecipeName + "_copy");
            File.WriteAllText(path, _recipeText);
            LoadRecipes(path);
            Status = $"Duplicated as {Path.GetFileName(path)}";
        }
        catch (Exception e) { Status = "Could not duplicate: " + e.Message; }
    }

    private static string UniquePath(string name)
    {
        var stem = string.Concat(name.Select(c => char.IsLetterOrDigit(c) || c is '_' or '-' ? c : '_'));
        var path = Path.Combine(RecipesFolder, stem + ".yaml");
        for (var k = 2; File.Exists(path); ++k) path = Path.Combine(RecipesFolder, $"{stem}_{k}.yaml");
        return path;
    }

    public async Task RunSelectedRecipe()
    {
        if (!RecipeOk) { Status = "Fix the recipe first: " + RecipeValid; return; }
        var dir = _recipe?.Path is { } p ? Path.GetDirectoryName(p)! : RecipesFolder;
        Directory.CreateDirectory(dir);
        await RunRecipeText(_recipeText, RecipeName, dir);
    }

    /// <summary>The recipe with its build stage replaced by the open structure (saved to a temporary file).</summary>
    public async Task RunRecipeOnDocument()
    {
        if (_doc == null || !RecipeOk) return;
        var j = JsonNode.Parse(CapsDocument.YamlToJson(_recipeText)) as JsonObject;
        if (j == null) return;
        var tmp = Path.Combine(Path.GetTempPath(), "caps-recipe-input.data");
        _doc.Save(tmp);
        j["build"] = new JsonObject { ["file"] = tmp };
        j.Remove("grow");   // the structure is already made
        var dir = _recipe?.Path is { } p ? Path.GetDirectoryName(p)! : RecipesFolder;
        Directory.CreateDirectory(dir);
        await RunRecipeText(j.ToJsonString(), RecipeName + "_on_" + Path.GetFileNameWithoutExtension(Title), dir);
    }
}
