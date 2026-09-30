using System;
using System.IO;
using System.Text.Json.Nodes;
using System.Threading;
using System.Threading.Tasks;
using CapsStudio.Interop;

namespace CapsStudio.ViewModels;

/// <summary>Start › From a recipe: a CAPS recipe (YAML or JSON, as `caps run` reads it) run off the UI thread — every
/// stage in the status bar — and the structure it made opened as the document. Files the recipe exports go beside it.</summary>
public partial class MainViewModel
{
    private bool _recipeRunning;
    private CancellationTokenSource? _recipeCancel;
    public bool RecipeRunning { get => _recipeRunning; private set { Set(ref _recipeRunning, value); Raise(nameof(Busy)); Raise(nameof(HasSelBar)); } }
    public string RecipeLog { get; private set; } = "";

    public void CancelRecipe() => _recipeCancel?.Cancel();

    public async Task RunRecipeFile(string path)
    {
        string text;
        try { text = File.ReadAllText(path); }
        catch (Exception e) { Status = "Could not read the recipe: " + e.Message; return; }
        await RunRecipeText(text, Path.GetFileNameWithoutExtension(path), Path.GetDirectoryName(Path.GetFullPath(path)) ?? ".");
    }

    /// <summary>Runs recipe text (YAML or JSON); relative paths and exported files are in `dir`.</summary>
    public async Task RunRecipeText(string text, string label, string dir)
    {
        if (Busy || _recipeRunning) { Status = "Wait for the run to finish (or cancel it) before running a recipe"; return; }
        var options = new JsonObject { ["base_dir"] = dir, ["out_dir"] = dir, ["forcefield_dir"] = Paths.ForceFields ?? "", ["threads"] = _settings.Threads }.ToJsonString();
        RecipeRunning = true;
        _recipeCancel = new CancellationTokenSource();
        var token = _recipeCancel.Token;
        var log = new System.Text.StringBuilder();
        Status = $"Recipe {label} · starting";
        var (doc, report) = await Task.Run(() => CapsDocument.RunRecipe(text, options, label, (k, n, name, st, detail, f) =>
        {
            var line = $"[{k}/{n}] {name} · {(detail.Length > 0 ? detail + " · " : "")}{st}";
            if (st != "running") lock (log) log.AppendLine(line);
            Avalonia.Threading.Dispatcher.UIThread.Post(() => Status = $"Recipe {label} " + line + (st == "running" ? $" {f * 100:0} %" : ""));
            return !token.IsCancellationRequested;
        }));
        RecipeRunning = false;
        RecipeLog = log.ToString();
        var rep = JsonNode.Parse(report.Length > 0 ? report : "{}");
        if (doc == null)
        {
            var code = rep?["exit"]?.GetValue<double>() ?? 4;
            var err = rep?["error"]?.GetValue<string>() ?? "failed";
            Status = token.IsCancellationRequested ? $"Recipe {label} cancelled" : $"Recipe {label} failed (exit {code:0}): {err}";
            return;
        }
        var files = rep?["files"] is JsonArray fa ? fa.Count : 0;
        Show(doc, label);
        GrownUnsaved = files == 0;
        SetModule(8);
        Status = $"Recipe {label} done" + (files > 0 ? $" · wrote {files} file{(files == 1 ? "" : "s")} beside it, with provenance" : " · nothing exported: save the cell to keep it");
    }
}
