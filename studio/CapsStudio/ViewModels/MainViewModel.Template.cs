using System.Collections.ObjectModel;
using System.Text.Json.Nodes;
using CapsStudio.Interop;

namespace CapsStudio.ViewModels;

/// <summary>A check of a reaction template.</summary>
public sealed record TemplateCheck(bool Ok, string Text);

/// <summary>Reaction template editor (design/boards/ReactionTemplate): the atom-mapped template as text, drawn before
/// and after the reaction with its map numbers, the bond changes, checks, a test on the open structure (reactive sites
/// and matches within the capture distance), saved to ~/.caps/templates and handed to React.</summary>
public sealed partial class MainViewModel
{
    public bool IsTemplate => _module == 45;
    public static string TemplateFolder => System.IO.Path.Combine(AppSettings.Folder, "templates");
    public ObservableCollection<string> TemplateNames { get; } = new();
    public ObservableCollection<ProvRow> TemplateChanges { get; } = new();
    public ObservableCollection<TemplateCheck> TemplateChecks { get; } = new();
    public ObservableCollection<string> TemplateTabs { get; } = new();
    private string _tplName = "", _tplText = "", _tplTitle = "", _tplError = "", _tplTest = "", _tplInfo = "";
    private int _tplIndex;
    private JsonArray? _tplViews;
    public event Action? TemplateDrawn;
    public JsonNode? TemplatePre { get; private set; }
    public JsonNode? TemplatePost { get; private set; }
    public List<(int, int)> TemplateFormed { get; } = new();
    public List<(int, int)> TemplateBroken { get; } = new();

    public string TemplateName { get => _tplName; set { if (value != null && Set(ref _tplName, value)) LoadTemplate(value); } }
    public string TemplateText { get => _tplText; set { if (Set(ref _tplText, value ?? "")) ParseTemplate(); } }
    public int TemplateIndex { get => _tplIndex; set { if (Set(ref _tplIndex, value)) ShowTemplate(); } }
    public string TemplateTitle { get => _tplTitle; private set => Set(ref _tplTitle, value); }
    public string TemplateError { get => _tplError; private set { if (Set(ref _tplError, value)) Raise(nameof(TemplateHasError)); } }
    public bool TemplateHasError => _tplError.Length > 0;
    public string TemplateTestText { get => _tplTest; private set => Set(ref _tplTest, value); }
    public string TemplateInfo { get => _tplInfo; private set => Set(ref _tplInfo, value); }
    public bool TemplateChecksPass => TemplateChecks.Count > 0 && TemplateChecks.All(c => c.Ok);
    public string TemplateChecksChip => TemplateChecksPass ? "pass" : $"{TemplateChecks.Count(c => !c.Ok)} to fix";
    public string TemplateFooter => $"Template {_tplIndex + 1} of {Math.Max(1, TemplateTabs.Count)}" + (TemplateTabs.Count > 1 ? " · " + string.Join(", ", TemplateTabs) : "");

    public void OpenTemplateEditor(string? name = null)
    {
        TemplateNames.Clear();
        try { foreach (var n in CapsDocument.ReactionTemplate("").Split('\n', StringSplitOptions.RemoveEmptyEntries | StringSplitOptions.TrimEntries)) TemplateNames.Add(n); } catch { }
        if (Directory.Exists(TemplateFolder))
            foreach (var f in Directory.EnumerateFiles(TemplateFolder, "*.txt").OrderBy(f => f)) TemplateNames.Add("my: " + System.IO.Path.GetFileNameWithoutExtension(f));
        SetModule(45);
        TemplateName = name ?? (_tplName.Length > 0 ? _tplName : TemplateNames.FirstOrDefault(n => n.StartsWith("epoxy_amine_primary", StringComparison.Ordinal)) ?? TemplateNames.FirstOrDefault() ?? "");
        if (_tplText.Length == 0 && _tplName.Length > 0) LoadTemplate(_tplName);
    }

    private void LoadTemplate(string name)
    {
        try
        {
            TemplateText = name.StartsWith("my: ", StringComparison.Ordinal)
                ? File.ReadAllText(System.IO.Path.Combine(TemplateFolder, name[4..] + ".txt"))
                : CapsDocument.ReactionTemplate(name);
            TemplateError = "";
        }
        catch (Exception e) { TemplateError = e.Message; }
    }

    private void ParseTemplate()
    {
        TemplateTabs.Clear();
        _tplViews = null;
        try
        {
            var j = JsonNode.Parse(CapsDocument.TemplateView(_tplText))!;
            if (j["ok"]?.GetValue<bool>() != true) { TemplateError = j["error"]?.GetValue<string>() ?? "cannot read the template"; ShowTemplate(); return; }
            TemplateError = "";
            _tplViews = (JsonArray)j["templates"]!;
            foreach (var v in _tplViews) TemplateTabs.Add(v!["name"]!.GetValue<string>());
        }
        catch (Exception e) { TemplateError = e.Message; }
        if (_tplIndex >= TemplateTabs.Count) _tplIndex = 0;
        Raise(nameof(TemplateIndex));
        ShowTemplate();
    }

    private void ShowTemplate()
    {
        TemplateChanges.Clear();
        TemplateChecks.Clear();
        TemplateFormed.Clear();
        TemplateBroken.Clear();
        var v = _tplViews != null && _tplIndex < _tplViews.Count ? _tplViews[_tplIndex] : null;
        TemplatePre = v?["pre"];
        TemplatePost = v?["post"];
        if (v != null)
        {
            TemplateTitle = "Reaction template · " + v["name"]!.GetValue<string>();
            foreach (var c in (JsonArray)v["changes"]!)
            {
                var kind = c!["kind"]!.GetValue<string>();
                TemplateChanges.Add(new ProvRow(kind == "formed" ? "bond formed" : kind == "broken" ? "bond broken" : kind == "moved" ? "atom moved" : "atom deleted",
                                                c["text"]!.GetValue<string>(), kind is "formed" ? "new" : kind is "broken" ? "breaks" : kind));
            }
            foreach (var c in (JsonArray)v["checks"]!) TemplateChecks.Add(new TemplateCheck(c!["ok"]!.GetValue<bool>(), c["text"]!.GetValue<string>()));
            HashSet<(int, int)> Bonds(JsonNode? d) => d?["bonds"] is JsonArray a ? a.Select(x => (Math.Min(x![0]!.GetValue<int>(), x[1]!.GetValue<int>()), Math.Max(x[0]!.GetValue<int>(), x[1]!.GetValue<int>()))).ToHashSet() : new();
            var pre = Bonds(TemplatePre);
            var post = Bonds(TemplatePost);
            TemplateFormed.AddRange(post.Except(pre));
            TemplateBroken.AddRange(pre.Except(post));
            var init = (JsonArray)v["initiators"]!;
            TemplateInfo = $"initiators {init[0]} and {init[1]} · capture {v["capture"]!.GetValue<double>():0.0#} Å · probability {v["probability"]!.GetValue<double>():0.##} · " +
                           (v["min_path"]!.GetValue<int>() is var mp && mp > 0 ? $"at least {mp} bonds apart or in different molecules" : "initiators in different molecules");
        }
        foreach (var n in new[] { nameof(TemplateChecksPass), nameof(TemplateChecksChip), nameof(TemplateFooter) }) Raise(n);
        TemplateDrawn?.Invoke();
    }

    /// <summary>Sites and matches of every template of the text on the open structure's current frame.</summary>
    public void TestTemplate()
    {
        if (_doc == null) { TemplateTestText = "Open a structure to test the template on"; return; }
        try
        {
            var j = JsonNode.Parse(_doc.TemplateTest(_tplText))!;
            if (j["ok"]?.GetValue<bool>() != true) { TemplateTestText = j["error"]?.GetValue<string>() ?? "failed"; return; }
            TemplateTestText = string.Join("\n", ((JsonArray)j["templates"]!).Select(t =>
                $"{t!["name"]!.GetValue<string>()}: {t["sites"]!.GetValue<int>()} reactive sites · {t["matches"]!.GetValue<double>():0} within capture" +
                (t["closest"]!.GetValue<double>() >= 0 ? $" (closest {t["closest"]!.GetValue<double>():0.00} Å)" : "")));
        }
        catch (Exception e) { TemplateTestText = e.Message; }
    }

    public string SaveTemplate()
    {
        Directory.CreateDirectory(TemplateFolder);
        var name = TemplateTabs.FirstOrDefault() ?? "template";
        var path = System.IO.Path.Combine(TemplateFolder, name + ".txt");
        File.WriteAllText(path, _tplText);
        if (!TemplateNames.Contains("my: " + name)) TemplateNames.Add("my: " + name);
        _tplName = "my: " + name;
        Raise(nameof(TemplateName));
        return $"Saved {RecentFiles.Tilde(path)}";
    }

    /// <summary>The template text as React's reactions (custom set) and the React page.</summary>
    public void UseTemplateInReact()
    {
        RxSet = ReactionSets.Length - 1;
        RxText = _tplText;
        SetModule(6);
        Status = "React uses the edited template" + (TemplateTabs.Count > 1 ? "s" : "");
    }

    public void ImportTemplate(string path)
    {
        try { TemplateText = File.ReadAllText(path); Status = $"Read {System.IO.Path.GetFileName(path)}"; }
        catch (Exception e) { TemplateError = e.Message; }
    }
}
