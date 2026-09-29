using System.Collections.ObjectModel;
using System.Text.Json.Nodes;
using CapsStudio.Interop;

namespace CapsStudio.ViewModels;

/// <summary>A page of the theory manual.</summary>
public sealed class ManualPage : ObservableObject
{
    public string Id { get; init; } = "";
    public string Group { get; init; } = "";
    public string Title { get; init; } = "";
    public string Summary { get; init; } = "";
    public string Equation { get; init; } = "";
    public string EquationNote { get; init; } = "";
    public string[][] Symbols { get; init; } = [];
    public string When { get; init; } = "";
    public string[] Cites { get; init; } = [];
    public string Source { get; init; } = "";
    public string Tested { get; init; } = "";
    public string Deviation { get; init; } = "";
    public string[] Engines { get; init; } = [];
    private bool _on;
    public bool On { get => _on; set => Set(ref _on, value); }
}

/// <summary>A group of the manual's navigation, or a page under it.</summary>
public sealed record ManualNavItem(string Text, ManualPage? Page)
{
    public bool IsGroup => Page == null;
}

/// <summary>Theory manual (design/boards/TheoryManual): one page per method CAPS implements — what it does, its equation
/// and symbols with CAPS's settings, when to use it, the reference, where it is implemented and tested, any deviation
/// from the paper — the steps of the open structure that used it, and its BibTeX.</summary>
public sealed partial class MainViewModel
{
    public bool IsManual => _module == 41;
    private List<ManualPage>? _manual;
    public ObservableCollection<ManualNavItem> ManualNav { get; } = new();
    public ObservableCollection<ProvRow> ManualSymbols { get; } = new();
    public ObservableCollection<string> ManualRefs { get; } = new();
    public ObservableCollection<string> ManualUsed { get; } = new();
    private ManualPage? _manualPage;
    private string _manualQuery = "";
    public ManualPage? ManualCurrent { get => _manualPage; private set { if (Set(ref _manualPage, value)) RaiseManual(); } }
    public string ManualQuery { get => _manualQuery; set { if (Set(ref _manualQuery, value ?? "")) FillManualNav(); } }
    public bool ManualHasDeviation => _manualPage != null && _manualPage.Deviation.Length > 0 && _manualPage.Deviation != "none";
    public string ManualDeviation => _manualPage == null ? "" : _manualPage.Deviation.Length == 0 ? "none" : _manualPage.Deviation;
    public bool ManualHasUsed => ManualUsed.Count > 0;
    public string ManualUsedTitle => _doc == null ? "Used in this structure" : "Used in " + Title;
    public int ManualCount => _manual?.Count ?? 0;

    private List<ManualPage> Manual
    {
        get
        {
            if (_manual != null) return _manual;
            _manual = new();
            try
            {
                if (Paths.Manual is { } path)
                    foreach (var p in (JsonArray)JsonNode.Parse(File.ReadAllText(path))!["pages"]!)
                        _manual.Add(new ManualPage
                        {
                            Id = S(p, "id"), Group = S(p, "group"), Title = S(p, "title"), Summary = S(p, "summary"), Equation = S(p, "equation"),
                            EquationNote = S(p, "equation_note"), When = S(p, "when"), Source = S(p, "source"), Tested = S(p, "tested"), Deviation = S(p, "deviation"),
                            Symbols = p!["symbols"] is JsonArray sy ? sy.Select(r => ((JsonArray)r!).Select(x => x?.GetValue<string>() ?? "").ToArray()).ToArray() : [],
                            Cites = p["cites"] is JsonArray c ? c.Select(x => x!.GetValue<string>()).ToArray() : [],
                            Engines = p["engines"] is JsonArray en ? en.Select(x => x!.GetValue<string>()).ToArray() : [],
                        });
            }
            catch (Exception e) { Status = "Cannot read the theory manual: " + e.Message; }
            return _manual;
        }
    }
    private static string S(JsonNode? n, string k) => n?[k]?.GetValue<string>() ?? "";

    public void OpenManual(string? id = null)
    {
        SetModule(41);
        FillManualNav();
        ShowManualPage(Manual.FirstOrDefault(p => p.Id == id) ?? _manualPage ?? Manual.FirstOrDefault(p => p.Id == "csvr") ?? Manual.FirstOrDefault());
    }

    private void FillManualNav()
    {
        ManualNav.Clear();
        var q = _manualQuery.Trim();
        var hits = Manual.Where(p => q.Length == 0 || p.Title.Contains(q, StringComparison.OrdinalIgnoreCase) || p.Summary.Contains(q, StringComparison.OrdinalIgnoreCase)
                                     || p.Group.Contains(q, StringComparison.OrdinalIgnoreCase) || p.Cites.Any(c => c.Contains(q, StringComparison.OrdinalIgnoreCase)));
        foreach (var g in hits.GroupBy(p => p.Group))
        {
            ManualNav.Add(new ManualNavItem(g.Key, null));
            foreach (var p in g) ManualNav.Add(new ManualNavItem(p.Title, p));
        }
    }

    public void ShowManualPage(ManualPage? p)
    {
        if (p == null) return;
        foreach (var x in Manual) x.On = x == p;
        ManualCurrent = p;
        ManualSymbols.Clear();
        foreach (var r in p.Symbols) ManualSymbols.Add(new ProvRow(r.ElementAtOrDefault(0) ?? "", r.ElementAtOrDefault(1) ?? "", r.ElementAtOrDefault(2) ?? ""));
        ManualRefs.Clear();
        foreach (var c in p.Cites) { var t = CapsDocument.CitationText(c); ManualRefs.Add(t.Length > 0 ? t : c); }
        ManualUsed.Clear();
        if (_doc != null && p.Engines.Length > 0)
            try
            {
                var steps = (JsonArray)JsonNode.Parse(_doc.Provenance())!["steps"]!;
                for (int k = 0; k < steps.Count; k++)
                {
                    var e = steps[k]!["engine"]!.GetValue<string>();
                    if (p.Engines.Contains(e)) ManualUsed.Add($"step {k + 1} · {e} · {steps[k]!["summary"]!.GetValue<string>()}");
                }
            }
            catch { }
        // the project: the folder's structures whose provenance used the method (with their steps), and this session's jobs
        ManualUsedProject.Clear();
        if (p.Engines.Length > 0)
        {
            foreach (var d in ProjectDocs)
            {
                if (d.Manifest?["steps"] is not JsonArray st) continue;
                var at = new List<int>();
                for (int k = 0; k < st.Count; k++)
                    if (st[k]?["engine"]?.GetValue<string>() is string e && p.Engines.Contains(e)) at.Add(k + 1);
                if (at.Count > 0) ManualUsedProject.Add($"{d.Name} · step{(at.Count > 1 ? "s" : "")} {string.Join(", ", at)}");
            }
            foreach (var j in Jobs)
                if (JobUsed(j, p.Engines)) ManualUsedProject.Add($"{j.Id} · {j.Title} · {j.Status}");
        }
        RaiseManual();
    }

    /// <summary>A job of this session ran one of the engines: its kind is the engine's module, and its title names the
    /// method where the engine does (dynamics.npt → an NPT run, relax.lbfgs → L-BFGS).</summary>
    private static bool JobUsed(Job j, string[] engines)
    {
        static string Norm(string x) => x.ToLowerInvariant().Replace("-", "").Replace(" ", "");
        foreach (var e in engines)
        {
            var dot = e.IndexOf('.');
            var (module, method) = dot < 0 ? (e, "") : (e[..dot], e[(dot + 1)..]);
            if (!string.Equals(module, j.Kind, StringComparison.OrdinalIgnoreCase)) continue;
            var t = Norm(j.Title);
            var hit = method switch
            {
                "" => true,
                "cg" => t.Contains("conjugategradient"),
                "sd" => t.Contains("steepestdescent"),
                "larsen21" => t.Contains("larsen"),
                "protocol" => true,
                "templates" => true,
                _ => t.Contains(Norm(method)),
            };
            if (hit) return true;
        }
        return false;
    }
    public ObservableCollection<string> ManualUsedProject { get; } = new();
    public bool ManualHasUsedProject => ManualUsedProject.Count > 0;
    public string ManualUsedProjectTitle => ProjectFolder.Length > 0 ? "Used in this project · " + System.IO.Path.GetFileName(ProjectFolder.TrimEnd('/', '\\')) : "Used in this project";

    private void RaiseManual()
    {
        foreach (var n in new[] { nameof(ManualHasDeviation), nameof(ManualDeviation), nameof(ManualHasUsed), nameof(ManualUsedTitle), nameof(ManualCount),
                                  nameof(ManualHasUsedProject), nameof(ManualUsedProjectTitle) }) Raise(n);
    }

    /// <summary>BibTeX of the page's references (from CAPS's built-in table).</summary>
    public string ManualBibtex()
    {
        if (_manualPage == null) return "";
        var m = new JsonObject { ["steps"] = new JsonArray(new JsonObject { ["engine"] = "manual", ["cites"] = new JsonArray(_manualPage.Cites.Select(c => (JsonNode)c).ToArray()) }) };
        return CapsDocument.ProvenanceBibtex(m.ToJsonString());
    }
}
