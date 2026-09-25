using System.Collections.ObjectModel;
using System.Globalization;
using System.Text;
using System.Text.Json.Nodes;
using CapsStudio.Interop;

namespace CapsStudio.ViewModels;

/// <summary>A row of the validation table: level, what, detail and its fix.</summary>
public sealed record CheckIssue(string Level, string Title, string Detail, string Fix, int[] Atoms)
{
    public bool IsError => Level == "error";
    public bool IsWarn => Level == "warn";
    public bool IsOk => Level == "ok";
    public string Icon => Level switch { "error" => "xcircle", "warn" => "alert", _ => "check" };
    public bool HasFix => Fix.Length > 0;
    public string FixText => Fix switch { "push_apart" => "Push apart", "wrap" => "Wrap", "add_h" => "Add H", _ => "" };
}

/// <summary>Interactions &amp; checks (design/boards/Interactions): hydrogen bonds (Luzar &amp; Chandler criterion), close
/// contacts and clashes drawn in the view, and the structure's problems with a fix each.</summary>
public sealed partial class MainViewModel
{
    private bool _ixOpen;
    public bool InteractionsOpen
    {
        get => _ixOpen;
        set
        {
            if (!Set(ref _ixOpen, value)) return;
            if (value) { AppearanceOpen = false; SelectionOpen = false; RunInteractions(); }
            else { _doc?.ClearChecks(); RenderRequested?.Invoke(); }
            RaiseIxPanel();
        }
    }
    public bool ShowInteractionsPanel => IsStudio && _ixOpen && _doc != null;
    private void RaiseIxPanel() { Raise(nameof(ShowInteractionsPanel)); Raise(nameof(ShowStudioTabs)); }

    private decimal _hbD = 3.5m, _hbA = 30;
    private bool _showContacts, _showClashes = true, _ixLive = true;
    public decimal HbDistance { get => _hbD; set { if (Set(ref _hbD, Math.Clamp(value, 2.0m, 5.0m))) RunInteractions(); } }
    public decimal HbAngle { get => _hbA; set { if (Set(ref _hbA, Math.Clamp(value, 5, 90))) RunInteractions(); } }
    public bool ShowContacts { get => _showContacts; set { if (Set(ref _showContacts, value)) RunInteractions(); } }
    public bool ShowClashes { get => _showClashes; set { if (Set(ref _showClashes, value)) RunInteractions(); } }
    /// <summary>Recompute after edits and on every frame.</summary>
    public bool InteractionsLive { get => _ixLive; set => Set(ref _ixLive, value); }

    public ObservableCollection<CheckIssue> CheckIssues { get; } = new();
    private string _hbChip = "", _molChip = "", _errChip = "", _ixStatus = "", _ixError = "";
    public string HbChip { get => _hbChip; private set => Set(ref _hbChip, value); }
    public string IxMoleculesChip { get => _molChip; private set => Set(ref _molChip, value); }
    public string IxErrorsChip { get => _errChip; private set => Set(ref _errChip, value); }
    public string IxStatus { get => _ixStatus; private set => Set(ref _ixStatus, value); }
    public string IxError { get => _ixError; private set { if (Set(ref _ixError, value)) Raise(nameof(IxHasError)); } }
    public bool IxHasError => _ixError.Length > 0;
    private JsonNode? _ix;

    public void RunInteractions()
    {
        if (_doc == null || !_ixOpen) return;
        try
        {
            var o = new JsonObject { ["hb_distance"] = (double)_hbD, ["hb_angle"] = (double)_hbA, ["show_hbonds"] = true, ["show_contacts"] = _showContacts, ["show_clashes"] = _showClashes };
            var j = JsonNode.Parse(_doc.Interactions(o.ToJsonString()))!;
            if (j["ok"]?.GetValue<bool>() != true) { IxError = j["error"]?.GetValue<string>() ?? "cannot check"; return; }
            IxError = "";
            _ix = j;
            var nh = (int)j["hbonds"]!.GetValue<double>();
            var nc = (int)j["clashes"]!.GetValue<double>();
            var mols = (int)j["molecules"]!.GetValue<double>();
            HbChip = $"{nh:N0} found";
            IxMoleculesChip = $"{mols:N0} molecule{(mols == 1 ? "" : "s")}";
            CheckIssues.Clear();
            foreach (var x in j["issues"]!.AsArray())
                CheckIssues.Add(new CheckIssue(x!["level"]!.GetValue<string>(), x["title"]!.GetValue<string>(), x["detail"]!.GetValue<string>(), x["fix"]!.GetValue<string>(),
                    x["atoms"]!.AsArray().Select(a => (int)a!.GetValue<double>()).ToArray()));
            var errors = CheckIssues.Count(c => c.IsError);
            IxErrorsChip = errors == 0 ? "no errors" : $"{errors} error{(errors == 1 ? "" : "s")}";
            IxStatus = $"{_doc.Summary().Atoms:N0} atoms · {nh:N0} H-bonds · {nc:N0} clashes";
            Raise(nameof(IxHasErrors));
        }
        catch (Exception e) { IxError = e.Message; }
        RenderRequested?.Invoke();
    }
    public bool IxHasErrors => CheckIssues.Any(c => c.IsError);

    /// <summary>Applies one issue's fix.</summary>
    public async Task FixIssue(CheckIssue issue)
    {
        if (_doc == null) return;
        switch (issue.Fix)
        {
            case "wrap":
                Wrap = true;
                break;
            case "add_h":
                RunEdit(new { op = "add_h", atoms = issue.Atoms });
                break;
            case "push_apart":
            {
                // the clashing atoms' molecules, relaxed with push-off (UFF), the rest held
                var doc = _doc;
                var mols = JsonNode.Parse(doc.Select(new JsonObject { ["mode"] = "molecule", ["atoms"] = new JsonArray(issue.Atoms.Select(a => (JsonNode)a).ToArray()), ["op"] = "replace" }.ToJsonString()))!;
                var sel = JsonNode.Parse(doc.SelectionJson())!["indices"]!.AsArray().Select(x => (int)x!.GetValue<double>()).ToArray();
                doc.Select("{\"mode\":\"none\"}");
                Status = $"Pushing {sel.Length:N0} atoms apart (UFF)…";
                var text = await Task.Run(() => doc.Edit(System.Text.Json.JsonSerializer.Serialize(new { op = "clean", atoms = sel })));
                var r = JsonNode.Parse(text)!;
                if (r["ok"]?.GetValue<bool>() != true) { IxError = r["error"]?.GetValue<string>() ?? "cannot push apart"; break; }
                AfterEdit("Pushed the clashing molecules apart (UFF)");
                _ = mols;
                break;
            }
        }
        RunInteractions();
    }

    /// <summary>Every fix that does not need a choice: wrap, hydrogens, push apart.</summary>
    public async Task FixAllSafe()
    {
        foreach (var fix in new[] { "wrap", "add_h", "push_apart" })
            if (CheckIssues.FirstOrDefault(c => c.Fix == fix) is { } issue) await FixIssue(issue);
    }

    /// <summary>The H-bonds and clashes as CSV.</summary>
    public string InteractionsCsv()
    {
        var sb = new StringBuilder("kind,atom_a,atom_b,atom_c,distance_A,angle_deg\n");
        if (_ix?["hbond_list"] is JsonArray hb)
            foreach (var h in hb)
                sb.AppendLine(string.Format(CultureInfo.InvariantCulture, "hbond,{0},{1},{2},{3:0.###},{4:0.#}", (int)h!["donor"]!.GetValue<double>() + 1, (int)h["hydrogen"]!.GetValue<double>() + 1,
                    (int)h["acceptor"]!.GetValue<double>() + 1, h["distance"]!.GetValue<double>(), h["angle"]!.GetValue<double>()));
        if (_ix?["clash_list"] is JsonArray cl)
            foreach (var c in cl)
                sb.AppendLine(string.Format(CultureInfo.InvariantCulture, "clash,{0},{1},,{2:0.###},", (int)c!["i"]!.GetValue<double>() + 1, (int)c["j"]!.GetValue<double>() + 1, c["distance"]!.GetValue<double>()));
        return sb.ToString();
    }
}
