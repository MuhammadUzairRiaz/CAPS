using System.Collections.ObjectModel;
using System.Globalization;
using System.Text.Json.Nodes;

namespace CapsStudio.ViewModels;

/// <summary>Atoms of one chemical environment the assigned force field has no type for.</summary>
public sealed record UntypedGroupRow(string Environment, int Count, int[] Atoms)
{
    public string Title => $"{Count} × {Environment}";
}

/// <summary>A library force field and how well it describes the open structure.</summary>
public sealed record CoverageRow(string Id, string Name, string Status, string Detail, bool Complete, bool Current)
{
    public bool NotCurrent => !Current;
}

/// <summary>Field › coverage: when the assigned force field cannot describe the structure (atoms it has no type for,
/// parameters it lacks, or its own charges not leaving every molecule neutral), which atoms, and which force fields of
/// the library describe the whole structure instead (caps_field_coverage, run in the background).</summary>
public sealed partial class FieldViewModel
{
    public ObservableCollection<UntypedGroupRow> UntypedGroups { get; } = new();
    public ObservableCollection<CoverageRow> Alternatives { get; } = new();
    private bool _coverageBusy;
    private string _coverageNote = "", _balanceNote = "";
    private int _coverageTicket;
    public bool CoverageBusy { get => _coverageBusy; private set { if (Set(ref _coverageBusy, value)) Raise(nameof(ShowCoverage)); } }
    public string CoverageNote { get => _coverageNote; private set => Set(ref _coverageNote, value); }
    /// <summary>The physics check: the force field's own charges leave the structure at its formal charge.</summary>
    public string BalanceNote { get => _balanceNote; private set { if (Set(ref _balanceNote, value)) Raise(nameof(HasBalanceNote)); } }
    public bool HasBalanceNote => _balanceNote.Length > 0;
    public bool HasUntypedGroups => UntypedGroups.Count > 0;
    public bool HasAlternatives => Alternatives.Count > 0;
    public bool ShowCoverage => _coverageBusy || HasUntypedGroups || HasBalanceNote || HasAlternatives || _coverageNote.Length > 0;

    /// <summary>Tries every library force field on the open structure. auto: only when the assigned one is incomplete.</summary>
    public async Task CheckCoverage(bool auto = false)
    {
        var doc = _doc();
        var dir = Paths.ForceFields;
        if (doc == null || dir == null) return;
        var current = Selected?.Id;
        if (auto && Complete) { ClearCoverage(); return; }
        var ticket = ++_coverageTicket;
        CoverageBusy = true;
        CoverageNote = "Trying the library's force fields on this structure…";
        try
        {
            var finished = false;   // progress posts that arrive after the result must not overwrite it
            var text = await Task.Run(() => doc.FieldCoverage(dir, (name, f) =>
            {
                Avalonia.Threading.Dispatcher.UIThread.Post(() => { if (ticket == _coverageTicket && !finished) CoverageNote = $"Trying {name}…"; });
                return ticket == _coverageTicket;
            }));
            finished = true;
            if (ticket != _coverageTicket) return;
            var r = JsonNode.Parse(text)!;
            if ((bool?)r["ok"] != true) { CoverageNote = "Could not check the library: " + (string?)r["error"]; return; }
            var list = (r["forcefields"] as JsonArray ?? []).OfType<JsonObject>().ToList();
            UntypedGroups.Clear();
            Alternatives.Clear();
            BalanceNote = "";
            var mine = list.FirstOrDefault(x => (string?)x["id"] == current);
            if (mine != null)
            {
                foreach (var g in (mine["untyped_groups"] as JsonArray ?? []).OfType<JsonObject>())
                    UntypedGroups.Add(new UntypedGroupRow((string?)g["environment"] ?? "", (int)(double)g["count"]!,
                        (g["atoms"] as JsonArray ?? []).Select(a => (int)(double)a!).ToArray()));
                if ((bool?)mine["balanced"] == false)
                    BalanceNote = string.Format(CultureInfo.InvariantCulture,
                        "{0}'s own charges leave the structure at {1:+0.000;−0.000} e, not its formal charge: its group charges do not balance for these atoms (a type is missing or a group is typed partly). CAPS will not run it until they do.",
                        _ffName, (double?)mine["net_charge"] ?? 0);
            }
            var ok = list.Where(x => (bool?)x["complete"] == true && (string?)x["id"] != current).ToList();
            // force fields whose own charges balance first, then the rest with Gasteiger charges, UFF last
            foreach (var x in ok.OrderBy(x => (string?)x["id"] == "uff" ? 2 : (string?)x["charges"] == "types" ? 0 : 1).ThenBy(x => (string?)x["name"]))
            {
                var ch = (string?)x["charges"] switch { "types" => "its own charges, neutral", "gasteiger" => "Gasteiger charges", _ => "no charges" };
                Alternatives.Add(new CoverageRow((string?)x["id"] ?? "", (string?)x["name"] ?? "", "describes every atom", ch, true, false));
            }
            var why = mine == null ? "" : (string?)mine["status"] switch
            {
                "untyped atoms" => $"{_ffName} has no type for {(int)(double)(mine["untyped"] ?? 0)} atoms of this structure.",
                "missing parameters" => $"{_ffName} types every atom but lacks {(int)(double)(mine["missing_count"] ?? 0)} parameters (listed below).",
                "charges do not balance" => $"{_ffName} types every atom, but its charges do not balance.",
                _ => "",
            };
            CoverageNote = (why.Length > 0 ? why + " " : "") + (Alternatives.Count == 0
                ? "No force field of the library describes the whole structure; enter or import the missing parameters, or use UFF for a clean-up."
                : $"{Alternatives.Count} force field{(Alternatives.Count == 1 ? "" : "s")} of the library describe{(Alternatives.Count == 1 ? "s" : "")} it completely:");
        }
        catch (Exception e) { CoverageNote = "Could not check the library: " + e.Message; }
        finally
        {
            if (ticket == _coverageTicket) CoverageBusy = false;
            Raise(nameof(HasUntypedGroups)); Raise(nameof(HasAlternatives)); Raise(nameof(ShowCoverage));
        }
    }

    private void ClearCoverage()
    {
        ++_coverageTicket;
        UntypedGroups.Clear();
        Alternatives.Clear();
        BalanceNote = "";
        CoverageNote = "";
        CoverageBusy = false;
        Raise(nameof(HasUntypedGroups)); Raise(nameof(HasAlternatives)); Raise(nameof(ShowCoverage));
    }

    /// <summary>Assigns another library force field (from the coverage list).</summary>
    public async Task UseForceField(string id)
    {
        var k = Library.ToList().FindIndex(e => e.Id == id);
        if (k < 0) return;
        FfIndex = k;
        await Assign();
    }

    /// <summary>Selects the atoms of one untyped environment in the view.</summary>
    public void SelectGroup(UntypedGroupRow g)
    {
        var doc = _doc();
        if (doc == null) return;
        doc.Select(new JsonObject { ["mode"] = "indices", ["atoms"] = new JsonArray(g.Atoms.Select(a => (JsonNode)a).ToArray()) }.ToJsonString());
        _changed();
        _status($"Selected {g.Count} atoms: {g.Environment}");
    }
}
