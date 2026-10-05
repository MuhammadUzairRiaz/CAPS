using System.Collections.ObjectModel;
using System.Globalization;
using System.Text.Json.Nodes;
using Avalonia.Media;

namespace CapsStudio.ViewModels;

/// <summary>A tag's chip in the strip over the view.</summary>
public sealed class TagChip : ObservableObject
{
    public string Name { get; init; } = "";
    public string Colour { get; init; } = "#F0A83C";
    public IBrush Brush => SolidColorBrush.Parse(Colour);
    public int Count { get; init; }
    public int Selected { get; init; }
    public string CountText => Count.ToString("N0", CultureInfo.InvariantCulture);
    public string Group { get; init; } = "";
    public bool Lit => Selected > 0 && Selected == Count;
    public string Tip => $"{Name} · {Count:N0} atoms" + (Selected > 0 ? $" ({Selected:N0} selected)" : "") +
                         $"\nClick selects it · ⇧ click adds it · ⌥ click shows only it\nLAMMPS group {Group} · GROMACS [ {Group} ] · Analyze tag:{Name}";
    private bool _renaming;
    public bool Renaming { get => _renaming; set { if (Set(ref _renaming, value)) Raise(nameof(NotRenaming)); } }
    public bool NotRenaming => !_renaming;
}

/// <summary>Tags (design/boards/Tags): any selection named and coloured, kept in the structure (NAME.tags.json beside it)
/// and a group everywhere — LAMMPS group lines, a GROMACS index group, Analyze of a tag, doc.tags in Python. An atom can
/// carry several. The strip over the view: click selects a tag, ⇧ click adds it, ⌥ click shows only it.</summary>
public sealed partial class MainViewModel
{
    public ObservableCollection<TagChip> TagChips { get; } = new();
    public bool HasTags => TagChips.Count > 0;
    public static readonly string[] TagColours = ["#F0A83C", "#5B8DEF", "#E07A5F", "#7CC784", "#9B7BD6", "#6CC4D8", "#E6C35C", "#D46BA3"];
    private string? _tagOnly;   // the tag shown alone (⌥ click), or none

    /// <summary>The strip from the structure's tags (after any change to them or to the selection).</summary>
    public void RefreshTags()
    {
        var keepRenaming = TagChips.FirstOrDefault(c => c.Renaming)?.Name;
        TagChips.Clear();
        if (_doc != null && !_doc.IsDisposed)
            try
            {
                foreach (var t in (JsonNode.Parse(_doc.TagsJson())?["tags"] as JsonArray ?? []).OfType<JsonObject>())
                    TagChips.Add(new TagChip
                    {
                        Name = (string?)t["name"] ?? "", Colour = (string?)t["colour"] ?? "#F0A83C", Count = (int)((double?)t["count"] ?? 0),
                        Selected = (int)((double?)t["selected"] ?? 0), Group = (string?)t["group"] ?? "", Renaming = (string?)t["name"] == keepRenaming,
                    });
            }
            catch { }
        Raise(nameof(HasTags));
    }

    private string TagEditJson(string op, string name, IEnumerable<int>? atoms = null, string colour = "", string newName = "")
    {
        var o = new JsonObject { ["op"] = op, ["name"] = name, ["colour"] = colour };
        if (newName.Length > 0) o["new_name"] = newName;
        if (atoms != null) o["atoms"] = new JsonArray(atoms.Select(a => (JsonNode)a).ToArray());
        return o.ToJsonString();
    }

    /// <summary>The selection tagged: a new tag, named "tag N" and ready to be renamed in the strip.</summary>
    public void TagSelection()
    {
        var a = SelectionAtoms();
        if (_doc == null || a.Length == 0) { Status = "Select atoms to tag"; return; }
        var used = TagChips.Select(c => c.Name).ToHashSet();
        var k = 1;
        while (used.Contains($"tag {k}")) ++k;
        var name = $"tag {k}";
        try { _doc.TagEdit(TagEditJson("set", name, a, TagColours[TagChips.Count % TagColours.Length])); }
        catch (Exception e) { Status = "Could not tag: " + e.Message; return; }
        RefreshTags();
        if (TagChips.FirstOrDefault(c => c.Name == name) is { } chip) chip.Renaming = true;
        TagsChanged();
        Status = $"{a.Length:N0} atoms tagged “{name}” · type its name in the strip over the view (Enter)";
    }

    /// <summary>A chip clicked: its atoms selected (add: added to the selection); only: shown alone (again: all shown).</summary>
    public void UseTag(TagChip chip, bool add, bool only)
    {
        if (_doc == null) return;
        var atoms = _doc.TagAtoms(chip.Name);
        if (only)
        {
            if (_tagOnly == chip.Name) { _tagOnly = null; ShowAllAtoms(); return; }
            _tagOnly = chip.Name;
            _doc.SetAtomState(null, 2);
            SetStates(atoms, 0, $"Only {chip.Name} shown ({atoms.Length:N0} atoms) · ⌥ click it again shows all");
            RefreshViewStates();
            return;
        }
        var r = JsonNode.Parse(_doc.Select(new JsonObject { ["mode"] = "indices", ["atoms"] = new JsonArray(atoms.Select(a => (JsonNode)a).ToArray()), ["op"] = add ? "add" : "replace" }.ToJsonString()))!;
        if (!add) { _selection.Clear(); RefreshSelection(); }
        SelectedCount = (int)((double?)r["count"] ?? 0);
        SelectHud = add ? SelectHud + " + " + chip.Name : chip.Name;
        RefreshSelBar();
        RefreshTags();
        Status = $"{SelectedCount:N0} atoms selected · tag {chip.Name}";
        RenderRequested?.Invoke();
    }

    public void RenameTag(TagChip chip, string to)
    {
        chip.Renaming = false;
        to = to.Trim();
        if (_doc == null || to.Length == 0 || to == chip.Name) { RefreshTags(); return; }
        try { _doc.TagEdit(TagEditJson("rename", chip.Name, newName: to)); Status = $"Tag {chip.Name} is now {to}"; }
        catch (Exception e) { Status = e.Message; }
        if (_tagOnly == chip.Name) _tagOnly = to;
        RefreshTags();
        TagsChanged();
    }

    public void RecolourTag(TagChip chip, string colour)
    {
        if (_doc == null) return;
        try { _doc.TagEdit(TagEditJson("colour", chip.Name, colour: colour)); } catch (Exception e) { Status = e.Message; }
        RefreshTags();
        TagsChanged();
    }

    /// <summary>The selection added to (or taken out of) a tag.</summary>
    public void TagWithSelection(TagChip chip, bool remove)
    {
        var a = SelectionAtoms();
        if (_doc == null || a.Length == 0) { Status = "Select the atoms first"; return; }
        try { _doc.TagEdit(TagEditJson(remove ? "remove" : "add", chip.Name, a)); } catch (Exception e) { Status = e.Message; return; }
        RefreshTags();
        TagsChanged();
        Status = $"{a.Length:N0} atoms {(remove ? "taken out of" : "added to")} {chip.Name}";
    }

    public void DeleteTag(TagChip chip)
    {
        if (_doc == null) return;
        try { _doc.TagEdit(TagEditJson("delete", chip.Name)); } catch (Exception e) { Status = e.Message; return; }
        if (_tagOnly == chip.Name) { _tagOnly = null; ShowAllAtoms(); }
        RefreshTags();
        TagsChanged();
        Status = $"Tag {chip.Name} deleted (its atoms stay)";
    }

    /// <summary>Tags changed: the structure is saved with them in a project (they live beside its file).</summary>
    private void TagsChanged() => QueueProjectSave();
}
