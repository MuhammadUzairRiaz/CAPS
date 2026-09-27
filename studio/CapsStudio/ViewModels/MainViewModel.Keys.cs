using System;
using System.Linq;
using System.Text.Json.Nodes;

namespace CapsStudio.ViewModels;

/// <summary>The keyboard map (design/boards/InteractionMap): what the keys do beyond the buttons.</summary>
public partial class MainViewModel
{
    /// <summary>The rail's pages in order, for ⌘1…⌘9.</summary>
    private static readonly int[] RailOrder = [8, 9, 0, 7, 5, 2, 4, 3, 1];
    public void GoRailPage(int k)
    {
        if (k < 1 || k > RailOrder.Length) return;
        GoModule(RailOrder[k - 1]);
    }

    /// <summary>⌘I: every atom not selected becomes selected, and the other way round.</summary>
    public void InvertSelectionKey()
    {
        if (_doc == null) return;
        var r = JsonNode.Parse(_doc.Select("{\"mode\":\"all\",\"op\":\"invert\"}"))!;
        if (r["ok"]?.GetValue<bool>() != true) return;
        SelectedCount = (int)r["count"]!.GetValue<double>();
        SelectHud = "inverted";
        Status = $"{SelectedCount:N0} atoms selected (inverted)";
        RenderRequested?.Invoke();
    }

    /// <summary>]: the selection (or the picked atoms) grows one bond outwards.</summary>
    public void GrowSelectionKey()
    {
        if (_doc == null) return;
        if (SelectedCount == 0 && _selection.Count > 0)
            _doc.Select(new JsonObject { ["mode"] = "indices", ["atoms"] = new JsonArray(_selection.Select(i => (JsonNode)i).ToArray()), ["op"] = "replace" }.ToJsonString());
        var r = JsonNode.Parse(_doc.Select("{\"mode\":\"grow\",\"steps\":1,\"op\":\"replace\"}"))!;
        if (r["ok"]?.GetValue<bool>() != true) { Status = (string?)r["error"] ?? "nothing selected to grow"; return; }
        SelectedCount = (int)r["count"]!.GetValue<double>();
        SelectHud = "grown along bonds";
        Status = $"{SelectedCount:N0} atoms selected (one bond further)";
        RenderRequested?.Invoke();
    }

    /// <summary>⌃Tab: the next structure of the project.</summary>
    public void NextStructure()
    {
        if (ProjectItems.Count < 2 || _activeItem == null) return;
        var k = ProjectItems.IndexOf(_activeItem);
        Activate(ProjectItems[(k + 1) % ProjectItems.Count]);
    }
}
