using System;
using System.Globalization;
using System.Linq;
using System.Text.Json.Nodes;

namespace CapsStudio.ViewModels;

// Unit cell page › Point group (molecule) (core pointgroup.hpp): the Schoenflies group of the open molecule within a
// tolerance, its symmetry elements, and the molecule made exactly symmetric (an undoable edit).
public partial class MainViewModel
{
    private decimal _pgTolerance = 0.1m;
    public decimal PgTolerance { get => _pgTolerance; set => Set(ref _pgTolerance, value); }
    private string _pgText = "", _pgNote = "";
    public string PgText { get => _pgText; private set => Set(ref _pgText, value); }
    public string PgNote { get => _pgNote; private set => Set(ref _pgNote, value); }

    private void ShowPointGroup(JsonNode r)
    {
        var inv = CultureInfo.InvariantCulture;
        var g = ((string?)r["group"] ?? "").Replace("*", "∞");
        var order = (int)(r["order"]?.GetValue<double>() ?? 0);
        var els = ((JsonArray)r["elements"]!).Select(e => (string?)e!["kind"] switch
        {
            "C" => $"C{((int)e["n"]!.GetValue<double>() is var n && n > 0 ? n.ToString(inv) : "∞")}",
            "S" => $"S{(int)e["n"]!.GetValue<double>()}",
            "sigma" => "σ",
            _ => "i",
        }).GroupBy(x => x).Select(x => x.Count() > 1 ? $"{x.Count()} {x.Key}" : x.Key);
        PgText = $"{g}{(order > 0 ? $" · order {order}" : "")} · {string.Join(", ", els)}";
        PgNote = $"Images meet to {r["max_deviation"]!.GetValue<double>().ToString("0.###e0", inv)} Å" +
                 string.Concat(((JsonArray)r["notes"]!).Select(x => " · " + (string?)x));
    }

    public void FindPointGroup()
    {
        if (_doc == null) return;
        try
        {
            var r = JsonNode.Parse(_doc.PointGroup(new JsonObject { ["tolerance"] = (double)_pgTolerance }.ToJsonString()))!;
            if (r["ok"]?.GetValue<bool>() != true) { PgText = ""; PgNote = (string?)r["error"] ?? "no group"; return; }
            ShowPointGroup(r);
            Status = "Point group " + PgText;
        }
        catch (Exception e) { PgNote = e.Message; }
    }

    public void SymmetrizeMolecule()
    {
        var r = RunEdit(new { op = "symmetrize_molecule", tolerance = (double)_pgTolerance });
        if (r?["point_group"] is JsonNode pg) ShowPointGroup(pg);
        else PgNote = EditError;
    }
}
