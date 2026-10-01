using System.Collections.ObjectModel;
using System.Globalization;
using System.Text.Json.Nodes;
using CapsStudio.Interop;

namespace CapsStudio.ViewModels;

/// <summary>A small molecule for Grow's cell (a solvent or a gas): what it is and how many go in.</summary>
public sealed record SmallChoice(string Name, string Smiles, string Kind)
{
    public override string ToString() => $"{Name} · {Kind}";
}

public sealed class GrowSmallRow : ObservableObject
{
    public string Name { get; init; } = "";
    public string Smiles { get; init; } = "";
    private decimal _count = 20;
    public decimal CountD { get => _count; set => Set(ref _count, Math.Clamp(Math.Round(value), 1, 100000)); }
}

/// <summary>Grow › Components: solvents and gases put into the grown cell's free space (caps_insert_molecules, CAPS Pack).</summary>
public sealed partial class MainViewModel
{
    private List<SmallChoice>? _smallChoices;
    /// <summary>The solvent library (Solvation) and common gases.</summary>
    public List<SmallChoice> GrowSmallChoices
    {
        get
        {
            if (_smallChoices != null) return _smallChoices;
            var list = new List<SmallChoice>();
            try
            {
                foreach (var s in JsonNode.Parse(CapsDocument.SolventLibrary())!["solvents"]!.AsArray())
                    list.Add(new SmallChoice(s!["name"]!.GetValue<string>(), s["smiles"]!.GetValue<string>(), "solvent"));
            }
            catch { }
            list.AddRange([new("Carbon dioxide", "O=C=O", "gas"), new("Nitrogen", "N#N", "gas"), new("Oxygen", "O=O", "gas"), new("Methane", "C", "gas")]);
            return _smallChoices = list;
        }
    }
    public ObservableCollection<GrowSmallRow> GrowSmall { get; } = new();
    public bool GrowHasSmall => GrowSmall.Count > 0;
    private bool _growSmallOpen;
    public bool GrowSmallOpen { get => _growSmallOpen; set => Set(ref _growSmallOpen, value); }
    private int _growSmallPick;
    public int GrowSmallPick { get => _growSmallPick; set => Set(ref _growSmallPick, value); }

    public void AddGrowSmall()
    {
        if (_growSmallPick < 0 || _growSmallPick >= GrowSmallChoices.Count) return;
        var c = GrowSmallChoices[_growSmallPick];
        if (GrowSmall.Any(r => r.Smiles == c.Smiles)) { Status = $"{c.Name} is in the cell's components already"; return; }
        GrowSmall.Add(new GrowSmallRow { Name = c.Name, Smiles = c.Smiles, CountD = c.Kind == "gas" ? 10 : 30 });
        GrowSmallOpen = false;
        Raise(nameof(GrowHasSmall));
        Status = $"{c.Name} goes into the free space after the chains are grown";
    }

    public void RemoveGrowSmall(GrowSmallRow r) { GrowSmall.Remove(r); Raise(nameof(GrowHasSmall)); }

    /// <summary>The groups of a grown cell with small molecules: every component by Grow's force field (null: one component).</summary>
    private string? GrowGroupSpec(CapsDocument doc)
    {
        if (!ReferenceEquals(_growParts.Doc, doc) || _growParts.Parts.Count < 2) return null;
        var idx = GrowFfIndex;
        if (idx < 0 || idx >= Field.Library.Count) return null;
        var groups = new System.Text.Json.Nodes.JsonArray();
        var used = new HashSet<string>(StringComparer.OrdinalIgnoreCase);
        foreach (var (name, mols) in _growParts.Parts)
        {
            var n = name;
            for (var k = 2; !used.Add(n); ++k) n = $"{name}_{k}";
            groups.Add(new System.Text.Json.Nodes.JsonObject { ["name"] = n, ["molecules"] = mols, ["forcefield"] = Field.Library[idx].File, ["charges"] = FieldViewModel.CoreChargesOf(_growCharges) });
        }
        return new System.Text.Json.Nodes.JsonObject { ["groups"] = groups }.ToJsonString();
    }

    /// <summary>After the chains: each small molecule inserted into the free space (2 Å between molecules); the lines
    /// for the report. Throws when a count does not fit.</summary>
    /// <summary>The grown cell's components (the chains, then each small molecule) and their molecules, for a force field
    /// per component: each its own atom types, and the LAMMPS input groups them (as Packing and the blend do).</summary>
    private (CapsDocument? Doc, List<(string Name, string Molecules)> Parts) _growParts = (null, new());

    private string InsertGrowSmall(CapsDocument doc, ulong seed)
    {
        var sb = new System.Text.StringBuilder();
        var k = 0UL;
        var parts = new List<(string, string)>();
        var before = (int)doc.Summary().Molecules;
        if (before > 0) parts.Add(("chains", $"1-{before}"));
        foreach (var r in GrowSmall.ToList())
        {
            var rep = doc.InsertMolecules(r.Smiles, (int)r.CountD, 2.0, seed + 17 * ++k);
            var now = (int)doc.Summary().Molecules;
            if (now > before) parts.Add((r.Name, $"{before + 1}-{now}"));
            before = now;
            sb.Append(CultureInfo.InvariantCulture, $"{r.CountD:0} × {r.Name}: ").Append(rep.Split('\n').FirstOrDefault() ?? "").Append('\n');
        }
        _growParts = (doc, parts);
        return sb.ToString();
    }
}
