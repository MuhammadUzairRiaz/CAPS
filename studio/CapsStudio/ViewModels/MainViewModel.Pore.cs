using System.Globalization;
using System.Text.Json.Nodes;
using CapsStudio.Interop;

namespace CapsStudio.ViewModels;

/// <summary>A fluid for a pore: its name and SMILES.</summary>
public sealed record PoreFluid(string Name, string Smiles, bool UnitedAtom = false)
{
    public override string ToString() => Name;
}

/// <summary>Pores in the nanostructure builder (design/boards/SlitPore): a slit between graphene walls, a cylindrical
/// channel carved from a crystal, or a framework's supercell, filled with a fluid packed inside the pore only.</summary>
public sealed partial class MainViewModel
{
    public bool NanoIsPore => _nanoKind == 3;
    public static readonly string[] PoreTypes = ["Slit", "Cylinder", "From CIF (zeolite, MOF)"];
    public static readonly string[] PoreWallChoices = ["1 graphene layer each", "2 layers each (AB)", "3 layers each (ABA)"];
    public static readonly PoreFluid[] PoreFluids =
    [
        new("Methane · all-atom", "C"), new("Methane · TraPPE-UA (one site)", "C", true), new("Ethane · TraPPE-UA", "CC", true),
        new("n-Butane · TraPPE-UA", "CCCC", true), new("n-Hexane · TraPPE-UA", "CCCCCC", true), new("Carbon dioxide", "O=C=O"), new("Water", "O"), new("Nitrogen", "N#N"),
        new("Ethanol", "CCO"), new("n-Hexane", "CCCCCC"), new("Toluene", "Cc1ccccc1"), new("Empty (no fluid)", ""),
    ];
    private int _poreType, _poreWalls, _poreFluid, _poreCrystal;
    private decimal _poreWidth = 10, _poreLx = 26, _poreLy = 22, _poreCount = 16, _poreWall = 6, _poreLength = 20, _poreRepeat = 2;
    private bool _poreVacuum, _porePassivate = true;
    public int PoreType { get => _poreType; set { if (Set(ref _poreType, value)) { RaisePore(); NanoPreview(); } } }
    public bool PoreIsSlit => _poreType == 0;
    public bool PoreIsCylinder => _poreType == 1;
    public bool PoreIsFramework => _poreType == 2;
    public bool PoreNeedsCrystal => _poreType != 0;
    public int PoreWalls { get => _poreWalls; set { if (Set(ref _poreWalls, value)) NanoPreview(); } }
    public decimal PoreWidth { get => _poreWidth; set { if (Set(ref _poreWidth, Math.Clamp(value, 3.4m, 200m))) { RaisePore(); NanoPreview(); } } }
    public decimal PoreLx { get => _poreLx; set { if (Set(ref _poreLx, Math.Clamp(value, 5, 300))) NanoPreview(); } }
    public decimal PoreLy { get => _poreLy; set { if (Set(ref _poreLy, Math.Clamp(value, 5, 300))) NanoPreview(); } }
    public bool PoreVacuum { get => _poreVacuum; set { if (Set(ref _poreVacuum, value)) NanoPreview(); } }
    public int PoreCrystal { get => _poreCrystal; set { if (Set(ref _poreCrystal, value)) { Raise(nameof(PoreCrystalItem)); NanoPreview(); } } }
    /// <summary>The crystal as an item (index bindings go blank when the list is filled after the page).</summary>
    public CrystalEntry? PoreCrystalItem
    {
        get => _poreCrystal >= 0 && _poreCrystal < Crystals.Count ? Crystals[_poreCrystal] : null;
        set { if (value != null && Crystals.IndexOf(value) is var k and >= 0) PoreCrystal = k; }
    }
    public decimal PoreWall { get => _poreWall; set { if (Set(ref _poreWall, Math.Clamp(value, 2, 50))) NanoPreview(); } }
    public decimal PoreLength { get => _poreLength; set { if (Set(ref _poreLength, Math.Clamp(value, 3, 300))) NanoPreview(); } }
    public decimal PoreRepeat { get => _poreRepeat; set { if (Set(ref _poreRepeat, Math.Clamp(Math.Round(value), 1, 8))) NanoPreview(); } }
    public bool PorePassivate { get => _porePassivate; set { if (Set(ref _porePassivate, value)) NanoPreview(); } }
    public int PoreFluidIndex { get => _poreFluid; set { if (Set(ref _poreFluid, value)) { RaisePore(); NanoPreview(); } } }
    public decimal PoreCount { get => _poreCount; set { if (Set(ref _poreCount, Math.Clamp(Math.Round(value), 0, 5000))) { RaisePore(); NanoPreview(); } } }
    public bool PoreHasFluid => PoreFluids[_poreFluid].Smiles.Length > 0;
    public string PoreWidthLabel => _poreType == 1 ? "Diameter D · Å" : "Width H · Å";
    public string PoreWidthNote => _poreType == 1 ? "channel along the crystal's c axis" : "carbon centre to centre";
    private string _poreChip = "";
    public string PoreChip { get => _poreChip; private set => Set(ref _poreChip, value); }

    private void RaisePore()
    {
        foreach (var n in new[] { nameof(PoreCrystalItem), nameof(NanoIsPore), nameof(PoreIsSlit), nameof(PoreIsCylinder), nameof(PoreIsFramework), nameof(PoreNeedsCrystal), nameof(PoreHasFluid),
                                  nameof(PoreWidthLabel), nameof(PoreWidthNote) }) Raise(n);
        RaiseNano();
    }

    private string PoreOptions()
    {
        var f = PoreFluids[_poreFluid];
        var o = new JsonObject
        {
            ["kind"] = _poreType switch { 1 => "cylinder", 2 => "framework", _ => "slit" },
            ["width"] = (double)_poreWidth, ["layers"] = _poreWalls + 1, ["lx"] = (double)_poreLx, ["ly"] = (double)_poreLy,
            ["vacuum"] = _poreVacuum, ["wall"] = (double)_poreWall, ["length"] = (double)_poreLength, ["passivate"] = _porePassivate,
            ["repeat"] = new JsonArray((int)_poreRepeat, (int)_poreRepeat, (int)_poreRepeat),
            ["fluid"] = f.Smiles, ["fluid_name"] = f.Name.Split(" · ")[0], ["count"] = f.Smiles.Length > 0 ? (int)_poreCount : 0, ["united_atom"] = f.UnitedAtom,
        };
        if (_poreType != 0) o["cif"] = _poreCrystal < Crystals.Count ? Crystals[_poreCrystal].File : "";
        return o.ToJsonString();
    }

    private string PoreTitle => _poreType switch
    {
        1 => $"{(_poreCrystal < Crystals.Count ? Crystals[_poreCrystal].Name : "crystal")} channel · D = {_poreWidth:0.#} Å",
        2 => $"{(_poreCrystal < Crystals.Count ? Crystals[_poreCrystal].Name : "crystal")} framework",
        _ => $"graphite slit · H = {_poreWidth:0.#} Å",
    };

    /// <summary>The pore's report (JSON) as the chips and notes under the preview.</summary>
    private string PoreReportText(string json)
    {
        try
        {
            var r = JsonNode.Parse(json)!;
            var inv = CultureInfo.InvariantCulture;
            var f = PoreFluids[_poreFluid];
            var wall = r["wall_atoms"]!.GetValue<int>();
            var fl = r["fluid_molecules"]!.GetValue<int>();
            PoreChip = fl > 0 ? $"{fl} × {f.Name.Split(" · ")[0]}" + (r["fluid_density"]!.GetValue<double>() > 0 ? $" · {r["fluid_density"]!.GetValue<double>().ToString("0.00", inv)} g/cm³ in the pore" : "") : "empty pore";
            return string.Join("\n", ((JsonArray)r["notes"]!).Select(n => n!.GetValue<string>())) + $"\n{wall:N0} wall atoms (molecule 1, held in Relax and Dynamics)";
        }
        catch { return json; }
    }
}
