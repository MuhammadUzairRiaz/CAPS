using System.Collections.ObjectModel;
using System.Globalization;
using System.Text.Json.Nodes;
using CapsStudio.Interop;

namespace CapsStudio.ViewModels;

public sealed record ColumnRow(string Column, string MapsTo, string Kind, bool Used);
public sealed record TypeRow(string Type, string Name, string Mass, string Element);
public sealed record FormatRow(string Name, string About, string Reads, string Writes);

/// <summary>Open file (design/boards/OpenLammps, OpenGromacs): the format found from the content, the first lines, how
/// dump columns map to properties, types and masses, where bonds come from and how many frames there are — then open.</summary>
public sealed partial class MainViewModel
{
    public bool IsOpenPreview => _module == 27;
    public ObservableCollection<PreviewLine> OpenHead { get; } = new();
    public ObservableCollection<ColumnRow> OpenColumns { get; } = new();
    public ObservableCollection<TypeRow> OpenTypes { get; } = new();
    public ObservableCollection<string> OpenNotes { get; } = new();
    public static readonly FormatRow[] FormatMatrix =
    [
        new("LAMMPS data", "atom_style full · molecular · charge · atomic; Masses, Bonds", "yes", "yes, with force-field sections"),
        new("LAMMPS dump", "text custom columns: id mol type element q x/xu/xs ix", "yes", "yes"),
        new("GROMACS .gro", "positions in nm, residues, box", "yes", "yes"),
        new("PDB", "ATOM/HETATM, CONECT, CRYST1", "yes", "yes"),
        new("XYZ · extended XYZ", "lattice in the comment line", "yes", "yes"),
        new("Tripos mol2", "bond orders, charges, atom types", "yes", "yes"),
        new("CIF", "crystals, symmetry expanded", "yes", "—"),
    ];
    private string _openPath = "", _openTopo = "", _openFormat = "", _openTitle = "", _openAtoms = "", _openFrames = "", _openBonds = "", _openUnits = "", _openMapped = "";
    private bool _openBusy, _openIsDump;
    private int _openGen;
    private int _returnModule = 8;

    public string OpenPath => _openPath;
    public string OpenTopology => _openTopo;
    public string OpenFileName => System.IO.Path.GetFileName(_openPath);
    public string OpenFormat { get => _openFormat; private set => Set(ref _openFormat, value); }
    public string OpenTitle { get => _openTitle; private set => Set(ref _openTitle, value); }
    public string OpenAtoms { get => _openAtoms; private set => Set(ref _openAtoms, value); }
    public string OpenFrames { get => _openFrames; private set => Set(ref _openFrames, value); }
    public string OpenBonds { get => _openBonds; private set => Set(ref _openBonds, value); }
    public string OpenUnits { get => _openUnits; private set => Set(ref _openUnits, value); }
    public string OpenMapped { get => _openMapped; private set => Set(ref _openMapped, value); }
    public bool OpenIsDump { get => _openIsDump; private set => Set(ref _openIsDump, value); }
    public bool OpenBusy { get => _openBusy; private set { if (Set(ref _openBusy, value)) Raise(nameof(OpenIdle)); } }
    public bool OpenIdle => !_openBusy;
    public string OpenButton => _openFrames.Length > 0 && _openFrames != "1" ? $"Open {_openFrames} frames" : "Open";

    /// <summary>Shows what the file holds (read off the UI thread); a dump picks up a data file beside it as its topology.</summary>
    public void PreviewOpen(string path, string? topology = null)
    {
        if (_module != 27) _returnModule = _module;
        _openPath = path;
        _openTopo = topology ?? TopologyFor(path);
        Raise(nameof(OpenPath)); Raise(nameof(OpenTopology)); Raise(nameof(OpenFileName));
        SetModule(27);
        var gen = ++_openGen;
        OpenBusy = true;
        OpenTitle = "Reading " + OpenFileName;
        var (p, t) = (_openPath, _openTopo);
        Task.Run(() =>
        {
            string? json = null, error = null;
            try { json = CapsDocument.InspectFile(p, t.Length > 0 ? t : null); } catch (Exception e) { error = e.Message; }
            Avalonia.Threading.Dispatcher.UIThread.Post(() =>
            {
                if (gen != _openGen) return;
                OpenBusy = false;
                FillOpenPreview(json, error);
            });
        });
    }

    private void FillOpenPreview(string? json, string? error)
    {
        OpenHead.Clear(); OpenColumns.Clear(); OpenTypes.Clear(); OpenNotes.Clear();
        if (json == null) { OpenTitle = "Cannot read " + OpenFileName; OpenFormat = "unknown"; OpenNotes.Add(error ?? ""); Raise(nameof(OpenButton)); return; }
        var j = JsonNode.Parse(json)!;
        var inv = CultureInfo.InvariantCulture;
        OpenFormat = (string?)j["format_name"] ?? "";
        OpenTitle = "Open · " + OpenFormat;
        OpenIsDump = (string?)j["format"] == "lammps-dump";
        int k = 0;
        foreach (var l in (JsonArray)j["head"]!) OpenHead.Add(new PreviewLine((++k).ToString(inv), (string?)l ?? ""));
        int used = 0, all = 0;
        foreach (var c in (JsonArray)j["columns"]!)
        {
            var u = (bool?)c!["used"] ?? true;
            used += u ? 1 : 0; ++all;
            OpenColumns.Add(new ColumnRow((string?)c["name"] ?? "", (string?)c["maps_to"] ?? "", (string?)c["kind"] ?? "", u));
        }
        OpenMapped = all > 0 ? $"{used} of {all} mapped" : "";
        foreach (var t in (JsonArray)j["types"]!)
            OpenTypes.Add(new TypeRow(((double?)t!["type"] ?? 0) > 0 ? ((double?)t["type"] ?? 0).ToString("0", inv) : "—", (string?)t["label"] ?? "",
                                      ((double?)t["mass"] ?? 0).ToString("F3", inv), (string?)t["element"] ?? ""));
        OpenAtoms = ((double?)j["atoms"] ?? 0).ToString("N0", inv);
        OpenFrames = ((double?)j["frames"] ?? 0).ToString("0", inv);
        OpenBonds = (string?)j["bonds_from"] ?? "";
        OpenUnits = (string?)j["units"] ?? "";
        foreach (var n in (JsonArray)j["notes"]!) OpenNotes.Add((string?)n ?? "");
        Raise(nameof(OpenButton));
    }

    public void SetOpenTopology(string? topology)
    {
        PreviewOpen(_openPath, topology ?? "");
    }

    public void CancelOpen() => SetModule(_returnModule is 27 ? 8 : _returnModule);

    /// <summary>Opens the previewed file (with its topology) and goes back to where it was asked for.</summary>
    public void ConfirmOpen()
    {
        var back = _returnModule is 27 ? 8 : _returnModule;
        SetModule(back);
        OpenRequested?.Invoke(_openPath + (_openTopo.Length > 0 ? "\n" + _openTopo : ""));
    }
}
