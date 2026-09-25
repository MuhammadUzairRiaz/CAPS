using System.Collections.ObjectModel;
using System.Globalization;
using System.Text.Json.Nodes;

namespace CapsStudio.ViewModels;

public sealed record ExportFormat(string Id, string Name, string About, string Extension);
public sealed record PreviewLine(string Number, string Text);

/// <summary>Export › Data (design/boards/ExportData): the structure (or the Visualize result) in a chosen format, with a
/// preview of the real file, the topology it holds and its size.</summary>
public sealed partial class MainViewModel
{
    public bool IsExport => _module == 21;

    public static readonly ExportFormat[] ExportFormats =
    [
        new("lammps-data", "LAMMPS data", "atom_style full · with force-field sections", "data"),
        new("lammps-dump", "LAMMPS dump", "every frame · id mol type xu yu zu", "lammpstrj"),
        new("gro", "GROMACS .gro", "nm, one residue per molecule", "gro"),
        new("pdb", "PDB", "CONECT, residues", "pdb"),
        new("xyz", "XYZ · extended XYZ", "lattice in the comment line", "xyz"),
        new("mol2", "Tripos mol2", "bond orders, charges, atom types", "mol2"),
    ];

    private int _exportFormat;
    private bool _exportCoeffs = true, _exportWrap, _exportPipeline, _exportBusy;
    private string _exportName = "", _exportSize = "", _exportLineCount = "", _exportError = "";
    private int _exportGen;
    public ObservableCollection<PreviewLine> ExportPreviewLines { get; } = new();
    public ObservableCollection<Row> ExportSections { get; } = new();
    public ObservableCollection<string> ExportNotes { get; } = new();

    public void OpenExport()
    {
        if (_doc == null) { Status = "Open or build a structure first"; return; }
        _exportPipeline = _module == 20 && PipelineRows.Count > 0;
        SetModule(21);
        if (_exportPipeline) ApplyPipelineForExport();
        Raise(nameof(ExportPipeline)); Raise(nameof(ExportPipelineAvailable));
        RefreshExport();
    }

    public int ExportFormatIndex { get => _exportFormat; set { if (Set(ref _exportFormat, Math.Clamp(value, 0, ExportFormats.Length - 1))) { Raise(nameof(ExportIsData)); RefreshExport(); } } }
    public bool ExportIsData => ExportFormats[_exportFormat].Id == "lammps-data";
    public bool ExportCoeffs { get => _exportCoeffs; set { if (Set(ref _exportCoeffs, value)) RefreshExport(); } }
    public bool ExportWrap { get => _exportWrap; set { if (Set(ref _exportWrap, value)) RefreshExport(); } }
    public bool ExportPipeline { get => _exportPipeline; set { if (Set(ref _exportPipeline, value)) { if (value) ApplyPipelineForExport(); RefreshExport(); } } }
    public bool ExportPipelineAvailable => PipelineRows.Count > 0;
    public string ExportName { get => _exportName; private set => Set(ref _exportName, value); }
    public string ExportFileSize { get => _exportSize; private set => Set(ref _exportSize, value); }
    public string ExportLineCount { get => _exportLineCount; private set => Set(ref _exportLineCount, value); }
    public string ExportError { get => _exportError; private set { if (Set(ref _exportError, value)) Raise(nameof(ExportHasError)); } }
    public bool ExportHasError => _exportError.Length > 0;
    public bool ExportBusy { get => _exportBusy; private set { if (Set(ref _exportBusy, value)) Raise(nameof(ExportIdle)); } }
    public bool ExportIdle => !_exportBusy;

    private string ExportOptionsJson() => new JsonObject { ["pipeline"] = _exportPipeline, ["wrap"] = _exportWrap, ["coeffs"] = _exportCoeffs }.ToJsonString();
    public string ExportSuggestedName => $"{System.IO.Path.GetFileNameWithoutExtension(_doc?.Path ?? "structure")}.{ExportFormats[_exportFormat].Extension}";

    // the core keeps the pipeline only while Visualize is open; exporting its result sets it again
    private void ApplyPipelineForExport()
    {
        if (_doc == null || PipelineRows.Count == 0) return;
        try { _doc.SetPipeline(PipelineJson()); } catch { }
    }

    /// <summary>Writes the chosen export to a scratch file off the UI thread and shows its head, size and sections.</summary>
    public void RefreshExport()
    {
        if (_doc == null || !IsExport) return;
        var gen = ++_exportGen;
        var doc = _doc;
        var fmt = ExportFormats[_exportFormat];
        var opts = ExportOptionsJson();
        ExportName = ExportSuggestedName;
        ExportBusy = true;
        Task.Run(() =>
        {
            string? json = null, error = null;
            try { json = doc.ExportPreview(fmt.Id, opts, 80); } catch (Exception e) { error = e.Message; }
            Avalonia.Threading.Dispatcher.UIThread.Post(() =>
            {
                if (gen != _exportGen) return;
                ExportBusy = false;
                ExportError = error ?? "";
                ExportPreviewLines.Clear();
                ExportSections.Clear();
                ExportNotes.Clear();
                if (json == null) return;
                var j = JsonNode.Parse(json)!;
                var inv = CultureInfo.InvariantCulture;
                int k = 0;
                foreach (var l in (JsonArray)j["lines"]!) ExportPreviewLines.Add(new PreviewLine((++k).ToString(inv), (string?)l ?? ""));
                var bytes = (double?)j["bytes"] ?? 0;
                ExportFileSize = bytes > 1 << 20 ? (bytes / 1048576).ToString("F1", inv) + " MB" : (bytes / 1024).ToString("F0", inv) + " KB";
                ExportLineCount = ((double?)j["line_count"] ?? 0).ToString("N0", inv) + " lines";
                string C(string key) => j[key] is JsonNode n ? ((double?)n ?? 0).ToString("N0", inv) : "–";
                ExportSections.Add(new Row("Atoms", C("atoms") + (j["atom_types"] != null ? " · " + C("atom_types") + " types" : "")));
                ExportSections.Add(new Row("Bonds", C("bonds") + (j["bond_types"] != null ? " · " + C("bond_types") + " types" : "")));
                if (j["angles"] != null) ExportSections.Add(new Row("Angles", C("angles") + (j["angle_types"] != null ? " · " + C("angle_types") + " types" : "")));
                if (j["dihedrals"] != null) ExportSections.Add(new Row("Dihedrals", C("dihedrals") + (j["dihedral_types"] != null ? " · " + C("dihedral_types") + " types" : "")));
                if (j["impropers"] != null) ExportSections.Add(new Row("Impropers", C("impropers") + (j["improper_types"] != null ? " · " + C("improper_types") + " types" : "")));
                foreach (var n in (JsonArray)j["notes"]!) ExportNotes.Add((string?)n ?? "");
            });
        });
    }

    public async Task ExportNow(string path)
    {
        if (_doc == null) return;
        var doc = _doc;
        var fmt = ExportFormats[_exportFormat];
        var opts = ExportOptionsJson();
        ExportBusy = true;
        try
        {
            await Task.Run(() => doc.ExportData(path, fmt.Id, opts));
            Status = $"Wrote {fmt.Name} · {path}";
            Remember(path, null);
        }
        catch (Exception e) { Status = "Export failed: " + e.Message; ExportError = e.Message; }
        finally { ExportBusy = false; }
    }
}
