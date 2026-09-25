using System.Collections.ObjectModel;
using System.Globalization;
using System.Text.Json;
using Avalonia.Threading;
using CapsStudio.Interop;

namespace CapsStudio.ViewModels;

/// <summary>One conformer in the builder's table.</summary>
public sealed record ConformerRow(int Number, string Energy, string Method, bool Minimised);

/// <summary>A force field that can clean a built molecule up (it has typing rules).</summary>
public sealed record CleanChoice(string Name, string? File);

/// <summary>Molecule (design/boards/Sketch): a 2D sketch or a SMILES to a 3D molecule with conformers.</summary>
public sealed partial class MainViewModel
{
    private static readonly CultureInfo Inv = CultureInfo.InvariantCulture;
    public bool IsMolecule => _module == 9;
    /// <summary>The rail's Studio entry covers the molecule builder too.</summary>
    public bool IsStudioRail => _module is 8 or 9 or 13 or 17 or 18 or 19;

    /// <summary>Opens the builder with a SMILES (Start's quick-start box, the Molecule builder tile).</summary>
    public void OpenBuilder(string? smiles = null)
    {
        SetModule(9);
        if (smiles != null) MolSmiles = smiles;
    }

    // ---------------------------------------------------------------- SMILES and what it says
    private string _molSmiles = "";
    private string _molFormula = "", _molMass = "", _molStereo = "", _molAtoms = "", _molValence = "", _molError = "";
    private bool _molOk;
    /// <summary>The SMILES being built (typed, or written from the sketch).</summary>
    public string MolSmiles
    {
        get => _molSmiles;
        set { if (Set(ref _molSmiles, value)) OnSmilesTyped(); }
    }
    public string MolFormula { get => _molFormula; private set => Set(ref _molFormula, value); }
    public string MolMass { get => _molMass; private set => Set(ref _molMass, value); }
    public string MolStereo { get => _molStereo; private set => Set(ref _molStereo, value); }
    public string MolAtoms { get => _molAtoms; private set => Set(ref _molAtoms, value); }
    public string MolValence { get => _molValence; private set => Set(ref _molValence, value); }
    public string MolError { get => _molError; private set { if (Set(ref _molError, value)) Raise(nameof(MolHasError)); } }
    public bool MolHasError => _molError.Length > 0;
    public bool MolOk { get => _molOk; private set { if (Set(ref _molOk, value)) Raise(nameof(MolCanUse)); } }

    /// <summary>The sketch's graph (caps_smiles_depict JSON); the canvas redraws when it changes from outside.</summary>
    public event Action<string>? SketchGraphChanged;
    private bool _fromSketch;

    private void OnSmilesTyped()
    {
        Describe(_molSmiles);
        if (!_fromSketch && _molOk)
        {
            var json = CapsDocument.SmilesDepict(_molSmiles);
            SketchGraphChanged?.Invoke(json);
        }
        else if (!_fromSketch && _molSmiles.Trim().Length == 0) SketchGraphChanged?.Invoke("");
        ScheduleBuild();
    }

    /// <summary>The sketch was edited: its SMILES replaces the text (the drawing keeps its layout).</summary>
    public void SmilesFromSketch(string smiles)
    {
        _fromSketch = true;
        try { MolSmiles = smiles; }
        finally { _fromSketch = false; }
    }

    private void Describe(string smiles)
    {
        if (smiles.Trim().Length == 0)
        {
            MolOk = false;
            MolError = "";
            MolFormula = MolMass = MolStereo = MolAtoms = MolValence = "";
            return;
        }
        try
        {
            using var js = JsonDocument.Parse(CapsDocument.SmilesInfo(smiles));
            var r = js.RootElement;
            var ok = r.TryGetProperty("ok", out var o) && o.GetBoolean();
            var err = r.TryGetProperty("error", out var e) ? e.GetString() ?? "" : "";
            if (r.TryGetProperty("formula", out var f))
            {
                MolFormula = Subscript(f.GetString() ?? "");
                MolMass = r.GetProperty("mass").GetDouble().ToString("F2", Inv) + " g/mol";
                var sc = r.GetProperty("stereocentres").GetInt32();
                var sb = r.GetProperty("stereo_bonds").GetInt32();
                MolStereo = sb > 0 ? $"{sc} · {sb} E/Z" : sc.ToString(Inv);
                MolAtoms = $"{r.GetProperty("atoms").GetInt32()} atoms after hydrogens";
            }
            else MolFormula = MolMass = MolStereo = MolAtoms = "";
            MolValence = ok ? "ok" : err;
            MolError = ok ? "" : err;
            MolOk = ok;
        }
        catch (Exception ex) { MolOk = false; MolError = ex.Message; }
    }

    private static string Subscript(string formula)
    {
        var sb = new System.Text.StringBuilder();
        foreach (var c in formula) sb.Append(c is >= '0' and <= '9' ? (char)('₀' + (c - '0')) : c);
        return sb.ToString();
    }

    // ---------------------------------------------------------------- 3D embedding
    public CleanChoice[] CleanChoices { get; } = LoadCleanChoices();
    private int _molClean;
    public int MolClean { get => _molClean; set { if (Set(ref _molClean, value)) ScheduleBuild(); } }
    private decimal _molConfCount = 5;
    public decimal MolConfCount { get => _molConfCount; set { if (Set(ref _molConfCount, Math.Clamp(value, 1, 50))) ScheduleBuild(); } }
    public string MolMethod => "Distance bounds (CAPS)";

    private static CleanChoice[] LoadCleanChoices()
    {
        var list = new List<CleanChoice>();
        var dir = Paths.ForceFields;
        if (dir != null)
        {
            try
            {
                using var js = JsonDocument.Parse(File.ReadAllText(Path.Combine(dir, "catalogue.json")));
                var order = new[] { "gaff-amber25-dlfield", "gaff-amber16-dlfield", "opls2005-dlfield", "cgenff-dlfield", "pcff-dlfield", "cvff-dlfield" };
                var found = new Dictionary<string, CleanChoice>();
                foreach (var e in js.RootElement.GetProperty("forcefields").EnumerateArray())
                {
                    if (!e.TryGetProperty("typing", out var t) || t.ValueKind != JsonValueKind.Object) continue;
                    var id = e.GetProperty("id").GetString() ?? "";
                    var name = id == "gaff-amber25-dlfield" ? "CAPS Field · GAFF2" : "CAPS Field · " + (e.GetProperty("name").GetString() ?? id);
                    found[id] = new CleanChoice(name, Path.Combine(dir, e.GetProperty("file").GetString()!));
                }
                foreach (var id in order) if (found.Remove(id, out var c)) list.Add(c);
                list.AddRange(found.Values);
                // UFF types every element from its bonds alone (silicones, phosphazenes, halogens, metal complexes)
                list.Insert(Math.Min(1, list.Count), new CleanChoice("UFF · every element", "uff"));
            }
            catch { /* no library: embedding only */ }
        }
        list.Add(new CleanChoice("None (embedding only)", null));
        return list.ToArray();
    }

    public ObservableCollection<ConformerRow> MolConformers { get; } = new();
    private CapsDocument? _molDoc;
    public CapsDocument? MolDoc { get => _molDoc; private set { if (Set(ref _molDoc, value)) { Raise(nameof(MolHasDoc)); Raise(nameof(MolCanUse)); } } }
    public bool MolHasDoc => _molDoc != null;
    public bool MolCanUse => _molDoc != null && !_molBuilding;
    private int _molConf;
    public int MolConf
    {
        get => _molConf;
        set
        {
            if (!Set(ref _molConf, value) || _molDoc == null || value < 0) return;
            try { _molDoc.SetFrame(value); } catch { }
            MolViewChanged?.Invoke();
        }
    }
    public event Action? MolViewChanged;
    private bool _molBuilding;
    public bool MolBuilding { get => _molBuilding; private set { if (Set(ref _molBuilding, value)) Raise(nameof(MolCanUse)); } }
    private string _molNotes = "";
    public string MolNotes { get => _molNotes; private set => Set(ref _molNotes, value); }
    private string _molMethodUsed = "";
    public string MolMethodUsed { get => _molMethodUsed; private set => Set(ref _molMethodUsed, value); }

    private DispatcherTimer? _buildTimer;
    private int _buildTicket;

    private void ScheduleBuild()
    {
        if (_buildTimer == null)
        {
            _buildTimer = new DispatcherTimer { Interval = TimeSpan.FromMilliseconds(450) };
            _buildTimer.Tick += async (_, _) => { _buildTimer!.Stop(); await BuildMolecule(); };
        }
        _buildTimer.Stop();
        if (_molOk) _buildTimer.Start();
    }

    /// <summary>Builds the current SMILES (the preview updates as the sketch changes; a newer request wins).</summary>
    public async Task BuildMolecule()
    {
        if (!_molOk) return;
        var ticket = ++_buildTicket;
        var smiles = _molSmiles;
        var ff = _molClean >= 0 && _molClean < CleanChoices.Length ? CleanChoices[_molClean].File : null;
        var n = (int)_molConfCount;
        MolBuilding = true;
        try
        {
            var (doc, report) = await Task.Run(() => CapsDocument.BuildSmiles(smiles, ff, n, 1, "molecule"));
            if (ticket != _buildTicket) { doc.Dispose(); return; }
            var old = _molDoc;
            MolDoc = doc;
            old?.Dispose();
            using var js = JsonDocument.Parse(report);
            var r = js.RootElement;
            MolMethodUsed = r.GetProperty("method").GetString() ?? "";
            var shortMethod = ff == null ? "embedding" : "embed + " + CleanChoices[_molClean].Name.Replace("CAPS Field · ", "") + " min.";
            MolConformers.Clear();
            var k = 0;
            foreach (var c in r.GetProperty("conformers").EnumerateArray())
            {
                var min = c.GetProperty("minimised").GetBoolean();
                MolConformers.Add(new ConformerRow(++k, min ? c.GetProperty("rel").GetDouble().ToString("F2", Inv) : "—", min ? shortMethod : "embedding", min));
            }
            MolNotes = string.Join("\n", r.GetProperty("notes").EnumerateArray().Select(x => x.GetString()));
            _molConf = -1;
            MolConf = 0;
            Status = $"Built {MolFormula} · {MolConformers.Count} conformer{(MolConformers.Count == 1 ? "" : "s")} · {MolMethodUsed}";
        }
        catch (Exception ex)
        {
            if (ticket == _buildTicket) { MolNotes = ex.Message; Status = "Could not build: " + ex.Message; }
        }
        finally
        {
            if (ticket == _buildTicket) MolBuilding = false;
        }
    }

    /// <summary>The built molecule (the chosen conformer) becomes the Studio document.</summary>
    public void OpenMoleculeInStudio()
    {
        if (_molDoc == null) return;
        if (Busy) { Status = "Wait for the run to finish (or cancel it) first"; return; }
        var doc = _molDoc;
        var conf = Math.Max(0, _molConf);
        _molDoc = null;
        Raise(nameof(MolDoc)); Raise(nameof(MolHasDoc)); Raise(nameof(MolCanUse));
        MolConformers.Clear();
        Show(doc, (_molFormula.Length > 0 ? _molFormula : "molecule") + " · built");
        if (conf > 0) { Frame = conf; }
        GrownUnsaved = true;
        SetModule(8);
    }

    public void SaveMolecule(string path)
    {
        if (_molDoc == null) return;
        _molDoc.Save(path);
        Status = $"Saved conformer {Math.Max(0, _molConf) + 1} to {path}";
    }
}
