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
    public bool IsStudioRail => _module is 8 or 9 or 17 or 19 or 25;   // Build and Export have their own rail entries

    /// <summary>Opens the builder with a SMILES (Start's quick-start box, the Molecule builder tile).</summary>
    public void OpenBuilder(string? smiles = null)
    {
        SetModule(9);
        if (smiles != null) MolSmiles = smiles;
        if (!_cgLoaded) LoadCgTemplates();
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
                var order = new[] { "gaff-amber25", "gaff-amber16", "opls2005", "cgenff", "pcff-frc", "cvff-frc" };
                var found = new Dictionary<string, CleanChoice>();
                foreach (var e in js.RootElement.GetProperty("forcefields").EnumerateArray())
                {
                    if (!e.TryGetProperty("typing", out var t) || t.ValueKind != JsonValueKind.Object) continue;
                    var id = e.GetProperty("id").GetString() ?? "";
                    var name = id == "gaff-amber25" ? "CAPS Field · GAFF2" : "CAPS Field · " + (e.GetProperty("name").GetString() ?? id);
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

    // ---------------------------------------------------------------- coarse-grained molecules (beads)
    /// <summary>Coarse-grained force fields with bead structures: MARTINI (templates from its sources), Dry MARTINI and
    /// SDK (bead SMILES; SDK also maps all-atom structures in the Force field step), Martini 3 (its molecules with their
    /// own topology: solvents, ions, small molecules, sugars, nucleobases, lipids).</summary>
    public static readonly (string Name, string File)[] CgForceFields =
        [("MARTINI · lipids", "martini-moltemplate.json"), ("MARTINI · polymers", "martini-polymers.json"),
         ("MARTINI · solvents", "martini-solvents.json"), ("MARTINI · surfactants", "martini-surfactants.json"),
         ("MARTINI · sugars", "martini-sugars.json"), ("MARTINI · amino acids", "martini-aminoacids.json"), ("Dry MARTINI", "drymartini-moltemplate.json"), ("SDK", "sdk-moltemplate.json"),
         ("Cooke–Deserno", "cooke-deserno-moltemplate.json"), ("Martini 3", "martini3.json")];
    public string[] CgFfNames => CgForceFields.Select(f => f.Name).ToArray();
    private int _cgFf;
    private string _cgBeadText = "", _cgTemplate = "";
    private readonly ObservableCollection<string> _cgTemplates = new();
    /// <summary>The chosen force field's bead templates (read the first time the list is shown).</summary>
    public ObservableCollection<string> CgTemplates { get { if (!_cgLoaded) LoadCgTemplates(); return _cgTemplates; } }
    private Dictionary<string, string> _cgTemplateText = new();
    private string? CgFfPath => Paths.ForceFields is { } dir ? System.IO.Path.Combine(dir, CgForceFields[_cgFf].File) : null;
    public int CgFf
    {
        get => _cgFf;
        set
        {
            if (Set(ref _cgFf, value)) LoadCgTemplates();
        }
    }
    private bool _cgLoaded;
    private void LoadCgTemplates()
    {
        _cgLoaded = true;
        _cgTemplates.Clear();
        _cgTemplateText = new();
        try { if (CgFfPath is { } p) _cgTemplateText = CapsDocument.BeadTemplates(p); } catch (Exception) { }
        foreach (var k in _cgTemplateText.Keys.OrderBy(k => k)) _cgTemplates.Add(k);
        Raise(nameof(CgHasTemplates));
        Raise(nameof(CgTemplateHint));
    }
    public bool CgHasTemplates => CgTemplates.Count > 0;   // (loads the list the first time)
    public string CgTemplateHint => CgHasTemplates ? "choose a template" : "none in this force field";
    public string CgTemplate
    {
        get => _cgTemplate;
        // bead SMILES are shown to edit; a molecule given term by term (Martini 3) is built by its name
        set { if (Set(ref _cgTemplate, value ?? "") && _cgTemplateText.TryGetValue(_cgTemplate, out var t)) CgBeadText = t.StartsWith('[') ? t : _cgTemplate; }
    }
    public string CgBeadText { get => _cgBeadText; set => Set(ref _cgBeadText, value ?? ""); }

    /// <summary>Builds the bead SMILES (or the chosen template) into the preview; Open in Studio takes it from there.</summary>
    public async Task BuildBeadsMolecule()
    {
        if (!_cgLoaded) LoadCgTemplates();
        var text = _cgBeadText.Trim();
        if (text.Length == 0) { MolNotes = "Choose a template or type bead SMILES, e.g. [Q0+1][Qa-1][Na]([Na][C1][C1])[C1][C1]"; return; }
        if (_cgTemplateText.TryGetValue(_cgTemplate, out var t) && (t == text || text == _cgTemplate)) text = _cgTemplate;
        var ff = CgFfPath;
        var ticket = ++_buildTicket;
        MolBuilding = true;
        try
        {
            var (doc, report) = await Task.Run(() => CapsDocument.BuildBeads(text, ff, 1, "beads"));
            if (ticket != _buildTicket) { doc.Dispose(); return; }
            var old = _molDoc;
            MolDoc = doc;
            old?.Dispose();
            using var js = JsonDocument.Parse(report);
            var r = js.RootElement;
            var beads = r.GetProperty("beads").GetDouble();
            var q = r.GetProperty("charge").GetDouble();
            MolConformers.Clear();
            MolFormula = r.GetProperty("template").GetString() is { Length: > 0 } tn ? tn : "beads";
            MolAtoms = beads.ToString("0", Inv) + " beads";
            MolNotes = $"{beads:0} beads, {r.GetProperty("bonds").GetDouble():0} bonds, charge {q:+0;−0;0} e · {CgForceFields[_cgFf].Name}: bond lengths from the force field; relax it with the force field before a run";
            Status = $"Built {MolFormula} · {beads:0} beads ({CgForceFields[_cgFf].Name})";
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
