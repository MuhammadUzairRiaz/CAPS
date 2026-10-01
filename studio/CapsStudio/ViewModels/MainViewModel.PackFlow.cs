using System;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Threading.Tasks;

namespace CapsStudio.ViewModels;

/// <summary>Pack as its own workflow (as Materials Studio keeps Amorphous Cell's Packing apart from Construction):
/// molecules packed into an empty box, or around the current structure (a grown cell, a surface), which stays fixed;
/// its own force field, assigned when packing finishes; then the cell goes straight to LAMMPS or GROMACS, with
/// minimisation and dynamics in CAPS as optional steps.</summary>
public partial class MainViewModel
{
    public static readonly string[] PackStartModes = ["An empty box", "Around the current structure (kept fixed; molecules fill its cell)"];
    private int _packStart;
    public int PackStart
    {
        get => _packStart;
        set
        {
            var v = value == 1 && !PackCanUseCurrent ? 0 : Math.Clamp(value, 0, 1);
            if (!Set(ref _packStart, v)) return;
            if (v == 1 && PackHostDoc() is { } host)
            {
                // the box is the structure's cell: new structure blocks go inside it
                var s = host.Summary();
                SetPackBox(s.CellA, s.CellB, s.CellC);
            }
            Raise(nameof(PackStartNote));
            Raise(nameof(PackIntoCurrent));
        }
    }
    /// <summary>Packing around a structure needs its cell (orthorhombic: packmol's boxes are axis-aligned).</summary>
    public bool PackCanUseCurrent
    {
        get
        {
            if (_doc == null) return false;
            var s = _doc.Summary();
            return s.CellValid != 0 && Math.Abs(s.Volume - s.CellA * s.CellB * s.CellC) < 1e-6 * s.Volume;   // V = abc only with right angles
        }
    }
    public string PackStartNote => _packStart == 1 && _doc != null
        ? PackRepacking
            ? string.Format(CultureInfo.InvariantCulture, "Packing again starts from {0} (as before the last packing) and replaces the packed cell: change the molecules or counts below, then Pack",
                _packHostTitle)
            : string.Format(CultureInfo.InvariantCulture, "{0} stays where it is; the molecules below fill the free space of its {1:0.#} × {2:0.#} × {3:0.#} Å cell (periodic)",
                Title.Replace(" (unsaved)", ""), _packX, _packY, _packZ)
        : "The molecules below are packed into the box of the input";

    // ---- packing again: the structure packed around is remembered, and the new packed cell replaces the last one
    private CapsStudio.Interop.CapsDocument? _packHost, _packResult, _packRunHost;
    private string _packHostTitle = "";
    /// <summary>The packed cell is open and its host still is: Pack starts from the host again.</summary>
    public bool PackRepacking => _packResult != null && ReferenceEquals(_doc, _packResult) && _packHost is { IsDisposed: false };
    /// <summary>What Pack packs around: the open structure, or (the last packed cell open) the structure it was packed around.</summary>
    private CapsStudio.Interop.CapsDocument? PackHostDoc() => PackRepacking ? _packHost : _doc;

    /// <summary>The box (and the molecule blocks that filled the old box) follow a new cell.</summary>
    private void SetPackBox(double a, double b, double c)
    {
        var inv = CultureInfo.InvariantCulture;
        var old = string.Format(inv, "{0:0.###} {1:0.###} {2:0.###}", _packX, _packY, _packZ);
        PackXD = (decimal)a; PackYD = (decimal)b; PackZD = (decimal)c;
        var now = string.Format(inv, "{0:0.###} {1:0.###} {2:0.###}", _packX, _packY, _packZ);
        if (old != now && _packText.Contains(old, StringComparison.Ordinal))
            PackText = _packText.Replace("to " + old, "to " + now).Replace("0. 0. 0. " + old, "0. 0. 0. " + now);
    }

    /// <summary>Another structure is open: the last run's result (banner, guarantee, convergence, log) is not its, so the
    /// page starts fresh for it; packed around a structure, the box follows its cell. The packed cell itself keeps them.</summary>
    private void PackDocumentChanged()
    {
        if (_doc != null && ReferenceEquals(_doc, _packResult)) { Raise(nameof(PackStartNote)); return; }
        PackDone = false;
        PackCurve.Clear();
        PackCurveChanged?.Invoke();
        PackLoop = 0; PackBad = 0; PackDmin = "—";
        PackLog = "Add molecules (or open an input: CAPS or packmol syntax), then Pack.";
        if (_packStart == 1)
        {
            if (!PackCanUseCurrent) PackStart = 0;
            else if (_doc != null) { var s = _doc.Summary(); SetPackBox(s.CellA, s.CellB, s.CellC); }
        }
        Raise(nameof(PackStartNote));
    }

    // ---- the Pack's own force field (Grow keeps its own too)
    private int _packFf = -1, _packCharges;
    private bool _packAssign = true;
    public int PackFfIndex { get => _packFf < 0 ? Field.FfIndex : _packFf; set { _packFf = value; Raise(); } }
    public int PackChargeMode { get => _packCharges; set => Set(ref _packCharges, value); }
    public bool PackAssignField { get => _packAssign; set => Set(ref _packAssign, value); }

    // ---- after packing: export straight away, or minimise / run dynamics first
    private bool _packDone;
    public bool PackDone { get => _packDone; private set => Set(ref _packDone, value); }

    /// <summary>The packmol text actually run: with the current structure as a fixed block and its cell as the periodic box.</summary>
    // stage 3: pack loosely, then compress the cell to a density (push-off and minimisation between affine steps)
    /// <summary>Packing around the current structure: it stays fixed, so the cell cannot be compressed.</summary>
    public bool PackIntoCurrent => _packStart == 1;
    private bool _packCompress;
    private decimal _packCompressTo = 0.9m;
    public bool PackCompress { get => _packCompress; set => Set(ref _packCompress, value); }
    public decimal PackCompressTo { get => _packCompressTo; set => Set(ref _packCompressTo, Math.Clamp(value, 0.05m, 5m)); }
    private string PackTextToRun()
    {
        var text = PackTextToRunCore();
        if (_packCompress && _packStart != 1 && !text.Split('\n').Any(l => l.TrimStart().StartsWith("compress ", StringComparison.OrdinalIgnoreCase)))
            text = string.Format(CultureInfo.InvariantCulture, "compress {0:0.###}\n", _packCompressTo) + text;
        return text;
    }
    private string PackTextToRunCore()
    {
        _packRunHost = null;
        if (_packStart != 1 || _doc == null) return _packText;
        var hostDoc = PackHostDoc()!;
        var hs = hostDoc.Summary();
        if (hs.CellValid == 0 || Math.Abs(hs.Volume - hs.CellA * hs.CellB * hs.CellC) > 1e-6 * hs.Volume)
            throw new InvalidOperationException("the current structure has no orthorhombic cell to pack into");
        _packRunHost = hostDoc;
        var dir = Path.Combine(Path.GetTempPath(), "caps-pack");
        Directory.CreateDirectory(dir);
        var host = Path.Combine(dir, "host_" + Environment.ProcessId + ".data");
        hostDoc.Save(host);
        var s = hs;
        var inv = CultureInfo.InvariantCulture;
        var caps = IsCapsPack(_packText);
        var cell = caps ? string.Format(inv, "cell      {0:0.####} {1:0.####} {2:0.####}", s.CellA, s.CellB, s.CellC)
                        : string.Format(inv, "pbc 0. 0. 0. {0:0.####} {1:0.####} {2:0.####}", s.CellA, s.CellB, s.CellC);
        var lines = _packText.Split('\n').Where(l => !l.TrimStart().StartsWith(caps ? "cell " : "pbc ", StringComparison.OrdinalIgnoreCase)).ToList();
        var first = lines.FindIndex(l => l.TrimStart().StartsWith(caps ? "molecule " : "structure ", StringComparison.OrdinalIgnoreCase));
        if (first < 0) first = lines.Count;
        var title = ReferenceEquals(hostDoc, _doc) ? Title.Replace(" (unsaved)", "") : _packHostTitle;
        lines.Insert(first, caps ? $"{cell}\n\nmolecule  {host}   # {title}, kept where it is\n  count   1\n  fixed   at 0 0 0\nend\n"
                                 : $"{cell}\n\nstructure {host}   # {title}, kept where it is\n  number 1\n  fixed 0. 0. 0. 0. 0. 0.\nend structure\n");
        return string.Join('\n', lines);
    }

    // ---- fill to a density: the number of the last molecule block from the target density of the whole cell
    private decimal _fillDensity = 0.9m;
    private string _fillText = "";
    public decimal? FillDensity { get => _fillDensity; set { if (value != null) Set(ref _fillDensity, Math.Clamp(value.Value, 0.01m, 5m)); } }
    public string FillText { get => _fillText; private set => Set(ref _fillText, value); }
    /// <summary>Sets the last structure block's number so the cell (the current structure plus the molecules) reaches the
    /// density: N = (ρ V N_A − m_host) / m_molecule, rounded down (ρ in g/cm³, V the cell's volume).</summary>
    public void FillToDensity()
    {
        var inv = CultureInfo.InvariantCulture;
        if (_doc == null || _packStart != 1) { FillText = "Pack around the current structure first (Pack into › Around the current structure)"; return; }
        var lines = _packText.Split('\n').ToList();
        var caps = IsCapsPack(_packText);
        string open = caps ? "molecule " : "structure ", countKey = caps ? "count " : "number ", endKey = caps ? "end" : "end structure";
        var k = lines.FindLastIndex(l => l.TrimStart().StartsWith(open, StringComparison.OrdinalIgnoreCase));
        if (k < 0) { FillText = "Add the molecule to fill with first"; return; }
        var path = lines[k].Trim()[open.Length..].Split('#')[0].Trim();
        if (!Path.IsPathRooted(path)) path = Path.Combine(_packBaseDir, path);
        double mMol;
        try { using var mol = Interop.CapsDocument.Open(path); mMol = mol.Summary().TotalMass; }
        catch (Exception e) { FillText = "Cannot read the molecule: " + e.Message; return; }
        var s = _doc.Summary();
        const double avogadroPerA3 = 0.602214076;   // g/cm³ × Å³ → g/mol
        var target = (double)_fillDensity * s.Volume * avogadroPerA3;
        var room = target - s.TotalMass;
        if (mMol <= 0) { FillText = "The molecule has no mass"; return; }
        if (room < mMol) { FillText = string.Format(inv, "The cell is already at {0:0.000} g/cm³: nothing to add for {1:0.000}", s.Density, _fillDensity); return; }
        var n = (int)Math.Floor(room / mMol);
        var j = k + 1;
        while (j < lines.Count && !lines[j].TrimStart().StartsWith(endKey, StringComparison.OrdinalIgnoreCase) && !lines[j].TrimStart().StartsWith(countKey, StringComparison.OrdinalIgnoreCase)) ++j;
        if (j < lines.Count && lines[j].TrimStart().StartsWith(countKey, StringComparison.OrdinalIgnoreCase)) lines[j] = $"  {countKey}{n}";
        else lines.Insert(k + 1, $"  {countKey}{n}");
        PackText = string.Join('\n', lines);
        var reached = (s.TotalMass + n * mMol) / (s.Volume * avogadroPerA3);
        FillText = string.Format(inv, "{0} × {1} ({2:0.0} g/mol): {3:0.000} → {4:0.000} g/cm³ (target {5:0.000}; whole molecules)", n, Path.GetFileNameWithoutExtension(path), mMol, s.Density, reached, _fillDensity);
    }

    /// <summary>After Grow or Pack made a cell: that builder's own force field, typed and checked.</summary>
    private async Task AssignForBuilder(bool pack)
    {
        if (_doc == null) return;
        if (pack ? !_packAssign : !_growAssignFf) return;
        var idx = pack ? PackFfIndex : GrowFfIndex;
        if (idx < 0 || idx >= Field.Library.Count) return;
        if (pack && PackGroupSpec(_doc) is { } spec)   // rows with their own force fields: a group each, the rest by the Pack's
        {
            // a water model row: its geometry, charges and (four-site) M sites on the waters first
            var wm = System.Text.RegularExpressions.Regex.Match(spec, "\"water\":\"([^\"]+)\"");
            if (wm.Success) RunEdit(new { op = "water_model", model = wm.Groups[1].Value });
            Field.FfIndex = idx;
            await Field.AssignGroupsJson(spec, "Assigned by molecule (Pack)");
            if (!Field.Assigned && !spec.Contains("\"water\"", StringComparison.Ordinal))   // not by component (united-atom …): the whole cell
            {
                var why = Field.Log;
                Field.ChargeMode = _packCharges;
                await Field.Assign();
                Status = $"Assigned to the whole cell (not by component: {why.Split('\n')[0]})";
            }
            _pipeAutoFf = true;
            RaiseGrowField();
            RefreshSteps();
            return;
        }
        if (!pack && GrowGroupSpec(_doc) is { } gspec)   // chains and small molecules: a group each, Grow's force field
        {
            Field.FfIndex = idx;
            await Field.AssignGroupsJson(gspec, "Assigned by component (Polymer cell)");
            if (Field.Assigned)
            {
                _pipeAutoFf = true;
                RaiseGrowField();
                RefreshSteps();
                return;
            }
            // a force field that cannot type the components apart (united-atom: CH groups become one site) — the whole
            // cell as one, and why
            var why = Field.Log;
            Field.ChargeMode = _growCharges;
            await Field.Assign();
            Status = $"Assigned to the whole cell (not by component: {why.Split('\n')[0]})";
            _pipeAutoFf = true;
            RaiseGrowField();
            RefreshSteps();
            return;
        }
        Field.FfIndex = idx;
        Field.ChargeMode = pack ? _packCharges : _growCharges;
        await Field.Assign();
        _pipeAutoFf = true;
        RaiseGrowField();
        RefreshSteps();
    }

    // ---- Grow's own force field
    private int _growFf = -1, _growCharges;
    public int GrowFfIndex { get => _growFf < 0 ? Field.FfIndex : _growFf; set { _growFf = value; Raise(); } }
    public int GrowChargeMode { get => _growCharges; set => Set(ref _growCharges, value); }

    /// <summary>Grow › Assign now: Grow's own force field on the structure in the viewer.</summary>
    public async Task AssignGrowFieldNow()
    {
        if (_doc == null) return;
        var idx = GrowFfIndex;
        if (idx >= 0 && idx < Field.Library.Count) Field.FfIndex = idx;
        Field.ChargeMode = _growCharges;
        await Field.Assign();
        RaiseGrowField();
        RefreshSteps();
    }

    // ---- a force field per molecule: each row may take its own (a library force field or a water model); the rest
    // of the cell (a polymer packed around, the rows left to it) takes the Pack's force field above
    public List<string> PackFfChoices => new[] { "The cell's force field (above)" }
        .Concat(Field.Library.Select(e => e.Label))
        .Concat(FieldViewModel.Waters.Select(w => "Water model · " + w.Name)).ToList();
    private int PackFfChoiceOf(string key, string value)
    {
        if (key == "water")
        {
            var j = FieldViewModel.Waters.FindIndex(w => w.Id.Equals(value, StringComparison.OrdinalIgnoreCase) || w.Name.Equals(value, StringComparison.OrdinalIgnoreCase));
            return j < 0 ? 0 : 1 + Field.Library.Count + j;
        }
        var lib = Field.Library.ToList();
        var i = lib.FindIndex(e => e.Id.Equals(value, StringComparison.OrdinalIgnoreCase) || e.File == value);
        return i < 0 ? 0 : 1 + i;
    }
    /// <summary>Row `row`'s force field: 0 the cell's, then the library's, then the water models (written into its block).</summary>
    public void SetPackRowForceField(int row, int choice)
    {
        var lines = _packText.Split('\n').ToList();
        int seen = -1, head = -1, end = -1;
        for (int i = 0; i < lines.Count; ++i)
        {
            var w = lines[i].Split('#')[0].Trim().Split((char[]?)null, StringSplitOptions.RemoveEmptyEntries);
            if (w.Length == 0) continue;
            var k = w[0].ToLowerInvariant();
            if (head < 0 && k is "molecule" or "structure" && ++seen == row) head = i;
            else if (head >= 0 && k == "end" && !(w.Length > 1 && w[1].Equals("atoms", StringComparison.OrdinalIgnoreCase))) { end = i; break; }
        }
        if (head < 0 || end < 0) return;
        for (int i = end - 1; i > head; --i)
        {
            var k = lines[i].Split('#')[0].Trim().Split(' ', StringSplitOptions.RemoveEmptyEntries).FirstOrDefault()?.ToLowerInvariant();
            if (k is "forcefield" or "water") lines.RemoveAt(i);
        }
        var nl = Field.Library.Count;
        string? add = choice <= 0 ? null
            : choice <= nl ? $"  forcefield {Field.Library[choice - 1].Id}"
            : choice - 1 - nl < FieldViewModel.Waters.Count ? $"  water      {FieldViewModel.Waters[choice - 1 - nl].Id}" : null;
        if (add != null) lines.Insert(head + 1, add);
        PackText = string.Join('\n', lines);
    }

    /// <summary>Row `row`'s number of copies (its count / number line; a fixed molecule keeps one).</summary>
    public void SetPackRowCount(int row, int count)
    {
        if (count < 1) return;
        var lines = _packText.Split('\n').ToList();
        int seen = -1, head = -1;
        for (int i = 0; i < lines.Count; ++i)
        {
            var w = lines[i].Split('#')[0].Trim().Split((char[]?)null, StringSplitOptions.RemoveEmptyEntries);
            if (w.Length == 0) continue;
            var k = w[0].ToLowerInvariant();
            if (head < 0)
            {
                if (k is "molecule" or "structure" && ++seen == row) head = i;
                continue;
            }
            if (k is "count" or "number") { lines[i] = $"  {k,-7} {count}"; PackText = string.Join('\n', lines); return; }
            if (k == "fixed") return;
            if (k == "end" && !(w.Length > 1 && w[1].Equals("atoms", StringComparison.OrdinalIgnoreCase)))
            {
                lines.Insert(head + 1, IsCapsPack(_packText) ? $"  count   {count}" : $"  number {count}");
                PackText = string.Join('\n', lines);
                return;
            }
        }
    }

    /// <summary>After packing: the groups by the rows' own force fields (null when no row has one).</summary>
    /// <summary>The packed cell's components as groups (as the blend's): the structure packed around, then each molecule row,
    /// in that order — each with its row's own force field or the Pack's, so every component has its own atom types and the
    /// LAMMPS input names it (group NAME type …), the pairs between components by the force fields' own mixing rule.</summary>
    private string? PackGroupSpec(Interop.CapsDocument doc)
    {
        System.Text.Json.Nodes.JsonArray items;
        try { items = System.Text.Json.Nodes.JsonNode.Parse(doc.PackItemsJson())!.AsArray(); } catch { return null; }
        var rows = items.Where(i => ((string?)i!["molecules"] ?? "").Length > 0).ToList();
        if (rows.Count < 2 && !rows.Any(i => ((string?)i!["forcefield"] ?? "").Length > 0)) return null;   // one component: one force field
        var idx = PackFfIndex;
        var packFf = idx >= 0 && idx < Field.Library.Count ? Field.Library[idx].File : null;
        var groups = new System.Text.Json.Nodes.JsonArray();
        var charges = FieldViewModel.CoreChargesOf(_packCharges);
        var used = new HashSet<string>(StringComparer.OrdinalIgnoreCase);
        foreach (var i in rows)
        {
            var ff = (string?)i!["forcefield"] ?? "";
            // the name: a row's comment (the structure packed around: "TITLE, kept where it is"), else its file's name
            var raw = (string?)i["name"] ?? "group";
            var hash = raw.IndexOf('#');
            var nm = hash >= 0 && raw[(hash + 1)..].Trim() is { Length: > 0 } c ? c.Split(", kept")[0].Trim() : Path.GetFileNameWithoutExtension(raw.Trim());
            if (nm.StartsWith("host_", StringComparison.Ordinal)) nm = _packHostTitle.Length > 0 ? _packHostTitle : "host";   // the structure packed around
            if (nm.Length == 0) nm = "group";
            var unique = nm;
            for (var k = 2; !used.Add(unique); ++k) unique = $"{nm}_{k}";
            var g = new System.Text.Json.Nodes.JsonObject { ["name"] = unique, ["molecules"] = (string)i["molecules"]! };
            if (ff.StartsWith("water:", StringComparison.Ordinal)) g["water"] = ff[6..];
            else
            {
                var file = ff.Length > 0 ? Field.Library.FirstOrDefault(x => x.Id.Equals(ff, StringComparison.OrdinalIgnoreCase))?.File ?? ff : packFf;
                if (file == null) return null;
                g["forcefield"] = file;
                g["charges"] = charges;
            }
            groups.Add(g);
        }
        if (packFf != null && rows.Count < items.Count)   // anything not in a row (should not happen): the Pack's force field
            groups.Add(new System.Text.Json.Nodes.JsonObject { ["name"] = "rest", ["molecules"] = "rest", ["forcefield"] = packFf, ["charges"] = charges });
        return new System.Text.Json.Nodes.JsonObject { ["groups"] = groups }.ToJsonString();
    }

    // ---- molecules to pack from their SMILES: rubber curatives and additives, solvents, or any SMILES
    /// <summary>Curatives and additives of the fragment library (whole molecules): sulfur, accelerators, antidegradants …</summary>
    public List<FragmentItem> PackAdditives => Fragments.Where(f => f.Category == "Rubber additives" && !f.Smiles.Contains('*')
                                                                    && !f.Name.Contains("fragment", StringComparison.OrdinalIgnoreCase) && !f.Name.Contains(" unit", StringComparison.OrdinalIgnoreCase)).ToList();
    public List<FragmentItem> PackSolvents => Fragments.Where(f => f.Category == "Solvents" && !f.Smiles.Contains('*')).ToList();
    private string _packSmiles = "";
    public string PackSmiles { get => _packSmiles; set => Set(ref _packSmiles, value ?? ""); }
    /// <summary>Where molecules built for Pack are kept (their inputs name them).</summary>
    public static string PackMoleculeFolder => Path.Combine(AppSettings.Override != null ? Path.GetDirectoryName(AppSettings.Override)! : AppSettings.Folder, "molecules");

    /// <summary>A molecule from its SMILES (hydrogens added, 3D, UFF-cleaned), kept as a .mol2 (bond orders kept) and added as a row.</summary>
    public async Task AddPackMolecule(string smiles, string? name = null)
    {
        smiles = smiles.Trim();
        if (smiles.Length == 0) { Status = "Type a SMILES (S1SSSSSSS1 for sulfur S8) or pick a molecule"; return; }
        // a readable name: the library's for the same SMILES, else (after building) the formula — never the SMILES mangled
        if (string.IsNullOrWhiteSpace(name)) name = Fragments.FirstOrDefault(f => f.Smiles == smiles)?.Name;
        var dir = PackMoleculeFolder;
        Directory.CreateDirectory(dir);
        try
        {
            var tmp = Path.Combine(dir, $"_building_{Environment.ProcessId}.mol2");
            var (atoms, formula) = await Task.Run(() =>
            {
                var (doc, _) = Interop.CapsDocument.BuildSmiles(smiles, "uff", 1, 1, name ?? smiles);
                using (doc)
                {
                    doc.Save(tmp);
                    var n = (int)doc.Summary().Atoms;
                    var el = new SortedDictionary<string, int>(StringComparer.Ordinal);
                    for (var i = 0; i < n; ++i) { var e = doc.Atom(i).ElementSymbol; el[e] = el.GetValueOrDefault(e) + 1; }
                    string Part(string e) => el.TryGetValue(e, out var c) ? e + (c > 1 ? c.ToString(CultureInfo.InvariantCulture) : "") : "";
                    var f = el.ContainsKey("C") ? Part("C") + Part("H") + string.Concat(el.Keys.Where(k => k is not ("C" or "H")).Select(Part))
                                                : string.Concat(el.Keys.Select(Part));   // Hill order
                    return (n, f);
                }
            });
            name = string.IsNullOrWhiteSpace(name) ? formula : name!;
            var safe = System.Text.RegularExpressions.Regex.Replace(new string(name.Select(c => char.IsLetterOrDigit(c) || c is '-' or '_' ? c : '_').ToArray()), "_+", "_").Trim('_');
            if (safe.Length == 0) safe = "molecule";
            if (safe.Length > 40) safe = safe[..40];
            var path = Path.Combine(dir, safe + ".mol2");
            File.Move(tmp, path, true);
            AddPackStructure(path);
            Status = $"{name} ({formula}): built from {smiles} ({atoms} atoms) and added — set its count and region below";
        }
        catch (Exception e) { Status = $"Could not build {name ?? smiles}: {e.Message}"; }
    }

    /// <summary>Pack › Export: the packed cell to LAMMPS or GROMACS (the Export center, with its checks).</summary>
    public void PackExport() => OpenExportCenter();
}
