using System.Collections.ObjectModel;
using System.Globalization;
using System.Text.Json.Nodes;
using CapsStudio.Interop;

namespace CapsStudio.ViewModels;

/// <summary>A run of one secondary structure for the HUD.</summary>
public sealed record BioRun(string Text, Avalonia.Media.IBrush Colour);

/// <summary>A residue of the sequence grid: its code and secondary structure (H E P C), selectable.</summary>
public sealed class ResidueCell : ObservableObject
{
    private char _ss;
    private bool _selected;
    public ResidueCell(int index, char code, char ss) { Index = index; Code = code.ToString(); _ss = ss; }
    public int Index { get; }
    public string Code { get; }
    public char Ss { get => _ss; set { if (Set(ref _ss, value)) { Raise(nameof(IsHelix)); Raise(nameof(IsStrand)); Raise(nameof(IsPpii)); Raise(nameof(IsCoil)); } } }
    public bool IsHelix => _ss == 'H';
    public bool IsStrand => _ss == 'E';
    public bool IsPpii => _ss == 'P';
    public bool IsCoil => _ss is not ('H' or 'E' or 'P');
    public bool Selected { get => _selected; set => Set(ref _selected, value); }
    public string Tip => $"{Code}{Index + 1}";
}

/// <summary>Biomolecule builder (design/boards/BioBuilder): a peptide from its sequence, a secondary structure per
/// residue (φ, ψ, ω), termini and protonation; all-atom with hydrogens, cleaned up with UFF.</summary>
public sealed partial class MainViewModel
{
    public bool IsBio => _module == 30;
    public ObservableCollection<ResidueCell> BioCells { get; } = new();
    public static readonly string[] BioNTerms = ["NH₃⁺", "NH₂", "Acetyl (ACE)"];
    public static readonly string[] BioCTerms = ["COO⁻", "COOH", "N-methyl amide (NME)"];
    private static readonly string[] NTermCodes = ["NH3+", "NH2", "ACE"], CTermCodes = ["COO-", "COOH", "NME"];
    public static readonly string[] BioProtonations = ["pH 7.0 defaults", "pH 4.0 (acidic)", "pH 10.0 (basic)", "Neutral side chains"];
    private static readonly double[] BioPhValues = [7.0, 4.0, 10.0, 7.0];
    public static readonly string[] BioCleanups = ["UFF clean-up", "None (placed geometry)"];

    private bool _bioLoaded;
    public void OpenBio()
    {
        if (!_bioLoaded)
        {
            _bioLoaded = true;
            _bioSeq = "AEAAAKEAAAKEAAAKAGGSPGGSPGSGSG";
            Raise(nameof(BioSequence));
            RebuildCells(new string('H', 22) + new string('C', 8));
        }
        SetModule(30);
        BioPreview();
    }

    // ---------------------------------------------------------------- sequence

    private string _bioSeq = "";
    public string BioSequence
    {
        get => _bioSeq;
        set
        {
            var clean = new string((value ?? "").Where(char.IsLetter).Select(char.ToUpperInvariant).ToArray());
            if (!Set(ref _bioSeq, clean)) return;
            RebuildCells(BioStructure);
            BioPreview();
        }
    }
    private string BioStructure => string.Concat(BioCells.Select(c => c.Ss));
    public string BioCountText => $"{_bioSeq.Length} residue{(_bioSeq.Length == 1 ? "" : "s")}";

    private void RebuildCells(string ss)
    {
        BioCells.Clear();
        for (var i = 0; i < _bioSeq.Length; i++) BioCells.Add(new ResidueCell(i, _bioSeq[i], i < ss.Length ? ss[i] : 'C'));
        _anchor = -1;
        Raise(nameof(BioCountText));
        RaiseBioSegments();
    }

    public void ImportFasta(string path)
    {
        try
        {
            var seq = CapsDocument.FastaSequence(File.ReadAllText(path));
            if (seq.Length == 0) { BioError = $"{Path.GetFileName(path)} holds no sequence"; return; }
            BioTitle = Path.GetFileNameWithoutExtension(path);
            BioSequence = seq;
            BioLog = $"Imported {seq.Length} residues from {Path.GetFileName(path)}";
        }
        catch (Exception e) { BioError = e.Message; }
    }

    // ---------------------------------------------------------------- selection and structure

    private int _anchor = -1;
    /// <summary>Click selects one residue; with Shift the range from the last click.</summary>
    public void SelectResidue(ResidueCell cell, bool extend)
    {
        if (!extend || _anchor < 0)
        {
            var only = cell.Selected && BioCells.Count(c => c.Selected) == 1;
            foreach (var c in BioCells) c.Selected = false;
            if (!only) { cell.Selected = true; _anchor = cell.Index; } else _anchor = -1;
        }
        else
        {
            int a = Math.Min(_anchor, cell.Index), b = Math.Max(_anchor, cell.Index);
            foreach (var c in BioCells) c.Selected = c.Index >= a && c.Index <= b;
        }
        Raise(nameof(BioSelectionText));
        var sel = BioCells.Where(c => c.Selected).ToList();
        if (sel.Count > 0 && sel.All(c => c.Ss == sel[0].Ss)) { _bioType = TypeOf(sel[0].Ss); RaiseBioType(); }
    }
    public string BioSelectionText
    {
        get
        {
            var sel = BioCells.Where(c => c.Selected).ToList();
            return sel.Count == 0 ? "No residues selected: a structure applies to all" : sel.Count == 1 ? $"Residue {sel[0].Tip} selected" : $"Residues {sel[0].Index + 1}–{sel[^1].Index + 1} selected";
        }
    }

    private static int TypeOf(char ss) => ss switch { 'H' => 0, 'E' => 1, 'P' => 2, _ => 3 };
    private static readonly char[] TypeCodes = ['H', 'E', 'P', 'C'];
    // the structure the φ ψ ω fields show; choosing one applies it to the selection (or to every residue)
    private int _bioType;
    public int BioType
    {
        get => _bioType;
        set
        {
            _bioType = Math.Clamp(value, 0, 3);
            var sel = BioCells.Where(c => c.Selected).ToList();
            foreach (var c in sel.Count > 0 ? sel : BioCells.ToList()) c.Ss = TypeCodes[_bioType];
            RaiseBioType();
            RaiseBioSegments();
            BioPreview();
        }
    }
    public bool BioHelix { get => _bioType == 0; set { if (value) BioType = 0; } }
    public bool BioStrand { get => _bioType == 1; set { if (value) BioType = 1; } }
    public bool BioPpii { get => _bioType == 2; set { if (value) BioType = 2; } }
    public bool BioCoil { get => _bioType == 3; set { if (value) BioType = 3; } }

    private decimal[] _helixT = [-57, -47, 180], _strandT = [-139, 135, 180], _ppiiT = [-75, 145, 180];
    private decimal[]? Torsions => _bioType switch { 0 => _helixT, 1 => _strandT, 2 => _ppiiT, _ => null };
    public bool BioTorsionsEditable => _bioType != 3;
    public decimal BioPhi { get => Torsions?[0] ?? -70; set => SetTorsion(0, value); }
    public decimal BioPsi { get => Torsions?[1] ?? 140; set => SetTorsion(1, value); }
    public decimal BioOmega { get => Torsions?[2] ?? 180; set => SetTorsion(2, value); }
    private void SetTorsion(int k, decimal v)
    {
        var t = Torsions;
        if (t == null || t[k] == v) return;
        t[k] = Math.Clamp(v, -180, 180);
        RaiseBioType();
        BioPreview();
    }
    public string BioCite => _bioType switch
    {
        0 => "Ideal α-helix φ/ψ — Pauling, Corey & Branson, PNAS 37, 205 (1951)",
        1 => "Antiparallel β-sheet φ/ψ — Pauling & Corey, PNAS 37, 729 (1951)",
        2 => "Polyproline II φ/ψ — Adzhubei, Sternberg & Makarov, J. Mol. Biol. 425, 2100 (2013)",
        _ => "Coil: φ/ψ drawn per residue from the αR, β and PPII basins (seeded); proline keeps φ −65°",
    };
    private void RaiseBioType()
    {
        foreach (var n in new[] { nameof(BioHelix), nameof(BioStrand), nameof(BioPpii), nameof(BioCoil), nameof(BioPhi), nameof(BioPsi), nameof(BioOmega), nameof(BioTorsionsEditable), nameof(BioCite) })
            Raise(n);
    }

    /// <summary>Runs of one structure: "α-helix 1–22", "coil 23–30", with the structure's colour.</summary>
    public List<BioRun> BioRuns
    {
        get
        {
            var runs = new List<BioRun>();
            var ss = BioStructure;
            for (var i = 0; i < ss.Length;)
            {
                var j = i;
                while (j + 1 < ss.Length && ss[j + 1] == ss[i]) j++;
                var (name, colour) = ss[i] switch { 'H' => ("α-helix", "#E07A5F"), 'E' => ("β-strand", "#6FA8DC"), 'P' => ("PPII", "#9B7BD6"), _ => ("coil", "#8A9097") };
                runs.Add(new BioRun(i == j ? $"{name} {i + 1}" : $"{name} {i + 1}–{j + 1}", Avalonia.Media.Brush.Parse(colour)));
                i = j + 1;
            }
            return runs;
        }
    }
    public List<string> BioSegments => BioRuns.Select(r => r.Text).ToList();
    /// <summary>The HUD: the first four runs and how many more.</summary>
    public List<BioRun> BioHudRuns
    {
        get
        {
            var all = BioRuns;
            if (all.Count <= 4) return all;
            var head = all.Take(3).ToList();
            head.Add(new BioRun($"+{all.Count - 3} more", Avalonia.Media.Brushes.Transparent));
            return head;
        }
    }
    public string BioSegmentsText => string.Join(" · ", BioSegments);
    private void RaiseBioSegments() { Raise(nameof(BioRuns)); Raise(nameof(BioSegments)); Raise(nameof(BioHudRuns)); Raise(nameof(BioSegmentsText)); }

    // ---------------------------------------------------------------- termini and state

    private int _nTerm, _cTerm, _prot, _cleanup;
    public int BioNTerm { get => _nTerm; set { if (Set(ref _nTerm, Math.Clamp(value, 0, 2))) BioPreview(); } }
    public int BioCTerm { get => _cTerm; set { if (Set(ref _cTerm, Math.Clamp(value, 0, 2))) BioPreview(); } }
    public int BioProtonation { get => _prot; set { if (Set(ref _prot, Math.Clamp(value, 0, 3))) BioPreview(); } }
    public int BioCleanup { get => _cleanup; set => Set(ref _cleanup, Math.Clamp(value, 0, 1)); }
    private string _bioTitle = "peptide";
    public string BioTitle { get => _bioTitle; set => Set(ref _bioTitle, value); }

    private string BioOptions(bool preview) => new JsonObject
    {
        ["sequence"] = _bioSeq, ["structure"] = BioStructure,
        ["helix"] = new JsonArray(_helixT.Select(v => (JsonNode)(double)v).ToArray()),
        ["strand"] = new JsonArray(_strandT.Select(v => (JsonNode)(double)v).ToArray()),
        ["ppii"] = new JsonArray(_ppiiT.Select(v => (JsonNode)(double)v).ToArray()),
        ["n_term"] = NTermCodes[_nTerm], ["c_term"] = CTermCodes[_cTerm],
        ["ph"] = BioPhValues[_prot], ["neutral"] = _prot == 3,
        ["cleanup"] = !preview && _cleanup == 0, ["ribbon"] = preview, ["seed"] = 1,
    }.ToJsonString();

    // ---------------------------------------------------------------- preview and build

    private CapsDocument? _bioDoc;
    public CapsDocument? BioDoc { get => _bioDoc; private set { if (Set(ref _bioDoc, value)) BioViewChanged?.Invoke(); } }
    public event Action? BioViewChanged;
    private int _bioTicket;
    private string _bioLog = "", _bioError = "", _bioSummary = "";
    public string BioLog { get => _bioLog; private set => Set(ref _bioLog, value); }
    public string BioError { get => _bioError; private set { if (Set(ref _bioError, value)) Raise(nameof(BioHasError)); } }
    public bool BioHasError => _bioError.Length > 0;
    public string BioSummary { get => _bioSummary; private set => Set(ref _bioSummary, value); }
    private bool _bioBuilding;
    public bool BioBuilding { get => _bioBuilding; private set { if (Set(ref _bioBuilding, value)) Raise(nameof(BioIdle)); } }
    public bool BioIdle => !_bioBuilding;
    public string BioHudTitle => $"Peptide · {BioCountText}";

    public void BioPreview()
    {
        Raise(nameof(BioHudTitle));
        if (_module != 30) return;
        if (_bioSeq.Length == 0) { BioError = "Type a sequence (one-letter codes)"; return; }
        var ticket = ++_bioTicket;
        var opts = BioOptions(true);
        var title = _bioTitle;
        Task.Run(() =>
        {
            try
            {
                var (d, r) = CapsDocument.PeptideBuild(opts, title);
                return (Doc: (CapsDocument?)d, Report: r, Error: (string?)null);
            }
            catch (Exception e) { return (Doc: (CapsDocument?)null, Report: "", Error: (string?)e.Message); }
        }).ContinueWith(t => Avalonia.Threading.Dispatcher.UIThread.Post(() =>
        {
            var (doc, rep, err) = t.Result;
            if (ticket != _bioTicket) { doc?.Dispose(); return; }
            if (err != null || doc == null) { BioError = err ?? "cannot build"; return; }
            BioError = "";
            var old = _bioDoc;
            BioDoc = doc;
            old?.Dispose();
            BioSummary = rep.Split('\n').FirstOrDefault() ?? "";
        }));
    }

    /// <summary>Builds the peptide (with the clean-up) and opens it as the Studio document.</summary>
    public async Task BuildPeptide()
    {
        if (BioBuilding || _bioSeq.Length == 0) return;
        BioBuilding = true;
        var opts = BioOptions(false);
        var title = _bioTitle;
        try
        {
            var (doc, rep) = await Task.Run(() => CapsDocument.PeptideBuild(opts, title));
            Show(doc, title);
            Record($"doc = caps.build.peptide({PyStr(_bioSeq)}, structure={PyStr(BioStructure)}, n_term={PyStr(NTermCodes[_nTerm])}, c_term={PyStr(CTermCodes[_cTerm])}, ph={BioPhValues[_prot].ToString(CultureInfo.InvariantCulture)}, cleanup={(_cleanup == 0 ? "True" : "False")})");
            GrownUnsaved = true;
            BioLog = rep;
            Status = "Built · " + (rep.Split('\n').FirstOrDefault() ?? "");
            SetModule(8);
        }
        catch (Exception e) { BioError = e.Message; Status = "Could not build: " + e.Message; }
        finally { BioBuilding = false; }
    }
}
