using System.Globalization;
using System.Text.Json.Nodes;
using CapsStudio.Interop;

namespace CapsStudio.ViewModels;

/// <summary>Biomolecule builder › Nucleic acid: a single DNA or RNA strand from its sequence, grown as a coil from
/// nucleotide units (D-sugars kept, phosphodiester links, 5′-OH and 3′-phosphate ends). No double helix: that needs
/// fibre coordinates CAPS does not ship.</summary>
public sealed partial class MainViewModel
{
    // head on O5′, tail on P: unit k's phosphate bonds to unit k+1's O5′ (O3′–P(=O)(O⁻)–O5′)
    private static readonly Dictionary<char, string> NaBases = new()
    {
        ['A'] = "n2cnc3c(N)ncnc32", ['C'] = "N2C=CC(N)=NC2=O", ['G'] = "N2C=NC3=C2N=C(N)NC3=O",
        ['T'] = "N2C=C(C)C(=O)NC2=O", ['U'] = "N2C=CC(=O)NC2=O",
    };
    private static string NaUnit(char b, bool rna) => rna
        ? $"*OC[C@H]1O[C@@H]({NaBases[b]})[C@H](O)[C@@H]1OP(=O)([O-])*"
        : $"*OC[C@H]1O[C@@H]({NaBases[b]})C[C@@H]1OP(=O)([O-])*";

    private bool _bioNucleic, _naRna;
    private CapsDocument? _naDoc;
    /// <summary>A copy of the last strand built, for the builder's preview.</summary>
    public CapsDocument? NaDoc { get => _naDoc; private set { var old = _naDoc; if (Set(ref _naDoc, value)) { old?.Dispose(); Raise(nameof(NaHud)); } } }
    public string NaHud => _naDoc == null ? "Build a strand to see it here" : $"{(_naRna ? "RNA" : "DNA")} · {_naDoc.Summary().Atoms:N0} atoms · single-stranded coil";
    private string _naSeq = "ACGTACGT";
    public bool BioNucleic { get => _bioNucleic; set { if (Set(ref _bioNucleic, value)) { Raise(nameof(BioPeptide)); Raise(nameof(BioBuildText)); Raise(nameof(NaInfo)); } } }
    public bool BioPeptide { get => !_bioNucleic; set => BioNucleic = !value; }
    public bool NaRna { get => _naRna; set { if (Set(ref _naRna, value)) { Raise(nameof(NaDna)); NaSequence = _naSeq; } } }
    public bool NaDna { get => !_naRna; set => NaRna = !value; }
    public string BioBuildText => _bioNucleic ? "Build strand" : "Build peptide";
    /// <summary>A, C, G and T (DNA) or U (RNA); spaces, digits and line breaks are ignored.</summary>
    public string NaSequence
    {
        get => _naSeq;
        set
        {
            var s = new string((value ?? "").ToUpperInvariant().Where(char.IsLetter).ToArray());
            s = _naRna ? s.Replace('T', 'U') : s.Replace('U', 'T');
            if (Set(ref _naSeq, s)) Raise(nameof(NaInfo));
            else Raise(nameof(NaInfo));
        }
    }
    public string NaInfo
    {
        get
        {
            var bad = _naSeq.Where(c => !(_naRna ? "ACGU" : "ACGT").Contains(c)).Distinct().ToArray();
            if (bad.Length > 0) return "Not a base here: " + string.Join(" ", bad);
            if (_naSeq.Length == 0) return "Type a sequence";
            var gc = _naSeq.Count(c => c is 'G' or 'C');
            return string.Format(CultureInfo.InvariantCulture, "{0} nt · GC {1:F0} % · net charge −{2} (one per phosphate; no counterions) · 5′-OH, 3′-phosphate",
                _naSeq.Length, 100.0 * gc / _naSeq.Length, _naSeq.Length);
        }
    }

    /// <summary>The strand as one chain of nucleotide units in sequence order, grown in a roomy cell and opened.</summary>
    public async Task BuildNucleic()
    {
        if (BioBuilding || _naSeq.Length == 0) return;
        var alphabet = _naRna ? "ACGU" : "ACGT";
        if (_naSeq.Any(c => !alphabet.Contains(c))) { BioError = NaInfo; return; }
        BioBuilding = true;
        var letters = _naSeq.Distinct().OrderBy(c => alphabet.IndexOf(c)).ToList();
        var spec = new JsonObject
        {
            ["units"] = new JsonArray(letters.Select(b => (JsonNode)new JsonObject { ["name"] = (_naRna ? "r" : "d") + b, ["smiles"] = NaUnit(b, _naRna) }).ToArray()),
            ["sequence"] = "pattern",
            ["pattern"] = new string(_naSeq.Select(c => (char)('A' + letters.IndexOf(c))).ToArray()),
            ["dp"] = _naSeq.Length,
            ["keep_configuration"] = 1,
        };
        var title = (_naRna ? "RNA " : "DNA ") + (_naSeq.Length <= 16 ? _naSeq : _naSeq[..16] + "…");
        try
        {
            var (doc, rep) = await Task.Run(() =>
            {
                var r = CapsDocument.GrowChains(spec.ToJsonString(), new CapsGrowOpts { Chains = 1, Dp = 0, Seed = 1, Density = 0.02, ContactScale = -0.9, Curve = 1 }, null, title);
                r.Doc.Edit("{\"op\":\"phosphate_ends\"}");   // the 3′ P–H cap → P–OH
                return r;
            });
            NaDoc = doc.Copy("strand preview");
            Show(doc, title);
            GrownUnsaved = true;
            BioLog = rep + "\nsingle strand, a coil (no double helix); D-sugars kept; 3′ end a phosphate monoester";
            Status = $"Built {title}: {doc.Summary().Atoms} atoms · a single-stranded coil";
            SetModule(8);
        }
        catch (Exception e) { BioError = e.Message; Status = "Could not build: " + e.Message; }
        finally { BioBuilding = false; }
    }
}
