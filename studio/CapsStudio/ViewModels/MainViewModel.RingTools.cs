namespace CapsStudio.ViewModels;

// The ring tool's panel (sketcher): rings of 3 to 8 fused onto a picked bond or spiro on a picked atom, an atom put into a
// picked bond, another structure file joined at a picked leaving atom (core ring_edit.hpp through caps_edit).
public partial class MainViewModel
{
    public int[] RingSizes { get; } = [3, 4, 5, 6, 7, 8];
    private int _ringSize = 6;
    public int RingSize { get => _ringSize; set { if (Set(ref _ringSize, value)) Raise(nameof(RingMembersHint)); } }
    private string _ringMembers = "";
    /// <summary>The new ring atoms in order, e.g. "O" or "N C=O C C"; empty: carbons.</summary>
    public string RingMembers { get => _ringMembers; set => Set(ref _ringMembers, value ?? ""); }
    private bool _ringAromatic;
    public bool RingAromatic { get => _ringAromatic; set => Set(ref _ringAromatic, value); }
    public string RingMembersHint => $"Fused: {_ringSize - 2} new atom(s), from the second picked atom round · spiro: {_ringSize - 1}. Element symbols, C=O for a carbonyl; empty for carbons.";

    public void RingFusePicked()
    {
        if (_selection.Count != 2) { Status = "Pick the two atoms of a bond (⇧ click), then Fuse"; return; }
        RunEdit(new { op = "fuse_ring", i = _selection[0], j = _selection[1], size = _ringSize, members = _ringMembers, aromatic = _ringAromatic });
    }
    public void RingSpiroPicked()
    {
        if (_selection.Count != 1) { Status = "Pick the atom the two rings share, then Spiro"; return; }
        RunEdit(new { op = "spiro_ring", atom = _selection[0], size = _ringSize, members = _ringMembers });
    }
    public void InsertAtomPicked()
    {
        if (_selection.Count != 2) { Status = "Pick the two atoms of a bond (⇧ click), then Insert"; return; }
        RunEdit(new { op = "insert_atom", i = _selection[0], j = _selection[1], element = _buildElement });
    }

    private string _joinPieceAtom = "";
    /// <summary>The piece's leaving atom, numbered from 1; empty: CAPS picks the first hydrogen of its first heavy atom.</summary>
    public string JoinPieceAtom { get => _joinPieceAtom; set => Set(ref _joinPieceAtom, value ?? ""); }
    private string _joinTorsion = "";
    /// <summary>Degrees about the new bond; empty: the turn that keeps the piece clearest of the structure.</summary>
    public string JoinTorsion { get => _joinTorsion; set => Set(ref _joinTorsion, value ?? ""); }
    public void JoinFile(string path)
    {
        if (_selection.Count != 1) { Status = "Pick the leaving atom here (a hydrogen, or any atom with one bond), then Join"; return; }
        var inv = System.Globalization.CultureInfo.InvariantCulture;
        var pieceAtom = int.TryParse(_joinPieceAtom.Trim(), out var pa) && pa >= 1 ? pa - 1 : -1;
        if (double.TryParse(_joinTorsion.Replace("°", "").Trim(), System.Globalization.NumberStyles.Float, inv, out var tor))
            RunEdit(new { op = "join", file = path, at = _selection[0], piece_atom = pieceAtom, torsion = tor });
        else RunEdit(new { op = "join", file = path, at = _selection[0], piece_atom = pieceAtom });
    }

    public long AtomCountForShot() => _doc == null ? 0 : (long)_doc.Summary().Atoms;

    /// <summary>Headless check (Screenshot ringop=OP|ATOMS|SIZE|MEMBERS): the atoms picked, the tool run.</summary>
    public void RingOpForShot(string op, int[] atoms, int size, string members)
    {
        for (var k = 0; k < atoms.Length; ++k) Pick(atoms[k], k > 0);
        RingSize = size;
        if (op == "join") { JoinFile(members); return; }
        RingMembers = members;
        if (op == "fuse") RingFusePicked();
        else if (op == "aromatic") { RingAromatic = true; RingFusePicked(); }
        else if (op == "spiro") RingSpiroPicked();
        else if (op == "insert") { BuildElement = members.Length > 0 ? members : "C"; InsertAtomPicked(); }
    }
}
