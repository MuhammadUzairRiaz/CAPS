using System.Collections.ObjectModel;
using System.Globalization;
using CapsStudio.Interop;

namespace CapsStudio.ViewModels;

public sealed record FocusRow(string Id, string Name, string Molecule, string Distance);

/// <summary>Keyboard and screen reader (design/boards/VisAccess): the 3D view is one focusable region; keys walk the
/// structure the way a reader walks text — atom, bond, molecule, frame — and each step is announced (a live region
/// for screen readers, and a bar in the view for sighted keyboard users).</summary>
public sealed partial class MainViewModel
{
    private int _focusAtom = -1, _focusPrev = -1;
    private int[]? _molOf;           // molecule of each atom (cached per document and topology)
    private CapsDocument? _molOfDoc;
    private long _molOfBonds = -1;
    private string _announcement = "";

    public static readonly string[] Verbosities = ["Brief · name and molecule", "Full · name, place, charge, nearest"];

    public int FocusAtom
    {
        get => _focusAtom;
        private set
        {
            if (!Set(ref _focusAtom, value)) return;
            Raise(nameof(HasFocusAtom));
            Raise(nameof(ShowAnnouncementBar));
            Raise(nameof(FocusTitle));
            RefreshFocusRows();
            RenderRequested?.Invoke();
        }
    }
    public bool HasFocusAtom => _focusAtom >= 0 && _doc != null;
    public string Announcement { get => _announcement; private set => Set(ref _announcement, value); }
    public bool ShowAnnouncementBar => HasFocusAtom && _settings.ShowAnnouncement;
    public ObservableCollection<FocusRow> FocusRows { get; } = new();
    public string FocusTitle => HasFocusAtom ? $"focus: atom {SafeAtom(_focusAtom)?.Id} · molecule {SafeAtom(_focusAtom)?.Mol}" : "no focus · click the view, then ↑ ↓";

    public int ReaderVerbosity
    {
        get => _settings.ReaderVerbosity;
        set { if (_settings.ReaderVerbosity == value) return; _settings.ReaderVerbosity = Math.Clamp(value, 0, 1); Raise(); Changed("Reader verbosity"); if (HasFocusAtom) AnnounceFocus(); }
    }
    public bool SelectionSound { get => _settings.SelectionSound; set { if (_settings.SelectionSound == value) return; _settings.SelectionSound = value; Raise(); Changed("Sound on selection"); if (value) SoundCue.Selection(); } }
    public bool AnnounceFrames { get => _settings.AnnounceFrames; set { if (_settings.AnnounceFrames == value) return; _settings.AnnounceFrames = value; Raise(); Changed("Announce frame changes"); } }
    public bool ShowAnnouncement
    {
        get => _settings.ShowAnnouncement;
        set { if (_settings.ShowAnnouncement == value) return; _settings.ShowAnnouncement = value; Raise(); Raise(nameof(ShowAnnouncementBar)); Changed("Show announcements"); }
    }

    private CapsAtomInfo? SafeAtom(int i)
    {
        if (_doc == null || i < 0) return null;
        try { return _doc.Atom(i); } catch { return null; }
    }

    private int[] MolOf()
    {
        var s = _doc!.Summary();
        if (_molOf == null || _molOfDoc != _doc || _molOfBonds != s.Bonds || _molOf.Length != s.Atoms)
        {
            _molOf = _doc.MoleculeIndex((int)s.Atoms).Mol;
            _molOfDoc = _doc;
            _molOfBonds = s.Bonds;
        }
        return _molOf;
    }

    /// <summary>Where a walk starts: the focused atom, else the picked one, else the first.</summary>
    private int FocusStart() => _focusAtom >= 0 ? _focusAtom : Picked >= 0 ? Picked : 0;

    private void MoveFocus(int i)
    {
        if (_doc == null) return;
        if (i != _focusAtom) _focusPrev = _focusAtom;
        FocusAtom = i;
        AnnounceFocus();
    }

    /// <summary>Puts the focus ring on atom i (tests, screenshots, a picked atom).</summary>
    public void FocusOn(int i) => MoveFocus(i);

    /// <summary>↑ ↓: the previous or next atom in the focused atom's molecule.</summary>
    public void FocusStep(int d)
    {
        if (_doc == null) return;
        if (_focusAtom < 0) { MoveFocus(FocusStart()); return; }
        var mol = MolOf();
        var m = mol[_focusAtom];
        var i = _focusAtom;
        for (int k = 0; k < mol.Length; ++k)
        {
            i = ((i + d) % mol.Length + mol.Length) % mol.Length;
            if (mol[i] == m) break;
        }
        MoveFocus(i);
    }

    /// <summary>[ ]: the first atom of the previous or next molecule.</summary>
    public void FocusMolecule(int d)
    {
        if (_doc == null) return;
        var mol = MolOf();
        var count = mol.Length == 0 ? 0 : mol.Max() + 1;
        if (count == 0) return;
        var cur = _focusAtom >= 0 ? mol[_focusAtom] : (d > 0 ? -1 : count);
        var target = ((cur + d) % count + count) % count;
        MoveFocus(Array.IndexOf(mol, target));
    }

    /// <summary>B: along a bond — onward (not back where it came from), heavy atoms first; at a dead end, back.</summary>
    public void FocusBond()
    {
        if (_doc == null) return;
        if (_focusAtom < 0) { MoveFocus(FocusStart()); return; }
        var bonded = _doc.Bonded(_focusAtom);
        if (bonded.Length == 0) { Announcement = $"Atom {SafeAtom(_focusAtom)?.Id} has no bonds."; return; }
        var onward = bonded.Where(b => b != _focusPrev)
            .OrderBy(b => SafeAtom(b)?.ElementSymbol == "H" ? 1 : 0).ThenBy(b => b).ToList();
        MoveFocus(onward.Count > 0 ? onward[0] : bonded[0]);
    }

    /// <summary>Space: select the focused atom (replacing the selection).</summary>
    public void FocusSelect()
    {
        if (!HasFocusAtom) return;
        Pick(_focusAtom);
        Announcement = "Selected. " + Describe(_focusAtom, 0);
    }

    /// <summary>M: add the focused atom to the selection, measuring from the atoms selected before it.</summary>
    public void FocusMeasure()
    {
        if (!HasFocusAtom) return;
        Pick(_focusAtom, true);
        Announcement = HasMeasure ? MeasureText.Replace("\n", ". ") + "." : "Selected. " + Describe(_focusAtom, 0);
    }

    /// <summary>⌘⇧A: the selection read out.</summary>
    public void AnnounceSelection()
    {
        if (_doc == null) return;
        if (_selection.Count == 0) { Announcement = "Nothing selected."; return; }
        var parts = _selection.Select(i => Describe(i, 0)).ToList();
        Announcement = $"{_selection.Count} selected: " + string.Join(" ", parts) + (HasMeasure ? " " + MeasureText.Replace("\n", ". ") + "." : "");
    }

    public void ClearFocus()
    {
        _focusPrev = -1;
        FocusAtom = -1;
        Announcement = "";
    }

    private void AnnounceFocus() => Announcement = HasFocusAtom ? Describe(_focusAtom, _settings.ReaderVerbosity) : "";

    /// <summary>Called when the frame changes: the announcement follows the focused atom.</summary>
    private void FocusOnFrame()
    {
        if (!HasFocusAtom) return;
        RefreshFocusRows();
        if (_settings.AnnounceFrames) Announcement = $"Frame {_frame + 1} of {_frames}. " + Describe(_focusAtom, _settings.ReaderVerbosity);
    }

    private static readonly Dictionary<string, string> ElementNames = new()
    {
        ["H"] = "hydrogen", ["C"] = "carbon", ["N"] = "nitrogen", ["O"] = "oxygen", ["F"] = "fluorine", ["Si"] = "silicon", ["P"] = "phosphorus",
        ["S"] = "sulfur", ["Cl"] = "chlorine", ["Br"] = "bromine", ["I"] = "iodine", ["Na"] = "sodium", ["K"] = "potassium", ["Zn"] = "zinc",
        ["Ca"] = "calcium", ["Mg"] = "magnesium", ["Al"] = "aluminium", ["Ti"] = "titanium", ["Fe"] = "iron", ["Cu"] = "copper", ["B"] = "boron",
    };

    /// <summary>The words for one atom: brief (0) — name and molecule; full (1) — also what it is bonded to, charge, place, nearest.</summary>
    public string Describe(int i, int verbosity)
    {
        var a = SafeAtom(i);
        if (a is not { } at) return "";
        var inv = CultureInfo.InvariantCulture;
        var sym = at.ElementSymbol;
        var elem = ElementNames.TryGetValue(sym, out var nm) ? nm : sym;
        var name = at.Name.Length > 0 && at.Name != sym ? $"{at.Name}, {elem}" : elem;
        if (verbosity == 0) return string.Format(inv, "Atom {0}, {1}, molecule {2}.", at.Id, name, at.Mol);
        var bonded = _doc!.Bonded(i).Select(b => SafeAtom(b)?.ElementSymbol ?? "?").OrderBy(x => x == "H" ? 1 : 0).ThenBy(x => x).ToList();
        var bondText = bonded.Count == 0 ? "no bonds" : "bonded to " + string.Join(", ", bonded);
        var sb = new System.Text.StringBuilder();
        sb.Append(string.Format(inv, "Atom {0}, {1}, {2}. Molecule {3}. ", at.Id, name, bondText, at.Mol));
        if (Math.Abs(at.Charge) > 0) sb.Append(string.Format(inv, "Charge {0:+0.0000;−0.0000;0} e. ", at.Charge));
        sb.Append(string.Format(inv, "Position {0:F2}, {1:F2}, {2:F2} Å.", at.X, at.Y, at.Z));
        var nb = _doc.Neighbours(i, 1);
        if (nb.Length > 0 && SafeAtom(nb[0].Index) is { } n) sb.Append(string.Format(inv, " Nearest: atom {0} at {1:F2} Å.", n.Id, nb[0].Distance));
        return sb.ToString();
    }

    private void RefreshFocusRows()
    {
        FocusRows.Clear();
        if (!HasFocusAtom) return;
        var inv = CultureInfo.InvariantCulture;
        foreach (var (idx, d) in _doc!.Neighbours(_focusAtom, 4))
            if (SafeAtom(idx) is { } n)
                FocusRows.Add(new(n.Id.ToString(inv), n.Name.Length > 0 ? n.Name : n.ElementSymbol, n.Mol.ToString(inv), d.ToString("F2", inv)));
    }
}
