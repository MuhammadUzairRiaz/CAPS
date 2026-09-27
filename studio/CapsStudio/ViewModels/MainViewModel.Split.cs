using System.Collections.ObjectModel;
using System.Globalization;
using System.Text.Json.Nodes;
using CapsStudio.Interop;

namespace CapsStudio.ViewModels;

/// <summary>A row of the split view's comparison: property, the left document's value, the right one's.</summary>
public sealed record SplitRow(string Property, string A, string B);

/// <summary>Split view (design/boards/LightSplit): the open document beside a second one, cameras synced, selections
/// and frames linked by index, their structure compared; the Paper (light) theme one click away.</summary>
public sealed partial class MainViewModel
{
    public bool IsSplit => _module == 34;
    public ObservableCollection<SplitRow> SplitRows { get; } = new();
    public event Action? SplitChanged;

    private CapsDocument? _splitB;
    private string _splitPathB = "";
    private bool _syncCamera = true, _syncSelection, _syncFrame = true;

    public CapsDocument? SplitDocB => _splitB;
    public string SplitTitleA => _doc != null ? Path.GetFileName(_doc.Path) : "No document";
    public string SplitTitleB => _splitB != null ? Path.GetFileName(_splitPathB) : "Open a second structure";
    public bool SplitHasB => _splitB != null;
    public bool SyncCamera { get => _syncCamera; set { if (Set(ref _syncCamera, value)) { Raise(nameof(SplitStatus)); SplitChanged?.Invoke(); } } }
    public bool SyncSelection { get => _syncSelection; set { if (Set(ref _syncSelection, value)) SplitChanged?.Invoke(); } }
    public bool SyncFrame { get => _syncFrame; set { if (Set(ref _syncFrame, value)) { if (value) SplitFrameB(); SplitChanged?.Invoke(); } } }
    public string SplitStatus => $"{(_splitB != null ? 2 : 1)} document{(_splitB != null ? "s" : "")}{(_syncCamera ? " · camera synced" : "")}";
    public string SplitThemeText => ThemeLight ? "Theme: Paper (light)" : ThemeDark ? "Theme: Graphite (dark)" : "Theme: follows the system";
    private string _splitTagA = "", _splitTagB = "", _splitError = "";
    public string SplitTagA { get => _splitTagA; private set => Set(ref _splitTagA, value); }
    public string SplitTagB { get => _splitTagB; private set => Set(ref _splitTagB, value); }
    public string SplitError { get => _splitError; private set => Set(ref _splitError, value); }
    /// <summary>Recent structures to put on the right (the Start page's list).</summary>
    public IEnumerable<RecentItem> SplitRecent => Recent.Where(r => _doc == null || r.Path != _doc.Path);

    public void OpenSplit()
    {
        if (_doc == null) return;
        SetModule(34);
        foreach (var n in new[] { nameof(SplitTitleA), nameof(SplitRecent), nameof(SplitThemeText), nameof(SplitStatus) }) Raise(n);
        RefreshSplit();
    }

    public async Task SetSplitB(string path, string? topology = null)
    {
        try
        {
            var topo = topology ?? TopologyFor(path);
            var doc = await Task.Run(() => CapsDocument.Open(path, topo.Length > 0 ? topo : null));
            _splitB?.Dispose();
            _splitB = doc;
            _splitPathB = path;
            SplitError = "";
            if (_syncFrame) SplitFrameB();
            foreach (var n in new[] { nameof(SplitDocB), nameof(SplitTitleB), nameof(SplitHasB), nameof(SplitStatus) }) Raise(n);
            RefreshSplit();
        }
        catch (Exception e) { SplitError = $"Could not open {Path.GetFileName(path)}: {e.Message}"; }
    }

    /// <summary>A document already open (a state of the left one) on the right.</summary>
    public void SetSplitBDocument(CapsDocument doc, string label)
    {
        _splitB?.Dispose();
        _splitB = doc;
        _splitPathB = label;
        SplitError = "";
        foreach (var n in new[] { nameof(SplitDocB), nameof(SplitTitleB), nameof(SplitHasB), nameof(SplitStatus) }) Raise(n);
        RefreshSplit();
    }

    /// <summary>The right document at the left one's frame (when it has that many).</summary>
    public void SplitFrameB()
    {
        if (_splitB == null) return;
        var n = _splitB.Summary().Frames;
        if (n > 1) { try { _splitB.SetFrame((int)Math.Min(_frame, n - 1)); } catch { } }
        SplitChanged?.Invoke();
    }

    public void ToggleSplitTheme()
    {
        SetTheme = ThemeLight ? "dark" : "light";
        RaiseTheme();
        Raise(nameof(SplitThemeText));
        SplitChanged?.Invoke();
    }

    private static (string Tag, List<string> Values) Describe(CapsDocument d)
    {
        var s = d.Summary();
        var tac = JsonNode.Parse(d.Tacticity())!;
        var label = tac["label"]!.GetValue<string>();
        var mols = d.Molecules();
        var rg = mols.Length > 0 ? mols.Average(m => m.Rg) : 0;
        var inv = CultureInfo.InvariantCulture;
        return (label.Length > 0 ? label : $"{s.Molecules:N0} molecules", new List<string>
        {
            s.Atoms.ToString("N0", inv), s.Molecules.ToString("N0", inv), s.Frames.ToString("N0", inv),
            s.CellValid != 0 ? s.Density.ToString("0.000", inv) + " g/cm³" : "no cell",
            tac["centres"]!.GetValue<double>().ToString("0", inv), label.Length > 0 ? label : "—",
            rg.ToString("0.00", inv) + " Å",
        });
    }

    public void RefreshSplit()
    {
        SplitRows.Clear();
        if (_doc == null || !IsSplit) return;
        try
        {
            var names = new[] { "Atoms", "Molecules", "Frames", "Density", "Stereocentres", "Tacticity", "Rg (mean molecule)" };
            var (tagA, a) = Describe(_doc);
            SplitTagA = tagA;
            List<string>? b = null;
            if (_splitB != null) { var (tagB, vb) = Describe(_splitB); SplitTagB = tagB; b = vb; }
            for (var k = 0; k < names.Length; k++) SplitRows.Add(new SplitRow(names[k], a[k], b?[k] ?? "—"));
        }
        catch (Exception e) { SplitError = e.Message; }
        Raise(nameof(SplitTitleA));
        SplitChanged?.Invoke();
    }

    /// <summary>The Studio's picked atoms (drawn in both panes when the selection is synced).</summary>
    public int[] PickedAtoms => _selection.ToArray();
    public void RefreshSplitSelection() => SplitChanged?.Invoke();

    public void CloseSplitB()
    {
        _splitB?.Dispose();
        _splitB = null;
        _splitPathB = "";
        foreach (var n in new[] { nameof(SplitDocB), nameof(SplitTitleB), nameof(SplitHasB), nameof(SplitStatus) }) Raise(n);
        RefreshSplit();
    }
}
