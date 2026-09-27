using System;
using System.Globalization;
using System.Linq;
using System.Text.Json.Nodes;

namespace CapsStudio.ViewModels;

public sealed record PeriodicRow(string What, string Wrapped, string MinImage, string Whole);

/// <summary>Studio › Periodic box (design/boards/PeriodicBox): the cell wrapped or with molecules whole, periodic images
/// around it, how many molecules cross a face and how many pieces wrapping makes, one chain measured wrapped, by minimum
/// image and whole, wrap-on-save, and centring on the selection.</summary>
public partial class MainViewModel
{
    public bool IsPeriodic => _module == 51;
    private int _pbShow, _pbSaveWrap, _pbImgA = 3, _pbImgB = 3, _pbImgC = 1, _pbFade = 30;
    private decimal _pbChain = 1;
    private string _pbChains = "—", _pbCross = "—", _pbPieces = "—", _pbBox = "—", _pbError = "", _pbNote = "";
    /// <summary>0 wrapped (atoms folded into the cell), 1 whole molecules, 2 whole with images around the cell.</summary>
    public int PbShow { get => _pbShow; set { if (Set(ref _pbShow, Math.Clamp(value, 0, 2))) ApplyPeriodicView(); } }
    public decimal PbImagesA { get => _pbImgA; set { if (Set(ref _pbImgA, (int)Math.Clamp(value, 1, 5))) ApplyPeriodicView(); } }
    public decimal PbImagesB { get => _pbImgB; set { if (Set(ref _pbImgB, (int)Math.Clamp(value, 1, 5))) ApplyPeriodicView(); } }
    public decimal PbImagesC { get => _pbImgC; set { if (Set(ref _pbImgC, (int)Math.Clamp(value, 1, 5))) ApplyPeriodicView(); } }
    public decimal PbFade { get => _pbFade; set { if (Set(ref _pbFade, (int)Math.Clamp(value, 0, 90))) ApplyPeriodicView(); } }
    public int PbSaveWrap { get => _pbSaveWrap; set { if (Set(ref _pbSaveWrap, Math.Clamp(value, 0, 2))) { try { _doc?.SetSaveWrap(value == 0 ? 0 : value == 1 ? 1 : 2); } catch { } } } }
    public decimal PbChain { get => _pbChain; set { if (Set(ref _pbChain, Math.Max(1, value))) RefreshPeriodic(); } }
    public string PbChains { get => _pbChains; private set => Set(ref _pbChains, value); }
    public string PbCross { get => _pbCross; private set => Set(ref _pbCross, value); }
    public string PbPieces { get => _pbPieces; private set => Set(ref _pbPieces, value); }
    public string PbBox { get => _pbBox; private set => Set(ref _pbBox, value); }
    public string PbError { get => _pbError; private set { Set(ref _pbError, value); Raise(nameof(PbHasError)); } }
    public bool PbHasError => _pbError.Length > 0;
    public string PbNote { get => _pbNote; private set => Set(ref _pbNote, value); }
    public System.Collections.ObjectModel.ObservableCollection<PeriodicRow> PbRows { get; } = new();
    public event Action? PeriodicChanged;

    public void OpenPeriodic()
    {
        SetModule(51);
        _pbShow = _wrap ? 0 : 1;
        Raise(nameof(PbShow));
        ApplyPeriodicView();
        RefreshPeriodic();
    }

    private void ApplyPeriodicView()
    {
        if (_doc == null) return;
        Wrap = _pbShow == 0;
        try
        {
            if (_pbShow == 2) _doc.SetImages(_pbImgA, _pbImgB, _pbImgC, 1 - _pbFade / 100.0);
            else _doc.SetImages(1, 1, 1, 0);
        }
        catch { /* an older core */ }
        RenderRequested?.Invoke();
        PeriodicChanged?.Invoke();
    }

    /// <summary>Images are a view of this page; they go when the page does.</summary>
    // during a run the document is held: the images go once the run ends
    private bool _leavePeriodicPending;
    private void LeavePeriodic()
    {
        if (_doc?.LongRunning == true) { _leavePeriodicPending = true; return; }
        _leavePeriodicPending = false;
        try { _doc?.SetImages(1, 1, 1, 0); } catch { }
    }

    public void RefreshPeriodic()
    {
        PbRows.Clear();
        if (_doc == null) return;
        var inv = CultureInfo.InvariantCulture;
        var r = JsonNode.Parse(_doc.Periodic(new JsonObject { ["molecule"] = (int)_pbChain }.ToJsonString()))!;
        if (r["ok"]?.GetValue<bool>() != true) { PbError = r["error"]?.GetValue<string>() ?? "no periodic cell"; return; }
        PbError = "";
        PbChains = ((int)r["molecules"]!.GetValue<double>()).ToString(inv);
        PbCross = ((int)r["crossing"]!.GetValue<double>()).ToString(inv);
        PbPieces = ((int)r["pieces"]!.GetValue<double>()).ToString(inv);
        var box = ((JsonArray)r["box"]!).Select(x => x!.GetValue<double>()).ToArray();
        PbBox = r["cubic"]?.GetValue<bool>() == true ? $"{box[0].ToString("0.0", inv)} Å cubic" : string.Join(" × ", box.Select(x => x.ToString("0.0", inv))) + " Å";
        _pbChain = (decimal)r["molecule"]!.GetValue<double>();
        Raise(nameof(PbChain));
        string D(JsonNode? n, string k) => n?[k]?.GetValue<double>().ToString("0.00", inv) ?? "—";
        if (r["bond"] is JsonNode b)
            PbRows.Add(new PeriodicRow($"Bond {b["i"]}–{b["j"]} (crosses)", D(b, "wrapped"), D(b, "min_image"), D(b, "whole")));
        else PbRows.Add(new PeriodicRow("No bond of it crosses a face", "—", "—", "—"));
        if (r["ends"] is JsonNode e)
            PbRows.Add(new PeriodicRow($"End to end ({e["i"]}–{e["j"]})", D(e, "wrapped"), D(e, "min_image"), D(e, "whole")));
        var ends = r["ends"];
        var halfBox = box.Min() / 2;
        PbNote = ends != null && ends["whole"]!.GetValue<double>() > halfBox
            ? "The minimum image fixes the bond but shortens this chain's end-to-end distance: it is longer than half the box. CAPS measures chains on the whole (unwrapped) copy."
            : "The minimum image fixes the bond; for a chain shorter than half the box it also gives the end-to-end distance. CAPS measures chains on the whole copy.";
        PeriodicChanged?.Invoke();
    }

    public void CentreOnSelection()
    {
        if (_doc == null) return;
        int[] idx = [];
        try { if (JsonNode.Parse(_doc.SelectionJson())?["indices"] is JsonArray a) idx = a.Select(x => (int)x!).ToArray(); } catch { }
        if (idx.Length == 0) idx = _selection.ToArray();
        if (idx.Length == 0) { PbError = "Select atoms first (click, or ⌘F for a query)"; return; }
        try { _doc.CentreOn(idx); }
        catch (Exception e) { PbError = e.Message; return; }
        AfterEdit($"Centred on {idx.Length} atoms (undo with ⌘Z)");
        RefreshPeriodic();
    }
}
