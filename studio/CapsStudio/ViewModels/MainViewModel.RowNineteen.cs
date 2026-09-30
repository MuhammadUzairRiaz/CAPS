using System;
using System.Collections.ObjectModel;
using System.Globalization;
using System.Linq;
using System.Text.Json.Nodes;
using System.Threading.Tasks;
using CapsStudio.Interop;

namespace CapsStudio.ViewModels;

public sealed record AutoStyleRow(string Atoms, string Default, string Detail, bool Now);
public sealed record ResolutionRow(string Quantity, string AllAtom, string UnitedAtom, string Coarse);

/// <summary>Row 19 of the design (DisplayStyles, GrowAllAtom, ModelResolution, AddHydrogens, LensView): how much of an
/// all-atom model the view draws, and conversions that make new documents; the model itself keeps every atom.</summary>
public partial class MainViewModel
{
    // ---------------------------------------------------------------- display styles (design/boards/DisplayStyles)
    public bool IsDisplayStyles => _module == 65;
    public static readonly string[] DisplayNames = ["All atoms", "No H", "Backbone"];
    private bool _dsPolarOnly, _dsSelectionFull = true;
    private string _dsAllCaption = "", _dsNoHCaption = "", _dsBackboneCaption = "", _dsModel = "";
    public string DsAllCaption { get => _dsAllCaption; private set => Set(ref _dsAllCaption, value); }
    public string DsNoHCaption { get => _dsNoHCaption; private set => Set(ref _dsNoHCaption, value); }
    public string DsBackboneCaption { get => _dsBackboneCaption; private set => Set(ref _dsBackboneCaption, value); }
    public string DsModel { get => _dsModel; private set => Set(ref _dsModel, value); }
    public ObservableCollection<AutoStyleRow> DsAuto { get; } = new();
    /// <summary>0 all atoms, 1 no H, 2 backbone (the main view's style).</summary>
    public int DsStyle
    {
        get => _style == 3 ? 1 : _style == 4 ? 2 : 0;
        set
        {
            var v = Math.Clamp(value, 0, 2);
            StyleIndex = v == 1 ? 3 : v == 2 ? 4 : DsAtoms;
            RaiseDisplay();
        }
    }
    private int _dsAtoms;
    /// <summary>How "All atoms" draws them: 0 ball and stick, 1 space filling, 2 sticks.</summary>
    public int DsAtoms
    {
        get => _dsAtoms;
        set { if (Set(ref _dsAtoms, Math.Clamp(value, 0, 2)) && DsStyle == 0) { StyleIndex = _dsAtoms; RaiseDisplay(); } }
    }
    public bool DsShowH { get => DsStyle == 0; set { DsStyle = value ? 0 : 1; } }
    public bool DsPolarOnly { get => _dsPolarOnly; set { if (Set(ref _dsPolarOnly, value)) ApplyDisplay(); } }
    public bool DsSelectionFull { get => _dsSelectionFull; set { if (Set(ref _dsSelectionFull, value)) ApplyDisplay(); } }
    public bool AutoStyleOn { get => _settings.AutoStyle; set { if (_settings.AutoStyle == value) return; _settings.AutoStyle = value; Raise(); Changed("Automatic display by size"); } }
    /// <summary>The status bar's reminder: which style is on, and what the model holds.</summary>
    public string DisplayStatus => _doc == null ? "" : "Display: " + DisplayStatusStyle() + (_lensOn ? " + lens" : "") + (_dsPolarOnly && DsStyle == 0 ? " · polar H" : "");
    // an Appearance style for all atoms is what the view draws; styles for some atoms are drawn over the display
    private string DisplayStatusStyle()
    {
        var all = AppLayers.LastOrDefault(l => l.Expression.Length == 0);
        var some = AppLayers.Count(l => l.Expression.Length > 0);
        var main = all != null ? AppStyleNames[all.Style] : DsStyle == 0 && _style is 1 or 2 ? Styles[_style] : DisplayNames[DsStyle];
        return some == 0 ? main : main + $" · {some} styled selection{(some == 1 ? "" : "s")}";
    }

    private void RaiseDisplay()
    {
        Raise(nameof(DsStyle)); Raise(nameof(DsShowH)); Raise(nameof(DisplayStatus)); Raise(nameof(DsAtoms));
    }

    public void OpenDisplayStyles() { SetModule(65); RefreshDisplayStyles(); }

    public void RefreshDisplayStyles()
    {
        DsAuto.Clear();
        long atoms = 0;
        if (_doc != null)
            try
            {
                var j = JsonNode.Parse(_doc.DisplayCounts())!;
                double D(string k) => (double?)j[k] ?? 0;
                atoms = (long)D("atoms");
                DsAllCaption = $"{Thin(D("atoms"))} atoms · {Thin(D("h"))} H";
                DsNoHCaption = $"{Thin(D("heavy"))} heavy atoms drawn · {Thin(D("h"))} H hidden";
                DsBackboneCaption = D("chains") > 0 ? $"{D("chains"):0} tubes through {Thin(D("backbone_atoms"))} backbone atoms" : "no chains: heavy atoms without H";
                DsModel = $"model: {Thin(D("atoms"))} atoms";
            }
            catch { }
        int Row(long n) => n <= _settings.AutoNoH ? 0 : n <= _settings.AutoBackbone ? 1 : n <= _settings.AutoLodAtoms ? 2 : 3;
        var now = _doc == null ? -1 : Row(atoms);
        DsAuto.Add(new($"≤ {Short(_settings.AutoNoH)}", "All atoms", "H shown", now == 0));
        DsAuto.Add(new($"{Short(_settings.AutoNoH)} – {Short(_settings.AutoBackbone)}", "No H", "H return in a lens", now == 1));
        DsAuto.Add(new($"{Short(_settings.AutoBackbone)} – {Short(_settings.AutoLodAtoms)}", "Backbone", "atoms in a 30 Å lens", now == 2));
        DsAuto.Add(new($"> {Short(_settings.AutoLodAtoms)}", "Backbone + LOD", "points far from the focus", now == 3));
        RaiseDisplay();
    }
    private static string Short(long n) => n >= 1_000_000 ? $"{n / 1_000_000.0:0.#} M" : n >= 1000 ? $"{n / 1000.0:0.#} k" : n.ToString(CultureInfo.InvariantCulture);

    /// <summary>The display options on the document (view only).</summary>
    private void ApplyDisplay()
    {
        if (_doc == null) return;
        try
        {
            var lens = new JsonObject
            {
                ["on"] = _lensOn, ["radius"] = (double)_lensRadius, ["inside"] = StyleOf(_lensInside), ["outside"] = StyleOf(_lensOutside), ["dim"] = _lensDim,
            };
            if (_lensCentre >= 0) lens["centre"] = _lensCentre;
            var clip = new JsonObject { ["on"] = _clipOn, ["axis"] = _clipAxis, ["from"] = (double)_clipFrom, ["to"] = (double)_clipTo, ["invert"] = _clipInvert };
            _doc.SetDisplay(new JsonObject { ["polar_h_only"] = _dsPolarOnly, ["selection_full"] = _dsSelectionFull, ["lens"] = lens, ["clip"] = clip }.ToJsonString());
            if (_lensOn) RefreshLensCounts();
        }
        catch (Exception e) { Status = "Display: " + e.Message; }
        Raise(nameof(DisplayStatus));
        RenderRequested?.Invoke();
    }
    private static int StyleOf(int display) => display == 1 ? 3 : display == 2 ? 4 : 0;

    /// <summary>Opening a structure: the style by its size (Settings thresholds), so a large model opens readable.</summary>
    private void AutoStyle(CapsDocument doc)
    {
        _lensCentre = -1;
        _lensOn = false;
        Raise(nameof(LensOn));
        if (_settings.AutoStyle)
        {
            var n = doc.Summary().Atoms;
            if (n > _settings.AutoNoH) StyleIndex = n <= _settings.AutoBackbone ? 3 : 4;
            if (n > _settings.AutoBackbone && n <= _settings.AutoLodAtoms) { _lensOn = true; _lensRadius = 15; Raise(nameof(LensOn)); Raise(nameof(LensRadius)); }
        }
        ApplyDisplay();
        RaiseDisplay();
    }

    public void ResetDisplay()
    {
        _dsAtoms = 0;
        StyleIndex = _settings.Style;
        _dsPolarOnly = false;
        _dsSelectionFull = true;
        Raise(nameof(DsPolarOnly)); Raise(nameof(DsSelectionFull));
        ApplyDisplay();
        RaiseDisplay();
    }
    public void SaveDisplayDefault()
    {
        _settings.Style = StyleIndex;
        Changed("Default display style");
        Status = $"New documents open as {DisplayNames[DsStyle]} (below the automatic thresholds)";
    }

    // ---------------------------------------------------------------- clip slab (design/boards/Appearance "Clip planes")
    // only the atoms between two planes normal to x, y or z are drawn (fractions of the cell); view only
    private bool _clipOn, _clipInvert;
    private int _clipAxis = 2;
    private decimal _clipFrom, _clipTo = 0.5m;
    public static readonly string[] ClipAxes = ["x", "y", "z"];
    public bool ClipOn { get => _clipOn; set { if (Set(ref _clipOn, value)) ApplyDisplay(); } }
    public int ClipAxis { get => _clipAxis; set { if (Set(ref _clipAxis, Math.Clamp(value, 0, 2))) ApplyDisplay(); } }
    public decimal ClipFrom { get => _clipFrom; set { if (Set(ref _clipFrom, Math.Clamp(value, 0m, 1m))) ApplyDisplay(); } }
    public decimal ClipTo { get => _clipTo; set { if (Set(ref _clipTo, Math.Clamp(value, 0m, 1m))) ApplyDisplay(); } }
    public bool ClipInvert { get => _clipInvert; set { if (Set(ref _clipInvert, value)) ApplyDisplay(); } }

    // ---------------------------------------------------------------- all-atom lens (design/boards/LensView)
    private bool _lensOpen, _lensOn, _lensDim, _lensMeasureInside, _lensHold;
    private int _lensInside, _lensOutside = 2, _lensCentre = -1, _lensFollows;
    private decimal _lensRadius = 10;
    public bool LensOpen
    {
        get => _lensOpen;
        set
        {
            if (!Set(ref _lensOpen, value)) return;
            if (value) { AppearanceOpen = false; SelectionOpen = false; InteractionsOpen = false; HistoryOpen = false; StatesOpen = false; LodOpen = false; if (!_lensOn) LensOn = true; else RefreshLensCounts(); }
            Raise(nameof(ShowLensPanel)); Raise(nameof(ShowStudioTabs));
        }
    }
    public bool ShowLensPanel => IsStudio && _lensOpen && _doc != null;
    public bool LensOn { get => _lensOn; set { if (Set(ref _lensOn, value)) { if (value && _lensFollows == 1) FollowSelection(); ApplyDisplay(); } } }
    public decimal LensRadius { get => _lensRadius; set { if (Set(ref _lensRadius, Math.Clamp(value, 2, 200))) ApplyDisplay(); } }
    /// <summary>0 all atoms, 1 no H, 2 backbone.</summary>
    public int LensInside { get => _lensInside; set { if (Set(ref _lensInside, Math.Clamp(value, 0, 2))) ApplyDisplay(); } }
    public int LensOutside { get => _lensOutside; set { if (Set(ref _lensOutside, Math.Clamp(value, 0, 2))) ApplyDisplay(); } }
    public bool LensDim { get => _lensDim; set { if (Set(ref _lensDim, value)) ApplyDisplay(); } }
    /// <summary>0 the cursor (hold L), 1 the selection.</summary>
    public int LensFollows { get => _lensFollows; set { if (Set(ref _lensFollows, Math.Clamp(value, 0, 1)) && value == 1) { FollowSelection(); ApplyDisplay(); } } }
    public string[] LensFollowNames { get; } = ["cursor (hold L)", "selection"];
    public bool LensMeasureInside { get => _lensMeasureInside; set => Set(ref _lensMeasureInside, value); }
    /// <summary>L held in the view: the lens follows the atom under the cursor.</summary>
    public bool LensHold { get => _lensHold; set => _lensHold = value && _lensOn && _lensFollows == 0; }
    private string _lensInsideCount = "—", _lensInsideH = "—", _lensChains = "—", _lensAtoms = "—", _lensChip = "";
    public string LensInsideCount { get => _lensInsideCount; private set => Set(ref _lensInsideCount, value); }
    public string LensInsideH { get => _lensInsideH; private set => Set(ref _lensInsideH, value); }
    public string LensChains { get => _lensChains; private set => Set(ref _lensChains, value); }
    public string LensModelAtoms { get => _lensAtoms; private set => Set(ref _lensAtoms, value); }
    public string LensChip { get => _lensChip; private set => Set(ref _lensChip, value); }
    public string LensHud => $"{DisplayNames[_lensOutside]} outside · {DisplayNames[_lensInside].ToLowerInvariant()} inside the lens · {LensModelAtoms} atoms in model";

    public void MoveLens(int atom)
    {
        if (atom < 0 || atom == _lensCentre) return;
        _lensCentre = atom;
        ApplyDisplay();
    }
    private void FollowSelection() { if (_selection.Count > 0) _lensCentre = _selection[0]; }

    private void RefreshLensCounts()
    {
        if (_doc == null) return;
        try
        {
            var j = JsonNode.Parse(_doc.DisplayCounts())!;
            double D(string k) => (double?)j[k] ?? double.NaN;
            if (double.IsFinite(D("lens_centre"))) _lensCentre = (int)D("lens_centre");
            LensInsideCount = double.IsFinite(D("lens_atoms")) ? Thin(D("lens_atoms")) : "—";
            LensInsideH = double.IsFinite(D("lens_h")) ? Thin(D("lens_h")) : "—";
            LensChains = Thin(D("chains"));
            LensModelAtoms = Thin(D("atoms"));
            LensChip = $"lens {LensInsideCount} atoms";
            Raise(nameof(LensHud));
        }
        catch { }
    }

    /// <summary>The lens outline for the view: its centre atom's position and radius in view points, and a label.</summary>
    public (double X, double Y, double R, string Label)? LensCircle(CapsCamera cam, CapsRenderOpts opt, double scaling)
    {
        if (_doc == null || !_lensOn || _lensCentre < 0 || !IsStudio) return null;
        try
        {
            var p = _doc.ProjectIndices(cam, opt, [_lensCentre]);   // the lens's atom only
            var scale = _doc.ViewScale(cam, opt);
            return (p[0] / scaling, p[1] / scaling, (double)_lensRadius * scale / scaling, $"lens · {_lensRadius:0.#} Å");
        }
        catch { return null; }
    }
    /// <summary>The lens circle from the GPU view's camera fit (one atom looked up; nothing over the whole structure).</summary>
    public (double X, double Y, double R, string Label)? LensCircleFit(in CapsViewFit fit, double scaling)
    {
        if (_doc == null || !_lensOn || _lensCentre < 0 || !IsStudio) return null;
        try
        {
            if (AtomXyz(_lensCentre) is not { } p) return null;
            var (x, y, k) = ProjectFit(fit, p);
            return (x / scaling, y / scaling, (double)_lensRadius * fit.Scale * k / scaling, $"lens · {_lensRadius:0.#} Å");
        }
        catch { return null; }
    }

    /// <summary>One atom's position (the current frame), or null.</summary>
    public (double X, double Y, double Z)? AtomXyz(int atom)
    {
        if (_doc == null || atom < 0) return null;
        if (System.Text.Json.Nodes.JsonNode.Parse(_doc.AtomProperties(atom))?["xyz"] is System.Text.Json.Nodes.JsonArray a && a.Count == 3)
            return ((double)a[0]!, (double)a[1]!, (double)a[2]!);
        return null;
    }

    /// <summary>A world point on screen for a camera fit (the core's projection): x, y in the fit's pixels, the perspective factor.</summary>
    public static (double X, double Y, double K) ProjectFit(in CapsViewFit f, (double X, double Y, double Z) p)
    {
        double dx = p.X - f.Cx, dy = p.Y - f.Cy, dz = p.Z - f.Cz;
        var x = dx * f.CosYaw + dz * f.SinYaw + f.PanX;
        var z = -dx * f.SinYaw + dz * f.CosYaw;
        var y = dy * f.CosPitch - z * f.SinPitch + f.PanY;
        var z2 = dy * f.SinPitch + z * f.CosPitch;
        var k = f.Perspective != 0 ? f.Dist / Math.Max(1e-3, f.Dist - z2) : 1.0;
        return (f.W / 2 + x * f.Scale * k, f.H / 2 - y * f.Scale * k, k);
    }

    public bool PickAllowed(int atom) => atom < 0 || !_lensOn || !_lensMeasureInside || _doc == null || _doc.LensInside(atom);

    // ---------------------------------------------------------------- add hydrogens (design/boards/AddHydrogens)
    public bool IsAddHydrogens => _module == 66;
    public ObservableCollection<HydrogenRow> AhRows { get; } = new();
    private bool _ahKeep = true, _ahRelax = true;
    public bool AhKeep { get => _ahKeep; set { if (Set(ref _ahKeep, value)) _ = AhPreview(); } }
    // protonation by pH (amino-acid residues, model pKa values): set on the document, used by the plan and add_h
    private bool _ahPhOn;
    private decimal _ahPh = 7.0m;
    public bool AhPhOn { get => _ahPhOn; set { if (Set(ref _ahPhOn, value)) ApplyAhPh(); } }
    public decimal AhPhD { get => _ahPh; set { if (Set(ref _ahPh, Math.Clamp(value, 0m, 14m)) && _ahPhOn) ApplyAhPh(); } }
    private void ApplyAhPh()
    {
        try { _doc?.SetPh(_ahPhOn ? (double)_ahPh : -1); } catch { }
        _ = AhPreview();
    }
    public bool AhRelax { get => _ahRelax; set => Set(ref _ahRelax, value); }
    private string _ahSource = "", _ahOrders = "", _ahHeavy = "—", _ahAdded = "—", _ahCharge = "—", _ahBefore = "", _ahAfter = "", _ahStatus = "", _ahTotal = "";
    public string AhSource { get => _ahSource; private set => Set(ref _ahSource, value); }
    public string AhOrders { get => _ahOrders; private set => Set(ref _ahOrders, value); }
    public string AhHeavy { get => _ahHeavy; private set => Set(ref _ahHeavy, value); }
    public string AhAdded { get => _ahAdded; private set => Set(ref _ahAdded, value); }
    private bool _ahCanAdd;
    /// <summary>The plan adds at least one hydrogen (else the button is off: nothing is missing).</summary>
    public bool AhCanAdd { get => _ahCanAdd; private set => Set(ref _ahCanAdd, value); }
    public string AhCharge { get => _ahCharge; private set => Set(ref _ahCharge, value); }
    public string AhBeforeCaption { get => _ahBefore; private set => Set(ref _ahBefore, value); }
    public string AhAfterCaption { get => _ahAfter; private set => Set(ref _ahAfter, value); }
    public string AhStatus { get => _ahStatus; private set => Set(ref _ahStatus, value); }
    public string AhTotal { get => _ahTotal; private set => Set(ref _ahTotal, value); }
    private CapsDocument? _ahBeforeDoc, _ahAfterDoc;
    public CapsDocument? AhBeforeDoc { get => _ahBeforeDoc; private set => Set(ref _ahBeforeDoc, value); }
    public CapsDocument? AhAfterDoc { get => _ahAfterDoc; private set => Set(ref _ahAfterDoc, value); }
    private bool _ahApplied;
    public bool AhCanUndo => _ahApplied;

    public async Task OpenAddHydrogens()
    {
        SetModule(66);
        _ahApplied = false;
        Raise(nameof(AhCanUndo));
        await AhPreview();
    }

    /// <summary>The structure as it is and a copy with the hydrogens added (new ones selected), and the rule table.</summary>
    public async Task AhPreview()
    {
        if (_doc == null) return;
        try
        {
            var doc = _doc;
            var plan = JsonNode.Parse(doc.HydrogenPlan())!;
            double D(string k) => (double?)plan[k] ?? 0;
            AhRows.Clear();
            foreach (var r in (plan["rows"] as JsonArray ?? []).OfType<JsonObject>())
                AhRows.Add(new((string?)r["label"] ?? "", ((double?)r["atoms"] ?? 0).ToString("0", Inv), ((double?)r["hydrogens"] ?? 0).ToString("0", Inv)));
            var s = doc.Summary();
            AhSource = System.IO.Path.GetFileName(Title.Replace(" (unsaved)", "")) + (D("h") == 0 ? " · heavy atoms only" : $" · {D("h"):0} H present");
            var geo = (bool?)plan["orders_from_geometry"] == true;
            AhOrders = (geo ? "from geometry" : "from the file") + (D("aromatic_bonds") > 0 ? $" · {D("aromatic_bonds"):0} aromatic bonds" : "");
            AhHeavy = D("heavy").ToString("0", Inv);
            AhCharge = D("net_charge").ToString("0.##", Inv);
            AhBeforeCaption = $"{D("heavy"):0} heavy atoms · {D("h"):0} H";
            // the preview: a copy with the hydrogens added (and the existing ones removed first when not kept)
            var before = doc.Copy("as it is");
            var after = doc.Copy("with hydrogens");
            int n0 = (int)after.Summary().Atoms;
            await Task.Run(() =>
            {
                if (!_ahKeep && D("h") > 0)   // replace the hydrogens present: remove them first
                {
                    after.Select("{\"mode\":\"element\",\"pattern\":\"H\"}");
                    after.Edit("{\"op\":\"delete\",\"atoms\":\"selection\"}");
                }
                after.Edit(_ahPhOn ? new JsonObject { ["op"] = "add_h", ["atoms"] = "", ["ph"] = (double)_ahPh }.ToJsonString() : "{\"op\":\"add_h\",\"atoms\":\"\"}");
            });
            var n1 = (int)after.Summary().Atoms;
            var added = n1 - n0;
            var first = (int)after.Summary().Atoms - Math.Max(0, added);
            try { after.Select(new JsonObject { ["mode"] = "indices", ["atoms"] = new JsonArray(Enumerable.Range(first, Math.Max(0, added)).Select(i => (JsonNode)i).ToArray()) }.ToJsonString()); } catch { }
            AhAdded = added.ToString(Inv);
            AhCanAdd = added > 0;
            AhAfterCaption = $"{n1} atoms · +{added} H";
            AhTotal = $"Rule total {D("add"):0} H" + (D("h") > 0 && _ahKeep ? $" (the {D("h"):0} present kept)" : "") + ".";
            AhStatus = $"{AhSource} · {n1} atoms after · {n1 - (int)D("heavy")} H";
            var ob = AhBeforeDoc; var oa = AhAfterDoc;
            AhBeforeDoc = before; AhAfterDoc = after;
            ob?.Dispose(); oa?.Dispose();
        }
        catch (Exception e) { AhStatus = "Add hydrogens: " + e.Message; }
    }

    /// <summary>Adds the hydrogens to the document (one undoable edit), then relaxes only the new ones when asked.</summary>
    public async Task AhApply()
    {
        if (_doc == null || !Idle) return;
        var n0 = (int)_doc.Summary().Atoms;
        AddHydrogensAll();
        var n1 = (int)_doc.Summary().Atoms;
        if (_ahRelax && n1 > n0)
        {
            var doc = _doc;
            var ids = new JsonArray(Enumerable.Range(n0, n1 - n0).Select(i => (JsonNode)i).ToArray());
            try { await Task.Run(() => doc.Edit(new JsonObject { ["op"] = "clean", ["atoms"] = ids, ["ftol"] = 0.5 }.ToJsonString())); } catch (Exception e) { Status = "Relax of the new H: " + e.Message; }
        }
        _ahApplied = true;
        Raise(nameof(AhCanUndo));
        RefreshSummary();
        RenderRequested?.Invoke();
        Status = $"Added {n1 - n0} hydrogens" + (_ahRelax ? ", only they were relaxed (UFF)" : "") + " · undo with ⌘Z";
        await AhPreview();
    }
    public async Task AhUndo()
    {
        if (!_ahApplied) return;
        UndoEdit(false);
        if (_ahRelax) UndoEdit(false);
        _ahApplied = false;
        Raise(nameof(AhCanUndo));
        await AhPreview();
    }

    // ---------------------------------------------------------------- model resolution (design/boards/ModelResolution)
    public bool IsModelResolution => _module == 67;
    public ObservableCollection<ResolutionRow> MrRows { get; } = new();
    private CapsDocument? _mrUa, _mrCg;
    public CapsDocument? MrUaDoc { get => _mrUa; private set => Set(ref _mrUa, value); }
    public CapsDocument? MrCgDoc { get => _mrCg; private set => Set(ref _mrCg, value); }
    private string _mrAaCaption = "", _mrUaCaption = "", _mrCgCaption = "", _mrTitle = "", _mrStatus = "", _mrKind = "all-atom model";
    public string MrAaCaption { get => _mrAaCaption; private set => Set(ref _mrAaCaption, value); }
    public string MrUaCaption { get => _mrUaCaption; private set => Set(ref _mrUaCaption, value); }
    public string MrCgCaption { get => _mrCgCaption; private set => Set(ref _mrCgCaption, value); }
    public string MrTitle { get => _mrTitle; private set => Set(ref _mrTitle, value); }
    public string MrStatus { get => _mrStatus; private set => Set(ref _mrStatus, value); }
    public string MrKind { get => _mrKind; private set => Set(ref _mrKind, value); }
    private decimal _mrPerBead = 5;
    public decimal MrPerBead { get => _mrPerBead; set { if (Set(ref _mrPerBead, Math.Clamp(value, 1, 50))) RefreshResolution(); } }
    /// <summary>What this document is (0 all-atom, 1 united-atom, 2 coarse-grained): picking another converts to a new document.</summary>
    public int MrResolution { get => _mrRes; set { var v = Math.Clamp(value, 0, 2); if (v == _mrRes) return; _ = ConvertResolution(v); } }
    private int _mrRes;

    public void OpenModelResolution() { SetModule(67); RefreshResolution(); }

    public void RefreshResolution()
    {
        if (_doc == null) return;
        try
        {
            var r = JsonNode.Parse(_doc.ResolutionSummary(new JsonObject { ["per_bead"] = (int)_mrPerBead }.ToJsonString()))!;
            var aa = r["all_atom"]!; var ua = r["united_atom"]!; var cg = r["coarse_grained"]!;
            double D(JsonNode n, string k) => (double?)n[k] ?? 0;
            // what the open document is: a model without hydrogens on carbon reads as united-atom, beads as coarse-grained
            var info = JsonNode.Parse(_doc.DisplayCounts())!;
            var hasH = ((double?)info["h"] ?? 0) > 0;
            _mrRes = Title.Contains("coarse-grained", StringComparison.Ordinal) ? 2 : !hasH && D(aa, "sites") == D(ua, "sites") ? 1 : 0;
            Raise(nameof(MrResolution));
            MrKind = _mrRes == 0 ? "all-atom model" : _mrRes == 1 ? "united-atom model" : "coarse-grained model";
            MrAaCaption = $"{D(aa, "sites"):0} sites · {D(aa, "hydrogens"):0} H";
            MrUaCaption = $"{D(ua, "sites"):0} sites · H folded in";
            MrCgCaption = $"{D(cg, "sites"):0} beads · {(int)_mrPerBead} backbone atoms each";
            MrRows.Clear();
            MrRows.Add(new("Sites", D(aa, "sites").ToString("0", Inv), D(ua, "sites").ToString("0", Inv), D(cg, "sites").ToString("0", Inv)));
            MrRows.Add(new("Hydrogens", $"{D(aa, "hydrogens"):0} explicit", "inside CH₂/CH₃" + (D(ua, "hydrogens") > 0 ? $" · {D(ua, "hydrogens"):0} polar kept" : ""), "inside beads"));
            MrRows.Add(new("Mass (g/mol)", D(aa, "mass").ToString("0.000", Inv), D(ua, "mass").ToString("0.000", Inv), D(cg, "mass").ToString("0.000", Inv)));
            MrRows.Add(new("Example force field", "GAFF2 · OPLS-AA", "TraPPE-UA", "Kremer–Grest (generic)"));
            MrRows.Add(new("Charges", "on every atom incl. H", "on united sites", "usually none"));
            MrRows.Add(new("Use it for", "chemistry, H-bonds, spectra", "alkanes, melts, faster", "long times, entanglement"));
            MrTitle = "Same structure, three resolutions · " + Title.Replace(" (unsaved)", "");
            MrStatus = $"{Title.Replace(" (unsaved)", "")} · {D(aa, "sites"):0} atoms · {D(aa, "mass"):0.000} g/mol";
            // the two conversions, as previews (not opened)
            var doc = _doc;
            var per = (int)_mrPerBead;
            var (u, _) = doc.ResolutionConvert("{\"to\":\"united-atom\"}", "united-atom");
            var (c, _) = doc.ResolutionConvert($"{{\"to\":\"coarse-grained\",\"per_bead\":{per}}}", "coarse-grained");
            var ou = MrUaDoc; var oc = MrCgDoc;
            MrUaDoc = u; MrCgDoc = c;
            ou?.Dispose(); oc?.Dispose();
        }
        catch (Exception e) { MrStatus = "Resolution: " + e.Message; }
    }

    /// <summary>Converting makes a new document; the original stays open in its own tab (Recent).</summary>
    public async Task ConvertResolution(int to)
    {
        if (_doc == null || !Idle) return;
        if (to == 0)
        {
            if (_mrRes != 1) { Status = "Coarse-grained → all-atom needs a backmap: see Builders › Coarse-grained"; Raise(nameof(MrResolution)); return; }
            // united-atom → all-atom: hydrogens on every carbon from valence, then only they are relaxed
            var copy = _doc.Copy("all-atom");
            var n0 = (int)copy.Summary().Atoms;
            await Task.Run(() => copy.Edit("{\"op\":\"add_h\",\"atoms\":\"\"}"));
            var n1 = (int)copy.Summary().Atoms;
            try { await Task.Run(() => copy.Edit(new JsonObject { ["op"] = "clean", ["atoms"] = new JsonArray(Enumerable.Range(n0, n1 - n0).Select(i => (JsonNode)i).ToArray()), ["ftol"] = 0.5 }.ToJsonString())); } catch { }
            Show(copy, Title.Replace(" (united-atom)", "").Replace(" (unsaved)", "") + " (all-atom, unsaved)");
            GrownUnsaved = true;
            SetModule(67);
            RefreshResolution();
            Status = $"All-atom again: {n1 - n0} H added at 1.09 Å and relaxed (UFF); the united-atom document is unchanged";
            return;
        }
        try
        {
            var json = to == 1 ? "{\"to\":\"united-atom\"}" : $"{{\"to\":\"coarse-grained\",\"per_bead\":{(int)_mrPerBead}}}";
            var name = Title.Replace(" (unsaved)", "") + (to == 1 ? " (united-atom" : " (coarse-grained") + ", unsaved)";
            var (d, report) = _doc.ResolutionConvert(json, name);
            Show(d, name);
            GrownUnsaved = true;
            SetModule(67);
            RefreshResolution();
            Status = "Converted: " + report.Replace("\n", " · ") + " · the original is unchanged (reopen it from Recent)";
        }
        catch (Exception e) { Status = "Conversion: " + e.Message; }
    }
}

public partial class MainViewModel
{
    /// <summary>Coarse-grained → all-atom: the open all-atom structure's beads (MrPerBead each) moved to the file's
    /// beads, the atoms carried with them, then relaxed (a new document).</summary>
    public async Task BackmapFrom(string beadsPath)
    {
        if (_doc == null || !Idle) return;
        var doc = _doc;
        var per = (int)_mrPerBead;
        var name = Title.Replace(" (unsaved)", "") + " (backmapped, unsaved)";
        Status = "Backmapping onto " + System.IO.Path.GetFileName(beadsPath) + " and relaxing…";
        try
        {
            var (d, report) = await Task.Run(() => doc.Backmap(beadsPath, per, true, name));
            Show(d, name);
            GrownUnsaved = true;
            SetModule(67);
            RefreshResolution();
            Status = "Backmapped: " + report.Split('\n').FirstOrDefault();
        }
        catch (Exception e) { Status = "Backmap: " + e.Message; }
    }
}

public sealed record HydrogenRow(string Atom, string Count, string Added);
