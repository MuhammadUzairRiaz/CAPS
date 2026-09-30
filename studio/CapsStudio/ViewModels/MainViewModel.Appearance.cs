using System.Collections.ObjectModel;
using System.Globalization;
using System.Text.Json.Nodes;
using CapsStudio.Interop;

namespace CapsStudio.ViewModels;

/// <summary>A style layer: the atoms an expression picks (empty: all) drawn in one style.</summary>
public sealed record StyleLayer(string Expression, int Style)
{
    public string Title => (Expression.Length == 0 ? "All atoms" : Expression) + " · " + MainViewModel.AppStyleNames[Style];
}

/// <summary>An atom label of the view: where it is drawn and its text.</summary>
/// <summary>One label kind (atom or bond property) and whether it is shown.</summary>
public sealed class LabelKindRow : ObservableObject
{
    private bool _on;
    private readonly Action _changed;
    public LabelKindRow(string id, string title, string group, bool bond, Action changed) { Id = id; Title = title; Group = group; Bond = bond; _changed = changed; }
    public string Id { get; }
    public string Title { get; }
    public string Group { get; }
    public bool Bond { get; }
    public bool On { get => _on; set { if (Set(ref _on, value)) _changed(); } }
}

public readonly record struct ViewLabel(double X, double Y, string Text, uint Argb = 0, bool Bond = false);
/// <summary>How labels look: font, size (pt), bold, a plate behind them, and the atom and bond colours (ARGB; null: the theme's text).</summary>
public sealed record LabelStyle(string Font, double Size, bool Bold, bool Plate, uint? AtomArgb, uint? BondArgb);

/// <summary>Appearance (design/boards/Appearance): per-selection styles (ball and stick, stick, wireframe, space
/// filling, coordination polyhedra, ribbon), colour by partial charge with a colour-blind-safe diverging ramp, atom
/// labels (element, R/S, type, charge), molecular surfaces coloured by the electrostatic potential, rendering switches.</summary>
public sealed partial class MainViewModel
{
    public static readonly string[] AppStyleNames = ["Ball & stick", "Stick", "Wireframe", "Space filling", "Polyhedra", "Ribbon"];
    private static readonly string[] AppStyleCodes = ["ball_and_stick", "sticks", "wireframe", "space_filling", "polyhedra", "ribbon"];
    public static readonly string[] AppColourNames = ["Element", "Molecule", "Type", "Distance to molecule centre", "Partial charge"];
    public static readonly string[] AppRampNames = ["Blue–orange (CB-safe)", "Viridis", "Red–white–blue"];
    private static readonly string[] AppRampCodes = ["blue_orange", "viridis", "red_white_blue"];
    public static readonly string[] AppSurfaceNames = ["None", "Solvent-accessible", "Solvent-excluded", "van der Waals"];
    private static readonly string[] AppSurfaceCodes = ["none", "accessible", "excluded", "vdw"];
    public static readonly string[] AppMapNames = ["Electrostatic potential", "Nearest atom (element)", "One colour"];
    private static readonly string[] AppMapCodes = ["potential", "atom", "uniform"];
    public static readonly string[] AppOpacities = ["20 %", "40 %", "60 %", "80 %", "100 %"];

    private bool _appOpen;
    /// <summary>The Appearance panel in place of the Studio inspector.</summary>
    public bool AppearanceOpen { get => _appOpen; set { if (Set(ref _appOpen, value)) { if (value) { SelectionOpen = false; InteractionsOpen = false; LodOpen = false; HistoryOpen = false; StatesOpen = false; LensOpen = false; } RaiseAppearanceVisibility(); if (value) AppRefreshInfo(); } } }
    public bool ShowStudioTabs => IsStudio && !_appOpen && !_selOpen && !_ixOpen && !_lodOpen && !_histOpen && !_lensOpen && !_stOpen;
    public bool ShowAppearance => IsStudio && _appOpen && _doc != null;
    private void RaiseAppearanceVisibility() { Raise(nameof(ShowStudioTabs)); Raise(nameof(ShowAppearance)); Raise(nameof(ShowAppLegend)); Raise(nameof(ShowSelectionPanel)); Raise(nameof(ShowInteractionsPanel)); }

    // ---------------------------------------------------------------- styles

    public ObservableCollection<StyleLayer> AppLayers { get; } = new();
    private int _appStyle;
    /// <summary>The style tile; choosing one styles the target atoms at once.</summary>
    public int AppStyle
    {
        get => _appStyle;
        set
        {
            _appStyle = Math.Clamp(value, 0, AppStyleNames.Length - 1);
            Raise();
            RaiseStyleTiles();
            var expr = _appTarget == 0 ? "" : _appExpr.Trim();
            if (expr.Length == 0) AppLayers.Clear();
            else foreach (var l in AppLayers.Where(l => l.Expression == expr).ToList()) AppLayers.Remove(l);
            AppLayers.Add(new StyleLayer(expr, _appStyle));
            Raise(nameof(AppChip));
            ApplyAppearance();
        }
    }
    public bool AppBallStick { get => _appStyle == 0; set { if (value) AppStyle = 0; } }
    public bool AppStick { get => _appStyle == 1; set { if (value) AppStyle = 1; } }
    public bool AppWire { get => _appStyle == 2; set { if (value) AppStyle = 2; } }
    public bool AppSpace { get => _appStyle == 3; set { if (value) AppStyle = 3; } }
    public bool AppPoly { get => _appStyle == 4; set { if (value) AppStyle = 4; } }
    public bool AppRibbon { get => _appStyle == 5; set { if (value) AppStyle = 5; } }
    private void RaiseStyleTiles() { foreach (var n in new[] { nameof(AppBallStick), nameof(AppStick), nameof(AppWire), nameof(AppSpace), nameof(AppPoly), nameof(AppRibbon) }) Raise(n); }

    private int _appTarget;
    /// <summary>0 all atoms, 1 the atoms an expression picks.</summary>
    public int AppTarget { get => _appTarget; set { if (Set(ref _appTarget, Math.Clamp(value, 0, 1))) Raise(nameof(AppTargetIsExpression)); } }
    public bool AppTargetIsExpression => _appTarget == 1;
    private string _appExpr = "Molecule <= 3";
    public string AppExpression { get => _appExpr; set => Set(ref _appExpr, value ?? ""); }
    public string AppChip => AppLayers.Count <= 1 ? (AppLayers.FirstOrDefault()?.Expression is { Length: > 0 } e ? "Selection: " + e : "All atoms") : $"Mixed styles: {AppLayers.Count} selections";

    public void RemoveStyleLayer(StyleLayer l)
    {
        AppLayers.Remove(l);
        Raise(nameof(AppChip));
        ApplyAppearance();
    }

    // ---------------------------------------------------------------- colour

    private int _appColour = -1, _appRamp;
    /// <summary>Colour by: element, molecule, type, distance (the view's colour modes) or partial charge.</summary>
    public int AppColour
    {
        get => _appColour >= 0 ? _appColour : ColourIndex;
        set
        {
            var v = Math.Clamp(value, 0, AppColourNames.Length - 1);
            _appColour = v;
            if (v < 4) ColourIndex = v;
            Raise();
            Raise(nameof(ColourText));
            Raise(nameof(AppRampEnabled));
            Raise(nameof(ShowAppLegend));
            ApplyAppearance();
        }
    }
    // a per-atom column of the trajectory (a LAMMPS dump's c_pe, fx, |f| …): colour by it in place of the choice above
    public System.Collections.ObjectModel.ObservableCollection<string> AppColumns { get; } = new();
    public bool AppHasColumns => AppColumns.Count > 1;
    private int _appColumn;
    public int AppColumn
    {
        get => _appColumn;
        set
        {
            if (!Set(ref _appColumn, Math.Clamp(value, 0, Math.Max(0, AppColumns.Count - 1)))) return;
            Raise(nameof(AppRampEnabled));
            ApplyAppearance();
        }
    }
    private string? ColumnName => _appColumn > 0 && _appColumn < AppColumns.Count ? AppColumns[_appColumn] : null;
    private void RefreshAppColumns()
    {
        AppColumns.Clear();
        _appColumn = 0;
        try
        {
            if (_doc != null && System.Text.Json.Nodes.JsonNode.Parse(_doc.TrajectoryColumns() is { Length: > 0 } t ? t : "{}")?["columns"] is JsonArray cols && cols.Count > 0)
            {
                AppColumns.Add("None (the colour above)");
                foreach (var c in cols) if ((string?)c?["name"] is { } n) AppColumns.Add(n);
            }
        }
        catch { }
        Raise(nameof(AppColumn));
        Raise(nameof(AppHasColumns));
    }
    public int AppRamp { get => _appRamp; set { if (Set(ref _appRamp, Math.Clamp(value, 0, AppRampNames.Length - 1))) { Raise(nameof(AppLegendBrush)); ApplyAppearance(); } } }
    public bool AppRampEnabled => AppColour == 4 || _appSurface > 0 || ColumnName != null;

    // ---------------------------------------------------------------- labels

    private List<LabelKindRow>? _atomKinds, _bondKinds;
    private void LoadLabelKinds()
    {
        if (_atomKinds != null) return;
        _atomKinds = new();
        _bondKinds = new();
        try
        {
            var j = JsonNode.Parse(CapsDocument.LabelKinds())!;
            foreach (var (key, list, bond) in new[] { ("atom", _atomKinds, false), ("bond", _bondKinds, true) })
                foreach (var k in j[key]!.AsArray())
                    list.Add(new LabelKindRow((string)k!["id"]!, (string)k["title"]!, (string)k["group"]!, bond, LabelsChanged));
        }
        catch { }
    }
    /// <summary>The atom properties that can label atoms (element, name, charge, oxidation state, hybridisation, x y z …).</summary>
    public List<LabelKindRow> AtomLabelKinds { get { LoadLabelKinds(); return _atomKinds!; } }
    /// <summary>The bond properties that can label bonds, at their midpoints (order, chemical type, length, force-field r₀ …).</summary>
    public List<LabelKindRow> BondLabelKinds { get { LoadLabelKinds(); return _bondKinds!; } }
    private LabelKindRow? Kind(string id) => AtomLabelKinds.FirstOrDefault(k => k.Id == id);
    private bool KindOn(string id) => Kind(id)?.On ?? false;
    private void SetKind(string id, bool on) { if (Kind(id) is { } k) k.On = on; }

    // the four labels the page had first (sessions, scripts and screenshots set these)
    public bool LabelElement { get => KindOn("element"); set => SetKind("element", value); }
    public bool LabelRs { get => KindOn("stereo"); set => SetKind("stereo", value); }
    public bool LabelType { get => KindOn("type"); set => SetKind("type", value); }
    public bool LabelCharge { get => KindOn("charge"); set => SetKind("charge", value); }
    public bool AnyLabels => AtomLabelKinds.Any(k => k.On) || BondLabelKinds.Any(k => k.On);

    // which atoms: all, the selection, all but hydrogens
    public static readonly string[] LabelScopes = ["All atoms", "Selected atoms", "All but hydrogens"];
    private int _labelScope;
    public int LabelScope { get => _labelScope; set { if (Set(ref _labelScope, Math.Clamp(value, 0, 2))) LabelsChanged(); } }

    // how labels look
    private static string[]? _labelFonts;
    /// <summary>The fonts labels can use: the Studio's own, then the system's.</summary>
    public static string[] LabelFonts => _labelFonts ??= new[] { "IBM Plex Mono", "IBM Plex Sans", "Inter" }
        .Concat(SystemFontNames()).Distinct().ToArray();
    private static IEnumerable<string> SystemFontNames()
    {
        try { return Avalonia.Media.FontManager.Current.SystemFonts.Select(f => f.Name).Where(n => !n.StartsWith('.')).OrderBy(n => n, StringComparer.OrdinalIgnoreCase).ToList(); }
        catch { return Array.Empty<string>(); }
    }
    private int _labelFont;
    private decimal _labelSize = 10.5m;
    private bool _labelBold, _labelPlate = true;
    public int LabelFont { get => _labelFont; set { if (Set(ref _labelFont, Math.Clamp(value, 0, LabelFonts.Length - 1))) RenderRequested?.Invoke(); } }
    public decimal? LabelSize { get => _labelSize; set { if (value != null && Set(ref _labelSize, Math.Clamp(value.Value, 5, 48))) RenderRequested?.Invoke(); } }
    public bool LabelBold { get => _labelBold; set { if (Set(ref _labelBold, value)) RenderRequested?.Invoke(); } }
    /// <summary>A translucent plate behind each label, so it reads over any colour.</summary>
    public bool LabelPlate { get => _labelPlate; set { if (Set(ref _labelPlate, value)) RenderRequested?.Invoke(); } }
    public static readonly string[] LabelColours = ["Theme text", "By element", "Accent", "White", "Black", "Red", "Blue", "Green", "Yellow", "Custom"];
    private static readonly uint[] LabelColourRgb = [0, 0, 0, 0xFFFFFF, 0x000000, 0xE5484D, 0x3E7BFA, 0x2FA36B, 0xF0C419, 0];
    private int _atomLabelColour, _bondLabelColour = 1;   // bonds in the accent colour
    private string _atomLabelHex = "#FF7A00", _bondLabelHex = "#3E7BFA";
    public int AtomLabelColour { get => _atomLabelColour; set { if (Set(ref _atomLabelColour, Math.Clamp(value, 0, LabelColours.Length - 1))) { Raise(nameof(AtomLabelCustom)); LabelsChanged(); } } }
    public int BondLabelColour { get => _bondLabelColour; set { if (Set(ref _bondLabelColour, Math.Clamp(value, 0, LabelColours.Length - 2))) { Raise(nameof(BondLabelCustom)); RenderRequested?.Invoke(); } } }
    public bool AtomLabelCustom => _atomLabelColour == LabelColours.Length - 1;
    public bool BondLabelCustom => _bondLabelColour == LabelColours.Length - 1;
    /// <summary>The custom colours as #RRGGBB.</summary>
    public string AtomLabelHex { get => _atomLabelHex; set { if (Set(ref _atomLabelHex, value ?? "")) RenderRequested?.Invoke(); } }
    public string BondLabelHex { get => _bondLabelHex; set { if (Set(ref _bondLabelHex, value ?? "")) RenderRequested?.Invoke(); } }
    public static readonly string[] BondLabelColours = LabelColours.Where(c => c != "By element").ToArray();

    /// <summary>How the labels are drawn (the overlay reads it).</summary>
    public LabelStyle LabelLook => new(LabelFonts[Math.Clamp(_labelFont, 0, LabelFonts.Length - 1)], (double)_labelSize, _labelBold, _labelPlate,
        ColourOf(_atomLabelColour, _atomLabelHex, LabelColours), ColourOf(_bondLabelColour, _bondLabelHex, BondLabelColours));
    // null: the theme's text colour; 0xFF000000 | rgb otherwise; "By element" is resolved per label
    private static uint? ColourOf(int pick, string hex, string[] names)
    {
        var name = names[Math.Clamp(pick, 0, names.Length - 1)];
        if (name == "Theme text" || name == "By element") return null;
        if (name == "Accent") return 0xFFF0A83C;
        if (name == "Custom") return Avalonia.Media.Color.TryParse(hex.Trim(), out var c) ? c.ToUInt32() : null;
        return 0xFF000000 | LabelColourRgb[Array.IndexOf(LabelColours, name)];
    }

    private void LabelsChanged() { _labelTexts = null; _bondTexts = null; Raise(nameof(AnyLabels)); RenderRequested?.Invoke(); }
    private string[]? _labelTexts;
    private uint[]? _labelRgb;
    private (int[] Pairs, string[] Texts, bool[] Crossing)? _bondTexts;
    private object? _labelDoc;
    private int _labelFrame = -1;
    /// <summary>Most labels drawn at once; above this the labels wait for a smaller view.</summary>
    public const int LabelLimit = 4000;
    private string _labelNote = "";
    public string LabelNote { get => _labelNote; private set => Set(ref _labelNote, value); }

    private bool LabelCacheValid => ReferenceEquals(_labelDoc, _doc) && _labelFrame == _frame;

    /// <summary>The label texts per atom (the chosen kinds joined), fetched once per document and frame.</summary>
    private string[]? LabelTexts()
    {
        if (_doc == null) return null;
        if (_labelTexts != null && LabelCacheValid) return _labelTexts;
        if (!LabelCacheValid) _bondTexts = null;
        var kinds = AtomLabelKinds.Where(k => k.On).Select(k => k.Id).ToList();
        var cols = new List<string[]>();
        foreach (var k in kinds)
        {
            try { cols.Add(JsonNode.Parse(_doc.AtomLabels(k))!.AsArray().Select(x => x?.GetValue<string>() ?? "").ToArray()); }
            catch { }
        }
        var n = cols.Count > 0 ? cols[0].Length : 0;
        var texts = new string[n];
        for (var i = 0; i < n; i++) texts[i] = string.Join(" ", cols.Select(c => i < c.Length ? c[i] : "").Where(x => x.Length > 0));
        // by element: the element's display colour
        _labelRgb = null;
        if (_atomLabelColour == 1 && n > 0)
        {
            try
            {
                var hexes = JsonNode.Parse(_doc.AtomLabels("colour"))!.AsArray().Select(x => x?.GetValue<string>() ?? "").ToArray();
                _labelRgb = hexes.Select(h => h.Length == 7 && uint.TryParse(h[1..], System.Globalization.NumberStyles.HexNumber, null, out var v) ? 0xFF000000 | v : 0u).ToArray();
            }
            catch { }
        }
        _labelTexts = texts;
        _labelDoc = _doc;
        _labelFrame = _frame;
        return texts;
    }

    private (int[] Pairs, string[] Texts, bool[] Crossing)? BondTexts()
    {
        if (_doc == null) return null;
        if (_bondTexts != null && LabelCacheValid) return _bondTexts;
        var kinds = BondLabelKinds.Where(k => k.On).Select(k => k.Id).ToList();
        if (kinds.Count == 0) return null;
        int[]? pairs = null;
        bool[]? crossing = null;
        var cols = new List<string[]>();
        foreach (var k in kinds)
        {
            try
            {
                var j = JsonNode.Parse(_doc.BondLabels(k))!;
                pairs ??= j["pairs"]!.AsArray().Select(x => (int)x!.GetValue<double>()).ToArray();
                crossing ??= j["crossing"]!.AsArray().Select(x => x!.GetValue<bool>()).ToArray();
                cols.Add(j["labels"]!.AsArray().Select(x => x?.GetValue<string>() ?? "").ToArray());
            }
            catch { }
        }
        if (pairs == null || crossing == null) return null;
        var n = pairs.Length / 2;
        var texts = new string[n];
        for (var b = 0; b < n; b++) texts[b] = string.Join(" ", cols.Select(c => b < c.Length ? c[b] : "").Where(x => x.Length > 0));
        _bondTexts = (pairs, texts, crossing);
        _labelDoc = _doc;
        _labelFrame = _frame;
        return _bondTexts;
    }

    /// <summary>Labels to draw after a render: visible atoms (and bonds, at their midpoints), in view pixels (points).</summary>
    public List<ViewLabel> ViewLabels(CapsCamera cam, CapsRenderOpts opt, double scaling)
    {
        if (_doc == null || !AnyLabels || !IsStudio) { LabelNote = ""; return new List<ViewLabel>(); }
        return LabelsFor(_doc, cam, opt, scaling, true);
    }

    // the view window's document (a copy of the structure, or a run's live frames): its own texts, fetched per document
    private CapsDocument? _otherDoc;
    private (string Kinds, string[]? Atoms, uint[]? Rgb, (int[] Pairs, string[] Texts, bool[] Crossing)? Bonds) _otherLabels;
    /// <summary>Labels for another view of the structure (the view window) with the same settings.</summary>
    public List<ViewLabel> ViewLabelsFor(CapsDocument doc, CapsCamera cam, CapsRenderOpts opt, double scaling) =>
        AnyLabels ? LabelsFor(doc, cam, opt, scaling, false) : new List<ViewLabel>();

    private List<ViewLabel> LabelsFor(CapsDocument doc, CapsCamera cam, CapsRenderOpts opt, double scaling, bool main)
    {
        var list = new List<ViewLabel>();
        string[]? texts;
        uint[]? rgb;
        (int[] Pairs, string[] Texts, bool[] Crossing)? bonds;
        if (main)
        {
            texts = AtomLabelKinds.Any(k => k.On) ? LabelTexts() : null;
            rgb = _labelRgb;
            bonds = BondTexts();
        }
        else
        {
            var key = string.Join(",", AtomLabelKinds.Concat(BondLabelKinds).Where(k => k.On).Select(k => (k.Bond ? "b:" : "a:") + k.Id)) + "|" + _atomLabelColour;
            if (!ReferenceEquals(_otherDoc, doc) || _otherLabels.Kinds != key)
            {
                var (d0, t0, r0, b0) = (_doc, _labelTexts, _labelRgb, _bondTexts);
                var (dd, ff) = (_labelDoc, _labelFrame);
                // the same code as the main view's, on the other document
                _doc = doc; _labelTexts = null; _bondTexts = null; _labelDoc = null;
                try { _otherLabels = (key, AtomLabelKinds.Any(k => k.On) ? LabelTexts() : null, _labelRgb, BondTexts()); }
                finally { (_doc, _labelTexts, _labelRgb, _bondTexts, _labelDoc, _labelFrame) = (d0, t0, r0, b0, dd, ff); }
                _otherDoc = doc;
            }
            (texts, rgb, bonds) = (_otherLabels.Atoms, _otherLabels.Rgb, _otherLabels.Bonds);
        }
        var atoms = (int)doc.Summary().Atoms;
        var p = doc.ProjectAtoms(cam, opt, atoms);
        // which atoms may carry labels
        bool[]? allow = null;
        if (_labelScope == 1)
        {
            allow = new bool[atoms];
            try { foreach (var x in JsonNode.Parse((_doc ?? doc).SelectionJson())!["indices"]!.AsArray()) { var i = (int)x!.GetValue<double>(); if (i >= 0 && i < atoms) allow[i] = true; } } catch { }
        }
        else if (_labelScope == 2)
        {
            try { allow = JsonNode.Parse(doc.AtomLabels("element"))!.AsArray().Select(x => (x?.GetValue<string>() ?? "") != "H").ToArray(); } catch { }
        }
        bool Ok(int i) => i >= 0 && i < atoms && p[3 * i + 2] > 0 && (allow == null || (i < allow.Length && allow[i]));
        if (texts != null)
            for (var i = 0; i < Math.Min(atoms, texts.Length); i++)
                if (texts[i].Length > 0 && Ok(i)) list.Add(new ViewLabel(p[3 * i] / scaling, p[3 * i + 1] / scaling, texts[i], rgb != null && i < rgb.Length ? rgb[i] : 0));
        if (bonds is { } bt)
            for (var b = 0; b < bt.Texts.Length; b++)
            {
                var (i, j) = (bt.Pairs[2 * b], bt.Pairs[2 * b + 1]);
                if (bt.Texts[b].Length == 0 || bt.Crossing[b] || !Ok(i) || !Ok(j)) continue;   // a bond across the cell has no midpoint on screen
                list.Add(new ViewLabel((p[3 * i] + p[3 * j]) / 2 / scaling, (p[3 * i + 1] + p[3 * j + 1]) / 2 / scaling, bt.Texts[b], 0, true));
            }
        if (list.Count > LabelLimit) { if (main) LabelNote = $"{list.Count:N0} labels: zoom in, label fewer kinds or only the selection (up to {LabelLimit:N0})"; return new List<ViewLabel>(); }
        if (main) LabelNote = "";
        return list;
    }

    // ---------------------------------------------------------------- surfaces

    private int _appSurface, _appMap, _appOpacity = 2;
    private decimal _appProbe = 1.4m;
    private string _appSurfaceExpr = "";
    public int AppSurface { get => _appSurface; set { if (Set(ref _appSurface, Math.Clamp(value, 0, 3))) { Raise(nameof(AppHasSurface)); Raise(nameof(AppRampEnabled)); Raise(nameof(ShowAppLegend)); ApplyAppearance(); } } }
    public bool AppHasSurface => _appSurface > 0;
    public decimal AppProbe { get => _appProbe; set { if (Set(ref _appProbe, Math.Clamp(value, 0.5m, 5m))) ApplyAppearance(); } }
    public int AppSurfaceMap { get => _appMap; set { if (Set(ref _appMap, Math.Clamp(value, 0, 2))) { Raise(nameof(ShowAppLegend)); ApplyAppearance(); } } }
    public int AppOpacity { get => _appOpacity; set { if (Set(ref _appOpacity, Math.Clamp(value, 0, 4))) ApplyAppearance(); } }
    /// <summary>The atoms the surface wraps (an expression; empty: all).</summary>
    public string AppSurfaceAtoms { get => _appSurfaceExpr; set { if (Set(ref _appSurfaceExpr, value ?? "")) ApplyAppearance(); } }

    // ---------------------------------------------------------------- rendering

    private bool _viewAo;
    public bool ViewAo { get => _viewAo; set { if (Set(ref _viewAo, value)) RenderRequested?.Invoke(); } }

    // ---------------------------------------------------------------- legend and info

    private string _appInfo = "", _appLegendTitle = "", _appLegendLo = "", _appLegendHi = "", _appError = "";
    public string AppInfo { get => _appInfo; private set => Set(ref _appInfo, value); }
    public string AppError { get => _appError; private set { if (Set(ref _appError, value)) Raise(nameof(AppHasError)); } }
    public bool AppHasError => _appError.Length > 0;
    public bool ShowAppLegend => IsStudio && _doc != null && (AppColour == 4 || (_appSurface > 0 && _appMap == 0));
    public string AppLegendTitle { get => _appLegendTitle; private set => Set(ref _appLegendTitle, value); }
    public string AppLegendLo { get => _appLegendLo; private set => Set(ref _appLegendLo, value); }
    public string AppLegendHi { get => _appLegendHi; private set => Set(ref _appLegendHi, value); }
    public Avalonia.Media.IBrush AppLegendBrush
    {
        get
        {
            var stops = _appRamp switch
            {
                1 => new[] { "#440154", "#3B528B", "#21918C", "#5EC962", "#FDE725" },
                2 => new[] { "#D6604D", "#F4F4F4", "#4393C3" },
                _ => new[] { "#3B7DD8", "#D9DCDF", "#E8893A" },
            };
            var g = new Avalonia.Media.LinearGradientBrush { StartPoint = new Avalonia.RelativePoint(0, 0, Avalonia.RelativeUnit.Relative), EndPoint = new Avalonia.RelativePoint(1, 0, Avalonia.RelativeUnit.Relative) };
            for (var k = 0; k < stops.Length; k++) g.GradientStops.Add(new Avalonia.Media.GradientStop(Avalonia.Media.Color.Parse(stops[k]), k / (double)(stops.Length - 1)));
            return g;
        }
    }

    private int _appTicket;
    private readonly object _appLock = new();
    private Task<(string Info, string? Error, bool Stale)>? _appTask;
    /// <summary>Waits until the last appearance change has reached the core (tests, exports).</summary>
    public void WaitAppearance() { try { _appTask?.Wait(); } catch { } }
    /// <summary>Sends the appearance to the core (surfaces are computed in the background) and redraws.</summary>
    public void ApplyAppearance()
    {
        Raise(nameof(DisplayStatus));
        if (_doc == null) return;
        var doc = _doc;
        var json = new JsonObject
        {
            ["active"] = AppLayers.Count > 0 || AppColour == 4 || _appSurface > 0 || ColumnName != null,
            ["layers"] = new JsonArray(AppLayers.Select(l => (JsonNode)new JsonObject { ["expression"] = l.Expression, ["style"] = AppStyleCodes[l.Style] }).ToArray()),
            ["colour"] = ColumnName is { } col ? "column:" + col : AppColour == 4 ? "charge" : "",
            ["charges"] = AppColour == 4 && _chgPreview is { } pq ? new JsonArray(pq.Select(x => (JsonNode)x).ToArray()) : null,
            ["ramp"] = AppRampCodes[_appRamp],
            ["surface"] = new JsonObject
            {
                ["kind"] = AppSurfaceCodes[_appSurface], ["probe"] = (double)_appProbe, ["opacity"] = (_appOpacity + 1) * 0.2,
                ["expression"] = _appSurfaceExpr.Trim(), ["colour"] = AppMapCodes[_appMap], ["spacing"] = 0.6,
            },
        }.ToJsonString();
        var ticket = ++_appTicket;
        AppInfo = _appSurface > 0 ? "Computing the surface…" : AppInfo;
        // one apply at a time, and a stale one (a newer change was made meanwhile) is skipped inside the lock: without
        // this, two quick changes could reach the core out of order and leave the older settings drawn
        _appTask = Task.Run(() =>
        {
            lock (_appLock)
            {
                if (ticket != Volatile.Read(ref _appTicket)) return (Info: "", Error: (string?)null, Stale: true);
                try { doc.SetAppearance(json); return (Info: doc.AppearanceInfo(), Error: (string?)null, Stale: false); }
                catch (Exception e) { return (Info: "", Error: (string?)e.Message, Stale: false); }
            }
        });
        _appTask.ContinueWith(t => Avalonia.Threading.Dispatcher.UIThread.Post(() =>
        {
            if (ticket != _appTicket || t.Result.Stale) return;
            AppError = t.Result.Error ?? "";
            if (t.Result.Error == null) ShowAppearanceInfo(t.Result.Info);
            RenderRequested?.Invoke();
        }));
    }

    private void AppRefreshInfo()
    {
        if (_doc == null) return;
        try { ShowAppearanceInfo(_doc.AppearanceInfo()); } catch { }
    }

    private void ShowAppearanceInfo(string text)
    {
        var j = JsonNode.Parse(text)!;
        var parts = new List<string>();
        var n = _doc?.Summary().Atoms ?? 0;
        parts.Add($"{n:N0} atoms");
        if (j["styles"] is JsonObject st)
            foreach (var (k, v) in st)
                if (k != "view" && v != null) parts.Add($"{v.GetValue<double>():N0} {k.Replace('_', ' ').Replace("sticks", "stick")}");
        if (j["surface"] is JsonObject sf)
        {
            parts.Add(string.Format(CultureInfo.InvariantCulture, "{0} {1:N0} Å²", AppSurfaceNames[_appSurface].ToLowerInvariant(), sf["area"]!.GetValue<double>()));
            if (_appMap == 0 && sf["potential"] is JsonArray pot)
            {
                AppLegendTitle = "Electrostatic potential (kcal/mol/e, 95th percentile) · colour-blind-safe diverging";
                AppLegendLo = pot[0]!.GetValue<double>().ToString("+0.0;−0.0;0", CultureInfo.InvariantCulture);
                AppLegendHi = pot[1]!.GetValue<double>().ToString("+0.0;−0.0;0", CultureInfo.InvariantCulture);
            }
        }
        if (j["polyhedra"] != null) parts.Add($"{j["polyhedra"]!.GetValue<double>():N0} polyhedron faces");
        if (AppColour == 4 && !(_appSurface > 0 && _appMap == 0) && j["charge"] is JsonArray q)
        {
            var m = Math.Max(Math.Abs(q[0]!.GetValue<double>()), Math.Abs(q[1]!.GetValue<double>()));
            var preview = j["preview"]?.GetValue<bool>() == true;
            AppLegendTitle = m < 1e-9 ? "Partial charge (e) · no charges in the file: assign them in Field"
                           : preview ? "Partial charge (e) · computed here, not applied yet · colour-blind-safe diverging" : "Partial charge (e) · colour-blind-safe diverging";
            AppLegendLo = (-m).ToString("+0.00;−0.00;0", CultureInfo.InvariantCulture);
            AppLegendHi = m.ToString("+0.00;−0.00;0", CultureInfo.InvariantCulture);
        }
        AppInfo = string.Join(" · ", parts);
        Raise(nameof(ShowAppLegend));
    }

    /// <summary>Clears every layer, surface and label (back to the view's style).</summary>
    public void ResetAppearance()
    {
        AppLayers.Clear();
        _appColour = -1;
        _appSurface = 0;
        foreach (var k in AtomLabelKinds.Concat(BondLabelKinds)) k.On = false;
        foreach (var n in new[] { nameof(AppColour), nameof(AppSurface), nameof(AppHasSurface), nameof(LabelElement), nameof(LabelRs), nameof(LabelType), nameof(LabelCharge), nameof(AppChip) }) Raise(n);
        _labelTexts = null;
        _bondTexts = null;
        ApplyAppearance();
    }
}
