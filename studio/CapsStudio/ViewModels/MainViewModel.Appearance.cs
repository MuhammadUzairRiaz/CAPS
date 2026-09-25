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
public readonly record struct ViewLabel(double X, double Y, string Text);

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
    public bool AppearanceOpen { get => _appOpen; set { if (Set(ref _appOpen, value)) { if (value) SelectionOpen = false; RaiseAppearanceVisibility(); if (value) AppRefreshInfo(); } } }
    public bool ShowStudioTabs => IsStudio && !_appOpen && !_selOpen;
    public bool ShowAppearance => IsStudio && _appOpen && _doc != null;
    private void RaiseAppearanceVisibility() { Raise(nameof(ShowStudioTabs)); Raise(nameof(ShowAppearance)); Raise(nameof(ShowAppLegend)); Raise(nameof(ShowSelectionPanel)); }

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
    public int AppRamp { get => _appRamp; set { if (Set(ref _appRamp, Math.Clamp(value, 0, AppRampNames.Length - 1))) { Raise(nameof(AppLegendBrush)); ApplyAppearance(); } } }
    public bool AppRampEnabled => AppColour == 4 || _appSurface > 0;

    // ---------------------------------------------------------------- labels

    private bool _lblElement, _lblRs, _lblType, _lblCharge;
    public bool LabelElement { get => _lblElement; set { if (Set(ref _lblElement, value)) LabelsChanged(); } }
    public bool LabelRs { get => _lblRs; set { if (Set(ref _lblRs, value)) LabelsChanged(); } }
    public bool LabelType { get => _lblType; set { if (Set(ref _lblType, value)) LabelsChanged(); } }
    public bool LabelCharge { get => _lblCharge; set { if (Set(ref _lblCharge, value)) LabelsChanged(); } }
    public bool AnyLabels => _lblElement || _lblRs || _lblType || _lblCharge;
    private void LabelsChanged() { _labelTexts = null; RenderRequested?.Invoke(); }
    private string[]? _labelTexts;
    private object? _labelDoc;
    private int _labelFrame = -1;
    /// <summary>Most atoms labelled at once; above this the labels wait for a smaller view.</summary>
    public const int LabelLimit = 4000;
    private string _labelNote = "";
    public string LabelNote { get => _labelNote; private set => Set(ref _labelNote, value); }

    /// <summary>The label texts per atom (joined kinds), fetched once per document and frame.</summary>
    private string[]? LabelTexts()
    {
        if (_doc == null || !AnyLabels) return null;
        if (_labelTexts != null && ReferenceEquals(_labelDoc, _doc) && _labelFrame == _frame) return _labelTexts;
        var kinds = new List<string>();
        if (_lblElement) kinds.Add("element");
        if (_lblRs) kinds.Add("rs");
        if (_lblType) kinds.Add("type");
        if (_lblCharge) kinds.Add("charge");
        var cols = kinds.Select(k => JsonNode.Parse(_doc.AtomLabels(k))!.AsArray().Select(x => x?.GetValue<string>() ?? "").ToArray()).ToList();
        var n = cols.Count > 0 ? cols[0].Length : 0;
        var texts = new string[n];
        for (var i = 0; i < n; i++)
        {
            var parts = cols.Select(c => i < c.Length ? c[i] : "").Where(x => x.Length > 0).ToList();
            // R/S alone labels only the stereocentres
            texts[i] = kinds.Count == 1 && kinds[0] == "rs" ? parts.FirstOrDefault() ?? "" : parts.Count == 0 ? "" : string.Join(" ", parts);
        }
        _labelTexts = texts;
        _labelDoc = _doc;
        _labelFrame = _frame;
        return texts;
    }

    /// <summary>Labels to draw after a render: visible atoms only, in view pixels (points).</summary>
    public List<ViewLabel> ViewLabels(CapsCamera cam, CapsRenderOpts opt, double scaling)
    {
        var list = new List<ViewLabel>();
        if (_doc == null || !AnyLabels || !IsStudio) { LabelNote = ""; return list; }
        var texts = LabelTexts();
        if (texts == null) return list;
        var atoms = texts.Length;
        var p = _doc.ProjectAtoms(cam, opt, atoms);
        var visible = 0;
        for (var i = 0; i < atoms; i++) if (p[3 * i + 2] > 0 && texts[i].Length > 0) visible++;
        if (visible > LabelLimit) { LabelNote = $"{visible:N0} labels: zoom in or label fewer kinds (up to {LabelLimit:N0})"; return list; }
        LabelNote = "";
        for (var i = 0; i < atoms; i++)
            if (p[3 * i + 2] > 0 && texts[i].Length > 0) list.Add(new ViewLabel(p[3 * i] / scaling, p[3 * i + 1] / scaling, texts[i]));
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
    /// <summary>Sends the appearance to the core (surfaces are computed in the background) and redraws.</summary>
    public void ApplyAppearance()
    {
        if (_doc == null) return;
        var doc = _doc;
        var json = new JsonObject
        {
            ["active"] = AppLayers.Count > 0 || AppColour == 4 || _appSurface > 0,
            ["layers"] = new JsonArray(AppLayers.Select(l => (JsonNode)new JsonObject { ["expression"] = l.Expression, ["style"] = AppStyleCodes[l.Style] }).ToArray()),
            ["colour"] = AppColour == 4 ? "charge" : "",
            ["ramp"] = AppRampCodes[_appRamp],
            ["surface"] = new JsonObject
            {
                ["kind"] = AppSurfaceCodes[_appSurface], ["probe"] = (double)_appProbe, ["opacity"] = (_appOpacity + 1) * 0.2,
                ["expression"] = _appSurfaceExpr.Trim(), ["colour"] = AppMapCodes[_appMap], ["spacing"] = 0.6,
            },
        }.ToJsonString();
        var ticket = ++_appTicket;
        AppInfo = _appSurface > 0 ? "Computing the surface…" : AppInfo;
        Task.Run(() =>
        {
            try { doc.SetAppearance(json); return (Info: doc.AppearanceInfo(), Error: (string?)null); }
            catch (Exception e) { return (Info: "", Error: (string?)e.Message); }
        }).ContinueWith(t => Avalonia.Threading.Dispatcher.UIThread.Post(() =>
        {
            if (ticket != _appTicket) return;
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
            AppLegendTitle = m < 1e-9 ? "Partial charge (e) · no charges in the file: assign them in Field" : "Partial charge (e) · colour-blind-safe diverging";
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
        _lblElement = _lblRs = _lblType = _lblCharge = false;
        foreach (var n in new[] { nameof(AppColour), nameof(AppSurface), nameof(AppHasSurface), nameof(LabelElement), nameof(LabelRs), nameof(LabelType), nameof(LabelCharge), nameof(AppChip) }) Raise(n);
        _labelTexts = null;
        ApplyAppearance();
    }
}
