using System;
using System.Globalization;
using System.Linq;
using System.Text.Json.Nodes;
using System.Threading.Tasks;
using CapsStudio.Interop;

namespace CapsStudio.ViewModels;

/// <summary>Analyze › Orientation (design/boards/Orientation): the nematic order parameter S of backbone chords (i → i+2),
/// its director, Herman's f along z, local crystallinity, S per frame and P₂ along the cell's z — beside an estimator
/// check on Bunn's orthorhombic polyethylene, where every chord lies on c and S must be exactly 1.</summary>
public sealed partial class MainViewModel
{
    public bool IsOrientation => _module == 52;
    public ResultCell OrS { get; } = new("S · this structure", "largest eigenvalue of Q");
    public ResultCell OrCheck { get; } = new("S · PE crystal (Bunn), 3 × 4 × 8 cells", "estimator check");
    public ResultCell OrHerman { get; } = new("Herman's f along z", "⟨P₂(cos θ)⟩ against z");
    public ResultCell OrCryst { get; } = new("Local crystallinity", "chords with ≥ 8 neighbours within 5 Å aligned within 10°");
    public ResultCell OrDirector { get; } = new("Director", "eigenvector of λ max");
    private string _orStatus = "";
    public string OrStatus { get => _orStatus; private set => Set(ref _orStatus, value); }
    public (double X, double Y)[] OrPerFrame { get; private set; } = [];
    public (double X, double Y)[] OrAlongZ { get; private set; } = [];
    public bool OrHasFrames => OrPerFrame.Length > 1;
    public bool OrHasZ => OrAlongZ.Length > 1;
    public event Action? OrChanged;

    public void OpenOrientation()
    {
        SetModule(52);
        if (!OrCheck.HasValue) _ = Task.Run(OrientationCheck).ContinueWith(t => Avalonia.Threading.Dispatcher.UIThread.Post(() =>
        {
            var (s, chords) = t.IsCompletedSuccessfully ? t.Result : (double.NaN, 0);
            OrCheck.Value = double.IsFinite(s) ? s.ToString("0.000", CultureInfo.InvariantCulture) : "—";
            OrCheck.Caption = $"S from {chords} chord vectors · chords i → i+2 lie on the c axis, so S = 1 exactly: a check that the estimator is right";
        }));
        OrChanged?.Invoke();
    }

    /// <summary>Bunn's polyethylene built and analysed off the UI thread (a few milliseconds).</summary>
    public static (double S, int Chords) OrientationCheck()
    {
        var spec = new JsonObject
        {
            ["space_group"] = "Pnam", ["a"] = 7.40, ["b"] = 4.93, ["c"] = 2.534, ["supercell"] = new JsonArray(3, 4, 8),
            ["sites"] = new JsonArray(
                new JsonObject { ["label"] = "C1", ["element"] = "C", ["x"] = 0.0380, ["y"] = 0.0650, ["z"] = 0.25 },
                new JsonObject { ["label"] = "H1", ["element"] = "H", ["x"] = 0.1848, ["y"] = 0.0466, ["z"] = 0.25 },
                new JsonObject { ["label"] = "H2", ["element"] = "H", ["x"] = 0.0068, ["y"] = 0.2811, ["z"] = 0.25 }),
        };
        var (pe, _) = CapsDocument.CrystalBuild(spec.ToJsonString(), "PE");
        using (pe)
        {
            var o = new CapsAnalyzeOpts { First = 0, Last = -1, Stride = 1, Blocks = 5 };
            var j = JsonNode.Parse(pe.Analyze("orientation", o, null))!;
            var p = ((JsonArray)j["properties"]!)[0]!;
            return (p["value"]!.GetValue<double>(), (int)(p["extra"]?["chord vectors per frame"]?.GetValue<double>() ?? 0));
        }
    }

    public async Task RunOrientation()
    {
        if (_doc == null || Analyze.Working) return;
        var inv = CultureInfo.InvariantCulture;
        OrStatus = "Chord vectors…";
        await RunChips("orientation");
        var card = Analyze.Results.FirstOrDefault(r => r.Id == "orientation");
        if (card == null || !double.IsFinite(card.Value)) { OrStatus = card?.Notes.FirstOrDefault() ?? Analyze.Log; OrChanged?.Invoke(); return; }
        double X(string k) => card.Extra.Where(e => e.Key.StartsWith(k, StringComparison.Ordinal)).Select(e => e.Value).DefaultIfEmpty(double.NaN).First();
        OrS.Value = card.Value.ToString("0.000", inv) + (double.IsFinite(card.Error) && card.Error > 0 ? " ± " + card.Error.ToString("0.000", inv) : "");
        OrS.Caption = $"S from {X("chord vectors").ToString("0", inv)} chord vectors per frame" + (card.Notes.Length > 0 ? " · " + card.Notes[0] : "");
        OrHerman.Value = X("Herman f").ToString("0.000", inv);
        OrCryst.Value = (X("local crystallinity") * 100).ToString("0.0", inv) + " %";
        OrDirector.Value = $"({X("director x").ToString("0.00", inv)}, {X("director y").ToString("0.00", inv)}, {X("director z").ToString("0.00", inv)})";
        var curves = Analyze.Curves.Where(c => c.Property == card.Name).ToList();
        OrPerFrame = curves.FirstOrDefault(c => c.Label == "S per frame") is { } f ? f.X.Zip(f.Y).ToArray() : [];
        OrAlongZ = curves.FirstOrDefault(c => c.Label.StartsWith("P₂")) is { } z ? z.X.Zip(z.Y).ToArray() : [];
        OrStatus = Analyze.Log;
        Raise(nameof(OrHasFrames)); Raise(nameof(OrHasZ));
        OrChanged?.Invoke();
    }
}
