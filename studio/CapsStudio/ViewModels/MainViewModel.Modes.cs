using System.Globalization;
using System.Text.Json.Nodes;

namespace CapsStudio.ViewModels;

/// <summary>Normal modes (Analyze › Normal modes): a mode played in the view — one period in a copy of the structure, so
/// the structure and its own frames stay as they are.</summary>
public sealed partial class MainViewModel
{
    public async Task PlayMode(int mode, double amplitude = 0.3)
    {
        if (_doc == null || !Idle) return;
        DuplicateStructure($" · mode {mode}");
        var doc = _doc!;
        Status = $"Normal mode {mode}: the Hessian of {Title.Replace(" (unsaved)", "")}…";
        try
        {
            var json = await Task.Run(() => doc.ModeAnimate(new JsonObject { ["mode"] = mode, ["amplitude"] = amplitude, ["frames"] = 24 }.ToJsonString()));
            var r = JsonNode.Parse(json)!;
            AfterRun(doc, "");
            IsPlaying = true;
            var inv = CultureInfo.InvariantCulture;
            var nu = (double?)r["wavenumber"] ?? double.NaN;
            Status = string.Format(inv, "Mode {0} of {1:0}: {2:0.0} cm⁻¹{3} · reduced mass {4:0.00} g/mol · one period in {5:0} frames, ±{6:0.##} Å",
                mode, (double?)r["modes"] ?? 0, Math.Abs(nu), nu < 0 ? " (imaginary)" : "", (double?)r["reduced_mass"] ?? 0, (double?)r["frames"] ?? 0, amplitude);
        }
        catch (Exception e) { Status = "Could not play the mode: " + e.Message; }
    }
}
