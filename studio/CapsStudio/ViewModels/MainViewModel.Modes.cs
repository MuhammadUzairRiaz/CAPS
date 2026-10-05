using System.Globalization;
using System.Text.Json.Nodes;

namespace CapsStudio.ViewModels;

/// <summary>Normal modes (Analyze › Normal modes): a mode played in the view — one period in a copy of the structure, so
/// the structure and its own frames stay as they are.</summary>
public sealed partial class MainViewModel
{
    /// <summary>The conformer search again (the same settings and seed) on a copy, whose frames become the conformers.</summary>
    public async Task OpenConformers()
    {
        if (_doc == null || !Idle) return;
        var temperature = (double)Analyze.FluctTD;
        var json = Analyze.ConfJson(temperature);
        DuplicateStructure(" · conformers");
        var doc = _doc!;
        Status = "Conformer search…";
        try
        {
            var r = JsonNode.Parse(await Task.Run(() => doc.ConformerFrames(json)))!;
            AfterRun(doc, "");
            var cs = r["conformers"] as JsonArray ?? [];
            var inv = CultureInfo.InvariantCulture;
            Status = string.Format(inv, "{0} conformers from {1:0} starts ({2:0} rotatable bonds) · frame 1 the lowest{3}", cs.Count, (double?)r["trials"] ?? 0, (double?)r["rotors"] ?? 0,
                cs.Count > 1 ? string.Format(inv, ", frame 2 +{0:0.00} kcal/mol", (double?)cs[1]?["relative"] ?? 0) : "");
            Frame = 0;
        }
        catch (Exception e) { Status = "Could not search conformers: " + e.Message; }
    }

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
