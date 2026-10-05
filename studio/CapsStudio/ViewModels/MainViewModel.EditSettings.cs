using System.Globalization;
using System.Text.Json.Nodes;
using CapsStudio.Interop;

namespace CapsStudio.ViewModels;

/// <summary>Edit in its builder: the page comes back with the settings that made the structure — kept with the structure
/// when it was built (and in the project file), else read from its provenance (structures grown before they were kept).</summary>
public sealed partial class MainViewModel
{
    private static readonly System.Text.Json.JsonSerializerOptions OneLine = new() { WriteIndented = false };

    /// <summary>Everything the Grow page holds, as the structure's build settings.</summary>
    private JsonObject GrowSettingsJson()
    {
        var o = new JsonObject
        {
            ["page"] = "grow", ["spec"] = _growSpec, ["spec_name"] = _growSpecName, ["spec_atoms"] = _growSpecAtoms, ["spec_mass"] = _growSpecMass,
            ["chains"] = _growChains, ["dp"] = _growDp, ["tacticity"] = _growTact, ["seed"] = _growSeed,
            ["density"] = _growDensity, ["use_box"] = _growUseBox, ["box"] = _growBox, ["shape"] = _growShape,
            ["slab_h"] = _growSlabH, ["slab_vacuum"] = _growSlabVac, ["cyl_r"] = _growCylR, ["cyl_len"] = _growCylLen,
            ["scale"] = _growScale, ["auto_scale"] = _growAutoScale, ["curve"] = _growCurve,
            ["method"] = _growMethod, ["method_t"] = _growMethodT, ["lookahead"] = _growLookahead, ["orient"] = _growOrient, ["orient_s"] = _growOrientS,
            ["ff"] = GrowFfIndex >= 0 && GrowFfIndex < Field.Library.Count ? Field.Library[GrowFfIndex].Id : null,
            ["assign"] = _growAssignFf, ["charges"] = _growCharges,
        };
        return o;
    }

    /// <summary>The Grow page set from a structure's build settings.</summary>
    private void ApplyGrowSettings(JsonObject o)
    {
        int I(string k, int d) => o[k] is JsonValue v && v.TryGetValue<int>(out var x) ? x : o[k] is JsonValue w && w.TryGetValue<double>(out var y) ? (int)y : d;
        double D(string k, double d) => o[k] is JsonValue v && v.TryGetValue<double>(out var x) ? x : d;
        bool B(string k, bool d) => o[k] is JsonValue v && v.TryGetValue<bool>(out var x) ? x : d;
        _growSpec = (string?)o["spec"];
        _growSpecName = (string?)o["spec_name"] ?? (_growSpec == null ? "PS" : "polymer");
        _growSpecAtoms = I("spec_atoms", 0);
        _growSpecMass = D("spec_mass", 0);
        GrowChains = I("chains", _growChains);
        GrowDp = I("dp", _growDp);
        GrowTacticity = I("tacticity", _growTact);
        GrowSeed = I("seed", _growSeed);
        GrowDensity = D("density", _growDensity);
        GrowShape = I("shape", 0);
        GrowUseBox = B("use_box", false);
        GrowBox = D("box", _growBox);
        _growSlabH = D("slab_h", _growSlabH); _growSlabVac = D("slab_vacuum", _growSlabVac); _growCylR = D("cyl_r", _growCylR); _growCylLen = D("cyl_len", _growCylLen);
        GrowScale = D("scale", _growScale);
        GrowAutoScale = B("auto_scale", true);
        GrowCurve = B("curve", true);
        GrowMethod = I("method", 0);
        _growMethodT = D("method_t", _growMethodT);
        _growLookahead = I("lookahead", 1);
        _growOrient = I("orient", 0);
        _growOrientS = D("orient_s", _growOrientS);
        if ((string?)o["ff"] is { Length: > 0 } ff && Field.Library.ToList().FindIndex(x => x.Id == ff) is var k and >= 0) GrowFfIndex = k;
        GrowAssignField = B("assign", _growAssignFf);
        GrowChargeMode = I("charges", _growCharges);
        foreach (var n in new[] { nameof(GrowHasSpec), nameof(GrowComponentName), nameof(GrowAtomsText), nameof(GrowEstimate), nameof(GrowCommand), nameof(GrowChainsD),
                                  nameof(GrowDpD), nameof(GrowSeedD), nameof(GrowDensityD), nameof(GrowBoxD), nameof(GrowScaleD), nameof(GrowRegionAD), nameof(GrowRegionBD),
                                  nameof(GrowMethodTempD), nameof(GrowOrientStrengthD), nameof(GrowLookaheadD) }) Raise(n);
    }

    /// <summary>A grown structure's Grow settings from its provenance (its grow.trials step and the force field assigned
    /// after it); what the provenance does not record — tacticity and seed — is read from the name CAPS gave the cell.</summary>
    private static JsonObject? GrowSettingsFromProvenance(CapsDocument doc, string name, IReadOnlyList<FfEntry> library)
    {
        JsonArray? steps;
        try { steps = JsonNode.Parse(doc.Provenance())?["steps"] as JsonArray; } catch { return null; }
        // a step's params: an object, or (as some steps write them) a list of [name, value] pairs
        static JsonObject? Params(JsonNode? step) => step?["params"] switch
        {
            JsonObject o => o,
            JsonArray a => new JsonObject(a.OfType<JsonArray>().Where(x => x.Count >= 2 && x[0] is JsonValue)
                                           .Select(x => new KeyValuePair<string, JsonNode?>((string?)x[0] ?? "", x[1]?.DeepClone()))
                                           .DistinctBy(kv => kv.Key)),
            _ => null,
        };
        var grow = steps?.FirstOrDefault(s => (string?)s?["engine"] == "grow.trials");
        if (Params(grow) is not JsonObject p) return null;
        string? S(string k) => p[k] is JsonValue v && v.TryGetValue<string>(out var s) ? s : p[k]?.ToJsonString();
        double? Num(string k) => S(k) is { } s && double.TryParse(s.Split(' ')[0], NumberStyles.Float, CultureInfo.InvariantCulture, out var x) ? x : null;
        JsonNode? Json(string k) { try { return S(k) is { } s ? JsonNode.Parse(s) : null; } catch { return null; } }
        var o = new JsonObject { ["page"] = "grow" };
        if (Json("units") is JsonArray units)
        {
            var spec = new JsonObject { ["units"] = units.DeepClone(), ["sequence"] = S("sequence") ?? "homopolymer" };
            if (Num("dp") is { } dp) spec["dp"] = (int)dp;
            foreach (var k in new[] { "weights", "blocks" }) if (Json(k) is JsonArray a) spec[k] = a.DeepClone();
            foreach (var k in new[] { "pattern", "head_cap", "tail_cap", "linkage", "architecture" }) if (S(k) is { } v) spec[k] = v;
            o["spec"] = spec.ToJsonString(OneLine);
            o["spec_name"] = units.Count == 1 ? (string?)units[0]?["name"] ?? "polymer" : "copolymer";
        }
        if (Num("dp") is { } d) o["dp"] = (int)d;
        if (Num("chains") is { } c) o["chains"] = (int)c;
        if (S("density") is { } dens && Num("density") is { } dv) { o["density"] = dv; o["use_box"] = !dens.Contains("g/cm", StringComparison.Ordinal) && dens.Contains('Å'); if ((bool)o["use_box"]!) o["box"] = dv; }
        if (S("box") is not null && Num("box") is { } bx) { o["use_box"] = true; o["box"] = bx; }
        if (S("contact scale") is { } cs) { o["auto_scale"] = cs.Contains("lowered", StringComparison.Ordinal) || cs.Contains("auto", StringComparison.Ordinal); if (Num("contact scale") is { } sc) o["scale"] = sc; }
        var lower = name.ToLowerInvariant();
        o["tacticity"] = lower.Contains("isotactic") ? 1 : lower.Contains("syndiotactic") ? 2 : 0;
        var m = System.Text.RegularExpressions.Regex.Match(name, @"seed(\d+)");
        if (m.Success) o["seed"] = int.Parse(m.Groups[1].Value, CultureInfo.InvariantCulture);
        // the force field: the one assigned after growing, else the one Grow typed with
        // (a force field read back whole from a saved .ff.json names no library file: the assignment before it does)
        FfEntry? Lib(string? f) => f == null ? null : library.FirstOrDefault(e => string.Equals(Path.GetFileName(e.File), Path.GetFileName(f), StringComparison.OrdinalIgnoreCase));
        var assigned = steps!.Where(s => (string?)s?["engine"] == "field.assign").Select(s => Params(s)?["force field"])
                             .Select(x => x is JsonValue v && v.TryGetValue<string>(out var f) ? Lib(f) : null).LastOrDefault(e => e != null);
        if ((assigned ?? Lib(S("forcefield"))) is { } fe) { o["ff"] = fe.Id; o["assign"] = assigned != null || (bool?)o["assign"] == true; }
        return o;
    }
}
