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

    // ---------------------------------------------------------------- Packing

    /// <summary>The Pack page as the structure's build settings: its input (molecules, counts, regions, each molecule's
    /// force field), where relative paths are read from, what it packed around, and the cell's force field.</summary>
    private JsonObject PackSettingsJson(string hostTitle) => new()
    {
        ["page"] = "pack", ["text"] = _packText, ["base"] = _packBaseDir, ["start"] = _packStart, ["host"] = hostTitle,
        ["ff"] = PackFfIndex >= 0 && PackFfIndex < Field.Library.Count ? Field.Library[PackFfIndex].Id : null,
        ["charges"] = _packCharges, ["assign"] = _packAssign, ["compress"] = _packCompress, ["compress_to"] = (double)_packCompressTo,
    };

    /// <summary>The Pack page set from a packed structure's settings; packed around a structure still in the project,
    /// Pack starts from it again and the new packing replaces this one.</summary>
    private void ApplyPackSettings(JsonObject o, ProjectItem packed)
    {
        if ((string?)o["base"] is { Length: > 0 } b && Directory.Exists(b)) PackBaseDir = b;
        if ((string?)o["text"] is { Length: > 0 } t) PackText = t;
        if ((string?)o["ff"] is { Length: > 0 } ff && Field.Library.ToList().FindIndex(x => x.Id == ff) is var k and >= 0) PackFfIndex = k;
        if (o["charges"] is JsonValue cv && cv.TryGetValue<int>(out var ch)) PackChargeMode = ch;
        if (o["assign"] is JsonValue av && av.TryGetValue<bool>(out var a)) PackAssignField = a;
        if (o["compress"] is JsonValue pv && pv.TryGetValue<bool>(out var c)) PackCompress = c;
        if (o["compress_to"] is JsonValue tv && tv.TryGetValue<double>(out var to)) PackCompressTo = (decimal)to;
        var host = (string?)o["host"] is { Length: > 0 } h ? ProjectItems.FirstOrDefault(x => x != packed && CleanName(x.Name) == h) : null;
        if ((int?)o["start"] == 1 && host != null)
        {
            _packResult = packed.Doc;
            _packHost = host.Doc;
            _packHostTitle = CleanName(host.Name);
            PackStart = 1;
        }
        else PackStart = 0;
        Raise(nameof(PackStartNote));
    }

    /// <summary>A packed structure's Pack input from its provenance (CAPS records it from this version on).</summary>
    private static JsonObject? PackSettingsFromProvenance(CapsDocument doc)
    {
        JsonArray? steps;
        try { steps = JsonNode.Parse(doc.Provenance())?["steps"] as JsonArray; } catch { return null; }
        var p = StepParams(steps?.LastOrDefault(s => (string?)s?["engine"] == "pack.lbfgs"));
        if (p?["input"] is not JsonValue v || !v.TryGetValue<string>(out var text) || text.Length == 0) return null;
        // packed around a structure: the run's input holds it as a fixed block (a temporary copy) and its cell; the page's
        // input is without them, and the structure is named in the block's comment
        var m = System.Text.RegularExpressions.Regex.Match(text,
            @"(?ms)^(?:cell\s[^\n]*|pbc\s[^\n]*)\n\n(?:molecule|structure)\s+\S*caps-pack\S*host_\d+\.data\s+#\s*(?<title>.*?), kept where it is\n.*?^end(?: structure)?\n");
        if (!m.Success) return new JsonObject { ["page"] = "pack", ["text"] = text };
        return new JsonObject { ["page"] = "pack", ["text"] = text.Remove(m.Index, m.Length), ["start"] = 1, ["host"] = m.Groups["title"].Value.Trim() };
    }

    // ---------------------------------------------------------------- Blend

    /// <summary>The Blend page as the structure's build settings: each component (its polymer, or the builder's chain
    /// itself), share, DP, density, chains and force field, and the cell's chains, density, start and growth.</summary>
    private JsonObject BlendSettingsJson() => new()
    {
        ["page"] = "blend", ["mode"] = _blendMode, ["chains"] = (double)_blendChains, ["density"] = (double)_blendDensity, ["morphology"] = _blendMorph,
        ["method"] = _blendGrowMethod, ["temperature"] = (double)_blendGrowTemp,
        ["rows"] = new JsonArray(BlendRows.Where(r => r.Polymer != null).Select(r => (JsonNode)new JsonObject
        {
            ["id"] = r.Polymer!.Id, ["name"] = r.Polymer.Name, ["spec"] = JsonNode.Parse(SpecOf(r.Polymer, (int)r.Dp)),
            ["weight"] = (double)r.Weight, ["dp"] = (double)r.Dp, ["density"] = (double)r.Density, ["count"] = (double)r.Count,
            ["ff"] = r.ForceField is { Id.Length: > 0 } f ? f.Id : null,
        }).ToArray()),
    };

    /// <summary>The Blend page set from a blend's settings: library polymers by their id, any other chain (the builder's,
    /// or one no longer in the library) as a chain of its own.</summary>
    private void ApplyBlendSettings(JsonObject o)
    {
        double D(JsonNode? n, double d) => n is JsonValue v && v.TryGetValue<double>(out var x) ? x : d;
        OpenBlend();
        BlendRows.Clear();
        BlendMode = (int)D(o["mode"], 0);
        foreach (var r in (o["rows"] as JsonArray ?? []).OfType<JsonObject>().Take(5))
        {
            var id = (string?)r["id"] ?? "";
            var entry = id.StartsWith("builder:", StringComparison.Ordinal) ? null : BlendLibrary.FirstOrDefault(e => e.Id == id);
            if (entry == null && r["spec"] is JsonObject spec)
            {
                var own = (JsonObject)spec.DeepClone();
                own.Remove("dp");
                var bid = $"builder:{_builtChains.Count + 1}";
                _builtChains[bid] = own;
                var units = own["units"] as JsonArray;
                entry = new LibraryEntry(bid, (string?)r["name"] ?? "chain", (string?)units?[0]?["smiles"] ?? "", false, (units?.Count ?? 1) > 1, null, true);
                BlendLibrary = [.. BlendLibrary, entry];
                Raise(nameof(BlendLibrary));
            }
            if (entry == null) continue;
            AddBlendRow(entry, (decimal)D(r["weight"], 50));
            var row = BlendRows[^1];
            row.Dp = (decimal)D(r["dp"], 20);
            row.Density = (decimal)D(r["density"], 1.0);
            row.Count = (decimal)D(r["count"], 8);
            if ((string?)r["ff"] is { Length: > 0 } ff && BlendForceFields.FirstOrDefault(e => e.Id == ff) is { } fe) row.ForceField = fe;
        }
        BlendChains = (decimal)D(o["chains"], 8);
        BlendDensity = (decimal)D(o["density"], 0.5);
        BlendMorph = (int)D(o["morphology"], 0);
        BlendGrowMethod = (int)D(o["method"], 0);
        BlendGrowTemp = (decimal)D(o["temperature"], 450);
        BlendRecount();
    }

    /// <summary>A blend's settings from its provenance (grow.blend records every component's chain and share).</summary>
    private static JsonObject? BlendSettingsFromProvenance(CapsDocument doc)
    {
        JsonArray? steps;
        try { steps = JsonNode.Parse(doc.Provenance())?["steps"] as JsonArray; } catch { return null; }
        var p = StepParams(steps?.FirstOrDefault(s => (string?)s?["engine"] == "grow.blend"));
        if (p == null) return null;
        JsonNode? J(string k) { try { return p[k] is JsonValue v && v.TryGetValue<string>(out var s) ? JsonNode.Parse(s) : p[k]?.DeepClone(); } catch { return p[k]?.DeepClone(); } }
        double? N(string k) => J(k) is JsonValue v && v.TryGetValue<double>(out var x) ? x : null;
        if (J("components") is not JsonArray comps) return null;
        var morph = (string?)J("morphology") ?? "mixed";
        var method = Array.IndexOf(GrowMethodIds, (string?)J("method") ?? "");
        var byCount = comps.OfType<JsonObject>().Any(c => (double?)c["chains"] > 0);
        var o = new JsonObject
        {
            ["page"] = "blend", ["mode"] = byCount ? 2 : 0, ["chains"] = N("chains") ?? 8, ["density"] = N("density") ?? 0.5,
            ["morphology"] = morph == "slabs" ? 1 : morph == "droplet" ? 2 : 0, ["method"] = Math.Max(0, method), ["temperature"] = N("temperature") ?? 450,
        };
        var rows = new JsonArray();
        foreach (var c in comps.OfType<JsonObject>())
        {
            var spec = c["spec"] as JsonObject;
            var units = spec?["units"] as JsonArray;
            var name = units is { Count: 1 } ? (string?)units[0]?["name"] ?? "chain" : string.Join("-", units?.Select(u => (string?)u?["name"]) ?? []);
            rows.Add(new JsonObject
            {
                ["id"] = "", ["name"] = name, ["spec"] = spec?.DeepClone(), ["weight"] = (double?)c["weight"] ?? 50, ["dp"] = (double?)spec?["dp"] ?? 20,
                ["count"] = (double?)c["chains"] is > 0 and var n ? n : 8,
            });
        }
        o["rows"] = rows;
        return o;
    }

    /// <summary>A provenance step's params: an object, or a list of [name, value] pairs.</summary>
    private static JsonObject? StepParams(JsonNode? step) => step?["params"] switch
    {
        JsonObject o => o,
        JsonArray a => new JsonObject(a.OfType<JsonArray>().Where(x => x.Count >= 2 && x[0] is JsonValue)
                                       .Select(x => new KeyValuePair<string, JsonNode?>((string?)x[0] ?? "", x[1]?.DeepClone())).DistinctBy(kv => kv.Key)),
        _ => null,
    };
}
