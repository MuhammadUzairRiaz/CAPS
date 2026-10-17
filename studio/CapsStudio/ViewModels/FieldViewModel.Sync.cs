using System.Collections.ObjectModel;
using System.Linq;
using System.Text.Json;
using System.Text.Json.Nodes;
using CapsStudio.Interop;

namespace CapsStudio.ViewModels;

/// <summary>An overlay file of the library (caps ff overlays): a parameter refit merged onto its force field before typing.</summary>
public sealed class FieldOverlayRow : ObservableObject
{
    public string File { get; init; } = "";
    public string Path { get; init; } = "";
    public string Name { get; init; } = "";
    public string Kinds { get; init; } = "";
    public string Tip { get; init; } = "";
    private bool _on;
    public bool On { get => _on; set { if (Set(ref _on, value)) Toggled?.Invoke(); } }
    internal Action? Toggled;
    internal void SetQuietly(bool on) { _on = on; Raise(nameof(On)); }
}

/// <summary>The command-line options of the Field step (caps ff apply / relax / md --overlay, caps ff complete) in the
/// Studio, through the same C API entry points: overlays on the assignment, what a force field cannot type (once per
/// motif, and why), and the missing bonded terms fitted with CAPS's machine-learned model when the build has it.</summary>
public sealed partial class FieldViewModel
{
    public ObservableCollection<FieldOverlayRow> Overlays { get; } = new();
    public bool HasOverlays => Overlays.Count > 0;
    private string _overlaysFor = "";
    private bool _loadingOverlays;

    protected override void OnChanged(string? name)
    {
        if (name == nameof(FfIndex)) LoadOverlays();
    }

    /// <summary>The library's overlays for the chosen force field (none for UFF or a force field without any).</summary>
    public void LoadOverlays()
    {
        var id = Selected?.Id ?? "";
        if (id == _overlaysFor) return;
        _overlaysFor = id;
        Overlays.Clear();
        if (id.Length > 0 && id != "uff" && Paths.ForceFields is { } dir)
        {
            try
            {
                var j = JsonNode.Parse(CapsDocument.ForcefieldOverlays(dir, id));
                foreach (var o in j?["overlays"]?.AsArray() ?? [])
                {
                    if (o == null) continue;
                    var k = o["kinds"];
                    var parts = new[] { "pairs", "bonds", "angles", "dihedrals", "impropers" }
                        .Select(x => ((int?)k?[x] ?? 0) is var n && n > 0 ? $"{n} {x}" : "").Where(x => x.Length > 0);
                    var notes = string.Join("\n", (o["notes"]?.AsArray() ?? []).Select(x => (string?)x ?? ""));
                    var row = new FieldOverlayRow
                    {
                        File = (string?)o["file"] ?? "", Path = (string?)o["path"] ?? "", Name = (string?)o["name"] ?? "",
                        Kinds = string.Join(", ", parts), Tip = notes.Length > 0 ? notes : (string?)o["name"] ?? "",
                    };
                    row.Toggled = OverlayToggled;
                    Overlays.Add(row);
                }
            }
            catch (Exception e) { Log = "Cannot list the overlays: " + e.Message; }
        }
        Raise(nameof(HasOverlays));
        Raise(nameof(ApplyCommand));
    }

    /// <summary>The overlay files switched on, in the library's order (each later file over the earlier).</summary>
    public IReadOnlyList<string> ChosenOverlays => Overlays.Where(o => o.On).Select(o => o.Path).ToList();

    private void OverlayToggled()
    {
        Raise(nameof(ApplyCommand));
        if (_loadingOverlays || !Assigned || IsGrouped || _assignedId != (Selected?.Id ?? "")) return;
        var chosen = ChosenOverlays;
        Recorder?.Invoke($"doc.field.set_overlays([{string.Join(", ", chosen.Select(p => "\"" + System.IO.Path.GetFileNameWithoutExtension(p) + "\""))}])");
        _ = Do(chosen.Count == 0 ? "Overlays removed" : $"Overlays: {string.Join(", ", chosen.Select(System.IO.Path.GetFileName))}", d => d.FieldSetOverlays(chosen));
    }

    /// <summary>After an assignment: the overlays switched on applied to it (the same merge as --overlay).</summary>
    private bool ApplyChosenOverlays(CapsDocument d, bool complete)
    {
        var chosen = ChosenOverlays;
        if (Selected?.Id == "uff") return complete;
        if (chosen.Count == 0 && d.FieldOverlays().Trim() is "[]" or "") return complete;
        return d.FieldSetOverlays(chosen);
    }

    /// <summary>The overlays of the report (what the assignment holds) shown as switched on.</summary>
    private void ReadOverlays(JsonElement r)
    {
        if (!r.TryGetProperty("overlays", out var ov) || ov.ValueKind != JsonValueKind.Array) return;
        var held = ov.EnumerateArray().Select(x => System.IO.Path.GetFileName(x.GetString() ?? "")).ToHashSet();
        _loadingOverlays = true;
        foreach (var row in Overlays) row.SetQuietly(held.Contains(row.File));
        _loadingOverlays = false;
        ExtraOverlays = ov.EnumerateArray().Select(x => x.GetString() ?? "").Where(p => Overlays.All(o => o.File != System.IO.Path.GetFileName(p))).ToList();
        Raise(nameof(ExtraOverlayText));
        Raise(nameof(HasExtraOverlays));
        Raise(nameof(ApplyCommand));
    }

    /// <summary>Overlays outside the library (a caps ff complete result) the assignment holds.</summary>
    public List<string> ExtraOverlays { get; private set; } = new();
    public bool HasExtraOverlays => ExtraOverlays.Count > 0;
    public string ExtraOverlayText => ExtraOverlays.Count == 0 ? "" : "Also merged: " + string.Join(", ", ExtraOverlays.Select(System.IO.Path.GetFileName));

    // ---- what the force field cannot type: once per chemical motif, and what that means for it (ffdef untyped_motifs)
    private string _motifs = "", _outOfDomain = "";
    public string UntypedMotifs { get => _motifs; private set { if (Set(ref _motifs, value)) Raise(nameof(HasMotifs)); } }
    public string OutOfDomain { get => _outOfDomain; private set => Set(ref _outOfDomain, value); }
    public bool HasMotifs => _motifs.Length > 0;

    private void ReadMotifs(JsonElement r)
    {
        UntypedMotifs = r.TryGetProperty("untyped_motifs", out var m) && m.ValueKind == JsonValueKind.Array
            ? string.Join("\n", m.EnumerateArray().Select(x => "· " + x.GetString())) : "";
        OutOfDomain = r.TryGetProperty("out_of_domain", out var d) && d.ValueKind == JsonValueKind.String ? d.GetString() ?? "" : "";
    }

    // ---- CAPS ML, two ways (a libtorch build only; caps_ml_available, the same check as caps_build_info's "ml_engine"):
    // "Fit missing terms" fits now and saves an overlay file (caps ff complete, caps_ff_complete); the checkbox "Fill
    // missing terms with CAPS ML at every assignment" fits whatever is missing each time the force field is assigned
    // (--complete ml, caps_field_set_options {"complete": "ml"}). Both are off with MlMissingNote in a build without it.
    public static bool MlAvailable { get; } = SafeMl();
    private static bool SafeMl() { try { return CapsDocument.MlAvailable; } catch { return false; } }
    public const string MlMissingNote = "This build of CAPS has no ML engine (it was built without CAPS_WITH_TORCH), so missing terms cannot be fitted here. Build CAPS with -DCAPS_WITH_TORCH=ON, or run caps ff complete from such a build and add its overlay with Import.";
    public static string MlNoteStatic => MlAvailable
        ? "Fit missing terms: CAPS finds every bond, angle and torsion the force field lacks for this structure, scans each on a fragment with its machine-learned potential and fits the term, now, and saves the fit as an overlay file that is merged (caps ff complete). To have this done at every assignment instead, tick \"Fill missing terms with CAPS ML at every assignment\" under Missing before run."
        : MlMissingNote;
    public bool CanComplete => MlAvailable && Assigned && HasMissing && !IsGrouped && Selected is { Id: not "uff" } && NotWorking;
    private string _completeText = "";
    public string CompleteText { get => _completeText; private set => Set(ref _completeText, value); }

    /// <summary>caps ff complete on the open structure; the overlay goes beside it (or in ~/CAPS/overlays) and is merged.</summary>
    public async Task CompleteWithModel(string? structurePath)
    {
        if (!CanComplete || Selected is not { } e) return;
        var doc = _doc();
        if (doc == null) return;
        var dir = structurePath is { Length: > 0 } sp && System.IO.Path.GetDirectoryName(sp) is { Length: > 0 } d0 ? d0
            : System.IO.Path.Combine(Environment.GetFolderPath(Environment.SpecialFolder.UserProfile), "CAPS", "overlays");
        Directory.CreateDirectory(dir);
        var stem = structurePath is { Length: > 0 } ? System.IO.Path.GetFileNameWithoutExtension(structurePath) : "structure";
        var overlay = System.IO.Path.Combine(dir, $"{stem}.{e.Id}-complete.json");
        var opts = new JsonObject { ["overlay"] = overlay, ["apply"] = true }.ToJsonString();
        Recorder?.Invoke($"doc.field.complete(\"{e.File.Replace("\\", "/")}\", \"{overlay.Replace("\\", "/")}\")");
        CompleteText = "Fitting the missing terms with CAPS's model…";
        await Do("Missing terms fitted", d =>
        {
            var rep = JsonNode.Parse(d.FfComplete(e.File, opts, line => { Avalonia.Threading.Dispatcher.UIThread.Post(() => CompleteText = line); return true; }));
            var terms = rep?["terms"]?.AsArray() ?? [];
            var fitted = terms.Count(t => (string?)t?["status"] == "fitted");
            var left = rep?["missing_after"]?.AsArray().Count ?? 0;
            Avalonia.Threading.Dispatcher.UIThread.Post(() => CompleteText = $"{fitted} terms fitted, {left} still missing · overlay {System.IO.Path.GetFileName(overlay)}");
            return rep?["complete"]?.GetValue<bool>() == true;
        });
    }

    /// <summary>--ff, --overlay and --charges of the assigned force field for caps relax / md / equilibrate / react; "" when
    /// none is assigned (the command then takes the data file's own); null when the assignment has no command-line form (by
    /// group, a water model, the topology file's own: save it whole and pass that file with --ff).</summary>
    public string? CliArgs()
    {
        if (!Assigned) return "";
        if (IsGrouped || _assignedId is "" or "file") return null;
        var e = Library.FirstOrDefault(x => x.Id == _assignedId);
        if (e == null) return null;
        var overlays = ChosenOverlays.Concat(ExtraOverlays).ToList();
        var charges = _chargeMode switch { 1 => "types", 2 => "gasteiger", 3 => "keep", 4 => "qeq", 5 => "increments", _ => "auto" };
        return $" --ff {Quote(e.Id == "uff" ? "uff" : e.File)}" + (overlays.Count > 0 ? " --overlay " + Quote(string.Join(",", overlays)) : "") + $" --charges {charges}" +
               (CompleteMl ? " --complete ml" : "");
    }

    /// <summary>The assignment as the command line runs it (caps ff apply; the same --ff, --overlay and --charges take
    /// relax, md, equilibrate and react).</summary>
    public string ApplyCommand
    {
        get
        {
            if (Selected is not { } e) return "";
            var ff = e.Id == "uff" ? "uff" : e.File;
            var overlays = ChosenOverlays.Concat(ExtraOverlays).ToList();
            var charges = _chargeMode switch { 1 => "types", 2 => "gasteiger", 3 => "keep", 4 => "qeq", 5 => "increments", _ => "auto" };
            return $"caps ff apply STRUCTURE --ff {Quote(ff)}" + (overlays.Count > 0 ? " --overlay " + Quote(string.Join(",", overlays)) : "") +
                   $" --charges {charges}" + (CompleteMl ? " --complete ml" : "") + " -o STRUCTURE.data";
        }
    }
    internal static string Quote(string s) => s.Length > 0 && s.All(c => char.IsLetterOrDigit(c) || "/._-,:=+@".Contains(c)) ? s : "'" + s.Replace("'", "'\\''") + "'";
}
