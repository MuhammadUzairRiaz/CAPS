using System.Globalization;
using System.Text.Json.Nodes;
using CapsStudio.Interop;

namespace CapsStudio.ViewModels;

/// <summary>χ from pair contacts (Fan, Olafson, Blanco &amp; Hsu, Macromolecules 1992; core chipair.hpp) on the Solvent screen
/// and the Blend phase diagram: pair energies of rigid molecules at van der Waals contact, Boltzmann-averaged, with
/// coordination numbers from packing, GAFF2 for every molecule. A predictive screen: the pages show it beside known
/// behaviour, failures included.</summary>
public sealed partial class MainViewModel
{
    private static string? Gaff2Path => Paths.ForceFields is { } d && File.Exists(Path.Combine(d, "gaff-amber25.json")) ? Path.Combine(d, "gaff-amber25.json") : null;

    private sealed record ContactFit(double A, double B, double Err, string Forcefield, string Notes);
    private readonly Dictionary<string, ContactFit> _ssContact = new();
    private bool _ssContactBusy;
    private string _ssContactProgress = "";
    private CancellationTokenSource? _ssContactCts;
    public bool SsContactBusy { get => _ssContactBusy; private set { if (Set(ref _ssContactBusy, value)) { Raise(nameof(SsCanContact)); Raise(nameof(SsContactIdle)); } } }
    public bool SsContactIdle => !_ssContactBusy;
    public string SsContactProgress { get => _ssContactProgress; private set => Set(ref _ssContactProgress, value); }
    /// <summary>A polymer from the file (its repeat unit is known) and GAFF2 in the library.</summary>
    public bool SsCanContact => !_ssContactBusy && SsPolymerJson()?["unit"] is JsonValue && Gaff2Path != null;
    public bool SsHasContact => SsRows.Any(r => r.ChiC != "—");

    private string SsKey(string solvent) => ((string?)SsPolymerJson()?["id"] ?? "?") + "|" + solvent;

    private (string Chi, string Pred, bool Agrees, bool Mismatch) SsContactFor(string solvent, double t, bool knownGood, bool knownBad, string known)
    {
        if (!_ssContact.TryGetValue(SsKey(solvent), out var f)) return ("—", "", false, false);
        var chi = f.A + f.B / t;
        var pred = chi < 0.45 ? "solvent" : chi <= 0.55 ? "borderline" : "non-solvent";   // the same rule as the Hildebrand column
        var mismatch = (pred == "solvent" && knownBad) || (pred == "non-solvent" && knownGood);
        var agrees = (pred == "solvent" && knownGood) || (pred == "non-solvent" && knownBad) || (pred == "borderline" && known.Contains('Θ'));
        return (chi.ToString("0.00", Inv), pred, agrees, mismatch);
    }

    private string SsContactSummary(bool known)
    {
        var rows = SsRows.Where(r => r.ChiC != "—").ToArray();
        if (rows.Length == 0) return "";
        var hild = SsRows.Where(r => r.ChiC != "—").Count(r => r.Agrees);
        var cont = rows.Count(r => r.AgreesC);
        var bad = rows.Where(r => r.MismatchC).Select(r => r.Name + (r.PredictedC == "solvent" ? " a solvent" : " a non-solvent")).ToArray();
        return known
            ? $"\n\nPair contacts (GAFF2, no measured input) agree with the known behaviour for {cont} of {rows.Length} solvents here (Hildebrand: {hild} of {rows.Length})" +
              (bad.Length > 0 ? $"; they call {string.Join(", ", bad)}." : ".") +
              " Averaging isolated pairs at contact favours the strongest single contacts (aromatic stacking, dipoles), which a liquid's crowded shell cannot all take: ethers and aromatics are where it fails most."
            : $"\n\nPair contacts (GAFF2): computed for {rows.Length} solvents; no known behaviour on file to check them against.";
    }

    /// <summary>χ(T) by pair contacts for the selected polymer's repeat unit against every solvent of the file.</summary>
    public async Task SsComputeContacts()
    {
        if (!SsCanContact || _solventData == null) return;
        var unit = (string)SsPolymerJson()!["unit"]!;
        var ff = Gaff2Path!;
        var solvents = (_solventData["solvents"] as JsonArray ?? []).OfType<JsonObject>().Where(s => s["smiles"] is JsonValue).ToArray();
        SsContactBusy = true;
        _ssContactCts = new CancellationTokenSource();
        var ct = _ssContactCts.Token;
        var failed = new List<string>();
        try
        {
            for (var i = 0; i < solvents.Length; ++i)
            {
                var name = (string?)solvents[i]["name"] ?? "";
                if (_ssContact.ContainsKey(SsKey(name))) continue;
                var smiles = (string)solvents[i]["smiles"]!;
                var k = i;
                SsContactProgress = $"{name} ({k + 1} of {solvents.Length})";
                var json = new JsonObject { ["a"] = unit, ["b"] = smiles, ["forcefield"] = ff, ["t"] = (double)_ssT }.ToJsonString();
                var res = await Task.Run(() => CapsDocument.ChiContacts(json, (st, f) =>
                {
                    Avalonia.Threading.Dispatcher.UIThread.Post(() => SsContactProgress = $"{name} ({k + 1} of {solvents.Length}) · {st}");
                    return !ct.IsCancellationRequested;
                }));
                if (ct.IsCancellationRequested) break;
                var r = JsonNode.Parse(res)!;
                if ((bool?)r["ok"] != true) { failed.Add($"{name}: {(string?)r["error"]}"); continue; }
                var notes = string.Join("\n", (r["notes"] as JsonArray ?? []).Select(x => (string?)x));
                _ssContact[SsKey(name)] = new ContactFit((double)r["fit_a"]!, (double)r["fit_b"]!, (double)r["chi_error"]!, (string?)r["forcefield"] ?? "", notes);
                SsRecompute();
                Raise(nameof(SsHasContact));
            }
            SsContactProgress = ct.IsCancellationRequested ? "stopped" : failed.Count == 0 ? "" : "not computed — " + string.Join("; ", failed);
            Status = ct.IsCancellationRequested ? "χ by pair contacts stopped" : $"χ by pair contacts for {solvents.Length - failed.Count} solvents";
        }
        catch (Exception e) { SsContactProgress = "failed: " + e.Message; }
        finally { SsContactBusy = false; SsRecompute(); Raise(nameof(SsHasContact)); }
    }
    public void SsStopContacts() => _ssContactCts?.Cancel();

    // ---------------------------------------------------------------- Blend phase: χ(T) = A + B/T fitted from pair contacts
    public string[] BpUnitNames { get { LoadPolymerLibrary(); return PolymerLibrary.Where(p => !p.Copolymer).Select(p => p.Name).ToArray(); } }
    private int _bpUnitA, _bpUnitB = 1, _bpSource;
    private bool _bpFitBusy;
    private string _bpFitNote = "";
    public int BpUnitA { get => _bpUnitA; set => Set(ref _bpUnitA, value); }
    public int BpUnitB { get => _bpUnitB; set => Set(ref _bpUnitB, value); }
    /// <summary>0: A and B entered; 1: fitted from pair contacts.</summary>
    public int BpSource { get => _bpSource; set { if (Set(ref _bpSource, value)) { Raise(nameof(BpInputsChip)); Raise(nameof(BpFootnote)); } } }
    public bool BpFitBusy { get => _bpFitBusy; private set { if (Set(ref _bpFitBusy, value)) Raise(nameof(BpCanFit)); } }
    public bool BpCanFit => !_bpFitBusy && Gaff2Path != null;
    public string BpFitNote { get => _bpFitNote; private set => Set(ref _bpFitNote, value); }

    public async Task BpFitContacts()
    {
        var lib = PolymerLibrary.Where(p => !p.Copolymer).ToList();
        if (!BpCanFit || _bpUnitA < 0 || _bpUnitB < 0 || _bpUnitA >= lib.Count || _bpUnitB >= lib.Count) return;
        var (pa, pb) = (lib[_bpUnitA], lib[_bpUnitB]);
        BpFitBusy = true;
        BpFitNote = $"{pa.Name} against {pb.Name}: pair contacts…";
        try
        {
            var json = new JsonObject { ["a"] = pa.Smiles, ["b"] = pb.Smiles, ["forcefield"] = Gaff2Path, ["t"] = (double)_bpT }.ToJsonString();
            var res = await Task.Run(() => CapsDocument.ChiContacts(json, (st, f) => { Avalonia.Threading.Dispatcher.UIThread.Post(() => BpFitNote = $"{pa.Name} against {pb.Name}: {st}"); return true; }));
            var r = JsonNode.Parse(res)!;
            if ((bool?)r["ok"] != true) { BpFitNote = "Could not fit: " + (string?)r["error"]; return; }
            double a = (double)r["fit_a"]!, b = (double)r["fit_b"]!, e = (double)r["chi_error"]!;
            _bpA = (decimal)Math.Round(a, 5);
            _bpB = (decimal)Math.Round(b, 3);
            Raise(nameof(BpA)); Raise(nameof(BpB));
            BpSource = 1;
            BpRecompute();
            var inv = CultureInfo.InvariantCulture;
            BpFitNote = string.Format(inv, "{0} against {1}, per repeat unit (capped with H): χ(T) = {2:0.0000} + {3:0.0}/T, ± {4:0.000} at {5:0} K · {6}. " +
                "A screen from pair contacts: in CAPS's checks it gets polyisoprene/polybutadiene (small χ) and polystyrene/polyisoprene (large) right but calls " +
                "polystyrene/PVME immiscible, which is miscible — treat the diagram as a first guess.", pa.Name, pb.Name, a, b, e, (double)_bpT, (string?)r["forcefield"]);
            Status = string.Format(inv, "χ(T) fitted from pair contacts: {0:0.0000} + {1:0.0}/T", a, b);
        }
        catch (Exception ex) { BpFitNote = "Could not fit: " + ex.Message; }
        finally { BpFitBusy = false; }
    }
}
