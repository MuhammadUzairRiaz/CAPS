using System.Collections.ObjectModel;
using System.Globalization;
using System.Text;
using System.Text.Json.Nodes;
using CapsStudio.Interop;

namespace CapsStudio.ViewModels;

/// <summary>A conformer of a torsion scan: rank, state, φ and energy above the lowest.</summary>
public sealed record TorsionRow(int Rank, string State, double Phi, double Energy, int Point)
{
    public string RankText => Rank.ToString(CultureInfo.InvariantCulture);
    public string PhiText => Phi.ToString("0.0", CultureInfo.InvariantCulture) + "°";
    public string EnergyText => Energy.ToString("0.00", CultureInfo.InvariantCulture);
}

/// <summary>Torsion scan (design/boards/TorsionScan): the energy along a dihedral of the open structure, rigid or with
/// the rest relaxed at each step (CAPS Field's force field, or UFF), the conformers it finds and the barrier.</summary>
public sealed partial class MainViewModel
{
    public bool IsTorsion => _module == 33;
    public ObservableCollection<TorsionRow> TorsionRows { get; } = new();
    public static readonly string[] TorsionModes = ["Rigid rotation", "Relax others · L-BFGS"];
    public static readonly string[] TorsionSides = ["d side (the last atom)", "a side (the first atom)"];
    public event Action? TorsionChanged;

    private int[] _torAtoms = [-1, -1, -1, -1];
    private string _torAtomsText = "";
    private decimal _torFrom = -180, _torTo = 180, _torStep = 15;
    private int _torMode = 1, _torSide;
    private bool _torRunning;
    private JsonNode? _tor;

    public void OpenTorsion()
    {
        if (_doc == null) return;
        // the four picked atoms, else a backbone torsion
        _torAtoms = _selection.Count == 4 ? _selection.ToArray() : _doc.DefaultTorsion() ?? [-1, -1, -1, -1];
        _torAtomsText = string.Join(" ", _torAtoms.Select(a => a >= 0 ? AtomName(a) : "?"));
        _tor = null;
        TorsionRows.Clear();
        Raise(nameof(TorsionAtomsText));
        SetModule(33);
        TorsionError = _torAtoms.Contains(-1) ? "Pick four bonded atoms in the Studio (⇧ click) or type their numbers, e.g. 1 2 3 4" : "";
        RaiseTorsion();
    }

    private string AtomName(int i)
    {
        try { var a = _doc!.Atom(i); return $"{a.ElementSymbol}{i + 1}"; }
        catch { return (i + 1).ToString(CultureInfo.InvariantCulture); }
    }

    /// <summary>The four atoms as typed: numbers (1-based) with or without element letters, "C1 C2 C3 C4" or "1 2 3 4".</summary>
    public string TorsionAtomsText
    {
        get => _torAtomsText;
        set
        {
            if (!Set(ref _torAtomsText, value ?? "")) return;
            var nums = _torAtomsText.Split([' ', ',', '-', '–'], StringSplitOptions.RemoveEmptyEntries)
                .Select(t => new string(t.Where(char.IsDigit).ToArray())).Where(t => t.Length > 0).Select(t => int.Parse(t, CultureInfo.InvariantCulture) - 1).ToArray();
            _torAtoms = nums.Length == 4 ? nums : [-1, -1, -1, -1];
            TorsionError = nums.Length == 4 ? "" : "Four atom numbers, e.g. 1 2 3 4";
            TorsionChanged?.Invoke();
        }
    }
    public int[] TorsionAtoms => _torSide == 0 ? _torAtoms : _torAtoms.Reverse().ToArray();
    public decimal TorsionFrom { get => _torFrom; set => Set(ref _torFrom, Math.Clamp(value, -180, 180)); }
    public decimal TorsionTo { get => _torTo; set => Set(ref _torTo, Math.Clamp(value, -180, 180)); }
    public decimal TorsionStep { get => _torStep; set => Set(ref _torStep, Math.Clamp(value, 1, 90)); }
    public int TorsionMode { get => _torMode; set => Set(ref _torMode, Math.Clamp(value, 0, 1)); }
    public int TorsionSide { get => _torSide; set { if (Set(ref _torSide, Math.Clamp(value, 0, 1))) TorsionChanged?.Invoke(); } }
    public bool TorsionIdle => !_torRunning;
    private string _torError = "", _torLog = "", _torFf = "";
    public string TorsionError { get => _torError; private set { if (Set(ref _torError, value)) Raise(nameof(TorsionHasError)); } }
    public bool TorsionHasError => _torError.Length > 0;
    public string TorsionLog { get => _torLog; private set => Set(ref _torLog, value); }
    public string TorsionBarrierText => _tor?["barrier"] is JsonNode b ? $"barrier {b.GetValue<double>():0.00} kcal/mol" : "barrier —";
    public string TorsionStepChip => $"scan points · {_torStep:0.#}° step";
    public string TorsionMinChip
    {
        get
        {
            var g = TorsionRows.FirstOrDefault(r => r.State.StartsWith("gauche"));
            return g != null ? $"gauche at ±{Math.Abs(g.Phi):0.0}°" : TorsionRows.FirstOrDefault() is { } t ? $"lowest at {t.Phi:0.0}°" : "";
        }
    }
    public string TorsionHud => $"{Title} · {_torFf}";
    public string TorsionPhiText { get; private set; } = "";
    public string TorsionStatus => $"{_doc?.Summary().Atoms ?? 0:N0} atoms · {_doc?.Summary().Bonds ?? 0:N0} bonds";
    public string TorsionCite => _torFf == "UFF" ? "Force field: UFF, Rappé et al., J. Am. Chem. Soc. 114, 10024 (1992)" : _torFf.Length > 0 ? "Force field: " + _torFf + " (CAPS Field)" : "Uses the Field assignment when complete, else UFF";

    private void RaiseTorsion()
    {
        foreach (var n in new[] { nameof(TorsionBarrierText), nameof(TorsionStepChip), nameof(TorsionMinChip), nameof(TorsionHud), nameof(TorsionPhiText), nameof(TorsionStatus), nameof(TorsionCite), nameof(TorsionIdle) })
            Raise(n);
        TorsionChanged?.Invoke();
    }

    /// <summary>The scan curve: (φ, E above the lowest) sorted by φ.</summary>
    public (double X, double Y)[] TorsionCurve() =>
        _tor?["points"] is JsonArray p ? p.Select(x => (x!["phi"]!.GetValue<double>(), x["energy"]!.GetValue<double>())).OrderBy(q => q.Item1).ToArray() : [];

    public async Task RunTorsionScan()
    {
        if (_doc == null || _torRunning) return;
        if (_torAtoms.Contains(-1)) { TorsionError = "Pick four bonded atoms first"; return; }
        _torRunning = true;
        Raise(nameof(TorsionIdle));
        TorsionError = "";
        var doc = _doc;
        var opts = new JsonObject
        {
            ["atoms"] = new JsonArray(TorsionAtoms.Select(a => (JsonNode)a).ToArray()), ["from"] = (double)_torFrom, ["to"] = (double)(_torTo <= _torFrom ? _torFrom + 360 : _torTo),
            ["step"] = (double)_torStep, ["relax"] = _torMode == 1, ["forcefield"] = "auto",
        }.ToJsonString();
        try
        {
            doc.TorsionShow(-1);   // from the structure as it was
        }
        catch { }
        try
        {
            var text = await Task.Run(() => doc.TorsionScan(opts, (d, t) =>
            {
                Avalonia.Threading.Dispatcher.UIThread.Post(() => Status = $"Torsion scan · point {d} of {t}");
                return true;
            }));
            var j = JsonNode.Parse(text)!;
            if (j["ok"]?.GetValue<bool>() != true) { TorsionError = j["error"]?.GetValue<string>() ?? "the scan failed"; return; }
            _tor = j;
            _torFf = j["forcefield"]!.GetValue<string>();
            TorsionRows.Clear();
            var points = j["points"]!.AsArray();
            var k = 0;
            foreach (var c in j["conformers"]!.AsArray())
            {
                var phi = c!["phi"]!.GetValue<double>();
                // the scan point nearest the conformer, to show it
                var nearest = Enumerable.Range(0, points.Count).OrderBy(i => Math.Abs(Math.IEEERemainder(points[i]!["phi"]!.GetValue<double>() - phi, 360))).First();
                TorsionRows.Add(new TorsionRow(++k, c["state"]!.GetValue<string>(), phi, c["energy"]!.GetValue<double>(), nearest));
            }
            TorsionPhiText = $"φ = {j["phi_start"]!.GetValue<double>():0.0}°";
            var notes = j["notes"]!.AsArray().Select(x => x!.GetValue<string>()).ToList();
            TorsionLog = $"{points.Count} points · {j["moving"]!.GetValue<double>():0} atoms rotated · {_torFf}" + (notes.Count > 0 ? "\n" + string.Join("\n", notes) : "");
            Status = $"Torsion scan done · {TorsionBarrierText}";
            if (TorsionRows.FirstOrDefault() is { } best) ShowTorsionRow(best);
        }
        catch (Exception e) { TorsionError = e.Message; }
        finally
        {
            _torRunning = false;
            RaiseTorsion();
        }
    }

    /// <summary>Shows a conformer's scan geometry in the structure.</summary>
    public void ShowTorsionRow(TorsionRow r)
    {
        if (_doc == null) return;
        try
        {
            _doc.TorsionShow(r.Point);
            TorsionPhiText = $"φ = {r.Phi:0.0}° · {r.State}";
            Raise(nameof(TorsionPhiText));
            RefreshSummary();
            RenderRequested?.Invoke();
            TorsionChanged?.Invoke();
        }
        catch (Exception e) { TorsionError = e.Message; }
    }

    /// <summary>φ, energies by term (kcal/mol) as CSV.</summary>
    public string TorsionCsv()
    {
        var sb = new StringBuilder("phi_deg,energy_kcal_mol,dihedral,vdw,coulomb,angle,bond,total\n");
        if (_tor?["points"] is JsonArray p)
            foreach (var x in p)
                sb.AppendLine(string.Join(",", new[] { "phi", "energy", "dihedral", "vdw", "coulomb", "angle", "bond", "total" }
                    .Select(k => x![k]!.GetValue<double>().ToString("0.#####", CultureInfo.InvariantCulture))));
        return sb.ToString();
    }
}
