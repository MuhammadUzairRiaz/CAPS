using System;
using System.Globalization;
using System.Linq;
using System.Threading.Tasks;

namespace CapsStudio.ViewModels;

public sealed record InterfaceRow(string Quantity, string ThisFrame, string Trajectory);

/// <summary>Analyze › Interface (design/boards/Interface): mass density along the surface normal for the surface's
/// elements and the film, the gap and adsorbed layer, the film's plateau density and the work of adhesion — for the frame
/// on screen and, with a trajectory, averaged over the chosen frames.</summary>
public sealed partial class MainViewModel
{
    public bool IsInterfacePage => _module == 48;
    public System.Collections.ObjectModel.ObservableCollection<InterfaceRow> IfRows { get; } = new();
    private string _ifStatus = "", _ifNote = "";
    public string IfStatus { get => _ifStatus; private set => Set(ref _ifStatus, value); }
    public string IfNote { get => _ifNote; private set => Set(ref _ifNote, value); }
    public (string Label, (double X, double Y)[] Data)[] IfSeries { get; private set; } = [];
    public (double From, double To)? IfGapBand { get; private set; }
    public event Action? IfChanged;

    public void OpenInterface() { SetModule(48); IfChanged?.Invoke(); }

    public async Task RunInterface()
    {
        if (_doc == null || Analyze.Working) return;
        var inv = CultureInfo.InvariantCulture;
        var frames = (int)Math.Max(1, _doc.Summary().Frames);
        // this frame
        var (f0, l0) = (Analyze.FirstD, Analyze.LastD);
        Analyze.FirstD = _frame; Analyze.LastD = _frame;
        IfStatus = "Density profile of this frame…";
        await RunChips("zprofile", "adhesion");
        var one = Snapshot();
        Analyze.FirstD = f0; Analyze.LastD = l0;
        // the trajectory
        (double[] Z, string[] Keys, System.Collections.Generic.Dictionary<string, double> V)? traj = null;
        if (frames > 1)
        {
            IfStatus = $"Averaging over {frames} frames…";
            await RunChips("zprofile", "adhesion");
            traj = Snapshot();
        }
        string Val(System.Collections.Generic.Dictionary<string, double>? v, string key, string fmt, string unit) =>
            v != null && v.TryGetValue(key, out var x) && double.IsFinite(x) ? x.ToString(fmt, inv) + unit : "—";
        IfRows.Clear();
        foreach (var (q, key, fmt, unit) in new[]
        {
            ("Surface top", "surface top (Å)", "0.00", " Å"),
            ("Film plateau ρ", "plateau", "0.000", " g/cm³"),
            ("Film reaches ½ plateau", "film reaches half its plateau at z (Å)", "0.0", " Å"),
            ("Gap to the surface", "gap to the surface (Å)", "0.0", " Å"),
            ("Adsorbed layer thickness", "adsorbed layer thickness (Å)", "0.0", " Å"),
            ("First-layer peak", "first-layer peak (g/cm³)", "0.000", " g/cm³"),
            ("Work of adhesion", "adhesion", "0", " mJ/m²"),
        })
            IfRows.Add(new InterfaceRow(q, Val(one.V, key, fmt, unit), frames > 1 ? Val(traj?.V, key, fmt, unit) : "needs a trajectory"));
        IfNote = frames > 1
            ? $"This frame: frame {_frame + 1}. Trajectory: frames {Analyze.FirstD + 1}–{(Analyze.LastD < 0 ? frames : Analyze.LastD + 1)}."
            : "One frame: a single snapshot, not an equilibrated average. Equilibrate, then run on the trajectory to fill the right-hand column.";
        IfStatus = Analyze.Log;
        IfChanged?.Invoke();
    }

    /// <summary>The profile curves and the numbers of the last run (zprofile + adhesion).</summary>
    private (double[] Z, string[] Keys, System.Collections.Generic.Dictionary<string, double> V) Snapshot()
    {
        var v = new System.Collections.Generic.Dictionary<string, double>();
        var z = Analyze.Results.FirstOrDefault(r => r.Id == "zprofile");
        if (z != null)
        {
            v["plateau"] = z.Value;
            foreach (var (k, x) in z.Extra) v[k] = x;
            var curves = Analyze.Curves.Where(c => c.Property == z.Name).ToList();
            // the film, then the surface's elements by the mass they carry (Si and O before H on silica)
            IfSeries = curves.Where(c => c.Label != "all atoms" && c.Label != "molecule 1 (surface)")
                             .OrderByDescending(c => c.Label == "other molecules (film)").ThenByDescending(c => c.Y.Where(double.IsFinite).Sum())
                             .Select(c => (c.Label == "other molecules (film)" ? "film" : c.Label, c.X.Zip(c.Y).ToArray())).ToArray();
            IfGapBand = v.TryGetValue("surface top (Å)", out var top) && v.TryGetValue("film reaches half its plateau at z (Å)", out var half) ? (top, half) : null;
        }
        var a = Analyze.Results.FirstOrDefault(r => r.Id == "adhesion");
        if (a != null) v["adhesion"] = a.Value;
        return ([], [], v);
    }
}
