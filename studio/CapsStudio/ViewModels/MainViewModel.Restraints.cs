using System.Collections.ObjectModel;
using System.Globalization;
using System.Text.Json.Nodes;
using CapsStudio.Interop;

namespace CapsStudio.ViewModels;

/// <summary>A restraint for Relax: k (r − r0)² between two atoms, or k (φ − φ0)² on the dihedral of four (indices
/// from 0; K and L −1 for a distance).</summary>
public sealed class RestraintRow : ObservableObject
{
    public int I { get; init; }
    public int J { get; init; }
    public int K { get; init; } = -1;
    public int L { get; init; } = -1;
    public bool IsDihedral => L >= 0;
    public string Atoms { get; init; } = "";
    public double Measured { get; init; }
    public string MeasuredText => IsDihedral ? Measured.ToString("0.0", CultureInfo.InvariantCulture) + "° now" : Measured.ToString("0.000", CultureInfo.InvariantCulture) + " Å now";
    public string TargetLabel => IsDihedral ? "φ₀ °" : "r₀ Å";
    public string KTip => IsDihedral ? "kcal/mol/rad²" : "kcal/mol/Å²";
    private decimal _r0 = 3, _k = 10;
    public Action? Changed { get; set; }
    public decimal R0D { get => _r0; set { if (Set(ref _r0, IsDihedral ? Math.Clamp(value, -180m, 180m) : Math.Clamp(value, 0m, 100m))) Changed?.Invoke(); } }
    public decimal KD { get => _k; set { if (Set(ref _k, Math.Clamp(value, 0m, 10000m))) Changed?.Invoke(); } }
}

public partial class MainViewModel
{
    // ---------------------------------------------------------------- Relax › distance restraints
    public ObservableCollection<RestraintRow> RelaxRestraints { get; } = new();
    public bool HasRestraints => RelaxRestraints.Count > 0;
    /// <summary>Two atoms picked in the view (Shift-click) give a distance restraint, four a dihedral one.</summary>
    public bool CanAddRestraint => _doc != null && _selection.Count is 2 or 4;

    public void AddMeasuredRestraint()
    {
        if (_doc == null || _selection.Count is not (2 or 4)) return;
        if (_selection.Count == 4)
        {
            var q = _selection.ToArray();
            if (RelaxRestraints.Any(r => r.IsDihedral && ((r.I == q[0] && r.J == q[1] && r.K == q[2] && r.L == q[3]) || (r.I == q[3] && r.J == q[2] && r.K == q[1] && r.L == q[0]))))
            { Status = "That dihedral is restrained already"; return; }
            var phi = _doc.Measure(q);
            var drow = new RestraintRow { I = q[0], J = q[1], K = q[2], L = q[3], Atoms = RestraintAtomsLabel(_doc, q), Measured = phi, R0D = Math.Round((decimal)phi, 1), KD = 50 };
            drow.Changed = ApplyRestraints;
            _restraintAtomCount = _doc.Summary().Atoms;
            RelaxRestraints.Add(drow);
            Raise(nameof(HasRestraints));
            ApplyRestraints();
            Status = string.Format(CultureInfo.InvariantCulture, "Dihedral restraint {0}: {1:F1}° now; set φ₀ (180 trans, ±60 gauche) and k, then Relax", drow.Atoms, phi);
            return;
        }
        var (i, j) = (_selection[0], _selection[1]);
        if (RelaxRestraints.Any(r => (r.I == i && r.J == j) || (r.I == j && r.J == i))) { Status = "That pair is restrained already"; return; }
        var d = _doc.Measure([i, j]);
        var row = new RestraintRow { I = i, J = j, Atoms = RestraintAtomsLabel(_doc, [i, j]), Measured = d, R0D = Math.Round((decimal)d, 2), KD = 10 };
        row.Changed = ApplyRestraints;
        _restraintAtomCount = _doc.Summary().Atoms;
        RelaxRestraints.Add(row);
        Raise(nameof(HasRestraints));
        ApplyRestraints();
        Status = string.Format(CultureInfo.InvariantCulture, "Restraint {0}: {1:F2} Å now; set the target and k, then Relax", row.Atoms, d);
    }

    public void RemoveRestraint(RestraintRow r)
    {
        RelaxRestraints.Remove(r);
        Raise(nameof(HasRestraints));
        ApplyRestraints();
    }

    /// <summary>The list as the core reads it (caps_set_restraints), on the open document.</summary>
    internal void ApplyRestraints()
    {
        if (_doc == null) return;
        var a = new JsonArray(RelaxRestraints.Select(r => (JsonNode)(r.IsDihedral
            ? new JsonObject { ["i"] = r.I, ["j"] = r.J, ["k"] = r.K, ["l"] = r.L, ["phi0"] = (double)r.R0D, ["kphi"] = (double)r.KD }
            : new JsonObject { ["i"] = r.I, ["j"] = r.J, ["r0"] = (double)r.R0D, ["k"] = (double)r.KD })).ToArray());
        try { _doc.SetRestraints(a.ToJsonString()); } catch (Exception e) { Status = "Restraints: " + e.Message; }
    }

    /// <summary>Atom count of the structure the restraints were picked on.</summary>
    private long _restraintAtomCount = -1;

    /// <summary>"id element – id element …" for the restrained atoms: the row label, and the identity a new document must match.</summary>
    private static string RestraintAtomsLabel(CapsDocument doc, IReadOnlyList<int> atoms) =>
        string.Join(" – ", atoms.Select(a => { var x = doc.Atom(a); return $"{x.Id} {x.ElementSymbol}"; }));

    /// <summary>A new document keeps the restraints only when it has the same atoms (a relaxed copy): the same atom count
    /// and the same id and element at every restrained index. Otherwise they go.</summary>
    private void RestraintsFollow(CapsDocument doc)
    {
        if (RelaxRestraints.Count == 0) return;
        var n = doc.Summary().Atoms;
        bool same = n == _restraintAtomCount && RelaxRestraints.All(r =>
            r.I < n && r.J < n && r.K < n && r.L < n
            && r.Atoms == RestraintAtomsLabel(doc, r.IsDihedral ? new[] { r.I, r.J, r.K, r.L } : new[] { r.I, r.J }));
        if (!same) { RelaxRestraints.Clear(); _restraintAtomCount = -1; Raise(nameof(HasRestraints)); return; }
        ApplyRestraints();
    }
}
