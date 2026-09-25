using System.Collections.ObjectModel;
using System.Globalization;
using System.Text.Json.Nodes;
using CapsStudio.Interop;

namespace CapsStudio.ViewModels;

/// <summary>A distance restraint for Relax: k (r − r0)² between two atoms (indices from 0).</summary>
public sealed class RestraintRow : ObservableObject
{
    public int I { get; init; }
    public int J { get; init; }
    public string Atoms { get; init; } = "";
    public double Measured { get; init; }
    public string MeasuredText => Measured.ToString("0.000", CultureInfo.InvariantCulture) + " Å now";
    private decimal _r0 = 3, _k = 10;
    public Action? Changed { get; set; }
    public decimal R0D { get => _r0; set { if (Set(ref _r0, Math.Clamp(value, 0m, 100m))) Changed?.Invoke(); } }
    public decimal KD { get => _k; set { if (Set(ref _k, Math.Clamp(value, 0m, 10000m))) Changed?.Invoke(); } }
}

public partial class MainViewModel
{
    // ---------------------------------------------------------------- Relax › distance restraints
    public ObservableCollection<RestraintRow> RelaxRestraints { get; } = new();
    public bool HasRestraints => RelaxRestraints.Count > 0;
    /// <summary>Two atoms picked in the view (Shift-click): their distance can become a restraint.</summary>
    public bool CanAddRestraint => _doc != null && _selection.Count == 2;

    public void AddMeasuredRestraint()
    {
        if (_doc == null || _selection.Count != 2) return;
        var (i, j) = (_selection[0], _selection[1]);
        if (RelaxRestraints.Any(r => (r.I == i && r.J == j) || (r.I == j && r.J == i))) { Status = "That pair is restrained already"; return; }
        var d = _doc.Measure([i, j]);
        var row = new RestraintRow { I = i, J = j, Atoms = $"{_doc.Atom(i).Id} {_doc.Atom(i).ElementSymbol} – {_doc.Atom(j).Id} {_doc.Atom(j).ElementSymbol}", Measured = d, R0D = Math.Round((decimal)d, 2), KD = 10 };
        row.Changed = ApplyRestraints;
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
        var a = new JsonArray(RelaxRestraints.Select(r => (JsonNode)new JsonObject { ["i"] = r.I, ["j"] = r.J, ["r0"] = (double)r.R0D, ["k"] = (double)r.KD }).ToArray());
        try { _doc.SetRestraints(a.ToJsonString()); } catch (Exception e) { Status = "Restraints: " + e.Message; }
    }

    /// <summary>A new document keeps the restraints only when it has the same atoms (a relaxed copy); otherwise they go.</summary>
    private void RestraintsFollow(CapsDocument doc)
    {
        if (RelaxRestraints.Count == 0) return;
        var n = doc.Summary().Atoms;
        if (RelaxRestraints.Any(r => r.I >= n || r.J >= n)) { RelaxRestraints.Clear(); Raise(nameof(HasRestraints)); return; }
        ApplyRestraints();
    }
}
