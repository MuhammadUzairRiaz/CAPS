using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using System.Text.Json.Nodes;

namespace CapsStudio.ViewModels;

// Minimise page › Cell with held parameters (core cell_relax.hpp through the caps_edit op relax_cell): the cell and the
// atoms relaxed with chosen lattice parameters kept; undoable like any edit.
public partial class MainViewModel
{
    private bool _holdA, _holdB, _holdC, _holdAlpha, _holdBeta, _holdGamma;
    public bool HoldA { get => _holdA; set => Set(ref _holdA, value); }
    public bool HoldB { get => _holdB; set => Set(ref _holdB, value); }
    public bool HoldC { get => _holdC; set => Set(ref _holdC, value); }
    public bool HoldAlpha { get => _holdAlpha; set => Set(ref _holdAlpha, value); }
    public bool HoldBeta { get => _holdBeta; set => Set(ref _holdBeta, value); }
    public bool HoldGamma { get => _holdGamma; set => Set(ref _holdGamma, value); }
    private decimal _cellPressure;
    public decimal CellPressure { get => _cellPressure; set => Set(ref _cellPressure, value); }
    private string _cellRelaxText = "";
    public string CellRelaxText { get => _cellRelaxText; private set { if (Set(ref _cellRelaxText, value)) Raise(nameof(CellRelaxHasText)); } }
    public bool CellRelaxHasText => _cellRelaxText.Length > 0;

    public void RunCellRelax()
    {
        if (_doc == null) return;
        var hold = new List<string>();
        if (_holdA) hold.Add("a");
        if (_holdB) hold.Add("b");
        if (_holdC) hold.Add("c");
        if (_holdAlpha) hold.Add("alpha");
        if (_holdBeta) hold.Add("beta");
        if (_holdGamma) hold.Add("gamma");
        Status = "Relaxing the cell…";
        var r = RunEdit(new { op = "relax_cell", hold = hold.ToArray(), pressure = (double)_cellPressure, ftol = 0.01 });
        if (r?["cell_relax"] is not JsonObject c) { CellRelaxText = EditError; return; }
        var inv = CultureInfo.InvariantCulture;
        double[] A(string k) => ((JsonArray)c[k]!).Select(x => x!.GetValue<double>()).ToArray();
        var b = A("before");
        var a = A("after");
        var st = A("stress");
        string[] names = ["a", "b", "c", "α", "β", "γ"];
        var parts = Enumerable.Range(0, 6).Select(k => $"{names[k]} {b[k].ToString(k < 3 ? "0.000" : "0.00", inv)} → {a[k].ToString(k < 3 ? "0.000" : "0.00", inv)}{(hold.Count > 0 && hold.Contains(new[] { "a", "b", "c", "alpha", "beta", "gamma" }[k]) ? " (held)" : "")}");
        var notes = ((JsonArray)c["notes"]!).Select(x => (string?)x ?? "").Where(x => x.Length > 0);
        CellRelaxText = string.Join(" · ", parts)
                        + $"\nV {c["volume_before"]!.GetValue<double>().ToString("0.00", inv)} → {c["volume_after"]!.GetValue<double>().ToString("0.00", inv)} Å³ · E {c["energy_before"]!.GetValue<double>().ToString("0.000", inv)} → {c["energy_after"]!.GetValue<double>().ToString("0.000", inv)} kcal/mol · {(int)c["steps"]!.GetValue<double>()} steps"
                        + $"\nstress GPa (xx yy zz yz xz xy) {string.Join(" ", st.Select(x => x.ToString("0.000", inv)))}"
                        + (notes.Any() ? "\n" + string.Join(" · ", notes) : "");
    }
}
