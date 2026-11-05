using System.Globalization;
using System.Text.Json.Nodes;

namespace CapsStudio.ViewModels;

/// <summary>Charges › QEq options (core qeq.hpp): the shielding kernel (Ohno–Klopman or the ns Slater densities of the
/// original method), the Ewald sum for periodic fillers, hydrogen's charge-dependent hardness and charge bounds per
/// element ("O:-1.2:0, Si:0:2.4").</summary>
public partial class MainViewModel
{
    public static readonly string[] QeqKernels = ["Ohno–Klopman shielding", "Slater ns densities (original QEq)"];
    private int _qeqKernel;
    private bool _qeqEwald, _qeqHydrogen;
    private string _qeqBounds = "";
    private decimal _qeqCutoff = 10;
    public bool ChargeIsQeq => _chgMethod == 1;
    public int QeqKernel { get => _qeqKernel; set { if (Set(ref _qeqKernel, Math.Clamp(value, 0, 1))) RecomputeQeq(); } }
    public bool QeqEwald { get => _qeqEwald; set { if (Set(ref _qeqEwald, value)) RecomputeQeq(); } }
    public bool QeqHydrogen { get => _qeqHydrogen; set { if (Set(ref _qeqHydrogen, value)) RecomputeQeq(); } }
    public string QeqBounds { get => _qeqBounds; set { if (Set(ref _qeqBounds, (value ?? "").Trim())) RecomputeQeq(); } }
    public decimal QeqCutoff { get => _qeqCutoff; set { if (Set(ref _qeqCutoff, Math.Clamp(value, 3, 30))) RecomputeQeq(); } }
    private void RecomputeQeq() { if (_chgMethod == 1) ComputeCharges(); }

    /// <summary>"O:-1.2:0, Si:0:2.4" → {"O": [-1.2, 0], "Si": [0, 2.4]}; a malformed entry throws with the reason.</summary>
    private JsonObject QeqJson()
    {
        var inv = CultureInfo.InvariantCulture;
        var bounds = new JsonObject();
        foreach (var item in _qeqBounds.Split(',', StringSplitOptions.RemoveEmptyEntries | StringSplitOptions.TrimEntries))
        {
            var p = item.Split(':');
            if (p.Length != 3 || !double.TryParse(p[1], NumberStyles.Float, inv, out var lo) || !double.TryParse(p[2], NumberStyles.Float, inv, out var hi))
                throw new FormatException($"bounds: '{item}' is not ELEMENT:min:max");
            bounds[p[0].Trim()] = new JsonArray(lo, hi);
        }
        return new JsonObject
        {
            ["kernel"] = _qeqKernel == 1 ? "slater" : "ohno", ["ewald"] = _qeqEwald, ["hydrogen"] = _qeqHydrogen,
            ["cutoff"] = (double)_qeqCutoff, ["bounds"] = bounds,
        };
    }
}
