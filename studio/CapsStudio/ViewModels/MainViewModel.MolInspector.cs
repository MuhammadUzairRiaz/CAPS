using System;
using System.Collections.ObjectModel;
using System.Globalization;
using System.Linq;
using System.Text.Json.Nodes;

namespace CapsStudio.ViewModels;

public sealed record MolInfoRow(string Section, string Key, string Value);

/// <summary>The inspector's molecule (design/boards/MoleculeInspector): the molecule of the picked atom (else the first) —
/// identity, masses, double-bond equivalents, shape from this geometry, and what needs a calculation.</summary>
public partial class MainViewModel
{
    public ObservableCollection<MolInfoRow> MolInfoRows { get; } = new();
    private string _miTitle = "", _miSmiles = "", _miFormula = "";
    public string MolInspectorTitle { get => _miTitle; private set => Set(ref _miTitle, value); }
    public bool HasMolInfo => MolInfoRows.Count > 0;

    public void RefreshMolInspector()
    {
        MolInfoRows.Clear();
        if (_doc == null) { Raise(nameof(HasMolInfo)); return; }
        var inv = CultureInfo.InvariantCulture;
        JsonNode r;
        try { r = JsonNode.Parse(_doc.MoleculeInfo(Math.Max(0, Picked)))!; }
        catch { Raise(nameof(HasMolInfo)); return; }
        if (r["ok"]?.GetValue<bool>() != true) { Raise(nameof(HasMolInfo)); return; }
        double D(string k) => r[k]?.GetValue<double>() ?? double.NaN;
        _miFormula = (string?)r["formula"] ?? "";
        _miSmiles = (string?)r["smiles"] ?? "";
        MolInspectorTitle = $"Molecule {D("molecule"):0} · {_miFormula}";
        void Row(string s, string k, string v) => MolInfoRows.Add(new MolInfoRow(s, k, v));
        Row("Identity", "Formula", _miFormula);
        Row("Identity", "SMILES", _miSmiles.Length > 0 ? _miSmiles : "—");
        Row("Identity", "Atoms · bonds · rings", $"{D("atoms"):0} · {D("bonds"):0} · {D("rings"):0}");
        Row("Mass", "Molecular weight", D("mass").ToString("0.000", inv) + " g/mol");
        Row("Mass", "Monoisotopic", double.IsFinite(D("monoisotopic")) ? D("monoisotopic").ToString("0.0000", inv) + " u" : "— (an element without a tabulated isotope)");
        Row("Mass", "Double-bond equivalents", D("dbe").ToString("0.#", inv));
        var I = ((JsonArray?)r["inertia"])?.Select(x => x!.GetValue<double>()).ToArray() ?? [];
        if (I.Length == 3) Row("Shape", "I_A · I_B · I_C", $"{I[0]:0.0} · {I[1]:0.0} · {I[2]:0.0} amu·Å²");
        Row("Shape", "I_C − I_A − I_B", D("inertia_defect").ToString("0.000", inv) + " amu·Å²");
        Row("Shape", "Radius of gyration (mass)", D("rg").ToString("0.000", inv) + " Å");
        Row("Shape", "Rotatable bonds", D("rotatable").ToString("0", inv));
        var charged = r["has_charges"]?.GetValue<bool>() == true;
        Row("Charges", "Net charge", charged ? D("net_charge").ToString("+0.000;−0.000;0.000", inv) + " e" : "no charges on the atoms");
        Row("Charges", "Dipole moment", charged && double.IsFinite(D("dipole")) ? D("dipole").ToString("0.00", inv) + " D" + (Math.Abs(D("net_charge")) > 0.01 ? " (about the centre of mass: the molecule is charged)" : "") : "needs charges (Studio › Charges)");
        Raise(nameof(HasMolInfo));
    }

    public string MolInfoTable => string.Join("\n", MolInfoRows.Select(r => $"{r.Key}\t{r.Value}"));

    public void AddMoleculeToLibrary()
    {
        if (_miSmiles.Length == 0) { Status = "No SMILES for this molecule"; return; }
        AddMyFragment(_miFormula, _miSmiles);
        Status = $"Added {_miFormula} to My fragments";
    }
}
