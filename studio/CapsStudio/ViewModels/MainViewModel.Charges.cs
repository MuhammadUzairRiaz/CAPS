using System;
using System.Collections.ObjectModel;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Text.Json.Nodes;

namespace CapsStudio.ViewModels;

public sealed record ChargeGroupRow(string Group, string N, string Mean, string Range);

/// <summary>Studio › Charges (design/boards/Charges): partial charges computed without touching the structure — Gasteiger,
/// QEq, the Field assignment's, or imported (.chg from RESP or AM1-BCC) — their distribution, net charge and groups, the
/// view coloured by charge; applied to the topology or exported as mol2 on request.</summary>
public partial class MainViewModel
{
    public bool IsCharges => _module == 50;
    public static readonly string[] ChargeMethods = ["gasteiger", "qeq", "forcefield", "file"];
    private int _chgMethod;
    private string _chgFile = "", _chgNet = "—", _chgMax = "—", _chgStatus = "", _chgNote = "", _chgError = "";
    public int ChargeMethod { get => _chgMethod; set { if (Set(ref _chgMethod, Math.Clamp(value, 0, 3)) && (value != 3 || _chgFile.Length > 0)) ComputeCharges(); } }
    public string ChargeFile { get => _chgFile; private set { if (Set(ref _chgFile, value)) Raise(nameof(HasChargeFile)); } }
    public bool HasChargeFile => _chgFile.Length > 0;
    public string ChargeNet { get => _chgNet; private set => Set(ref _chgNet, value); }
    public string ChargeMax { get => _chgMax; private set => Set(ref _chgMax, value); }
    public string ChargeStatus { get => _chgStatus; private set => Set(ref _chgStatus, value); }
    public string ChargeNote { get => _chgNote; private set => Set(ref _chgNote, value); }
    public string ChargeError { get => _chgError; private set { Set(ref _chgError, value); Raise(nameof(ChargeHasError)); } }
    public bool ChargeHasError => _chgError.Length > 0;
    private int _chgFormal;
    public string ChargeTarget => _chgFormal == 0 ? "0 (neutral)" : _chgFormal.ToString("+0;−0", CultureInfo.InvariantCulture) + " (formal charges)";
    public ObservableCollection<ChargeGroupRow> ChargeGroups { get; } = new();
    /// <summary>The computed charges the view shows while the Charges page is open (none once applied or left).</summary>
    private double[]? _chgPreview;
    public bool ChargePreviewing => _chgPreview != null;
    private void SetChargePreview(double[]? q)
    {
        _chgPreview = q;
        Raise(nameof(ChargePreviewing));
        ApplyAppearance();
    }
    /// <summary>Leaving the Charges page: the view shows the structure's own charges again.</summary>
    private void EndChargePreview() { if (_chgPreview != null) SetChargePreview(null); }
    public (double X, double Y)[] ChargeHistogram { get; private set; } = [];
    public string ChargeRange { get; private set; } = "";
    public event Action? ChargesChanged;

    public void OpenCharges()
    {
        SetModule(50);
        AppColour = 4;   // the view coloured by partial charge
        ComputeCharges();
    }

    public void ChooseChargeFile(string path) { ChargeFile = path; _chgMethod = 3; Raise(nameof(ChargeMethod)); ComputeCharges(); }

    private JsonObject ChargeRequest(bool apply) => new()
    {
        ["method"] = ChargeMethods[_chgMethod], ["path"] = _chgFile, ["apply"] = apply,
    };

    public void ComputeCharges()
    {
        if (_doc == null) return;
        var inv = CultureInfo.InvariantCulture;
        var r = JsonNode.Parse(_doc.Charges(ChargeRequest(false).ToJsonString()))!;
        ChargeGroups.Clear();
        if (r["ok"]?.GetValue<bool>() != true)
        {
            ChargeError = r["error"]?.GetValue<string>() ?? "cannot compute charges";
            ChargeNet = ChargeMax = "—";
            ChargeHistogram = [];
            EndChargePreview();
            ChargesChanged?.Invoke();
            return;
        }
        ChargeError = "";
        _chgFormal = (int)(r["formal"]?.GetValue<double>() ?? 0);
        var net = r["net"]!.GetValue<double>();
        var max = r["max_abs"]!.GetValue<double>();
        ChargeNet = net.ToString("+0.000000;−0.000000;0.000000", inv) + " e";
        ChargeMax = max.ToString("0.0000", inv) + " e";
        ChargeRange = $"−{max.ToString("0.000", inv)}   0   +{max.ToString("0.000", inv)}";
        foreach (var g in (JsonArray)r["groups"]!)
            ChargeGroups.Add(new ChargeGroupRow((string)g!["name"]!, ((int)g["n"]!.GetValue<double>()).ToString(inv),
                g["mean"]!.GetValue<double>().ToString("+0.0000;−0.0000;0.0000", inv),
                $"{g["lo"]!.GetValue<double>().ToString("+0.000;−0.000;0.000", inv)} … {g["hi"]!.GetValue<double>().ToString("+0.000;−0.000;0.000", inv)}"));
        var e = ((JsonArray)r["edges"]!).Select(x => x!.GetValue<double>()).ToArray();
        var c = ((JsonArray)r["counts"]!).Select(x => x!.GetValue<double>()).ToArray();
        ChargeHistogram = c.Select((n, i) => ((e[i] + e[i + 1]) / 2, n)).ToArray();
        ChargeNote = string.Join(" · ", ((JsonArray)r["notes"]!).Select(x => (string?)x ?? ""));
        ChargeStatus = $"{_doc.Summary().Atoms:N0} atoms · net {net.ToString("0.0e0", inv)} e · {MethodName(_chgMethod)} · shown in the view, not applied";
        if (IsCharges && r["q"] is JsonArray qa) SetChargePreview(qa.Select(x => x!.GetValue<double>()).ToArray());
        Raise(nameof(ChargeTarget)); Raise(nameof(ChargeRange));
        ChargesChanged?.Invoke();
    }

    private static string MethodName(int m) => m switch { 0 => "Gasteiger–Marsili · 6 iterations", 1 => "QEq", 2 => "from the force field", _ => "imported" };

    /// <summary>Sets the charges on the structure (undoable; the Field assignment is cleared).</summary>
    public void ApplyCharges()
    {
        if (_doc == null) return;
        var r = JsonNode.Parse(_doc.Charges(ChargeRequest(true).ToJsonString()))!;
        if (r["ok"]?.GetValue<bool>() != true) { ChargeError = r["error"]?.GetValue<string>() ?? "cannot apply"; return; }
        _chgPreview = null;
        Raise(nameof(ChargePreviewing));
        AfterEdit($"Charges applied · {MethodName(_chgMethod)} (undo with ⌘Z)");
        ApplyAppearance();
        ChargeStatus = $"Applied · {MethodName(_chgMethod)}";
    }

    /// <summary>A mol2 with these charges (the structure itself is left as it was).</summary>
    public void ExportChargesMol2(string path)
    {
        if (_doc == null) return;
        var r = JsonNode.Parse(_doc.Charges(ChargeRequest(true).ToJsonString()))!;
        if (r["ok"]?.GetValue<bool>() != true) { ChargeError = r["error"]?.GetValue<string>() ?? "cannot export"; return; }
        try { _doc.Save(path); }
        finally { _doc.Undo(false); }
        Status = $"Wrote {Path.GetFileName(path)} with {MethodName(_chgMethod)} charges";
    }
}
