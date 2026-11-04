using System.Globalization;
using System.Text;
using System.Text.Json.Nodes;

namespace CapsStudio.ViewModels;

/// <summary>Relax page › restraints in every run (core restraints.hpp): distances, angles, torsions and positions held near a value or
/// between bounds in Relax, Dynamics and Equilibrate, and written into the LAMMPS and GROMACS exports.</summary>
public sealed partial class MainViewModel
{
    private string _rsText = "", _rsReport = "", _rsError = "";
    public string RunRestraintText { get => _rsText; set => Set(ref _rsText, value); }
    public string RunRestraintReport { get => _rsReport; private set { if (Set(ref _rsReport, value)) Raise(nameof(HasRunRestraints)); } }
    public bool HasRunRestraints => _rsReport.Length > 0;
    public string RunRestraintError { get => _rsError; private set { if (Set(ref _rsError, value)) Raise(nameof(RunRestraintHasError)); } }
    public bool RunRestraintHasError => _rsError.Length > 0;

    /// <summary>Sets the restraints from the text box.</summary>
    public void ApplyRunRestraints()
    {
        if (_doc == null) return;
        var (n, err) = _doc.SetRestraintText(_rsText);
        RunRestraintError = n < 0 ? err : "";
        if (n >= 0) Status = n == 0 ? "No restraints" : $"{n} restraint(s) act in Relax, Dynamics and Equilibrate and go into the LAMMPS and GROMACS files";
        SyncRunRestraints(keepText: n < 0);
    }

    /// <summary>A restraint from the selection at its present value: 2 atoms a distance, 3 an angle, 4 a torsion, else positions.</summary>
    public void RestrainSelectionForRuns()
    {
        if (_doc == null) return;
        var sel = new List<int>(_selection);
        try
        {
            if (JsonNode.Parse(_doc.SelectionJson())?["indices"] is JsonArray a)
                foreach (var x in a) { var i = (int)((double?)x ?? -1); if (i >= 0 && !sel.Contains(i)) sel.Add(i); }
        }
        catch { }
        if (sel.Count == 0) { Status = "Select atoms first: two for a distance, three for an angle, four for a torsion, any other number to tether them"; return; }
        var inv = CultureInfo.InvariantCulture;
        string line;
        var atoms = string.Join(" ", sel.Select(i => (i + 1).ToString(inv)));
        if (sel.Count >= 2 && sel.Count <= 4)
        {
            var kind = sel.Count == 2 ? "distance" : sel.Count == 3 ? "angle" : "torsion";
            // the present value, from the report of a trial restraint
            var probe = $"{kind} {atoms} 0";
            var (n, err) = _doc.SetRestraintText((_rsText.Length > 0 ? _rsText.TrimEnd() + "\n" : "") + probe);
            double v = 0;
            if (n > 0 && JsonNode.Parse(_doc.RestraintReport())?["restraints"] is JsonArray rr && rr.Count > 0) v = (double?)rr[^1]!["value"] ?? 0;
            _doc.SetRestraintText(_rsText);
            if (n < 0) { RunRestraintError = err; return; }
            line = $"{kind} {atoms} {v.ToString(sel.Count == 2 ? "0.000" : "0.0", inv)} {(sel.Count == 2 ? 10 : 50)}";
        }
        else line = $"position {RangesOf(sel)} 10";
        RunRestraintText = (_rsText.Length > 0 ? _rsText.TrimEnd() + "\n" : "") + line;
        ApplyRunRestraints();
    }

    private static string RangesOf(List<int> sel)
    {
        var a = sel.Distinct().OrderBy(i => i).ToList();
        var parts = new List<string>();
        for (var k = 0; k < a.Count;)
        {
            var j = k;
            while (j + 1 < a.Count && a[j + 1] == a[j] + 1) ++j;
            parts.Add(j == k ? $"{a[k] + 1}" : $"{a[k] + 1}-{a[j] + 1}");
            k = j + 1;
        }
        return string.Join(",", parts);
    }

    private void SyncRunRestraints(bool keepText = false)
    {
        if (_doc == null) { RunRestraintReport = ""; return; }
        try
        {
            var r = JsonNode.Parse(_doc.RestraintReport())!;
            if (!keepText) RunRestraintText = ((string?)r["text"] ?? "").TrimEnd();
            var list = r["restraints"] as JsonArray ?? [];
            if (list.Count == 0) { RunRestraintReport = ""; return; }
            var inv = CultureInfo.InvariantCulture;
            var sb = new StringBuilder();
            var shown = 0;
            foreach (var x in list)
            {
                if (++shown > 12) break;
                var kind = (string?)x!["kind"] ?? "";
                var at = string.Join("–", ((JsonArray)x["atoms"]!).Select(q => ((int)((double?)q ?? 0) + 1).ToString(inv)));
                var lo = (double?)x["lower"] ?? 0; var hi = (double?)x["upper"] ?? 0;
                var unit = kind == "angle" || kind == "torsion" ? "°" : " Å";
                var target = kind == "position" ? (hi > 0 ? $"within {hi.ToString("0.##", inv)} Å" : "at its anchor") : lo == hi ? lo.ToString("0.###", inv) + unit : $"{lo.ToString("0.###", inv)}…{hi.ToString("0.###", inv)}{unit}";
                var val = x["value"] is JsonNode vn ? ((double?)vn ?? 0).ToString("0.###", inv) + unit : "–";
                sb.AppendLine($"{kind} {at}: {target}, now {val}, {((double?)x["energy"] ?? 0).ToString("0.###", inv)} kcal/mol");
            }
            if (list.Count > 12) sb.AppendLine($"… {list.Count - 12} more");
            sb.Append($"total {((double?)r["energy"] ?? 0).ToString("0.000", inv)} kcal/mol");
            RunRestraintReport = sb.ToString();
        }
        catch (Exception e) { RunRestraintError = e.Message; }
    }
}
