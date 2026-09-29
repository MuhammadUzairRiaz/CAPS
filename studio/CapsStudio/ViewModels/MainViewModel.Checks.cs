using System.Collections.ObjectModel;
using System.Text.Json.Nodes;
using Avalonia.Media;
using CapsStudio.Interop;

namespace CapsStudio.ViewModels;

/// <summary>One file check: what was found, what CAPS did, and one action.</summary>
public sealed class FileCheckRow
{
    public string Level { get; init; } = "note";
    public string Title { get; init; } = "";
    public string Detail { get; init; } = "";
    public string Action { get; init; } = "";
    public string ActionText => Action switch { "wrap" => "Wrap view", "relax" => "Relax", "field" => "Field", _ => "" };
    public bool HasAction => ActionText.Length > 0;
    public string Icon => Level switch { "pass" => "check", "error" => "xcircle", _ => "alert" };
    public IBrush Colour => Tokens.Brush(Level switch { "pass" => "OkB", "note" => "SelB", "warn" => "WarnB", _ => "ErrB" });
    public IBrush Border => Level is "warn" or "error" ? Colour : Tokens.Brush("LineB");
    public bool NeedsLook => Level is "warn" or "error";
}

/// <summary>One atom type of the open file, its element and mass editable (File checks › Atom types).</summary>
public sealed class AtomTypeRow : ObservableObject
{
    public int Type { get; init; }
    public string Label { get; init; } = "";
    public string Count { get; init; } = "";
    public string FromMass { get; init; } = "";
    private string _element = "";
    private decimal? _mass;
    public string Element { get => _element; set { if (Set(ref _element, value ?? "")) Raise(nameof(Mismatch)); } }
    public decimal? Mass { get => _mass; set => Set(ref _mass, value); }
    /// <summary>The element differs from the one the mass points to (a guess worth checking).</summary>
    public bool Mismatch => FromMass.Length > 0 && !string.Equals(FromMass, _element, StringComparison.OrdinalIgnoreCase);
    public string MismatchTip => Mismatch ? $"The mass points to {FromMass}" : "";
}

/// <summary>File checks (design/boards/VisProblems): every opened file is checked; each finding says what was found,
/// what was done and what can be changed. The Studio's Validation panel shows the same list.</summary>
public sealed partial class MainViewModel
{
    public bool IsChecks => _module == 17;
    public ObservableCollection<FileCheckRow> FileChecks { get; } = new();
    public IEnumerable<FileCheckRow> ChecksFine => FileChecks.Where(c => !c.NeedsLook);
    public IEnumerable<FileCheckRow> ChecksLook => FileChecks.Where(c => c.NeedsLook);
    public bool HasChecksLook => FileChecks.Any(c => c.NeedsLook);
    public bool NoChecksLook => !HasChecksLook;
    public string ChecksSummary
    {
        get
        {
            int pass = FileChecks.Count(c => c.Level == "pass"), note = FileChecks.Count(c => c.Level == "note"), look = FileChecks.Count(c => c.NeedsLook);
            return $"{pass} passed · {note} note{(note == 1 ? "" : "s")}" + (look > 0 ? $" · {look} need a look" : "");
        }
    }

    private void LoadFileChecks()
    {
        FileChecks.Clear();
        if (_doc != null)
        {
            try
            {
                foreach (var c in JsonNode.Parse(_doc.FileChecks())!.AsArray())
                    FileChecks.Add(new FileCheckRow
                    {
                        Level = (string?)c!["level"] ?? "note", Title = (string?)c["title"] ?? "", Detail = (string?)c["detail"] ?? "", Action = (string?)c["action"] ?? "",
                    });
            }
            catch (Exception e) { FileChecks.Add(new FileCheckRow { Level = "error", Title = "Checks could not run", Detail = e.Message }); }
        }
        foreach (var n in new[] { nameof(ChecksFine), nameof(ChecksLook), nameof(HasChecksLook), nameof(NoChecksLook), nameof(ChecksSummary) }) Raise(n);
    }

    public void OpenChecks() { LoadFileChecks(); LoadTypeRows(); SetModule(17); }

    // ---------------------------------------------------------------- atom types: element and mass by hand
    public ObservableCollection<AtomTypeRow> TypeRows { get; } = new();
    public bool HasTypeRows => TypeRows.Count > 0;
    private void LoadTypeRows()
    {
        TypeRows.Clear();
        try
        {
            if (_doc != null)
                foreach (var t in JsonNode.Parse(_doc.TypeTable())!.AsArray())
                {
                    if ((int)t!["type"]!.GetValue<double>() == 0) continue;   // untyped atoms (a file without types)
                    TypeRows.Add(new AtomTypeRow
                    {
                        Type = (int)t["type"]!.GetValue<double>(), Label = (string?)t["label"] ?? "", Count = ((int)t["count"]!.GetValue<double>()).ToString("N0"),
                        FromMass = (string?)t["from_mass"] ?? "", Element = (string?)t["element"] ?? "",
                        Mass = t["mass"]!.GetValue<double>() > 0 ? Math.Round((decimal)t["mass"]!.GetValue<double>(), 4) : null,
                    });
                }
        }
        catch { }
        Raise(nameof(HasTypeRows));
    }
    /// <summary>Every atom of the row's type made its element, with its mass (one undoable edit; the type is kept).</summary>
    public void ApplyTypeRow(AtomTypeRow r)
    {
        if (_doc == null) return;
        var ok = r.Mass is { } m && m > 0
            ? RunEdit(new { op = "type_element", type = r.Type, element = r.Element.Trim(), mass = (double)m })
            : RunEdit(new { op = "type_element", type = r.Type, element = r.Element.Trim() });
        if (ok == null) return;
        LoadFileChecks();
        LoadTypeRows();
    }

    public void RunCheckAction(FileCheckRow r)
    {
        switch (r.Action)
        {
            case "wrap": Wrap = true; SetModule(8); Status = "Wrapped into the box (view only; the coordinates are kept)"; break;
            case "relax": SetModule(2); break;
            case "field": SetModule(7); break;
        }
    }

    public void ExportChecks(string path)
    {
        var lines = new List<string> { $"# File checks · {Title}", "" };
        lines.AddRange(FileChecks.Select(c => $"- [{c.Level}] {c.Title} — {c.Detail}"));
        File.WriteAllLines(path, lines);
        Status = $"Check report written to {path}";
    }
}
