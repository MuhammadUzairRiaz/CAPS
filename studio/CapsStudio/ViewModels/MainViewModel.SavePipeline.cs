using System.Security.Cryptography;
using System.Text;
using CapsStudio.Interop;

namespace CapsStudio.ViewModels;

/// <summary>One file a saved pipeline writes (the YAML's outputs: block): a table as CSV, a plot as SVG, the global
/// attributes, a render of the view or the grid.</summary>
public sealed class PipelineOutputRow : System.ComponentModel.INotifyPropertyChanged
{
    public event System.ComponentModel.PropertyChangedEventHandler? PropertyChanged;
    public string Kind { get; init; } = "table";
    public Action? Changed;
    private string _what = "", _path = "";
    private int _w = 1920, _h = 1080;
    public string What { get => _what; set { if (_what == value) return; _what = value; PropertyChanged?.Invoke(this, new(nameof(What))); Changed?.Invoke(); } }
    public string Path { get => _path; set { if (_path == value) return; _path = value; PropertyChanged?.Invoke(this, new(nameof(Path))); Changed?.Invoke(); } }
    public decimal Width { get => _w; set { _w = (int)Math.Clamp(value, 64, 8192); Changed?.Invoke(); } }
    public decimal Height { get => _h; set { _h = (int)Math.Clamp(value, 64, 8192); Changed?.Invoke(); } }
    public bool HasWhat => Kind is "table" or "plot";
    public bool IsRender => Kind == "render";
    public string KindText => Kind switch { "table" => "Table · CSV", "plot" => "Plot · SVG", "attributes" => "Attributes · CSV", "render" => "Render", "grid" => "Grid", _ => Kind };
    public System.Text.Json.Nodes.JsonObject Json()
    {
        var o = new System.Text.Json.Nodes.JsonObject { ["kind"] = Kind, ["path"] = _path };
        if (HasWhat) o["what"] = _what;
        if (IsRender) { o["width"] = _w; o["height"] = _h; }
        return o;
    }
    public static PipelineOutputRow From(System.Text.Json.Nodes.JsonObject o) => new()
    {
        Kind = (string?)o["kind"] ?? "table", _what = (string?)o["what"] ?? "", _path = (string?)o["path"] ?? "",
        _w = (int?)o["width"] ?? 1920, _h = (int?)o["height"] ?? 1080,
    };
}

/// <summary>Save pipeline (design/boards/SavePipeline): the steps as YAML to diff, review and version, its sha256 (the
/// hash results record), and the command that runs it elsewhere.</summary>
public sealed partial class MainViewModel
{
    public bool IsSavePipeline => _module == 28;
    private string _pipeName = "structure report", _pipeYaml = "", _pipeYamlHash = "";
    private bool _pipeWithSource = true;

    public void OpenSavePipeline()
    {
        SetModule(28);
        RefreshPipelineYaml();
    }

    public string PipelineName { get => _pipeName; set { if (Set(ref _pipeName, value)) RefreshPipelineYaml(); } }
    public bool PipelineWithSource { get => _pipeWithSource; set { if (Set(ref _pipeWithSource, value)) RefreshPipelineYaml(); } }
    public string PipelineYaml { get => _pipeYaml; private set => Set(ref _pipeYaml, value); }
    public string PipelineYamlHash { get => _pipeYamlHash; private set { if (Set(ref _pipeYamlHash, value)) Raise(nameof(PipelineYamlShort)); } }
    public string PipelineYamlShort => _pipeYamlHash.Length > 12 ? _pipeYamlHash[..12] + "…" : _pipeYamlHash;
    public string PipelineFileName => Slug(_pipeName) + ".caps-pipeline.yaml";
    public string PipelineRunCommand => $"caps run {PipelineFileName} --input 'runs/*/{System.IO.Path.GetFileName(_doc?.Path is { Length: > 0 } p ? p : "structure.data")}'" +
                                        (PipelineOutputs.Count > 0 ? " --out outputs" : "");
    public string PipelineStepsText => PipelineRows.Count == 0 ? "No steps yet: add them in Visualize." :
        string.Join("\n", PipelineRows.Reverse().Select((r, k) => $"{k + 1}. {r.Title}{(r.Enabled ? "" : " (off)")} · {r.Summary}"));

    private static string Slug(string s)
    {
        var b = new StringBuilder();
        foreach (var c in s.ToLowerInvariant()) b.Append(char.IsLetterOrDigit(c) ? c : '_');
        var t = b.ToString().Trim('_');
        while (t.Contains("__")) t = t.Replace("__", "_");
        return t.Length == 0 ? "pipeline" : t;
    }

    // scope (design/boards/SavePipeline): this project's pipelines/ folder, or the shared library in ~/CAPS/pipelines
    public static readonly string[] PipelineScopes = ["This project", "Shared library"];
    public static string PipelineLibrary => Environment.GetEnvironmentVariable("CAPS_PIPELINES") is { Length: > 0 } d ? d
        : System.IO.Path.Combine(Environment.GetFolderPath(Environment.SpecialFolder.UserProfile), "CAPS", "pipelines");
    private int _pipeScope = 1;
    public int PipelineScope { get => _projFile == null ? 1 : _pipeScope; set { if (Set(ref _pipeScope, Math.Clamp(value, 0, 1))) Raise(nameof(PipelineScopeFolder)); } }
    public string PipelineScopeFolder => PipelineScope == 0 && _projFile != null ? System.IO.Path.Combine(CapsProjectFile.FolderOf(_projFile), "pipelines") : PipelineLibrary;
    public string PipelineScopeNote => _projFile == null ? "no project open: the shared library" : RecentFiles.Tilde(PipelineScopeFolder);
    public bool PipelineStoreWithResults { get => _settings.PipelineStoreWithResults; set { if (_settings.PipelineStoreWithResults == value) return; _settings.PipelineStoreWithResults = value; Raise(); Changed("Store the pipeline with its results"); } }
    public bool PipelineAskPython { get => _settings.PipelineAskPython; set { if (_settings.PipelineAskPython == value) return; _settings.PipelineAskPython = value; Raise(); Changed("Ask before running a Python file"); RefreshPipelineYaml(); } }

    /// <summary>Save into the scope's folder under the pipeline's file name (no picker).</summary>
    public string? SavePipelineToScope()
    {
        try
        {
            var dir = PipelineScopeFolder;
            Directory.CreateDirectory(dir);
            var path = System.IO.Path.Combine(dir, PipelineFileName);
            SavePipelineYaml(path);
            return path;
        }
        catch (Exception e) { Status = "Could not save the pipeline: " + e.Message; return null; }
    }

    /// <summary>A pipeline file in your own library or the open project's pipelines/ folder (it may keep its Python steps on).</summary>
    private bool IsOwnPipelineFile(string path)
    {
        var full = System.IO.Path.GetFullPath(path);
        bool Under(string dir) => full.StartsWith(System.IO.Path.GetFullPath(dir) + System.IO.Path.DirectorySeparatorChar, OperatingSystem.IsLinux() ? StringComparison.Ordinal : StringComparison.OrdinalIgnoreCase);
        return Under(PipelineLibrary) || (_projFile != null && Under(System.IO.Path.Combine(CapsProjectFile.FolderOf(_projFile), "pipelines")));
    }

    private void RefreshPipelineYaml()
    {
        var file = _pipeWithSource && _doc?.Path is { Length: > 0 } p ? System.IO.Path.GetFileName(p) : null;
        var topo = file != null ? TopologyFor(_doc!.Path) : "";
        try
        {
            PipelineYaml = CapsDocument.PipelineToYaml(PipelineJson(), _pipeName, file, topo.Length > 0 ? System.IO.Path.GetFileName(topo) : null);
            // Ask before running a Python file off: said in the file (honoured only for your own library and project)
            if (!_settings.PipelineAskPython && PipelineRows.Any(r => r.Type == "python"))
            {
                var lines = _pipeYaml.Split('\n').ToList();
                var at = lines.FindIndex(l => l.StartsWith("name:", StringComparison.Ordinal));
                lines.Insert(at >= 0 ? at + 1 : Math.Min(1, lines.Count), "ask_before_python: false");
                PipelineYaml = string.Join("\n", lines);
            }
        }
        catch (Exception e) { PipelineYaml = "# " + e.Message; }
        PipelineYamlHash = Convert.ToHexString(SHA256.HashData(Encoding.UTF8.GetBytes(_pipeYaml))).ToLowerInvariant();
        Raise(nameof(PipelineFileName)); Raise(nameof(PipelineRunCommand)); Raise(nameof(PipelineStepsText)); Raise(nameof(PipelineTableNames));
        Raise(nameof(PipelineScope)); Raise(nameof(PipelineScopeFolder)); Raise(nameof(PipelineScopeNote));
    }

    // ---------------------------------------------------------------- outputs (the YAML's outputs: block)
    public System.Collections.ObjectModel.ObservableCollection<PipelineOutputRow> PipelineOutputs { get; } = new();
    public bool HasPipelineOutputs => PipelineOutputs.Count > 0;
    /// <summary>The tables the current result has (for an output's "what").</summary>
    public string PipelineTableNames => _pipeResult?["tables"] is System.Text.Json.Nodes.JsonArray ts && ts.Count > 0
        ? "tables: " + string.Join(", ", ts.Select(t => (string?)t?["name"] ?? "")) : "no tables in the result yet";

    public void AddPipelineOutput(string kind)
    {
        var first = _pipeResult?["tables"] is System.Text.Json.Nodes.JsonArray ts && ts.Count > 0 ? (string?)ts[^1]?["name"] ?? "" : "";
        var row = new PipelineOutputRow { Kind = kind };
        row.What = first;
        row.Path = kind switch
        {
            "table" => (first.Length > 0 ? first : "table") + ".csv",
            "plot" => (first.Length > 0 ? first : "plot") + ".svg",
            "attributes" => "attributes.csv",
            "render" => "view.png",
            _ => "density.cube",
        };
        WireOutput(row);
        PipelineOutputs.Add(row);
        Raise(nameof(HasPipelineOutputs));
        RefreshPipelineYaml();
    }
    public void RemovePipelineOutput(PipelineOutputRow row)
    {
        PipelineOutputs.Remove(row);
        Raise(nameof(HasPipelineOutputs));
        RefreshPipelineYaml();
    }
    private void WireOutput(PipelineOutputRow row) => row.Changed = () => { if (IsSavePipeline) RefreshPipelineYaml(); };
    private System.Text.Json.Nodes.JsonArray PipelineOutputsJson() => new(PipelineOutputs.Select(o => (System.Text.Json.Nodes.JsonNode)o.Json()).ToArray());

    /// <summary>Writes the outputs now, from the shown frame's result, under dir.</summary>
    public void WritePipelineOutputs(string dir)
    {
        if (_doc == null || PipelineOutputs.Count == 0) return;
        try
        {
            var r = System.Text.Json.Nodes.JsonNode.Parse(_doc.PipelineWriteOutputs(PipelineJson(), dir));
            var lines = (r?["lines"] as System.Text.Json.Nodes.JsonArray ?? []).Select(x => (string?)x ?? "").ToList();
            var wrote = lines.Count(l => l.StartsWith("wrote"));
            if (_settings.PipelineStoreWithResults && (string?)r?["error"] == null)
            {
                // the pipeline beside its results, so they say how they were made
                RefreshPipelineYaml();
                File.WriteAllText(System.IO.Path.Combine(dir, PipelineFileName), _pipeYaml);
                lines.Add("wrote " + PipelineFileName);
                ++wrote;
            }
            Status = (string?)r?["error"] ?? $"Wrote {wrote} of {lines.Count} outputs to {dir}" + (wrote < lines.Count ? " · " + lines.First(l => !l.StartsWith("wrote")) : "");
        }
        catch (Exception e) { Status = "Could not write the outputs: " + e.Message; }
    }

    public void SavePipelineYaml(string path)
    {
        RefreshPipelineYaml();
        File.WriteAllText(path, _pipeYaml);
        Status = $"Saved {System.IO.Path.GetFileName(path)} · sha256 {_pipeYamlHash[..12]}";
    }
}
