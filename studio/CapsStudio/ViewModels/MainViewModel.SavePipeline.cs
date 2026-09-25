using System.Security.Cryptography;
using System.Text;
using CapsStudio.Interop;

namespace CapsStudio.ViewModels;

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
    public string PipelineRunCommand => $"caps run {PipelineFileName} --input 'runs/*/{System.IO.Path.GetFileName(_doc?.Path is { Length: > 0 } p ? p : "structure.data")}'";
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

    private void RefreshPipelineYaml()
    {
        var file = _pipeWithSource && _doc?.Path is { Length: > 0 } p ? System.IO.Path.GetFileName(p) : null;
        var topo = file != null ? TopologyFor(_doc!.Path) : "";
        try { PipelineYaml = CapsDocument.PipelineToYaml(PipelineJson(), _pipeName, file, topo.Length > 0 ? System.IO.Path.GetFileName(topo) : null); }
        catch (Exception e) { PipelineYaml = "# " + e.Message; }
        PipelineYamlHash = Convert.ToHexString(SHA256.HashData(Encoding.UTF8.GetBytes(_pipeYaml))).ToLowerInvariant();
        Raise(nameof(PipelineFileName)); Raise(nameof(PipelineRunCommand)); Raise(nameof(PipelineStepsText));
    }

    public void SavePipelineYaml(string path)
    {
        RefreshPipelineYaml();
        File.WriteAllText(path, _pipeYaml);
        Status = $"Saved {System.IO.Path.GetFileName(path)} · sha256 {_pipeYamlHash[..12]}";
    }
}
