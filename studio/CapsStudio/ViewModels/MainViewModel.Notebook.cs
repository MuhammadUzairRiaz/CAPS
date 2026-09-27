using System;
using System.IO;
using System.Text.Json.Nodes;

namespace CapsStudio.ViewModels;

/// <summary>Open in notebook (design/boards/JupyterNotebook): a Jupyter notebook beside the structure that opens it with
/// the caps package, its force field assigned, a first analysis, and the 3D view — then the system opens it.</summary>
public partial class MainViewModel
{
    public string? WriteNotebook()
    {
        if (_doc == null) { Status = "Open or build a structure first"; return null; }
        // the structure as a file: its own when saved, else a copy in Documents/CAPS
        var src = _doc.Path;
        var saved = src is { Length: > 0 } && File.Exists(src) && !GrownUnsaved;
        var stem = string.Concat(Title.Replace(" (unsaved)", "").Split(Path.GetInvalidFileNameChars())).Replace(' ', '_');
        if (stem.Length == 0) stem = "structure";
        var dir = saved ? Path.GetDirectoryName(src)! : Path.Combine(Environment.GetFolderPath(Environment.SpecialFolder.MyDocuments), "CAPS");
        Directory.CreateDirectory(dir);
        var file = saved ? src! : Path.Combine(dir, Path.GetFileNameWithoutExtension(stem) + ".data");
        if (!saved) _doc.Save(file);
        var nb = Path.Combine(dir, Path.GetFileNameWithoutExtension(file) + ".ipynb");
        for (var k = 2; File.Exists(nb); k++) nb = Path.Combine(dir, $"{Path.GetFileNameWithoutExtension(file)}_{k}.ipynb");
        var data = Path.GetDirectoryName(Paths.ForceFields ?? "") ?? "";   // …/data (its forcefields folder's parent)
        var py = Path.Combine(data, "python");
        string R(string p) => "r\"" + p.Replace("\"", "\\\"") + "\"";
        var cells = new JsonArray();
        int id = 0;
        void Md(string text) => cells.Add(new JsonObject { ["cell_type"] = "markdown", ["id"] = $"cell-{id++}", ["metadata"] = new JsonObject(), ["source"] = text });
        void Code(string text) => cells.Add(new JsonObject { ["cell_type"] = "code", ["id"] = $"cell-{id++}", ["execution_count"] = null, ["metadata"] = new JsonObject(), ["outputs"] = new JsonArray(), ["source"] = text });
        Md($"# {Title.Replace(" (unsaved)", "")} in CAPS\n\nWritten by CAPS Studio: the structure opened with the `caps` package, its force field, a first analysis and the 3D view. Every Studio module has a Python call (`help(caps)`).");
        Code($"import sys\nsys.path.insert(0, {R(py)})   # the caps package shipped with the Studio\nimport caps\n\ndoc = caps.open({R(file)})\ndoc.summary()");
        var fieldLine = Field.Assigned && Field.Selected is { } fe && fe.File.Length > 0
            ? $"report = doc.field.assign({R(fe.File)})\nprint(report[\"complete\"], report.get(\"missing\", [])[:5])"
            : "# a force field from the library: pcff-frc, compass-frc, oplsaa2024-moltemplate, gaff2, uff …\nreport = doc.field.assign(\"uff\")\nprint(report[\"complete\"])";
        Code(fieldLine);
        Code("props = doc.analyze(\"density,rdf\")\nprops");
        Code("doc.view(style=\"ball-and-stick\")");
        var nbook = new JsonObject
        {
            ["cells"] = cells, ["nbformat"] = 4, ["nbformat_minor"] = 5,
            ["metadata"] = new JsonObject
            {
                ["kernelspec"] = new JsonObject { ["display_name"] = "Python 3", ["language"] = "python", ["name"] = "python3" },
                ["language_info"] = new JsonObject { ["name"] = "python" },
            },
        };
        File.WriteAllText(nb, nbook.ToJsonString(new System.Text.Json.JsonSerializerOptions { WriteIndented = true }));
        Status = $"Wrote {Path.GetFileName(nb)} beside {Path.GetFileName(file)}";
        return nb;
    }

    /// <summary>Writes the notebook, then opens it with the system's handler for .ipynb (Jupyter, VS Code …).</summary>
    public void OpenInNotebook()
    {
        var nb = WriteNotebook();
        if (nb == null) return;
        try
        {
            System.Diagnostics.Process.Start(new System.Diagnostics.ProcessStartInfo(nb) { UseShellExecute = true });
        }
        catch (Exception e) { Status = $"Wrote {Path.GetFileName(nb)}; no application opens .ipynb files here ({e.Message})"; }
    }
}
