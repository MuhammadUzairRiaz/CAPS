using System.Collections.ObjectModel;
using System.Diagnostics;
using System.Globalization;
using System.Text;
using System.Text.RegularExpressions;

namespace CapsStudio.ViewModels;

/// <summary>A recorded Studio action: its number and the Python line it became.</summary>
public sealed record RecordedCommand(int Number, string Text);

/// <summary>A parameter of a macro: name, type, default (promoted from a literal in the script).</summary>
public sealed record MacroParameter(string Name, string Type, string Default);

/// <summary>Macro recorder (design/boards/MacroRecorder): every Studio action records as one line of Python using the
/// caps package; the script is edited, literals are promoted to parameters, and it runs with python3 on this machine
/// with its output streamed. Scripts live in ~/.caps/macros.</summary>
public sealed partial class MainViewModel
{
    public bool IsMacro => _module == 36;
    public ObservableCollection<RecordedCommand> RecordedCommands { get; } = new();
    public ObservableCollection<MacroParameter> MacroParameters { get; } = new();
    public ObservableCollection<string> MacroFiles { get; } = new();
    public static string? MacroFolderOverride { get; set; }
    public static string MacroFolder => MacroFolderOverride ?? Path.Combine(AppSettings.Folder, "macros");

    private bool _recording;
    public bool Recording { get => _recording; set { if (Set(ref _recording, value)) { Raise(nameof(RecordingText)); if (value) Status = "Recording: every action becomes a line of Python"; } } }
    public string RecordingText => _recording ? "Recording" : "Paused";

    /// <summary>Records one action (a line of Python acting on `doc`).</summary>
    public void Record(string line)
    {
        if (!_recording) return;
        RecordedCommands.Add(new RecordedCommand(RecordedCommands.Count + 1, line));
        var body = _macroText.TrimEnd();
        // lines go inside the function when the script has one, else at the end
        var marker = "    return doc";
        if (line.StartsWith("doc = ") && body.Contains("    doc = None\n")) _macroText = body.Replace("    doc = None\n", "    " + line + "\n") + "\n";
        else if (body.Contains(marker)) _macroText = body.Replace(marker, "    " + line + "\n" + marker) + "\n";
        else _macroText = body + "\n" + line + "\n";
        Raise(nameof(MacroText));
    }

    private static string PyStr(string s) => "\"" + s.Replace("\\", "\\\\").Replace("\"", "\\\"") + "\"";
    private static string PyJsonArgs(string json)
    {
        // {"op":"add_h","atoms":[1,2]} → op="add_h", atoms=[1, 2]
        try
        {
            var o = System.Text.Json.Nodes.JsonNode.Parse(json)!.AsObject();
            return string.Join(", ", o.Select(kv => kv.Key + "=" + PyValue(kv.Value)));
        }
        catch { return "**" + json; }
    }
    private static string PyValue(System.Text.Json.Nodes.JsonNode? v) => v switch
    {
        null => "None",
        System.Text.Json.Nodes.JsonArray a => "[" + string.Join(", ", a.Select(PyValue)) + "]",
        System.Text.Json.Nodes.JsonValue x when x.TryGetValue<string>(out var s) => PyStr(s),
        System.Text.Json.Nodes.JsonValue x when x.TryGetValue<bool>(out var b) => b ? "True" : "False",
        _ => v.ToJsonString(),
    };

    public void RecordOpen(string path, string? topology)
    {
        // absolute paths: the script runs from the macros folder
        path = Path.GetFullPath(path);
        Record(!string.IsNullOrEmpty(topology) ? $"doc = caps.open({PyStr(path)}, {PyStr(Path.GetFullPath(topology))})" : $"doc = caps.open({PyStr(path)})");
    }
    public void RecordEdit(string json) => Record($"doc.edit({PyJsonArgs(json)})");

    private string _macroText = "# Recorded in CAPS Studio · replayable with python3\nimport caps\n\n\ndef macro():\n    doc = None\n    return doc\n\n\nif __name__ == \"__main__\":\n    macro()\n";
    public string MacroText { get => _macroText; set => Set(ref _macroText, value ?? ""); }
    private string _macroName = "macro.py";
    public string MacroName { get => _macroName; set { if (value != null && Set(ref _macroName, value)) { LoadMacro(value); Raise(nameof(MacroFileIndex)); } } }
    public int MacroFileIndex { get => MacroFiles.IndexOf(_macroName); set { if (value >= 0 && value < MacroFiles.Count) MacroName = MacroFiles[value]; } }

    public void OpenMacro()
    {
        Directory.CreateDirectory(MacroFolder);
        MacroFiles.Clear();
        foreach (var f in Directory.EnumerateFiles(MacroFolder, "*.py").OrderBy(f => f)) MacroFiles.Add(Path.GetFileName(f));
        if (!MacroFiles.Contains(_macroName)) MacroFiles.Insert(0, _macroName);
        SetModule(36);
        Raise(nameof(MacroName)); Raise(nameof(MacroFileIndex));
    }

    private void LoadMacro(string name)
    {
        var p = Path.Combine(MacroFolder, name);
        if (File.Exists(p)) { _macroText = File.ReadAllText(p); Raise(nameof(MacroText)); RefreshParameters(); }
    }

    public void SaveMacro(string? name = null)
    {
        Directory.CreateDirectory(MacroFolder);
        name ??= _macroName;
        if (!name.EndsWith(".py")) name += ".py";
        File.WriteAllText(Path.Combine(MacroFolder, name), _macroText);
        if (!MacroFiles.Contains(name)) MacroFiles.Add(name);
        _macroName = name;
        Raise(nameof(MacroName)); Raise(nameof(MacroFileIndex));
        Status = $"Saved {Path.Combine(MacroFolder, name)}";
    }

    public void NewMacro()
    {
        var k = 1;
        while (MacroFiles.Contains($"macro_{k}.py")) k++;
        _macroName = $"macro_{k}.py";
        _macroText = "# Recorded in CAPS Studio · replayable with python3\nimport caps\n\n\ndef macro():\n    doc = None\n    return doc\n\n\nif __name__ == \"__main__\":\n    macro()\n";
        RecordedCommands.Clear();
        MacroFiles.Add(_macroName);
        Raise(nameof(MacroName)); Raise(nameof(MacroFileIndex));
        Raise(nameof(MacroText));
        RefreshParameters();
    }

    /// <summary>The literal selected in the editor becomes a parameter of the macro's function (its default).</summary>
    public string? PromoteToParameter(string literal, string? name = null)
    {
        literal = literal.Trim();
        if (literal.Length == 0) return "Select a number or a quoted string in the script first";
        var isString = literal.Length >= 2 && (literal[0] == '"' && literal[^1] == '"' || literal[0] == '\'' && literal[^1] == '\'');
        var isNumber = double.TryParse(literal, NumberStyles.Float, CultureInfo.InvariantCulture, out _);
        if (!isString && !isNumber) return "Only a literal (a number or a quoted string) can become a parameter";
        name ??= isString ? $"text_{MacroParameters.Count + 1}" : $"value_{MacroParameters.Count + 1}";
        var type = isString ? "str" : literal.Contains('.') || literal.Contains('e') ? "float" : "int";
        var body = _macroText;
        var idx = body.IndexOf(literal, StringComparison.Ordinal);
        if (idx < 0) return "That text is not in the script";
        body = body[..idx] + name + body[(idx + literal.Length)..];
        // add to the signature: def macro(…):
        var m = Regex.Match(body, @"def macro\(([^)]*)\):");
        if (m.Success)
        {
            var args = m.Groups[1].Value.Trim();
            var add = $"{name}: {type} = {literal}";
            body = body[..m.Index] + $"def macro({(args.Length > 0 ? args + ", " : "")}{add}):" + body[(m.Index + m.Length)..];
        }
        _macroText = body;
        Raise(nameof(MacroText));
        RefreshParameters();
        return null;
    }

    private void RefreshParameters()
    {
        MacroParameters.Clear();
        var m = Regex.Match(_macroText, @"def macro\(([^)]*)\):");
        if (!m.Success) return;
        foreach (var a in m.Groups[1].Value.Split(',', StringSplitOptions.RemoveEmptyEntries | StringSplitOptions.TrimEntries))
        {
            var pm = Regex.Match(a, @"^(\w+)\s*(?::\s*(\w+))?\s*(?:=\s*(.+))?$");
            if (pm.Success) MacroParameters.Add(new MacroParameter(pm.Groups[1].Value, pm.Groups[2].Value, pm.Groups[3].Value));
        }
    }

    private string _macroOutput = "";
    public string MacroOutput { get => _macroOutput; private set => Set(ref _macroOutput, value); }
    private bool _macroRunning, _stopOnError = true;
    public bool MacroRunning { get => _macroRunning; private set { if (Set(ref _macroRunning, value)) Raise(nameof(MacroIdle)); } }
    public bool MacroIdle => !_macroRunning;
    public bool MacroStopOnError { get => _stopOnError; set => Set(ref _stopOnError, value); }
    private Process? _macroProc;

    /// <summary>The native library this Studio loaded (for the script's caps package).</summary>
    private static string NativeLibraryPath =>
        Path.Combine(AppContext.BaseDirectory, OperatingSystem.IsMacOS() ? "libcaps.dylib" : OperatingSystem.IsWindows() ? "caps.dll" : "libcaps.so");

    /// <summary>Runs the script with python3 (CAPS_PYTHON overrides), the caps package and this Studio's library.</summary>
    public async Task RunMacro()
    {
        if (_macroRunning) return;
        SaveMacro();
        var script = Path.Combine(MacroFolder, _macroName);
        MacroRunning = true;
        MacroOutput = $">>> run {_macroName}\n";
        try
        {
            var python = Environment.GetEnvironmentVariable("CAPS_PYTHON") is { Length: > 0 } p ? p : "python3";
            var psi = new ProcessStartInfo(python) { RedirectStandardOutput = true, RedirectStandardError = true, UseShellExecute = false, WorkingDirectory = MacroFolder };
            psi.ArgumentList.Add("-u");
            if (_stopOnError) { psi.ArgumentList.Add("-X"); psi.ArgumentList.Add("faulthandler"); }
            psi.ArgumentList.Add(script);
            if (Paths.Python is { } pkg) psi.Environment["PYTHONPATH"] = pkg + (Environment.GetEnvironmentVariable("PYTHONPATH") is { Length: > 0 } pp ? Path.PathSeparator + pp : "");
            psi.Environment["CAPS_LIB"] = NativeLibraryPath;
            using var proc = Process.Start(psi) ?? throw new InvalidOperationException("cannot start " + python);
            _macroProc = proc;
            var sb = new StringBuilder(MacroOutput);
            void Pump(StreamReader r) => Task.Run(async () =>
            {
                string? line;
                while ((line = await r.ReadLineAsync()) != null)
                {
                    var l = line;
                    Avalonia.Threading.Dispatcher.UIThread.Post(() => { sb.AppendLine(l); MacroOutput = sb.ToString(); });
                }
            });
            Pump(proc.StandardOutput);
            Pump(proc.StandardError);
            await proc.WaitForExitAsync();
            await Task.Delay(100);
            var code = proc.ExitCode;
            Avalonia.Threading.Dispatcher.UIThread.Post(() => { sb.AppendLine(code == 0 ? "done" : $"exit code {code}"); MacroOutput = sb.ToString(); });
            Status = code == 0 ? $"{_macroName} finished" : $"{_macroName} stopped with exit code {code}";
        }
        catch (Exception e) { MacroOutput += e.Message + "\n(set CAPS_PYTHON to a Python 3 interpreter)\n"; }
        finally { MacroRunning = false; _macroProc = null; }
    }

    public void StopMacro() { try { _macroProc?.Kill(true); } catch { } }
}
