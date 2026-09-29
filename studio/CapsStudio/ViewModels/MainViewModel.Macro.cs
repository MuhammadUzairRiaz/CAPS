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

    /// <summary>A run recorded from its Copy-as-Python script: the call itself (the script's opening — import, open,
    /// force field — and its closing prints and saves are left out, as the macro has its own).</summary>
    private void RecordScript(string script)
    {
        if (!_recording) return;
        var lines = script.Replace("\r", "").Split('\n');
        var body = new List<string>();
        var inString = false;
        foreach (var l in lines)
        {
            var t = l.TrimStart();
            if (!inString && (t.StartsWith("#") || t.StartsWith("import caps") || t.StartsWith("doc = caps.open(") || t.StartsWith("doc.field.assign(") || t.Length == 0)) continue;
            if (!inString && (t.StartsWith("print(") || t.StartsWith("doc.save(") || t.StartsWith("doc.save_trajectory("))) continue;
            body.Add(l);
            if (l.Split("\"\"\"").Length % 2 == 0) inString = !inString;   // inside a triple-quoted block
        }
        foreach (var l in body) Record(l);
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
    internal static string NativeLibraryPath =>
        Path.Combine(AppContext.BaseDirectory, OperatingSystem.IsMacOS() ? "libcaps.dylib" : OperatingSystem.IsWindows() ? "caps.dll" : "libcaps.so");

    // Target: 0 new documents (the script makes its own), 1 the open structure (caps.current(); caps.hand_back(doc) returns
    // the result, which opens when the macro ends). Where: this machine or a host (the script and the structure go up with
    // scp, python3 runs there with the host's caps package, the output streams back).
    public static readonly string[] MacroTargets = ["New documents", "The open structure (caps.current())"];
    private int _macroTarget, _macroWhere;
    public int MacroTarget { get => _macroTarget; set => Set(ref _macroTarget, Math.Clamp(value, 0, 1)); }
    public int MacroWhere { get => Math.Min(_macroWhere, _settings.Hosts.Count); set => Set(ref _macroWhere, Math.Clamp(value, 0, _settings.Hosts.Count)); }

    /// <summary>Runs the script with python3 (CAPS_PYTHON overrides), the caps package and this Studio's library — here or
    /// on a host.</summary>
    public async Task RunMacro()
    {
        if (_macroRunning) return;
        if (_macroTarget == 1 && _doc == null) { Status = "Target is the open structure: open or build one first"; return; }
        SaveMacro();
        var script = Path.Combine(MacroFolder, _macroName);
        var host = MacroWhere > 0 ? _settings.Hosts[MacroWhere - 1] : null;
        MacroRunning = true;
        MacroOutput = $">>> run {_macroName}{(host != null ? " on " + host.Name : "")}\n";
        var runDir = Path.Combine(MacroFolder, ".run", Path.GetFileNameWithoutExtension(_macroName));
        var result = Path.Combine(runDir, "result.data");
        string remoteDir = "";
        try
        {
            Directory.CreateDirectory(runDir);
            if (File.Exists(result)) File.Delete(result);
            if (_macroTarget == 1) _doc!.Save(Path.Combine(runDir, "current.data"));
            var python = PythonExe;   // Settings › Python & scripting, else CAPS_PYTHON, else python3
            ProcessStartInfo psi;
            if (host == null)
            {
                psi = new ProcessStartInfo(python) { RedirectStandardOutput = true, RedirectStandardError = true, UseShellExecute = false, WorkingDirectory = MacroFolder };
                psi.ArgumentList.Add("-u");
                if (_stopOnError) { psi.ArgumentList.Add("-X"); psi.ArgumentList.Add("faulthandler"); }
                psi.ArgumentList.Add(script);
                if (Paths.Python is { } pkg) psi.Environment["PYTHONPATH"] = pkg + (Environment.GetEnvironmentVariable("PYTHONPATH") is { Length: > 0 } pp ? Path.PathSeparator + pp : "");
                psi.Environment["CAPS_LIB"] = NativeLibraryPath;
                if (_macroTarget == 1) psi.Environment["CAPS_DOC"] = Path.Combine(runDir, "current.data");
                psi.Environment["CAPS_OUT"] = result;
            }
            else
            {
                var id = $"macro-{DateTime.Now:yyyyMMdd-HHmmss}";
                var mk = await Tool("ssh", SshArgs(host, $"mkdir -p \"{host.WorkDir}/{id}\" && cd \"{host.WorkDir}/{id}\" && pwd"), 30000);
                if (mk.Code != 0) throw new InvalidOperationException("ssh: " + (mk.Err.Length > 0 ? mk.Err.Split('\n')[0] : $"exit {mk.Code}"));
                remoteDir = mk.Out.Split('\n').Last().Trim();
                var up = new List<string> { script };
                if (_macroTarget == 1) up.Add(Path.Combine(runDir, "current.data"));
                var sent = await Tool("scp", ScpArgs(host, up, $"{Target(host)}:{remoteDir}/"), 120000);
                if (sent.Code != 0) throw new InvalidOperationException("scp: " + (sent.Err.Length > 0 ? sent.Err.Split('\n')[0] : $"exit {sent.Code}"));
                MacroOutput += $"sent to {host.Name}:{remoteDir}\n";
                var env = (_macroTarget == 1 ? "CAPS_DOC=current.data " : "") + "CAPS_OUT=result.data";
                psi = new ProcessStartInfo("ssh") { RedirectStandardOutput = true, RedirectStandardError = true, UseShellExecute = false };
                foreach (var a2 in SshArgs(host, $"cd {Q(remoteDir)} && {env} python3 -u {Q(Path.GetFileName(script))}")) psi.ArgumentList.Add(a2);
            }
            using var proc = Process.Start(psi) ?? throw new InvalidOperationException("cannot start " + psi.FileName);
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
            if (host != null && code == 0)   // the result, when the script handed one back
            {
                var back = await Tool("scp", ScpArgs(host, [$"{Target(host)}:{remoteDir}/result.data"], result), 120000);
                if (back.Code != 0 && File.Exists(result)) File.Delete(result);
            }
            var opened = code == 0 && File.Exists(result);
            if (opened) Open(result);
            Avalonia.Threading.Dispatcher.UIThread.Post(() =>
            {
                sb.AppendLine(code == 0 ? (opened ? "done · the result is open" : "done") : $"exit code {code}");
                MacroOutput = sb.ToString();
            });
            Status = code == 0 ? $"{_macroName} finished{(opened ? " · its result is open" : "")}" : $"{_macroName} stopped with exit code {code}";
        }
        catch (Exception e) { MacroOutput += e.Message + "\n(choose a Python 3 interpreter in Settings › Python & scripting; hosts: Settings › Compute & remote)\n"; }
        finally { MacroRunning = false; _macroProc = null; }
    }

    public void StopMacro() { try { _macroProc?.Kill(true); } catch { } }
}
