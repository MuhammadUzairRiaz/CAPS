using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Text;
using CapsStudio.Interop;

namespace CapsStudio.ViewModels;

/// <summary>Python console (the analysis dock's Python tab): one interpreter that keeps running, with this Studio's caps
/// package and library; the open structure is `doc` (saved for it, as a macro gets it), and "Open doc in Studio"
/// brings `doc` back as a new structure. Lines run as typed; ↑/↓ recall earlier ones.</summary>
public sealed partial class MainViewModel
{
    private Process? _py;
    private readonly StringBuilder _pyOut = new();
    private string _pyInput = "", _pyDir = "";
    private readonly List<string> _pyHistory = new();
    private int _pyHistoryAt = -1;
    public string ConsoleText => _pyOut.ToString();
    public string ConsoleInput { get => _pyInput; set => Set(ref _pyInput, value ?? ""); }
    public bool ConsoleRunning => _py is { HasExited: false };
    public event Action? ConsoleChanged;

    private void ConsoleWrite(string text)
    {
        lock (_pyOut)
        {
            _pyOut.Append(text);
            if (_pyOut.Length > 200_000) _pyOut.Remove(0, _pyOut.Length - 150_000);   // a long session keeps its tail
        }
        Avalonia.Threading.Dispatcher.UIThread.Post(() => { Raise(nameof(ConsoleText)); ConsoleChanged?.Invoke(); });
    }

    /// <summary>Starts the interpreter (again): `import caps`, and `doc` = the open structure when there is one.</summary>
    public void StartConsole()
    {
        StopConsole();
        _pyDir = Path.Combine(Path.GetTempPath(), "caps-console-" + Environment.ProcessId);
        Directory.CreateDirectory(_pyDir);
        var cur = Path.Combine(_pyDir, "current.data");
        var hasDoc = false;
        if (_doc != null) { try { _doc.Save(cur); hasDoc = true; } catch { } }
        var psi = new ProcessStartInfo(PythonExe) { RedirectStandardInput = true, RedirectStandardOutput = true, RedirectStandardError = true, UseShellExecute = false, WorkingDirectory = _pyDir };
        PythonProcess.Utf8Io(psi);
        foreach (var a in new[] { "-u", "-i", "-q" }) psi.ArgumentList.Add(a);
        if (Paths.Python is { } pkg) psi.Environment["PYTHONPATH"] = pkg + (Environment.GetEnvironmentVariable("PYTHONPATH") is { Length: > 0 } pp ? Path.PathSeparator + pp : "");
        if (Native.LibraryPath is { } lib) psi.Environment["CAPS_LIB"] = lib;
        if (hasDoc) psi.Environment["CAPS_DOC"] = cur;
        psi.Environment["CAPS_OUT"] = Path.Combine(_pyDir, "result.data");
        try { _py = Process.Start(psi); }
        catch (Exception e) { ConsoleWrite($"cannot start {PythonExe}: {e.Message}\n(choose the interpreter in Settings › Python & scripting)\n"); return; }
        if (_py == null) return;
        void Pump(StreamReader r) => Task.Run(async () =>
        {
            var buf = new char[4096];
            int n;
            while ((n = await r.ReadAsync(buf, 0, buf.Length)) > 0) ConsoleWrite(new string(buf, 0, n));
        });
        Pump(_py.StandardOutput);
        Pump(_py.StandardError);
        _py.EnableRaisingEvents = true;
        _py.Exited += (_, _) => { ConsoleWrite("\n[the interpreter stopped]\n"); Avalonia.Threading.Dispatcher.UIThread.Post(() => Raise(nameof(ConsoleRunning))); };
        // quiet prompts (the dock echoes each line itself), then caps and the structure
        _py.StandardInput.WriteLine("import sys; sys.ps1 = ''; sys.ps2 = ''");
        _py.StandardInput.WriteLine("import caps");
        _py.StandardInput.WriteLine(hasDoc
            ? "doc = caps.current(); print('caps ABI %d · doc = %s (%d atoms)' % (caps.abi_version(), " + PyString(Title.Replace(" (unsaved)", "")) + ", doc.summary()['atoms']))"
            : "doc = None; print('caps ABI %d · no structure open: doc = None' % caps.abi_version())");
        _py.StandardInput.Flush();
        Raise(nameof(ConsoleRunning));
    }

    private static string PyString(string s) => "'" + s.Replace("\\", "\\\\").Replace("'", "\\'") + "'";

    public void StopConsole()
    {
        if (_py is { HasExited: false }) { try { _py.Kill(true); } catch { } }
        _py?.Dispose();
        _py = null;
        Raise(nameof(ConsoleRunning));
    }

    /// <summary>Runs what is typed (several lines run as a block; a block is closed with an empty line).</summary>
    public void ConsoleRun()
    {
        var text = _pyInput;
        if (text.Trim().Length == 0 && !ConsoleRunning) return;
        if (!ConsoleRunning) StartConsole();
        if (_py == null) return;
        ConsoleWrite(string.Join("\n", text.Split('\n').Select(l => ">>> " + l)) + "\n");
        if (text.Trim().Length > 0) { _pyHistory.Add(text); _pyHistoryAt = -1; }
        foreach (var l in text.Split('\n')) _py.StandardInput.WriteLine(l);
        if (text.Contains('\n')) _py.StandardInput.WriteLine();   // close a block
        _py.StandardInput.Flush();
        ConsoleInput = "";
    }

    /// <summary>↑ (−1) and ↓ (+1) through the lines run before.</summary>
    public void ConsoleHistory(int step)
    {
        if (_pyHistory.Count == 0) return;
        _pyHistoryAt = _pyHistoryAt < 0 ? (step < 0 ? _pyHistory.Count - 1 : -1) : Math.Clamp(_pyHistoryAt + step, 0, _pyHistory.Count);
        ConsoleInput = _pyHistoryAt >= 0 && _pyHistoryAt < _pyHistory.Count ? _pyHistory[_pyHistoryAt] : "";
    }

    public void ClearConsole() { lock (_pyOut) _pyOut.Clear(); Raise(nameof(ConsoleText)); }

    /// <summary>`doc` back into the Studio as a new structure (caps.hand_back writes it; the Studio opens it).</summary>
    public async Task ConsoleOpenDoc()
    {
        if (!ConsoleRunning || _py == null) { Status = "Start the console first"; return; }
        var result = Path.Combine(_pyDir, "result.data");
        try { File.Delete(result); } catch { }
        _py.StandardInput.WriteLine("print('handed back: ' + caps.hand_back(doc))");
        _py.StandardInput.Flush();
        for (var k = 0; k < 100 && !File.Exists(result); ++k) await Task.Delay(100);
        await Task.Delay(150);   // written whole
        if (!File.Exists(result)) { Status = "doc was not written (is doc a structure?)"; return; }
        var keep = Path.Combine(_pyDir, $"console-{DateTime.Now:HHmmss}.data");
        File.Copy(result, keep, true);
        Show(CapsDocument.Open(keep), "doc from the console");
        Status = "The console's doc is open as a new structure";
    }
}
