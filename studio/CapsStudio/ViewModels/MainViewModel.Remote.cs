using System.Collections.ObjectModel;
using System.Diagnostics;
using System.Globalization;

namespace CapsStudio.ViewModels;

/// <summary>A row of the host list: this machine or a remote host, with its last known state.</summary>
public sealed class HostRow : ObservableObject
{
    public HostRow(RemoteHost? host, string name, string detail) { Host = host; _name = name; _detail = detail; }
    public RemoteHost? Host { get; }
    private string _name, _detail, _state = "not tested";
    private int _level;   // 0 unknown, 1 ok, 2 error
    public string Name { get => _name; set => Set(ref _name, value); }
    public string Detail { get => _detail; set => Set(ref _detail, value); }
    public string State { get => _state; set => Set(ref _state, value); }
    public int Level { get => _level; set { if (Set(ref _level, value)) { Raise(nameof(IsOk)); Raise(nameof(IsError)); } } }
    public bool IsOk => _level == 1;
    public bool IsError => _level == 2;
    public bool IsLocal => Host == null;
}

/// <summary>Compute &amp; remote (design/boards/RemoteCompute): this machine and SSH hosts (credentials stay in the SSH
/// agent), a connection test that runs caps on the host, the batch-job template and what moves with a job.</summary>
public sealed partial class MainViewModel
{
    public static readonly string[] Schedulers = ["SLURM", "PBS", "none"];
    public ObservableCollection<HostRow> Hosts { get; } = new();
    private HostRow? _host;

    private void LoadHosts()
    {
        Hosts.Clear();
        Hosts.Add(new HostRow(null, "This machine", $"local · {Environment.ProcessorCount} threads · CPU") { State = "ready", Level = 1 });
        foreach (var h in _settings.Hosts) Hosts.Add(new HostRow(h, h.Name, HostDetail(h)));
        SelectedHost = Hosts.Count > 1 ? Hosts[1] : Hosts[0];
        Raise(nameof(HostsStatus));
    }

    private static string HostDetail(RemoteHost h) =>
        $"{(h.User.Length > 0 ? h.User + "@" : "")}{(h.Hostname.Length > 0 ? h.Hostname : "no hostname")}:{h.Port} · SSH · {(h.Scheduler == "none" ? "no scheduler" : h.Scheduler)}";

    public HostRow? SelectedHost
    {
        get => _host;
        set
        {
            if (!Set(ref _host, value)) return;
            foreach (var n in new[] { nameof(HostEditable), nameof(HostName), nameof(HostHostname), nameof(HostUser), nameof(HostPort), nameof(HostScheduler),
                                      nameof(HostPartition), nameof(HostWorkDir), nameof(HostTitle), nameof(HostTestText) }) Raise(n);
        }
    }
    public bool HostEditable => _host?.Host != null;
    public string HostTitle => _host?.Name ?? "";
    private RemoteHost H => _host?.Host ?? new RemoteHost();
    private void HostChanged()
    {
        if (_host?.Host is not { } h) return;
        _host.Name = h.Name;
        _host.Detail = HostDetail(h);
        _host.State = "not tested";
        _host.Level = 0;
        foreach (var n in new[] { nameof(HostTitle), nameof(HostName), nameof(HostHostname), nameof(HostUser), nameof(HostPort), nameof(HostScheduler), nameof(HostPartition), nameof(HostWorkDir) })
            Raise(n);
        Changed("Host");
    }
    public string HostName { get => H.Name; set { if (_host?.Host is { } h && value != null) { h.Name = value; HostChanged(); } } }
    public string HostHostname { get => H.Hostname; set { if (_host?.Host is { } h && value != null) { h.Hostname = value.Trim(); HostChanged(); } } }
    public string HostUser { get => H.User; set { if (_host?.Host is { } h && value != null) { h.User = value.Trim(); HostChanged(); } } }
    public decimal HostPort { get => H.Port; set { if (_host?.Host is { } h) { h.Port = (int)Math.Clamp(value, 1, 65535); HostChanged(); } } }
    public string HostScheduler { get => H.Scheduler; set { if (_host?.Host is { } h && value != null) { h.Scheduler = value; HostChanged(); } } }
    public string HostPartition { get => H.Partition; set { if (_host?.Host is { } h && value != null) { h.Partition = value.Trim(); HostChanged(); } } }
    public string HostWorkDir { get => H.WorkDir; set { if (_host?.Host is { } h && value != null) { h.WorkDir = value.Trim(); HostChanged(); } } }
    public string HostsStatus => $"{Hosts.Count(h => h.IsOk)} of {Hosts.Count} hosts reachable";

    public void AddHost()
    {
        var h = new RemoteHost { Name = $"host-{_settings.Hosts.Count + 1}" };
        _settings.Hosts.Add(h);
        var row = new HostRow(h, h.Name, HostDetail(h));
        Hosts.Add(row);
        SelectedHost = row;
        Changed("Host added");
    }

    public void RemoveHost()
    {
        if (_host?.Host is not { } h) return;
        _settings.Hosts.Remove(h);
        Hosts.Remove(_host);
        SelectedHost = Hosts.LastOrDefault();
        Changed("Host removed");
    }

    private string _hostTest = "";
    public string HostTestText { get => _hostTest; private set => Set(ref _hostTest, value); }
    private bool _hostTesting;
    public bool HostIdle { get => !_hostTesting; private set => Set(ref _hostTesting, !value); }

    /// <summary>ssh in batch mode (the agent's keys, no password prompt) and run `caps --version` on the host.</summary>
    public async Task TestHost()
    {
        if (_host?.Host is not { } h || _hostTesting) return;
        if (h.Hostname.Length == 0) { HostTestText = "Enter the hostname first"; return; }
        HostIdle = false;
        HostTestText = "Connecting…";
        var row = _host;
        try
        {
            var psi = new ProcessStartInfo("ssh") { RedirectStandardOutput = true, RedirectStandardError = true, UseShellExecute = false };
            foreach (var a in new[] { "-o", "BatchMode=yes", "-o", "ConnectTimeout=8", "-p", h.Port.ToString(CultureInfo.InvariantCulture),
                                      h.User.Length > 0 ? $"{h.User}@{h.Hostname}" : h.Hostname, "caps --version || echo caps-not-found" })
                psi.ArgumentList.Add(a);
            using var p = Process.Start(psi) ?? throw new InvalidOperationException("cannot start ssh");
            var outTask = p.StandardOutput.ReadToEndAsync();
            var errTask = p.StandardError.ReadToEndAsync();
            var done = await Task.Run(() => p.WaitForExit(15000));
            if (!done) { try { p.Kill(); } catch { } throw new TimeoutException("no answer in 15 s"); }
            var output = (await outTask).Trim();
            var error = (await errTask).Trim();
            if (p.ExitCode == 0 && output.StartsWith("caps "))
            {
                row.State = "connected";
                row.Level = 1;
                HostTestText = $"Reachable · {output} on the host";
            }
            else if (p.ExitCode == 0)
            {
                row.State = "no caps";
                row.Level = 2;
                HostTestText = "Reachable, but caps is not on the host's PATH: install the CAPS command-line tools there";
            }
            else
            {
                row.State = "unreachable";
                row.Level = 2;
                HostTestText = "ssh failed: " + (error.Length > 0 ? error.Split('\n')[0] : $"exit {p.ExitCode}") + " (keys come from your SSH agent; no passwords are asked or stored)";
            }
        }
        catch (Exception e)
        {
            row.State = "unreachable";
            row.Level = 2;
            HostTestText = e.Message;
        }
        finally
        {
            HostIdle = true;
            Raise(nameof(HostsStatus));
        }
    }

    /// <summary>The batch-job template ({job}, {partition}, {workdir}, {recipe} are filled per job).</summary>
    public string JobTemplate { get => _settings.JobTemplate; set { if (value != null && value != _settings.JobTemplate) { _settings.JobTemplate = value; Raise(); Changed("Job template"); } } }
    public void ResetJobTemplate() { JobTemplate = RemoteHost.DefaultTemplate; }
}
