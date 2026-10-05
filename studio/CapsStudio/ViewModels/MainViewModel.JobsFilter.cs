using System.Collections.ObjectModel;
using System.Collections.Specialized;
using System.ComponentModel;

namespace CapsStudio.ViewModels;

/// <summary>The Jobs list filtered and sorted (D12): a search over the job's title and where it runs, its state, and
/// newest, oldest, longest or by kind; the list follows new jobs and state changes as they happen.</summary>
public sealed partial class MainViewModel
{
    public ObservableCollection<Job> JobsShown { get; } = new();
    public static readonly string[] JobStates = ["Every job", "Running or queued", "Finished", "Failed", "Stopped or cancelled"];
    public static readonly string[] JobSorts = ["Newest first", "Oldest first", "Longest first", "By kind"];
    private string _jobQuery = "";
    private int _jobState, _jobSort;
    private bool _jobsViewHooked;
    public string JobQuery { get => _jobQuery; set { if (Set(ref _jobQuery, value ?? "")) RefreshJobsShown(); } }
    public int JobState { get => _jobState; set { if (Set(ref _jobState, Math.Clamp(value, 0, JobStates.Length - 1))) RefreshJobsShown(); } }
    public int JobSort { get => _jobSort; set { if (Set(ref _jobSort, Math.Clamp(value, 0, JobSorts.Length - 1))) RefreshJobsShown(); } }
    public string JobsShownText => JobsShown.Count == Jobs.Count ? "" : $"{JobsShown.Count} of {Jobs.Count} shown";

    private void HookJobsView()
    {
        if (_jobsViewHooked) return;
        _jobsViewHooked = true;
        foreach (var j in Jobs) j.PropertyChanged += OnJobStateChanged;
        Jobs.CollectionChanged += (_, e) =>
        {
            foreach (var j in e.NewItems?.OfType<Job>() ?? []) j.PropertyChanged += OnJobStateChanged;
            foreach (var j in e.OldItems?.OfType<Job>() ?? []) j.PropertyChanged -= OnJobStateChanged;
            if (e.Action == NotifyCollectionChangedAction.Reset) foreach (var j in Jobs) { j.PropertyChanged -= OnJobStateChanged; j.PropertyChanged += OnJobStateChanged; }
            RefreshJobsShown();
        };
    }

    private void OnJobStateChanged(object? s, PropertyChangedEventArgs e)
    {
        if (e.PropertyName is nameof(Job.Status) or nameof(Job.Ended)) RefreshJobsShown();
    }

    public void RefreshJobsShown()
    {
        HookJobsView();
        IEnumerable<Job> q = Jobs;
        var text = _jobQuery.Trim();
        if (text.Length > 0)
            q = q.Where(j => j.Title.Contains(text, StringComparison.OrdinalIgnoreCase) || j.Where.Contains(text, StringComparison.OrdinalIgnoreCase) ||
                             j.Kind.Contains(text, StringComparison.OrdinalIgnoreCase) || j.Document.Contains(text, StringComparison.OrdinalIgnoreCase));
        q = _jobState switch
        {
            1 => q.Where(j => j.IsRunning || j.IsQueued),
            2 => q.Where(j => j.IsDone),
            3 => q.Where(j => j.IsFailed),
            4 => q.Where(j => j.IsQuiet),
            _ => q,
        };
        TimeSpan Length(Job j) => (j.Ended ?? DateTime.Now) - j.Started;
        q = _jobSort switch
        {
            1 => q.OrderBy(j => j.Started),
            2 => q.OrderByDescending(Length),
            3 => q.OrderBy(j => j.Kind).ThenByDescending(j => j.Started),
            _ => q.OrderByDescending(j => j.Started),
        };
        var list = q.ToList();
        if (!list.SequenceEqual(JobsShown))
        {
            var keep = SelectedJob;
            JobsShown.Clear();
            foreach (var j in list) JobsShown.Add(j);
            if (keep != null && list.Contains(keep) && SelectedJob != keep) SelectedJob = keep;
        }
        Raise(nameof(JobsShownText));
    }
}
