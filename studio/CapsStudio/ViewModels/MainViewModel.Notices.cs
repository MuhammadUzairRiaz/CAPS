using System;
using System.Collections.ObjectModel;
using System.Linq;

namespace CapsStudio.ViewModels;

/// <summary>A system message (design/boards/SystemStates): what happened, in plain words; then what still works; then one
/// primary action. Severity: info (nothing to do), check (worth a look), warning (may affect work), blocks (must be
/// fixed first — a run refused).</summary>
public sealed class Notice
{
    public required string Key { get; init; }
    public required string Severity { get; init; }      // info | check | warning | blocks
    public required string Title { get; init; }
    public required string Body { get; init; }
    public string Icon { get; init; } = "alert";
    public string Primary { get; init; } = "";
    public string Secondary { get; init; } = "";
    public Action? OnPrimary { get; init; }
    public Action? OnSecondary { get; init; }
    public string SeverityText => Severity == "blocks" ? "blocks run" : Severity;
    public bool HasPrimary => Primary.Length > 0;
    public bool HasSecondary => Secondary.Length > 0;
    public bool IsInfo => Severity == "info";
    public bool IsCheck => Severity == "check";
    public bool IsWarning => Severity == "warning";
    public bool IsBlocks => Severity == "blocks";
}

public partial class MainViewModel
{
    /// <summary>The messages shown over the view, newest first (at most three).</summary>
    public ObservableCollection<Notice> Notices { get; } = new();
    public bool HasNotices => Notices.Count > 0;

    /// <summary>Shows a notice; one with the same key replaces the old one.</summary>
    public void Notify(Notice n)
    {
        foreach (var old in Notices.Where(x => x.Key == n.Key).ToList()) Notices.Remove(old);
        Notices.Insert(0, n);
        while (Notices.Count > 3) Notices.RemoveAt(Notices.Count - 1);
        Raise(nameof(HasNotices));
        Announcement = $"{n.Title}. {n.Body}";
    }

    public void Dismiss(string key)
    {
        foreach (var old in Notices.Where(x => x.Key == key).ToList()) Notices.Remove(old);
        Raise(nameof(HasNotices));
    }

    public void DismissNotice(Notice n) { Notices.Remove(n); Raise(nameof(HasNotices)); }

    public void NoticePrimary(Notice n) { DismissNotice(n); n.OnPrimary?.Invoke(); }
    public void NoticeSecondary(Notice n) { DismissNotice(n); n.OnSecondary?.Invoke(); }

    /// <summary>Runs refuse an incomplete Field assignment: says which terms are missing and where to fix them.</summary>
    private bool BlockedByField(string run)
    {
        if (!Field.Incomplete) { Dismiss("field.missing"); return false; }
        var missing = Field.Missing.Count;
        var what = missing > 0 ? string.Join(" and ", Field.Missing.GroupBy(m => m.Kind).Select(g => $"{g.Count()} {g.Key}{(g.Count() == 1 ? "" : "s")}")) : "some atoms";
        Notify(new Notice
        {
            Key = "field.missing", Severity = "blocks", Icon = "tag",
            Title = missing > 0 ? $"{missing} missing parameter{(missing == 1 ? "" : "s")}" : "Atoms without a type",
            Body = $"{run} cannot start: {what} have no {Field.ForceFieldName} parameters. CAPS never guesses them.",
            Primary = "Open Field report", OnPrimary = () => SetModule(7),
            Secondary = "Use UFF", OnSecondary = () => { Field.UseUff(); },
        });
        Status = $"{run} blocked: the Field assignment is incomplete";
        return true;
    }
}
