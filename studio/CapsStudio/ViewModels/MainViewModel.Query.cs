using System;
using System.Collections.ObjectModel;
using System.Diagnostics;
using System.Globalization;
using System.Linq;
using System.Text.Json.Nodes;

namespace CapsStudio.ViewModels;

public sealed class SavedQueryRow : ObservableObject
{
    public required string Name { get; init; }
    public required string Query { get; init; }
    private string _count = "—";
    public string Count { get => _count; set => Set(ref _count, value); }
}

/// <summary>Select by query (design/boards/SmartSelect, ⌘F): one grammar — smarts "…", element, type, chain, index,
/// ring, stereo R|S|*, within D of …, sel, and / or / not / ( ) — counted as you type, applied on Enter, saved by name,
/// and kept as a named group.</summary>
public partial class MainViewModel
{
    private bool _queryOpen;
    private string _queryText = "", _queryInfo = "", _queryTiming = "";
    private bool _queryBad;
    public bool QueryOpen { get => _queryOpen; set { if (Set(ref _queryOpen, value) && value) { RefreshSavedQueries(); PreviewQuery(); } } }
    public string QueryText { get => _queryText; set { if (Set(ref _queryText, value ?? "")) PreviewQuery(); } }
    public string QueryInfo { get => _queryInfo; private set => Set(ref _queryInfo, value); }
    public string QueryTiming { get => _queryTiming; private set => Set(ref _queryTiming, value); }
    public bool QueryBad { get => _queryBad; private set => Set(ref _queryBad, value); }
    public ObservableCollection<SavedQueryRow> SavedQueries { get; } = new();

    private JsonNode? QueryCall(string query, string op)
    {
        if (_doc == null) return null;
        return JsonNode.Parse(_doc.Select(new JsonObject { ["mode"] = "query", ["pattern"] = query, ["op"] = op }.ToJsonString()));
    }

    /// <summary>Counts the atoms (and whole rings) the query matches, without selecting them.</summary>
    public void PreviewQuery()
    {
        if (_doc == null || _queryText.Trim().Length == 0) { QueryInfo = _doc == null ? "No structure open" : "Type a query"; QueryBad = false; QueryTiming = ""; return; }
        var sw = Stopwatch.StartNew();
        var r = QueryCall(_queryText, "preview");
        var ms = sw.Elapsed.TotalMilliseconds;
        if (r?["ok"]?.GetValue<bool>() != true) { QueryInfo = r?["error"]?.GetValue<string>() ?? "cannot parse"; QueryBad = true; QueryTiming = ""; return; }
        QueryBad = false;
        var n = (int)r["count"]!.GetValue<double>();
        var rings = (int)(r["rings"]?.GetValue<double>() ?? 0);
        QueryInfo = $"{n:N0} atom{(n == 1 ? "" : "s")}" + (rings > 0 ? $" · {rings} ring{(rings == 1 ? "" : "s")}" : "");
        QueryTiming = $"query parsed · {ms.ToString("0.0", CultureInfo.InvariantCulture)} ms";
    }

    /// <summary>Enter: the query becomes the selection (op: replace, add, subtract, intersect).</summary>
    public void ApplyQuery(string op = "replace")
    {
        var r = QueryCall(_queryText, op);
        if (r?["ok"]?.GetValue<bool>() != true) { QueryInfo = r?["error"]?.GetValue<string>() ?? "cannot parse"; QueryBad = true; return; }
        SelectedCount = (int)r["count"]!.GetValue<double>();
        SelectHud = _queryText.Length > 40 ? _queryText[..40] + "…" : _queryText;
        Record($"doc.select(\"query\", {PyStr(_queryText)}{(op == "replace" ? "" : $", op={PyStr(op)}")})");
        Status = $"Selected {SelectedCount:N0} atoms · {_queryText}";
        RenderRequested?.Invoke();
    }

    public void RefreshSavedQueries()
    {
        SavedQueries.Clear();
        foreach (var q in _settings.SavedQueries)
        {
            var row = new SavedQueryRow { Name = q.Name, Query = q.Query };
            try
            {
                var r = QueryCall(q.Query, "preview");
                row.Count = r?["ok"]?.GetValue<bool>() == true ? ((int)r["count"]!.GetValue<double>()).ToString("N0", CultureInfo.InvariantCulture) : "error";
            }
            catch { row.Count = "—"; }
            SavedQueries.Add(row);
        }
    }

    public void UseSavedQuery(SavedQueryRow r) => QueryText = r.Query;

    public void SaveQuery(string? name = null)
    {
        if (_queryText.Trim().Length == 0 || QueryBad) return;
        name = string.IsNullOrWhiteSpace(name) ? $"Query {_settings.SavedQueries.Count + 1}" : name.Trim();
        _settings.SavedQueries.RemoveAll(q => q.Name == name);
        _settings.SavedQueries.Add(new SavedQueryData { Name = name, Query = _queryText.Trim() });
        Changed("Saved query");
        RefreshSavedQueries();
    }

    public void DeleteSavedQuery(SavedQueryRow r)
    {
        _settings.SavedQueries.RemoveAll(q => q.Name == r.Name);
        Changed("Saved query removed");
        RefreshSavedQueries();
    }

    /// <summary>Selects and keeps the result as a named group (Selection › sets).</summary>
    public void QueryAsGroup()
    {
        ApplyQuery();
        if (!QueryBad) SaveSelectionAsSet(SavedQueries.FirstOrDefault(q => q.Query == _queryText.Trim())?.Name ?? _queryText.Trim());
    }
}
