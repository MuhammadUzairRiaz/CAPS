namespace CapsStudio.ViewModels;

/// <summary>A step of the Studio tour: the region it lights and what it says.</summary>
public sealed record TourStepInfo(string Region, string Name, string Outline, string Title, string Body, string[] Keys);

/// <summary>First-run tour (design/boards/FirstRunTour): five steps over the Studio — modules, tools, the 3D view, the
/// inspector, curves and runs — each lighting its region; shown once when the first structure opens, and from ⌘K.</summary>
public sealed partial class MainViewModel
{
    public static readonly TourStepInfo[] TourSteps =
    [
        new("Rail", "Modules", "Studio, then Grow, Pack, Relax … in workflow order", "The workflow runs down the rail",
            "Studio holds the structure; Grow, Pack, Relax, Dynamics, Equilibrate and React change it; Analyze measures it; Field types its atoms. Jobs keeps every run, Bench checks this machine.",
            ["⌘O open", "⌘W close"]),
        new("Toolbar", "Tools", "select, build, measure, clean", "Tools act on what you select",
            "Pick Select or a builder tool — place an atom, draw a bond, delete, add hydrogens — then click in the 3D view. Every edit can be undone and recorded as Python in the Macro recorder.",
            ["⌘Z undo", "⇧E element", "⌘K commands"]),
        new("ViewHost", "3D view", "orbit, pan, zoom · live monitors", "The view is yours to turn",
            "Drag to rotate, Shift-drag to pan, scroll to zoom. Click an atom to inspect it; Shift-click up to four to measure a distance, angle or dihedral. R resets the view.",
            ["R reset", "⌘E export"]),
        new("InspectorPanel", "Inspector", "atom, molecule, view, file", "What you picked, and how it is drawn",
            "The inspector shows the picked atom and its molecule, the view options and the file's own settings. Selection, Appearance and Checks open here when you choose them in the toolbar.",
            []),
        new("AnalysisDock", "Curves", "g(r), minimisation, dynamics, chains", "Runs report as they go",
            "Curves fill in here while a run works: g(r), energy and force during Relax, temperature and density during Dynamics. Every run is also kept in Jobs with its log and provenance.",
            []),
    ];

    private bool _tourOpen;
    private int _tourStep;
    public bool TourOpen { get => _tourOpen; private set => Set(ref _tourOpen, value); }
    public int TourStep { get => _tourStep; private set { if (Set(ref _tourStep, value)) RaiseTour(); } }
    public TourStepInfo TourCurrent => TourSteps[_tourStep];
    public string TourCount => $"{_tourStep + 1} of {TourSteps.Length}";
    public bool TourFirst => _tourStep == 0;
    public bool TourLast => _tourStep == TourSteps.Length - 1;
    public string TourNextText => TourLast ? "Done" : "Next";
    public event Action? TourChanged;

    private void RaiseTour()
    {
        foreach (var n in new[] { nameof(TourCurrent), nameof(TourCount), nameof(TourFirst), nameof(TourLast), nameof(TourNextText) }) Raise(n);
        TourChanged?.Invoke();
    }

    public void StartTour()
    {
        if (!IsStudio) SetModule(8);
        _tourStep = 0;
        TourOpen = true;
        RaiseTour();
    }

    /// <summary>The tour starts by itself once, with the first structure in the Studio.</summary>
    public void MaybeStartTour()
    {
        if (!_settings.TourDone && HasDocument && IsStudio && !TourOpen && Environment.GetEnvironmentVariable("CAPS_NO_TOUR") is not { Length: > 0 }) StartTour();
    }

    public void TourNext() { if (TourLast) EndTour(); else TourStep++; }
    public void TourBack() { if (!TourFirst) TourStep--; }
    public void EndTour()
    {
        TourOpen = false;
        if (!_settings.TourDone) { _settings.TourDone = true; Changed("Tour"); }
        TourChanged?.Invoke();
    }
}
