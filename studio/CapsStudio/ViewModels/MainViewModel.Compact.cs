namespace CapsStudio.ViewModels;

/// <summary>Compact layout (design/boards/CompactLayout): below 1440 px the side panels become drawers over the view,
/// the rail and the toolbar's view options show icons only, and the curves dock folds to its tabs.</summary>
public sealed partial class MainViewModel
{
    private bool _compact, _inspectorDrawer = true, _projectDrawer, _dockOpen;
    public bool Compact { get => _compact; set { if (Set(ref _compact, value)) RaiseCompact(); } }
    public bool InspectorDrawer { get => _inspectorDrawer; set { if (Set(ref _inspectorDrawer, value)) { if (value) _projectDrawer = false; RaiseCompact(); } } }
    public bool ProjectDrawer { get => _projectDrawer; set { if (Set(ref _projectDrawer, value)) { if (value) _inspectorDrawer = false; RaiseCompact(); } } }
    public bool DockOpen { get => _dockOpen; set { if (Set(ref _dockOpen, value)) RaiseCompact(); } }
    // Fragments: the project drawer on its fragments tab; Monitors: the pinned measurements over the view, shown or hidden
    public bool FragmentsDrawer
    {
        get => _projectDrawer && LeftTab == 1;
        set { if (value) { LeftTab = 1; ProjectDrawer = true; } else if (LeftTab == 1) ProjectDrawer = false; RaiseCompact(); }
    }
    public bool ProjectOnlyDrawer => _projectDrawer && LeftTab != 1;
    private bool _monitorsShown = true;
    public bool MonitorsShown { get => _monitorsShown; set { if (Set(ref _monitorsShown, value)) Raise(nameof(MonitorsVisible)); } }
    public bool MonitorsVisible => HasMonitors && _monitorsShown;
    public double DockChevronAngle => _dockOpen ? 0 : 180;
    /// <summary>The Project Explorer: on the left with every module (a drawer in the compact layout).</summary>
    public bool ProjectPanelShown => !_compact || _projectDrawer;
    private double _projectWidth = 300;
    /// <summary>The explorer's width (its right edge drags it), kept in the settings.</summary>
    public double ProjectWidth { get => _projectWidth; set => Set(ref _projectWidth, Math.Clamp(value, 220, 560)); }
    public bool InspectorShown => !_compact || _inspectorDrawer;
    public event Action? CompactChanged;

    private void RaiseCompact()
    {
        foreach (var n in new[] { nameof(InspectorDrawer), nameof(ProjectDrawer), nameof(FragmentsDrawer), nameof(ProjectOnlyDrawer), nameof(DockOpen), nameof(DockChevronAngle), nameof(ProjectPanelShown), nameof(InspectorShown) }) Raise(n);
        CompactChanged?.Invoke();
    }
}
