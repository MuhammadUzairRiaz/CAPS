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
    public double DockChevronAngle => _dockOpen ? 0 : 180;
    public bool ProjectPanelShown => IsStudio && (!_compact || _projectDrawer);
    public bool InspectorShown => !_compact || _inspectorDrawer;
    public event Action? CompactChanged;

    private void RaiseCompact()
    {
        foreach (var n in new[] { nameof(InspectorDrawer), nameof(ProjectDrawer), nameof(DockOpen), nameof(DockChevronAngle), nameof(ProjectPanelShown), nameof(InspectorShown) }) Raise(n);
        CompactChanged?.Invoke();
    }
}
