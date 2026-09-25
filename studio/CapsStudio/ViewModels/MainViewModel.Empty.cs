namespace CapsStudio.ViewModels;

/// <summary>Nothing open (design/boards/VisEmpty): the pages that work on a structure show where to get one.</summary>
public sealed partial class MainViewModel
{
    public string EmptyHeading => "Open a structure or trajectory";
    public bool ShowEmpty => _doc == null && _module is 1 or 20;
    private string _emptyError = "";
    public string EmptyError { get => _emptyError; set { if (Set(ref _emptyError, value)) Raise(nameof(EmptyHasError)); } }
    public bool EmptyHasError => _emptyError.Length > 0;
    public List<RecentItem> RecentFew => Recent.Take(4).ToList();
    public string RecentCount => Recent.Count.ToString(System.Globalization.CultureInfo.InvariantCulture);

    private static string Mod => OperatingSystem.IsMacOS() ? "⌘" : "Ctrl+";
    public string KeyOpen => Mod + "O";
    public string KeyOpenRecent => OperatingSystem.IsMacOS() ? "⌘⇧O" : "Ctrl+Shift+O";
    public string KeyPalette => Mod + "K";
    public string KeyAnnounce => OperatingSystem.IsMacOS() ? "⌘⇧A" : "Ctrl+Shift+A";
}
