namespace CapsStudio.ViewModels;

/// <summary>The status bar's mouse hint, with the modifier keys named as words for the platform (Option on macOS, Alt elsewhere).</summary>
public sealed partial class MainViewModel
{
    private static string AltKeyName => OperatingSystem.IsMacOS() ? "Option" : "Alt";
    public string MouseHint => $"dbl-click molecule · {AltKeyName} drag box · Shift pan";
    public string MouseHintTip => $"Click picks (Shift click adds) · double-click: the molecule · {AltKeyName} double-click: the same type · {AltKeyName} drag: a box · "
                                + $"drag turns · Shift or right drag pans · hold the right button: the ring · H hides · {AltKeyName}+H shows all · Esc clears";
}
