using Avalonia;
using Avalonia.Controls;
using Avalonia.Media;
using Avalonia.Styling;

namespace CapsStudio;

/// <summary>Design tokens from code (Themes/Tokens.axaml): the brush of the current theme (Graphite or Paper).</summary>
public static class Tokens
{
    public static IBrush Brush(string key)
    {
        var app = Application.Current;
        if (app != null && app.TryFindResource(key, app.ActualThemeVariant, out var r) && r is IBrush b) return b;
        return Brushes.Gray;
    }

    public static FontFamily Sans => Application.Current?.TryFindResource("Sans", out var r) == true && r is FontFamily f ? f : FontFamily.Default;

    public static FontFamily Mono => Application.Current?.TryFindResource("Mono", out var r) == true && r is FontFamily f ? f : FontFamily.Default;

    /// <summary>Switches between Graphite (dark) and Paper (light).</summary>
    public static void Use(bool light)
    {
        if (Application.Current != null) Application.Current.RequestedThemeVariant = light ? ThemeVariant.Light : ThemeVariant.Dark;
    }
}
