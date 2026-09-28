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

    /// <summary>Graphite (dark), Paper (light) or the system's choice.</summary>
    public static void UseTheme(string theme)
    {
        if (Application.Current != null)
            Application.Current.RequestedThemeVariant = theme switch { "light" => ThemeVariant.Light, "system" => ThemeVariant.Default, _ => ThemeVariant.Dark };
        Reapply();
    }

    private static bool _highContrast, _focusRings;
    /// <summary>High contrast is on (the GPU view draws its outlines twice as wide and near-full ink).</summary>
    public static bool HighContrast => _highContrast;
    private static Styles? _focusStyle;

    /// <summary>High contrast (Settings › Accessibility): firmer borders and secondary text for the theme in use; call again
    /// after a theme change.</summary>
    public static void UseHighContrast(bool on)
    {
        _highContrast = on;
        var app = Application.Current;
        if (app == null) return;
        foreach (var k in new[] { "LineC", "MutedC", "DimC" }) app.Resources.Remove(k);
        if (!on) return;
        var light = app.ActualThemeVariant == ThemeVariant.Light;
        app.Resources["LineC"] = Color.Parse(light ? "#8A857D" : "#6A737C");
        app.Resources["MutedC"] = Color.Parse(light ? "#2E3338" : "#D4D8DC");
        app.Resources["DimC"] = Color.Parse(light ? "#40454B" : "#BCC2C8");
    }

    /// <summary>Focus rings on every focused control (not only after Tab): an accent border on the focused control.</summary>
    public static void UseFocusRings(bool on)
    {
        _focusRings = on;
        var app = Application.Current;
        if (app == null) return;
        if (_focusStyle != null) { app.Styles.Remove(_focusStyle); _focusStyle = null; }
        if (!on) return;
        var accent = Brush("AccB");
        var st = new Style(x => x.Is<Avalonia.Controls.Primitives.TemplatedControl>().Class(":focus"));
        st.Setters.Add(new Setter(Avalonia.Controls.Primitives.TemplatedControl.BorderBrushProperty, accent));
        st.Setters.Add(new Setter(Avalonia.Controls.Primitives.TemplatedControl.BorderThicknessProperty, new Thickness(2)));
        _focusStyle = new Styles { st };
        app.Styles.Add(_focusStyle);
    }

    /// <summary>After a theme change: the accessibility overrides follow the new theme.</summary>
    public static void Reapply()
    {
        if (_highContrast) UseHighContrast(true);
        if (_focusRings) UseFocusRings(true);
    }

    /// <summary>Switches between Graphite (dark) and Paper (light).</summary>
    public static void Use(bool light)
    {
        if (Application.Current != null) Application.Current.RequestedThemeVariant = light ? ThemeVariant.Light : ThemeVariant.Dark;
    }
}
