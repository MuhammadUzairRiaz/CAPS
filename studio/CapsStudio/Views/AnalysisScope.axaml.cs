using Avalonia;
using Avalonia.Controls;
using Avalonia.Markup.Xaml;

namespace CapsStudio.Views;

/// <summary>Frames, atoms and (optionally) the axis an analysis page runs on — bound to Analyze's own settings, so the
/// page and Analyze's side panel always agree.</summary>
public partial class AnalysisScope : UserControl
{
    public static readonly StyledProperty<bool> ShowGroupProperty = AvaloniaProperty.Register<AnalysisScope, bool>(nameof(ShowGroup), true);
    public static readonly StyledProperty<bool> ShowAxisProperty = AvaloniaProperty.Register<AnalysisScope, bool>(nameof(ShowAxis));
    public bool ShowGroup { get => GetValue(ShowGroupProperty); set => SetValue(ShowGroupProperty, value); }
    public bool ShowAxis { get => GetValue(ShowAxisProperty); set => SetValue(ShowAxisProperty, value); }

    public AnalysisScope()
    {
        AvaloniaXamlLoader.Load(this);
        Apply();
        PropertyChanged += (_, e) => { if (e.Property == ShowGroupProperty || e.Property == ShowAxisProperty) Apply(); };
    }

    private void Apply()
    {
        this.FindControl<Grid>("GroupRow")!.IsVisible = ShowGroup;
        this.FindControl<StackPanel>("AxisRow")!.IsVisible = ShowAxis;
    }
}
