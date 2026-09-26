using Avalonia;
using Avalonia.Controls;
using Avalonia.Markup.Xaml;

namespace CapsStudio.Views;

/// <summary>The "Structure" box of a module page: the project's structures, the active one selected.</summary>
public partial class StructurePicker : UserControl
{
    public static readonly StyledProperty<bool> ShowResultOptionProperty = AvaloniaProperty.Register<StructurePicker, bool>(nameof(ShowResultOption));
    /// <summary>Run modules (Minimise, Equilibrate, Dynamics) offer to keep the original structure.</summary>
    public bool ShowResultOption { get => GetValue(ShowResultOptionProperty); set => SetValue(ShowResultOptionProperty, value); }

    public StructurePicker()
    {
        AvaloniaXamlLoader.Load(this);
        this.FindControl<CheckBox>("ResultBox")!.IsVisible = false;
        PropertyChanged += (_, e) => { if (e.Property == ShowResultOptionProperty) this.FindControl<CheckBox>("ResultBox")!.IsVisible = ShowResultOption; };
    }
}
