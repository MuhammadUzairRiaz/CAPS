using Avalonia.Controls;
using Avalonia.Interactivity;
using Avalonia.Markup.Xaml;

namespace CapsStudio.Views;

/// <summary>A run's seed (SeedChoice): fixed for reproducibility, or new each time with the seed used shown.</summary>
public partial class SeedPicker : UserControl
{
    public static readonly Avalonia.StyledProperty<bool> ShowLabelProperty = Avalonia.AvaloniaProperty.Register<SeedPicker, bool>(nameof(ShowLabel), true);
    /// <summary>False where the page puts its own "Seed" label above (a grid cell).</summary>
    public bool ShowLabel { get => GetValue(ShowLabelProperty); set => SetValue(ShowLabelProperty, value); }

    public SeedPicker()
    {
        AvaloniaXamlLoader.Load(this);
        PropertyChanged += (_, e) => { if (e.Property == ShowLabelProperty) this.FindControl<TextBlock>("Label")!.IsVisible = ShowLabel; };
    }
    private void OnKeep(object? s, RoutedEventArgs e) => (DataContext as ViewModels.SeedChoice)?.Keep();
}
