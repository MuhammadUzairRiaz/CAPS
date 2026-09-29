using Avalonia.Controls;
using Avalonia.Markup.Xaml;

namespace CapsStudio.Views;

/// <summary>The label settings (MainViewModel's label kinds, scope, font, size and colours).</summary>
public partial class LabelSettings : UserControl
{
    public LabelSettings() => AvaloniaXamlLoader.Load(this);
}
