using Avalonia.Controls;
using Avalonia.Markup.Xaml;

namespace CapsStudio.Views;

/// <summary>The cluster options of a run page, shown when Run where is a host.</summary>
public partial class RemoteRunPanel : UserControl
{
    public RemoteRunPanel() { AvaloniaXamlLoader.Load(this); }
}
