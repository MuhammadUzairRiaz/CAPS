using Avalonia.Controls;
using Avalonia.Interactivity;
using Avalonia.Markup.Xaml;

namespace CapsStudio.Views.Pages;

public partial class GrowPage : PageBase
{
    public GrowPage() => AvaloniaXamlLoader.Load(this);
    /// <summary>Where the 3D view goes while this page shows.</summary>
    public Decorator Slot => this.FindControl<Decorator>("ViewSlot")!;

    private async void OnGrow(object? s, RoutedEventArgs e) => await Vm.Grow();
    private void OnCancel(object? s, RoutedEventArgs e) => Vm.CancelGrow();
    private void OnDensityMode(object? s, RoutedEventArgs e) => Vm.GrowUseBox = false;
    private void OnDisplay(object? s, RoutedEventArgs e) { if (s is Control { Tag: string t }) Vm.StyleIndex = int.Parse(t); }
    private async void OnSave(object? s, RoutedEventArgs e) { if (Window != null) await Window.SaveAsAsync("data", "LAMMPS data"); }
    private async void OnCopyCommand(object? s, RoutedEventArgs e)
    {
        var clip = TopLevel.GetTopLevel(this)?.Clipboard;
        if (clip != null) await clip.SetTextAsync(Vm.GrowCommand);
        Vm.Status = "Copied: " + Vm.GrowCommand;
    }
}
