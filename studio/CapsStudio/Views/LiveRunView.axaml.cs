using Avalonia.Controls;
using Avalonia.Markup.Xaml;

namespace CapsStudio.Views;

public partial class LiveRunView : UserControl
{
    private ViewModels.MainViewModel? _hooked;

    public LiveRunView()
    {
        AvaloniaXamlLoader.Load(this);
        DataContextChanged += (_, _) =>
        {
            if (DataContext is not ViewModels.MainViewModel vm || vm == _hooked) return;
            _hooked = vm;
            vm.PropertyChanged += (_, e) =>
            {
                if (e.PropertyName is nameof(ViewModels.MainViewModel.RunLiveDoc) or nameof(ViewModels.MainViewModel.Document)
                    or nameof(ViewModels.MainViewModel.Idle)) Show(vm);
            };
            Show(vm);
        };
    }

    /// <summary>The run's latest snapshot; with none, the structure itself while nothing runs (a running MD holds the
    /// document, so it is never drawn from here then).</summary>
    private void Show(ViewModels.MainViewModel vm)
    {
        var v = this.FindControl<MolView>("View")!;
        var hint = this.FindControl<TextBlock>("Hint")!;
        var doc = vm.RunLiveDoc ?? (vm.Idle ? vm.Document : null);
        hint.IsVisible = doc == null;
        v.IsVisible = doc != null;
        if (ReferenceEquals(v.Document, doc)) return;
        var first = v.Document == null;
        v.ColourMode = 1; v.ShowCell = true; v.DrawStyle = vm.StyleIndex;
        v.Document = doc;   // the camera stays where the user turned it between snapshots
        if (first && doc != null) v.Reset();
    }
}
