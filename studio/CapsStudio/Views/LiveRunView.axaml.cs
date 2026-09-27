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
                if (e.PropertyName != nameof(ViewModels.MainViewModel.RunLiveDoc)) return;
                var v = this.FindControl<MolView>("View")!;
                var first = v.Document == null;
                v.ColourMode = 1; v.ShowCell = true; v.DrawStyle = vm.StyleIndex;
                v.Document = vm.RunLiveDoc;   // the camera stays where the user turned it between snapshots
                if (first && vm.RunLiveDoc != null) v.Reset();
            };
        };
    }
}
