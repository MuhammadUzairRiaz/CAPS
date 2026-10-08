namespace CapsStudio.ViewModels;

/// <summary>One stage of a packing (design/boards/Pack): pending, running, done or failed.</summary>
public sealed class PackStageRow : ObservableObject
{
    public string Title { get; init; } = "";
    private string _state = "pending";
    public string State { get => _state; set { if (Set(ref _state, value)) Raise(nameof(Brush)); } }
    public Avalonia.Media.IBrush Brush => Tokens.Brush(_state switch { "running" => "AccB", "done" => "OkB", "failed" => "ErrB", _ => "DimB" });
}
