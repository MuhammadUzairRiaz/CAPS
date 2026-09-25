using System.Collections.ObjectModel;
using System.ComponentModel;
using System.Globalization;
using System.Text.Json.Nodes;
using Avalonia.Media;
using Avalonia.Media.Imaging;

namespace CapsStudio.ViewModels;

/// <summary>One colouring in the gallery: its preview and legend.</summary>
public sealed class ColourTile : INotifyPropertyChanged
{
    public event PropertyChangedEventHandler? PropertyChanged;
    private void Raise(string n) => PropertyChanged?.Invoke(this, new PropertyChangedEventArgs(n));
    public string Title { get; init; } = "";
    public string Property { get; init; } = "";
    public string Mode { get; init; } = "auto";
    public string Map { get; init; } = "viridis";
    private string _sub = "", _lo = "", _hi = "";
    private Bitmap? _image;
    private bool _continuous;
    public string Sub { get => _sub; set { _sub = value; Raise(nameof(Sub)); } }
    public Bitmap? Image { get => _image; set { _image = value; Raise(nameof(Image)); } }
    public bool Continuous { get => _continuous; set { _continuous = value; Raise(nameof(Continuous)); Raise(nameof(Categorical)); } }
    public bool Categorical => !_continuous;
    public string Lo { get => _lo; set { _lo = value; Raise(nameof(Lo)); } }
    public string Hi { get => _hi; set { _hi = value; Raise(nameof(Hi)); } }
    public bool Diverging => Map == "diverging";
    public bool Viridis => Map != "diverging";
    public ObservableCollection<LegendSwatch> Entries { get; } = new();
}

/// <summary>Analyze › Visualize › Colour by (design/boards/ColourBy): one cell, six colourings; choosing one sets the
/// Visualize pipeline's colour coding.</summary>
public sealed partial class MainViewModel
{
    public bool IsColourBy => _module == 24;
    public const int ColourTileW = 432, ColourTileH = 250;

    public ObservableCollection<ColourTile> ColourTiles { get; } =
    [
        new ColourTile { Title = "Element", Property = "Element", Mode = "categorical" },
        new ColourTile { Title = "Molecule (chain)", Property = "Molecule", Mode = "categorical" },
        new ColourTile { Title = "Atom type", Property = "Type", Mode = "categorical" },
        new ColourTile { Title = "Distance to chain centre", Property = "DistanceToCOM", Mode = "continuous" },
        new ColourTile { Title = "Position z", Property = "Position.Z", Mode = "continuous" },
        new ColourTile { Title = "Partial charge", Property = "Charge", Mode = "continuous", Map = "diverging" },
    ];
    private int _colourGen;

    public void OpenColourBy()
    {
        if (_doc == null) { Status = "Open a structure first"; return; }
        SetModule(24);
        RenderColourTiles();
    }

    /// <summary>Renders each colouring with a one-step pipeline (off the UI thread), then puts the Visualize pipeline back.</summary>
    private void RenderColourTiles()
    {
        if (_doc == null) return;
        var doc = _doc;
        var gen = ++_colourGen;
        var cam = Camera;
        cam.Zoom *= 1.2;   // the tiles are small: a little closer, as on the board
        var opt = ViewOptions(ColourTileW, ColourTileH, 2);
        opt.Highlight0 = opt.Highlight1 = opt.Highlight2 = opt.Highlight3 = -1;
        opt.Focus = 0;
        opt.ShowCell = 1;
        var tiles = ColourTiles.ToList();
        Task.Run(() =>
        {
            foreach (var t in tiles)
            {
                var step = new JsonObject { ["type"] = "colour_coding", ["property"] = t.Property, ["mode"] = t.Mode, ["map"] = t.Map, ["lighten_h"] = t.Property is "Molecule" or "Type" };
                var rgba = new byte[ColourTileW * ColourTileH * 4];
                string result;
                try
                {
                    doc.SetPipeline(new JsonObject { ["steps"] = new JsonArray(step) }.ToJsonString());
                    doc.Render(cam, opt, rgba);
                    result = doc.PipelineResult();
                }
                catch { continue; }
                Avalonia.Threading.Dispatcher.UIThread.Post(() =>
                {
                    if (gen != _colourGen || !IsColourBy) return;
                    t.Image = ToBitmap(rgba, ColourTileW, ColourTileH);
                    FillColourLegend(t, result);
                });
            }
            try { doc.SetPipeline(""); } catch { }   // the gallery leaves no pipeline behind; Visualize sets its own
        });
    }

    private static void FillColourLegend(ColourTile t, string result)
    {
        t.Entries.Clear();
        if (result.Length == 0 || JsonNode.Parse(result)?["legend"] is not JsonObject l) { t.Sub = "no legend"; return; }
        var inv = CultureInfo.InvariantCulture;
        t.Continuous = (bool?)l["continuous"] ?? false;
        var lo = (double?)l["lo"] ?? 0;
        var hi = (double?)l["hi"] ?? 0;
        t.Lo = lo.ToString("G3", inv);
        t.Hi = hi.ToString("G3", inv);
        var e = l["entries"] as JsonArray;
        if (e != null)
            foreach (var x in e.Take(10))
                t.Entries.Add(new LegendSwatch((t.Property == "Molecule" ? "mol " : "") + ((string?)x?["label"] ?? ""), new SolidColorBrush(Color.Parse((string?)x?["colour"] ?? "#888888"))));
        t.Sub = t.Continuous ? $"{t.Lo} → {t.Hi}" + (t.Map == "diverging" ? " · diverging" : " · viridis") : $"{(e?.Count ?? 0)} values · categorical";
    }

    /// <summary>Uses the chosen colouring: the Visualize pipeline's colour coding takes its property (added if missing).</summary>
    public void ApplyColourTile(ColourTile t)
    {
        var row = PipelineRows.FirstOrDefault(r => r.Type == "colour_coding");
        if (row == null)
        {
            PipeSelected = PipelineRows.FirstOrDefault();
            AddStep("colour_coding");
            row = PipelineRows.First(r => r.Type == "colour_coding");
        }
        row.Params["property"] = t.Property;
        row.Params["mode"] = t.Mode;
        row.Params["map"] = t.Map;
        row.Params["lighten_h"] = t.Property == "Molecule";
        PipeSelected = null;
        PipeSelected = row;
        SetModule(20);
        Status = $"Colour by {t.Title.ToLowerInvariant()}";
    }
}
