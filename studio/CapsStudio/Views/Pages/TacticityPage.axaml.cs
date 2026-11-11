using System.Linq;
using System.Text.Json.Nodes;
using Avalonia;
using Avalonia.Controls;
using Avalonia.Interactivity;
using Avalonia.Markup.Xaml;
using CapsStudio.Interop;
using CapsStudio.ViewModels;

namespace CapsStudio.Views.Pages;

public partial class TacticityPage : PageBase
{
    // the dyad colours: the strip, its legend dots and the units in the view use these same two
    private const string MesoColour = "#F0A83C", RacemoColour = "#2271DB";
    private const int Units = 8;

    private MainViewModel? _hooked;
    public TacticityPage()
    {
        AvaloniaXamlLoader.Load(this);
        this.FindControl<SequenceStrip>("Dyads")!.Palette = [MesoColour, RacemoColour];
        // the view's empty state: said in words, so a dark viewport is never left unexplained
        _empty = new TextBlock
        {
            Text = NoChainText,
            TextWrapping = Avalonia.Media.TextWrapping.Wrap, TextAlignment = Avalonia.Media.TextAlignment.Center, MaxWidth = 320,
            HorizontalAlignment = Avalonia.Layout.HorizontalAlignment.Center, VerticalAlignment = Avalonia.Layout.VerticalAlignment.Center,
            Margin = new Thickness(16), IsHitTestVisible = false,
        };
        _empty.Classes.Add("dim");
        var view = this.FindControl<MolView>("View")!;
        if (view.Parent is Panel host) host.Children.Add(_empty);
        ToolTip.SetTip(view, "Repeat units 1 to 8 of the chain, made whole across the cell boundary. Unit 1 keeps its element colours; "
                             + "each later unit is amber when its dyad with the unit before is meso and blue when it is racemo, as in the dyad strip.");
        DataContextChanged += (_, _) =>
        {
            if (DataContext is not MainViewModel vm || vm == _hooked) return;
            _hooked = vm;
            vm.TsChanged += () => Update(vm);
            vm.PropertyChanged += (_, e) =>
            {
                // a new chain is cut at once; its dyad colours follow when the statistics are in (TsChanged)
                if (e.PropertyName == nameof(MainViewModel.TsDoc)) SyncUnits(vm, false);
                else if (e.PropertyName == nameof(MainViewModel.TsBuilding)) ShowEmpty(vm);
            };
            SyncUnits(vm, false);
        };
    }

    private const string NoChainText = "No chain to show yet. A chain is grown from the stereo model as you change it, and Redraw grows another.";
    private readonly TextBlock _empty;
    private void ShowEmpty(MainViewModel vm)
    {
        _empty.Text = vm.TsDoc == null ? NoChainText : "The first eight repeat units could not be cut from this chain: its atoms carry no unit numbers.";
        _empty.IsVisible = _units == null && !vm.TsBuilding;
    }

    private void Update(MainViewModel vm)
    {
        this.FindControl<SequenceStrip>("Dyads")!.SetUnits(vm.TsDyads);
        this.FindControl<SciPlot>("Pentads")!.Figure = ViewModels.PageFigures.Pentads(vm.TsPentadNames, vm.TsModelPentads, vm.TsMarkov ? "Markov model" : "Bernoulli model", vm.TsChainPentads);
        SyncUnits(vm, true);
    }

    // The panel's own document: only the atoms of repeat units 1 to 8 of the grown chain, unwrapped and without a cell,
    // so the renderer's fit frames those atoms and nothing else of the chain is in the scene.
    private CapsDocument? _units;
    private CapsDocument? _unitsSource;   // the chain (vm.TsDoc) _units was cut from; compared by reference only

    private void SyncUnits(MainViewModel vm, bool colour)
    {
        var v = this.FindControl<MolView>("View")!;
        v.DrawStyle = 0;
        v.ShowCell = false;
        if (!ReferenceEquals(vm.TsDoc, _unitsSource))
        {
            _unitsSource = vm.TsDoc;
            var old = _units;
            _units = _unitsSource == null ? null : CutFirstUnits(_unitsSource);
            v.Document = _units;
            v.Reset();   // the default camera: the renderer fits the view to this document's atoms
            old?.Dispose();
        }
        if (colour && _units != null && vm.TsDoc != null)
        {
            ColourDyads(vm.TsDoc, _units);
            v.Refresh();
        }
        ShowEmpty(vm);
    }

    /// <summary>
    /// A new document holding the atoms of repeat units 1 to <see cref="Units"/> of the chain: a copy of its current
    /// state, every other atom deleted, then the segment made whole across the periodic boundary and the cell dropped
    /// (the builder's cluster edit with radius 0). Null when the chain has no unit numbers or the cut fails.
    /// </summary>
    private static CapsDocument? CutFirstUnits(CapsDocument chain)
    {
        CapsDocument? doc = null;
        try
        {
            var res = chain.AtomResidues();
            if (!res.Any(r => r >= 1 && r <= Units)) return null;
            var drop = Enumerable.Range(0, res.Length).Where(i => res[i] < 1 || res[i] > Units).ToArray();
            doc = chain.StateDocument("{}", $"first {Units} units");
            if (drop.Length > 0)
            {
                var r = JsonNode.Parse(doc.Edit("{\"op\":\"delete\",\"atoms\":[" + string.Join(",", drop) + "]}"));
                if ((bool?)r?["ok"] != true) { doc.Dispose(); return null; }
            }
            // whole across the boundary, no cell; a chain grown without a cell is whole already, so a refusal is fine
            try { doc.Edit("{\"op\":\"cluster\",\"radius\":0}"); } catch { }
            return doc;
        }
        catch
        {
            doc?.Dispose();
            return null;
        }
    }

    /// <summary>
    /// Colours the cut units by the chain's dyads: the dyad between stereocentres k and k+1 colours the unit that holds
    /// centre k+1, amber for meso and blue for racemo (the strip's palette). The first unit has no dyad before it and keeps
    /// its element colours. The centres and dyads are the chain's own (first chain of its tacticity report), the same
    /// source as the strip.
    /// </summary>
    private static void ColourDyads(CapsDocument chain, CapsDocument units)
    {
        try
        {
            var tac = JsonNode.Parse(chain.Tacticity());
            var first = (tac?["chains"] as JsonArray)?.OfType<JsonObject>().FirstOrDefault();
            var centres = (first?["centres"] as JsonArray ?? []).Select(x => (int)(double)x!).ToArray();
            var dyads = (string?)first?["dyads"] ?? "";
            var res = chain.AtomResidues();
            var meso = new HashSet<int>();
            var racemo = new HashSet<int>();
            for (var k = 0; k < dyads.Length && k + 1 < centres.Length; k++)
            {
                int a = centres[k], b = centres[k + 1];
                if (a < 0 || b < 0 || a >= res.Length || b >= res.Length) continue;
                int ra = res[a], rb = res[b];
                if (ra < 1 || ra > Units || rb < 1 || rb > Units) continue;
                (dyads[k] == 'm' ? meso : racemo).Add(rb);
            }
            var ures = units.AtomResidues();
            var steps = new List<string>();
            void Paint(HashSet<int> unitSet, string colour)
            {
                var atoms = Enumerable.Range(0, ures.Length).Where(i => unitSet.Contains(ures[i])).ToArray();
                if (atoms.Length == 0) return;
                steps.Add("{\"type\":\"manual_selection\",\"mode\":\"replace\",\"atoms\":[" + string.Join(",", atoms) + "]}");
                steps.Add("{\"type\":\"assign_colour\",\"colour\":\"" + colour + "\"}");
            }
            Paint(meso, MesoColour);
            Paint(racemo, RacemoColour);
            units.SetPipeline(steps.Count == 0 ? null : "{\"steps\":[" + string.Join(",", steps) + "]}");
        }
        catch
        {
            // colours are a reading aid: without them the units show in element colours
            try { units.SetPipeline(null); } catch { }
        }
    }

    protected override void OnPropertyChanged(AvaloniaPropertyChangedEventArgs change)
    {
        base.OnPropertyChanged(change);
        if (change.Property == IsVisibleProperty && IsVisible && DataContext is MainViewModel vm) Update(vm);
    }
    private void OnRedraw(object? s, RoutedEventArgs e) => Vm.TsRedraw();
    private void OnApply(object? s, RoutedEventArgs e) => Vm.TsApply();
    private void OnFit(object? s, RoutedEventArgs e) => Vm.TsFitMeasured();
}
