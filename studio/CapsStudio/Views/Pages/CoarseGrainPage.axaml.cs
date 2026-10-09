using Avalonia;
using Avalonia.Controls;
using Avalonia.Interactivity;
using Avalonia.Markup.Xaml;
using Avalonia.Platform.Storage;
using CapsStudio.ViewModels;

namespace CapsStudio.Views.Pages;

public partial class CoarseGrainPage : PageBase
{
    private MainViewModel? _hooked;
    private static readonly string[] Palette = ["AccB", "SelB", "OkB", "WarnB", "ErrB", "MutedB", "TextB"];

    public CoarseGrainPage()
    {
        AvaloniaXamlLoader.Load(this);
        DataContextChanged += (_, _) =>
        {
            if (DataContext is not MainViewModel vm || vm == _hooked) return;
            _hooked = vm;
            vm.CgwChanged += () => { if (vm.IsCgw) Draw(vm); };
        };
    }

    protected override void OnPropertyChanged(AvaloniaPropertyChangedEventArgs change)
    {
        base.OnPropertyChanged(change);
        if (change.Property == IsVisibleProperty && IsVisible && DataContext is MainViewModel vm) Draw(vm);
    }

    private void Draw(MainViewModel vm)
    {
        var kinds = this.FindControl<XyChart>("KindChart")!;
        kinds.Categories = vm.CgwKindCats;
        kinds.Set(vm.CgwKindBars.Select((s, i) => new ChartSeries(s.Name, s.Pts, "bars", Palette[i % Palette.Length])).ToArray());

        var p = this.FindControl<XyChart>("PChart")!;
        var rg = vm.CgwPlotRange;
        p.XLabel = vm.CgwPlotUnit;
        p.XMin = rg?.Lo; p.XMax = rg?.Hi;
        p.Set(vm.CgwPlotP.Length > 0 ? [new ChartSeries("P", vm.CgwPlotP, "line", "AccB")] : []);
        var u = this.FindControl<XyChart>("UChart")!;
        u.XLabel = vm.CgwPlotUnit;
        u.XMin = rg?.Lo; u.XMax = rg?.Hi; u.YMax = rg?.UMax;
        var uPts = rg is { } r0 ? vm.CgwPlotU.Where(q => q.X >= r0.Lo && q.X <= r0.Hi && q.Y <= r0.UMax * 1.5 + 1).ToArray() : vm.CgwPlotU;
        u.Set(uPts.Length > 0 ? [new ChartSeries("U", uPts, "line", "SelB")] : []);

        var gr = this.FindControl<XyChart>("GrChart")!;
        gr.Set(vm.CgwGr.Where(s => s.Value.ContainsKey(vm.CgwPairKey)).Select((s, i) => new ChartSeries(s.Key, s.Value[vm.CgwPairKey], "line", Palette[i % Palette.Length])).ToArray(),
               guides: vm.CgwGr.Count > 0 ? [new ChartGuide(1, false)] : null);
        var pair = this.FindControl<XyChart>("PairChart")!;
        var ps = new List<ChartSeries>();
        if (vm.CgwPairU.Length > 0) ps.Add(new ChartSeries("table " + Path.GetFileName(Path.GetDirectoryName(vm.CgwPairsPath) ?? ""), vm.CgwPairU, "line", "AccB"));
        if (vm.CgwFitU.Length > 0) ps.Add(new ChartSeries("fit " + vm.CgwForm, vm.CgwFitU, "dash", "WarnB"));
        // the repulsive core reaches thousands of kcal/mol: show the well
        var ys = ps.SelectMany(s => s.Points.Select(q => q.Y)).Where(double.IsFinite).ToArray();
        pair.YMin = ys.Length > 0 ? Math.Min(-0.1, ys.Min() * 1.3) : null;
        pair.YMax = ys.Length > 0 ? Math.Max(0.5, -3 * ys.Min()) : null;
        pair.Set(ps.ToArray(), guides: ps.Count > 0 ? [new ChartGuide(0, false)] : null);

        var r2 = this.FindControl<XyChart>("R2Chart")!;
        var rs = new List<ChartSeries>();
        if (vm.CgwR2Built.Length > 0) rs.Add(new ChartSeries("built CG chains", vm.CgwR2Built, "line", "AccB"));
        if (vm.CgwR2Aa.Length > 0) rs.Add(new ChartSeries("all-atom (mapped)", vm.CgwR2Aa, "dash", "TextB"));
        r2.Set(rs.ToArray());

        var st = this.FindControl<XyChart>("StressChart")!;
        st.Set(vm.CgwStress.Select((s, i) => new ChartSeries(s.Name, s.Pts, i == 0 ? "line" : "dash", Palette[i % Palette.Length], i == 0 ? 2.4 : 1.4)).ToArray(),
               guides: vm.CgwStress.Count > 0 ? [new ChartGuide(0, false)] : null);

        var dyn = this.FindControl<XyChart>("DynChart")!;
        dyn.Set(vm.CgwDyn.Where(s => s.Pts.Length > 0).Select((s, i) => new ChartSeries(s.Name, s.Pts, "line", Palette[i % Palette.Length])).ToArray());
    }

    // ---------------------------------------------------------------- pickers (start in the last folder used)
    private async Task<IStorageFolder?> Start(TopLevel top) => Vm.CgwPickDir.Length > 0 && Directory.Exists(Vm.CgwPickDir) ? await top.StorageProvider.TryGetFolderFromPathAsync(Vm.CgwPickDir) : null;
    private async Task<string[]> Files(string title, bool many)
    {
        if (TopLevel.GetTopLevel(this) is not { } top) return [];
        var f = await top.StorageProvider.OpenFilePickerAsync(new FilePickerOpenOptions { Title = title, AllowMultiple = many, SuggestedStartLocation = await Start(top) });
        var paths = f.Select(x => x.TryGetLocalPath()).OfType<string>().ToArray();
        if (paths.Length > 0) Vm.CgwPickDir = Path.GetDirectoryName(paths[0]) ?? "";
        return paths;
    }
    private async Task<string?> Folder(string title)
    {
        if (TopLevel.GetTopLevel(this) is not { } top) return null;
        var d = await top.StorageProvider.OpenFolderPickerAsync(new FolderPickerOpenOptions { Title = title, SuggestedStartLocation = await Start(top) });
        var p = d.Count > 0 ? d[0].TryGetLocalPath() : null;
        if (p != null) Vm.CgwPickDir = p;
        return p;
    }
    private async Task<string?> SaveAs(string title, string name)
    {
        if (TopLevel.GetTopLevel(this) is not { } top) return null;
        var f = await top.StorageProvider.SaveFilePickerAsync(new FilePickerSaveOptions { Title = title, SuggestedFileName = name, SuggestedStartLocation = await Start(top) });
        return f?.TryGetLocalPath();
    }
    /// <summary>Tag "file:Prop", "folder:Prop", "save:Prop" or "clear:Prop": picks a path into that view-model property.</summary>
    private async void OnPickFile(object? s, RoutedEventArgs e)
    {
        if ((s as Control)?.Tag is not string tag || tag.Split(':', 2) is not [var kind, var prop]) return;
        var pi = typeof(MainViewModel).GetProperty(prop);
        if (pi == null) return;
        string? path = kind switch
        {
            "file" => (await Files(Label(prop), false)).FirstOrDefault(),
            "folder" => await Folder(Label(prop)),
            "save" => await SaveAs(Label(prop), Path.GetFileName(pi.GetValue(Vm) as string ?? "calibration.json")),
            "clear" => "",
            _ => null,
        };
        if (path != null) pi.SetValue(Vm, path);
    }
    private static string Label(string prop) => prop switch
    {
        "CgwMapDump" => "An all-atom LAMMPS dump of the structure",
        "CgwTypes" => "The shared type list (types.json)",
        "CgwTargets" => "The IBI targets (targets.json)",
        "CgwPairsPath" => "The current tables (pairs.json)",
        "CgwFits" => "The fit (fits.json)",
        "CgwHistory" => "The calibration history (calibration.json)",
        "CgwStressFile" => "A tension run's STEM.stress_strain.dat",
        "CgwMechMap" or "CgwDynMap" or "CgwBackCg" => "The CG system's map.json",
        "CgwMechDump" or "CgwDynDump" => "The CG run's dump",
        "CgwAaCsv" => "dynamics.csv of the mapped all-atom run",
        "CgwCgCsv" => "dynamics.csv of the CG run",
        "CgwBackRefMap" => "The reference's map.json",
        "CgwBackRefData" => "The reference all-atom LAMMPS data",
        "CgwBackInput" => "The reference LAMMPS input",
        "CgwBackFrame" => "A CG frame (data file or dump)",
        "CgwBackCheckData" => "An all-atom LAMMPS data file to check",
        _ => "Choose a folder",
    };

    private async void OnRoot(object? s, RoutedEventArgs e) { if (await Folder("The coarse-grain working folder") is { } p) Vm.CgwRoot = p; }
    private async void OnAddInputs(object? s, RoutedEventArgs e) => Vm.CgwAddInputs(await Files("All-atom structures (LAMMPS data with bonds, PDB, mol2 …)", true));
    private void OnRemoveInput(object? s, RoutedEventArgs e) { if ((s as Control)?.Tag is string p) Vm.CgwInputs.Remove(p); }
    private void OnRemoveSys(object? s, RoutedEventArgs e)
    {
        if ((s as Control)?.Tag is not CgwSysRow r) return;
        Vm.CgwSystems.Remove(r); Vm.CgwStepRows.Remove(r); Vm.CgwPpaRows.Remove(r);
    }
    /// <summary>A system row: its map.json, then its frames (sys, ppa) or its dump and log (step).</summary>
    private async void OnAddSystem(object? s, RoutedEventArgs e)
    {
        var which = (s as Control)?.Tag as string ?? "sys";
        var map = (await Files("The system's map.json", false)).FirstOrDefault();
        if (map == null) return;
        string[] rest;
        if (which == "step")
        {
            var dump = (await Files("Its CG dump from the last iteration", false)).FirstOrDefault();
            if (dump == null) return;
            var log = (await Files("Its LAMMPS log (cancel: none)", false)).FirstOrDefault();
            rest = log == null ? [dump] : [dump, log];
        }
        else
        {
            rest = await Files("Its frames (CG dumps or STEM.cg.data)", true);
            if (rest.Length == 0) return;
        }
        var row = new CgwSysRow(map, rest);
        (which == "step" ? Vm.CgwStepRows : which == "ppa" ? Vm.CgwPpaRows : Vm.CgwSystems).Add(row);
    }
    private async void OnTgScan(object? s, RoutedEventArgs e)
    {
        var f = (await Files("A cooling scan: 'T density' per line", false)).FirstOrDefault();
        if (f == null) return;
        Vm.CgwTgFile = f;
        await Vm.CgwRunTg();
    }

    private void OnStage(object? s, RoutedEventArgs e) { if ((s as Control)?.Tag is CgwStageItem st) Vm.CgwStage = st.Index; }
    private async void OnPrimary(object? s, RoutedEventArgs e)
    {
        if (Vm.CgwStage == 4) { await Copy(Vm.CgwEquilCommand); return; }
        await Vm.CgwPrimary();
    }
    private async void OnMap(object? s, RoutedEventArgs e) => await Vm.CgwRunMap();
    private async void OnBonded(object? s, RoutedEventArgs e) => await Vm.CgwRunBonded();
    private async void OnTargets(object? s, RoutedEventArgs e) => await Vm.CgwRunTargets();
    private async void OnIbiStart(object? s, RoutedEventArgs e) => await Vm.CgwIbiStart();
    private async void OnIbiStep(object? s, RoutedEventArgs e) => await Vm.CgwIbiStep();
    private async void OnFit(object? s, RoutedEventArgs e) => await Vm.CgwFit();
    private async void OnCalibrate(object? s, RoutedEventArgs e) => await Vm.CgwCalibrate();
    private async void OnBuild(object? s, RoutedEventArgs e) => await Vm.CgwBuild();
    private async void OnPpa(object? s, RoutedEventArgs e) => await Vm.CgwRunPpa();
    private async void OnMechDecks(object? s, RoutedEventArgs e) => await Vm.CgwMechDecks();
    private async void OnMechAnalyze(object? s, RoutedEventArgs e) => await Vm.CgwMechAnalyze();
    private async void OnDyn(object? s, RoutedEventArgs e) => await Vm.CgwRunDyn();
    private async void OnTimeMap(object? s, RoutedEventArgs e) => await Vm.CgwRunTimeMap();
    private async void OnBackmap(object? s, RoutedEventArgs e) => await Vm.CgwBackmap();
    private async void OnBackCheck(object? s, RoutedEventArgs e) => await Vm.CgwBackmapCheck();

    private async Task Copy(string text) { if (text.Length > 0 && TopLevel.GetTopLevel(this)?.Clipboard is { } c) await c.SetTextAsync(text); }
    private async void OnCopyCli(object? s, RoutedEventArgs e) => await Copy(Vm.CgwCli);
    private async void OnCopyText(object? s, RoutedEventArgs e) => await Copy((s as Control)?.Tag as string == "ibi" ? Vm.CgwIbiLoop : Vm.CgwEquilCommand);
}
