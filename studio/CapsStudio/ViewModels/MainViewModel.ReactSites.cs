using System.Globalization;
using System.Text.Json.Nodes;
using CapsStudio.Interop;

namespace CapsStudio.ViewModels;

/// <summary>React › reactive sites and crosslinking: each chain's sites for the chosen reaction (the epoxides of an ENR
/// chain, the CH₂ of a PE chain), how many of them may react per chain, and the crosslinking asked for — as a degree of
/// crosslinking or as phr of a crosslinker — with everything that follows from it: links, crosslinker molecules, phr.
///   DC = 2 N_links / N_units × 100 %  (a link joins two repeat units; Vasilev et al. 2021)
///   phr = 100 N_x M_x / m_rubber       (one crosslinker molecule per link, as MAH bridging two epoxides)</summary>
public sealed partial class MainViewModel
{
    private string _rxSitesText = "";
    private int _rxSitesMin, _rxChains, _rxUnits;
    private double _rxRubberMass, _rxUnitMass;
    private decimal _rxSitesPer;   // 0: all of each chain's sites may react
    private decimal _rxDc = 20, _rxPhrX;
    private string _rxXSmiles = "", _rxXName = "";
    private double _rxXMass;
    private bool _rxCalcBusy;

    /// <summary>"4 chains · 5 reactive sites each (20) · 40 repeat units of 76.32 g/mol".</summary>
    public string RxSitesText { get => _rxSitesText; private set => Set(ref _rxSitesText, value); }
    public bool RxHasSites => _rxChains > 0;
    /// <summary>Sites per chain that may react (0: all); never more than the chain with the fewest has.</summary>
    public decimal RxSitesPerChainD
    {
        get => _rxSitesPer;
        set { if (Set(ref _rxSitesPer, Math.Clamp(Math.Round(value), 0, Math.Max(0, _rxSitesMin)))) RxRecalc(); }
    }
    public int RxSitesMax => _rxSitesMin;
    /// <summary>Degree of crosslinking asked for, %.</summary>
    public decimal RxDcD { get => _rxDc; set { if (Set(ref _rxDc, Math.Clamp(value, 0, 100))) { _rxFromPhr = false; _rxAimDc = true; RxRecalc(); } } }
    // the run aims at the degree shown here whenever it can be counted; choosing another target under Network › Stop at
    // takes over (giving a degree or phr here again hands it back)
    private bool _rxAimDc = true, _rxTargetFromPanel;
    /// <summary>The crosslinker in phr (parts per hundred rubber): sets the degree of crosslinking.</summary>
    public decimal RxPhrXD { get => _rxPhrX; set { if (Set(ref _rxPhrX, Math.Max(0, value))) { _rxFromPhr = true; _rxAimDc = true; RxRecalc(); } } }
    private bool _rxFromPhr;
    /// <summary>The crosslinker molecule (one per link), as SMILES; empty: the reaction joins the chains directly (C–C, peroxide).</summary>
    public string RxCrosslinkerSmiles
    {
        get => _rxXSmiles;
        set { if (Set(ref _rxXSmiles, (value ?? "").Trim())) { Raise(nameof(RxHasCrosslinker)); _ = RxCrosslinkerMassNow(); } }
    }
    public bool RxHasCrosslinker => _rxXSmiles.Length > 0;
    private string _rxCalcText = "";
    public string RxCalcText { get => _rxCalcText; private set => Set(ref _rxCalcText, value); }
    private int _rxLinksNeeded, _rxXNeeded;
    public bool RxCanInsertX => RxHasCrosslinker && _rxXNeeded > 0 && _doc != null && !Busy;
    public string RxInsertXText => _rxXNeeded > 0 ? $"Insert {_rxXNeeded:N0} {(_rxXName.Length > 0 ? _rxXName : "crosslinker")}" : "Insert the crosslinker";
    public const string RxCalcTex = @"\mathrm{DC}=\frac{2\,N_{\mathrm{links}}}{N_{\mathrm{units}}}\times 100\,\%\qquad \mathrm{phr}=\frac{100\,N_{x}\,M_{x}}{m_{\mathrm{rubber}}}";

    /// <summary>The crosslinker a cure brings (MAH for ENR + MAH; maleic acid for ENR + acid), else none.</summary>
    private void RxCrosslinkerDefault()
    {
        (_rxXName, var smi) = _rxSet switch
        {
            6 => ("MAH", "O=C1OC(=O)C=C1"),             // maleic anhydride: one per ENR–MAH–ENR bridge
            5 => ("maleic acid", "OC(=O)/C=C\\C(=O)O"),  // a diacid bridging two epoxides
            _ => ("", ""),
        };
        RxCrosslinkerSmiles = smi;
        Raise(nameof(RxInsertXText));
    }

    /// <summary>The crosslinker's molar mass from its SMILES (built once, hydrogens added).</summary>
    private async Task RxCrosslinkerMassNow()
    {
        _rxXMass = 0;
        var smi = _rxXSmiles;
        if (smi.Length > 0)
            try
            {
                _rxXMass = await Task.Run(() =>
                {
                    var (d, _) = CapsDocument.BuildSmiles(smi, "uff", 1, 1, "crosslinker");
                    using (d) return d.Summary().TotalMass;
                });
                if (smi != _rxXSmiles) return;   // changed meanwhile
            }
            catch (Exception) { _rxXMass = 0; }
        RxRecalc();
    }

    private bool _rxSitesQueued;
    /// <summary>The chains' sites for the reactions in the text, counted again (the text, the structure or a run changed).</summary>
    public void QueueRxSites()
    {
        if (_rxSitesQueued) return;
        _rxSitesQueued = true;
        Avalonia.Threading.Dispatcher.UIThread.Post(() => { _rxSitesQueued = false; RefreshRxSites(); }, Avalonia.Threading.DispatcherPriority.Background);
    }

    public void RefreshRxSites()
    {
        _rxChains = 0;
        if (_doc != null && _rxText.Trim().Length > 0 && !Busy)
            try
            {
                var j = JsonNode.Parse(_doc.ReactSites(_rxText))!;
                var chains = j["chains"]!.AsArray();
                _rxChains = chains.Count;
                _rxUnits = (int)(double)j["units"]!;
                _rxRubberMass = (double)j["mass"]!;
                _rxUnitMass = (double)j["repeat_unit_mass"]!;
                var sites = chains.Select(c => (int)(double)c!["sites"]!).ToList();
                _rxSitesMin = sites.Count > 0 ? sites.Min() : 0;
                var total = sites.Sum();
                var inv = CultureInfo.InvariantCulture;
                var each = sites.Count > 0 && sites.Min() == sites.Max() ? $"{sites[0]:N0} reactive sites each" : $"{sites.Min():N0}–{sites.Max():N0} reactive sites per chain";
                RxSitesText = _rxChains == 0 ? "No chains (molecules of 30 atoms or more) to react"
                    : string.Format(inv, "{0:N0} chains · {1} ({2:N0} in all) · {3:N0} repeat units of {4:0.00} g/mol ({5:N0} g/mol of rubber)",
                                    _rxChains, each, total, _rxUnits, _rxUnitMass, _rxRubberMass);
            }
            catch (Exception e) { RxSitesText = "Reactive sites: " + e.Message; }
        else RxSitesText = _doc == null ? "" : "Choose a reaction to count each chain's reactive sites";
        if (_rxSitesPer > _rxSitesMin) { _rxSitesPer = _rxSitesMin; Raise(nameof(RxSitesPerChainD)); }
        Raise(nameof(RxHasSites));
        Raise(nameof(RxSitesMax));
        RxRecalc();
    }

    /// <summary>Links, crosslinker molecules and phr from the degree of crosslinking (or the degree from the phr), and the run's
    /// target set to it.</summary>
    private void RxRecalc()
    {
        if (_rxCalcBusy) return;
        _rxCalcBusy = true;
        try
        {
            var inv = CultureInfo.InvariantCulture;
            if (_rxChains == 0) { RxCalcText = ""; _rxXNeeded = 0; RxPanelTargetOff(); return; }
            if (_rxUnits == 0)
            {
                RxCalcText = "The chains carry no repeat-unit numbers (residues, as Polymer cell writes them): a degree of crosslinking cannot be counted here — stop at a number of links or links per chain instead (Network › Stop at)";
                _rxXNeeded = 0;
                RxPanelTargetOff();
                return;
            }
            var x = RxHasCrosslinker && _rxXMass > 0;
            if (_rxFromPhr && x)
            {
                // phr → molecules → links (one per molecule) → DC
                var molecules = (double)_rxPhrX * _rxRubberMass / (100 * _rxXMass);
                var dc = Math.Min(100, 200 * molecules / _rxUnits);
                _rxDc = (decimal)Math.Round(dc, 3);
                Raise(nameof(RxDcD));
            }
            var linksExact = (double)_rxDc / 100 * _rxUnits / 2;
            _rxLinksNeeded = (int)Math.Round(linksExact);
            var cap = _rxSitesPer > 0 ? _rxChains * (int)_rxSitesPer / 2 : int.MaxValue;
            _rxXNeeded = x ? _rxLinksNeeded : 0;
            if (x && !_rxFromPhr) { _rxPhrX = (decimal)Math.Round(100 * _rxXNeeded * _rxXMass / _rxRubberMass, 3); Raise(nameof(RxPhrXD)); }
            var lines = new List<string>
            {
                string.Format(inv, "DC {0:0.##} % of {1:N0} repeat units → {2:0.##} links → {3:N0} in this cell", _rxDc, _rxUnits, linksExact, _rxLinksNeeded),
            };
            if (x)
                lines.Add(string.Format(inv, "{0:N0} {1} ({2:0.00} g/mol, one per link) = {3:0.###} phr", _rxXNeeded, _rxXName.Length > 0 ? _rxXName : "crosslinker", _rxXMass,
                                        100 * _rxXNeeded * _rxXMass / _rxRubberMass));
            if (_rxSitesPer > 0)
                lines.Add(_rxLinksNeeded > cap
                    ? string.Format(inv, "⚠ {0} sites per chain allow at most {1:N0} links ({2} chains × {0} ÷ 2): raise it, or ask for less", _rxSitesPer, cap, _rxChains)
                    : string.Format(inv, "{0} of each chain's {1} sites may react (up to {2:N0} links)", _rxSitesPer, _rxSitesMin, cap));
            if (_rxLinksNeeded < 1) lines.Add("⚠ fewer than one link in this cell: grow more chains or ask for more");
            RxCalcText = string.Join("\n", lines);
            // the run aims at it: the degree of crosslinking
            if (_rxAimDc && _rxLinksNeeded >= 1)
            {
                _rxTargetKind = 5; _rxTargetValue = (double)_rxDc; _rxTargetFromPanel = true;
                Raise(nameof(RxTargetKind)); Raise(nameof(RxTargetValueD)); Raise(nameof(RxTargetIsConversion)); Raise(nameof(RxTargetHelp));
            }
            else RxPanelTargetOff();
        }
        finally
        {
            _rxCalcBusy = false;
            Raise(nameof(RxCanInsertX));
            Raise(nameof(RxInsertXText));
        }
    }

    /// <summary>A degree the panel set as the target, taken back when it can no longer be counted (no chains, no repeat units).</summary>
    private void RxPanelTargetOff()
    {
        if (!_rxTargetFromPanel) return;
        _rxTargetFromPanel = false;
        if (_rxTargetKind != 5) return;
        _rxTargetKind = 0; _rxTargetValue = 1;
        Raise(nameof(RxTargetKind)); Raise(nameof(RxTargetValueD)); Raise(nameof(RxTargetIsConversion)); Raise(nameof(RxTargetHelp));
    }
    /// <summary>Another target chosen under Network › Stop at: the panel no longer sets it.</summary>
    internal void RxTargetChosen(int kind) { if (kind != 5) { _rxAimDc = false; _rxTargetFromPanel = false; } }

    /// <summary>The crosslinker molecules the calculation asks for, into the cell's free space.</summary>
    public async Task InsertRxCrosslinker()
    {
        if (_doc == null || !RxCanInsertX) return;
        var doc = _doc;
        var smi = _rxXSmiles;
        var n = _rxXNeeded;
        try
        {
            Status = $"Inserting {n} {_rxXName}…";
            var rep = await Task.Run(() => doc.InsertMolecules(smi, n, 2.0, (ulong)Environment.TickCount64));
            Field.Reset();   // new molecules: the force field is assigned again before the run
            AfterEdit($"{n} {(_rxXName.Length > 0 ? _rxXName : "crosslinker")} inserted · {rep.Split('\n').FirstOrDefault()}");
            RefreshRxSites();
        }
        catch (Exception e) { Status = "Could not insert the crosslinker: " + e.Message; }
    }

    // ---- the force field the reaction runs with (assigned before the run, re-typed after every cycle)
    private int _rxFf = -1;
    /// <summary>The Field library's force fields, first "the assigned one".</summary>
    public List<string> RxForceFields => ["The assigned force field", .. Field.Library.Select(e => e.Name)];
    /// <summary>0: the assigned force field; k: the library's k-1-th, assigned before the run.</summary>
    public int RxFfChoice { get => _rxFf + 1; set { if (Set(ref _rxFf, Math.Max(-1, value - 1))) Raise(nameof(RxFfNote)); } }
    public string RxFfNote => _rxFf < 0
        ? "The assigned force field types the product after every cycle: new bonds change atom types (CH₂ → CH, epoxide O → ester), their charges and parameters follow"
        : $"{Field.Library[_rxFf].Name} is assigned before the run and types the product after every cycle (the new types' parameters from it; a missing one stops the run, the last good cycle kept)";

    /// <summary>The reaction's force field assigned (when another than the assigned one is chosen, or none is assigned yet).</summary>
    private async Task<bool> RxAssignField()
    {
        if (!_rxUseField) return true;
        if (_rxFf >= 0 && (!Field.Assigned || Field.FfIndex != _rxFf))
        {
            Field.FfIndex = _rxFf;
            await Field.Assign();
        }
        if (!Field.Assigned) return true;   // none assigned and none chosen: the built-in default types the run, as before
        if (!Field.Complete)
        {
            RxLog = $"{Field.ForceFieldName} does not cover this structure (untyped atoms or missing terms): see the Force field step, or choose another force field for the reaction";
            Status = "React: the force field is incomplete";
            return false;
        }
        return true;
    }
}
