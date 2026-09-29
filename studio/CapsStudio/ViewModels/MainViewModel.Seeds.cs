using System.Text.RegularExpressions;

namespace CapsStudio.ViewModels;

/// <summary>The random seed of every kind of run (SeedPicker on each page): fixed to reproduce a run, or new each time
/// with the seed used shown. Each run takes its seed once, when it starts (a queued run when it is queued).</summary>
public sealed partial class MainViewModel
{
    private SeedChoice? _polySeedC, _growSeedC, _mdSeedC, _packSeedC, _rxSeedC, _cgSeedC, _mtSeedC, _adsSeedC, _sorbSeedC,
        _dpdSeedC, _blendSeedC, _nanoSeedC, _filmSeedC, _solvSeedC;
    private int _polyDraw = 1;   // the seed of the chain shown on the Polymer builder (strip, 3D preview, Build in document)

    /// <summary>Polymer builder: the sequence and the chain drawn; new each time redraws on every change.</summary>
    public SeedChoice PolySeed => _polySeedC ??= new SeedChoice(1, false, () => PolyChanged());
    public SeedChoice GrowSeedChoice => _growSeedC ??= new SeedChoice(System.Math.Max(1, _growSeed), false, () => { _growSeed = (int)(_growSeedC!.Value ?? 1); Raise(nameof(GrowSeedD)); });
    /// <summary>Dynamics, Equilibrate and CBMC: velocities and thermostat noise.</summary>
    public SeedChoice MdSeedChoice => _mdSeedC ??= new SeedChoice(System.Math.Max(1, _mdSeed), false, () => { _mdSeed = (int)(_mdSeedC!.Value ?? 1); Raise(nameof(MdSeedD)); });
    /// <summary>Pack: the input's seed line follows it.</summary>
    public SeedChoice PackSeedChoice => _packSeedC ??= new SeedChoice(System.Math.Max(1, _packSeed), false, () => { _packSeed = (int)(_packSeedC!.Value ?? 1); SetPackSeedLine(_packSeed); });
    public SeedChoice RxSeedChoice => _rxSeedC ??= new SeedChoice(System.Math.Max(1, _rxSeed), false, () => { _rxSeed = (int)(_rxSeedC!.Value ?? 1); Raise(nameof(RxSeedD)); });
    public SeedChoice CgSeedChoice => _cgSeedC ??= new SeedChoice((int)_cgSeed, false, () => { _cgSeed = _cgSeedC!.Value ?? 1; CgPreview(); });
    public SeedChoice MtSeedChoice => _mtSeedC ??= new SeedChoice((int)_mtSeed, false, () => _mtSeed = _mtSeedC!.Value ?? 1);
    public SeedChoice AdsSeed => _adsSeedC ??= new SeedChoice();
    public SeedChoice SorbSeed => _sorbSeedC ??= new SeedChoice();
    public SeedChoice DpdSeed => _dpdSeedC ??= new SeedChoice();
    public SeedChoice BlendSeed => _blendSeedC ??= new SeedChoice();
    /// <summary>Nanostructures: where groups are grafted, and the polymer matrix grown around them.</summary>
    public SeedChoice NanoSeed => _nanoSeedC ??= new SeedChoice();
    public SeedChoice FilmSeed => _filmSeedC ??= new SeedChoice();
    public SeedChoice SolvSeed => _solvSeedC ??= new SeedChoice();

    // the Pack input's seed line (added when the input has none)
    private static readonly Regex PackSeedLine = new(@"^[ \t]*seed[ \t]+\S+[ \t]*$", RegexOptions.Multiline | RegexOptions.IgnoreCase);
    private static string WithSeedLine(string text, int seed) =>
        PackSeedLine.IsMatch(text) ? PackSeedLine.Replace(text, $"seed {seed}", 1) : $"seed {seed}\n" + text;
    private void SetPackSeedLine(int seed)
    {
        if (_packText.Trim().Length > 0) PackText = WithSeedLine(_packText, seed);
    }

    /// <summary>Polymer builder: another draw of the sequence and chain (new each run).</summary>
    public void PolyRedraw() => PolyChanged();
}
