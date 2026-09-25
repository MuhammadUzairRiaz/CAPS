using System.Runtime.InteropServices;

namespace CapsStudio.Interop;

// Mirrors capi/include/caps_c.h (ABI v7). Keep field order and types identical.

[StructLayout(LayoutKind.Sequential)]
public struct CapsCamera
{
    public double Yaw, Pitch, Zoom, PanX, PanY;
    public int Perspective;
}

[StructLayout(LayoutKind.Sequential)]
public struct CapsRenderOpts
{
    public int Width, Height, Supersample;
    public int Background;      // 0 dark, 1 white, 2 transparent, 3 custom
    public uint CustomRgb;
    public int ColourBy;        // 0 element, 1 molecule, 2 type, 3 distance to own COM
    public int Style;           // 0 ball & stick, 1 space filling, 2 sticks, 3 no H, 4 backbone
    public int Outlines, DepthCue, ShowCell;
    public int Highlight0, Highlight1, Highlight2, Highlight3;   // selected atoms, -1 unused
    public int Focus;           // atom index + 1 with the keyboard-focus ring, 0 none (ABI 18)
    public int AmbientOcclusion;   // darken atoms by the open sky they see (ABI 19)
}

[StructLayout(LayoutKind.Sequential, CharSet = CharSet.Ansi)]
public struct CapsAtomInfo
{
    public long Id, Mol;
    public int Type, Element, Index;
    public double Charge, X, Y, Z;
    [MarshalAs(UnmanagedType.ByValTStr, SizeConst = 4)] public string ElementSymbol;
    [MarshalAs(UnmanagedType.ByValTStr, SizeConst = 16)] public string Name;
}

[StructLayout(LayoutKind.Sequential)]
public struct CapsMolecule
{
    public int Molecule, Atoms;
    public double Mass, Rg, Kappa2, ComX, ComY, ComZ;
}

[StructLayout(LayoutKind.Sequential, CharSet = CharSet.Ansi)]
public struct CapsSummary
{
    public long Atoms, Bonds, Molecules, Frames;
    public int BondsFromFile, HasCharges, CellValid, Unwrapped;
    public double CellA, CellB, CellC, Volume, Density, TotalMass, TotalCharge;
    [MarshalAs(UnmanagedType.ByValTStr, SizeConst = 24)] public string Format;
}

[StructLayout(LayoutKind.Sequential)]
public struct CapsBuildOpts
{
    public int Conformers;
    public ulong Seed;
}

[StructLayout(LayoutKind.Sequential)]
public struct CapsGrowOpts
{
    public int Chains, Dp, Tacticity;
    public ulong Seed;
    public double Box, Density, ContactScale;
    public int Curve;
}

[UnmanagedFunctionPointer(CallingConvention.Cdecl)]
public delegate int CapsProgress(int done, int total, int restarts, IntPtr user);

[StructLayout(LayoutKind.Sequential)]
public struct CapsRelaxOpts
{
    public int Method;            // 0 steepest descent, 1 CG, 2 L-BFGS, 3 FIRE
    public double Ftol;
    public int MaxIterations;
    public double TargetDensity, CompressStep;
    public int Pushoff, RelaxBox;
    public double Pressure, Cutoff;
    public int Coulomb, Threads;
}

[StructLayout(LayoutKind.Sequential)]
public struct CapsMdOpts
{
    public double Dt;
    public long Steps;
    public double Temperature;
    public int Thermostat;        // 0 none, 1 Bussi, 2 Langevin
    public double TauT;
    public int Barostat;          // 0 none, 1 C-rescale, 2 Berendsen
    public double Pressure, TauP;
    public int NewVelocities;
    public ulong Seed;
    public int ThermoEvery, FrameEvery;
    public double Cutoff;
    public int Coulomb, Tail, Threads;
}

[StructLayout(LayoutKind.Sequential)]
public struct CapsThermo
{
    public long Step;
    public double TimePs, Temperature, Potential, Kinetic, Total, Conserved, Pressure, Volume, Density;
}

[StructLayout(LayoutKind.Sequential)]
public struct CapsProtocolParams
{
    public double TFinal, TMax, PFinal, PMax, TimeScale;
    public int Cycles;
    public double TLow, THigh, RampPs, HoldPs;
}

[StructLayout(LayoutKind.Sequential)]
public struct CapsEquilOpts
{
    public double Dt;
    public int Thermostat, Barostat;
    public double TauT, TauP;
    public ulong Seed;
    public double Cutoff;
    public int Coulomb, Tail, Threads;
    public double FramePs, ThermoPs;
    public int UntilConverged;
    public double BlockPs;
    public int MaxBlocks;
    public double TolDensity, TolEnergy, TolRg;
}

[UnmanagedFunctionPointer(CallingConvention.Cdecl)]
public delegate int CapsEquilProgress(int stage, int stages, IntPtr label, in CapsThermo row, IntPtr user);

[StructLayout(LayoutKind.Sequential)]
public struct CapsReactOpts
{
    public ulong Seed;
    public int MaxCycles, MaxPerCycle;
    public double TargetConversion, Capture;
    public int Relax, RelaxIterations;
    public double MdPs, Temperature, Cutoff;
    public int Coulomb;
}

[StructLayout(LayoutKind.Sequential)]
public struct CapsReactCycle
{
    public int Cycle, Reactions, Total, Clusters, Atoms;
    public double Conversion, LargestFraction, ReducedMw, Energy;
}

[StructLayout(LayoutKind.Sequential)]
public struct CapsAnalyzeOpts
{
    public long First, Last, Stride;
    public double FramePs, TimestepFs;
    public int Blocks;
    public int ElemA, ElemB, InterOnly;
    public double RdfRmax, RdfDr, Qmax, Dq, QDirect;
    public double FitFrom, FitTo;
    public double Probe, Grid;
    public double Cutoff;
    public int Threads;
}

[StructLayout(LayoutKind.Sequential)]
public struct CapsMechOpts
{
    public int Configurations;
    public double Strain;
    public double Temperature;
    public int Axis;
    public double Rate, MaxStrain, FitStrain;
    public int LateralFixed;
    public double TStart, TEnd, TStep, PsPerStep;
    public double Dt, Pressure;
    public ulong Seed;
    public double RunPs;
    public double EquilibratePs;
}

[UnmanagedFunctionPointer(CallingConvention.Cdecl)]
public delegate int CapsOpenProgress(int stage, double fraction, IntPtr detail, IntPtr user);

[UnmanagedFunctionPointer(CallingConvention.Cdecl)]
public delegate int CapsAnalyzeProgress(IntPtr what, double fraction, IntPtr user);

[UnmanagedFunctionPointer(CallingConvention.Cdecl)]
public delegate int CapsReactProgress(in CapsReactCycle row, IntPtr user);

[UnmanagedFunctionPointer(CallingConvention.Cdecl)]
public delegate int CapsPackProgress(int loop, int loops, double penalty, int bad, IntPtr user);
[UnmanagedFunctionPointer(CallingConvention.Cdecl)]
public delegate int CapsSeriesProgress(int done, int total, IntPtr user);
[UnmanagedFunctionPointer(CallingConvention.Cdecl)]
public delegate int CapsStageProgress(int stage, int loop, int loops, double dmin, int bad, IntPtr user);

[UnmanagedFunctionPointer(CallingConvention.Cdecl)]
public delegate int CapsMdProgress(in CapsThermo row, long steps, IntPtr user);

[UnmanagedFunctionPointer(CallingConvention.Cdecl)]
public delegate int CapsRelaxProgress(int stage, int stages, int iteration, double energy, double fmax, double density, IntPtr user);

internal static class Native
{
    private const string Lib = "caps";

    [DllImport(Lib, EntryPoint = "caps_abi_version")] public static extern int AbiVersion();
    [DllImport(Lib, EntryPoint = "caps_last_error")] private static extern IntPtr LastErrorPtr();
    public static string LastError() => Marshal.PtrToStringUTF8(LastErrorPtr()) ?? "";

    [DllImport(Lib, EntryPoint = "caps_inspect_file")] public static extern int InspectFile([MarshalAs(UnmanagedType.LPUTF8Str)] string path, [MarshalAs(UnmanagedType.LPUTF8Str)] string? topology, byte[]? json, int cap);
    [DllImport(Lib, EntryPoint = "caps_open_staged")] public static extern IntPtr OpenStaged([MarshalAs(UnmanagedType.LPUTF8Str)] string path, [MarshalAs(UnmanagedType.LPUTF8Str)] string? topology, int maxFrames, CapsOpenProgress? progress, IntPtr user);
    [DllImport(Lib, EntryPoint = "caps_adopt_frames")] public static extern int AdoptFrames(IntPtr dst, IntPtr src);
    [DllImport(Lib, EntryPoint = "caps_import")] public static extern IntPtr Import([MarshalAs(UnmanagedType.LPUTF8Str)] string path, [MarshalAs(UnmanagedType.LPUTF8Str)] string? topology, [MarshalAs(UnmanagedType.LPUTF8Str)] string options);
    [DllImport(Lib, EntryPoint = "caps_import_preview")] public static extern int ImportPreview([MarshalAs(UnmanagedType.LPUTF8Str)] string path, [MarshalAs(UnmanagedType.LPUTF8Str)] string options, byte[]? json, int cap);
    [DllImport(Lib, EntryPoint = "caps_import_fragment")] public static extern IntPtr ImportFragment([MarshalAs(UnmanagedType.LPUTF8Str)] string path, [MarshalAs(UnmanagedType.LPUTF8Str)] string options);
    [DllImport(Lib, EntryPoint = "caps_export_image")] public static extern int ExportImage(IntPtr doc, in CapsCamera cam, in CapsRenderOpts opt, [MarshalAs(UnmanagedType.LPUTF8Str)] string path, [MarshalAs(UnmanagedType.LPUTF8Str)] string options, byte[]? overlay);
    [DllImport(Lib, EntryPoint = "caps_png_text")] public static extern int PngText([MarshalAs(UnmanagedType.LPUTF8Str)] string path, byte[]? json, int cap);
    [DllImport(Lib, EntryPoint = "caps_export_movie")] public static extern int ExportMovie(IntPtr doc, in CapsCamera cam, in CapsRenderOpts opt, [MarshalAs(UnmanagedType.LPUTF8Str)] string path, [MarshalAs(UnmanagedType.LPUTF8Str)] string options, CapsSeriesProgress? progress, IntPtr user);
    [DllImport(Lib, EntryPoint = "caps_provenance")] public static extern int Provenance(IntPtr doc, byte[]? json, int cap);
    [DllImport(Lib, EntryPoint = "caps_provenance_file")] public static extern int ProvenanceFile([MarshalAs(UnmanagedType.LPUTF8Str)] string path, byte[]? json, int cap);
    [DllImport(Lib, EntryPoint = "caps_provenance_compare")] public static extern int ProvenanceCompare([MarshalAs(UnmanagedType.LPUTF8Str)] string a, [MarshalAs(UnmanagedType.LPUTF8Str)] string b, byte[]? json, int cap);
    [DllImport(Lib, EntryPoint = "caps_provenance_bibtex")] public static extern int ProvenanceBibtex([MarshalAs(UnmanagedType.LPUTF8Str)] string manifest, byte[]? text, int cap);
    [DllImport(Lib, EntryPoint = "caps_open")] public static extern IntPtr Open([MarshalAs(UnmanagedType.LPUTF8Str)] string path, [MarshalAs(UnmanagedType.LPUTF8Str)] string? topology);
    [DllImport(Lib, EntryPoint = "caps_grow")] public static extern IntPtr Grow(in CapsGrowOpts o, CapsProgress? progress, IntPtr user, byte[] report, int cap);
    [DllImport(Lib, EntryPoint = "caps_relax")] public static extern int Relax(IntPtr doc, in CapsRelaxOpts o, CapsRelaxProgress? progress, IntPtr user, byte[] report, int cap);
    [DllImport(Lib, EntryPoint = "caps_md")] public static extern int Md(IntPtr doc, in CapsMdOpts o, CapsMdProgress? progress, IntPtr user, byte[] report, int cap);
    [DllImport(Lib, EntryPoint = "caps_save_trajectory")] public static extern int SaveTrajectory(IntPtr doc, [MarshalAs(UnmanagedType.LPUTF8Str)] string path);
    [DllImport(Lib, EntryPoint = "caps_protocol_text")] public static extern int ProtocolText([MarshalAs(UnmanagedType.LPUTF8Str)] string name, in CapsProtocolParams p, byte[] text, int cap);
    [DllImport(Lib, EntryPoint = "caps_equilibrate")] public static extern int Equilibrate(IntPtr doc, byte[] protocol, in CapsEquilOpts o, CapsEquilProgress? progress, IntPtr user, byte[] report, int cap);
    [DllImport(Lib, EntryPoint = "caps_internal_distances")] public static extern int InternalDistances(IntPtr doc, [Out] int[] n, [Out] double[] ratio, int cap, out int chains, out double b2);
    [UnmanagedFunctionPointer(CallingConvention.Cdecl)]
    public delegate int BenchProgressFn([MarshalAs(UnmanagedType.LPUTF8Str)] string table, [MarshalAs(UnmanagedType.LPUTF8Str)] string what, double fraction, IntPtr user);
    [DllImport(Lib, EntryPoint = "caps_bench_list")] public static extern int BenchList(byte[]? json, int cap);
    [DllImport(Lib, EntryPoint = "caps_bench_run")] private static extern int BenchRun([MarshalAs(UnmanagedType.LPUTF8Str)] string id, [MarshalAs(UnmanagedType.LPUTF8Str)] string samples,
        [MarshalAs(UnmanagedType.LPUTF8Str)] string? forcefields, int repeats, int quick, BenchProgressFn? progress, IntPtr user, byte[]? json, int cap);
    [DllImport(Lib, EntryPoint = "caps_bench_write")] public static extern int BenchWrite([MarshalAs(UnmanagedType.LPUTF8Str)] string tables, [MarshalAs(UnmanagedType.LPUTF8Str)] string dir);

    /// <summary>Runs one bench table (the progress callback runs on the calling thread; false cancels) and returns its JSON.</summary>
    public static string BenchRunJson(string id, string samples, string? forcefields, int repeats, bool quick, Func<string, string, double, bool>? progress)
    {
        BenchProgressFn? cb = progress == null ? null : (t, w, f, _) => progress(t, w, f) ? 0 : 1;
        var buf = new byte[1 << 20];
        var n = BenchRun(id, samples, forcefields, repeats, quick ? 1 : 0, cb, IntPtr.Zero, buf, buf.Length);
        GC.KeepAlive(cb);
        if (n < 0) throw new InvalidOperationException(LastError());
        return System.Text.Encoding.UTF8.GetString(buf, 0, Math.Min(n, buf.Length) - 1);
    }

    [DllImport(Lib, EntryPoint = "caps_unit_info")] public static extern int UnitInfo([MarshalAs(UnmanagedType.LPUTF8Str)] string smiles, byte[]? json, int cap);
    [DllImport(Lib, EntryPoint = "caps_chain_preview")] public static extern int ChainPreview([MarshalAs(UnmanagedType.LPUTF8Str)] string spec, ulong seed, byte[]? json, int cap);
    [DllImport(Lib, EntryPoint = "caps_grow_chains")] public static extern IntPtr GrowChains([MarshalAs(UnmanagedType.LPUTF8Str)] string spec, in CapsGrowOpts o, CapsProgress? progress, IntPtr user, byte[] report, int cap);
    [DllImport(Lib, EntryPoint = "caps_surface_terminations")] public static extern int SurfaceTerminations([MarshalAs(UnmanagedType.LPUTF8Str)] string cif, int h, int k, int l, byte[]? json, int cap);
    [DllImport(Lib, EntryPoint = "caps_surface_build")] public static extern IntPtr SurfaceBuild([MarshalAs(UnmanagedType.LPUTF8Str)] string cif, [MarshalAs(UnmanagedType.LPUTF8Str)] string options, byte[] report, int cap);
    [DllImport(Lib, EntryPoint = "caps_interface_build")] public static extern IntPtr InterfaceBuild([MarshalAs(UnmanagedType.LPUTF8Str)] string options, [MarshalAs(UnmanagedType.LPUTF8Str)] string spec, in CapsGrowOpts o, CapsProgress? progress, IntPtr user, byte[] report, int cap);
    [DllImport(Lib, EntryPoint = "caps_nano_build")] public static extern IntPtr NanoBuild([MarshalAs(UnmanagedType.LPUTF8Str)] string options, byte[] report, int cap);
    [DllImport(Lib, EntryPoint = "caps_nano_embed")] public static extern IntPtr NanoEmbed([MarshalAs(UnmanagedType.LPUTF8Str)] string options, [MarshalAs(UnmanagedType.LPUTF8Str)] string spec, in CapsGrowOpts o, CapsProgress? progress, IntPtr user, byte[] report, int cap);
    [DllImport(Lib, EntryPoint = "caps_grow_blend")] public static extern IntPtr GrowBlend([MarshalAs(UnmanagedType.LPUTF8Str)] string options, in CapsGrowOpts o, CapsProgress? progress, IntPtr user, byte[] report, int cap);
    [DllImport(Lib, EntryPoint = "caps_file_checks")] public static extern int FileChecks(IntPtr doc, byte[]? json, int cap);
    [DllImport(Lib, EntryPoint = "caps_insert_molecules")] public static extern int InsertMolecules(IntPtr doc, [MarshalAs(UnmanagedType.LPUTF8Str)] string smiles, int count, double tolerance, ulong seed, byte[] report, int cap);
    [DllImport(Lib, EntryPoint = "caps_set_held_molecule")] public static extern void SetHeldMolecule(IntPtr doc, long mol);
    [DllImport(Lib, EntryPoint = "caps_held_molecule")] public static extern long HeldMolecule(IntPtr doc);
    [DllImport(Lib, EntryPoint = "caps_peptide_info")] public static extern int PeptideInfo([MarshalAs(UnmanagedType.LPUTF8Str)] string options, byte[]? json, int cap);
    [DllImport(Lib, EntryPoint = "caps_peptide_build")] public static extern IntPtr PeptideBuild([MarshalAs(UnmanagedType.LPUTF8Str)] string options, byte[] report, int cap);
    [DllImport(Lib, EntryPoint = "caps_fasta_sequence")] public static extern int FastaSequence([MarshalAs(UnmanagedType.LPUTF8Str)] string text, byte[]? seq, int cap);
    [DllImport(Lib, EntryPoint = "caps_element_number")] public static extern int ElementNumber([MarshalAs(UnmanagedType.LPUTF8Str)] string symbol);
    [DllImport(Lib, EntryPoint = "caps_element_info")] private static extern int ElementInfo(int z, out double mass, out double covalent, out double vdw, out uint rgb);
    public static double ElementMass(int z) => ElementInfo(z, out var m, out _, out _, out _) == 0 ? m : 0;
    public static double ElementCovalent(int z) => ElementInfo(z, out _, out var c, out _, out _) == 0 ? c : 0;
    public static double ElementVdw(int z) => ElementInfo(z, out _, out _, out var v, out _) == 0 ? v : 0;
    public static uint ElementColour(int z) => ElementInfo(z, out _, out _, out _, out var c) == 0 ? c : 0x9AA1A8;
    [DllImport(Lib, EntryPoint = "caps_interactions")] public static extern int Interactions(IntPtr doc, [MarshalAs(UnmanagedType.LPUTF8Str)] string options, byte[]? json, int cap);
    [DllImport(Lib, EntryPoint = "caps_clear_checks")] public static extern void ClearChecks(IntPtr doc);
    [DllImport(Lib, EntryPoint = "caps_edit")] public static extern int Edit(IntPtr doc, [MarshalAs(UnmanagedType.LPUTF8Str)] string json, byte[]? outp, int cap);
    [DllImport(Lib, EntryPoint = "caps_undo")] public static extern int Undo(IntPtr doc, int redo);
    [DllImport(Lib, EntryPoint = "caps_history")] public static extern int History(IntPtr doc, byte[]? json, int cap);
    [DllImport(Lib, EntryPoint = "caps_select")] public static extern int Select(IntPtr doc, [MarshalAs(UnmanagedType.LPUTF8Str)] string json, byte[]? outp, int cap);
    [DllImport(Lib, EntryPoint = "caps_selection")] public static extern int Selection(IntPtr doc, byte[]? json, int cap);
    [DllImport(Lib, EntryPoint = "caps_tacticity")] public static extern int Tacticity(IntPtr doc, byte[]? json, int cap);
    [DllImport(Lib, EntryPoint = "caps_torsion_scan")] public static extern int TorsionScan(IntPtr doc, [MarshalAs(UnmanagedType.LPUTF8Str)] string options, CapsSeriesProgress? progress, IntPtr user, byte[]? json, int cap);
    [DllImport(Lib, EntryPoint = "caps_torsion_show")] public static extern int TorsionShow(IntPtr doc, int index);
    [DllImport(Lib, EntryPoint = "caps_default_torsion")] public static extern int DefaultTorsion(IntPtr doc, int[] atoms);
    [DllImport(Lib, EntryPoint = "caps_trajectory_series")] public static extern int TrajectorySeries(IntPtr doc, [MarshalAs(UnmanagedType.LPUTF8Str)] string options, IntPtr progress, IntPtr user, byte[]? json, int cap);
    [DllImport(Lib, EntryPoint = "caps_set_smoothing")] public static extern void SetSmoothing(IntPtr doc, int window);
    [DllImport(Lib, EntryPoint = "caps_set_appearance")] public static extern int SetAppearance(IntPtr doc, [MarshalAs(UnmanagedType.LPUTF8Str)] string json);
    [DllImport(Lib, EntryPoint = "caps_appearance_info")] public static extern int AppearanceInfo(IntPtr doc, byte[]? json, int cap);
    [DllImport(Lib, EntryPoint = "caps_atom_labels")] public static extern int AtomLabels(IntPtr doc, [MarshalAs(UnmanagedType.LPUTF8Str)] string kind, byte[]? json, int cap);
    [DllImport(Lib, EntryPoint = "caps_project_atoms")] public static extern int ProjectAtoms(IntPtr doc, in CapsCamera cam, in CapsRenderOpts opt, float[] xyv, int count);
    [DllImport(Lib, EntryPoint = "caps_solvent_library")] public static extern int SolventLibrary(byte[]? json, int cap);
    [DllImport(Lib, EntryPoint = "caps_solvate_plan")] public static extern int SolvatePlan(IntPtr solute, [MarshalAs(UnmanagedType.LPUTF8Str)] string options, byte[]? json, int cap);
    [DllImport(Lib, EntryPoint = "caps_solvate")] public static extern IntPtr Solvate(IntPtr solute, [MarshalAs(UnmanagedType.LPUTF8Str)] string options, CapsStageProgress? progress, IntPtr user, byte[] report, int cap);
    [DllImport(Lib, EntryPoint = "caps_space_groups")] public static extern int SpaceGroups(byte[]? json, int cap);
    [DllImport(Lib, EntryPoint = "caps_crystal_info")] public static extern int CrystalInfo([MarshalAs(UnmanagedType.LPUTF8Str)] string spec, byte[]? json, int cap);
    [DllImport(Lib, EntryPoint = "caps_crystal_build")] public static extern IntPtr CrystalBuild([MarshalAs(UnmanagedType.LPUTF8Str)] string spec, byte[] report, int cap);
    [DllImport(Lib, EntryPoint = "caps_crystal_symmetrize")] public static extern int CrystalSymmetrize([MarshalAs(UnmanagedType.LPUTF8Str)] string spec, double snap, byte[]? json, int cap);
    [DllImport(Lib, EntryPoint = "caps_crystal_find_symmetry")] public static extern int CrystalFindSymmetry([MarshalAs(UnmanagedType.LPUTF8Str)] string? spec, [MarshalAs(UnmanagedType.LPUTF8Str)] string? cif, double tolerance, byte[]? json, int cap);
    [DllImport(Lib, EntryPoint = "caps_set_palette")] public static extern void SetPalette(int palette);
    [DllImport(Lib, EntryPoint = "caps_set_threads")] public static extern void SetThreads(int threads);
    [DllImport(Lib, EntryPoint = "caps_set_electrostatics")] public static extern void SetElectrostatics(int mode, double ewaldRtol, double pmeSpacing, int pmeOrder);
    [DllImport(Lib, EntryPoint = "caps_smiles_info")] public static extern int SmilesInfo([MarshalAs(UnmanagedType.LPUTF8Str)] string smiles, byte[]? json, int cap);
    [DllImport(Lib, EntryPoint = "caps_smiles_depict")] public static extern int SmilesDepict([MarshalAs(UnmanagedType.LPUTF8Str)] string smiles, byte[]? json, int cap);
    [DllImport(Lib, EntryPoint = "caps_smiles_write")] public static extern int SmilesWrite([MarshalAs(UnmanagedType.LPUTF8Str)] string graph, byte[]? smiles, int cap);
    [DllImport(Lib, EntryPoint = "caps_build_smiles")] public static extern IntPtr BuildSmiles([MarshalAs(UnmanagedType.LPUTF8Str)] string smiles, [MarshalAs(UnmanagedType.LPUTF8Str)] string? ff, in CapsBuildOpts o, byte[] report, int cap);
    [DllImport(Lib, EntryPoint = "caps_pack")] public static extern IntPtr Pack(byte[] text, [MarshalAs(UnmanagedType.LPUTF8Str)] string baseDir, int threads, CapsPackProgress? progress, IntPtr user, byte[] report, int cap);
    [DllImport(Lib, EntryPoint = "caps_reaction_template")] public static extern int ReactionTemplate([MarshalAs(UnmanagedType.LPUTF8Str)] string name, byte[] text, int cap);
    [DllImport(Lib, EntryPoint = "caps_react")] public static extern int React(IntPtr doc, byte[] templates, in CapsReactOpts o, CapsReactProgress? progress, IntPtr user, byte[] report, int cap);
    [DllImport(Lib, EntryPoint = "caps_field_info")] public static extern int FieldInfo(IntPtr doc, byte[] text, int cap);
    [DllImport(Lib, EntryPoint = "caps_field_assign")] public static extern int FieldAssign(IntPtr doc, [MarshalAs(UnmanagedType.LPUTF8Str)] string ff, [MarshalAs(UnmanagedType.LPUTF8Str)] string? rules, int charges);
    [DllImport(Lib, EntryPoint = "caps_analyze")] public static extern int Analyze(IntPtr doc, [MarshalAs(UnmanagedType.LPUTF8Str)] string props, in CapsAnalyzeOpts o, CapsAnalyzeProgress? progress, IntPtr user);
    [DllImport(Lib, EntryPoint = "caps_analyze_ex")] public static extern int AnalyzeEx(IntPtr doc, [MarshalAs(UnmanagedType.LPUTF8Str)] string props, in CapsAnalyzeOpts o, in CapsMechOpts m, CapsAnalyzeProgress? progress, IntPtr user);
    [DllImport(Lib, EntryPoint = "caps_lammps_input")] public static extern int LammpsInput(IntPtr doc, [MarshalAs(UnmanagedType.LPUTF8Str)] string dataName, byte[]? text, int cap);
    [DllImport(Lib, EntryPoint = "caps_equilibrate_checks")] public static extern int EquilibrateChecks(IntPtr doc, byte[]? json, int cap);
    [DllImport(Lib, EntryPoint = "caps_analyze_report")] public static extern int AnalyzeReport(IntPtr doc, byte[]? json, int cap);
    [DllImport(Lib, EntryPoint = "caps_field_report")] public static extern int FieldReport(IntPtr doc, byte[]? json, int cap);
    [DllImport(Lib, EntryPoint = "caps_field_override")] public static extern int FieldOverride(IntPtr doc, int index, [MarshalAs(UnmanagedType.LPUTF8Str)] string? type);
    [DllImport(Lib, EntryPoint = "caps_field_add_rule")] public static extern int FieldAddRule(IntPtr doc, [MarshalAs(UnmanagedType.LPUTF8Str)] string kind, [MarshalAs(UnmanagedType.LPUTF8Str)] string types,
        [MarshalAs(UnmanagedType.LPUTF8Str)] string style, [MarshalAs(UnmanagedType.LPUTF8Str)] string pars);
    [DllImport(Lib, EntryPoint = "caps_field_import")] public static extern int FieldImport(IntPtr doc, [MarshalAs(UnmanagedType.LPUTF8Str)] string path);
    [DllImport(Lib, EntryPoint = "caps_field_remove_rules")] public static extern int FieldRemoveRules(IntPtr doc);
    [DllImport(Lib, EntryPoint = "caps_field_clear")] public static extern int FieldClear(IntPtr doc);
    [DllImport(Lib, EntryPoint = "caps_field_types_file")] public static extern int FieldTypesFile(IntPtr doc, [MarshalAs(UnmanagedType.LPUTF8Str)] string path);
    [DllImport(Lib, EntryPoint = "caps_save")] public static extern int Save(IntPtr doc, [MarshalAs(UnmanagedType.LPUTF8Str)] string path);
    [DllImport(Lib, EntryPoint = "caps_close")] public static extern void Close(IntPtr doc);
    [DllImport(Lib, EntryPoint = "caps_summary_get")] public static extern int Summary(IntPtr doc, out CapsSummary s);
    [DllImport(Lib, EntryPoint = "caps_set_frame")] public static extern int SetFrame(IntPtr doc, long frame);
    [DllImport(Lib, EntryPoint = "caps_set_wrap")] public static extern int SetWrap(IntPtr doc, int wrap);
    [DllImport(Lib, EntryPoint = "caps_atom")] public static extern int Atom(IntPtr doc, int index, out CapsAtomInfo a);
    [DllImport(Lib, EntryPoint = "caps_note_count")] public static extern int NoteCount(IntPtr doc);
    [DllImport(Lib, EntryPoint = "caps_note")] private static extern IntPtr NotePtr(IntPtr doc, int k);
    public static string Note(IntPtr doc, int k) => Marshal.PtrToStringUTF8(NotePtr(doc, k)) ?? "";

    [DllImport(Lib, EntryPoint = "caps_render")] public static extern unsafe int Render(IntPtr doc, in CapsCamera cam, in CapsRenderOpts opt, byte* rgba);
    [DllImport(Lib, EntryPoint = "caps_pick")] public static extern int Pick(IntPtr doc, int x, int y);
    [DllImport(Lib, EntryPoint = "caps_export_png")] public static extern int ExportPng(IntPtr doc, in CapsCamera cam, in CapsRenderOpts opt, [MarshalAs(UnmanagedType.LPUTF8Str)] string path);
    [DllImport(Lib, EntryPoint = "caps_export_svg")] public static extern int ExportSvg(IntPtr doc, in CapsCamera cam, in CapsRenderOpts opt, [MarshalAs(UnmanagedType.LPUTF8Str)] string path);
    [DllImport(Lib, EntryPoint = "caps_measure")] public static extern int Measure(IntPtr doc, int[] idx, int n, out double value);
    [DllImport(Lib, EntryPoint = "caps_rdf")] public static extern int Rdf(IntPtr doc, int ea, int eb, double rmax, double dr, int inter, [Out] double[] r, [Out] double[] g, int cap);
    [DllImport(Lib, EntryPoint = "caps_molecules")] public static extern int Molecules(IntPtr doc, [Out] CapsMolecule[] out_, int cap);
    [DllImport(Lib, EntryPoint = "caps_property_range")] public static extern int PropertyRange(IntPtr doc, out double lo, out double hi);
    [DllImport(Lib, EntryPoint = "caps_pipeline_set")] public static extern int PipelineSet(IntPtr doc, [MarshalAs(UnmanagedType.LPUTF8Str)] string? json);
    [DllImport(Lib, EntryPoint = "caps_pipeline_result")] public static extern int PipelineResult(IntPtr doc, byte[]? json, int cap);
    [DllImport(Lib, EntryPoint = "caps_pipeline_particles")] public static extern int PipelineParticles(IntPtr doc, [MarshalAs(UnmanagedType.LPUTF8Str)] string? filter, int offset, int count, byte[]? json, int cap);
    [DllImport(Lib, EntryPoint = "caps_pipeline_bonds")] public static extern int PipelineBonds(IntPtr doc, int offset, int count, byte[]? json, int cap);
    [DllImport(Lib, EntryPoint = "caps_pipeline_series")] public static extern int PipelineSeries(IntPtr doc, int stride, CapsAnalyzeProgress? progress, IntPtr user, byte[]? json, int cap);
    [DllImport(Lib, EntryPoint = "caps_export_data")] public static extern int ExportData(IntPtr doc, [MarshalAs(UnmanagedType.LPUTF8Str)] string path, [MarshalAs(UnmanagedType.LPUTF8Str)] string format, [MarshalAs(UnmanagedType.LPUTF8Str)] string options);
    [DllImport(Lib, EntryPoint = "caps_export_preview")] public static extern int ExportPreview(IntPtr doc, [MarshalAs(UnmanagedType.LPUTF8Str)] string format, [MarshalAs(UnmanagedType.LPUTF8Str)] string options, int lines, byte[]? json, int cap);
    [DllImport(Lib, EntryPoint = "caps_bundle_preview")] public static extern int BundlePreview(IntPtr doc, [MarshalAs(UnmanagedType.LPUTF8Str)] string options, byte[]? json, int cap);
    [DllImport(Lib, EntryPoint = "caps_bundle_write")] public static extern int BundleWrite(IntPtr doc, [MarshalAs(UnmanagedType.LPUTF8Str)] string path, [MarshalAs(UnmanagedType.LPUTF8Str)] string options, in CapsCamera cam, in CapsRenderOpts opt);
    [DllImport(Lib, EntryPoint = "caps_pipeline_to_yaml")] public static extern int PipelineToYaml([MarshalAs(UnmanagedType.LPUTF8Str)] string json, [MarshalAs(UnmanagedType.LPUTF8Str)] string? name, [MarshalAs(UnmanagedType.LPUTF8Str)] string? file, [MarshalAs(UnmanagedType.LPUTF8Str)] string? topology, byte[]? yaml, int cap);
    [DllImport(Lib, EntryPoint = "caps_pipeline_from_yaml")] public static extern int PipelineFromYaml([MarshalAs(UnmanagedType.LPUTF8Str)] string yaml, byte[]? json, int cap);
    [DllImport(Lib, EntryPoint = "caps_set_python")] public static extern void SetPython([MarshalAs(UnmanagedType.LPUTF8Str)] string? packageDir, [MarshalAs(UnmanagedType.LPUTF8Str)] string? interpreter);
    [DllImport(Lib, EntryPoint = "caps_pipeline_catalogue")] public static extern int PipelineCatalogue(byte[]? json, int cap);
    [DllImport(Lib, EntryPoint = "caps_view_scale")] public static extern double ViewScale(IntPtr doc, in CapsCamera cam, in CapsRenderOpts opt);
    [DllImport(Lib, EntryPoint = "caps_bonded")] public static extern int Bonded(IntPtr doc, int index, [Out] int[]? idx, int cap);
    [DllImport(Lib, EntryPoint = "caps_molecule_index")] public static extern int MoleculeIndex(IntPtr doc, [Out] int[] mol, int cap);
    [DllImport(Lib, EntryPoint = "caps_neighbours")] public static extern int Neighbours(IntPtr doc, int index, int k, [Out] int[] idx, [Out] double[] dist);
}

/// <summary>An opened structure or trajectory. All calls are serialised; the core is not re-entrant per document.</summary>
public sealed class CapsDocument : IDisposable
{
    private IntPtr _h;
    private readonly object _lock = new();

    public string Path { get; }

    private CapsDocument(IntPtr h, string path) { _h = h; Path = path; }

    public static CapsDocument Open(string path, string? topology = null)
    {
        var h = Native.Open(path, topology);
        if (h == IntPtr.Zero) throw new InvalidOperationException(Native.LastError());
        return new CapsDocument(h, path);
    }

    /// <summary>A file's saved provenance (its .provenance.json), JSON with ok.</summary>
    public static string ProvenanceFile(string path) => Sized((b, c) => Native.ProvenanceFile(path, b, c));
    public static string ProvenanceCompare(string a, string b) => Sized((x, c) => Native.ProvenanceCompare(a, b, x, c));
    public static string ProvenanceBibtex(string manifest) => Sized((b, c) => Native.ProvenanceBibtex(manifest, b, c));
    /// <summary>The PNG's text chunks (a provenance manifest among them) as a JSON object.</summary>
    public static string PngText(string path) => Sized((b, c) => Native.PngText(path, b, c));

    /// <summary>Opens with the import choices (caps_import): bonds perceived / from the file / none, tolerance, bond orders,
    /// molecules, unwrap, cell.</summary>
    public static CapsDocument Import(string path, string? topology, string options)
    {
        var h = Native.Import(path, topology, options);
        if (h == IntPtr.Zero) throw new InvalidOperationException(Native.LastError());
        return new CapsDocument(h, path);
    }
    /// <summary>The import choices applied to frame 0 (counts, cell, first lines) as JSON.</summary>
    public static string ImportPreview(string path, string options) => Sized((b, c) => Native.ImportPreview(path, options, b, c));
    /// <summary>The preview's fragment (ten connected heavy atoms with their hydrogens) as a document.</summary>
    public static CapsDocument ImportFragment(string path, string options)
    {
        var h = Native.ImportFragment(path, options);
        if (h == IntPtr.Zero) throw new InvalidOperationException(Native.LastError());
        return new CapsDocument(h, path);
    }

    /// <summary>Staged open: progress(stage, fraction, detail) runs on the calling thread for each stage (0 format, 1 frame 0,
    /// 2 topology, 3 frames read); return false to stop reading, keeping the frames so far. maxFrames 1 reads frame 0 only.</summary>
    public static CapsDocument OpenStaged(string path, string? topology, int maxFrames, Func<int, double, string, bool>? progress)
    {
        CapsOpenProgress? cb = progress == null ? null : (st, f, d, _) => progress(st, f, Marshal.PtrToStringUTF8(d) ?? "") ? 0 : 1;
        var h = Native.OpenStaged(path, topology, maxFrames, cb, IntPtr.Zero);
        GC.KeepAlive(cb);
        if (h == IntPtr.Zero) throw new InvalidOperationException(Native.LastError());
        return new CapsDocument(h, path);
    }

    /// <summary>Moves the frames of other (the same file read in full) into this document; returns the frame count, or -1.</summary>
    public int AdoptFrames(CapsDocument other) => Native.AdoptFrames(_h, other._h);

    /// <summary>Grow a cell. The progress callback runs on the calling (worker) thread; return false to cancel.</summary>
    public static (CapsDocument Doc, string Report) Grow(CapsGrowOpts o, Func<int, int, int, bool>? progress, string label)
    {
        var report = new byte[2048];
        CapsProgress? cb = progress == null ? null : (d, t, r, _) => progress(d, t, r) ? 0 : 1;
        var h = Native.Grow(o, cb, IntPtr.Zero, report, report.Length);
        GC.KeepAlive(cb);
        if (h == IntPtr.Zero) throw new InvalidOperationException(Native.LastError());
        var text = System.Text.Encoding.UTF8.GetString(report).TrimEnd('\0').Trim();
        return (new CapsDocument(h, label), text);
    }

    /// <summary>Parses a SMILES without building: JSON with formula, mass, counts, or ok = false and the error.</summary>
    public static string SmilesInfo(string smiles)
    {
        var n = Native.SmilesInfo(smiles, null, 0);
        var buf = new byte[Math.Max(1, n)];
        Native.SmilesInfo(smiles, buf, buf.Length);
        return System.Text.Encoding.UTF8.GetString(buf, 0, Math.Max(0, n - 1));
    }

    private static string JsonCall(Func<byte[]?, int, int> f)
    {
        var n = f(null, 0);
        var buf = new byte[Math.Max(1, n)];
        f(buf, buf.Length);
        return System.Text.Encoding.UTF8.GetString(buf, 0, Math.Max(0, n - 1));
    }
    public static string UnitInfo(string smiles) => JsonCall((b, c) => Native.UnitInfo(smiles, b, c));
    public static string ChainPreview(string spec, ulong seed) => JsonCall((b, c) => Native.ChainPreview(spec, seed, b, c));

    /// <summary>Grows chains of a polymer spec (caps_grow_chains).</summary>
    public static (CapsDocument Doc, string Report) GrowChains(string spec, CapsGrowOpts o, Func<int, int, int, bool>? progress, string label)
    {
        var report = new byte[8192];
        CapsProgress? cb = progress == null ? null : (d, t, r, _) => progress(d, t, r) ? 0 : 1;
        var h = Native.GrowChains(spec, o, cb, IntPtr.Zero, report, report.Length);
        GC.KeepAlive(cb);
        if (h == IntPtr.Zero) throw new InvalidOperationException(Native.LastError());
        return (new CapsDocument(h, label), System.Text.Encoding.UTF8.GetString(report).TrimEnd('\0').Trim());
    }

    /// <summary>The terminations of a (hkl) plane of a CIF crystal (caps_surface_terminations, JSON).</summary>
    public static string SurfaceTerminations(string cif, int h, int k, int l) => JsonCall((b, c) => Native.SurfaceTerminations(cif, h, k, l, b, c));

    /// <summary>A slab cleaved from a CIF crystal (caps_surface_build).</summary>
    public static (CapsDocument Doc, string Report) SurfaceBuild(string cif, string options, string label)
    {
        var report = new byte[4096];
        var h = Native.SurfaceBuild(cif, options, report, report.Length);
        if (h == IntPtr.Zero) throw new InvalidOperationException(Native.LastError());
        return (new CapsDocument(h, label), System.Text.Encoding.UTF8.GetString(report).TrimEnd('\0').Trim());
    }

    /// <summary>A polymer film grown onto a slab (caps_interface_build); the slab is molecule 1 and held in Relax.</summary>
    public static (CapsDocument Doc, string Report) InterfaceBuild(string options, string spec, CapsGrowOpts o, Func<int, int, int, bool>? progress, string label)
    {
        var report = new byte[8192];
        CapsProgress? cb = progress == null ? null : (d, t, r, _) => progress(d, t, r) ? 0 : 1;
        var h = Native.InterfaceBuild(options, spec, o, cb, IntPtr.Zero, report, report.Length);
        GC.KeepAlive(cb);
        if (h == IntPtr.Zero) throw new InvalidOperationException(Native.LastError());
        return (new CapsDocument(h, label), System.Text.Encoding.UTF8.GetString(report).TrimEnd('\0').Trim());
    }

    /// <summary>Calls that do real work: one call into a large buffer, a second only when the answer did not fit.</summary>
    private static string JsonCallOnce(Func<byte[]?, int, int> f)
    {
        var buf = new byte[1 << 16];
        var n = f(buf, buf.Length);
        if (n > buf.Length) { buf = new byte[n]; f(buf, buf.Length); }
        return System.Text.Encoding.UTF8.GetString(buf, 0, Math.Max(0, Math.Min(n, buf.Length) - 1));
    }
    /// <summary>What a peptide spec builds, without the clean-up (caps_peptide_info).</summary>
    public static string PeptideInfo(string options) => JsonCallOnce((b, c) => Native.PeptideInfo(options, b, c));
    /// <summary>The sequence of the first record of a FASTA text (caps_fasta_sequence).</summary>
    public static string FastaSequence(string text) => JsonCall((b, c) => Native.FastaSequence(text, b, c));
    /// <summary>An all-atom peptide (caps_peptide_build).</summary>
    public static (CapsDocument Doc, string Report) PeptideBuild(string options, string label)
    {
        var report = new byte[4096];
        var h = Native.PeptideBuild(options, report, report.Length);
        if (h == IntPtr.Zero) throw new InvalidOperationException(Native.LastError());
        return (new CapsDocument(h, label), System.Text.Encoding.UTF8.GetString(report).TrimEnd('\0').Trim());
    }
    /// <summary>Solvents and salts for the solvation builder (caps_solvent_library).</summary>
    public static string SolventLibrary() => JsonCall(Native.SolventLibrary);
    /// <summary>The box and counts of a solvation, without packing (caps_solvate_plan); solute may be null.</summary>
    public static string SolvatePlan(CapsDocument? solute, string options)
    {
        if (solute == null) return JsonCallOnce((b, c) => Native.SolvatePlan(IntPtr.Zero, options, b, c));
        lock (solute._lock) return JsonCallOnce((b, c) => Native.SolvatePlan(solute._h, options, b, c));
    }
    /// <summary>Solvent and ions packed around the solute (caps_solvate): a new document.</summary>
    public static (CapsDocument Doc, string Report) Solvate(CapsDocument? solute, string options, Func<int, int, int, double, int, bool>? progress, string label)
    {
        var report = new byte[4096];
        CapsStageProgress? cb = progress == null ? null : (st, l, n, d, b, _) => progress(st, l, n, d, b) ? 0 : 1;
        IntPtr h;
        if (solute == null) h = Native.Solvate(IntPtr.Zero, options, cb, IntPtr.Zero, report, report.Length);
        else lock (solute._lock) h = Native.Solvate(solute._h, options, cb, IntPtr.Zero, report, report.Length);
        GC.KeepAlive(cb);
        GC.KeepAlive(solute);
        if (h == IntPtr.Zero) throw new InvalidOperationException(Native.LastError());
        return (new CapsDocument(h, label), System.Text.Encoding.UTF8.GetString(report).TrimEnd('\0').Trim());
    }
    /// <summary>The 530 space-group settings (caps_space_groups).</summary>
    public static string SpaceGroups() => JsonCall(Native.SpaceGroups);
    /// <summary>What a crystal spec builds (caps_crystal_info).</summary>
    public static string CrystalInfo(string spec) => JsonCallOnce((b, c) => Native.CrystalInfo(spec, b, c));
    /// <summary>The spec with its sites on their special positions (caps_crystal_symmetrize).</summary>
    public static string CrystalSymmetrize(string spec, double snap) => JsonCallOnce((b, c) => Native.CrystalSymmetrize(spec, snap, b, c));
    /// <summary>The space group and asymmetric unit of a spec's crystal or of a CIF file (caps_crystal_find_symmetry).</summary>
    public static string CrystalFindSymmetry(string? spec, string? cif, double tolerance) => JsonCallOnce((b, c) => Native.CrystalFindSymmetry(spec, cif, tolerance, b, c));
    /// <summary>A crystal from a space group, a lattice and an asymmetric unit (caps_crystal_build).</summary>
    public static (CapsDocument Doc, string Report) CrystalBuild(string spec, string label)
    {
        var report = new byte[4096];
        var h = Native.CrystalBuild(spec, report, report.Length);
        if (h == IntPtr.Zero) throw new InvalidOperationException(Native.LastError());
        return (new CapsDocument(h, label), System.Text.Encoding.UTF8.GetString(report).TrimEnd('\0').Trim());
    }

    /// <summary>A graphene sheet, nanotube or nanoparticle (caps_nano_build).</summary>
    public static (CapsDocument Doc, string Report) NanoBuild(string options, string label)
    {
        var report = new byte[4096];
        var h = Native.NanoBuild(options, report, report.Length);
        if (h == IntPtr.Zero) throw new InvalidOperationException(Native.LastError());
        return (new CapsDocument(h, label), System.Text.Encoding.UTF8.GetString(report).TrimEnd('\0').Trim());
    }

    /// <summary>The filler held at the centre of a cell with polymer chains grown around it (caps_nano_embed).</summary>
    public static (CapsDocument Doc, string Report) NanoEmbed(string options, string spec, CapsGrowOpts o, Func<int, int, int, bool>? progress, string label)
    {
        var report = new byte[8192];
        CapsProgress? cb = progress == null ? null : (d, t, r, _) => progress(d, t, r) ? 0 : 1;
        var h = Native.NanoEmbed(options, spec, o, cb, IntPtr.Zero, report, report.Length);
        GC.KeepAlive(cb);
        if (h == IntPtr.Zero) throw new InvalidOperationException(Native.LastError());
        return (new CapsDocument(h, label), System.Text.Encoding.UTF8.GetString(report).TrimEnd('\0').Trim());
    }

    /// <summary>A polymer blend grown component after component (caps_grow_blend).</summary>
    public static (CapsDocument Doc, string Report) GrowBlend(string options, CapsGrowOpts o, Func<int, int, int, bool>? progress, string label)
    {
        var report = new byte[8192];
        CapsProgress? cb = progress == null ? null : (d, t, r, _) => progress(d, t, r) ? 0 : 1;
        var h = Native.GrowBlend(options, o, cb, IntPtr.Zero, report, report.Length);
        GC.KeepAlive(cb);
        if (h == IntPtr.Zero) throw new InvalidOperationException(Native.LastError());
        return (new CapsDocument(h, label), System.Text.Encoding.UTF8.GetString(report).TrimEnd('\0').Trim());
    }

    /// <summary>Holds molecule `mol` in place in Relax (0: none).</summary>
    public void SetHeldMolecule(long mol) { lock (_lock) Native.SetHeldMolecule(_h, mol); }
    public long HeldMolecule() { lock (_lock) return Native.HeldMolecule(_h); }

    /// <summary>The file checks of this document as JSON (caps_file_checks).</summary>
    private static string Sized(Func<byte[]?, int, int> call)
    {
        var n = call(null, 0);
        if (n < 0) throw new InvalidOperationException(Native.LastError());
        var buf = new byte[Math.Max(1, n)];
        call(buf, buf.Length);
        return System.Text.Encoding.UTF8.GetString(buf, 0, Math.Max(0, n - 1));
    }

    /// <summary>Sets the visualize pipeline (JSON steps; null or "" clears it) and runs it on the shown frame.</summary>
    public void SetPipeline(string? json) { lock (_lock) { if (Native.PipelineSet(_h, json) != 0) throw new InvalidOperationException(Native.LastError()); } }
    public string PipelineResult() { lock (_lock) return Sized((b, c) => Native.PipelineResult(_h, b, c)); }
    public string PipelineParticles(string filter, int offset, int count) { lock (_lock) return Sized((b, c) => Native.PipelineParticles(_h, filter, offset, count, b, c)); }
    public string PipelineBonds(int offset, int count) { lock (_lock) return Sized((b, c) => Native.PipelineBonds(_h, offset, count, b, c)); }
    public static string PipelineCatalogue() => Sized(Native.PipelineCatalogue);
    public static string PipelineToYaml(string json, string? name, string? file, string? topology) => Sized((b, c) => Native.PipelineToYaml(json, name, file, topology, b, c));
    public static string PipelineFromYaml(string yaml) => Sized((b, c) => Native.PipelineFromYaml(yaml, b, c));
    /// <summary>What a file holds before opening it (scans a dump for its frames: call off the UI thread).</summary>
    public static string InspectFile(string path, string? topology) => Sized((b, c) => Native.InspectFile(path, topology, b, c));
    /// <summary>The bundle's files with sizes and hashes (figures are made on write). Runs the pipeline twice: off the UI thread.</summary>
    public string BundlePreview(string options) { lock (_lock) return Sized((b, c) => Native.BundlePreview(_h, options, b, c)); }
    public int BundleWrite(string path, string options, in CapsCamera cam, in CapsRenderOpts opt)
    {
        lock (_lock)
        {
            var n = Native.BundleWrite(_h, path, options, cam, opt);
            if (n < 0) throw new InvalidOperationException(Native.LastError());
            return n;
        }
    }
    public void ExportData(string path, string format, string options) { lock (_lock) { if (Native.ExportData(_h, path, format, options) != 0) throw new InvalidOperationException(Native.LastError()); } }
    /// <summary>Writes the export to a scratch file and returns its first lines, size and section counts (JSON). Writes twice; call off the UI thread.</summary>
    public string ExportPreview(string format, string options, int lines)
    {
        lock (_lock)
        {
            var buf = new byte[1 << 16];
            var n = Native.ExportPreview(_h, format, options, lines, buf, buf.Length);
            if (n < 0) throw new InvalidOperationException(Native.LastError());
            if (n > buf.Length) { buf = new byte[n]; Native.ExportPreview(_h, format, options, lines, buf, buf.Length); }
            return System.Text.Encoding.UTF8.GetString(buf, 0, Math.Max(0, Math.Min(n, buf.Length) - 1));
        }
    }
    /// <summary>The pipeline's attributes on every stride-th frame (runs the pipeline once to size, once to fill: call off the UI thread).</summary>
    public string PipelineSeries(int stride, Func<double, bool>? progress)
    {
        CapsAnalyzeProgress? cb = progress == null ? null : (_, f, _) => progress(f) ? 0 : 1;
        lock (_lock)
        {
            // one run fills a generous buffer; a second only if it was too small
            var buf = new byte[1 << 20];
            var n = Native.PipelineSeries(_h, stride, cb, IntPtr.Zero, buf, buf.Length);
            if (n < 0) throw new InvalidOperationException(Native.LastError());
            if (n > buf.Length) { buf = new byte[n]; Native.PipelineSeries(_h, stride, null, IntPtr.Zero, buf, buf.Length); }
            GC.KeepAlive(cb);
            return System.Text.Encoding.UTF8.GetString(buf, 0, Math.Max(0, Math.Min(n, buf.Length) - 1));
        }
    }

    public string FileChecks()
    {
        lock (_lock)
        {
            var n = Native.FileChecks(_h, null, 0);
            if (n < 0) throw new InvalidOperationException(Native.LastError());
            var buf = new byte[Math.Max(1, n)];
            Native.FileChecks(_h, buf, buf.Length);
            return System.Text.Encoding.UTF8.GetString(buf, 0, Math.Max(0, n - 1));
        }
    }

    /// <summary>Inserts copies of a molecule (SMILES) into the free space of the current frame (caps_insert_molecules).</summary>
    public string InsertMolecules(string smiles, int count, double tolerance, ulong seed)
    {
        var report = new byte[4096];
        lock (_lock) Check(Native.InsertMolecules(_h, smiles, count, tolerance, seed, report, report.Length));
        return System.Text.Encoding.UTF8.GetString(report).TrimEnd('\0').Trim();
    }

    /// <summary>The 2D drawing of a SMILES (JSON graph with coordinates, bond length 1).</summary>
    public static string SmilesDepict(string smiles)
    {
        var n = Native.SmilesDepict(smiles, null, 0);
        var buf = new byte[Math.Max(1, n)];
        Native.SmilesDepict(smiles, buf, buf.Length);
        return System.Text.Encoding.UTF8.GetString(buf, 0, Math.Max(0, n - 1));
    }

    /// <summary>SMILES from a graph in the caps_smiles_depict JSON form.</summary>
    public static string SmilesWrite(string graphJson)
    {
        var n = Native.SmilesWrite(graphJson, null, 0);
        if (n < 0) throw new InvalidOperationException(Native.LastError());
        var buf = new byte[Math.Max(1, n)];
        Native.SmilesWrite(graphJson, buf, buf.Length);
        return System.Text.Encoding.UTF8.GetString(buf, 0, Math.Max(0, n - 1));
    }

    /// <summary>A 3D molecule from SMILES, one frame per conformer (lowest energy first); ff: a caps-forcefield JSON
    /// with typing rules for the clean-up, or null. Returns the document and the JSON report.</summary>
    public static (CapsDocument Doc, string Report) BuildSmiles(string smiles, string? ff, int conformers, ulong seed, string label)
    {
        var report = new byte[65536];
        var h = Native.BuildSmiles(smiles, ff, new CapsBuildOpts { Conformers = conformers, Seed = seed }, report, report.Length);
        if (h == IntPtr.Zero) throw new InvalidOperationException(Native.LastError());
        return (new CapsDocument(h, label), System.Text.Encoding.UTF8.GetString(report).TrimEnd('\0'));
    }

    public void Save(string path) { lock (_lock) Check(Native.Save(_h, path)); }

    /// <summary>Pack from packmol-format text. Throws with the report attached when the tolerance cannot be met.</summary>
    public static (CapsDocument Doc, string Report) Pack(string text, string baseDir, Func<int, int, double, int, bool>? progress, string label)
    {
        var report = new byte[4096];
        CapsPackProgress? cb = progress == null ? null : (l, n, f, b, _) => progress(l, n, f, b) ? 0 : 1;
        var h = Native.Pack(System.Text.Encoding.UTF8.GetBytes(text + "\0"), baseDir, 0, cb, IntPtr.Zero, report, report.Length);
        GC.KeepAlive(cb);
        var rep = System.Text.Encoding.UTF8.GetString(report).TrimEnd('\0').Trim();
        if (h == IntPtr.Zero) throw new InvalidOperationException(Native.LastError() + (rep.Length > 0 ? "\n" + rep : ""));
        return (new CapsDocument(h, label), rep);
    }

    /// <summary>Relax the current frame in place (the document gains one frame per stage). The progress callback runs on
    /// the calling (worker) thread; return false to cancel. Returns whether the force tolerance was met, and the report.</summary>
    public (bool Converged, string Report) Relax(CapsRelaxOpts o, Func<int, int, int, double, double, double, bool>? progress)
    {
        lock (_lock)
        {
            var report = new byte[8192];
            CapsRelaxProgress? cb = progress == null ? null : (st, n, it, e, f, d, _) => progress(st, n, it, e, f, d) ? 0 : 1;
            var rc = Native.Relax(_h, o, cb, IntPtr.Zero, report, report.Length);
            GC.KeepAlive(cb);
            Check(rc);
            return (rc == 0, System.Text.Encoding.UTF8.GetString(report).TrimEnd('\0').Trim());
        }
    }

    /// <summary>Molecular dynamics from the current frame (the document becomes the recorded trajectory). The progress
    /// callback runs on the calling (worker) thread with each thermo row; return false to cancel.</summary>
    public string Md(CapsMdOpts o, Func<CapsThermo, long, bool>? progress)
    {
        lock (_lock)
        {
            var report = new byte[8192];
            CapsMdProgress? cb = progress == null ? null : (in CapsThermo r, long n, IntPtr _) => progress(r, n) ? 0 : 1;
            var rc = Native.Md(_h, o, cb, IntPtr.Zero, report, report.Length);
            GC.KeepAlive(cb);
            Check(rc);
            return System.Text.Encoding.UTF8.GetString(report).TrimEnd('\0').Trim();
        }
    }

    /// <summary>Text of a named protocol (larsen21, annealing, pushoff).</summary>
    public static string ProtocolText(string name, CapsProtocolParams p)
    {
        var buf = new byte[16384];
        if (Native.ProtocolText(name, p, buf, buf.Length) < 0) throw new InvalidOperationException(Native.LastError());
        return System.Text.Encoding.UTF8.GetString(buf).TrimEnd('\0');
    }

    /// <summary>Run a protocol from the current frame; the document becomes the recorded trajectory.
    /// Returns whether the convergence checks passed (true when they were not asked for) and the report.</summary>
    public (bool Converged, string Report) Equilibrate(string protocol, CapsEquilOpts o, Func<int, int, string, CapsThermo, bool>? progress)
    {
        lock (_lock)
        {
            var report = new byte[16384];
            var text = System.Text.Encoding.UTF8.GetBytes(protocol + "\0");
            CapsEquilProgress? cb = progress == null ? null
                : (int st, int n, IntPtr label, in CapsThermo r, IntPtr _) => progress(st, n, Marshal.PtrToStringUTF8(label) ?? "", r) ? 0 : 1;
            var rc = Native.Equilibrate(_h, text, o, cb, IntPtr.Zero, report, report.Length);
            GC.KeepAlive(cb);
            Check(rc);
            return (rc == 0, System.Text.Encoding.UTF8.GetString(report).TrimEnd('\0').Trim());
        }
    }

    public (int[] N, double[] Ratio, int Chains, double B2) InternalDistances()
    {
        lock (_lock)
        {
            var n = new int[20000];
            var r = new double[20000];
            var m = Native.InternalDistances(_h, n, r, n.Length, out var chains, out var b2);
            Check(m);
            return (n[..m], r[..m], chains, b2);
        }
    }

    public static string ReactionTemplate(string name)
    {
        var buf = new byte[8192];
        if (Native.ReactionTemplate(name, buf, buf.Length) < 0) throw new InvalidOperationException(Native.LastError());
        return System.Text.Encoding.UTF8.GetString(buf).TrimEnd('\0');
    }

    /// <summary>React the current frame; progress runs on the worker thread with each cycle; return false to cancel.</summary>
    public string React(string templates, CapsReactOpts o, Func<CapsReactCycle, bool>? progress)
    {
        lock (_lock)
        {
            var report = new byte[8192];
            CapsReactProgress? cb = progress == null ? null : (in CapsReactCycle r, IntPtr _) => progress(r) ? 0 : 1;
            var rc = Native.React(_h, System.Text.Encoding.UTF8.GetBytes(templates + "\0"), o, cb, IntPtr.Zero, report, report.Length);
            GC.KeepAlive(cb);
            Check(rc);
            return System.Text.Encoding.UTF8.GetString(report).TrimEnd('\0').Trim();
        }
    }

    public void SaveTrajectory(string path) { lock (_lock) Check(Native.SaveTrajectory(_h, path)); }

    /// <summary>GAFF typing, term counts and energy of the current frame, as text.</summary>
    public string FieldInfo()
    {
        lock (_lock)
        {
            var text = new byte[4096];
            Check(Native.FieldInfo(_h, text, text.Length));
            return System.Text.Encoding.UTF8.GetString(text).TrimEnd('\0').Trim();
        }
    }

    // ---- CAPS Field: each call returns true when the assignment is complete (every atom typed, every parameter found)

    /// <summary>Types and parameterises the structure with a library force field. charges: 0 force field, 1 Gasteiger, 2 file.</summary>
    public bool FieldAssign(string ffPath, string? rulesPath, int charges) { lock (_lock) return CheckField(Native.FieldAssign(_h, ffPath, rulesPath, charges)); }
    public bool FieldOverride(int index, string? type) { lock (_lock) return CheckField(Native.FieldOverride(_h, index, type)); }
    public bool FieldAddRule(string kind, string types, string style, string pars) { lock (_lock) return CheckField(Native.FieldAddRule(_h, kind, types, style, pars)); }
    public bool FieldImport(string path) { lock (_lock) return CheckField(Native.FieldImport(_h, path)); }
    public bool FieldRemoveRules() { lock (_lock) return CheckField(Native.FieldRemoveRules(_h)); }
    public void FieldClear() { lock (_lock) Check(Native.FieldClear(_h)); }
    public void FieldTypesFile(string path) { lock (_lock) Check(Native.FieldTypesFile(_h, path)); }

    /// <summary>The assignment as JSON, or "" when there is none.</summary>
    public string FieldReport()
    {
        lock (_lock)
        {
            var n = Native.FieldReport(_h, null, 0);
            if (n <= 1) return "";
            var buf = new byte[n];
            Native.FieldReport(_h, buf, n);
            return System.Text.Encoding.UTF8.GetString(buf, 0, n - 1);
        }
    }

    /// <summary>Analyze properties over the frames (ids comma-separated, as caps analyze). The progress callback runs on the
    /// calling (worker) thread with what is being computed and the fraction done; return false to cancel. Returns the JSON
    /// report: {frames, of, atoms, properties: [...]}.</summary>
    public string Analyze(string props, CapsAnalyzeOpts o, Func<string, double, bool>? progress) => Analyze(props, o, default, progress);

    /// <summary>Analyze with the mechanics and Tg protocols (cij_strain, cij_fluct, tensile, tg); the document is not changed.</summary>
    public string Analyze(string props, CapsAnalyzeOpts o, CapsMechOpts m, Func<string, double, bool>? progress)
    {
        lock (_lock)
        {
            CapsAnalyzeProgress? cb = progress == null ? null : (w, f, _) => progress(Marshal.PtrToStringUTF8(w) ?? "", f) ? 1 : 0;
            var rc = Native.AnalyzeEx(_h, props, o, m, cb, IntPtr.Zero);
            GC.KeepAlive(cb);
            Check(rc);
            var n = Native.AnalyzeReport(_h, null, 0);
            if (n <= 1) return "";
            var buf = new byte[n];
            Native.AnalyzeReport(_h, buf, n);
            return System.Text.Encoding.UTF8.GetString(buf, 0, n - 1);
        }
    }

    /// <summary>Convergence checks of the last Equilibrate run (JSON), or "".</summary>
    public string EquilibrateChecks()
    {
        lock (_lock)
        {
            var n = Native.EquilibrateChecks(_h, null, 0);
            if (n <= 1) return "";
            var buf = new byte[n];
            Native.EquilibrateChecks(_h, buf, n);
            return System.Text.Encoding.UTF8.GetString(buf, 0, n - 1);
        }
    }

    /// <summary>The LAMMPS input setup (styles, read_data, neighbour list) for the data file Save writes.</summary>
    public string LammpsInput(string dataName)
    {
        lock (_lock)
        {
            var n = Native.LammpsInput(_h, dataName, null, 0);
            if (n < 0) throw new InvalidOperationException(Native.LastError());
            var buf = new byte[n];
            Native.LammpsInput(_h, dataName, buf, n);
            return System.Text.Encoding.UTF8.GetString(buf, 0, n - 1);
        }
    }

    private static bool CheckField(int rc) { Check(rc); return rc == 0; }

    private static void Check(int rc) { if (rc < 0) throw new InvalidOperationException(Native.LastError()); }

    public CapsSummary Summary() { lock (_lock) { Check(Native.Summary(_h, out var s)); return s; } }
    public void SetFrame(long f) { lock (_lock) { Alive(); Check(Native.SetFrame(_h, f)); } }
    public void SetWrap(bool wrap) { lock (_lock) Check(Native.SetWrap(_h, wrap ? 1 : 0)); }
    public CapsAtomInfo Atom(int i) { lock (_lock) { Check(Native.Atom(_h, i, out var a)); return a; } }

    public IReadOnlyList<string> Notes()
    {
        lock (_lock)
        {
            var n = Native.NoteCount(_h);
            var list = new List<string>(n);
            for (var k = 0; k < n; k++) list.Add(Native.Note(_h, k));
            return list;
        }
    }

    public unsafe void Render(in CapsCamera cam, in CapsRenderOpts opt, byte[] rgba)
    {
        lock (_lock)
        {
            Alive();
            fixed (byte* p = rgba) Check(Native.Render(_h, cam, opt, p));
        }
    }

    public int Pick(int x, int y) { lock (_lock) return Native.Pick(_h, x, y); }
    /// <summary>H-bonds, contacts, clashes and checks of the frame (caps_interactions); drawn until ClearChecks.</summary>
    public string Interactions(string options) { lock (_lock) { Alive(); return JsonCallOnce((b, c) => Native.Interactions(_h, options, b, c)); } }
    public void ClearChecks() { lock (_lock) { Alive(); Native.ClearChecks(_h); } }
    /// <summary>One structure edit (caps_edit), JSON {ok, error, what, atoms, added}.</summary>
    public string Edit(string json) { lock (_lock) { Alive(); return JsonCallOnce((b, c) => Native.Edit(_h, json, b, c)); } }
    /// <summary>Undo (redo = false) or redo the last edit; false when there is none.</summary>
    public bool Undo(bool redo) { lock (_lock) { Alive(); return Native.Undo(_h, redo ? 1 : 0) == 0; } }
    public string History() { lock (_lock) return JsonCall((b, c) => Native.History(_h, b, c)); }
    /// <summary>Changes the selection (caps_select), JSON {ok, error, count, matched}.</summary>
    public string Select(string json) { lock (_lock) { Alive(); return JsonCallOnce((b, c) => Native.Select(_h, json, b, c)); } }
    public string SelectionJson() { lock (_lock) return JsonCallOnce((b, c) => Native.Selection(_h, b, c)); }
    public string Tacticity() { lock (_lock) return JsonCallOnce((b, c) => Native.Tacticity(_h, b, c)); }
    /// <summary>A torsion scan of the current frame (caps_torsion_scan), JSON; progress(done, total) → false stops.</summary>
    public string TorsionScan(string options, Func<int, int, bool>? progress)
    {
        CapsSeriesProgress? cb = progress == null ? null : (d, t, _) => progress(d, t) ? 0 : 1;
        string r;
        lock (_lock) { Alive(); r = JsonCallOnce((b, c) => Native.TorsionScan(_h, options, cb, IntPtr.Zero, b, c)); }
        GC.KeepAlive(cb);
        return r;
    }
    /// <summary>Puts scan point k into the current frame (−1: the frame before the scan).</summary>
    public void TorsionShow(int k) { lock (_lock) { Alive(); Check(Native.TorsionShow(_h, k)); } }
    public int[]? DefaultTorsion() { var a = new int[4]; lock (_lock) { Alive(); return Native.DefaultTorsion(_h, a) == 0 ? a : null; } }
    /// <summary>Per-frame series for the trajectory player (caps_trajectory_series), JSON.</summary>
    public string TrajectorySeries(string options) { lock (_lock) { Alive(); return JsonCallOnce((b, c) => Native.TrajectorySeries(_h, options, IntPtr.Zero, IntPtr.Zero, b, c)); } }
    /// <summary>Positions shown averaged over `window` frames (1: off).</summary>
    public void SetSmoothing(int window) { lock (_lock) { Alive(); Native.SetSmoothing(_h, window); } }
    /// <summary>Styles, colours, surfaces and polyhedra of the view (caps_set_appearance).</summary>
    public void SetAppearance(string json) { lock (_lock) { Alive(); if (Native.SetAppearance(_h, json) != 0) throw new InvalidOperationException(Native.LastError()); } }
    public string AppearanceInfo() { lock (_lock) return JsonCall((b, c) => Native.AppearanceInfo(_h, b, c)); }
    /// <summary>One label per atom of the current frame: element, rs, type, charge or name (caps_atom_labels).</summary>
    public string AtomLabels(string kind) { lock (_lock) return JsonCallOnce((b, c) => Native.AtomLabels(_h, kind, b, c)); }
    /// <summary>x, y (pixels) and visibility of every atom after the last render with the same camera and options.</summary>
    public float[] ProjectAtoms(in CapsCamera cam, in CapsRenderOpts opt, int atoms)
    {
        var buf = new float[Math.Max(1, atoms) * 3];
        lock (_lock) { Alive(); var n = Native.ProjectAtoms(_h, cam, opt, buf, atoms); if (n < 0) throw new InvalidOperationException(Native.LastError()); }
        return buf;
    }
    /// <summary>The steps that produced this structure (caps-manifest/1.0).</summary>
    public string Provenance() { lock (_lock) { Alive(); return Sized((b, c) => Native.Provenance(_h, b, c)); } }
    public void ExportPng(in CapsCamera cam, in CapsRenderOpts opt, string path) { lock (_lock) { Alive(); Check(Native.ExportPng(_h, cam, opt, path)); } }
    /// <summary>PNG with 8 or 16 bits, dpi, colour profile and the provenance manifest (caps_export_image); overlay is
    /// a width × height straight-alpha RGBA layer (labels, measurements) or null.</summary>
    public void ExportImage(in CapsCamera cam, in CapsRenderOpts opt, string path, string options, byte[]? overlay)
    {
        lock (_lock) { Alive(); Check(Native.ExportImage(_h, cam, opt, path, options, overlay)); }
    }
    /// <summary>An animated PNG or a PNG sequence of the frames, or a turntable (caps_export_movie); returns the frames
    /// written. progress(done, total) → false stops.</summary>
    public int ExportMovie(in CapsCamera cam, in CapsRenderOpts opt, string path, string options, Func<int, int, bool>? progress)
    {
        CapsSeriesProgress? cb = progress == null ? null : (d, t, _) => progress(d, t) ? 0 : 1;
        int n;
        lock (_lock) { Alive(); n = Native.ExportMovie(_h, cam, opt, path, options, cb, IntPtr.Zero); }
        GC.KeepAlive(cb);
        if (n < 0) throw new InvalidOperationException(Native.LastError());
        return n;
    }
    /// <summary>Throws once the document is closed (a view may still hold it while a newer one replaces it).</summary>
    private void Alive() { if (_h == IntPtr.Zero) throw new ObjectDisposedException(nameof(CapsDocument)); }
    public void ExportSvg(in CapsCamera cam, in CapsRenderOpts opt, string path) { lock (_lock) Check(Native.ExportSvg(_h, cam, opt, path)); }

    public CapsMolecule[] Molecules()
    {
        lock (_lock)
        {
            var n = Native.Molecules(_h, [], 0);
            Check(n);
            var arr = new CapsMolecule[n];
            Check(Native.Molecules(_h, arr, n));
            return arr;
        }
    }

    public (double Lo, double Hi) PropertyRange()
    {
        lock (_lock) { Check(Native.PropertyRange(_h, out var lo, out var hi)); return (lo, hi); }
    }

    public double Measure(int[] idx)
    {
        lock (_lock) { Check(Native.Measure(_h, idx, idx.Length, out var v)); return v; }
    }

    public (double R, double G)[] Rdf(int elemA, int elemB, double rmax, double dr, bool interOnly)
    {
        lock (_lock)
        {
            var cap = (int)(rmax / dr) + 1;
            var r = new double[cap];
            var g = new double[cap];
            var n = Native.Rdf(_h, elemA, elemB, rmax, dr, interOnly ? 1 : 0, r, g, cap);
            Check(n);
            return Enumerable.Range(0, n).Select(k => (r[k], g[k])).ToArray();
        }
    }

    /// <summary>Pixels per Å at the focal plane for an image of opt's size (exact when orthographic): scale bars.</summary>
    public double ViewScale(in CapsCamera cam, in CapsRenderOpts opt) { lock (_lock) return Native.ViewScale(_h, cam, opt); }

    /// <summary>Atoms bonded to atom i.</summary>
    public int[] Bonded(int i)
    {
        lock (_lock)
        {
            var buf = new int[16];
            var n = Native.Bonded(_h, i, buf, buf.Length);
            if (n < 0) throw new InvalidOperationException(Native.LastError());
            if (n > buf.Length) { buf = new int[n]; Native.Bonded(_h, i, buf, n); }
            return buf[..n];
        }
    }

    /// <summary>The molecule (0-based, connected by bonds) of every atom, and the number of molecules.</summary>
    public (int[] Mol, int Count) MoleculeIndex(int atoms)
    {
        lock (_lock)
        {
            var m = new int[atoms];
            var n = Native.MoleculeIndex(_h, m, atoms);
            if (n < 0) throw new InvalidOperationException(Native.LastError());
            return (m, n);
        }
    }

    public (int Index, double Distance)[] Neighbours(int i, int k)
    {
        lock (_lock)
        {
            var idx = new int[k];
            var d = new double[k];
            var n = Native.Neighbours(_h, i, k, idx, d);
            Check(n);
            return Enumerable.Range(0, n).Select(q => (idx[q], d[q])).ToArray();
        }
    }

    public void Dispose()
    {
        lock (_lock)
        {
            if (_h != IntPtr.Zero) { Native.Close(_h); _h = IntPtr.Zero; }
        }
    }
}
