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
public delegate int CapsMdProgress(in CapsThermo row, long steps, IntPtr user);

[UnmanagedFunctionPointer(CallingConvention.Cdecl)]
public delegate int CapsRelaxProgress(int stage, int stages, int iteration, double energy, double fmax, double density, IntPtr user);

internal static class Native
{
    private const string Lib = "caps";

    [DllImport(Lib, EntryPoint = "caps_abi_version")] public static extern int AbiVersion();
    [DllImport(Lib, EntryPoint = "caps_last_error")] private static extern IntPtr LastErrorPtr();
    public static string LastError() => Marshal.PtrToStringUTF8(LastErrorPtr()) ?? "";

    [DllImport(Lib, EntryPoint = "caps_open_staged")] public static extern IntPtr OpenStaged([MarshalAs(UnmanagedType.LPUTF8Str)] string path, [MarshalAs(UnmanagedType.LPUTF8Str)] string? topology, int maxFrames, CapsOpenProgress? progress, IntPtr user);
    [DllImport(Lib, EntryPoint = "caps_adopt_frames")] public static extern int AdoptFrames(IntPtr dst, IntPtr src);
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
    public void ExportPng(in CapsCamera cam, in CapsRenderOpts opt, string path) { lock (_lock) { Alive(); Check(Native.ExportPng(_h, cam, opt, path)); } }
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
