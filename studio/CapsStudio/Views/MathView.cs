using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using System.Text;
using Avalonia;
using Avalonia.Controls;
using Avalonia.Controls.Documents;
using Avalonia.Input.Platform;
using Avalonia.Media;

namespace CapsStudio.Views;

/// <summary>Typeset mathematics for the theory manual: a small subset of LaTeX (fractions, sub- and superscripts, roots,
/// sums and integrals with limits, sized brackets, accents, upright function names and text) laid out the way TeX does
/// it, italic single-letter variables, upright numbers, names and upper-case Greek. Long equations wrap at their \quad
/// separators; \\ starts a new line. Symbol mode reads plain symbols (τ_P, k_max, N_f k_B T); Prose mode keeps text
/// upright and only lowers _sub and raises ^sup.</summary>
public sealed class MathView : Control
{
    public enum MathMode { Tex, Symbol, Prose }

    public static readonly StyledProperty<string?> TexProperty = AvaloniaProperty.Register<MathView, string?>(nameof(Tex));
    public static readonly StyledProperty<MathMode> ModeProperty = AvaloniaProperty.Register<MathView, MathMode>(nameof(Mode));
    public static readonly StyledProperty<double> FontSizeProperty = AvaloniaProperty.Register<MathView, double>(nameof(FontSize), 19);
    public static readonly StyledProperty<IBrush?> ForegroundProperty = TextElement.ForegroundProperty.AddOwner<MathView>();
    public static readonly StyledProperty<FontFamily> FontFamilyProperty = TextElement.FontFamilyProperty.AddOwner<MathView>();

    public string? Tex { get => GetValue(TexProperty); set => SetValue(TexProperty, value); }
    public MathMode Mode { get => GetValue(ModeProperty); set => SetValue(ModeProperty, value); }
    public double FontSize { get => GetValue(FontSizeProperty); set => SetValue(FontSizeProperty, value); }
    public IBrush? Foreground { get => GetValue(ForegroundProperty); set => SetValue(ForegroundProperty, value); }

    static MathView()
    {
        AffectsMeasure<MathView>(TexProperty, ModeProperty, FontSizeProperty);
        AffectsRender<MathView>(ForegroundProperty);
    }

    public MathView()
    {
        var copyTex = new MenuItem { Header = "Copy LaTeX" };
        copyTex.Click += async (_, _) => { if (TopLevel.GetTopLevel(this)?.Clipboard is { } c) await c.SetTextAsync(Tex ?? ""); };
        ContextMenu = new ContextMenu { ItemsSource = new[] { copyTex } };
    }

    private List<(Box line, double indent)> _lines = new();
    private Node? _tree;
    private string? _parsed;
    private MathMode _parsedMode;

    private Node Tree()
    {
        if (_tree == null || _parsed != Tex || _parsedMode != Mode)
        {
            _parsed = Tex;
            _parsedMode = Mode;
            try { _tree = MathParser.Parse(Tex ?? "", Mode); }
            catch (Exception e) { _tree = new Seq([new Sym(Tex ?? "", SymKind.Text), new Sym("  (" + e.Message + ")", SymKind.Text)]); }
        }
        return _tree;
    }

    protected override Size MeasureOverride(Size available)
    {
        // maths in a serif with italics and Greek; prose in the interface's own font
        var fam = Mode == MathMode.Prose ? GetValue(FontFamilyProperty) : new FontFamily("STIX Two Text, Cambria Math, Cambria, Times New Roman, DejaVu Serif, serif");
        var lay = new Layout(fam, FontSize, Foreground ?? Brushes.Black);
        _lines = lay.Lines(Tree(), double.IsInfinity(available.Width) ? double.MaxValue : Math.Max(40, available.Width));
        double w = 0, h = 0;
        for (var i = 0; i < _lines.Count; ++i)
        {
            var (b, ind) = _lines[i];
            w = Math.Max(w, ind + b.W);
            h += b.A + b.D + (i > 0 ? LineGap : 0);
        }
        return new Size(Math.Ceiling(w) + 1, Math.Ceiling(h) + 1);
    }

    private double LineGap => 0.55 * FontSize;

    public override void Render(DrawingContext ctx)
    {
        var brush = Foreground ?? Brushes.Black;
        double y = 0;
        for (var i = 0; i < _lines.Count; ++i)
        {
            var (b, ind) = _lines[i];
            if (i > 0) y += LineGap;
            y += b.A;
            b.Draw(ctx, ind, y, brush);
            y += b.D;
        }
    }

    // ================================================================ the parsed tree
    public enum SymKind { Italic, Upright, Number, Bin, Rel, Punct, Open, Close, Func, Text }
    public abstract record Node;
    public sealed record Seq(List<Node> Items) : Node;
    public sealed record Sym(string Text, SymKind Kind) : Node;
    public sealed record Frac(Node Num, Node Den, bool Small) : Node;
    public sealed record Sqrt(Node Body) : Node;
    public sealed record Scripts(Node Base, Node? Sup, Node? Sub) : Node;
    public sealed record BigOp(string Glyph, bool Integral) : Node;
    public sealed record LimFunc(string Name) : Node;
    public sealed record Fenced(string Open, string Close, Node Body) : Node;
    public sealed record Accent(string Kind, Node Body) : Node;
    public sealed record Space(double Em, bool Breakable) : Node;
    public sealed record LineBreak : Node;

    // ================================================================ parsing
    public static class MathParser
    {
        private static readonly Dictionary<string, string> Greek = new()
        {
            ["alpha"] = "α", ["beta"] = "β", ["gamma"] = "γ", ["delta"] = "δ", ["epsilon"] = "ε", ["varepsilon"] = "ε", ["zeta"] = "ζ",
            ["eta"] = "η", ["theta"] = "θ", ["iota"] = "ι", ["kappa"] = "κ", ["lambda"] = "λ", ["mu"] = "μ", ["nu"] = "ν", ["xi"] = "ξ",
            ["pi"] = "π", ["rho"] = "ρ", ["sigma"] = "σ", ["tau"] = "τ", ["upsilon"] = "υ", ["phi"] = "φ", ["varphi"] = "φ", ["chi"] = "χ",
            ["psi"] = "ψ", ["omega"] = "ω", ["ell"] = "ℓ",
            ["Gamma"] = "Γ", ["Delta"] = "Δ", ["Theta"] = "Θ", ["Lambda"] = "Λ", ["Xi"] = "Ξ", ["Pi"] = "Π", ["Sigma"] = "Σ",
            ["Phi"] = "Φ", ["Psi"] = "Ψ", ["Omega"] = "Ω",
        };
        private static readonly Dictionary<string, (string, SymKind)> Symbols = new()
        {
            ["cdot"] = ("·", SymKind.Bin), ["times"] = ("×", SymKind.Bin), ["pm"] = ("±", SymKind.Bin), ["mp"] = ("∓", SymKind.Bin),
            ["le"] = ("≤", SymKind.Rel), ["leq"] = ("≤", SymKind.Rel), ["ge"] = ("≥", SymKind.Rel), ["geq"] = ("≥", SymKind.Rel),
            ["ne"] = ("≠", SymKind.Rel), ["neq"] = ("≠", SymKind.Rel), ["approx"] = ("≈", SymKind.Rel), ["equiv"] = ("≡", SymKind.Rel),
            ["to"] = ("→", SymKind.Rel), ["rightarrow"] = ("→", SymKind.Rel), ["leftarrow"] = ("←", SymKind.Rel), ["in"] = ("∈", SymKind.Rel),
            ["propto"] = ("∝", SymKind.Rel), ["sim"] = ("∼", SymKind.Rel),
            ["infty"] = ("∞", SymKind.Upright), ["nabla"] = ("∇", SymKind.Upright), ["partial"] = ("∂", SymKind.Upright),
            ["dots"] = ("…", SymKind.Upright), ["ldots"] = ("…", SymKind.Upright), ["cdots"] = ("⋯", SymKind.Upright),
            ["circ"] = ("∘", SymKind.Upright), ["prime"] = ("′", SymKind.Upright), ["#"] = ("#", SymKind.Upright), ["%"] = ("%", SymKind.Upright),
            ["{"] = ("{", SymKind.Open), ["}"] = ("}", SymKind.Close), ["langle"] = ("⟨", SymKind.Open), ["rangle"] = ("⟩", SymKind.Close),
            ["|"] = ("‖", SymKind.Upright), ["AA"] = ("Å", SymKind.Upright),
        };
        private static readonly HashSet<string> Funcs = ["exp", "ln", "log", "cos", "sin", "tan", "erfc", "erf", "sign", "tr", "dev", "det", "arccos"];
        private static readonly HashSet<string> Limits = ["lim", "min", "max", "argmax", "argmin", "sup", "inf"];
        private static readonly Dictionary<string, double> Spaces = new() { [","] = 0.17, [":"] = 0.22, [";"] = 0.28, ["!"] = -0.17, [" "] = 0.3, ["quad"] = 1, ["qquad"] = 2 };

        /// <summary>The tree of a formula; throws on an unknown command or unbalanced braces.</summary>
        public static Node Parse(string s, MathMode mode = MathMode.Tex)
        {
            var p = new P(s, mode);
            var n = p.Row(top: true);
            if (p.I < s.Length) throw new FormatException($"unexpected '{s[p.I]}' at {p.I}");
            return n;
        }

        private sealed class P(string s, MathMode mode)
        {
            public int I;
            private bool More => I < s.Length;

            public Node Row(bool top = false, bool inText = false, string? until = null)
            {
                var items = new List<Node>();
                while (More)
                {
                    var c = s[I];
                    if (c == '}') { if (top) throw new FormatException($"unmatched '}}' at {I}"); break; }
                    if (until != null && s.AsSpan(I).StartsWith(until)) break;
                    if (c == '^' || c == '_')
                    {
                        ++I;
                        var arg = ScriptArg();
                        var last = items.Count > 0 ? items[^1] : new Sym("", SymKind.Upright);
                        if (items.Count > 0) items.RemoveAt(items.Count - 1);
                        items.Add(last is Scripts sc
                            ? (c == '^' ? sc with { Sup = arg } : sc with { Sub = arg })
                            : new Scripts(last, c == '^' ? arg : null, c == '_' ? arg : null));
                        continue;
                    }
                    if (c == '{') { ++I; items.Add(Group(inText)); continue; }
                    if (c == '\\') { items.Add(Command(inText || mode == MathMode.Prose)); continue; }
                    if (char.IsWhiteSpace(c))
                    {
                        ++I;
                        while (More && char.IsWhiteSpace(s[I])) ++I;
                        if (inText || mode != MathMode.Tex) items.Add(new Space(0.28, top));
                        continue;
                    }
                    items.Add(Atom(inText || mode == MathMode.Prose));
                }
                return items.Count == 1 ? items[0] : new Seq(items);
            }

            private Node Group(bool inText = false)
            {
                var n = Row(inText: inText);
                if (!More || s[I] != '}') throw new FormatException("missing '}'");
                ++I;
                return n;
            }

            private Node Arg(bool inText = false)
            {
                while (More && s[I] == ' ') ++I;
                if (!More) throw new FormatException("missing argument");
                if (s[I] == '{') { ++I; return Group(inText); }
                if (s[I] == '\\') return Command(inText);
                return Atom(inText);
            }

            // _x, _{…}; in symbol and prose modes _word takes the whole word (k_max, t_ramp)
            private Node ScriptArg()
            {
                if (More && s[I] != '{' && s[I] != '\\' && mode != MathMode.Tex && char.IsLetterOrDigit(s[I]))
                {
                    var j = I;
                    while (j < s.Length && (char.IsLetterOrDigit(s[j]) || s[j] == ',' && j + 1 < s.Length && char.IsLetterOrDigit(s[j + 1]) && s[j - 1] != ',')) ++j;
                    var w = s[I..j];
                    I = j;
                    return mode == MathMode.Prose ? new Sym(w, SymKind.Text) : WordOrLetters(w);
                }
                return Arg();
            }

            private static Node WordOrLetters(string w) =>
                w.Length > 1 && w.All(char.IsLetter) ? new Sym(w, SymKind.Upright)
                : w.All(char.IsDigit) ? new Sym(w, SymKind.Number)
                : w.Length == 1 ? new Sym(w, char.IsLetter(w[0]) && IsItalicLetter(w[0]) ? SymKind.Italic : SymKind.Upright)
                : new Seq(w.Select(ch => (Node)new Sym(ch.ToString(), char.IsDigit(ch) ? SymKind.Number : IsItalicLetter(ch) ? SymKind.Italic : SymKind.Upright)).ToList());

            private Node Atom(bool inText)
            {
                var c = s[I];
                if (inText)
                {
                    var j = I;
                    while (j < s.Length && s[j] is not ('\\' or '{' or '}' or '^' or '_') && !char.IsWhiteSpace(s[j])) ++j;
                    if (j == I) j = I + 1;
                    var t = s[I..j];
                    I = j;
                    return new Sym(t, SymKind.Text);
                }
                if (char.IsDigit(c) || c == '.' && I + 1 < s.Length && char.IsDigit(s[I + 1]))
                {
                    var j = I;
                    while (j < s.Length && (char.IsDigit(s[j]) || s[j] == '.' && j + 1 < s.Length && char.IsDigit(s[j + 1]))) ++j;
                    var t = s[I..j];
                    I = j;
                    return new Sym(t, SymKind.Number);
                }
                if (mode == MathMode.Symbol && char.IsLetter(c))
                {
                    // symbol mode: a word of two or more letters is a name, upright
                    var j = I;
                    while (j < s.Length && char.IsLetter(s[j])) ++j;
                    var w = s[I..j];
                    I = j;
                    return WordOrLetters(w);
                }
                ++I;
                if (char.IsSurrogate(c) && More) { ++I; return new Sym(s.Substring(I - 2, 2), SymKind.Upright); }
                var str = c.ToString();
                // a combining mark (K̄, ṗ written as Unicode) stays with its letter
                if (More && CharUnicodeInfo.GetUnicodeCategory(s[I]) == UnicodeCategory.NonSpacingMark) { str += s[I]; ++I; }
                return c switch
                {
                    '+' or '−' or '*' or '·' or '×' or '±' => new Sym(c == '*' ? "∗" : str, SymKind.Bin),
                    '-' => new Sym("−", SymKind.Bin),
                    '=' or '<' or '>' or '≤' or '≥' or '≠' or '≈' or '→' or '←' or '∈' or ':' => new Sym(str, SymKind.Rel),
                    ',' or ';' => new Sym(str, SymKind.Punct),
                    '(' or '[' or '⟨' => new Sym(str, SymKind.Open),
                    ')' or ']' or '⟩' => new Sym(str, SymKind.Close),
                    '\'' => new Scripts(new Sym("", SymKind.Upright), new Sym("′", SymKind.Upright), null),
                    _ when char.IsLetter(c) && IsItalicLetter(c) => new Sym(str, SymKind.Italic),
                    _ => new Sym(str, SymKind.Upright),
                };
            }

            private Node Command(bool inText)
            {
                ++I;   // the backslash
                if (!More) throw new FormatException("lone '\\'");
                string name;
                if (char.IsLetter(s[I]))
                {
                    var j = I;
                    while (j < s.Length && char.IsLetter(s[j])) ++j;
                    name = s[I..j];
                    I = j;
                    if (s.Length > I && s[I] == '*' && name == "operatorname") ++I;
                }
                else name = s[I++].ToString();
                if (name == "\\") return new LineBreak();
                if (Spaces.TryGetValue(name, out var em)) return new Space(em, name is "quad" or "qquad");
                if (Greek.TryGetValue(name, out var g))
                    return new Sym(g, char.IsUpper(g[0]) && name != "ell" ? SymKind.Upright : SymKind.Italic);
                if (Symbols.TryGetValue(name, out var sy)) return new Sym(sy.Item1, sy.Item2);
                if (Funcs.Contains(name)) return new Sym(name, SymKind.Func);
                if (Limits.Contains(name)) return new LimFunc(name == "argmax" ? "arg max" : name == "argmin" ? "arg min" : name);
                switch (name)
                {
                    case "frac": case "dfrac": return new Frac(Arg(), Arg(), false);
                    case "tfrac": return new Frac(Arg(), Arg(), true);
                    case "sqrt": return new Sqrt(Arg());
                    case "text": case "mathrm": case "textrm": case "mathsf":
                        {
                            var n = Arg(inText: name == "text");
                            return name == "text" ? n : Upright(n);
                        }
                    case "operatorname": return new Sym(Flatten(Arg(inText: true)), SymKind.Func);
                    case "sum": return new BigOp("∑", false);
                    case "prod": return new BigOp("∏", false);
                    case "int": return new BigOp("∫", true);
                    case "oint": return new BigOp("∮", true);
                    case "dot": case "ddot": case "hat": case "bar": case "tilde": case "vec": case "overline": return new Accent(name, Arg());
                    case "left":
                        {
                            var open = Delim();
                            var body = Row(until: "\\right");
                            if (!s.AsSpan(I).StartsWith("\\right")) throw new FormatException("\\left without \\right");
                            I += 6;
                            return new Fenced(open, Delim(), body);
                        }
                }
                throw new FormatException($"unknown command \\{name}");
            }

            private string Delim()
            {
                while (More && s[I] == ' ') ++I;
                if (!More) throw new FormatException("missing delimiter");
                if (s[I] == '\\')
                {
                    ++I;
                    var j = I;
                    while (j < s.Length && char.IsLetter(s[j])) ++j;
                    if (j == I) j = I + 1;
                    var n = s[I..j];
                    I = j;
                    return n switch { "langle" => "⟨", "rangle" => "⟩", "{" => "{", "}" => "}", "|" => "‖", "lbrace" => "{", "rbrace" => "}", _ => throw new FormatException($"unknown delimiter \\{n}") };
                }
                var d = s[I++];
                return d == '.' ? "" : d.ToString();
            }

            private static Node Upright(Node n) => n switch
            {
                Sym y when y.Kind == SymKind.Italic => y with { Kind = SymKind.Upright },
                Seq q => new Seq(q.Items.Select(Upright).ToList()),
                Scripts sc => new Scripts(Upright(sc.Base), sc.Sup, sc.Sub),
                _ => n,
            };

            private static string Flatten(Node n) => n switch
            {
                Sym y => y.Text,
                Seq q => string.Concat(q.Items.Select(Flatten)),
                Space => " ",
                _ => "",
            };
        }

        // upright by convention: the upper-case Greek, symbols; italic: Latin letters and lower-case Greek
        private static bool IsItalicLetter(char c) => c is >= 'a' and <= 'z' or >= 'A' and <= 'Z' || c is >= 'α' and <= 'ω';
    }

    // ================================================================ boxes
    private abstract class Box
    {
        public double W, A, D;
        public abstract void Draw(DrawingContext ctx, double x, double y, IBrush brush);
    }

    private sealed class GlyphBox : Box
    {
        private readonly FormattedText _ft;
        public GlyphBox(FormattedText ft, double w, double a, double d) { _ft = ft; W = w; A = a; D = d; }
        public override void Draw(DrawingContext ctx, double x, double y, IBrush brush) => ctx.DrawText(_ft, new Point(x, y - _ft.Baseline));
    }

    // an ink-measured glyph (big operators, delimiters) placed with its ink centred where asked
    private sealed class InkBox : Box
    {
        private readonly FormattedText _ft;
        private readonly double _dy;
        public InkBox(FormattedText ft, double w, double a, double d, double dy) { _ft = ft; W = w; A = a; D = d; _dy = dy; }
        public override void Draw(DrawingContext ctx, double x, double y, IBrush brush) => ctx.DrawText(_ft, new Point(x, y + _dy));
    }

    private sealed class HBox : Box
    {
        public readonly List<(Box b, double dx)> Items = new();
        public void Add(Box b)
        {
            Items.Add((b, W));
            W += b.W;
            A = Math.Max(A, b.A);
            D = Math.Max(D, b.D);
        }
        public override void Draw(DrawingContext ctx, double x, double y, IBrush brush)
        {
            foreach (var (b, dx) in Items) b.Draw(ctx, x + dx, y, brush);
        }
    }

    private sealed class Placed : Box   // children at offsets (dx, dy of their baselines)
    {
        public readonly List<(Box b, double dx, double dy)> Items = new();
        public readonly List<Action<DrawingContext, double, double, IBrush>> Paint = new();
        public override void Draw(DrawingContext ctx, double x, double y, IBrush brush)
        {
            foreach (var (b, dx, dy) in Items) b.Draw(ctx, x + dx, y + dy, brush);
            foreach (var p in Paint) p(ctx, x, y, brush);
        }
    }

    private sealed class SpaceBox : Box
    {
        public readonly bool Breakable;
        public SpaceBox(double w, bool breakable) { W = w; Breakable = breakable; }
        public override void Draw(DrawingContext ctx, double x, double y, IBrush brush) { }
    }

    // ================================================================ layout
    private sealed class Layout(FontFamily fam, double baseSize, IBrush brush)
    {
        private readonly Typeface _up = new(fam);
        private readonly Typeface _it = new(fam, FontStyle.Italic);

        private FormattedText Text(string t, double size, bool italic) =>
            new(t, CultureInfo.InvariantCulture, FlowDirection.LeftToRight, italic ? _it : _up, size, brush);

        private static double Axis(double s) => 0.25 * s;

        // nominal heights per character, so baselines and fraction gaps do not depend on the font's line box
        private static (double a, double d) Metrics(string t, double s)
        {
            double a = 0, d = 0;
            foreach (var c in t)
            {
                var (ca, cd) = c switch
                {
                    _ when "acemnorsuvwxzıαεικνοπστυω".Contains(c) => (0.47, 0.02),
                    _ when "gpqy".Contains(c) => (0.47, 0.22),
                    _ when "γημρφχψ".Contains(c) => (0.47, 0.22),
                    _ when "βζξ".Contains(c) => (0.72, 0.22),
                    'j' or 'f' => (0.72, 0.22),
                    '(' or ')' or '[' or ']' or '{' or '}' or '|' or '‖' or '⟨' or '⟩' => (0.75, 0.25),
                    '+' or '−' or '=' or '<' or '>' or '≤' or '≥' or '≠' or '≈' or '±' or '×' or '·' or '→' or '←' or '∼' or '∝' or '∈' => (0.55, 0.05),
                    ',' or ';' => (0.1, 0.18),
                    '.' or '…' or '⋯' => (0.1, 0.0),
                    '∘' or '′' => (0.6, 0.0),
                    ' ' => (0, 0),
                    _ => (0.70, 0.02),
                };
                a = Math.Max(a, ca);
                d = Math.Max(d, cd);
            }
            return (a * s, d * s);
        }

        private Box Glyph(string t, double s, bool italic, double padL = 0, double padR = 0)
        {
            if (t == "∘") s *= 1.6;   // the ring operator is drawn small in text fonts: P° needs a visible ring
            var ft = Text(t, s, italic);
            var (a, d) = Metrics(t, s);
            var w = ft.WidthIncludingTrailingWhitespace + padL + padR;
            var g = new GlyphBox(ft, w, a, d);
            if (padL == 0) return g;
            var p = new Placed { W = w, A = a, D = d };
            p.Items.Add((new GlyphBox(ft, ft.WidthIncludingTrailingWhitespace, a, d), padL, 0));
            return p;
        }

        public List<(Box, double)> Lines(Node root, double maxW)
        {
            var items = root is Seq q ? q.Items : [root];
            var boxes = Spaced(items, baseSize, 0, true);
            var lines = new List<(Box, double)>();
            var cur = new List<Box>();
            double indent = 0;
            void Flush(double nextIndent)
            {
                while (cur.Count > 0 && cur[^1] is SpaceBox) cur.RemoveAt(cur.Count - 1);
                var h = new HBox();
                foreach (var b in cur) h.Add(b);
                if (h.A == 0 && h.D == 0) { h.A = 0.72 * baseSize; h.D = 0.22 * baseSize; }
                lines.Add((h, indent));
                cur = new List<Box>();
                indent = nextIndent;
            }
            foreach (var b in boxes)
            {
                if (b == null) { Flush(0); continue; }   // \\
                if (cur.Count == 0 && b is SpaceBox) continue;
                cur.Add(b);
                if (indent + cur.Sum(x => x.W) > maxW)
                {
                    // wrap at the last breakable space of the line
                    var k = cur.FindLastIndex(x => x is SpaceBox { Breakable: true });
                    if (k > 0)
                    {
                        var rest = cur.Skip(k + 1).ToList();
                        cur = cur.Take(k).ToList();
                        Flush(1.5 * baseSize);
                        cur = rest;
                    }
                }
            }
            Flush(0);
            return lines;
        }

        // a row's boxes with TeX's spacing between operators, relations and punctuation (null marks \\)
        private List<Box?> Spaced(List<Node> items, double s, int level, bool top)
        {
            var outp = new List<Box?>();
            Node? prev = null;
            var script = level > 0;
            foreach (var n in items)
            {
                if (n is LineBreak) { outp.Add(null); prev = null; continue; }
                var kind = KindOf(n);
                var prevKind = prev == null ? (SymKind?)null : KindOf(prev);
                // a binary operator after nothing, an opening bracket, an operator or a relation is a sign
                var unary = kind == SymKind.Bin && (prevKind is null or SymKind.Bin or SymKind.Rel or SymKind.Open or SymKind.Punct);
                if (!script)
                {
                    if (kind == SymKind.Rel && prevKind is not null and not SymKind.Open) outp.Add(new SpaceBox(0.28 * s, false));
                    else if (kind == SymKind.Bin && !unary) outp.Add(new SpaceBox(0.22 * s, false));
                    else if (prevKind == SymKind.Func && (kind is SymKind.Italic or SymKind.Upright or SymKind.Number || n is Sym { Text: "⟨" })) outp.Add(new SpaceBox(0.17 * s, false));
                    else if (prev is LimFunc or BigOp or Scripts { Base: BigOp or LimFunc } && n is not Space) outp.Add(new SpaceBox(0.12 * s, false));
                }
                outp.Add(Build(n, s, level, top));
                if (!script)
                {
                    if (kind == SymKind.Rel) outp.Add(new SpaceBox(0.28 * s, false));
                    else if (kind == SymKind.Bin && !unary) outp.Add(new SpaceBox(0.22 * s, false));
                    else if (kind == SymKind.Punct) outp.Add(new SpaceBox(0.17 * s, false));
                }
                else if (kind == SymKind.Punct) outp.Add(new SpaceBox(0.1 * s, false));
                prev = n is Space ? prev : n;
            }
            return outp;
        }

        private static SymKind KindOf(Node n) => n switch
        {
            Sym y => y.Kind,
            Scripts sc => KindOf(sc.Base) is SymKind.Bin or SymKind.Rel ? SymKind.Upright : KindOf(sc.Base),
            Fenced => SymKind.Close,
            _ => SymKind.Upright,
        };

        private Box Row(Node n, double s, int level)
        {
            var items = n is Seq q ? q.Items : [n];
            var h = new HBox();
            foreach (var b in Spaced(items, s, level, false)) if (b != null) h.Add(b);
            return h;
        }

        private Box Build(Node n, double s, int level, bool top)
        {
            switch (n)
            {
                case Sym y:
                    return y.Kind switch
                    {
                        SymKind.Italic => Glyph(y.Text, s, true),
                        SymKind.Text => Glyph(y.Text, s, false),
                        _ => Glyph(y.Text, s, false),
                    };
                case Space sp: return new SpaceBox(sp.Em * s, sp.Breakable && top);
                case Seq q: return Row(q, s, level);
                case Frac f: return BuildFrac(f, s, level);
                case Sqrt r: return BuildSqrt(r, s, level);
                case Scripts sc: return BuildScripts(sc, s, level, top);
                case BigOp op: return BuildBigOp(op, s, level);
                case LimFunc lf: return Glyph(lf.Name, s, false);
                case Fenced fe: return BuildFenced(fe, s, level);
                case Accent ac: return BuildAccent(ac, s, level);
                default: return new SpaceBox(0, false);
            }
        }

        private Box BuildFrac(Frac f, double s, int level)
        {
            // display fractions keep the size at the top level; nested and \tfrac ones step down
            var cs = f.Small ? s * 0.72 : level == 0 ? s : s * 0.85;
            var num = Row(f.Num, cs, f.Small ? level + 1 : level);
            var den = Row(f.Den, cs, f.Small ? level + 1 : level);
            var t = Math.Max(1, 0.05 * s);
            var gap = (f.Small ? 0.1 : 0.14) * s;
            var axis = Axis(s);
            var up = axis + t / 2 + gap + num.D;
            var down = -axis + t / 2 + gap + den.A;
            var pad = 0.1 * s;
            var w = Math.Max(num.W, den.W) + 2 * pad;
            var p = new Placed { W = w + 0.1 * s, A = up + num.A, D = down + den.D };
            p.Items.Add((num, 0.05 * s + (w - num.W) / 2, -up));
            p.Items.Add((den, 0.05 * s + (w - den.W) / 2, down));
            p.Paint.Add((ctx, x, y, b) => ctx.FillRectangle(b, new Rect(x + 0.05 * s, y - axis - t / 2, w, t)));
            return p;
        }

        private Box BuildSqrt(Sqrt r, double s, int level)
        {
            var body = Row(r.Body, s, level);
            var t = Math.Max(1, 0.05 * s);
            var top = body.A + 0.14 * s;
            var bottom = body.D + 0.02 * s;
            var sign = 0.55 * s;
            var p = new Placed { W = sign + body.W + 0.12 * s, A = top + t, D = bottom };
            p.Items.Add((body, sign + 0.04 * s, 0));
            p.Paint.Add((ctx, x, y, b) =>
            {
                var pen = new Pen(b, t, lineCap: PenLineCap.Round, lineJoin: PenLineJoin.Round);
                var g = new StreamGeometry();
                using (var c = g.Open())
                {
                    var h = top + bottom;
                    c.BeginFigure(new Point(x + 0.04 * s, y - 0.45 * h + bottom + 0.05 * s), false);
                    c.LineTo(new Point(x + 0.15 * s, y - 0.5 * h + bottom));
                    c.LineTo(new Point(x + 0.3 * s, y + bottom));
                    c.LineTo(new Point(x + sign, y - top));
                    c.LineTo(new Point(x + sign + body.W + 0.12 * s, y - top));
                    c.EndFigure(false);
                }
                ctx.DrawGeometry(null, pen, g);
            });
            return p;
        }

        private Box BuildBigOp(BigOp op, double s, int level)
        {
            var size = level == 0 ? (op.Integral ? 1.9 : 1.6) * s : 1.15 * s;
            var ft = Text(op.Glyph, size, false);
            var ink = ft.BuildGeometry(new Point(0, 0))?.Bounds ?? new Rect(0, 0, ft.Width, ft.Height);
            // centre the ink on the maths axis
            var half = ink.Height / 2;
            var axis = Axis(s);
            var dy = -axis - half - ink.Top;   // text origin (top of line box) relative to the baseline
            return new InkBox(ft, ft.WidthIncludingTrailingWhitespace + 0.05 * s, axis + half, half - axis, dy);
        }

        private Box BuildScripts(Scripts sc, double s, int level, bool top)
        {
            var ss = level == 0 ? s * 0.7 : s * 0.8;   // scripts of scripts: about 0.56 of the text size
            var baseBox = Build(sc.Base, s, level, false);
            var sup = sc.Sup == null ? null : Row(sc.Sup, ss, level + 1);
            var sub = sc.Sub == null ? null : Row(sc.Sub, ss, level + 1);
            // limits above and below: sums and products, lim/min/max at the top level (display)
            if (level == 0 && (sc.Base is BigOp { Integral: false } || sc.Base is LimFunc))
            {
                var w = Math.Max(baseBox.W, Math.Max(sup?.W ?? 0, sub?.W ?? 0));
                var p = new Placed { W = w, A = baseBox.A, D = baseBox.D };
                p.Items.Add((baseBox, (w - baseBox.W) / 2, 0));
                if (sup != null)
                {
                    var dy = baseBox.A + 0.12 * s + sup.D;
                    p.Items.Add((sup, (w - sup.W) / 2, -dy));
                    p.A = dy + sup.A;
                }
                if (sub != null)
                {
                    var dy = baseBox.D + 0.1 * s + sub.A;
                    p.Items.Add((sub, (w - sub.W) / 2, dy));
                    p.D = dy + sub.D;
                }
                return p;
            }
            var italic = sc.Base is Sym { Kind: SymKind.Italic };
            var isInt = sc.Base is BigOp { Integral: true };
            var supUp = isInt ? baseBox.A - 0.3 * s : Math.Max(0.38 * s, baseBox.A - 0.28 * s);
            var subDown = isInt ? baseBox.D - 0.05 * s : Math.Max(0.17 * s, baseBox.D + 0.02 * s - (baseBox.D > 0.1 * s ? 0.12 * s : 0));
            if (sup != null && sub != null)
            {
                // keep a gap between them
                var gap = (supUp - sup.D) - (sub.A - subDown);
                if (gap < 0.16 * s) { var fix = (0.16 * s - gap) / 2; supUp += fix; subDown += fix; }
            }
            var q = new Placed { W = baseBox.W, A = baseBox.A, D = baseBox.D };
            q.Items.Add((baseBox, 0, 0));
            double w2 = 0;
            if (sup != null)
            {
                var dx = baseBox.W + (italic ? 0.06 * s : 0.02 * s) + (isInt ? -0.1 * s : 0);
                q.Items.Add((sup, dx, -supUp));
                q.A = Math.Max(q.A, supUp + sup.A);
                w2 = Math.Max(w2, dx + sup.W);
            }
            if (sub != null)
            {
                var dx = baseBox.W + (italic ? -0.02 * s : 0.01 * s) + (isInt ? -0.35 * s : 0);
                q.Items.Add((sub, dx, subDown));
                q.D = Math.Max(q.D, subDown + sub.D);
                w2 = Math.Max(w2, dx + sub.W);
            }
            q.W = Math.Max(baseBox.W, w2) + 0.04 * s;
            return q;
        }

        private Box BuildFenced(Fenced fe, double s, int level)
        {
            var body = Row(fe.Body, s, level);
            var axis = Axis(s);
            var half = Math.Max(body.A - axis, body.D + axis) + 0.08 * s;
            var small = half <= 0.52 * s;
            var p = new Placed { A = body.A, D = body.D };
            double x = 0;
            Box? Side(string d)
            {
                if (d.Length == 0) return null;
                if (small) return Glyph(d, s, false);
                var h = 2 * half;
                var w = d switch { "(" or ")" => 0.3 * s + 0.05 * h, "[" or "]" => 0.32 * s, "{" or "}" => 0.42 * s, "⟨" or "⟩" => 0.3 * s + 0.08 * h, _ => 0.25 * s };
                var box = new Placed { W = w, A = axis + half, D = half - axis };
                box.Paint.Add((ctx, bx, by, b) => DrawDelim(ctx, d, bx, by - axis - half, w, h, s, b));
                return box;
            }
            var l = Side(fe.Open);
            var r = Side(fe.Close);
            if (l != null) { p.Items.Add((l, 0, 0)); x += l.W; p.A = Math.Max(p.A, l.A); p.D = Math.Max(p.D, l.D); }
            p.Items.Add((body, x + 0.03 * s, 0));
            x += body.W + 0.06 * s;
            if (r != null) { p.Items.Add((r, x, 0)); x += r.W; p.A = Math.Max(p.A, r.A); p.D = Math.Max(p.D, r.D); }
            p.W = x;
            return p;
        }

        private static void DrawDelim(DrawingContext ctx, string d, double x, double top, double w, double h, double s, IBrush b)
        {
            var t = Math.Max(1, 0.055 * s);
            var pen = new Pen(b, t, lineCap: PenLineCap.Round, lineJoin: PenLineJoin.Round);
            var g = new StreamGeometry();
            var m = 0.1 * s;   // side margin
            using (var c = g.Open())
            {
                var bot = top + h;
                var mid = top + h / 2;
                switch (d)
                {
                    case "(":
                        c.BeginFigure(new Point(x + w - m, top), false);
                        c.CubicBezierTo(new Point(x + m - 0.02 * h, top + 0.25 * h), new Point(x + m - 0.02 * h, bot - 0.25 * h), new Point(x + w - m, bot));
                        break;
                    case ")":
                        c.BeginFigure(new Point(x + m, top), false);
                        c.CubicBezierTo(new Point(x + w - m + 0.02 * h, top + 0.25 * h), new Point(x + w - m + 0.02 * h, bot - 0.25 * h), new Point(x + m, bot));
                        break;
                    case "[":
                        c.BeginFigure(new Point(x + w - m, top), false); c.LineTo(new Point(x + m, top)); c.LineTo(new Point(x + m, bot)); c.LineTo(new Point(x + w - m, bot));
                        break;
                    case "]":
                        c.BeginFigure(new Point(x + m, top), false); c.LineTo(new Point(x + w - m, top)); c.LineTo(new Point(x + w - m, bot)); c.LineTo(new Point(x + m, bot));
                        break;
                    case "⟨":
                        c.BeginFigure(new Point(x + w - m, top), false); c.LineTo(new Point(x + m, mid)); c.LineTo(new Point(x + w - m, bot));
                        break;
                    case "⟩":
                        c.BeginFigure(new Point(x + m, top), false); c.LineTo(new Point(x + w - m, mid)); c.LineTo(new Point(x + m, bot));
                        break;
                    case "{":
                    case "}":
                        {
                            var xo = d == "{" ? x + w - m : x + m;
                            var xi = d == "{" ? x + m : x + w - m;
                            var xc = (xo + xi) / 2;
                            c.BeginFigure(new Point(xo, top), false);
                            c.CubicBezierTo(new Point(xc, top), new Point(xc, top), new Point(xc, top + 0.2 * h));
                            c.LineTo(new Point(xc, mid - 0.12 * h));
                            c.CubicBezierTo(new Point(xc, mid - 0.02 * h), new Point(xc, mid), new Point(xi, mid));
                            c.CubicBezierTo(new Point(xc, mid), new Point(xc, mid + 0.02 * h), new Point(xc, mid + 0.12 * h));
                            c.LineTo(new Point(xc, bot - 0.2 * h));
                            c.CubicBezierTo(new Point(xc, bot), new Point(xc, bot), new Point(xo, bot));
                            break;
                        }
                    case "‖":
                        c.BeginFigure(new Point(x + w / 2 - 0.06 * s, top), false); c.LineTo(new Point(x + w / 2 - 0.06 * s, bot)); c.EndFigure(false);
                        c.BeginFigure(new Point(x + w / 2 + 0.06 * s, top), false); c.LineTo(new Point(x + w / 2 + 0.06 * s, bot));
                        break;
                    default:   // |
                        c.BeginFigure(new Point(x + w / 2, top), false); c.LineTo(new Point(x + w / 2, bot));
                        break;
                }
                c.EndFigure(false);
            }
            ctx.DrawGeometry(null, pen, g);
        }

        private Box BuildAccent(Accent ac, double s, int level)
        {
            var body = Build(ac.Body, s, level, false);
            var italic = ac.Body is Sym { Kind: SymKind.Italic };
            var shift = italic ? 0.08 * s : 0;   // italic letters lean right: the accent follows
            var gap = 0.08 * s;
            var ah = ac.Kind switch { "dot" or "ddot" => 0.1 * s, "bar" or "overline" => 0.06 * s, _ => 0.16 * s };
            var p = new Placed { W = body.W, A = body.A + gap + ah, D = body.D };
            p.Items.Add((body, 0, 0));
            p.Paint.Add((ctx, x, y, b) =>
            {
                var cx = x + body.W / 2 + shift;
                var yb = y - body.A - gap;   // bottom of the accent
                var t = Math.Max(1, 0.05 * s);
                var pen = new Pen(b, t, lineCap: PenLineCap.Round, lineJoin: PenLineJoin.Round);
                switch (ac.Kind)
                {
                    case "dot":
                        ctx.DrawEllipse(b, null, new Point(cx, yb - 0.05 * s), 0.05 * s, 0.05 * s);
                        break;
                    case "ddot":
                        ctx.DrawEllipse(b, null, new Point(cx - 0.11 * s, yb - 0.05 * s), 0.05 * s, 0.05 * s);
                        ctx.DrawEllipse(b, null, new Point(cx + 0.11 * s, yb - 0.05 * s), 0.05 * s, 0.05 * s);
                        break;
                    case "bar":
                    case "overline":
                        {
                            var hw = ac.Kind == "overline" ? body.W / 2 : Math.Max(0.2 * s, body.W / 2 - 0.06 * s);
                            ctx.DrawLine(pen, new Point(cx - hw, yb - 0.03 * s), new Point(cx + hw, yb - 0.03 * s));
                            break;
                        }
                    case "vec":
                        ctx.DrawLine(pen, new Point(cx - 0.22 * s, yb - 0.08 * s), new Point(cx + 0.22 * s, yb - 0.08 * s));
                        ctx.DrawLine(pen, new Point(cx + 0.12 * s, yb - 0.15 * s), new Point(cx + 0.22 * s, yb - 0.08 * s));
                        ctx.DrawLine(pen, new Point(cx + 0.12 * s, yb - 0.01 * s), new Point(cx + 0.22 * s, yb - 0.08 * s));
                        break;
                    case "tilde":
                        {
                            var g = new StreamGeometry();
                            using (var c = g.Open())
                            {
                                c.BeginFigure(new Point(cx - 0.2 * s, yb - 0.03 * s), false);
                                c.CubicBezierTo(new Point(cx - 0.1 * s, yb - 0.2 * s), new Point(cx + 0.1 * s, yb + 0.06 * s), new Point(cx + 0.2 * s, yb - 0.12 * s));
                                c.EndFigure(false);
                            }
                            ctx.DrawGeometry(null, pen, g);
                            break;
                        }
                    default:   // hat
                        ctx.DrawLine(pen, new Point(cx - 0.18 * s, yb), new Point(cx, yb - 0.15 * s));
                        ctx.DrawLine(pen, new Point(cx, yb - 0.15 * s), new Point(cx + 0.18 * s, yb));
                        break;
                }
            });
            return p;
        }
    }
}
