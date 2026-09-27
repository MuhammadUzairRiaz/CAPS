using System.Runtime.InteropServices;
using Avalonia.OpenGL;
using Avalonia.OpenGL.Controls;
using CapsStudio.Interop;

namespace CapsStudio.Views;

/// <summary>The GPU 3D view: the scene of caps_render_scene (spheres, half-bond capsules, lines) drawn as ray-cast
/// impostors — per pixel the same sphere and capsule surfaces, Blinn–Phong shading, depth cue, edge ink and selection
/// rings as the CPU renderer — on the camera caps_view_fit gives, so it lines up with the CPU view to the pixel. The
/// scene is uploaded once; turning the view only changes uniforms, so large systems rotate at the display's rate.
/// OpenGL 3.3 core or OpenGL ES 3.0 (ANGLE on Windows); where neither is available Ready stays false and the CPU view
/// is used.</summary>
public sealed unsafe class GlMolView : OpenGlControlBase
{
    private const int GL_ARRAY_BUFFER = 0x8892, GL_STATIC_DRAW = 0x88E4, GL_FLOAT = 0x1406, GL_UNSIGNED_BYTE = 0x1401;
    private const int GL_TRIANGLE_STRIP = 0x0005, GL_DEPTH_TEST = 0x0B71, GL_LESS = 0x0201, GL_BLEND = 0x0BE2, GL_CULL_FACE = 0x0B44;
    private const int GL_COLOR_BUFFER_BIT = 0x4000, GL_DEPTH_BUFFER_BIT = 0x0100, GL_FRAMEBUFFER = 0x8D40;
    private const int GL_VERTEX_SHADER = 0x8B31, GL_FRAGMENT_SHADER = 0x8B30, GL_SCISSOR_TEST = 0x0C11;

    private delegate* unmanaged<int, int, void> _divisor;
    private delegate* unmanaged<int, int, int, int*, void> _attachParam;
    private const int GL_DEPTH_ATTACHMENT = 0x8D00;
    /// <summary>The size of the framebuffer drawn last (pixels).</summary>
    public (int W, int H) FramebufferSize { get; private set; }
    private delegate* unmanaged<int, int, int, int, void> _drawInstanced;
    private delegate* unmanaged<int, float, float, void> _u2f;
    private delegate* unmanaged<int, float, float, float, void> _u3f;
    private delegate* unmanaged<int, float, float, float, float, void> _u4f;

    private int _progSphere, _progCapsule, _progLine;
    private int _vaoSphere, _vaoCapsule, _vaoLine;
    private int _cornerSphere, _cornerCapsule;
    private int _bufSphere, _bufSphereRgb, _bufSphereRing, _bufCapsule, _bufCapsuleRgb, _bufLine, _bufLineRgb;
    private int _nSphere, _nCapsule, _nLine;
    private CapsSceneData? _pending, _scene;
    private CapsViewFit _fit;
    private bool _haveFit;
    private (double X, double Y, double Z, double R) _bounds;

    /// <summary>OpenGL came up and the shaders compiled: the view can draw.</summary>
    public bool Ready { get; private set; }
    /// <summary>Why the GPU view is not in use, when it is not (shown in Settings).</summary>
    public string Status { get; private set; } = "starting";
    public event Action? ReadyChanged;
    /// <summary>The framebuffer no longer has the size the camera was fitted for.</summary>
    public event Action? FitStale;
    /// <summary>Milliseconds the last frame took to submit (CPU side).</summary>
    public double LastFrameMs { get; private set; }

    /// <summary>A new scene (content changed): uploaded at the next frame.</summary>
    public void SetScene(CapsSceneData scene)
    {
        _pending = scene;
        RequestNextFrameRendering();
    }

    /// <summary>The camera for the next frame (turning the view: nothing is uploaded).</summary>
    public void SetView(CapsViewFit fit)
    {
        _fit = fit;
        _haveFit = true;
        RequestNextFrameRendering();
    }

    protected override void OnOpenGlInit(GlInterface gl)
    {
        base.OnOpenGlInit(gl);
        try
        {
            var v = GlVersion;
            var es = v.Type == GlProfileType.OpenGLES;
            if (es ? v.Major < 3 : (v.Major < 3 || (v.Major == 3 && v.Minor < 3)))
                throw new NotSupportedException($"OpenGL {(es ? "ES " : "")}{v.Major}.{v.Minor} (needs 3.3, or ES 3.0)");
            IntPtr P(string name)
            {
                var p = gl.GetProcAddress(name);
                if (p == IntPtr.Zero) p = gl.GetProcAddress(name + "ARB");
                if (p == IntPtr.Zero) throw new NotSupportedException(name + " missing");
                return p;
            }
            _divisor = (delegate* unmanaged<int, int, void>)P("glVertexAttribDivisor");
            _drawInstanced = (delegate* unmanaged<int, int, int, int, void>)P("glDrawArraysInstanced");
            _u2f = (delegate* unmanaged<int, float, float, void>)P("glUniform2f");
            _u3f = (delegate* unmanaged<int, float, float, float, void>)P("glUniform3f");
            _u4f = (delegate* unmanaged<int, float, float, float, float, void>)P("glUniform4f");
            var ap = gl.GetProcAddress("glGetFramebufferAttachmentParameteriv");
            if (ap != IntPtr.Zero) _attachParam = (delegate* unmanaged<int, int, int, int*, void>)ap;
            var head = es ? "#version 300 es\nprecision highp float;\nprecision highp int;\n" : "#version 330 core\n";
            _progSphere = Program(gl, head + Common + SphereVs, head + Common + Shade + SphereFs);
            _progCapsule = Program(gl, head + Common + CapsuleVs, head + Common + Shade + CapsuleFs);
            _progLine = Program(gl, head + Common + LineVs, head + Common + LineFs);
            _cornerSphere = Buffer(gl, [-1f, -1f, 1f, -1f, -1f, 1f, 1f, 1f]);
            _cornerCapsule = Buffer(gl, [0f, -1f, 1f, -1f, 0f, 1f, 1f, 1f]);
            _vaoSphere = gl.GenVertexArray();
            _vaoCapsule = gl.GenVertexArray();
            _vaoLine = gl.GenVertexArray();
            _bufSphere = gl.GenBuffer(); _bufSphereRgb = gl.GenBuffer(); _bufSphereRing = gl.GenBuffer();
            _bufCapsule = gl.GenBuffer(); _bufCapsuleRgb = gl.GenBuffer();
            _bufLine = gl.GenBuffer(); _bufLineRgb = gl.GenBuffer();
            Ready = true;
            Status = $"OpenGL {(es ? "ES " : "")}{v.Major}.{v.Minor} · {gl.Renderer}";
        }
        catch (Exception e)
        {
            Ready = false;
            Status = "not available: " + e.Message;
        }
        ReadyChanged?.Invoke();
    }

    protected override void OnOpenGlDeinit(GlInterface gl)
    {
        foreach (var p in new[] { _progSphere, _progCapsule, _progLine }) if (p != 0) gl.DeleteProgram(p);
        foreach (var b in new[] { _cornerSphere, _cornerCapsule, _bufSphere, _bufSphereRgb, _bufSphereRing, _bufCapsule, _bufCapsuleRgb, _bufLine, _bufLineRgb })
            if (b != 0) gl.DeleteBuffer(b);
        foreach (var a in new[] { _vaoSphere, _vaoCapsule, _vaoLine }) if (a != 0) gl.DeleteVertexArray(a);
        _progSphere = _progCapsule = _progLine = 0;
        Ready = false;
        base.OnOpenGlDeinit(gl);
    }

    protected override void OnOpenGlLost()
    {
        Ready = false;
        Status = "context lost";
        ReadyChanged?.Invoke();
        base.OnOpenGlLost();
    }

    protected override void OnOpenGlRender(GlInterface gl, int fb)
    {
        if (!Ready) return;
        var sw = System.Diagnostics.Stopwatch.StartNew();
        if (_pending != null) { Upload(gl, _pending); _scene = _pending; _pending = null; }
        var scale = VisualRoot?.RenderScaling ?? 1;
        var pw = Math.Max(1, (int)(Bounds.Width * scale));
        var ph = Math.Max(1, (int)(Bounds.Height * scale));
        gl.BindFramebuffer(GL_FRAMEBUFFER, fb);
        // the framebuffer's real size, from the depth renderbuffer Avalonia attaches to it
        if (_attachParam != null)
        {
            int name = 0, type = 0;
            _attachParam(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, 0x8CD0, &type);   // GL_FRAMEBUFFER_ATTACHMENT_OBJECT_TYPE
            if (type == 0x8D41)   // GL_RENDERBUFFER
            {
                _attachParam(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, 0x8CD1, &name);   // GL_FRAMEBUFFER_ATTACHMENT_OBJECT_NAME
                if (name != 0)
                {
                    gl.BindRenderbuffer(0x8D41, name);
                    gl.GetRenderbufferParameteriv(0x8D41, 0x8D42, out var rw);   // GL_RENDERBUFFER_WIDTH
                    gl.GetRenderbufferParameteriv(0x8D41, 0x8D43, out var rh);   // GL_RENDERBUFFER_HEIGHT
                    gl.BindRenderbuffer(0x8D41, 0);
                    if (rw > 0 && rh > 0) { pw = rw; ph = rh; }
                }
            }
        }
        FramebufferSize = (pw, ph);
        // standard-DPI screens: draw at twice the size into our own framebuffer and filter it down (the CPU view's
        // 2 × 2 supersampling); Retina screens already have two samples per point
        var ss = (scale < 1.5 || Environment.GetEnvironmentVariable("CAPS_GL_SS") == "2") && gl.IsBlitFramebufferAvailable ? 2 : 1;   // CAPS_GL_SS=2: also on Retina
        var target = fb;
        if (ss > 1 && EnsureSupersample(gl, pw * ss, ph * ss)) target = _ssFbo;
        else ss = 1;
        var outW = pw; var outH = ph;
        pw *= ss; ph *= ss;
        gl.BindFramebuffer(GL_FRAMEBUFFER, target);
        gl.Viewport(0, 0, pw, ph);
        gl.Disable(GL_SCISSOR_TEST);
        var sc = _scene;
        var bg = sc?.Background ?? 0x0F1113u;
        if (sc?.Transparent == true) gl.ClearColor(0, 0, 0, 0);
        else gl.ClearColor(((bg >> 16) & 255) / 255f, ((bg >> 8) & 255) / 255f, (bg & 255) / 255f, 1);
        gl.ClearDepth(1);
        gl.DepthMask(1);
        gl.Clear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        if (sc == null || !_haveFit) { Resolve(gl, target, fb, pw, ph, outW, outH); LastFrameMs = sw.Elapsed.TotalMilliseconds; return; }
        gl.Enable(GL_DEPTH_TEST);
        gl.DepthFunc(GL_LESS);
        gl.Disable(GL_BLEND);
        gl.Disable(GL_CULL_FACE);
        var f = _fit;
        // the fit was made for this framebuffer's size; if it lags a resize, scale it across and ask for a new one
        var sx = (float)Math.Min(pw / Math.Max(1, f.W), ph / Math.Max(1, f.H));
        if (Math.Abs(f.W - outW) > 1 || Math.Abs(f.H - outH) > 1) Avalonia.Threading.Dispatcher.UIThread.Post(() => FitStale?.Invoke());
        // the depth range covers the whole scene in its current orientation
        var cz = RotZ(f, _bounds.X, _bounds.Y, _bounds.Z);
        var pad = _bounds.R + (sc.MaxRadius * 2 + 2);
        var zlo = (float)(cz - pad);
        var zhi = (float)(cz + pad);
        void Common(int prog)
        {
            gl.UseProgram(prog);
            _u3f(gl.GetUniformLocationString(prog, "uCentre"), (float)f.Cx, (float)f.Cy, (float)f.Cz);
            _u4f(gl.GetUniformLocationString(prog, "uRot"), (float)f.CosYaw, (float)f.SinYaw, (float)f.CosPitch, (float)f.SinPitch);
            _u2f(gl.GetUniformLocationString(prog, "uPan"), (float)f.PanX, (float)f.PanY);
            gl.Uniform1f(gl.GetUniformLocationString(prog, "uScale"), (float)f.Scale * sx);
            _u2f(gl.GetUniformLocationString(prog, "uSize"), pw, ph);
            gl.Uniform1f(gl.GetUniformLocationString(prog, "uPersp"), f.Perspective != 0 ? 1f : 0f);
            gl.Uniform1f(gl.GetUniformLocationString(prog, "uDist"), (float)f.Dist);
            _u2f(gl.GetUniformLocationString(prog, "uDepth"), zlo, zhi);
            _u2f(gl.GetUniformLocationString(prog, "uZ"), (float)f.ZMin, (float)f.ZMax);
            _u3f(gl.GetUniformLocationString(prog, "uBg"), ((bg >> 16) & 255) / 255f, ((bg >> 8) & 255) / 255f, (bg & 255) / 255f);
            gl.Uniform1f(gl.GetUniformLocationString(prog, "uCue"), sc.DepthCue ? 1f : 0f);
            gl.Uniform1f(gl.GetUniformLocationString(prog, "uTransparent"), sc.Transparent ? 1f : 0f);
            gl.Uniform1f(gl.GetUniformLocationString(prog, "uOutline"), sc.Outlines ? 1f : 0f);
            gl.Uniform1f(gl.GetUniformLocationString(prog, "uPx"), ss);   // edge, ring and line widths in drawn pixels
            var ink = sc.Dark ? (0.04f, 0.045f, 0.05f) : (0.08f, 0.08f, 0.08f);
            _u3f(gl.GetUniformLocationString(prog, "uInk"), ink.Item1, ink.Item2, ink.Item3);
        }
        if (_nCapsule > 0)
        {
            Common(_progCapsule);
            gl.BindVertexArray(_vaoCapsule);
            _drawInstanced(GL_TRIANGLE_STRIP, 0, 4, _nCapsule);
        }
        if (_nSphere > 0)
        {
            Common(_progSphere);
            gl.BindVertexArray(_vaoSphere);
            _drawInstanced(GL_TRIANGLE_STRIP, 0, 4, _nSphere);
        }
        if (_nLine > 0)
        {
            Common(_progLine);
            gl.BindVertexArray(_vaoLine);
            _drawInstanced(GL_TRIANGLE_STRIP, 0, 4, _nLine);
        }
        gl.BindVertexArray(0);
        gl.UseProgram(0);
        gl.Disable(GL_DEPTH_TEST);
        Resolve(gl, target, fb, pw, ph, outW, outH);
        LastFrameMs = sw.Elapsed.TotalMilliseconds;
    }

    // ---- supersampling on standard-DPI screens
    private int _ssFbo, _ssColour, _ssDepth, _ssW, _ssH;

    private bool EnsureSupersample(GlInterface gl, int w, int h)
    {
        if (_ssFbo != 0 && _ssW == w && _ssH == h) return true;
        if (_ssFbo == 0) { _ssFbo = gl.GenFramebuffer(); _ssColour = gl.GenRenderbuffer(); _ssDepth = gl.GenRenderbuffer(); }
        gl.BindRenderbuffer(0x8D41, _ssColour);
        gl.RenderbufferStorage(0x8D41, 0x8058, w, h);   // GL_RGBA8
        gl.BindRenderbuffer(0x8D41, _ssDepth);
        gl.RenderbufferStorage(0x8D41, 0x81A6, w, h);   // GL_DEPTH_COMPONENT24
        gl.BindRenderbuffer(0x8D41, 0);
        gl.BindFramebuffer(GL_FRAMEBUFFER, _ssFbo);
        gl.FramebufferRenderbuffer(GL_FRAMEBUFFER, 0x8CE0, 0x8D41, _ssColour);   // GL_COLOR_ATTACHMENT0
        gl.FramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, 0x8D41, _ssDepth);
        var ok = gl.CheckFramebufferStatus(GL_FRAMEBUFFER) == 0x8CD5;   // GL_FRAMEBUFFER_COMPLETE
        _ssW = w; _ssH = h;
        if (!ok) { _ssW = _ssH = -1; }
        return ok;
    }

    private static void Resolve(GlInterface gl, int from, int to, int w, int h, int outW, int outH)
    {
        if (from == to) return;
        gl.BindFramebuffer(0x8CA8, from);   // GL_READ_FRAMEBUFFER
        gl.BindFramebuffer(0x8CA9, to);     // GL_DRAW_FRAMEBUFFER
        gl.BlitFramebuffer(0, 0, w, h, 0, 0, outW, outH, GL_COLOR_BUFFER_BIT, 0x2601);   // GL_LINEAR: a 2 × 2 box
        gl.BindFramebuffer(GL_FRAMEBUFFER, to);
    }

    private static double RotZ(in CapsViewFit f, double x, double y, double z)
    {
        double dx = x - f.Cx, dy = y - f.Cy, dz = z - f.Cz;
        var zr = -dx * f.SinYaw + dz * f.CosYaw;
        return dy * f.SinPitch + zr * f.CosPitch;
    }

    private void Upload(GlInterface gl, CapsSceneData d)
    {
        // a sphere around everything drawn (the depth range is taken from it each frame)
        double lo0 = double.MaxValue, lo1 = double.MaxValue, lo2 = double.MaxValue, hi0 = double.MinValue, hi1 = double.MinValue, hi2 = double.MinValue;
        void Grow(float x, float y, float z) { lo0 = Math.Min(lo0, x); lo1 = Math.Min(lo1, y); lo2 = Math.Min(lo2, z); hi0 = Math.Max(hi0, x); hi1 = Math.Max(hi1, y); hi2 = Math.Max(hi2, z); }
        for (var k = 0; k + 3 < d.Spheres.Length; k += 4) Grow(d.Spheres[k], d.Spheres[k + 1], d.Spheres[k + 2]);
        for (var k = 0; k + 6 < d.Capsules.Length; k += 7) { Grow(d.Capsules[k], d.Capsules[k + 1], d.Capsules[k + 2]); Grow(d.Capsules[k + 3], d.Capsules[k + 4], d.Capsules[k + 5]); }
        for (var k = 0; k + 5 < d.Lines.Length; k += 6) { Grow(d.Lines[k], d.Lines[k + 1], d.Lines[k + 2]); Grow(d.Lines[k + 3], d.Lines[k + 4], d.Lines[k + 5]); }
        if (lo0 > hi0) lo0 = lo1 = lo2 = hi0 = hi1 = hi2 = 0;
        _bounds = ((lo0 + hi0) / 2, (lo1 + hi1) / 2, (lo2 + hi2) / 2, 0.5 * Math.Sqrt((hi0 - lo0) * (hi0 - lo0) + (hi1 - lo1) * (hi1 - lo1) + (hi2 - lo2) * (hi2 - lo2)));

        _nSphere = d.SphereId.Length;
        _nCapsule = d.CapsuleRgb.Length;
        _nLine = d.LineRgb.Length;
        // spheres: corner (per vertex), x y z r, colour, ring (per instance)
        gl.BindVertexArray(_vaoSphere);
        Attrib(gl, _cornerSphere, 0, 2, GL_FLOAT, false, 8, 0, 0);
        Data(gl, _bufSphere, d.Spheres);
        Attrib(gl, _bufSphere, 1, 4, GL_FLOAT, false, 16, 0, 1);
        Data(gl, _bufSphereRgb, d.SphereRgb);
        Attrib(gl, _bufSphereRgb, 2, 4, GL_UNSIGNED_BYTE, true, 4, 0, 1);
        Data(gl, _bufSphereRing, d.SphereRing);
        Attrib(gl, _bufSphereRing, 3, 1, GL_UNSIGNED_BYTE, false, 1, 0, 1);
        // capsules: a, (b, r), colour
        gl.BindVertexArray(_vaoCapsule);
        Attrib(gl, _cornerCapsule, 0, 2, GL_FLOAT, false, 8, 0, 0);
        Data(gl, _bufCapsule, d.Capsules);
        Attrib(gl, _bufCapsule, 1, 3, GL_FLOAT, false, 28, 0, 1);
        Attrib(gl, _bufCapsule, 2, 4, GL_FLOAT, false, 28, 12, 1);
        Data(gl, _bufCapsuleRgb, d.CapsuleRgb);
        Attrib(gl, _bufCapsuleRgb, 3, 4, GL_UNSIGNED_BYTE, true, 4, 0, 1);
        // lines: a, (b, width), colour
        var lw = new float[7 * _nLine];
        for (var k = 0; k < _nLine; k++)
        {
            Array.Copy(d.Lines, 6 * k, lw, 7 * k, 6);
            lw[7 * k + 6] = d.LineWidth.Length > k ? d.LineWidth[k] : 1f;
        }
        gl.BindVertexArray(_vaoLine);
        Attrib(gl, _cornerCapsule, 0, 2, GL_FLOAT, false, 8, 0, 0);
        Data(gl, _bufLine, lw);
        Attrib(gl, _bufLine, 1, 3, GL_FLOAT, false, 28, 0, 1);
        Attrib(gl, _bufLine, 2, 4, GL_FLOAT, false, 28, 12, 1);
        Data(gl, _bufLineRgb, d.LineRgb);
        Attrib(gl, _bufLineRgb, 3, 4, GL_UNSIGNED_BYTE, true, 4, 0, 1);
        gl.BindVertexArray(0);
        gl.BindBuffer(GL_ARRAY_BUFFER, 0);
    }

    private static void Data<T>(GlInterface gl, int buf, T[] data) where T : unmanaged
    {
        gl.BindBuffer(GL_ARRAY_BUFFER, buf);
        fixed (T* p = data.Length > 0 ? data : new T[1])
            gl.BufferData(GL_ARRAY_BUFFER, (IntPtr)(Math.Max(1, data.Length) * sizeof(T)), (IntPtr)p, GL_STATIC_DRAW);
    }

    private void Attrib(GlInterface gl, int buf, int loc, int size, int type, bool normalized, int stride, int offset, int divisor)
    {
        gl.BindBuffer(GL_ARRAY_BUFFER, buf);
        gl.EnableVertexAttribArray(loc);
        gl.VertexAttribPointer(loc, size, type, normalized ? 1 : 0, stride, (IntPtr)offset);
        _divisor(loc, divisor);
    }

    private static int Buffer(GlInterface gl, float[] data)
    {
        var b = gl.GenBuffer();
        Data(gl, b, data);
        return b;
    }

    private static int Program(GlInterface gl, string vs, string fs)
    {
        var v = gl.CreateShader(GL_VERTEX_SHADER);
        var err = gl.CompileShaderAndGetError(v, vs);
        if (err != null) throw new InvalidOperationException("vertex shader: " + err);
        var f = gl.CreateShader(GL_FRAGMENT_SHADER);
        err = gl.CompileShaderAndGetError(f, fs);
        if (err != null) throw new InvalidOperationException("fragment shader: " + err);
        var p = gl.CreateProgram();
        gl.AttachShader(p, v);
        gl.AttachShader(p, f);
        err = gl.LinkProgramAndGetError(p);
        if (err != null) throw new InvalidOperationException("link: " + err);
        gl.DeleteShader(v);
        gl.DeleteShader(f);
        return p;
    }

    // ---------------------------------------------------------------- shaders
    // The projection of caps_view_fit: r = R(p − c) + pan; screen x = w/2 + r.x·s·k, y = h/2 − r.y·s·k (y down, as the
    // CPU image); view depth r.z (nearer is larger) mapped onto the depth range [uDepth.x, uDepth.y].
    private const string Common = """
        uniform vec3 uCentre; uniform vec4 uRot; uniform vec2 uPan; uniform float uScale; uniform vec2 uSize;
        uniform float uPersp; uniform float uDist; uniform vec2 uDepth; uniform vec2 uZ; uniform vec3 uBg;
        uniform float uCue; uniform float uTransparent; uniform float uOutline; uniform float uPx; uniform vec3 uInk;
        vec3 rotv(vec3 p) {
          vec3 d = p - uCentre;
          float x = d.x * uRot.x + d.z * uRot.y;
          float z = -d.x * uRot.y + d.z * uRot.x;
          float y2 = d.y * uRot.z - z * uRot.w;
          float z2 = d.y * uRot.w + z * uRot.z;
          return vec3(x + uPan.x, y2 + uPan.y, z2);
        }
        float kOf(float z) { return uPersp > 0.5 ? uDist / max(1e-3, uDist - z) : 1.0; }
        vec2 screenOf(vec3 r, float k) { return vec2(uSize.x * 0.5 + r.x * uScale * k, uSize.y * 0.5 - r.y * uScale * k); }
        float ndcDepth(float z) { return clamp(1.0 - 2.0 * (z - uDepth.x) / (uDepth.y - uDepth.x), -1.0, 1.0); }
        vec4 clipOf(vec2 px, float z) { return vec4(px.x / uSize.x * 2.0 - 1.0, 1.0 - px.y / uSize.y * 2.0, ndcDepth(z), 1.0); }

        """;

    // The CPU renderer's material and fog (render.cpp: shade, depth cue).
    private const string Shade = """
        vec3 shade(vec3 base, vec3 n) {
          vec3 L = normalize(vec3(-0.45, 0.6, 1.0));
          vec3 H = normalize(L + vec3(0.0, 0.0, 1.0));
          float diff = max(0.0, dot(n, L));
          float spec = pow(max(0.0, dot(n, H)), 40.0);
          return min(vec3(1.0), base * (0.32 + 0.72 * diff) + 0.5 * spec);
        }
        vec3 cue(vec3 c, float z) {
          if (uCue < 0.5 || uZ.y <= uZ.x) return c;
          float t = clamp(0.35 * (uZ.y - z) / (uZ.y - uZ.x), 0.0, 0.35);
          return uTransparent > 0.5 ? mix(c, vec3(0.0), t * 0.6) : mix(c, uBg, t);
        }
        float fragDepth(float z) { return ndcDepth(z) * 0.5 + 0.5; }

        """;

    private const string SphereVs = """
        layout(location = 0) in vec2 aCorner;
        layout(location = 1) in vec4 aSphere;
        layout(location = 2) in vec4 aColour;
        layout(location = 3) in float aRing;
        out vec2 vOff; out float vR; out float vZ; out float vRw; out vec3 vCol; flat out int vRing;
        void main() {
          vec3 r = rotv(aSphere.xyz);
          float k = kOf(r.z);
          vec2 c = screenOf(r, k);
          float R = aSphere.w * uScale * k;
          int ring = int(aRing + 0.5);
          float ext = R + (ring > 0 ? 9.0 * uPx : 1.5 * uPx);
          vOff = aCorner * ext; vR = R; vZ = r.z; vRw = aSphere.w; vCol = aColour.zyx; vRing = ring;
          gl_Position = clipOf(c + vOff, r.z + aSphere.w);
        }

        """;

    private const string SphereFs = """
        in vec2 vOff; in float vR; in float vZ; in float vRw; in vec3 vCol; flat in int vRing;
        out vec4 frag;
        void main() {
          float d = length(vOff);
          if (d > vR) {   // outside the atom: its selection and focus rings, drawn over everything
            if (vRing == 0) discard;
            bool on = false; vec3 rc = vec3(0.0);
            if ((vRing & 1) != 0 && abs(d - (vR + 3.0 * uPx)) <= 1.2 * uPx) { on = true; rc = vec3(0.4235, 0.7686, 0.8471); }
            if ((vRing & 2) != 0 && abs(d - (vR + 6.5 * uPx)) <= 1.6 * uPx && abs(vOff.y / d) > 0.26) { on = true; rc = vec3(0.9608, 0.6471, 0.1412); }
            if (!on) discard;
            frag = vec4(rc, 1.0);
            gl_FragDepth = 0.0;
            return;
          }
          float nz = sqrt(max(0.0, 1.0 - d * d / (vR * vR)));
          vec3 n = vec3(vOff.x / vR, -vOff.y / vR, nz);
          float z = vZ + vRw * nz;
          vec3 c = cue(shade(vCol, n), z);
          if (uOutline > 0.5 && vR - d < uPx) c = mix(c, uInk, 0.5);
          frag = vec4(c, 1.0);
          gl_FragDepth = fragDepth(z);
        }

        """;

    private const string CapsuleVs = """
        layout(location = 0) in vec2 aCorner;
        layout(location = 1) in vec3 aA;
        layout(location = 2) in vec4 aBR;
        layout(location = 3) in vec4 aColour;
        out vec2 vP; flat out vec2 vPa; flat out vec2 vPb; flat out vec2 vZab; flat out float vR; flat out float vRw; flat out vec3 vCol;
        void main() {
          vec3 ra = rotv(aA), rb = rotv(aBR.xyz);
          float ka = kOf(ra.z), kb = kOf(rb.z);
          vec2 pa = screenOf(ra, ka), pb = screenOf(rb, kb);
          float R = aBR.w * uScale * (ka + kb) * 0.5;
          vec2 e = pb - pa;
          float L = length(e);
          vec2 u = L > 1e-4 ? e / L : vec2(1.0, 0.0);
          vec2 nrm = vec2(-u.y, u.x);
          float ext = R + 1.5 * uPx;
          vec2 p = pa - u * ext + u * (L + 2.0 * ext) * aCorner.x + nrm * ext * aCorner.y;
          vP = p; vPa = pa; vPb = pb; vZab = vec2(ra.z, rb.z); vR = R; vRw = aBR.w; vCol = aColour.zyx;
          gl_Position = clipOf(p, max(ra.z, rb.z) + aBR.w);
        }

        """;

    private const string CapsuleFs = """
        in vec2 vP; flat in vec2 vPa; flat in vec2 vPb; flat in vec2 vZab; flat in float vR; flat in float vRw; flat in vec3 vCol;
        out vec4 frag;
        void main() {
          vec2 e = vPb - vPa;
          float L2 = dot(e, e);
          float t = L2 > 1e-9 ? clamp(dot(vP - vPa, e) / L2, 0.0, 1.0) : 0.0;
          vec2 dv = vP - vPa - t * e;
          float d2 = dot(dv, dv);
          if (d2 > vR * vR) discard;
          float nz = sqrt(max(0.0, 1.0 - d2 / (vR * vR)));
          float z = vZab.x + t * (vZab.y - vZab.x) + vRw * nz;
          vec3 c = cue(shade(vCol, vec3(dv.x / vR, -dv.y / vR, nz)), z);
          if (uOutline > 0.5 && vR - sqrt(d2) < uPx) c = mix(c, uInk, 0.5);
          frag = vec4(c, 1.0);
          gl_FragDepth = fragDepth(z);
        }

        """;

    private const string LineVs = """
        layout(location = 0) in vec2 aCorner;
        layout(location = 1) in vec3 aA;
        layout(location = 2) in vec4 aBW;
        layout(location = 3) in vec4 aColour;
        flat out vec3 vCol;
        void main() {
          vec3 ra = rotv(aA), rb = rotv(aBW.xyz);
          vec2 pa = screenOf(ra, kOf(ra.z)), pb = screenOf(rb, kOf(rb.z));
          vec2 e = pb - pa;
          float L = length(e);
          vec2 u = L > 1e-4 ? e / L : vec2(1.0, 0.0);
          vec2 nrm = vec2(-u.y, u.x);
          float hw = max(0.5, 0.5 * aBW.w * uPx);
          vec2 p = pa + e * aCorner.x + nrm * hw * aCorner.y;
          float z = mix(ra.z, rb.z, aCorner.x) + 0.05;   // a small bias: edges stay visible where they touch atoms
          vCol = aColour.zyx;
          gl_Position = clipOf(p, z);
        }

        """;

    private const string LineFs = """
        flat in vec3 vCol;
        out vec4 frag;
        void main() { frag = vec4(vCol, 1.0); }

        """;
}
