using System;
using System.IO;
using System.Runtime.InteropServices;
using System.Text.Json.Nodes;

namespace CapsStudio;

/// <summary>A problem report for the CAPS team: the version and commit, the system, and the newest entry of crash.log.
/// Help › Report a problem opens it as an email to <see cref="Email"/> in the user's mail program (and puts the whole
/// report on the clipboard, for when no mail program is set up). Nothing is sent unless the user sends that email.</summary>
public static class ProblemReport
{
    public const string Email = "support@caps-studio.org";
    // mail programs and Windows cut a mailto link near 2,000 characters; the clipboard copy is never shortened
    private const int MailBodyChars = 1600;

    public static string Text()
    {
        var w = new System.Text.StringBuilder();
        w.AppendLine("What I was doing when it happened:");
        w.AppendLine();
        w.AppendLine();
        w.AppendLine("---- CAPS (filled in by CAPS Studio)");
        w.AppendLine($"Studio {typeof(ProblemReport).Assembly.GetName().Version}");
        try
        {
            var j = JsonNode.Parse(Interop.CapsDocument.BuildInfo())!;
            w.AppendLine($"Core {j["version"]} · commit {j["commit"]}");
        }
        catch (Exception x) { w.AppendLine($"Core not loaded: {x.GetType().Name}: {x.Message}"); }
        w.AppendLine($"System {RuntimeInformation.OSDescription} · {RuntimeInformation.OSArchitecture} · .NET {Environment.Version}");
        if (CrashLog.Newest() is { Length: > 0 } last)
        {
            w.AppendLine();
            w.AppendLine("---- newest error (crash.log)");
            w.AppendLine(last);
        }
        return w.ToString();
    }

    public static Uri Mailto(string text)
    {
        var body = text.Length <= MailBodyChars ? text : text[..MailBodyChars] + "\n… (the whole report is on the clipboard: paste it here)";
        var subject = $"CAPS problem report ({typeof(ProblemReport).Assembly.GetName().Version}, {(OperatingSystem.IsMacOS() ? "macOS" : OperatingSystem.IsWindows() ? "Windows" : "Linux")})";
        return new Uri($"mailto:{Email}?subject={Uri.EscapeDataString(subject)}&body={Uri.EscapeDataString(body)}");
    }
}
