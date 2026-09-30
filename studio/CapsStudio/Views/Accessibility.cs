using Avalonia;
using Avalonia.Automation;
using Avalonia.Controls;
using Avalonia.Controls.Primitives;
using Avalonia.Interactivity;
using Avalonia.VisualTree;

namespace CapsStudio.Views;

/// <summary>Names for assistive technology (VoiceOver, Narrator, Orca) and UI automation: a button, toggle or combo box
/// whose content is not plain text (an icon, an icon with a label) is announced by its visible label, else by its
/// tooltip — not by its class name ("Avalonia.Controls.Panel") — and its tooltip becomes its help text.</summary>
public static class Accessibility
{
    private static bool _done;

    public static void Init()
    {
        if (_done) return;
        _done = true;
        Control.LoadedEvent.AddClassHandler<Button>((c, _) => Label(c), RoutingStrategies.Direct);
        Control.LoadedEvent.AddClassHandler<ToggleButton>((c, _) => Label(c), RoutingStrategies.Direct);
        Control.LoadedEvent.AddClassHandler<ComboBox>((c, _) => Label(c), RoutingStrategies.Direct);
    }

    private static void Label(Control c)
    {
        var tip = ToolTip.GetTip(c) as string;
        if (!string.IsNullOrWhiteSpace(tip) && string.IsNullOrEmpty(AutomationProperties.GetHelpText(c))) AutomationProperties.SetHelpText(c, tip);
        if (!string.IsNullOrEmpty(AutomationProperties.GetName(c))) return;
        if (c is ContentControl { Content: string s } && s.Length > 0) return;   // plain text: its own name
        var label = c.GetVisualDescendants().OfType<TextBlock>().Select(t => t.Text).FirstOrDefault(t => !string.IsNullOrWhiteSpace(t));
        // a combo box is named by what it is for (its tooltip), its selection is read out by the control itself
        var name = c is ComboBox ? tip : label ?? tip;
        if (!string.IsNullOrWhiteSpace(name)) AutomationProperties.SetName(c, name!.Length > 120 ? name[..120] : name);
    }
}
