using System;
using System.Globalization;

namespace CapsStudio.ViewModels;

/// <summary>The random seed of one kind of run: a fixed seed reproduces the run exactly; "new each time" draws a fresh
/// seed for every run and shows it, so any run can be repeated (Keep makes the shown seed the fixed one).</summary>
public sealed class SeedChoice : ObservableObject
{
    private static readonly Random Dice = new();
    private bool _fresh;
    private decimal? _value;
    private long _used;
    private readonly Action? _changed;

    public SeedChoice(int value = 1, bool fresh = false, Action? changed = null)
    {
        _value = Math.Max(1, value);
        _fresh = fresh;
        _changed = changed;
    }

    /// <summary>A new seed for every run.</summary>
    public bool Fresh { get => _fresh; set { if (Set(ref _fresh, value)) { Raise(nameof(Fixed)); Raise(nameof(UsedText)); Raise(nameof(HasUsed)); _changed?.Invoke(); } } }
    public bool Fixed => !_fresh;
    /// <summary>The fixed seed (1 … 2³¹ − 1). An emptied box keeps the last one.</summary>
    public decimal? Value
    {
        get => _value;
        set { if (value == null) return; if (Set(ref _value, Math.Clamp(Math.Round(value.Value), 1, int.MaxValue))) _changed?.Invoke(); }
    }
    /// <summary>The seed the last run used (0 before any).</summary>
    public long Used { get => _used; private set { if (Set(ref _used, value)) { Raise(nameof(UsedText)); Raise(nameof(HasUsed)); } } }
    public bool HasUsed => _used > 0 && _fresh;   // a fixed seed is the one in the box
    public string UsedText => _used > 0 ? "this run: " + _used.ToString(CultureInfo.InvariantCulture) : "";

    /// <summary>The seed for a run: the fixed one, or a fresh one; recorded as Used.</summary>
    public int Take()
    {
        int s;
        lock (Dice) s = _fresh ? Dice.Next(1, int.MaxValue) : (int)(_value ?? 1);
        Used = s;
        return s;
    }

    /// <summary>The seed of the last run becomes the fixed seed (to reproduce it).</summary>
    public void Keep()
    {
        if (_used <= 0) return;
        _value = _used;
        Raise(nameof(Value));
        Fresh = false;
    }

    /// <summary>For recipes and saved settings: the fixed seed, or -1 for "new each time".</summary>
    public long Setting => _fresh ? -1 : (long)(_value ?? 1);
}
