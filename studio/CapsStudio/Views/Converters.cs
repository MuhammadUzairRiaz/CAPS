using System.Globalization;
using Avalonia.Data;
using Avalonia.Data.Converters;

namespace CapsStudio.Views;

/// <summary>For option cards bound to one integer choice: checked when the value equals the parameter; checking sets it.</summary>
public sealed class IntEqualsConverter : IValueConverter
{
    public static readonly IntEqualsConverter Instance = new();
    public object? Convert(object? value, Type targetType, object? parameter, CultureInfo culture) =>
        value is int v && int.TryParse(parameter?.ToString(), NumberStyles.Integer, CultureInfo.InvariantCulture, out var p) && v == p;
    public object? ConvertBack(object? value, Type targetType, object? parameter, CultureInfo culture) =>
        value is true && int.TryParse(parameter?.ToString(), NumberStyles.Integer, CultureInfo.InvariantCulture, out var p) ? p : BindingOperations.DoNothing;
}

/// <summary>The same for a string choice.</summary>
public sealed class StringEqualsConverter : IValueConverter
{
    public static readonly StringEqualsConverter Instance = new();
    public object? Convert(object? value, Type targetType, object? parameter, CultureInfo culture) => value is string v && v == parameter?.ToString();
    public object? ConvertBack(object? value, Type targetType, object? parameter, CultureInfo culture) =>
        value is true && parameter != null ? parameter.ToString() : BindingOperations.DoNothing;
}

/// <summary>True when an integer is below the parameter (one field for several choices).</summary>
public sealed class IntLessConverter : IValueConverter
{
    public static readonly IntLessConverter Instance = new();
    public object? Convert(object? value, Type targetType, object? parameter, CultureInfo culture) =>
        value is int v && int.TryParse(parameter?.ToString(), NumberStyles.Integer, CultureInfo.InvariantCulture, out var p) && v < p;
    public object? ConvertBack(object? value, Type targetType, object? parameter, CultureInfo culture) => BindingOperations.DoNothing;
}
