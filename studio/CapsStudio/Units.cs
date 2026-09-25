using System;
using System.Collections.Generic;

namespace CapsStudio;

/// <summary>Display units (design/boards/Units): how quantities are shown — never how files are written (every format keeps
/// its own units and states them). Conversions from SI 2019 constants; the calorie is the thermochemical one (exact).</summary>
public static class DisplayUnits
{
    public static string Energy { get; set; } = "kcal/mol";     // kcal/mol | kJ/mol | eV
    public static string Length { get; set; } = "Å";            // Å | nm
    public static string Pressure { get; set; } = "atm";        // atm | bar | MPa | GPa
    public static string Time { get; set; } = "ps";             // fs | ps | ns
    public static string Density { get; set; } = "g/cm³";       // g/cm³ | kg/m³

    public static readonly string[] Systems = ["CAPS default · Å · ps · kcal/mol · K · atm · g/cm³", "SI-derived · nm · ps · kJ/mol · K · bar · kg/m³", "LAMMPS metal · Å · ps · eV · K · bar"];
    public static void UseSystem(int k)
    {
        (Energy, Length, Pressure, Time, Density) = k switch
        {
            1 => ("kJ/mol", "nm", "bar", "ps", "kg/m³"),
            2 => ("eV", "Å", "bar", "ps", "g/cm³"),
            _ => ("kcal/mol", "Å", "atm", "ps", "g/cm³"),
        };
    }

    private const double KcalToKj = 4.184, EvToKj = 96.48533212, AtmToBar = 1.01325;
    /// <summary>The factor and unit a value in `unit` is shown with (1 and the unit itself when nothing applies).</summary>
    public static (double Factor, string Unit) For(string unit)
    {
        switch (unit)
        {
            case "kcal/mol":
                return Energy switch { "kJ/mol" => (KcalToKj, "kJ/mol"), "eV" => (KcalToKj / EvToKj, "eV"), _ => (1, unit) };
            case "Å":
                return Length == "nm" ? (0.1, "nm") : (1, unit);
            case "Å²":
                return Length == "nm" ? (0.01, "nm²") : (1, unit);
            case "Å³":
                return Length == "nm" ? (0.001, "nm³") : (1, unit);
            case "atm":
                return Pressure switch { "bar" => (AtmToBar, "bar"), "MPa" => (0.101325, "MPa"), "GPa" => (1.01325e-4, "GPa"), _ => (1, unit) };
            case "bar":
                return Pressure switch { "atm" => (1 / AtmToBar, "atm"), "MPa" => (0.1, "MPa"), "GPa" => (1e-4, "GPa"), _ => (1, unit) };
            case "ps":
                return Time switch { "ns" => (1e-3, "ns"), "fs" => (1e3, "fs"), _ => (1, unit) };
            case "g/cm³":
                return Density == "kg/m³" ? (1000, "kg/m³") : (1, unit);
            default:
                return (1, unit);
        }
    }

    /// <summary>The conversion table of the Units board, computed from the constants.</summary>
    public static IEnumerable<(string Quantity, string From, string Value, string Basis)> Table()
    {
        const double kB = 1.380649e-23, NA = 6.02214076e23, e = 1.602176634e-19, u = 1.66053906660e-24;
        var kt = kB * 298.15 * NA / 1000;   // kJ/mol
        yield return ("Energy", "1 kcal/mol", "4.184 kJ/mol", "exact (thermochemical calorie)");
        yield return ("Energy", "1 eV per particle", $"{e * NA / 1000:0.000} kJ/mol · {e * NA / 4184:0.0000} kcal/mol", "e · N_A");
        yield return ("Thermal energy", "k_B T at 298.15 K", $"{kt:0.0000} kJ/mol · {kt / 4.184:0.0000} kcal/mol", "k_B exact");
        yield return ("Pressure", "1 atm", "1.01325 bar", "exact");
        yield return ("Pressure", "1 GPa", "10 000 bar", "exact");
        yield return ("Density", "1 amu/Å³", $"{u / 1e-24:0.00000} g/cm³", "1 u = 1.66053907 × 10⁻²⁴ g");
        yield return ("Diffusion", "1 Å²/ps", "1 × 10⁻⁸ m²/s", "exact");
        yield return ("Dipole", "1 e·Å", $"{e * 1e-10 / 3.33564e-30:0.0000} D", "1 D = 3.33564 × 10⁻³⁰ C·m");
    }
}
