#include "caps/elements.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <string>

namespace caps {
namespace {

// Masses: IUPAC 2021 conventional weights. Covalent radii: Cordero et al., Dalton Trans. 2008.
// vdW radii: Bondi 1964 where tabulated, otherwise 2.00. Colours: CAPS display palette.
const Element kTable[] = {
    {0, "X", 0.0, 0.77, 1.70, 0xB0B6BC},
    {1, "H", 1.008, 0.31, 1.20, 0xE9ECEF},
    {2, "He", 4.0026, 0.28, 1.40, 0xD9FFFF},
    {3, "Li", 6.94, 1.28, 1.82, 0xCC80FF},
    {4, "Be", 9.0122, 0.96, 2.00, 0xC2FF00},
    {5, "B", 10.81, 0.84, 2.00, 0xFFB5B5},
    {6, "C", 12.011, 0.76, 1.70, 0x8B969E},
    {7, "N", 14.007, 0.71, 1.55, 0x2271DB},
    {8, "O", 15.999, 0.66, 1.52, 0xE35049},
    {9, "F", 18.998, 0.57, 1.47, 0x90E050},
    {10, "Ne", 20.180, 0.58, 1.54, 0xB3E3F5},
    {11, "Na", 22.990, 1.66, 2.27, 0x9B7AD5},
    {12, "Mg", 24.305, 1.41, 1.73, 0x8AFF00},
    {13, "Al", 26.982, 1.21, 2.00, 0xBFA6A6},
    {14, "Si", 28.085, 1.11, 2.10, 0xD6A357},
    {15, "P", 30.974, 1.07, 1.80, 0xFF8000},
    {16, "S", 32.06, 1.05, 1.80, 0xE4C84C},
    {17, "Cl", 35.45, 1.02, 1.75, 0x58B573},
    {18, "Ar", 39.948, 1.06, 1.88, 0x80D1E3},
    {19, "K", 39.098, 2.03, 2.75, 0x8F40D4},
    {20, "Ca", 40.078, 1.76, 2.00, 0x3DFF00},
    {21, "Sc", 44.956, 1.70, 2.00, 0xE6E6E6},
    {22, "Ti", 47.867, 1.60, 2.00, 0xBFC2C7},
    {23, "V", 50.942, 1.53, 2.00, 0xA6A6AB},
    {24, "Cr", 51.996, 1.39, 2.00, 0x8A99C7},
    {25, "Mn", 54.938, 1.39, 2.00, 0x9C7AC7},
    {26, "Fe", 55.845, 1.32, 2.00, 0xE06633},
    {27, "Co", 58.933, 1.26, 2.00, 0xF090A0},
    {28, "Ni", 58.693, 1.24, 1.63, 0x50D050},
    {29, "Cu", 63.546, 1.32, 1.40, 0xC88033},
    {30, "Zn", 65.38, 1.22, 1.39, 0x7D80B0},
    {31, "Ga", 69.723, 1.22, 1.87, 0xC28F8F},
    {32, "Ge", 72.630, 1.20, 2.00, 0x668F8F},
    {33, "As", 74.922, 1.19, 1.85, 0xBD80E3},
    {34, "Se", 78.971, 1.20, 1.90, 0xFFA100},
    {35, "Br", 79.904, 1.20, 1.85, 0xA62929},
    {36, "Kr", 83.798, 1.16, 2.02, 0x5CB8D1},
};
constexpr int kCount = sizeof(kTable) / sizeof(kTable[0]);

// Heavier elements (Rb–Lr; Bk–Lr covalent radii from Pyykkö & Atsumi 2009): IUPAC weights, Cordero covalent radii, Bondi van der Waals radii where tabulated,
// Jmol colours.
const Element kExtra[] = {
    {37, "Rb", 85.468, 2.20, 1.70, 0x702EB0},
    {38, "Sr", 87.62, 1.95, 1.70, 0x00FF00},
    {39, "Y", 88.906, 1.90, 1.70, 0x94FFFF},
    {40, "Zr", 91.224, 1.75, 1.70, 0x94E0E0},
    {41, "Nb", 92.906, 1.64, 1.70, 0x73C2C9},
    {42, "Mo", 95.95, 1.54, 1.70, 0x54B5B5},
    {43, "Tc", 98.0, 1.47, 1.70, 0x3B9E9E},
    {44, "Ru", 101.07, 1.46, 1.70, 0x248F8F},
    {45, "Rh", 102.91, 1.42, 1.70, 0x0A7D8C},
    {46, "Pd", 106.42, 1.39, 1.63, 0x006985},
    {47, "Ag", 107.87, 1.45, 1.72, 0xC0C0C0},
    {48, "Cd", 112.41, 1.44, 1.58, 0xFFD98F},
    {49, "In", 114.82, 1.42, 1.93, 0xA67573},
    {50, "Sn", 118.71, 1.39, 2.17, 0x668080},
    {51, "Sb", 121.76, 1.39, 1.70, 0x9E63B5},
    {52, "Te", 127.60, 1.38, 2.06, 0xD47A00},
    {53, "I", 126.90, 1.39, 1.98, 0x940094},
    {54, "Xe", 131.29, 1.40, 2.16, 0x429EB0},
    {55, "Cs", 132.91, 2.44, 1.70, 0x57178F},
    {56, "Ba", 137.33, 2.15, 1.70, 0x00C900},
    {57, "La", 138.91, 2.07, 1.70, 0x70D4FF},
    {58, "Ce", 140.12, 2.04, 1.70, 0xFFFFC7},
    {59, "Pr", 140.91, 2.03, 1.70, 0xD9FFC7},
    {60, "Nd", 144.24, 2.01, 1.70, 0xC7FFC7},
    {61, "Pm", 145.0, 1.99, 1.70, 0xA3FFC7},
    {62, "Sm", 150.36, 1.98, 1.70, 0x8FFFC7},
    {63, "Eu", 151.96, 1.98, 1.70, 0x61FFC7},
    {64, "Gd", 157.25, 1.96, 1.70, 0x45FFC7},
    {65, "Tb", 158.93, 1.94, 1.70, 0x30FFC7},
    {66, "Dy", 162.50, 1.92, 1.70, 0x1FFFC7},
    {67, "Ho", 164.93, 1.92, 1.70, 0x00FF9C},
    {68, "Er", 167.26, 1.89, 1.70, 0x00E675},
    {69, "Tm", 168.93, 1.90, 1.70, 0x00D452},
    {70, "Yb", 173.05, 1.87, 1.70, 0x00BF38},
    {71, "Lu", 174.97, 1.87, 1.70, 0x00AB24},
    {72, "Hf", 178.49, 1.75, 1.70, 0x4DC2FF},
    {73, "Ta", 180.95, 1.70, 1.70, 0x4DA6FF},
    {74, "W", 183.84, 1.62, 1.70, 0x2194D6},
    {75, "Re", 186.21, 1.51, 1.70, 0x267DAB},
    {76, "Os", 190.23, 1.44, 1.70, 0x266696},
    {77, "Ir", 192.22, 1.41, 1.70, 0x175487},
    {78, "Pt", 195.08, 1.36, 1.75, 0xD0D0E0},
    {79, "Au", 196.97, 1.36, 1.66, 0xFFD123},
    {80, "Hg", 200.59, 1.32, 1.55, 0xB8B8D0},
    {81, "Tl", 204.38, 1.45, 1.96, 0xA6544D},
    {82, "Pb", 207.2, 1.46, 2.02, 0x575961},
    {83, "Bi", 208.98, 1.48, 1.70, 0x9E4FB5},
    {84, "Po", 209.0, 1.40, 1.70, 0xAB5C00},
    {85, "At", 210.0, 1.50, 1.70, 0x754F45},
    {86, "Rn", 222.0, 1.50, 1.70, 0x428296},
    {87, "Fr", 223.0, 2.60, 1.70, 0x420066},
    {88, "Ra", 226.0, 2.21, 1.70, 0x007D00},
    {89, "Ac", 227.0, 2.15, 1.70, 0x70ABFA},
    {90, "Th", 232.04, 2.06, 1.70, 0x00BAFF},
    {91, "Pa", 231.04, 2.00, 1.70, 0x00A1FF},
    {92, "U", 238.03, 1.96, 1.86, 0x008FFF},
    {93, "Np", 237.0, 1.90, 1.70, 0x0080FF},
    {94, "Pu", 244.0, 1.87, 1.70, 0x006BFF},
    {95, "Am", 243.0, 1.80, 1.70, 0x545CF2},
    {96, "Cm", 247.0, 1.69, 1.70, 0x785CE3},
    {97, "Bk", 247.0, 1.68, 1.70, 0x8A4FE3},
    {98, "Cf", 251.0, 1.68, 1.70, 0xA136D4},
    {99, "Es", 252.0, 1.65, 1.70, 0xB31FD4},
    {100, "Fm", 257.0, 1.67, 1.70, 0xB31FBA},
    {101, "Md", 258.0, 1.73, 1.70, 0xB30DA6},
    {102, "No", 259.0, 1.76, 1.70, 0xBD0D87},
    {103, "Lr", 266.0, 1.61, 1.70, 0xC70066},
};

std::string upper(std::string_view s) {
  std::string o(s);
  for (auto& ch : o) ch = static_cast<char>(std::toupper(static_cast<unsigned char>(ch)));
  return o;
}

}  // namespace

int max_element() { return 103; }

const Element& element(int z) {
  if (z > 0 && z < kCount) return kTable[z];
  for (const auto& e : kExtra)
    if (e.z == z) return e;
  return kTable[0];
}

int element_from_symbol(std::string_view s) {
  if (s.empty()) return 0;
  const std::string u = upper(s);
  for (int z = 1; z < kCount; ++z)
    if (upper(kTable[z].symbol) == u) return z;
  for (const auto& e : kExtra)
    if (upper(e.symbol) == u) return e.z;
  return 0;
}

int element_from_mass(double mass, double tol) {
  int best = 0;
  double bd = tol;
  for (int z = 1; z < kCount; ++z) {
    const double d = std::fabs(kTable[z].mass - mass);
    if (d <= bd) { bd = d; best = z; }
  }
  for (const auto& e : kExtra) {
    const double d = std::fabs(e.mass - mass);
    if (d <= bd) { bd = d; best = e.z; }
  }
  return best;
}

int element_from_name(std::string_view name) {
  // Strip leading digits (PDB-style "1HB"), then try two letters, then one.
  size_t i = 0;
  while (i < name.size() && std::isdigit(static_cast<unsigned char>(name[i]))) ++i;
  std::string_view n = name.substr(i);
  if (n.empty()) return 0;
  if (n.size() >= 2 && std::isalpha(static_cast<unsigned char>(n[1])) && std::islower(static_cast<unsigned char>(n[1]))) {
    // Mixed case like "Cl", "Na": trust it.
    if (int z = element_from_symbol(n.substr(0, 2))) return z;
  }
  // All-caps names (GROMACS/PDB): only a few two-letter symbols are common in soft matter.
  const std::string u = upper(n.substr(0, std::min<size_t>(2, n.size())));
  // "CA" stays carbon (alpha / aromatic carbon names); calcium must be written "Ca".
  if (n.size() == 2 && (u == "CL" || u == "NA" || u == "BR" || u == "SI" || u == "MG" || u == "ZN" || u == "FE"))
    if (int z = element_from_symbol(u)) return z;
  return element_from_symbol(n.substr(0, 1));
}

#include "element_extra.inc"

const char* element_name(int z) { return z > 0 && z <= 118 ? kElementExtra[z].name : ""; }
double pauling_electronegativity(int z) { return z > 0 && z <= 118 ? kElementExtra[z].pauling : 0.0; }
int most_common_mass_number(int z) { return z > 0 && z <= 118 ? kElementExtra[z].mass_number : 0; }
const char* electron_configuration(int z) { return z > 0 && z <= 118 ? kElementExtra[z].configuration : ""; }

}  // namespace caps
