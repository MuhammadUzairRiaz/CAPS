#pragma once
#include <string>
#include <string_view>

namespace caps {

struct Element {
  int z;
  const char* symbol;
  double mass;        // g/mol, IUPAC standard atomic weight
  double covalent;    // Å, Cordero et al. 2008
  double vdw;         // Å, Bondi 1964 (1.70 when not tabulated)
  unsigned rgb;       // display colour 0xRRGGBB
};

const Element& element(int z);                       // z = 0 gives the "unknown" entry
int element_from_symbol(std::string_view s);         // case-insensitive, 0 if unknown
int element_from_mass(double mass, double tol = 0.1);  // nearest standard weight within tol
int element_from_name(std::string_view atom_name);   // "CA1" -> C, "HW2" -> H, "Cl" -> Cl
int max_element();

// Reference data for labels (element_extra.inc: IUPAC/CIAAW, NIST, NUBASE2020 and the Blue Obelisk Data Repository):
// the IUPAC name, the Pauling electronegativity (0 when not tabulated), the mass number of the most abundant isotope
// (of the longest-lived one for elements with no natural composition), and the ground-state electron configuration
// ("[He] 2s2 2p2"; "" when not tabulated).
const char* element_name(int z);
double pauling_electronegativity(int z);
int most_common_mass_number(int z);
const char* electron_configuration(int z);

}  // namespace caps
