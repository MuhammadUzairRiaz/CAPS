// A force field as assigned to one structure, written whole and read back exactly: every atom's type, charge and mass,
// every term with its atoms and parameters, the pair table, the 1-4 and exclusion settings, many-body and hydrogen-bond
// data, virtual sites (numbers to 17 significant digits). A recipe run on another machine (a remote job, a saved Tg
// recipe) takes it as its force field instead of typing the structure again — the same parameters, whatever the
// force-field library there, and the assignment by group (a filler and its matrix, solvents, a water model) kept.
#pragma once

#include <string>

#include "caps/field.hpp"

namespace caps {

// {"format": "caps-assigned-forcefield", "version": 1, "atoms": N, …}
std::string forcefield_to_json(const ForceField& ff);
// Throws FieldError when the text is not such a file.
ForceField forcefield_from_json(const std::string& text);

}  // namespace caps
