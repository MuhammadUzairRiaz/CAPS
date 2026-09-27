// A moltemplate system (.lt) from CAPS's LAMMPS data and input files: the same styles, coefficients (class II cross
// terms included) and topology, the types named after the force field's (see caps/io.hpp write_moltemplate).
#include <algorithm>
#include <cctype>
#include <fstream>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <vector>

#include "caps/io.hpp"

namespace caps {
namespace {

std::vector<std::string> words(const std::string& l) {
  std::istringstream in(l);
  std::vector<std::string> w;
  for (std::string t; in >> t;) w.push_back(t);
  return w;
}

std::string clean(std::string s) {
  for (char& c : s)
    if (!(std::isalnum(static_cast<unsigned char>(c)) || c == '_' || c == '-' || c == '.' || c == '+')) c = '_';
  return s.empty() ? std::string("t") : s;
}

// "1 a b c  # comment" → number, values, comment
struct Row { std::string id; std::vector<std::string> v; std::string comment; };
Row row_of(const std::string& l) {
  Row r;
  const auto hash = l.find('#');
  const std::string body = hash == std::string::npos ? l : l.substr(0, hash);
  if (hash != std::string::npos) {
    r.comment = l.substr(hash + 1);
    const auto a = r.comment.find_first_not_of(' '), b = r.comment.find_last_not_of(" \r");
    r.comment = a == std::string::npos ? "" : r.comment.substr(a, b - a + 1);
  }
  auto w = words(body);
  if (!w.empty()) { r.id = w[0]; r.v.assign(w.begin() + 1, w.end()); }
  return r;
}

std::string join(const std::vector<std::string>& v, size_t from = 0) {
  std::string t;
  for (size_t k = from; k < v.size(); ++k) t += (k > from ? " " : "") + v[k];
  return t;
}

}  // namespace

std::string lammps_to_moltemplate(const std::string& data_path, const std::string& input_path, const std::string& title) {
  std::ifstream din(data_path), iin(input_path);
  if (!din) throw std::runtime_error("cannot read " + data_path);
  if (!iin) throw std::runtime_error("cannot read " + input_path);
  std::vector<std::string> data, input;
  for (std::string l; std::getline(din, l);) data.push_back(l);
  for (std::string l; std::getline(iin, l);) input.push_back(l);

  // the data file's sections
  std::map<std::string, std::vector<Row>> sec;
  std::vector<std::string> box;
  std::string cur;
  static const std::set<std::string> known = {"Masses", "Pair Coeffs", "PairIJ Coeffs", "Bond Coeffs", "Angle Coeffs", "Dihedral Coeffs", "Improper Coeffs",
                                              "BondBond Coeffs", "BondAngle Coeffs", "MiddleBondTorsion Coeffs", "EndBondTorsion Coeffs", "AngleTorsion Coeffs",
                                              "AngleAngleTorsion Coeffs", "BondBond13 Coeffs", "AngleAngle Coeffs", "Atoms", "Velocities", "Bonds", "Angles",
                                              "Dihedrals", "Impropers"};
  for (size_t k = 1; k < data.size(); ++k) {
    const std::string& l = data[k];
    const std::string head = l.substr(0, l.find('#'));
    std::string h = head;
    while (!h.empty() && (h.back() == ' ' || h.back() == '\r')) h.pop_back();
    if (known.count(h)) { cur = h; continue; }
    if (l.find("xlo xhi") != std::string::npos || l.find("ylo yhi") != std::string::npos || l.find("zlo zhi") != std::string::npos ||
        l.find("xy xz yz") != std::string::npos) {
      box.push_back(head);
      continue;
    }
    if (cur.empty()) continue;
    const Row r = row_of(l);
    if (!r.id.empty()) sec[cur].push_back(r);
  }
  if (!sec.count("Atoms")) throw std::runtime_error(data_path + ": no Atoms section");

  // type names: the force field's names from the comments, unique; interactions numbered with their names
  std::map<std::string, std::string> atom_name, bond_name, angle_name, dih_name, imp_name;
  std::set<std::string> used;
  for (const auto& r : sec["Masses"]) {
    std::string n = clean(r.comment.empty() ? "type" + r.id : r.comment);
    std::string u = n;
    for (int k = 2; used.count(u); ++k) u = n + "_" + std::to_string(k);
    used.insert(u);
    atom_name[r.id] = "@atom:" + u;
  }
  auto name_types = [&](const char* section, const char* kind, std::map<std::string, std::string>& out) {
    for (const auto& r : sec[section]) out[r.id] = std::string("@") + kind + ":" + kind[0] + r.id + (r.comment.empty() ? "" : "_" + clean(r.comment));
  };
  name_types("Bond Coeffs", "bond", bond_name);
  name_types("Angle Coeffs", "angle", angle_name);
  name_types("Dihedral Coeffs", "dihedral", dih_name);
  name_types("Improper Coeffs", "improper", imp_name);
  // interaction types without a Coeffs section (coefficients in the input): numbered only
  auto ensure = [&](std::map<std::string, std::string>& m, const char* kind, const std::string& id) {
    if (!m.count(id)) m[id] = std::string("@") + kind + ":" + kind[0] + id;
    return m[id];
  };
  auto atom_of = [&](const std::string& id) {
    if (!atom_name.count(id)) atom_name[id] = "@atom:type" + id;
    return atom_name[id];
  };

  std::ostringstream o;
  o << "# moltemplate system written by CAPS: " << title << "\n"
    << "# the same force field, coefficients and topology as the LAMMPS files CAPS writes for it\n"
    << "# run: moltemplate.sh -overlay-all system.lt   (-overlay-all keeps each interaction's atom order: class II cross terms depend on it)\n\n";
  // In Init: units, styles and settings from the input (up to read_data)
  std::vector<std::string> after;   // kspace and pair_coeff lines after read_data
  o << "write_once(\"In Init\") {\n";
  bool past = false;
  for (const auto& l : input) {
    const auto w = words(l);
    if (w.empty() || w[0][0] == '#') continue;
    if (w[0] == "read_data") { past = true; continue; }
    if (w[0] == "pair_coeff" || w[0] == "kspace_style" || w[0] == "kspace_modify") { after.push_back(l); continue; }
    if (past) continue;   // neighbour settings, thermo and run belong to the user's run script
    if (w[0] == "units" || w[0] == "atom_style" || w[0] == "boundary" || w[0].find("_style") != std::string::npos || w[0] == "special_bonds" ||
        w[0] == "pair_modify" || w[0] == "dielectric")
      o << "  " << join(w) << "\n";
  }
  o << "}\n\n";

  o << "write_once(\"Data Masses\") {\n";
  for (const auto& r : sec["Masses"]) o << "  " << atom_of(r.id) << " " << join(r.v) << "\n";
  o << "}\n\n";

  o << "write_once(\"In Settings\") {\n";
  for (const auto& l : after) {
    auto w = words(l.substr(0, l.find('#')));
    if (w.empty()) continue;
    if (w[0] == "pair_coeff" && w.size() >= 3) {
      if (w[1] != "*") w[1] = atom_of(w[1]);
      if (w[2] != "*") w[2] = atom_of(w[2]);
    }
    o << "  " << join(w) << "\n";
  }
  for (const auto& r : sec["Pair Coeffs"]) o << "  pair_coeff " << atom_of(r.id) << " " << atom_of(r.id) << " " << join(r.v) << "\n";
  for (const auto& r : sec["PairIJ Coeffs"])
    if (!r.v.empty()) o << "  pair_coeff " << atom_of(r.id) << " " << atom_of(r.v[0]) << " " << join(r.v, 1) << "\n";
  for (const auto& r : sec["Bond Coeffs"]) o << "  bond_coeff " << bond_name[r.id] << " " << join(r.v) << "\n";
  for (const auto& r : sec["Angle Coeffs"]) o << "  angle_coeff " << angle_name[r.id] << " " << join(r.v) << "\n";
  for (const auto& r : sec["BondBond Coeffs"]) o << "  angle_coeff " << ensure(angle_name, "angle", r.id) << " bb " << join(r.v) << "\n";
  for (const auto& r : sec["BondAngle Coeffs"]) o << "  angle_coeff " << ensure(angle_name, "angle", r.id) << " ba " << join(r.v) << "\n";
  for (const auto& r : sec["Dihedral Coeffs"]) o << "  dihedral_coeff " << dih_name[r.id] << " " << join(r.v) << "\n";
  for (const auto& r : sec["MiddleBondTorsion Coeffs"]) o << "  dihedral_coeff " << ensure(dih_name, "dihedral", r.id) << " mbt " << join(r.v) << "\n";
  for (const auto& r : sec["EndBondTorsion Coeffs"]) o << "  dihedral_coeff " << ensure(dih_name, "dihedral", r.id) << " ebt " << join(r.v) << "\n";
  for (const auto& r : sec["AngleTorsion Coeffs"]) o << "  dihedral_coeff " << ensure(dih_name, "dihedral", r.id) << " at " << join(r.v) << "\n";
  for (const auto& r : sec["AngleAngleTorsion Coeffs"]) o << "  dihedral_coeff " << ensure(dih_name, "dihedral", r.id) << " aat " << join(r.v) << "\n";
  for (const auto& r : sec["BondBond13 Coeffs"]) o << "  dihedral_coeff " << ensure(dih_name, "dihedral", r.id) << " bb13 " << join(r.v) << "\n";
  for (const auto& r : sec["Improper Coeffs"]) o << "  improper_coeff " << imp_name[r.id] << " " << join(r.v) << "\n";
  for (const auto& r : sec["AngleAngle Coeffs"]) o << "  improper_coeff " << ensure(imp_name, "improper", r.id) << " aa " << join(r.v) << "\n";
  o << "}\n\n";

  // the box
  o << "write_once(\"Data Boundary\") {\n";
  for (const auto& b : box) o << "  " << join(words(b)) << "\n";
  o << "}\n\n";
  // box lengths, for unwrapping with the image flags (orthorhombic and triclinic)
  double lo[3] = {0, 0, 0}, hi[3] = {0, 0, 0}, xy = 0, xz = 0, yz = 0;
  for (const auto& b : box) {
    const auto w = words(b);
    if (w.size() >= 4 && w[2] == "xlo") lo[0] = std::stod(w[0]), hi[0] = std::stod(w[1]);
    if (w.size() >= 4 && w[2] == "ylo") lo[1] = std::stod(w[0]), hi[1] = std::stod(w[1]);
    if (w.size() >= 4 && w[2] == "zlo") lo[2] = std::stod(w[0]), hi[2] = std::stod(w[1]);
    if (w.size() >= 6 && w[3] == "xy") xy = std::stod(w[0]), xz = std::stod(w[1]), yz = std::stod(w[2]);
  }
  const double lx = hi[0] - lo[0], ly = hi[1] - lo[1], lz = hi[2] - lo[2];

  // atoms (full: id mol type q x y z [ix iy iz]), unwrapped so no image flags are needed
  o << "write(\"Data Atoms\") {\n";
  char b[256];
  for (const auto& r : sec["Atoms"]) {
    if (r.v.size() < 6) throw std::runtime_error("an Atoms line with fewer than 7 columns (atom_style full expected)");
    double x = std::stod(r.v[3]), y = std::stod(r.v[4]), z = std::stod(r.v[5]);
    if (r.v.size() >= 9) {
      const double ix = std::stod(r.v[6]), iy = std::stod(r.v[7]), iz = std::stod(r.v[8]);
      x += ix * lx + iy * xy + iz * xz;
      y += iy * ly + iz * yz;
      z += iz * lz;
    }
    std::snprintf(b, sizeof b, "  $atom:a%s $mol:m%s %s %s %.8f %.8f %.8f\n", r.id.c_str(), r.v[0].c_str(), atom_of(r.v[1]).c_str(), r.v[2].c_str(), x, y, z);
    o << b;
  }
  o << "}\n\n";
  auto topo = [&](const char* section, const char* block, const char* kind, std::map<std::string, std::string>& names, size_t natoms) {
    if (!sec.count(section) || sec[section].empty()) return;
    o << "write(\"" << block << "\") {\n";
    for (const auto& r : sec[section]) {
      if (r.v.size() < 1 + natoms) continue;
      o << "  $" << kind << ":" << kind[0] << r.id << " " << ensure(names, kind, r.v[0]);
      for (size_t k = 1; k <= natoms; ++k) o << " $atom:a" << r.v[k];
      o << "\n";
    }
    o << "}\n\n";
  };
  topo("Bonds", "Data Bonds", "bond", bond_name, 2);
  topo("Angles", "Data Angles", "angle", angle_name, 3);
  topo("Dihedrals", "Data Dihedrals", "dihedral", dih_name, 4);
  topo("Impropers", "Data Impropers", "improper", imp_name, 4);
  return o.str();
}

}  // namespace caps
