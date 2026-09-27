// CAPS file checks (see caps/checks.hpp).
#include "caps/checks.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <sstream>

#include "caps/elements.hpp"
#include "caps/json.hpp"
#include "cell_list.hpp"

namespace caps {

namespace {

std::string fmt(const char* f, double a) {
  char b[64];
  std::snprintf(b, sizeof b, f, a);
  return b;
}

std::string num(size_t n) {   // 1 300 with a thin grouping space, as the boards write counts
  std::string s = std::to_string(n), o;
  const int len = int(s.size());
  for (int i = 0; i < len; ++i) {
    if (i && (len - i) % 3 == 0) o += " ";
    o += s[size_t(i)];
  }
  return o;
}

}  // namespace

std::vector<FileCheck> file_checks(const Trajectory& t) {
  std::vector<FileCheck> out;
  const System& s = t.topology;
  const size_t n = s.atoms.size();
  if (n == 0) {
    out.push_back({"error", "No atoms", "The file defines no atoms; nothing to show.", ""});
    return out;
  }
  // atom counts and identifiers
  {
    bool frames_ok = true;
    for (const auto& p : t.positions) frames_ok = frames_ok && p.size() == n;
    std::vector<int64_t> ids;
    for (const auto& a : s.atoms) ids.push_back(a.id);
    std::sort(ids.begin(), ids.end());
    const bool contiguous = ids.front() == 1 && ids.back() == int64_t(n) && std::adjacent_find(ids.begin(), ids.end()) == ids.end();
    if (!frames_ok) out.push_back({"error", "Frames disagree on the atom count", "Some frames hold a different number of atoms than the topology; those frames cannot be shown.", ""});
    else
      out.push_back({"pass", "Atom counts agree",
                     num(n) + " atoms" + (t.frames() > 1 ? " in each of " + num(t.frames()) + " frames" : "") +
                         (contiguous ? " with ids 1–" + num(n) + "." : "; ids are not 1…N (kept as read)."),
                     ""});
  }
  // elements
  {
    size_t none = 0;
    for (const auto& a : s.atoms) none += a.element == 0;
    if (none) out.push_back({"warn", num(none) + " atoms without an element", "No element, mass or recognisable name: they are drawn grey and left untyped. Give masses or element names in the file.", ""});
  }
  // the cell and atoms outside it
  const bool cell = s.cell.valid();
  if (!cell) out.push_back({"note", "No periodic cell", "The structure is treated as a molecule in vacuum: no periodic images, no density.", ""});
  else {
    size_t outside = 0;
    for (const auto& a : s.atoms) {
      const Vec3 f = s.cell.to_fractional(a.pos);
      for (int k = 0; k < 3; ++k)
        if (f[k] < -1e-9 || f[k] >= 1 + 1e-9) { ++outside; break; }
    }
    if (outside)
      out.push_back({"note", num(outside) + " atoms lie outside the box",
                     s.unwrapped ? "Coordinates are unwrapped (whole molecules). Kept as read; wrapping changes only the view."
                                 : "Kept as read; wrapping changes only the view.",
                     "wrap"});
    const double rho = s.density();
    if (rho > 0 && (rho < 0.05 || rho > 8)) out.push_back({"warn", "Unusual density " + fmt("%.3g", rho) + " g/cm³", "Check the box and the masses (a box in nm read as Å gives 1000× the density).", ""});
  }
  // bonds
  if (!s.bonds.empty()) {
    double longest = 0;
    size_t stretched = 0, across = 0;
    for (const auto& b : s.bonds) {
      const Vec3 raw = s.atoms[b.j].pos - s.atoms[b.i].pos;
      const Vec3 d = cell ? s.cell.minimum_image(raw) : raw;
      const double r = norm(d);
      longest = std::max(longest, r);
      const double ref = element(s.atoms[b.i].element).covalent + element(s.atoms[b.j].element).covalent;
      if (r > ref + 0.8) ++stretched;
      if (cell && norm(raw) > r + 1e-6) ++across;
    }
    std::string d = num(s.bonds.size()) + " bonds " + (s.bonds_from_file ? "from the file" : "perceived from distances") + "; the longest is " + fmt("%.2f", longest) + " Å.";
    // a bond longer than half the cell's narrowest width cannot be read by the minimum-image convention: the engines
    // (and CAPS) would join the wrong images of its atoms. Only unwrapped positions (whole molecules, from image flags
    // or a builder) show the bond's real length; the minimum image of a wrapped pair is never longer than half the box.
    if (cell && s.unwrapped) {
      const double v = s.cell.volume();
      const double w = std::min({v / norm(cross(s.cell.b, s.cell.c)), v / norm(cross(s.cell.c, s.cell.a)), v / norm(cross(s.cell.a, s.cell.b))});
      size_t ambiguous = 0;
      for (const auto& b : s.bonds)
      {
        // a plausible bond (within 0.8 Å of the covalent sum) longer than half the width: the cell is narrower than
        // twice the bond. A raw vector far longer than any bond is a bond across the boundary of a periodic network
        // (a sheet, a crystal): its image is the bond, and that is fine.
        const double raw = norm(s.atoms[b.j].pos - s.atoms[b.i].pos);
        const double ref = element(s.atoms[b.i].element).covalent + element(s.atoms[b.j].element).covalent;
        if (raw > 0.5 * w && raw <= ref + 0.8) ++ambiguous;
      }
      if (ambiguous)
        out.push_back({"error", num(ambiguous) + " bonds longer than half the box",
                       "The narrowest cell width is " + fmt("%.2f", w) + " Å: a bond longer than " + fmt("%.2f", 0.5 * w) +
                           " Å has no unique nearest image, so LAMMPS, GROMACS and CAPS join the wrong copies of its atoms. Make the cell larger (a supercell) or check the topology.",
                       ""});
    }
    if (across) d += " " + num(across) + " cross the box boundary (molecules are made whole for analysis).";
    if (stretched)
      out.push_back({"warn", num(stretched) + " bonds are stretched", d + " Bonds more than 0.8 Å beyond the covalent sum are kept; relax the structure or check the topology.", "relax"});
    else out.push_back({"pass", "Bonds are plausible", d, ""});
  } else if (n > 1) {
    out.push_back({"note", "No bonds", "Neither the file nor the distances give bonds (ions, atoms, or a coarse model).", ""});
  }
  // bond angles: squeezed ones (below 70° outside three-membered rings) and flattened tetrahedral centres (four
  // neighbours, an angle above 150°: its configuration is undefined); aromatic bonds that close no ring (a cut ring)
  if (!s.bonds.empty() && n < 2000000) {
    std::vector<std::vector<uint32_t>> nb(n);
    std::vector<int> arom(n, 0);
    for (const auto& b : s.bonds) {
      nb[b.i].push_back(b.j), nb[b.j].push_back(b.i);
      if (b.order == 4) ++arom[b.i], ++arom[b.j];
    }
    auto vec = [&](uint32_t a, uint32_t b) { const Vec3 d = s.atoms[b].pos - s.atoms[a].pos; return cell ? s.cell.minimum_image(d) : d; };
    size_t squeezed = 0, flat = 0;
    double worst_small = 180, worst_flat = 0;
    uint32_t at_small = 0, at_flat = 0;
    for (uint32_t j = 0; j < n; ++j) {
      const auto& L = nb[j];
      if (L.size() < 2 || L.size() > 8 || s.atoms[j].element == 0) continue;   // beads (a coarse model) have no chemical angles
      bool flattened = false;
      for (size_t a = 0; a < L.size(); ++a)
        for (size_t b = a + 1; b < L.size(); ++b) {
          const Vec3 u = vec(j, L[a]), w = vec(j, L[b]);
          const double nu = norm(u), nw = norm(w);
          if (nu < 1e-6 || nw < 1e-6) continue;
          const double ang = std::acos(std::clamp(dot(u, w) / (nu * nw), -1.0, 1.0)) * 180 / M_PI;
          const bool ring3 = std::find(nb[L[a]].begin(), nb[L[a]].end(), L[b]) != nb[L[a]].end();
          if (ang < 70 && !ring3) {
            ++squeezed;
            if (ang < worst_small) worst_small = ang, at_small = j;
          }
          if (L.size() == 4 && s.atoms[j].element == 6 && ang > 150) {
            flattened = true;
            if (ang > worst_flat) worst_flat = ang, at_flat = j;
          }
        }
      if (flattened) ++flat;
    }
    if (squeezed)
      out.push_back({"warn", num(squeezed) + " bond angles below 70°",
                     "The smallest, " + fmt("%.0f", worst_small) + "°, is at atom " + num(at_small + 1) +
                         " (outside a three-membered ring): overlapping groups or a wrong bond. Clean or relax the structure, or check the bonds there.",
                     "relax"});
    if (flat)
      out.push_back({"warn", num(flat) + " flattened tetrahedral carbons",
                     "Atom " + num(at_flat + 1) + " has a bond angle of " + fmt("%.0f", worst_flat) +
                         "°: a four-bonded carbon this flat has no defined configuration (R/S, tacticity). Clean up the geometry (⌘⇧C) before typing or runs.",
                     "relax"});
    size_t cut = 0;
    for (uint32_t i = 0; i < n; ++i) cut += arom[i] == 1;
    if (cut)
      out.push_back({"warn", num(cut) + " aromatic atoms with a single aromatic bond",
                     "An aromatic bond that closes no ring: a ring was cut (at the cell boundary or by a deletion) or the bond orders are wrong. Check the rings there; Field types these atoms as they stand.",
                     ""});
  }
  // close contacts between atoms that are not bonded
  {
    Grid g(s, 0.7);
    std::vector<std::vector<uint32_t>> nb(n);
    for (const auto& b : s.bonds) nb[b.i].push_back(b.j), nb[b.j].push_back(b.i);
    size_t close = 0;
    double dmin = 1e9;
    for (uint32_t i = 0; i < n; ++i) {
      const Vec3& f = g.frac[i];
      g.for_neighbour_bins(g.bin(f, 0), g.bin(f, 1), g.bin(f, 2), [&](const std::vector<uint32_t>& bin) {
        for (uint32_t j : bin) {
          if (j <= i || std::find(nb[i].begin(), nb[i].end(), j) != nb[i].end()) continue;
          const double r = norm(g.sep(i, j));
          if (r < 0.7) ++close, dmin = std::min(dmin, r);
        }
      });
    }
    if (close) out.push_back({"warn", num(close) + " atom pairs closer than 0.7 Å", "The closest unbonded pair is " + fmt("%.2f", dmin) + " Å apart. Relax with push-off before dynamics.", "relax"});
  }
  // charges
  bool all_zero = true;
  for (const auto& a : s.atoms) all_zero = all_zero && std::fabs(a.charge) < 1e-12;
  if (s.has_charges && all_zero) {
    out.push_back({"note", "All charges are zero", "The file carries a charge column, all zero. Assign charges in Field (from the force field, Gasteiger or QEq) if electrostatics matter.", "field"});
  } else if (s.has_charges) {
    double q = 0;
    for (const auto& a : s.atoms) q += a.charge;
    if (std::fabs(q) < 1e-3) out.push_back({"pass", "Charges are neutral", "Sum of q = " + fmt("%+.1e", q) + " e over " + num(n) + " atoms.", ""});
    else out.push_back({"warn", "Net charge " + fmt("%+.4f", q) + " e", "Electrostatics in a periodic cell assume a neutralising background. Check the charges or assign them in Field.", "field"});
  } else {
    out.push_back({"note", "No charges in the file", "Assign a force field in Field (charges from the force field, Gasteiger or QEq) before runs with electrostatics.", "field"});
  }
  // frame spacing
  if (t.timesteps.size() > 2) {
    const int64_t d0 = t.timesteps[1] - t.timesteps[0];
    for (size_t k = 2; k < t.timesteps.size(); ++k)
      if (t.timesteps[k] - t.timesteps[k - 1] != d0) {
        out.push_back({"warn", "Frames " + std::to_string(k - 1) + " → " + std::to_string(k) + " change the timestep spacing",
                       "Step " + std::to_string(t.timesteps[k - 1]) + " to " + std::to_string(t.timesteps[k]) + " after a spacing of " + std::to_string(d0) +
                           ". Time-based analysis uses the real timesteps, not frame numbers.",
                       ""});
        break;
      }
  }
  // what the reader inferred or skipped
  for (const auto& note : s.notes) {
    if (note.find("frame(s)") != std::string::npos && note.size() < 20) continue;   // the frame count, reported above
    const bool skipped = note.find("skip") != std::string::npos || note.find("ignored") != std::string::npos;
    std::string title = skipped ? "Lines skipped while reading" : "Read with an assumption";
    if (note.find("elements guessed") != std::string::npos) title = "Elements guessed from masses";
    else if (note.find("image flags") != std::string::npos) title = "Positions unwrapped with the image flags";
    else if (note.find("perceived") != std::string::npos) title = "Bonds perceived from distances";
    else if (note.find("occupancy") != std::string::npos) title = "Partly occupied sites left out";
    else if (note.find("symmetry operations") != std::string::npos) title = "Crystal expanded by its symmetry";
    std::string d = note;
    if (!d.empty()) d[0] = char(std::toupper(static_cast<unsigned char>(d[0])));
    if (d.back() != '.') d += ".";
    out.push_back({skipped ? "warn" : "note", title, d, ""});
  }
  return out;
}

std::string file_checks_json(const std::vector<FileCheck>& checks) {
  Json a = Json::array();
  for (const auto& c : checks) {
    Json o = Json::object();
    o["level"] = c.level;
    o["title"] = c.title;
    o["detail"] = c.detail;
    o["action"] = c.action;
    a.push_back(o);
  }
  return a.dump();
}

std::string file_checks_text(const std::vector<FileCheck>& checks, const std::string& title) {
  std::ostringstream o;
  o << "# File checks · " << title << "\n\n";
  for (const auto& c : checks) o << "- [" << c.level << "] " << c.title << " — " << c.detail << "\n";
  return o.str();
}

}  // namespace caps
