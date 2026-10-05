// CAPS interactions and checks (see caps/interactions.hpp).
#include "caps/interactions.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <set>

#include "caps/edit.hpp"
#include "caps/elements.hpp"

namespace caps {

namespace {

bool polar(int z) { return z == 7 || z == 8 || z == 9; }

// Pairs within rc, each once (i < j), with the minimum-image separation when the cell is periodic.
template <class F>
void pairs_within(const System& s, double rc, F&& visit) {
  const size_t n = s.atoms.size();
  if (n < 2) return;
  const bool periodic = s.cell.valid();
  Vec3 lo{1e300, 1e300, 1e300}, hi{-1e300, -1e300, -1e300};
  std::vector<Vec3> p(n);
  for (size_t i = 0; i < n; ++i) {
    p[i] = periodic ? s.cell.wrap(s.atoms[i].pos) : s.atoms[i].pos;
    for (int k = 0; k < 3; ++k) lo[k] = std::min(lo[k], p[i][k]), hi[k] = std::max(hi[k], p[i][k]);
  }
  int nb[3];
  for (int k = 0; k < 3; ++k) nb[k] = std::clamp(int((hi[k] - lo[k]) / rc), 1, 200);
  auto bin_of = [&](const Vec3& q, int k) { return std::clamp(int((q[k] - lo[k]) / ((hi[k] - lo[k]) / nb[k] + 1e-12)), 0, nb[k] - 1); };
  std::vector<std::vector<uint32_t>> bins(static_cast<size_t>(nb[0]) * nb[1] * nb[2]);
  for (size_t i = 0; i < n; ++i) bins[(size_t(bin_of(p[i], 0)) * nb[1] + bin_of(p[i], 1)) * nb[2] + bin_of(p[i], 2)].push_back(uint32_t(i));
  const double rc2 = rc * rc;
  for (int x = 0; x < nb[0]; ++x)
    for (int y = 0; y < nb[1]; ++y)
      for (int z = 0; z < nb[2]; ++z) {
        const auto& here = bins[(size_t(x) * nb[1] + y) * nb[2] + z];
        std::set<size_t> seen;
        for (int dx = -1; dx <= 1; ++dx)
          for (int dy = -1; dy <= 1; ++dy)
            for (int dz = -1; dz <= 1; ++dz) {
              int X = x + dx, Y = y + dy, Z = z + dz;
              if (periodic) X = (X + nb[0]) % nb[0], Y = (Y + nb[1]) % nb[1], Z = (Z + nb[2]) % nb[2];
              else if (X < 0 || Y < 0 || Z < 0 || X >= nb[0] || Y >= nb[1] || Z >= nb[2]) continue;
              const size_t key = (size_t(X) * nb[1] + Y) * nb[2] + Z;
              if (!seen.insert(key).second) continue;
              for (uint32_t i : here)
                for (uint32_t j : bins[key]) {
                  if (j <= i) continue;
                  Vec3 d = p[j] - p[i];
                  if (periodic) d = s.cell.minimum_image(d);
                  const double r2 = dot(d, d);
                  if (r2 <= rc2) visit(i, j, d, std::sqrt(r2));
                }
            }
      }
}

}  // namespace

int molecules_across_edge(const System& s) {
  if (!s.cell.valid()) return 0;
  const auto mol = s.molecules();
  std::set<int> cut;
  for (size_t i = 0; i < s.atoms.size(); ++i) {
    const Vec3 f = s.cell.to_fractional(s.atoms[i].pos - s.cell.origin);
    for (int k = 0; k < 3; ++k)
      if (f[k] < -1e-6 || f[k] >= 1 + 1e-6) { cut.insert(mol[i]); break; }
  }
  return int(cut.size());
}

InteractionReport find_interactions(const System& s, const InteractionOptions& o) {
  InteractionReport R;
  const size_t n = s.atoms.size();
  int nmol = 0;
  const auto mol = s.molecules(&nmol);
  R.molecules = nmol;
  // bonded neighbours and 1-3 partners
  std::vector<std::vector<uint32_t>> nb(n);
  for (const auto& b : s.bonds) nb[b.i].push_back(b.j), nb[b.j].push_back(b.i);
  auto close_in_graph = [&](uint32_t i, uint32_t j, int depth) {
    // is j within `depth` bonds of i?
    std::vector<uint32_t> frontier{i}, next;
    std::set<uint32_t> seen{i};
    for (int d = 0; d < depth; ++d) {
      next.clear();
      for (uint32_t a : frontier)
        for (uint32_t q : nb[a]) {
          if (q == j) return true;
          if (seen.insert(q).second) next.push_back(q);
        }
      frontier.swap(next);
    }
    return false;
  };
  double rmax = 0;
  for (const auto& a : s.atoms) rmax = std::max(rmax, element(a.element).vdw);
  const double rc = std::max(o.hb_distance, 2 * rmax);
  const double cos_limit = std::cos(o.hb_angle * 3.14159265358979323846 / 180.0);
  auto hbond = [&](uint32_t d, uint32_t a, const Vec3& da, double r) {
    // da: from donor d to acceptor a
    for (uint32_t h : nb[d]) {
      if (s.atoms[h].element != 1) continue;
      Vec3 dh = s.atoms[h].pos - s.atoms[d].pos;
      if (s.cell.valid()) dh = s.cell.minimum_image(dh);
      const double c = dot(dh, da) / (norm(dh) * r);
      if (c >= cos_limit) {
        ++R.n_hbonds;
        if (R.hbonds.size() < o.max_pairs) R.hbonds.push_back({d, h, a, r, std::acos(std::clamp(c, -1.0, 1.0)) * 180 / 3.14159265358979323846});
      }
    }
  };
  pairs_within(s, rc, [&](uint32_t i, uint32_t j, const Vec3& d, double r) {
    const int zi = s.atoms[i].element, zj = s.atoms[j].element;
    const bool same = mol[i] == mol[j];
    // hydrogen bonds between polar heavy atoms
    if (polar(zi) && polar(zj) && r <= o.hb_distance && (!same || !close_in_graph(i, j, 3))) {
      hbond(i, j, d, r);
      hbond(j, i, d * -1.0, r);
    }
    // contacts and clashes: not bonded, not 1-3
    const double sum = element(zi).vdw + element(zj).vdw;
    if (r >= sum - o.contact_margin) return;
    if (same && close_in_graph(i, j, 2)) return;
    // a donor hydrogen and its acceptor are an H-bond, not a contact
    if ((zi == 1 && polar(zj)) || (zj == 1 && polar(zi))) {
      const uint32_t h = zi == 1 ? i : j;
      bool on_polar = false;
      for (uint32_t q : nb[h]) on_polar |= polar(s.atoms[q].element);
      if (on_polar && r > 1.5) return;
    }
    if (r < o.clash_factor * sum) {
      ++R.n_clashes;
      if (R.clashes.size() < o.max_pairs) R.clashes.push_back({i, j, r, o.clash_factor * sum});
    } else {
      ++R.n_contacts;
      if (R.contacts.size() < o.max_pairs) R.contacts.push_back({i, j, r, sum - o.contact_margin});
    }
  });
  // checks
  char b[200];
  if (R.n_clashes > 0) {
    double dmin = 1e300;
    Issue is;
    for (const auto& c : R.clashes) {
      dmin = std::min(dmin, c.distance);
      if (is.atoms.size() < 2000) is.atoms.push_back(c.i), is.atoms.push_back(c.j);
    }
    std::snprintf(b, sizeof b, "closest %.2f Å", dmin);
    is.level = "error", is.title = std::to_string(R.n_clashes) + " clash" + (R.n_clashes == 1 ? "" : "es") + " (< 0.75 × vdW)", is.detail = b, is.fix = "push_apart";
    R.issues.push_back(is);
  } else {
    R.issues.push_back({"ok", "No clashes", "no pair closer than 0.75 × vdW", "", {}});
  }
  const int cut = molecules_across_edge(s);
  if (cut > 0) R.issues.push_back({"warn", "Cell edge cuts " + std::to_string(cut) + " molecule" + (cut == 1 ? "" : "s"), "atoms outside the cell: wrap them", "wrap", {}});
  // valences: over-bonded atoms; heavy atoms lacking hydrogens
  int over = 0, lacking = 0;
  std::vector<uint32_t> over_atoms, lack_atoms;
  for (size_t i = 0; i < n; ++i) {
    const int z = s.atoms[i].element;
    const int q = int(std::lround(s.atoms[i].charge));
    const int val = default_valence(z, std::fabs(s.atoms[i].charge - q) < 0.05 ? q : 0);
    if (val == 0) continue;
    double sum = 0;
    for (const auto& bd : s.bonds)
      if (bd.i == i || bd.j == i) sum += bd.order == 2 ? 2 : bd.order == 3 ? 3 : bd.order == 4 ? 1.5 : 1;
    if (sum > val + 0.6) { ++over; if (over_atoms.size() < 500) over_atoms.push_back(uint32_t(i)); }
    else if (z != 1 && sum < val - 0.6 && !s.bonds.empty()) { ++lacking; if (lack_atoms.size() < 500) lack_atoms.push_back(uint32_t(i)); }
  }
  if (over > 0) R.issues.push_back({"error", std::to_string(over) + " atom(s) over their valence", "more bonds than the element takes", "", over_atoms});
  if (lacking > 0) R.issues.push_back({"warn", std::to_string(lacking) + " atom(s) lack hydrogens", "open valences (united atoms or missing H)", "add_h", lack_atoms});
  if (over == 0 && lacking == 0) R.issues.push_back({"ok", "Valences", "every atom within its valence", "", {}});
  for (const auto& a : s.atoms) R.net_charge += a.charge;
  std::snprintf(b, sizeof b, "%.3f e", R.net_charge);
  const bool neutral = std::fabs(R.net_charge) < 0.01;
  R.issues.push_back({neutral ? "ok" : "warn", "Net charge", b, "", {}});
  return R;
}

}  // namespace caps
