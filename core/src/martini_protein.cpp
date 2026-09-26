// Martini 2.2 proteins and DSSP. The model's rules come from data/martini/martini22-protein.json (vermouth-martinize's
// martini22 force field and CHARMM36 mappings, converted by bench/ff/convert_vermouth_martini22.py); DSSP follows
// Kabsch & Sander (1983) as DSSP 2.0 (CMBI) implements it, the version vermouth's reference outputs were made with.
#include "caps/martini_protein.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <map>
#include <tuple>
#include <memory>
#include <set>
#include <sstream>
#include <stdexcept>

#include "caps/elements.hpp"
#include "caps/json.hpp"
#include "caps/resolution.hpp"

namespace caps {

namespace {

constexpr double kPi = 3.14159265358979323846;
constexpr double kDeg = kPi / 180;
constexpr double kKJ = 4.184;

const Json& model(const std::string& path) {
  static std::map<std::string, std::unique_ptr<Json>> cache;
  auto& slot = cache[path];
  if (!slot) {
    std::ifstream in(path);
    if (!in) throw std::invalid_argument("cannot open " + path);
    std::stringstream ss;
    ss << in.rdbuf();
    slot = std::make_unique<Json>(Json::parse(ss.str()));
  }
  return *slot;
}

Vec3 pos(const System& s, int i) { return s.atoms[size_t(i)].pos; }

}  // namespace

std::vector<ProteinResidue> protein_residues(const System& s) {
  std::vector<ProteinResidue> out;
  const std::vector<int> mol = s.molecules();
  // a residue: its molecule, its chain (a PDB's chain identifier is kept as the atom's mol) and its number; chains
  // joined by disulfides are one molecule but keep their own residues
  std::map<std::tuple<int, int64_t, int64_t>, size_t> index;
  for (uint32_t a = 0; a < s.atoms.size(); ++a) {
    const Atom& at = s.atoms[a];
    auto key = std::make_tuple(mol[a], at.mol, at.resid);
    auto it = index.find(key);
    if (it == index.end()) {
      it = index.emplace(key, out.size()).first;
      ProteinResidue r;
      r.name = at.resname;
      r.resid = at.resid;
      r.chain = mol[a];
      out.push_back(r);
    }
    ProteinResidue& r = out[it->second];
    r.atoms.push_back(a);
    const std::string& n = at.name;
    if (n == "N") r.n = int(a);
    else if (n == "CA") r.ca = int(a);
    else if (n == "C") r.c = int(a);
    else if (n == "O" || (n == "OT1" && r.o < 0)) r.o = int(a);
  }
  // chains: a residue continues the previous one when their C and N are bonded
  std::set<std::pair<uint32_t, uint32_t>> bonded;
  for (const auto& b : s.bonds) bonded.insert({std::min(b.i, b.j), std::max(b.i, b.j)});
  int chain = 0;
  for (size_t k = 0; k < out.size(); ++k) {
    bool cont = false;
    if (k > 0 && out[k - 1].c >= 0 && out[k].n >= 0) {
      const uint32_t a = uint32_t(out[k - 1].c), b = uint32_t(out[k].n);
      cont = bonded.count({std::min(a, b), std::max(a, b)}) > 0;
    }
    if (k > 0 && !cont) ++chain;
    out[k].chain = chain;
  }
  return out;
}

// ---------------------------------------------------------------------------------------------------- DSSP

std::string dssp(const System& s) {
  // residues with no backbone atom at all (water, ions, ligands) are not part of the protein
  std::vector<ProteinResidue> res;
  for (auto& r : protein_residues(s))
    if (r.n >= 0 || r.ca >= 0 || r.c >= 0) res.push_back(std::move(r));
  const size_t n = res.size();
  // N, CA and C are needed; a residue without its O (a terminal one written without it) accepts no H-bond
  for (const auto& r : res)
    if (r.n < 0 || r.ca < 0 || r.c < 0) throw std::invalid_argument("DSSP: residue " + r.name + std::to_string(r.resid) + " lacks backbone atoms");
  auto brk = [&](size_t a, size_t b) {   // a chain break between residues a < b
    for (size_t k = a; k < b; ++k)
      if (res[k].chain != res[k + 1].chain) return true;
    return false;
  };
  // amide H as DSSP places it: N plus the unit vector from the previous residue's O to its C
  std::vector<Vec3> H(n), N(n), CA(n), C(n), O(n);
  for (size_t k = 0; k < n; ++k) {
    N[k] = pos(s, res[k].n), CA[k] = pos(s, res[k].ca), C[k] = pos(s, res[k].c), O[k] = res[k].o >= 0 ? pos(s, res[k].o) : C[k];
    H[k] = N[k];
  }
  for (size_t k = 1; k < n; ++k)
    if (res[k].chain == res[k - 1].chain && res[k - 1].o >= 0) {
      const Vec3 co = C[k - 1] - O[k - 1];
      H[k] = N[k] + co * (1 / norm(co));
    }
  // H-bond energies (kcal/mol), the two best acceptors of each donor
  constexpr double kCoupling = -27.888, kMin = -9.9, kMax = -0.5;
  struct HB { int res = -1; double e = 0; };
  std::vector<std::array<HB, 2>> acc(n);   // acc[donor]: the CO partners of donor's NH
  auto energy = [&](size_t d, size_t a) {
    if (res[d].name == "PRO" || res[a].o < 0) return 0.0;
    const double ho = norm(H[d] - O[a]), hc = norm(H[d] - C[a]), nc = norm(N[d] - C[a]), no = norm(N[d] - O[a]);
    double e;
    if (ho < 0.5 || hc < 0.5 || nc < 0.5 || no < 0.5) e = kMin;
    else e = kCoupling / ho - kCoupling / hc + kCoupling / nc - kCoupling / no;
    e = std::round(e * 1000) / 1000;   // DSSP rounds the energy to 0.001
    return std::max(e, kMin);
  };
  auto record = [&](size_t d, size_t a) {
    const double e = energy(d, a);
    auto& h = acc[d];
    if (e < h[0].e) h[1] = h[0], h[0] = {int(a), e};
    else if (e < h[1].e) h[1] = {int(a), e};
  };
  for (size_t i = 0; i + 1 < n; ++i)
    for (size_t j = i + 1; j < n; ++j) {
      if (norm(CA[i] - CA[j]) >= 9.0) continue;
      record(i, j);
      if (j != i + 1) record(j, i);
    }
  auto bond = [&](long d, long a) {   // NH of d to CO of a
    if (d < 0 || a < 0 || d >= long(n) || a >= long(n)) return false;
    for (const auto& h : acc[size_t(d)])
      if (h.res == a && h.e < kMax) return true;
    return false;
  };
  enum SS { Loop, Alpha, Bridge, Strand, H3, H5, Turn, Bend };
  std::vector<SS> ss(n, Loop);
  // bends: CA(i−2) → CA(i) → CA(i+2) turning by more than 70°
  std::vector<char> bend(n, 0);
  for (size_t i = 2; i + 2 < n; ++i) {
    if (brk(i - 2, i + 2)) continue;
    const Vec3 u = CA[i] - CA[i - 2], v = CA[i + 2] - CA[i];
    const double c = dot(u, v) / (norm(u) * norm(v));
    bend[i] = std::acos(std::clamp(c, -1.0, 1.0)) / kDeg > 70;
  }
  // β bridges and ladders
  struct Ladder { bool parallel; std::vector<long> i, j; };
  std::vector<Ladder> ladders;
  auto test_bridge = [&](long i, long j) -> int {   // 1 parallel, 2 antiparallel, 0 none
    if (i < 1 || j < 1 || i + 1 >= long(n) || j + 1 >= long(n)) return 0;
    if (brk(size_t(i - 1), size_t(i + 1)) || brk(size_t(j - 1), size_t(j + 1))) return 0;
    const long a = i - 1, b = i, c = i + 1, d = j - 1, e = j, f = j + 1;
    if ((bond(c, e) && bond(e, a)) || (bond(f, b) && bond(b, d))) return 1;
    if ((bond(c, d) && bond(f, a)) || (bond(e, b) && bond(b, e))) return 2;
    return 0;
  };
  for (long i = 1; i + 4 < long(n); ++i)
    for (long j = i + 3; j + 1 < long(n); ++j) {
      const int t = test_bridge(i, j);
      if (!t) continue;
      bool found = false;
      for (auto& l : ladders) {
        if (l.parallel != (t == 1) || i != l.i.back() + 1) continue;
        if (t == 1 && l.j.back() + 1 == j) { l.i.push_back(i); l.j.push_back(j); found = true; break; }
        if (t == 2 && l.j.front() - 1 == j) { l.i.push_back(i); l.j.insert(l.j.begin(), j); found = true; break; }
      }
      if (!found) ladders.push_back({t == 1, {i}, {j}});
    }
  // ladders joined across bulges
  std::sort(ladders.begin(), ladders.end(), [](const Ladder& a, const Ladder& b) { return a.i.front() < b.i.front(); });
  for (bool merged = true; merged;) {
    merged = false;
    for (size_t x = 0; x < ladders.size() && !merged; ++x)
      for (size_t y = x + 1; y < ladders.size() && !merged; ++y) {
        Ladder& p = ladders[x];
        Ladder& q = ladders[y];
        const long ibi = p.i.front(), iei = p.i.back(), jbi = p.j.front(), jei = p.j.back();
        const long ibj = q.i.front(), iej = q.i.back(), jbj = q.j.front(), jej = q.j.back();
        if (p.parallel != q.parallel || brk(size_t(std::min(ibi, ibj)), size_t(std::max(iei, iej))) ||
            brk(size_t(std::min(jbi, jbj)), size_t(std::max(jei, jej))) || ibj - iei >= 6 || (iei >= ibj && ibi <= iej))
          continue;
        const bool bulge = p.parallel ? ((jbj - jei < 6 && ibj - iei < 3) || jbj - jei < 3) : ((jbi - jej < 6 && ibj - iei < 3) || jbi - jej < 3);
        if (!bulge) continue;
        p.i.insert(p.i.end(), q.i.begin(), q.i.end());
        if (p.parallel) p.j.insert(p.j.end(), q.j.begin(), q.j.end());
        else p.j.insert(p.j.begin(), q.j.begin(), q.j.end());
        ladders.erase(ladders.begin() + long(y));
        merged = true;
      }
  }
  for (const auto& l : ladders) {
    const SS t = l.i.size() > 1 ? Strand : Bridge;
    for (long k = l.i.front(); k <= l.i.back(); ++k)
      if (ss[size_t(k)] != Strand) ss[size_t(k)] = t;
    for (long k = *std::min_element(l.j.begin(), l.j.end()); k <= *std::max_element(l.j.begin(), l.j.end()); ++k)
      if (ss[size_t(k)] != Strand) ss[size_t(k)] = t;
  }
  // helices and turns: flags per stride (none, start, end, middle, start-and-end)
  enum Flag { None, Start, End, Middle, StartEnd };
  std::vector<std::array<Flag, 6>> flag(n);
  for (auto& f : flag) f.fill(None);
  for (int st = 3; st <= 5; ++st)
    for (size_t i = 0; i + st < n; ++i)
      if (!brk(i, i + size_t(st)) && bond(long(i) + st, long(i))) {
        flag[i + size_t(st)][size_t(st)] = End;
        for (size_t j = i + 1; j < i + size_t(st); ++j)
          if (flag[j][size_t(st)] == None) flag[j][size_t(st)] = Middle;
        flag[i][size_t(st)] = flag[i][size_t(st)] == End ? StartEnd : Start;
      }
  auto start = [&](long i, int st) { return i >= 0 && i < long(n) && (flag[size_t(i)][size_t(st)] == Start || flag[size_t(i)][size_t(st)] == StartEnd); };
  for (long i = 1; i + 4 < long(n); ++i)
    if (start(i, 4) && start(i - 1, 4))
      for (long j = i; j <= i + 3; ++j) ss[size_t(j)] = Alpha;
  for (long i = 1; i + 3 < long(n); ++i)
    if (start(i, 3) && start(i - 1, 3)) {
      bool empty = true;
      for (long j = i; j <= i + 2; ++j) empty = empty && (ss[size_t(j)] == Loop || ss[size_t(j)] == H3);
      if (empty)
        for (long j = i; j <= i + 2; ++j) ss[size_t(j)] = H3;
    }
  for (long i = 1; i + 5 < long(n); ++i)
    if (start(i, 5) && start(i - 1, 5)) {
      bool empty = true;
      for (long j = i; j <= i + 4; ++j) empty = empty && (ss[size_t(j)] == Loop || ss[size_t(j)] == H5);
      if (empty)
        for (long j = i; j <= i + 4; ++j) ss[size_t(j)] = H5;
    }
  for (long i = 1; i + 1 < long(n); ++i) {
    if (ss[size_t(i)] != Loop) continue;
    bool turn = false;
    for (int st = 3; st <= 5 && !turn; ++st)
      for (int k = 1; k < st && !turn; ++k) turn = start(i - k, st);
    if (turn) ss[size_t(i)] = Turn;
    else if (bend[size_t(i)]) ss[size_t(i)] = Bend;
  }
  std::string out(n, 'C');
  for (size_t k = 0; k < n; ++k) out[k] = "CHBEGITS"[ss[k]];
  return out;
}

std::string dssp_to_martini(const std::string& letters, const std::string& data_path) {
  const Json& m = model(data_path)["ss"];
  std::string cg;
  for (char c : letters) {
    const std::string key(1, c);
    if (!m["dssp_to_cg"].has(key)) throw std::invalid_argument(std::string("secondary structure letter '") + c + "' is not DSSP's");
    cg += m["dssp_to_cg"][key].str();
  }
  std::string wild = ".";
  for (char c : cg) wild += c == 'H' ? 'H' : '.';
  wild += '.';
  for (const auto& p : m["helix_patterns"].items()) {
    const std::string a = p[0].str(), b = p[1].str();
    for (size_t k; (k = wild.find(a)) != std::string::npos;) wild.replace(k, a.size(), b);
  }
  wild = wild.substr(1, wild.size() - 2);
  for (size_t k = 0; k < cg.size(); ++k)
    if (wild[k] != '.') cg[k] = wild[k];
  return cg;
}

// ---------------------------------------------------------------------------------------------------- Martini 2.2

System martini22_protein(const System& aa, const std::string& ss_in, const std::string& data_path, MartiniProteinReport* rep_out) {
  MartiniProteinReport rep;
  const Json& m = model(data_path);
  const Json& R = m["residues"];
  const Json& L = m["links"];
  auto res = protein_residues(aa);
  const size_t nres = res.size();
  if (!nres) throw std::invalid_argument("no residues");
  const std::string dssp_letters = ss_in.empty() ? dssp(aa) : ss_in;
  if (dssp_letters.size() != nres)
    throw std::invalid_argument("the secondary structure has " + std::to_string(dssp_letters.size()) + " letters for " + std::to_string(nres) + " residues");
  const std::string cg = dssp_to_martini(dssp_letters, data_path);
  rep.dssp = dssp_letters, rep.cg_ss = cg;
  const auto nb = aa.neighbours();
  auto has_h = [&](uint32_t a, const std::string& hname = "") {
    for (uint32_t v : nb[a])
      if (aa.atoms[v].element == 1 && (hname.empty() || aa.atoms[v].name == hname)) return true;
    return false;
  };
  auto atom_named = [&](const ProteinResidue& r, const std::string& name) -> int {
    for (uint32_t a : r.atoms)
      if (aa.atoms[a].name == name) return int(a);
    return -1;
  };
  // the Martini residue: histidine by its protonation, neutral acids / bases by their extra hydrogen
  std::vector<std::string> mres(nres);
  for (size_t k = 0; k < nres; ++k) {
    std::string n = res[k].name;
    auto prot = [&](const char* heavy) { const int a = atom_named(res[k], heavy); return a >= 0 && has_h(uint32_t(a)); };
    if (n == "HIS" || n == "HID" || n == "HIE" || n == "HIP" || n == "HSD" || n == "HSE" || n == "HSP") {
      const bool d = prot("ND1"), e = prot("NE2"), any = [&] { for (uint32_t a : res[k].atoms) if (aa.atoms[a].element == 1) return true; return false; }();
      n = !any || (d && e) ? "HSP" : d ? "HSD" : "HSE";   // no hydrogens at all: the mapping's HIS, which is HSP
    } else if (n == "ASP" && prot("OD1") + prot("OD2")) n = "ASP0";
    else if (n == "GLU" && prot("OE1") + prot("OE2")) n = "GLU0";
    else if (n == "LYS") {
      const int nz = atom_named(res[k], "NZ");
      int h = 0;
      if (nz >= 0) for (uint32_t v : nb[size_t(nz)]) h += aa.atoms[v].element == 1;
      if (nz >= 0 && h == 2) n = "LSN";
    } else if (n == "ARG") {
      int h = 0;
      for (const char* x : {"NE", "NH1", "NH2"}) {
        const int a = atom_named(res[k], x);
        if (a >= 0) for (uint32_t v : nb[size_t(a)]) h += aa.atoms[v].element == 1;
      }
      if (h == 4) n = "ARG0";
    }
    if (!R.has(n)) throw std::invalid_argument("residue " + res[k].name + std::to_string(res[k].resid) + ": Martini 2.2 has no " + n);
    mres[k] = n;
  }
  // beads: the listed heavy atoms of each bead (and the hydrogens on those the mapping lists hydrogens for), mass-weighted
  System out;
  out.title = aa.title;
  out.cell = aa.cell;
  out.source_format = "caps-martini22-protein";
  std::vector<std::map<std::string, uint32_t>> bead(nres);
  std::vector<char> nter(nres, 0), cter(nres, 0);
  for (size_t k = 0; k < nres; ++k) {
    nter[k] = k == 0 || res[k - 1].chain != res[k].chain;
    cter[k] = k + 1 == nres || res[k + 1].chain != res[k].chain;
  }
  auto alias = [&](const std::string& resn, const std::string& atom) {
    if ((resn == "ILE") && atom == "CD") return std::string("CD1");
    return atom;
  };
  for (size_t k = 0; k < nres; ++k) {
    const Json& t = R[mres[k]];
    for (const auto& b : t["beads"].items()) {
      const std::string bn = b["name"].str();
      std::vector<uint32_t> atoms;
      for (const auto& e : t["atoms"][bn].items()) {
        int a = atom_named(res[k], alias(mres[k], e[0].str()));
        if (a < 0) a = atom_named(res[k], e[0].str());
        if (a < 0) throw std::invalid_argument("residue " + res[k].name + std::to_string(res[k].resid) + " lacks atom " + e[0].str());
        atoms.push_back(uint32_t(a));
        if (e[1].number() > 0)
          for (uint32_t v : nb[size_t(a)])
            if (aa.atoms[v].element == 1) atoms.push_back(v);
      }
      if (bn == "BB") {   // termini: the C-terminal OXT (and its H) join the backbone bead; the N-terminal N's hydrogens are all on N
        const int oxt = atom_named(res[k], "OXT") >= 0 ? atom_named(res[k], "OXT") : atom_named(res[k], "OT2");
        if (oxt >= 0) {
          atoms.push_back(uint32_t(oxt));
          for (uint32_t v : nb[size_t(oxt)]) if (aa.atoms[v].element == 1) atoms.push_back(v);
        }
      }
      const Vec3 ref = aa.atoms[atoms[0]].pos;
      Vec3 c{0, 0, 0};
      double w = 0;
      for (uint32_t a : atoms) {
        const double ma = element(aa.atoms[a].element).mass;
        const Vec3 d = aa.cell.valid() ? aa.cell.minimum_image(aa.atoms[a].pos - ref) : aa.atoms[a].pos - ref;
        c = c + (ref + d) * ma;
        w += ma;
      }
      Atom at;
      at.id = int64_t(out.atoms.size() + 1);
      at.mol = res[k].chain + 1;
      at.resname = mres[k];
      at.resid = res[k].resid;
      at.element = 0;
      at.name = b["type"].str();
      at.charge = b["charge"].number();
      at.pos = c * (1 / w);
      bead[k][bn] = uint32_t(out.atoms.size());
      out.atoms.push_back(at);
    }
  }
  auto in = [](const std::string& list, const std::string& x) {   // "ALA PRO HYP" contains x
    std::istringstream is(list);
    for (std::string w; is >> w;) if (w == x) return true;
    return false;
  };
  auto has_code = [](const std::string& codes, char c) { return codes.find(c) != std::string::npos; };
  // backbone bead types by secondary structure (not on termini, whose modifications keep their own types)
  for (size_t k = 0; k < nres; ++k) {
    Atom& bb = out.atoms[bead[k].at("BB")];
    if (nter[k] || cter[k]) continue;
    for (const auto& r : L["bb_type"].items()) {
      const std::string rs = r["res"].str();
      if ((rs == "*" || in(rs, mres[k])) && has_code(r["ss"].str(), cg[k])) bb.name = r["type"].str();
    }
  }
  const Json& term = L["termini"];
  for (size_t k = 0; k < nres; ++k) {
    Atom& bb = out.atoms[bead[k].at("BB")];
    if (nter[k]) {   // N-ter (NH3+) or NH2-ter by the N's hydrogens
      const int nn = res[k].n;
      int h = 0;
      if (nn >= 0) for (uint32_t v : nb[size_t(nn)]) h += aa.atoms[v].element == 1;
      if (h == 2) bb.name = "P5", bb.charge = 0;   // NH2-ter
      else bb.name = term["N"]["type"].str(), bb.charge += term["N"]["charge"].number();
    }
    if (cter[k]) {   // C-ter (COO−) or COOH-ter by a hydrogen on the carboxyl oxygens
      const int oxt = atom_named(res[k], "OXT");
      const bool cooh = (oxt >= 0 && has_h(uint32_t(oxt))) || (res[k].o >= 0 && has_h(uint32_t(res[k].o)));
      if (cooh) bb.name = "P5", bb.charge = 0;
      else bb.name = term["C"]["type"].str(), bb.charge += term["C"]["charge"].number();
    }
  }
  // the topology
  auto topo = std::make_shared<ExplicitTopology>();
  topo->source = "Martini 2.2 protein";
  const double kc = L["constraint_k"].number();
  std::set<std::pair<uint32_t, uint32_t>> have;
  auto add_bond = [&](uint32_t i, uint32_t j, double r0_nm, double k, const std::string& g) {
    if (!have.insert({std::min(i, j), std::max(i, j)}).second) return;
    topo->bonds.push_back({i, j, k / (2 * kKJ * 100), r0_nm * 10, g});
  };
  auto add_angle = [&](uint32_t i, uint32_t j, uint32_t k, double th, double kk, const std::string& g) {
    topo->angles.push_back({i, j, k, 1, kk / (2 * kKJ), th * kDeg, g});
  };
  auto idx = [&](size_t k, const Json& v) { return bead[k].at(R[mres[k]]["beads"][size_t(v.number())]["name"].str()); };
  for (size_t k = 0; k < nres; ++k) {   // side chains
    const Json& t = R[mres[k]];
    for (const auto& b : t["bonds"].items()) add_bond(idx(k, b[0]), idx(k, b[1]), b[2].number(), b[3].number(), "side chain");
    for (const auto& b : t["constraints"].items()) add_bond(idx(k, b[0]), idx(k, b[1]), b[2].number(), kc, "side chain constraint");
    for (const auto& a : t["angles"].items()) {
      if (a[3].number() != 2) throw std::invalid_argument("angle function " + std::to_string(int(a[3].number())) + " in " + mres[k]);
      add_angle(idx(k, a[0]), idx(k, a[1]), idx(k, a[2]), a[4].number(), a[5].number(), "side chain");
    }
    for (const auto& d : t["dihedrals"].items()) {
      const int f = int(d[4].number());
      const uint32_t a = idx(k, d[0]), b = idx(k, d[1]), c = idx(k, d[2]), e = idx(k, d[3]);
      if (f == 2) topo->dihedrals.push_back({a, b, c, e, 2, d[6].number() / (2 * kKJ), d[5].number() * kDeg, 0, "side chain"});
      else if (f == 1) topo->dihedrals.push_back({a, b, c, e, 1, d[6].number() / kKJ, d[5].number() * kDeg, d.size() > 7 ? int(d[7].number()) : 1, "side chain"});
      else throw std::invalid_argument("dihedral function " + std::to_string(f) + " in " + mres[k]);
    }
  }
  const Json& bbb = L["bb_bond"];
  const std::string helix = bbb["helix"].str(), other = bbb["other"].str();
  auto same_chain = [&](size_t a, size_t b) { return res[a].chain == res[b].chain; };
  for (size_t k = 0; k + 1 < nres; ++k) {   // backbone bonds
    if (!same_chain(k, k + 1)) continue;
    const uint32_t a = bead[k].at("BB"), b = bead[k + 1].at("BB");
    const bool ha = has_code(helix, cg[k]), hb = has_code(helix, cg[k + 1]);
    if (ha && hb) add_bond(a, b, bbb["helix_helix_constraint"].number(), kc, "backbone constraint");
    else if ((ha && has_code(other, cg[k + 1])) || (hb && has_code(other, cg[k]))) add_bond(a, b, bbb["helix_other_constraint"].number(), kc, "backbone constraint");
    else add_bond(a, b, bbb["coil"][0].number(), bbb["coil"][1].number(), "backbone");
  }
  for (size_t k = 1; k + 1 < nres; ++k) {   // BBB angles: the last rule that applies
    if (!same_chain(k - 1, k + 1)) continue;
    double th = 0, kk = 0;
    bool any = false;
    for (const auto& r : L["bbb_angle"].items()) {
      const std::string rs = r["res"].str(), codes = r["ss"].str();
      auto ok = [&](size_t x) {
        const bool rn = rs == "*" || (rs[0] == '!' ? !in(rs.substr(1), mres[x]) : in(rs, mres[x]));
        return rn && has_code(codes, cg[x]);
      };
      const bool hit = r.has("middle_only") ? ok(k) : (ok(k - 1) || ok(k) || ok(k + 1));
      if (hit) th = r["theta"].number(), kk = r["k"].number(), any = true;
    }
    if (any) add_angle(bead[k - 1].at("BB"), bead[k].at("BB"), bead[k + 1].at("BB"), th, kk, "BBB");
  }
  const Json& bbs = L["bbs_angle"];
  const Json& sbb = L["first_sbb_angle"];
  for (size_t k = 0; k < nres; ++k) {   // BBS and the first residue's SBB
    if (!bead[k].count("SC1")) continue;
    if (k > 0 && same_chain(k - 1, k))
      add_angle(bead[k - 1].at("BB"), bead[k].at("BB"), bead[k].at("SC1"), bbs[0].number(), bbs[1].number(), "BBS");
    else if (k + 1 < nres && same_chain(k, k + 1))
      add_angle(bead[k].at("SC1"), bead[k].at("BB"), bead[k + 1].at("BB"), sbb[0].number(), sbb[1].number(), "first SBB");
  }
  const Json& hd = L["helix_dihedral"];
  for (size_t k = 0; k + 3 < nres; ++k) {   // helix dihedrals: four helical residues in a row
    if (!same_chain(k, k + 3)) continue;
    bool all = true;
    for (size_t x = k; x <= k + 3; ++x) all = all && has_code(hd["ss"].str(), cg[x]);
    if (all)
      topo->dihedrals.push_back({bead[k].at("BB"), bead[k + 1].at("BB"), bead[k + 2].at("BB"), bead[k + 3].at("BB"), 1,
                                 hd["k"].number() / kKJ, hd["phi"].number() * kDeg, int(hd["n"].number()), "backbone dihedral"});
  }
  const Json& el = L["extended_elastic"];
  for (size_t k = 0; k + 3 < nres; ++k) {   // elastic bonds of extended regions: four extended residues in a row
    if (!same_chain(k, k + 3)) continue;
    bool all = true;
    for (size_t x = k; x <= k + 3; ++x) all = all && has_code(el["ss"].str(), cg[x]);
    if (!all) continue;
    add_bond(bead[k].at("BB"), bead[k + 2].at("BB"), el["short"][0].number(), el["short"][1].number(), "elastic (extended)");
    add_bond(bead[k + 1].at("BB"), bead[k + 3].at("BB"), el["short"][0].number(), el["short"][1].number(), "elastic (extended)");
    add_bond(bead[k].at("BB"), bead[k + 3].at("BB"), el["long"][0].number(), el["long"][1].number(), "elastic (extended)");
  }
  // disulfides: bonded SG atoms
  for (size_t a = 0; a < nres; ++a)
    for (size_t b = a + 1; b < nres; ++b) {
      if (mres[a] != "CYS" || mres[b] != "CYS") continue;
      const int sa = atom_named(res[a], "SG"), sb = atom_named(res[b], "SG");
      if (sa < 0 || sb < 0) continue;
      if (std::find(nb[size_t(sa)].begin(), nb[size_t(sa)].end(), uint32_t(sb)) == nb[size_t(sa)].end()) continue;
      add_bond(bead[a].at("SC1"), bead[b].at("SC1"), L["disulfide_constraint"].number(), kc, "disulfide constraint");
      ++rep.disulfides;
    }
  for (const auto& b : topo->bonds) out.bonds.push_back({b.i, b.j, 1});
  topo->natoms = out.atoms.size();
  out.topology = topo;
  out.bonds_from_file = true;
  out.has_charges = true;
  out.has_mol = true;
  // one type entry per bead type (72 until typed: the force field gives each bead type its mass)
  std::map<std::string, int> tid;
  for (auto& a : out.atoms) {
    auto [it, fresh] = tid.emplace(a.name, int(tid.size()) + 1);
    a.type = it->second;
    if (fresh) {
      TypeInfo t;
      t.type = it->second;
      t.label = a.name;
      t.mass = 72.0;
      out.types.push_back(t);
    }
  }
  rep.residues = int(nres);
  rep.beads = int(out.atoms.size());
  rep.chains = nres ? res.back().chain + 1 : 0;
  rep.residue_names = mres;
  if (rep_out) *rep_out = std::move(rep);
  return out;
}

System gromacs_molecule(const Json& mol, const std::function<double(const std::string&)>& type_mass, double constraint_kj, uint64_t seed) {
  BeadMolecule m;
  std::vector<double> len;
  const auto& atoms = mol["atoms"].items();
  for (const auto& a : atoms) {
    m.type.push_back(a["type"].str());
    m.charge.push_back(a["charge"].number());
  }
  auto topo = std::make_shared<ExplicitTopology>();
  topo->source = "Martini 3 molecule";
  const size_t n = atoms.size();
  topo->masses.assign(n, std::nan(""));
  bool any_mass = false;
  for (size_t i = 0; i < n; ++i)
    if (atoms[i].has("mass") && !atoms[i]["mass"].is_null()) topo->masses[i] = atoms[i]["mass"].number(), any_mass = true;
  if (!any_mass) topo->masses.clear();
  auto u = [](const Json& v) { return uint32_t(v.number()); };
  for (const auto& b : mol["bonds"].items()) {
    if (int(b[2].number()) != 1) throw std::invalid_argument("bond function " + std::to_string(int(b[2].number())) + " is not handled");
    topo->bonds.push_back({u(b[0]), u(b[1]), b[4].number() / (2 * kKJ * 100), b[3].number() * 10, "bond"});
    m.bonds.push_back({int(b[0].number()), int(b[1].number())});
    len.push_back(b[3].number() * 10);
  }
  for (const auto& c : mol["constraints"].items()) {
    topo->bonds.push_back({u(c[0]), u(c[1]), constraint_kj / (2 * kKJ * 100), c[3].number() * 10, "constraint"});
    m.bonds.push_back({int(c[0].number()), int(c[1].number())});
    len.push_back(c[3].number() * 10);
  }
  for (const auto& a : mol["angles"].items()) {
    const int f = int(a[3].number());
    const double th = a[4].number() * kDeg, k = a[5].number();
    if (f == 1) topo->angles.push_back({u(a[0]), u(a[1]), u(a[2]), 0, k / (2 * kKJ), th, "angle"});
    else if (f == 2) topo->angles.push_back({u(a[0]), u(a[1]), u(a[2]), 1, k / (2 * kKJ), th, "angle"});
    else if (f == 10) topo->angles.push_back({u(a[0]), u(a[1]), u(a[2]), 5, k / (2 * kKJ), th, "angle"});
    else throw std::invalid_argument("angle function " + std::to_string(f) + " is not handled");
  }
  for (const auto& d : mol["dihedrals"].items()) {
    const int f = int(d[4].number());
    const uint32_t a = u(d[0]), b = u(d[1]), c = u(d[2]), e = u(d[3]);
    if (f == 1 || f == 9 || f == 4) {
      const int mult = d.size() > 7 ? int(d[7].number()) : 1;
      topo->dihedrals.push_back({a, b, c, e, f == 4 ? 4 : 1, d[6].number() / kKJ, d[5].number() * kDeg, mult, "dihedral"});
    } else if (f == 2) {
      topo->dihedrals.push_back({a, b, c, e, 2, d[6].number() / (2 * kKJ), d[5].number() * kDeg, 0, "improper"});
    } else {
      throw std::invalid_argument("dihedral function " + std::to_string(f) + " is not handled");
    }
  }
  for (const auto& x : mol["exclusions"].items())
    for (size_t k = 1; k < x.size(); ++k) topo->exclusions.push_back({u(x[0]), u(x[k])});
  for (const auto& v : mol["vsites"].items()) {
    ExplicitTopology::VSite vs;
    vs.site = u(v["site"]);
    for (const auto& a : v["from"].items()) vs.from.push_back(u(a));
    if (v.has("weights") && !v["weights"].is_null())
      for (const auto& w : v["weights"].items()) vs.w.push_back(w.number());
    topo->vsites.push_back(vs);
  }
  BeadBuildOptions o;
  o.seed = seed;
  o.mass = type_mass;
  System s = build_bead_graph(m, len, o);
  // virtual sites where their atoms put them (weights, or the centre of mass)
  for (const auto& v : topo->vsites) {
    Vec3 c{0, 0, 0};
    double wt = 0;
    for (size_t k = 0; k < v.from.size(); ++k) {
      const double w = !v.w.empty() ? v.w[k] : (!topo->masses.empty() && !std::isnan(topo->masses[v.from[k]]) ? topo->masses[v.from[k]] : type_mass(m.type[v.from[k]]));
      c = c + s.atoms[v.from[k]].pos * w;
      wt += w;
    }
    if (wt != 0) s.atoms[v.site].pos = c * (1 / wt);
  }
  topo->natoms = n;
  s.topology = topo;
  return s;
}

std::string martini_itp(const System& s, const std::string& data_path) {
  if (!s.topology) throw std::invalid_argument("the structure has no explicit topology");
  const double kc = model(data_path)["links"]["constraint_k"].number() / (2 * kKJ * 100);
  std::ostringstream o;
  char b[256];
  o << "[ moleculetype ]\nprotein 1\n\n[ atoms ]\n";
  for (size_t i = 0; i < s.atoms.size(); ++i) {
    const auto& a = s.atoms[i];
    std::snprintf(b, sizeof b, "%4zu %-4s %4lld %-4s %-4s %4zu %6.2f\n", i + 1, a.name.c_str(), static_cast<long long>(a.resid), a.resname.c_str(), "B",
                  i + 1, a.charge);
    o << b;
  }
  o << "\n[ bonds ]\n";
  for (const auto& x : s.topology->bonds)
    if (std::fabs(x.k - kc) > 1e-9) {
      std::snprintf(b, sizeof b, "%4u %4u 1 %.4f %.1f ; %s\n", x.i + 1, x.j + 1, x.r0 / 10, x.k * 2 * kKJ * 100, x.group.c_str());
      o << b;
    }
  o << "\n[ constraints ]\n";
  for (const auto& x : s.topology->bonds)
    if (std::fabs(x.k - kc) <= 1e-9) {
      std::snprintf(b, sizeof b, "%4u %4u 1 %.4f ; %s\n", x.i + 1, x.j + 1, x.r0 / 10, x.group.c_str());
      o << b;
    }
  o << "\n[ angles ]\n";
  for (const auto& x : s.topology->angles) {
    std::snprintf(b, sizeof b, "%4u %4u %4u 2 %.1f %.1f ; %s\n", x.i + 1, x.j + 1, x.k + 1, x.theta0 / kDeg, x.kt * 2 * kKJ, x.group.c_str());
    o << b;
  }
  o << "\n[ dihedrals ]\n";
  for (const auto& x : s.topology->dihedrals) {
    if (x.form == 1) std::snprintf(b, sizeof b, "%4u %4u %4u %4u 1 %.1f %.1f %d ; %s\n", x.i + 1, x.j + 1, x.k + 1, x.l + 1, x.phi0 / kDeg, x.kd * kKJ, x.n, x.group.c_str());
    else std::snprintf(b, sizeof b, "%4u %4u %4u %4u 2 %.1f %.1f ; %s\n", x.i + 1, x.j + 1, x.k + 1, x.l + 1, x.phi0 / kDeg, x.kd * 2 * kKJ, x.group.c_str());
    o << b;
  }
  return o.str();
}

}  // namespace caps
