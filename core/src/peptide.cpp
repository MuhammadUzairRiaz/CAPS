// CAPS peptide builder: sequence and secondary structure to an all-atom peptide (see caps/peptide.hpp).
#include "caps/peptide.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <map>
#include <random>
#include <sstream>
#include <stdexcept>

#include "caps/elements.hpp"
#include "caps/molecule.hpp"

namespace caps {

namespace {

constexpr double kDeg = 3.14159265358979323846 / 180.0;
Vec3 unitv(const Vec3& v) { const double n = norm(v); return n > 0 ? v * (1.0 / n) : Vec3{1, 0, 0}; }

// NeRF: d from a-b-c with |cd| = bond, angle bcd = ang, torsion abcd = tor (degrees, IUPAC sign).
Vec3 place(const Vec3& a, const Vec3& b, const Vec3& c, double bond, double ang, double tor) {
  const Vec3 bc = unitv(c - b);
  const Vec3 n = unitv(cross(b - a, bc));
  const Vec3 m = cross(n, bc);
  const double A = ang * kDeg, T = tor * kDeg;
  const Vec3 d2{-bond * std::cos(A), bond * std::sin(A) * std::cos(T), bond * std::sin(A) * std::sin(T)};
  return c + bc * d2[0] + m * d2[1] + n * d2[2];
}

// One residue: its codes and the side chain as SMILES written after N[C@@H]( for L, neutral and charged forms.
struct ResidueDef {
  char code;
  const char* name;
  const char* side;      // neutral
  const char* charged;   // the ionised form, or nullptr
  double pka;            // of the ionisable group (0: none)
  bool acid;             // ionised above pKa (acids) or below (bases)
};

// L-isoleucine and L-threonine carry their second centre (2S,3S and 2S,3R) as in PubChem's isomeric SMILES.
const ResidueDef kResidues[] = {
    {'A', "ALA", "C", nullptr, 0, false},
    {'R', "ARG", "CCCNC(=N)N", "CCCNC(=[NH2+])N", 12.5, false},
    {'N', "ASN", "CC(=O)N", nullptr, 0, false},
    {'D', "ASP", "CC(=O)O", "CC(=O)[O-]", 3.9, true},
    {'C', "CYS", "CS", "C[S-]", 8.3, true},
    {'Q', "GLN", "CCC(=O)N", nullptr, 0, false},
    {'E', "GLU", "CCC(=O)O", "CCC(=O)[O-]", 4.3, true},
    {'G', "GLY", "", nullptr, 0, false},
    {'H', "HIS", "Cc1c[nH]cn1", "Cc1c[nH]c[nH+]1", 6.0, false},
    {'I', "ILE", "[C@@H](C)CC", nullptr, 0, false},
    {'L', "LEU", "CC(C)C", nullptr, 0, false},
    {'K', "LYS", "CCCCN", "CCCC[NH3+]", 10.5, false},
    {'M', "MET", "CCSC", nullptr, 0, false},
    {'F', "PHE", "Cc1ccccc1", nullptr, 0, false},
    {'P', "PRO", "CCC1", nullptr, 0, false},
    {'S', "SER", "CO", nullptr, 0, false},
    {'T', "THR", "[C@H](O)C", nullptr, 0, false},
    {'W', "TRP", "Cc1c[nH]c2ccccc12", nullptr, 0, false},
    {'Y', "TYR", "Cc1ccc(O)cc1", "Cc1ccc([O-])cc1", 10.1, true},
    {'V', "VAL", "C(C)C", nullptr, 0, false},
};

const ResidueDef* residue_def(char c) {
  c = char(std::toupper(static_cast<unsigned char>(c)));
  for (const auto& r : kResidues) if (r.code == c) return &r;
  return nullptr;
}

std::string side_at(const ResidueDef& r, const PeptideOptions& o) {
  if (!r.charged || o.neutral) return r.side;
  const bool ionised = r.acid ? o.ph > r.pka : o.ph < r.pka;
  return ionised ? r.charged : r.side;
}

// The residue as written in the peptide: N, CA, side chain, C, O (Pro closes its ring onto N with digit 1).
std::string residue_smiles(const ResidueDef& r, const std::string& side, const std::string& n_atom) {
  if (r.code == 'G') return n_atom + "CC(=O)";
  if (r.code == 'P') return (n_atom == "N" ? std::string("N1") : n_atom.substr(0, n_atom.size() - 1) + "]1") + "[C@@H](CCC1)C(=O)";
  return n_atom + "[C@@H](" + side + ")C(=O)";
}

// Hydrogen bond lengths by the parent element.
double h_bond(int z) { return z == 6 ? 1.09 : z == 7 ? 1.01 : z == 8 ? 0.96 : z == 16 ? 1.34 : 1.0; }

// Superposes three source points on three target points (frames on a–b and b–c): the map from source to target space.
struct Frame {
  Vec3 o, e1, e2, e3;
  static Frame of(const Vec3& a, const Vec3& b, const Vec3& c) {
    Frame f;
    f.o = b;
    f.e1 = unitv(c - b);
    const Vec3 u = a - b;
    f.e2 = unitv(u - f.e1 * dot(u, f.e1));
    f.e3 = cross(f.e1, f.e2);
    return f;
  }
};
Vec3 carry(const Frame& from, const Frame& to, const Vec3& p) {
  const Vec3 d = p - from.o;
  const double x = dot(d, from.e1), y = dot(d, from.e2), z = dot(d, from.e3);
  return to.o + to.e1 * x + to.e2 * y + to.e3 * z;
}

struct Template {
  MolGraph g;
  std::vector<Vec3> pos;
  std::vector<std::vector<int>> h_of;   // hydrogens of each written atom
};

std::vector<std::vector<int>> hydrogens_of(const MolGraph& g) {
  std::vector<std::vector<int>> h(g.atoms.size());
  for (const auto& b : g.bonds) {
    if (g.atoms[size_t(b.b)].element == 1 && b.b >= g.heavy) h[size_t(b.a)].push_back(b.b);
    else if (g.atoms[size_t(b.a)].element == 1 && b.a >= g.heavy) h[size_t(b.b)].push_back(b.a);
  }
  return h;
}

const Template& residue_template(const std::string& smiles) {
  static std::map<std::string, Template> cache;
  auto it = cache.find(smiles);
  if (it != cache.end()) return it->second;
  Template t;
  t.g = parse_smiles(smiles);
  add_hydrogens(t.g);
  EmbedOptions eo;
  eo.seed = 7;
  t.pos = embed(t.g, eo);
  std::vector<std::string> notes;
  if (auto ff = molecule_forcefield(t.g, "uff", "gasteiger", notes, nullptr, &t.pos)) minimise_molecule(t.g, ff, 0.5, t.pos);
  t.h_of = hydrogens_of(t.g);
  return cache.emplace(smiles, std::move(t)).first->second;
}

}  // namespace

std::string residue_name(char code) {
  const auto* r = residue_def(code);
  return r ? r->name : "";
}

std::array<double, 3> structure_torsions(const PeptideOptions& o, char ss) {
  switch (std::toupper(static_cast<unsigned char>(ss))) {
    case 'H': return o.helix;
    case 'E': return o.strand;
    case 'P': return o.ppii;
    default: return {-70, 140, 180};
  }
}

std::string parse_fasta(const std::string& text) {
  std::istringstream in(text);
  std::string line, seq;
  bool header_seen = false;
  while (std::getline(in, line)) {
    if (!line.empty() && line[0] == '>') {
      if (header_seen && !seq.empty()) break;   // the first record only
      header_seen = true;
      continue;
    }
    if (!line.empty() && line[0] == ';') continue;
    for (char c : line)
      if (std::isalpha(static_cast<unsigned char>(c))) seq += char(std::toupper(static_cast<unsigned char>(c)));
  }
  return seq;
}

System build_peptide(const PeptideOptions& o, PeptideReport* report) {
  // the sequence
  std::string seq;
  for (char c : o.sequence) {
    if (std::isspace(static_cast<unsigned char>(c)) || std::isdigit(static_cast<unsigned char>(c)) || c == '-' || c == '*') continue;
    if (!residue_def(c)) throw std::invalid_argument(std::string("'") + c + "' is not one of the 20 amino-acid codes");
    seq += char(std::toupper(static_cast<unsigned char>(c)));
  }
  if (seq.empty()) throw std::invalid_argument("the sequence is empty");
  if (seq.size() > 2000) throw std::invalid_argument("more than 2000 residues");
  const size_t n = seq.size();
  std::string ss;
  for (size_t i = 0; i < n; ++i) {
    char c = i < o.structure.size() ? char(std::toupper(static_cast<unsigned char>(o.structure[i]))) : 'C';
    ss += (c == 'H' || c == 'E' || c == 'P') ? c : 'C';
  }
  const bool ace = o.n_term == "ACE", nh2 = o.n_term == "NH2", nme = o.c_term == "NME", cooh = o.c_term == "COOH";

  // the peptide as SMILES, remembering where each residue's atoms are written
  struct Res { const ResidueDef* def; std::string side; int N, CA, C, O, side0, nside; };
  std::vector<Res> res(n);
  std::string smi;
  int at = 0;
  int ace_c = -1, ace_o = -1, ace_me = -1, oxt = -1, nme_n = -1, nme_c = -1;
  if (ace) { smi = "CC(=O)"; ace_me = 0; ace_c = 1; ace_o = 2; at = 3; }
  for (size_t i = 0; i < n; ++i) {
    Res& r = res[i];
    r.def = residue_def(seq[i]);
    r.side = side_at(*r.def, o);
    std::string n_atom = "N";
    if (i == 0 && !ace && !nh2) n_atom = r.def->code == 'P' ? "[NH2+]" : "[NH3+]";
    smi += residue_smiles(*r.def, r.side, n_atom);
    // side-chain atoms: those of the template between CA and C
    const Template& t = residue_template(residue_smiles(*r.def, r.side, "N") + "O");
    r.nside = t.g.heavy - 5;
    r.N = at; r.CA = at + 1; r.side0 = at + 2; r.C = at + 2 + r.nside; r.O = r.C + 1;
    at = r.O + 1;
  }
  if (nme) { smi += "NC"; nme_n = at; nme_c = at + 1; at += 2; }
  else { smi += cooh ? "O" : "[O-]"; oxt = at; at += 1; }

  MolGraph g = parse_smiles(smi);
  if (g.heavy != at) throw std::logic_error("peptide SMILES bookkeeping: " + std::to_string(g.heavy) + " atoms written, " + std::to_string(at) + " counted");
  add_hydrogens(g);
  const auto h_of = hydrogens_of(g);
  std::vector<std::vector<int>> nb(g.atoms.size());
  for (const auto& b : g.bonds) nb[size_t(b.a)].push_back(b.b), nb[size_t(b.b)].push_back(b.a);
  std::vector<Vec3> pos(g.atoms.size(), Vec3{0, 0, 0});
  std::vector<char> placed(g.atoms.size(), 0);
  auto put = [&](int i, const Vec3& p) { pos[size_t(i)] = p; placed[size_t(i)] = 1; };

  // backbone torsions: coil residues drawn from the αR, β and PPII basins; proline keeps φ −65°
  std::mt19937_64 rng(o.seed);
  std::normal_distribution<double> jitter(0.0, 12.0);
  std::uniform_real_distribution<double> pick(0.0, 1.0);
  auto torsions_of = [&](size_t i) {
    std::array<double, 3> t = structure_torsions(o, ss[i]);
    if (ss[i] == 'C') {
      const double u = pick(rng);
      t = u < 0.35 ? std::array<double, 3>{-63, -43, 180} : u < 0.7 ? std::array<double, 3>{-120, 130, 180} : std::array<double, 3>{-75, 145, 180};
      t[0] += jitter(rng), t[1] += jitter(rng);
    }
    if (res[i].def->code == 'P') t[0] = -65;
    return t;
  };

  // N, CA, C, O by NeRF; a coil residue is redrawn (up to 30 times) when its CA comes within 3.6 Å of an earlier one
  std::vector<std::array<double, 3>> tor(n);
  const double bNCA = 1.458, bCAC = 1.525, bCN = 1.329, bCO = 1.231;
  const double aNCAC = 111.2, aCACN = 116.2, aCNCA = 121.7, aCACO = 120.5;
  for (size_t i = 0; i < n; ++i) {
    for (int attempt = 0;; ++attempt) {
      tor[i] = torsions_of(i);
      if (i == 0) {
        put(res[0].N, {0, 0, 0});
        put(res[0].CA, {bNCA, 0, 0});
        put(res[0].C, pos[size_t(res[0].CA)] + Vec3{-bCAC * std::cos(aNCAC * kDeg), bCAC * std::sin(aNCAC * kDeg), 0});
        break;
      }
      const Res &p = res[i - 1], &r = res[i];
      put(r.N, place(pos[size_t(p.N)], pos[size_t(p.CA)], pos[size_t(p.C)], bCN, aCACN, tor[i - 1][1]));
      put(r.CA, place(pos[size_t(p.CA)], pos[size_t(p.C)], pos[size_t(r.N)], bNCA, aCNCA, tor[i - 1][2]));
      put(r.C, place(pos[size_t(p.C)], pos[size_t(r.N)], pos[size_t(r.CA)], bCAC, aNCAC, tor[i][0]));
      bool clash = false;
      for (size_t j = 0; j + 2 < i && !clash; ++j) clash = norm(pos[size_t(res[j].CA)] - pos[size_t(r.CA)]) < 3.6;
      if (!clash || ss[i] != 'C' || attempt >= 30) break;
      // redraw this residue's φ and the previous residue's ψ
      auto t = torsions_of(i - 1);
      tor[i - 1][1] = t[1];
    }
  }
  for (size_t i = 0; i < n; ++i) {
    const Res& r = res[i];
    const Vec3 &N = pos[size_t(r.N)], &CA = pos[size_t(r.CA)], &C = pos[size_t(r.C)];
    put(r.O, place(N, CA, C, bCO, aCACO, tor[i][1] + 180));
    // side chain and hydrogens on CA and the side chain from the embedded free amino acid, superposed on N, CA, C
    const Template& t = residue_template(residue_smiles(*r.def, r.side, "N") + "O");
    const int tC = 2 + r.nside;
    const Frame from = Frame::of(t.pos[0], t.pos[1], t.pos[size_t(tC)]), to = Frame::of(N, CA, C);
    auto take = [&](int tatom, int atom) {
      if (tatom != 0 && tatom != 1 && tatom != tC) put(atom, carry(from, to, t.pos[size_t(tatom)]));
      const auto& th = t.h_of[size_t(tatom)];
      const auto& ph = h_of[size_t(atom)];
      for (size_t k = 0; k < th.size() && k < ph.size(); ++k) put(ph[k], carry(from, to, t.pos[size_t(th[k])]));
    };
    take(1, r.CA);
    for (int k = 0; k < r.nside; ++k) take(2 + k, r.side0 + k);
  }
  // termini and caps
  const Res &first = res.front(), &last = res.back();
  if (ace) {
    const Vec3 &N = pos[size_t(first.N)], &CA = pos[size_t(first.CA)], &C = pos[size_t(first.C)];
    put(ace_c, place(C, CA, N, bCN, aCNCA, tor[0][0]));
    put(ace_o, place(CA, N, pos[size_t(ace_c)], bCO, 123.0, 0));
    put(ace_me, place(CA, N, pos[size_t(ace_c)], 1.52, aCACN, 180));
  }
  {
    const Vec3 &N = pos[size_t(last.N)], &CA = pos[size_t(last.CA)], &C = pos[size_t(last.C)];
    if (nme) {
      put(nme_n, place(N, CA, C, bCN, aCACN, tor[n - 1][1]));
      put(nme_c, place(CA, C, pos[size_t(nme_n)], 1.46, aCNCA, 180));
    } else {
      put(oxt, place(N, CA, C, 1.25, 117.0, tor[n - 1][1]));
      if (!cooh) pos[size_t(last.O)] = place(N, CA, C, 1.25, 117.0, tor[n - 1][1] + 180);   // carboxylate: both C–O 1.25 Å
    }
  }
  // every hydrogen not yet placed, from its parent's geometry
  for (size_t p = 0; p < size_t(g.heavy); ++p) {
    std::vector<int> todo;
    for (int h : h_of[p]) if (!placed[size_t(h)]) todo.push_back(h);
    if (todo.empty()) continue;
    std::vector<int> heavy;
    for (int q : nb[p]) if (g.atoms[size_t(q)].element != 1) heavy.push_back(q);
    const Vec3 P = pos[p];
    const double b = h_bond(g.atoms[p].element);
    const size_t k = todo.size();
    if (heavy.size() >= 2) {
      Vec3 s{0, 0, 0};
      for (int q : heavy) s = s + unitv(P - pos[size_t(q)]);
      if (k == 1) { put(todo[0], P + unitv(s) * b); continue; }
      // two hydrogens on a centre with two heavy neighbours: tetrahedral, either side of their plane
      const Vec3 nrm = unitv(cross(pos[size_t(heavy[0])] - P, pos[size_t(heavy[1])] - P));
      const Vec3 bis = unitv(s);
      put(todo[0], P + (bis * std::cos(54.75 * kDeg) + nrm * std::sin(54.75 * kDeg)) * b);
      put(todo[1], P + (bis * std::cos(54.75 * kDeg) - nrm * std::sin(54.75 * kDeg)) * b);
      continue;
    }
    if (heavy.empty()) continue;
    const int A = heavy[0];
    int R = -1;   // a reference atom across the bond for the torsions
    for (int q : nb[size_t(A)]) if (size_t(q) != p && placed[size_t(q)] && g.atoms[size_t(q)].element != 1) { R = q; break; }
    const Vec3 Rp = R >= 0 ? pos[size_t(R)] : pos[size_t(A)] + Vec3{0, 1, 0.3};
    // sp3: staggered on the reference (NH3+, CH3: 60 180 300; NH2: 60 300, the lone pair anti; OH, SH: anti)
    static const double k1[] = {180}, k2[] = {60, 300}, k3[] = {60, 180, 300};
    const double* ts = k == 1 ? k1 : k == 2 ? k2 : k3;
    for (size_t m = 0; m < k && m < 3; ++m) put(todo[m], place(Rp, pos[size_t(A)], P, b, 109.5, ts[m]));
  }

  // clean-up
  std::vector<std::string> notes;
  if (o.cleanup) {
    std::string why;
    if (auto ff = molecule_forcefield(g, "uff", "gasteiger", notes, nullptr, &pos)) {
      std::vector<Vec3> p = pos;
      double e = 0;
      if (minimise_molecule(g, ff, o.ftol, p, &e, &why)) {
        pos = p;
        char b[120];
        std::snprintf(b, sizeof b, "UFF clean-up: E = %.1f kcal/mol (ftol %.2g)", e, o.ftol);
        notes.push_back(b);
      } else {
        notes.push_back("UFF clean-up skipped: " + why);
      }
    }
  }
  const auto chk = chirality_check(g, pos);
  const int wrong = int(std::count(chk.begin(), chk.end(), -1));
  if (wrong) notes.push_back(std::to_string(wrong) + " stereocentre(s) inverted");

  // the system: PDB-style atom names, three-letter residue names
  System s = molecule_system(g, pos);
  s.title = "peptide " + (seq.size() > 24 ? seq.substr(0, 24) + "…" : seq);
  s.source_format = "peptide";
  s.has_mol = true;
  s.unwrapped = true;
  static const char* kGreek = "ABGDEZH";
  std::vector<std::string> names(g.atoms.size());
  std::vector<std::string> resname(g.atoms.size(), "");
  std::vector<int64_t> resid(g.atoms.size(), 0);
  const int64_t first_res = ace ? 2 : 1;   // an acetyl cap is residue 1
  for (size_t i = 0; i < n; ++i) {
    const Res& r = res[i];
    names[size_t(r.N)] = "N"; names[size_t(r.CA)] = "CA"; names[size_t(r.C)] = "C"; names[size_t(r.O)] = "O";
    for (int a : {r.N, r.CA, r.C, r.O}) resname[size_t(a)] = r.def->name, resid[size_t(a)] = first_res + int64_t(i);
    // side chain by bonds from CA: B, G, D, E, Z, H; numbered where an element repeats at one depth
    std::map<int, int> depth{{r.CA, 0}};
    std::vector<int> frontier{r.CA};
    while (!frontier.empty()) {
      std::vector<int> next;
      for (int a : frontier)
        for (int q : nb[size_t(a)])
          if (q >= r.side0 && q < r.side0 + r.nside && !depth.count(q)) { depth[q] = depth[a] + 1; next.push_back(q); }
      frontier = next;
    }
    std::map<std::pair<int, int>, std::vector<int>> groups;
    for (int k = 0; k < r.nside; ++k) {
      const int a = r.side0 + k;
      groups[{depth.count(a) ? depth[a] : 6, g.atoms[size_t(a)].element}].push_back(a);
      resname[size_t(a)] = r.def->name;
      resid[size_t(a)] = first_res + int64_t(i);
    }
    for (const auto& [key, atoms] : groups)
      for (size_t k = 0; k < atoms.size(); ++k)
        names[size_t(atoms[k])] = std::string(element(key.second).symbol) + kGreek[std::min(key.first, 6)] + (atoms.size() > 1 ? std::to_string(k + 1) : "");
  }
  if (oxt >= 0) { names[size_t(oxt)] = "OXT"; resname[size_t(oxt)] = res.back().def->name; resid[size_t(oxt)] = first_res + int64_t(n) - 1; }
  if (ace) { names[size_t(ace_me)] = "CH3"; names[size_t(ace_c)] = "C"; names[size_t(ace_o)] = "O"; for (int a : {ace_me, ace_c, ace_o}) resname[size_t(a)] = "ACE", resid[size_t(a)] = 1; }
  if (nme) {
    names[size_t(nme_n)] = "N"; names[size_t(nme_c)] = "CH3";
    for (int a : {nme_n, nme_c}) resname[size_t(a)] = "NME", resid[size_t(a)] = first_res + int64_t(n);
  }
  for (size_t p = 0; p < size_t(g.heavy); ++p) {
    const auto& hs = h_of[p];
    std::string base = names[p];
    if (base == "N") base = "";
    else if (base.size() > 1) base = base.substr(1);
    else base = "";
    for (size_t k = 0; k < hs.size(); ++k) {
      names[size_t(hs[k])] = "H" + base + (hs.size() > 1 ? std::to_string(k + 1) : "");
      resname[size_t(hs[k])] = resname[p];
      resid[size_t(hs[k])] = resid[p];
    }
  }
  std::map<int, int> type_of;
  for (size_t i = 0; i < s.atoms.size(); ++i) {
    auto& a = s.atoms[i];
    if (!names[i].empty()) a.name = names[i];
    a.resname = resname[i].empty() ? "UNK" : resname[i];
    a.resid = resid[i];
    auto it = type_of.find(a.element);
    if (it == type_of.end()) {
      it = type_of.emplace(a.element, int(type_of.size()) + 1).first;
      TypeInfo ti;
      ti.type = it->second;
      ti.mass = element(a.element).mass;
      ti.label = element(a.element).symbol;
      s.types.push_back(ti);
    }
    a.type = it->second;
  }

  const MolInfo info = molecule_info(g);
  if (report) {
    report->smiles = smi;
    report->formula = info.formula;
    report->structure = ss;
    report->residues = int(n);
    report->charge = info.charge;
    report->atoms = s.atoms.size();
    report->mass = info.mass;
    report->notes = notes;
  }
  char b[200];
  std::snprintf(b, sizeof b, "%zu residues · %zu atoms · %s · %.1f g/mol · charge %+d", n, s.atoms.size(), info.formula.c_str(), info.mass, info.charge);
  s.notes.insert(s.notes.begin(), b);
  for (const auto& note : notes) s.notes.push_back(note);
  return s;
}

}  // namespace caps
