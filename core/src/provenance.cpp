#include "caps/provenance.hpp"

#include <algorithm>
#include <ctime>
#include <fstream>
#include <map>
#include <set>

namespace caps {

namespace {

Json kv_json(const KeyValues& kv) {
  Json o = Json::object();
  for (const auto& [k, v] : kv) o[k] = v;
  return o;
}
KeyValues kv_from(const Json& j, const char* key) {
  KeyValues out;
  if (!j.has(key) || !j[key].is_object()) return out;
  for (const auto& [k, v] : j[key].members()) out.push_back({k, v.is_string() ? v.str() : v.dump(0)});
  return out;
}

struct Cite { const char* key; const char* entry; };
// The methods CAPS implements, as published.
const Cite kCites[] = {
    {"rappe1992", "@article{rappe1992,\n  author = {Rapp{\\'e}, A. K. and Casewit, C. J. and Colwell, K. S. and Goddard, W. A. and Skiff, W. M.},\n  title = {{UFF}, a full periodic table force field for molecular mechanics and molecular dynamics simulations},\n  journal = {J. Am. Chem. Soc.}, volume = {114}, pages = {10024--10035}, year = {1992}, doi = {10.1021/ja00051a040}\n}"},
    {"rappe1991", "@article{rappe1991,\n  author = {Rapp{\\'e}, A. K. and Goddard, W. A.},\n  title = {Charge equilibration for molecular dynamics simulations},\n  journal = {J. Phys. Chem.}, volume = {95}, pages = {3358--3363}, year = {1991}, doi = {10.1021/j100161a070}\n}"},
    {"wang2004", "@article{wang2004,\n  author = {Wang, J. and Wolf, R. M. and Caldwell, J. W. and Kollman, P. A. and Case, D. A.},\n  title = {Development and testing of a general amber force field},\n  journal = {J. Comput. Chem.}, volume = {25}, pages = {1157--1174}, year = {2004}, doi = {10.1002/jcc.20035}\n}"},
    {"gasteiger1980", "@article{gasteiger1980,\n  author = {Gasteiger, J. and Marsili, M.},\n  title = {Iterative partial equalization of orbital electronegativity---a rapid access to atomic charges},\n  journal = {Tetrahedron}, volume = {36}, pages = {3219--3228}, year = {1980}, doi = {10.1016/0040-4020(80)80168-2}\n}"},
    {"cordero2008", "@article{cordero2008,\n  author = {Cordero, B. and G{\\'o}mez, V. and Platero-Prats, A. E. and Rev{\\'e}s, M. and Echeverr{\\'i}a, J. and Cremades, E. and Barrag{\\'a}n, F. and Alvarez, S.},\n  title = {Covalent radii revisited},\n  journal = {Dalton Trans.}, pages = {2832--2838}, year = {2008}, doi = {10.1039/b801115j}\n}"},
    {"liu1989", "@article{liu1989,\n  author = {Liu, D. C. and Nocedal, J.},\n  title = {On the limited memory {BFGS} method for large scale optimization},\n  journal = {Math. Program.}, volume = {45}, pages = {503--528}, year = {1989}, doi = {10.1007/BF01589116}\n}"},
    {"polak1969", "@article{polak1969,\n  author = {Polak, E. and Ribi{\\`e}re, G.},\n  title = {Note sur la convergence de m{\\'e}thodes de directions conjugu{\\'e}es},\n  journal = {Rev. Fr. Inform. Rech. Op{\\'e}r.}, volume = {3}, pages = {35--43}, year = {1969}\n}"},
    {"bitzek2006", "@article{bitzek2006,\n  author = {Bitzek, E. and Koskinen, P. and G{\\\"a}hler, F. and Moseler, M. and Gumbsch, P.},\n  title = {Structural relaxation made simple},\n  journal = {Phys. Rev. Lett.}, volume = {97}, pages = {170201}, year = {2006}, doi = {10.1103/PhysRevLett.97.170201}\n}"},
    {"auhl2003", "@article{auhl2003,\n  author = {Auhl, R. and Everaers, R. and Grest, G. S. and Kremer, K. and Plimpton, S. J.},\n  title = {Equilibration of long chain polymer melts in computer simulations},\n  journal = {J. Chem. Phys.}, volume = {119}, pages = {12718--12728}, year = {2003}, doi = {10.1063/1.1628670}\n}"},
    {"larsen2011", "@article{larsen2011,\n  author = {Larsen, G. S. and Lin, P. and Hart, K. E. and Colina, C. M.},\n  title = {Molecular simulations of {PIM-1}-like polymers of intrinsic microporosity},\n  journal = {Macromolecules}, volume = {44}, pages = {6944--6951}, year = {2011}, doi = {10.1021/ma200345v}\n}"},
    {"bussi2007", "@article{bussi2007,\n  author = {Bussi, G. and Donadio, D. and Parrinello, M.},\n  title = {Canonical sampling through velocity rescaling},\n  journal = {J. Chem. Phys.}, volume = {126}, pages = {014101}, year = {2007}, doi = {10.1063/1.2408420}\n}"},
    {"bernetti2020", "@article{bernetti2020,\n  author = {Bernetti, M. and Bussi, G.},\n  title = {Pressure control using stochastic cell rescaling},\n  journal = {J. Chem. Phys.}, volume = {153}, pages = {114107}, year = {2020}, doi = {10.1063/5.0020514}\n}"},
    {"berendsen1984", "@article{berendsen1984,\n  author = {Berendsen, H. J. C. and Postma, J. P. M. and van Gunsteren, W. F. and DiNola, A. and Haak, J. R.},\n  title = {Molecular dynamics with coupling to an external bath},\n  journal = {J. Chem. Phys.}, volume = {81}, pages = {3684--3690}, year = {1984}, doi = {10.1063/1.448118}\n}"},
    {"swope1982", "@article{swope1982,\n  author = {Swope, W. C. and Andersen, H. C. and Berens, P. H. and Wilson, K. R.},\n  title = {A computer simulation method for the calculation of equilibrium constants for the formation of physical clusters of molecules: application to small water clusters},\n  journal = {J. Chem. Phys.}, volume = {76}, pages = {637--649}, year = {1982}, doi = {10.1063/1.442716}\n}"},
    {"essmann1995", "@article{essmann1995,\n  author = {Essmann, U. and Perera, L. and Berkowitz, M. L. and Darden, T. and Lee, H. and Pedersen, L. G.},\n  title = {A smooth particle mesh {E}wald method},\n  journal = {J. Chem. Phys.}, volume = {103}, pages = {8577--8593}, year = {1995}, doi = {10.1063/1.470117}\n}"},
    {"fennell2006", "@article{fennell2006,\n  author = {Fennell, C. J. and Gezelter, J. D.},\n  title = {Is the {E}wald summation still necessary? {P}airwise alternatives to the accepted standard for long-range electrostatics},\n  journal = {J. Chem. Phys.}, volume = {124}, pages = {234104}, year = {2006}, doi = {10.1063/1.2206581}\n}"},
    {"martinez2009", "@article{martinez2009,\n  author = {Mart{\\'i}nez, L. and Andrade, R. and Birgin, E. G. and Mart{\\'i}nez, J. M.},\n  title = {{PACKMOL}: a package for building initial configurations for molecular dynamics simulations},\n  journal = {J. Comput. Chem.}, volume = {30}, pages = {2157--2164}, year = {2009}, doi = {10.1002/jcc.21224}\n}"},
    {"matsumoto1998", "@article{matsumoto1998,\n  author = {Matsumoto, M. and Nishimura, T.},\n  title = {Mersenne twister: a 623-dimensionally equidistributed uniform pseudo-random number generator},\n  journal = {ACM Trans. Model. Comput. Simul.}, volume = {8}, pages = {3--30}, year = {1998}, doi = {10.1145/272991.272995}\n}"},
    {"abascal2005", "@article{abascal2005,\n  author = {Abascal, J. L. F. and Vega, C.},\n  title = {A general purpose model for the condensed phases of water: {TIP4P/2005}},\n  journal = {J. Chem. Phys.}, volume = {123}, pages = {234505}, year = {2005}, doi = {10.1063/1.2121687}\n}"},
    {"jorgensen1983", "@article{jorgensen1983,\n  author = {Jorgensen, W. L. and Chandrasekhar, J. and Madura, J. D. and Impey, R. W. and Klein, M. L.},\n  title = {Comparison of simple potential functions for simulating liquid water},\n  journal = {J. Chem. Phys.}, volume = {79}, pages = {926--935}, year = {1983}, doi = {10.1063/1.445869}\n}"},
    {"berendsen1987", "@article{berendsen1987,\n  author = {Berendsen, H. J. C. and Grigera, J. R. and Straatsma, T. P.},\n  title = {The missing term in effective pair potentials},\n  journal = {J. Phys. Chem.}, volume = {91}, pages = {6269--6271}, year = {1987}, doi = {10.1021/j100308a038}\n}"},
    {"hall1981", "@article{hall1981,\n  author = {Hall, S. R.},\n  title = {Space-group notation with an explicit origin},\n  journal = {Acta Cryst. A}, volume = {37}, pages = {517--525}, year = {1981}, doi = {10.1107/S0567739481001228}\n}"},
    {"engh1991", "@article{engh1991,\n  author = {Engh, R. A. and Huber, R.},\n  title = {Accurate bond and angle parameters for {X}-ray protein structure refinement},\n  journal = {Acta Cryst. A}, volume = {47}, pages = {392--400}, year = {1991}, doi = {10.1107/S0108767391001071}\n}"},
    {"parsons2005", "@article{parsons2005,\n  author = {Parsons, J. and Holmes, J. B. and Rojas, J. M. and Tsai, J. and Strauss, C. E. M.},\n  title = {Practical conversion from torsion space to {C}artesian space for in silico protein synthesis},\n  journal = {J. Comput. Chem.}, volume = {26}, pages = {1063--1068}, year = {2005}, doi = {10.1002/jcc.20237}\n}"},
};

}  // namespace

std::string now_iso() {
  const std::time_t t = std::time(nullptr);
  char buf[32];
  std::strftime(buf, sizeof buf, "%Y-%m-%dT%H:%M:%SZ", std::gmtime(&t));
  return buf;
}

Json manifest_json(const Manifest& m) {
  Json j = Json::object();
  j["schema"] = "caps-manifest/1.0";
  j["generator"] = m.generator;
  j["deterministic"] = m.deterministic;
  Json in = Json::array();
  for (const auto& [name, sha] : m.inputs) {
    Json o = Json::object();
    o["name"] = name;
    o["sha256"] = sha;
    in.push_back(std::move(o));
  }
  j["inputs"] = std::move(in);
  Json steps = Json::array();
  for (const auto& s : m.steps) {
    Json o = Json::object();
    o["engine"] = s.engine;
    o["summary"] = s.summary;
    o["params"] = kv_json(s.params);
    if (!s.rng.empty()) o["rng"] = s.rng;
    Json c = Json::array();
    for (const auto& k : s.cites) c.push_back(k);
    o["cites"] = std::move(c);
    if (!s.approximations.empty()) o["approximations"] = kv_json(s.approximations);
    o["time"] = s.time;
    steps.push_back(std::move(o));
  }
  j["steps"] = std::move(steps);
  Json ap = Json::object();
  for (const auto& [k, v] : approximations(m)) ap[k] = v;
  j["approximations"] = std::move(ap);
  return j;
}

Manifest manifest_from_json(const Json& j) {
  Manifest m;
  m.generator = j.text("generator", m.generator);
  m.deterministic = !j.has("deterministic") || j["deterministic"].boolean();
  if (j.has("inputs") && j["inputs"].is_array())
    for (const auto& o : j["inputs"].items()) m.inputs.push_back({o.text("name"), o.text("sha256")});
  if (j.has("steps") && j["steps"].is_array())
    for (const auto& o : j["steps"].items()) {
      ProvStep s;
      s.engine = o.text("engine");
      s.summary = o.text("summary");
      s.params = kv_from(o, "params");
      s.rng = o.text("rng");
      if (o.has("cites") && o["cites"].is_array())
        for (const auto& c : o["cites"].items()) s.cites.push_back(c.str());
      s.approximations = kv_from(o, "approximations");
      s.time = o.text("time");
      m.steps.push_back(std::move(s));
    }
  return m;
}

std::string sidecar_path(const std::string& data_path) { return data_path + ".provenance.json"; }

void write_manifest(const Manifest& m, const std::string& data_path) {
  std::ofstream f(sidecar_path(data_path));
  if (!f) throw std::runtime_error("cannot write " + sidecar_path(data_path));
  f << manifest_json(m).dump(2) << "\n";
}

std::optional<Manifest> read_manifest(const std::string& data_path) {
  std::ifstream f(sidecar_path(data_path));
  if (!f) return std::nullopt;
  std::string text((std::istreambuf_iterator<char>(f)), std::istreambuf_iterator<char>());
  try {
    const Json j = Json::parse(text);
    if (j.text("schema") != "caps-manifest/1.0") return std::nullopt;
    return manifest_from_json(j);
  } catch (...) {
    return std::nullopt;
  }
}

KeyValues approximations(const Manifest& m) {
  static const char* order[] = {"van der Waals", "Electrostatics", "Constraints", "Precision", "Estimated parameters"};
  std::map<std::string, std::string> latest;
  for (const auto& s : m.steps)
    for (const auto& [k, v] : s.approximations) latest[k] = v;
  KeyValues out;
  for (const char* k : order)
    if (auto it = latest.find(k); it != latest.end()) out.push_back({k, it->second});
  for (const auto& [k, v] : latest)
    if (std::find(std::begin(order), std::end(order), k) == std::end(order)) out.push_back({k, v});
  return out;
}

ManifestDiff compare(const Manifest& a, const Manifest& b) {
  ManifestDiff d;
  d.steps_a = int(a.steps.size());
  d.steps_b = int(b.steps.size());
  d.same_inputs = a.inputs == b.inputs;
  d.same_generator = a.generator == b.generator;
  if (!d.same_generator) d.notes.push_back("generator: " + a.generator + " vs " + b.generator);
  if (!d.same_inputs) d.notes.push_back("the inputs differ (names or sha256)");
  const size_t n = std::min(a.steps.size(), b.steps.size());
  for (size_t k = 0; k < n; ++k) {
    const auto& x = a.steps[k];
    const auto& y = b.steps[k];
    bool differs = false;
    if (x.engine != y.engine) {
      d.notes.push_back("step " + std::to_string(k + 1) + ": " + x.engine + " vs " + y.engine);
      ++d.differing_steps;
      continue;
    }
    std::map<std::string, std::string> px(x.params.begin(), x.params.end()), py(y.params.begin(), y.params.end());
    std::set<std::string> keys;
    for (const auto& [kk, v] : x.params) keys.insert(kk);
    for (const auto& [kk, v] : y.params) keys.insert(kk);
    for (const auto& [kk, v] : x.params) {
      const auto va = v, vb = py.count(kk) ? py[kk] : "—";
      if (va != vb) { d.rows.push_back({int(k), x.engine, kk, va, vb}); differs = true; }
    }
    for (const auto& [kk, v] : y.params)
      if (!px.count(kk)) { d.rows.push_back({int(k), x.engine, kk, "—", v}); differs = true; }
    if (x.rng != y.rng) { d.rows.push_back({int(k), x.engine, "rng", x.rng, y.rng}); differs = true; }
    d.differing_steps += differs;
  }
  for (size_t k = n; k < a.steps.size(); ++k) d.notes.push_back("step " + std::to_string(k + 1) + " (" + a.steps[k].engine + ") only in the first");
  for (size_t k = n; k < b.steps.size(); ++k) d.notes.push_back("step " + std::to_string(k + 1) + " (" + b.steps[k].engine + ") only in the second");
  return d;
}

bool known_citation(const std::string& key) {
  return std::any_of(std::begin(kCites), std::end(kCites), [&](const Cite& c) { return key == c.key; });
}

std::string bibtex(const std::vector<std::string>& keys) {
  std::string out;
  for (const auto& k : keys)
    for (const auto& c : kCites)
      if (k == c.key) { out += c.entry; out += "\n\n"; }
  return out;
}

std::vector<std::string> all_cites(const Manifest& m) {
  std::vector<std::string> out;
  for (const auto& s : m.steps)
    for (const auto& c : s.cites)
      if (std::find(out.begin(), out.end(), c) == out.end()) out.push_back(c);
  return out;
}

}  // namespace caps
