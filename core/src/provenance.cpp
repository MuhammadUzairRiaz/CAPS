#include "caps/provenance.hpp"

#include <algorithm>
#include <ctime>
#include <fstream>
#include <cstring>
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
    {"soldera2006", "@article{soldera2006,\n  author = {Soldera, A. and Metatla, N.},\n  title = {Glass transition of polymers: atomistic simulation versus experiments},\n  journal = {Phys. Rev. E}, volume = {74}, pages = {061803}, year = {2006}, doi = {10.1103/PhysRevE.74.061803}\n}"},
    {"yeh2004", "@article{yeh2004,\n  author = {Yeh, I.-C. and Hummer, G.},\n  title = {System-size dependence of diffusion coefficients and viscosities from molecular dynamics simulations with periodic boundary conditions},\n  journal = {J. Phys. Chem. B}, volume = {108}, pages = {15873--15879}, year = {2004}, doi = {10.1021/jp0477147}\n}"},
    {"bernetti2020", "@article{bernetti2020,\n  author = {Bernetti, M. and Bussi, G.},\n  title = {Pressure control using stochastic cell rescaling},\n  journal = {J. Chem. Phys.}, volume = {153}, pages = {114107}, year = {2020}, doi = {10.1063/5.0020514}\n}"},
    {"berendsen1984", "@article{berendsen1984,\n  author = {Berendsen, H. J. C. and Postma, J. P. M. and van Gunsteren, W. F. and DiNola, A. and Haak, J. R.},\n  title = {Molecular dynamics with coupling to an external bath},\n  journal = {J. Chem. Phys.}, volume = {81}, pages = {3684--3690}, year = {1984}, doi = {10.1063/1.448118}\n}"},
    {"ryckaert1977", "@article{ryckaert1977,\n  author = {Ryckaert, J.-P. and Ciccotti, G. and Berendsen, H. J. C.},\n  title = {Numerical integration of the cartesian equations of motion of a system with constraints: molecular dynamics of n-alkanes},\n  journal = {J. Comput. Phys.}, volume = {23}, pages = {327--341}, year = {1977}, doi = {10.1016/0021-9991(77)90098-5}\n}"},
    {"andersen1983", "@article{andersen1983,\n  author = {Andersen, H. C.},\n  title = {Rattle: a ``velocity'' version of the shake algorithm for molecular dynamics calculations},\n  journal = {J. Comput. Phys.}, volume = {52}, pages = {24--34}, year = {1983}, doi = {10.1016/0021-9991(83)90014-1}\n}"},
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
    {"rosenbluth1955", "@article{rosenbluth1955,\n  author = {Rosenbluth, M. N. and Rosenbluth, A. W.},\n  title = {Monte {C}arlo calculation of the average extension of molecular chains},\n  journal = {J. Chem. Phys.}, volume = {23}, pages = {356--359}, year = {1955}, doi = {10.1063/1.1741967}\n}"},
    {"theodorou1985", "@article{theodorou1985,\n  author = {Theodorou, D. N. and Suter, U. W.},\n  title = {Detailed molecular structure of a vinyl polymer glass},\n  journal = {Macromolecules}, volume = {18}, pages = {1467--1478}, year = {1985}, doi = {10.1021/ma00149a018}\n}"},
    {"siepmann1992", "@article{siepmann1992,\n  author = {Siepmann, J. I. and Frenkel, D.},\n  title = {Configurational bias {M}onte {C}arlo: a new sampling scheme for flexible chains},\n  journal = {Mol. Phys.}, volume = {75}, pages = {59--70}, year = {1992}, doi = {10.1080/00268979200100061}\n}"},
    {"jorgensen1984", "@article{jorgensen1984,\n  author = {Jorgensen, W. L. and Madura, J. D. and Swenson, C. J.},\n  title = {Optimized intermolecular potential functions for liquid hydrocarbons},\n  journal = {J. Am. Chem. Soc.}, volume = {106}, pages = {6638--6646}, year = {1984}, doi = {10.1021/ja00334a030}\n}"},
    {"tuckerman1992", "@article{tuckerman1992,\n  author = {Tuckerman, M. and Berne, B. J. and Martyna, G. J.},\n  title = {Reversible multiple time scale molecular dynamics},\n  journal = {J. Chem. Phys.}, volume = {97}, pages = {1990--2001}, year = {1992}, doi = {10.1063/1.463137}\n}"},
    {"leimkuhler2013", "@article{leimkuhler2013,\n  author = {Leimkuhler, B. and Matthews, C.},\n  title = {Rational construction of stochastic numerical methods for molecular sampling},\n  journal = {Appl. Math. Res. Express}, volume = {2013}, pages = {34--56}, year = {2013}, doi = {10.1093/amrx/abs010}\n}"},
    {"ewald1921", "@article{ewald1921,\n  author = {Ewald, P. P.},\n  title = {Die Berechnung optischer und elektrostatischer Gitterpotentiale},\n  journal = {Ann. Phys.}, volume = {369}, pages = {253--287}, year = {1921}, doi = {10.1002/andp.19213690304}\n}"},
    {"faber1965", "@article{faber1965,\n  author = {Faber, T. E. and Ziman, J. M.},\n  title = {A theory of the electrical properties of liquid metals {III}. {T}he resistivity of binary alloys},\n  journal = {Phil. Mag.}, volume = {11}, pages = {153--173}, year = {1965}, doi = {10.1080/14786436508211931}\n}"},
    {"lorch1969", "@article{lorch1969,\n  author = {Lorch, E.},\n  title = {Neutron diffraction by germania, silica and radiation-damaged silica glasses},\n  journal = {J. Phys. C}, volume = {2}, pages = {229--237}, year = {1969}, doi = {10.1088/0022-3719/2/2/305}\n}"},
    {"sears1992", "@article{sears1992,\n  author = {Sears, V. F.},\n  title = {Neutron scattering lengths and cross sections},\n  journal = {Neutron News}, volume = {3}, number = {3}, pages = {26--37}, year = {1992}, doi = {10.1080/10448639208218770}\n}"},
    {"bondi1964", "@article{bondi1964,\n  author = {Bondi, A.},\n  title = {van der {W}aals volumes and radii},\n  journal = {J. Phys. Chem.}, volume = {68}, pages = {441--451}, year = {1964}, doi = {10.1021/j100785a001}\n}"},
    {"gelb1999", "@article{gelb1999,\n  author = {Gelb, L. D. and Gubbins, K. E.},\n  title = {Pore size distributions in porous glasses: a computer simulation study},\n  journal = {Langmuir}, volume = {15}, pages = {305--308}, year = {1999}, doi = {10.1021/la9808418}\n}"},
    {"einstein1905", "@article{einstein1905,\n  author = {Einstein, A.},\n  title = {{\\\"U}ber die von der molekularkinetischen {T}heorie der {W}{\\\"a}rme geforderte {B}ewegung von in ruhenden {F}l{\\\"u}ssigkeiten suspendierten {T}eilchen},\n  journal = {Ann. Phys.}, volume = {322}, pages = {549--560}, year = {1905}, doi = {10.1002/andp.19053220806}\n}"},
    {"theodorou1986", "@article{theodorou1986,\n  author = {Theodorou, D. N. and Suter, U. W.},\n  title = {Atomistic modeling of mechanical properties of polymeric glasses},\n  journal = {Macromolecules}, volume = {19}, pages = {139--154}, year = {1986}, doi = {10.1021/ma00155a022}\n}"},
    {"lutsko1989", "@article{lutsko1989,\n  author = {Lutsko, J. F.},\n  title = {Generalized expressions for the calculation of elastic constants by computer simulation},\n  journal = {J. Appl. Phys.}, volume = {65}, pages = {2991--2997}, year = {1989}, doi = {10.1063/1.342716}\n}"},
    {"clavier2017", "@article{clavier2017,\n  author = {Clavier, G. and Desbiens, N. and Bourasseau, E. and Lachet, V. and Brusselle-Dupend, N. and Rousseau, B.},\n  title = {Computation of elastic constants of solids using molecular simulation: comparison of constant volume and constant pressure ensemble methods},\n  journal = {Mol. Simul.}, volume = {43}, pages = {1413--1422}, year = {2017}, doi = {10.1080/08927022.2017.1313418}\n}"},
    {"prince2004", "@book{prince2004,\n  editor = {Prince, E.},\n  title = {International Tables for Crystallography, Volume C: Mathematical, Physical and Chemical Tables},\n  edition = {3rd}, publisher = {Kluwer}, year = {2004}, note = {Table 6.1.1.4, Cromer--Mann coefficients}\n}"},
    {"fan1992", "@article{fan1992,\n  author = {Fan, C. F. and Olafson, B. D. and Blanco, M. and Hsu, S. L.},\n  title = {Application of molecular simulation to derive phase diagrams of binary mixtures},\n  journal = {Macromolecules}, volume = {25}, pages = {3667--3676}, year = {1992}\n}"},
    {"blanco1991", "@article{blanco1991,\n  author = {Blanco, M.},\n  title = {Molecular silverware. I. General solutions to excluded volume constrained problems},\n  journal = {J. Comput. Chem.}, volume = {12}, pages = {237--247}, year = {1991}\n}"},
    {"kremer1990", "@article{kremer1990,\n  author = {Kremer, K. and Grest, G. S.},\n  title = {Dynamics of entangled linear polymer melts: a molecular-dynamics simulation},\n  journal = {J. Chem. Phys.}, volume = {92}, pages = {5057--5086}, year = {1990}, doi = {10.1063/1.458541}\n}"},
    {"parsons2005", "@article{parsons2005,\n  author = {Parsons, J. and Holmes, J. B. and Rojas, J. M. and Tsai, J. and Strauss, C. E. M.},\n  title = {Practical conversion from torsion space to {C}artesian space for in silico protein synthesis},\n  journal = {J. Comput. Chem.}, volume = {26}, pages = {1063--1068}, year = {2005}, doi = {10.1002/jcc.20237}\n}"},
    {"larsen2016", "@article{larsen2016,\n  author = {Larsen, P. M. and Schmidt, S. and Schi{\\o}tz, J.},\n  title = {Robust structural identification via polyhedral template matching},\n  journal = {Modelling Simul. Mater. Sci. Eng.}, volume = {24}, pages = {055007}, year = {2016}, doi = {10.1088/0965-0393/24/5/055007}\n}"},
    {"rycroft2009", "@article{rycroft2009,\n  author = {Rycroft, C. H.},\n  title = {{VORO++}: A three-dimensional {V}oronoi cell library in {C++}},\n  journal = {Chaos}, volume = {19}, pages = {041111}, year = {2009}, doi = {10.1063/1.3215722}\n}"},
    {"gellatly1982", "@article{gellatly1982,\n  author = {Gellatly, B. J. and Finney, J. L.},\n  title = {Characterisation of models of multicomponent amorphous metals: the radical alternative to the {V}oronoi polyhedron},\n  journal = {J. Non-Cryst. Solids}, volume = {50}, pages = {313--329}, year = {1982}, doi = {10.1016/0022-3093(82)90093-X}\n}"},
    {"love2005", "@article{love2005,\n  author = {Love, J. C. and Estroff, L. A. and Kriebel, J. K. and Nuzzo, R. G. and Whitesides, G. M.},\n  title = {Self-assembled monolayers of thiolates on metals as a form of nanotechnology},\n  journal = {Chem. Rev.}, volume = {105}, pages = {1103--1170}, year = {2005}, doi = {10.1021/cr0300789}\n}"},
    {"banks2005", "@article{banks2005,\n  author = {Banks, J. L. and Beard, H. S. and Cao, Y. and Cho, A. E. and Damm, W. and Farid, R. and Felts, A. K. and Halgren, T. A. and Mainz, D. T. and Maple, J. R. and others},\n  title = {Integrated Modeling Program, Applied Chemical Theory ({IMPACT})},\n  journal = {J. Comput. Chem.}, volume = {26}, pages = {1752--1780}, year = {2005}, doi = {10.1002/jcc.20292}\n}"},
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

std::string citation_text(const std::string& key) {
  for (const auto& c : kCites) {
    if (key != c.key) continue;
    const std::string e = c.entry;
    auto field = [&](const std::string& name) -> std::string {
      const size_t at = e.find(name + " = {");
      if (at == std::string::npos) return "";
      size_t i = at + name.size() + 4;
      int depth = 1;
      std::string out;
      for (; i < e.size() && depth > 0; ++i) {
        if (e[i] == '{') { ++depth; continue; }
        if (e[i] == '}') { if (--depth == 0) break; continue; }
        if (e[i] == '\\') {   // LaTeX accents: \'e → é and the like, kept simple
          // letters of their own: \o ø, \O Ø, \ss ß, \aa å (followed by a brace, a space or the end)
          static const std::pair<const char*, const char*> letters[] = {{"ss", "ß"}, {"aa", "å"}, {"AA", "Å"}, {"o", "ø"}, {"O", "Ø"}, {"l", "ł"}};
          bool done = false;
          for (const auto& [cmd, u] : letters) {
            const size_t k = std::strlen(cmd);
            if (e.compare(i + 1, k, cmd) == 0 && (i + 1 + k >= e.size() || e[i + 1 + k] == '}' || e[i + 1 + k] == ' ' || e[i + 1 + k] == '{')) {
              out += u;
              i += k + (i + 1 + k < e.size() && e[i + 1 + k] == ' ' ? 1 : 0);
              done = true;
              break;
            }
          }
          if (done) continue;
          if (i + 2 < e.size()) {
            const char a = e[i + 1], b = e[i + 2] == '{' && i + 3 < e.size() ? e[i + 3] : e[i + 2];
            static const std::map<std::string, std::string> acc = {{"'e", "é"}, {"'a", "á"}, {"'o", "ó"}, {"'i", "í"}, {"`e", "è"}, {"\"a", "ä"}, {"\"o", "ö"}, {"\"u", "ü"}, {"\"U", "Ü"}};
            if (auto it = acc.find(std::string{a, b}); it != acc.end()) {
              out += it->second;
              i += e[i + 2] == '{' ? 4 : 2;   // a closing brace after it is counted by the loop
              continue;
            }
          }
          continue;
        }
        out += e[i];
      }
      return out;
    };
    std::string authors = field("author");
    if (authors.empty()) authors = field("editor") + " (ed.)";
    // "A. and B. and C." → "A., B., C."
    for (size_t p = authors.find(" and "); p != std::string::npos; p = authors.find(" and ")) authors.replace(p, 5, ", ");
    std::string t = authors + ", \"" + field("title") + "\"";
    const std::string journal = field("journal");
    if (!journal.empty()) t += ", " + journal + " " + field("volume") + ", " + field("pages");
    else t += ", " + field("publisher");
    t += " (" + field("year") + ")";
    for (size_t p = t.find("--"); p != std::string::npos; p = t.find("--")) t.replace(p, 2, "–");
    const std::string doi = field("doi");
    if (!doi.empty()) t += ". doi:" + doi;
    return t;
  }
  return "";
}

namespace {
std::string param(const ProvStep& s, const std::string& key) {
  for (const auto& [k, v] : s.params) if (k == key) return v;
  return "";
}
}  // namespace

std::string methods_text(const Manifest& m, std::vector<std::string>* refs, const std::vector<Manifest>& replicas) {
  std::vector<std::string> order;
  auto cite = [&](const std::vector<std::string>& keys) {
    std::string out;
    for (const auto& k : keys) {
      if (!known_citation(k) || k == "matsumoto1998") continue;   // the generator is named, not cited
      auto it = std::find(order.begin(), order.end(), k);
      if (it == order.end()) { order.push_back(k); it = order.end() - 1; }
      out += (out.empty() ? "" : ", ") + std::to_string(it - order.begin() + 1);
    }
    return out.empty() ? std::string() : " [" + out + "]";
  };
  auto has = [](const ProvStep& s, const char* key) { return std::find(s.cites.begin(), s.cites.end(), std::string(key)) != s.cites.end(); };
  auto only = [&](const ProvStep& s, std::initializer_list<const char*> keys) {
    std::vector<std::string> v;
    for (const char* k : keys) if (has(s, k)) v.push_back(k);
    return v;
  };
  std::vector<std::string> sentences;
  std::string version = m.generator.empty() ? "CAPS" : m.generator;
  bool built = false;
  for (const auto& s : m.steps) {
    const std::string& e = s.engine;
    std::string t;
    if (e == "io.read" || e == "io.import") {
      if (!built) t = "The starting structure (" + param(s, "atoms") + " atoms) was read from " + param(s, "file") +
                      (e == "io.import" && param(s, "bonds").rfind("perceived", 0) == 0 ? ", with bonds perceived from covalent radii" + cite({"cordero2008"}) : std::string()) + ".";
      built = true;
    } else if (e == "chem.build") {
      t = "The molecule " + param(s, "smiles") + " was built from its SMILES string with " + version + (param(s, "clean-up") == "UFF" ? " and cleaned up with UFF" + cite({"rappe1992"}) : std::string()) + ".";
      built = true;
    } else if (e == "grow.trials" || e == "grow.blend" || e == "interface.build" || e == "nano.embed") {
      const std::string chains = param(s, "chains"), dp = param(s, "DP").empty() ? param(s, "dp") : param(s, "DP");
      const std::string where = !param(s, "density").empty() ? " at an initial density of " + param(s, "density") : !param(s, "box").empty() ? " in a " + param(s, "box") + " cubic cell" : std::string();
      t = (chains.empty() ? std::string("Polymer chains were") : chains + (chains == "1" ? " chain" : " chains") + (dp.empty() ? "" : " of " + dp + " repeat units") + (chains == "1" ? " was" : " were")) + " grown" +
          (e == "interface.build" ? " against the surface" : e == "nano.embed" ? " around the filler" : e == "grow.blend" ? " as a blend" : " in a periodic cell") + where +
          " with " + version + ", each unit placed from internal coordinates" + cite({"parsons2005"}) + " by choosing among trial torsions the one with the largest clearance from atoms already placed.";
      built = true;
    } else if (e == "field.assign") {
      const std::string ff = param(s, "force field"), ch = param(s, "charges");
      t = "Atom types were assigned from " + ff + cite(only(s, {"rappe1992", "wang2004"})) +
          (ch == "from the force field" ? ", with the force field's charges" : ", with " + ch + " charges" + cite(only(s, {"rappe1991", "gasteiger1980"}))) + ".";
    } else if (e.rfind("relax.", 0) == 0) {
      t = "The structure was minimised with " + param(s, "minimiser") + cite(only(s, {"liu1989", "polak1969", "bitzek2006"})) + " to a largest atomic force of " + param(s, "|F|max") +
          (param(s, "push-off") == "on" ? ", after capped-force push-off stages" + cite({"auhl2003"}) : std::string()) + ".";
    } else if (e.rfind("dynamics.", 0) == 0) {
      const bool npt = e == "dynamics.npt", nvt = e == "dynamics.nvt";
      std::string th = param(s, "thermostat");
      th = th.substr(0, th.find(" · "));
      const std::string len = param(s, "length");
      const size_t of = len.find(" of ");
      const std::string duration = len.substr(0, len.find(" · ")), step = of == std::string::npos ? std::string() : len.substr(of + 4);
      const std::string cons = param(s, "constraints");
      t = std::string(npt ? "NPT" : nvt ? "NVT" : "NVE") + " molecular dynamics (velocity Verlet" + cite({"swope1982"}) + (step.empty() ? "" : ", " + step + " time step") +
          (cons.empty() ? "" : ", " + std::string(cons.rfind("all bonds", 0) == 0 ? "all bonds" : "bonds to hydrogen") + " constrained with SHAKE/RATTLE" + cite({"ryckaert1977", "andersen1983"})) +
          ") was run for " + duration +
          " at " + param(s, "temperature") +
          (nvt || npt ? " with the " + th + " thermostat" + cite(only(s, {"bussi2007"})) : std::string()) +
          (npt ? " and " + param(s, "barostat").substr(0, param(s, "barostat").find(" · ")) + " pressure control" + cite(only(s, {"bernetti2020", "berendsen1984"})) : std::string()) + ".";
    } else if (e == "recipe.run") {
      t = "The structure was made by the CAPS recipe “" + param(s, "recipe") + "” (SHA-256 " + param(s, "sha256").substr(0, 12) + "…), whose steps follow.";
    } else if (e == "analysis.tg") {
      t = "The glass transition temperature was estimated from a stepwise NPT cooling scan from " + param(s, "from") + " to " + param(s, "to") + " in steps of " +
          param(s, "step") + ", holding " + param(s, "hold") + " at each temperature, by a two-line fit of the specific volume against temperature" +
          cite(only(s, {"soldera2006"})) + "; simulated cooling rates are many orders of magnitude faster than calorimetry, so Tg is expected above experiment.";
    } else if (e == "equilibrate.larsen21") {
      t = "The cell was equilibrated with the 21-step compression–decompression protocol (P_max = " + param(s, "Pmax") + ")" + cite({"larsen2011"}) + ", " + param(s, "length") + " of dynamics in all.";
    } else if (e == "equilibrate.protocol") {
      t = "The cell was equilibrated with a " + param(s, "stages") + "-stage protocol, " + param(s, "length") + " of dynamics in all.";
    } else if (e == "pack.lbfgs" || e == "pack.insert" || e == "solvate.pack") {
      t = std::string(e == "solvate.pack" ? "Solvent" : e == "pack.insert" ? "Additional molecules" : "Molecules") + " were packed without overlaps by minimising a pair-penalty function" + cite({"martinez2009"}) + ".";
      built = true;
    } else if (e == "nano.pore") {
      t = "The pore (" + s.summary + ") was built with " + version + cite(only(s, {"martinez2009"})) + "; its walls were held fixed.";
      built = true;
    } else if (e == "crystal.build") {
      t = "The crystal was built from its space group" + cite({"hall1981"}) + " and asymmetric unit.";
      built = true;
    } else if (e == "bio.peptide") {
      t = "The peptide " + param(s, "sequence") + " was built from its sequence with standard backbone geometry" + cite({"engh1991", "parsons2005"}) + ".";
      built = true;
    } else if (e == "edit.builder") {
      t = "The structure was edited by hand in the builder (" + param(s, "operations") + " operations).";
    } else if (e == "cg.kremer_grest") {
      t = param(s, "chains") + " Kremer–Grest bead-spring chains of " + param(s, "beads") + " beads" + cite({"kremer1990"}) + " were built as random walks at ρσ³ = " +
          param(s, "density").substr(0, param(s, "density").find(' ')) + (param(s, "k_theta").rfind("0 ", 0) == 0 ? "" : " with a bending stiffness k_θ = " + param(s, "k_theta")) + ".";
      built = true;
    } else if (e == "nano.build" || e == "surface.build") {
      t = std::string(e == "surface.build" ? "The surface slab" : "The nanostructure") + " was built with " + version + ".";
      built = true;
    }
    if (!t.empty()) sentences.push_back(t);
  }
  // approximations of the energy evaluations
  const auto ap = approximations(m);
  std::string vdw, elec;
  for (const auto& [k, v] : ap) { if (k == "van der Waals") vdw = v; if (k == "Electrostatics") elec = v; }
  if (!vdw.empty() || !elec.empty()) {
    // "cut-off 12 Å + tail correction" → "were truncated at 12 Å with analytic tail corrections"
    std::string t = "Van der Waals interactions";
    const size_t cp = vdw.find("cut-off ");
    if (cp != std::string::npos) {
      const std::string rest = vdw.substr(cp + 8);
      const size_t plus = rest.find(" + ");
      t += " were truncated at " + rest.substr(0, plus) + (plus != std::string::npos ? " with analytic tail corrections" : "");
    } else {
      t += ": " + vdw;
    }
    if (!elec.empty()) {
      const std::string sep = " · ";
      const size_t dot = elec.find(sep);
      const std::string detail = dot == std::string::npos ? std::string() : " (" + elec.substr(dot + sep.size()) + ")";
      if (elec.rfind("SPME", 0) == 0) t += "; electrostatics used smooth particle-mesh Ewald summation" + cite({"essmann1995"}) + detail;
      else if (elec.rfind("damped", 0) == 0) t += "; electrostatics used the damped shifted force method" + cite({"fennell2006"}) + detail;
      else if (elec == "off") t += "; electrostatics were switched off";
      else t += "; electrostatics: " + elec;
    }
    sentences.push_back(t + ".");
  }
  if (!replicas.empty())
    sentences.push_back(std::to_string(replicas.size() + 1) + " independent replicas were prepared the same way with different random seeds.");
  std::string out;
  for (const auto& s : sentences) out += (out.empty() ? "" : " ") + s;
  if (refs) {
    refs->clear();
    for (const auto& k : order) refs->push_back(citation_text(k));
  }
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
