#include "caps/ffio.hpp"

#include <functional>

#include "caps/json.hpp"

namespace caps {

namespace {

Json nums(std::initializer_list<double> v) {
  Json a = Json::array();
  for (double x : v) a.push_back(Json(x));
  return a;
}
template <class T, class F>
Json rows(const std::vector<T>& v, F f) {
  Json a = Json::array();
  for (const auto& x : v) a.push_back(f(x));
  return a;
}
Json strings(const std::vector<std::string>& v) {
  Json a = Json::array();
  for (const auto& s : v) a.push_back(Json(s));
  return a;
}
template <class T>
Json numbers(const std::vector<T>& v) {
  Json a = Json::array();
  for (const auto& x : v) a.push_back(Json(double(x)));
  return a;
}

// reading
const Json& at(const Json& j, const char* k) {
  static const Json empty = Json::array();
  return j.has(k) ? j[k] : empty;
}
double d(const Json& r, size_t i) { return r[i].number(); }
uint32_t u(const Json& r, size_t i) { return uint32_t(r[i].number()); }
int n(const Json& r, size_t i) { return int(r[i].number()); }
template <class T, class F>
void read_rows(const Json& j, const char* k, std::vector<T>& out, F f) {
  out.clear();
  for (const auto& r : at(j, k).items()) out.push_back(f(r));
}
std::vector<std::string> read_strings(const Json& j, const char* k) {
  std::vector<std::string> v;
  for (const auto& s : at(j, k).items()) v.push_back(s.str());
  return v;
}
template <class T>
std::vector<T> read_numbers(const Json& j, const char* k) {
  std::vector<T> v;
  for (const auto& x : at(j, k).items()) v.push_back(T(x.number()));
  return v;
}
bool flag(const Json& j, const char* k, bool def) { return j.has(k) ? j[k].boolean() : def; }

}  // namespace

std::string forcefield_to_json(const ForceField& ff) {
  Json j = Json::object();
  j["format"] = Json("caps-assigned-forcefield");
  j["version"] = Json(1);
  j["atoms"] = Json(double(ff.type_index.size()));
  j["name"] = Json(ff.name);
  j["atom_type"] = strings(ff.atom_type);
  j["why"] = strings(ff.why);
  j["type_index"] = numbers(ff.type_index);
  j["type_names"] = strings(ff.type_names);
  j["lj"] = rows(ff.lj, [](const PairType& p) { return nums({p.eps, p.sigma}); });
  j["charge"] = numbers(ff.charge);
  j["mass"] = numbers(ff.mass);
  j["bonds"] = rows(ff.bonds, [](const BondTerm& t) { return nums({double(t.i), double(t.j), t.k, t.r0}); });
  j["angles"] = rows(ff.angles, [](const AngleTerm& t) { return nums({double(t.i), double(t.j), double(t.k), t.kt, t.theta0}); });
  auto tors = [](const TorsionTerm& t) { return nums({double(t.i), double(t.j), double(t.k), double(t.l), t.v, double(t.n), t.delta}); };
  j["dihedrals"] = rows(ff.dihedrals, tors);
  j["impropers"] = rows(ff.impropers, tors);
  j["impropers_dlpoly"] = rows(ff.impropers_dlpoly, tors);
  j["impropers_harmonic"] = rows(ff.impropers_harmonic, [](const HarmonicTorsion& t) { return nums({double(t.i), double(t.j), double(t.k), double(t.l), t.k2, t.chi0}); });
  j["inversions"] = rows(ff.inversions, [](const InversionTerm& t) { return nums({double(t.c), double(t.a), double(t.b), double(t.d), t.kw, t.w0, double(t.form)}); });
  j["bonds_x"] = rows(ff.bonds_x, [](const BondX& t) { return nums({double(t.i), double(t.j), double(t.form), t.a, t.b, t.c, t.d}); });
  j["angles_x"] = rows(ff.angles_x, [](const AngleX& t) { return nums({double(t.i), double(t.j), double(t.k), double(t.form), t.a, t.b}); });
  j["urey_bradley"] = rows(ff.urey_bradley, [](const UreyBradley& t) { return nums({double(t.i), double(t.k), t.kub, t.r0}); });
  j["cbt"] = rows(ff.cbt, [](const CbtTorsion& t) { return nums({double(t.i), double(t.j), double(t.k), double(t.l), t.a[0], t.a[1], t.a[2], t.a[3], t.a[4]}); });
  j["lj_pairs"] = rows(ff.lj_pairs, [](const PairLJ& t) { return nums({double(t.i), double(t.j), t.eps, t.sigma}); });
  {
    Json a = Json::array();
    for (const auto& [k, f] : ff.pair_func) a.push_back(nums({double(k.first), double(k.second), double(f.form), f.a, f.b, f.c}));
    j["pair_func"] = a;
  }
  j["lj14_types"] = rows(ff.lj14_types, [](const PairType& p) { return nums({p.eps, p.sigma}); });
  j["bonds2"] = rows(ff.bonds2, [](const Class2Bond& t) { return nums({double(t.i), double(t.j), t.r0, t.k2, t.k3, t.k4}); });
  j["angles2"] = rows(ff.angles2, [](const Class2Angle& t) {
    return nums({double(t.i), double(t.j), double(t.k), t.theta0, t.k2, t.k3, t.k4, t.bb_m, t.bb_r1, t.bb_r2, t.ba_n1, t.ba_n2, t.ba_r1, t.ba_r2});
  });
  j["dihedrals2"] = rows(ff.dihedrals2, [](const Class2Dihedral& t) {
    return nums({double(t.i), double(t.j), double(t.k), double(t.l), t.k1, t.phi1, t.k2, t.phi2, t.k3, t.phi3, t.mbt[0], t.mbt[1], t.mbt[2], t.mbt_r2,
                 t.ebt_b[0], t.ebt_b[1], t.ebt_b[2], t.ebt_c[0], t.ebt_c[1], t.ebt_c[2], t.ebt_r1, t.ebt_r3, t.at_d[0], t.at_d[1], t.at_d[2], t.at_e[0],
                 t.at_e[1], t.at_e[2], t.at_theta1, t.at_theta2, t.aat_m, t.aat_theta1, t.aat_theta2, t.bb13_n, t.bb13_r1, t.bb13_r3});
  });
  j["impropers2"] = rows(ff.impropers2, [](const Class2Improper& t) {
    return nums({double(t.i), double(t.j), double(t.k), double(t.l), t.kchi, t.chi0, t.m1, t.m2, t.m3, t.theta1, t.theta2, t.theta3});
  });
  j["pair_form"] = Json(ff.pair_form);
  j["mixing"] = Json(ff.mixing);
  j["native_pair"] = Json(ff.native_pair), j["native_dihedral"] = Json(ff.native_dihedral), j["native_improper"] = Json(ff.native_improper);
  j["improper_written"] = Json(ff.improper_written);
  j["impropers_dlpoly_ok"] = Json(ff.impropers_dlpoly_ok), j["impropers_dlpoly_set"] = Json(ff.impropers_dlpoly_set);
  j["native_cutoff"] = Json(ff.native_cutoff);
  j["native_gromacs_lj"] = Json(ff.native_gromacs_lj);
  j["native_special"] = Json(ff.native_special);
  j["native_timestep"] = Json(ff.native_timestep);
  {
    Json a = Json::array();
    for (const auto& [k, p] : ff.pair_override) a.push_back(nums({double(k.first), double(k.second), p.eps, p.sigma}));
    j["pair_override"] = a;
  }
  j["pairs14"] = rows(ff.pairs14, [](const std::array<uint32_t, 2>& p) { return nums({double(p[0]), double(p[1])}); });
  j["pairs14_lj"] = numbers(ff.pairs14_lj);
  j["pairs14_coul"] = numbers(ff.pairs14_coul);
  j["type_part"] = numbers(ff.type_part);
  j["part14"] = rows(ff.part14, [](const std::array<double, 2>& p) { return nums({p[0], p[1]}); });
  j["excluded"] = rows(ff.excluded, [](const std::vector<uint32_t>& e) { return numbers(e); });
  {
    Json s = Json::object();
    s["on"] = Json(ff.sw.on);
    s["p"] = nums({ff.sw.eps, ff.sw.sigma, ff.sw.a, ff.sw.lambda, ff.sw.gamma, ff.sw.cos0, ff.sw.A, ff.sw.B, ff.sw.p, ff.sw.q});
    s["atom"] = numbers(ff.sw.atom);
    j["sw"] = s;
  }
  {
    const auto& m = ff.manybody;
    Json o = Json::object();
    o["style"] = Json(m.style), o["file"] = Json(m.file), o["units"] = Json(m.units), o["tagged"] = Json(m.tagged);
    o["element"] = strings(m.element), o["citation"] = Json(m.citation), o["args"] = Json(m.args), o["file2"] = Json(m.file2);
    o["entry"] = strings(m.entry), o["extract"] = strings(m.extract), o["metal_only"] = Json(m.metal_only);
    j["manybody"] = o;
  }
  j["lj14"] = Json(ff.lj14), j["coul14"] = Json(ff.coul14), j["keep13"] = Json(ff.keep13);
  {
    Json a = Json::array();
    for (const auto& [x, y] : ff.excluded_type_pairs) a.push_back(nums({double(x), double(y)}));
    j["excluded_type_pairs"] = a;
  }
  j["coul_gromacs"] = Json(ff.coul_gromacs);
  j["coul_inner"] = Json(ff.coul_inner), j["lj_inner"] = Json(ff.lj_inner), j["dielectric"] = Json(ff.dielectric), j["cutoff"] = Json(ff.cutoff);
  j["coul_rf"] = Json(ff.coul_rf), j["eps_rf"] = Json(ff.eps_rf), j["lj_shift"] = Json(ff.lj_shift), j["lj_fsw"] = Json(ff.lj_fsw);
  {
    const auto& h = ff.hbond;
    Json o = Json::object();
    o["inner"] = Json(h.inner), o["outer"] = Json(h.outer), o["cos_cut"] = Json(h.cos_cut), o["power"] = Json(h.power), o["angle_deg"] = Json(h.angle_deg);
    o["hyd"] = rows(h.hyd, [](const std::vector<uint32_t>& e) { return numbers(e); });
    o["acceptor"] = numbers(h.acceptor);
    Json p = Json::array();
    for (const auto& [k, v] : h.param) p.push_back(nums({double(k.first), double(k.second), v[0], v[1], v[2]}));
    o["param"] = p;
    Json t = Json::array();
    for (const auto& [k, v] : h.htype) t.push_back(nums({double(k.first), double(k.second), double(v)}));
    o["htype"] = t;
    j["hbond"] = o;
  }
  j["vsites"] = rows(ff.vsites, [](const VirtualSite& v) {
    Json o = Json::object();
    o["site"] = Json(double(v.site)), o["from"] = numbers(v.from), o["w"] = numbers(v.w);
    if (v.c != 0) o["c"] = Json(v.c);
    return o;
  });
  j["notes"] = strings(ff.notes);
  return j.dump_exact(0);
}

ForceField forcefield_from_json(const std::string& text) {
  Json j;
  try { j = Json::parse(text); } catch (const std::exception& e) { throw FieldError(std::string("not a force field file: ") + e.what()); }
  if (!j.is_object() || j.text("format", "") != "caps-assigned-forcefield") throw FieldError("not a CAPS assigned force field (format caps-assigned-forcefield)");
  if (j.num("version", 0) > 1) throw FieldError("a newer CAPS wrote this force field file (version " + std::to_string(int(j.num("version", 0))) + ")");
  ForceField f;
  f.name = j.text("name", "");
  f.atom_type = read_strings(j, "atom_type");
  f.why = read_strings(j, "why");
  f.type_index = read_numbers<int>(j, "type_index");
  f.type_names = read_strings(j, "type_names");
  read_rows(j, "lj", f.lj, [](const Json& r) { return PairType{d(r, 0), d(r, 1)}; });
  f.charge = read_numbers<double>(j, "charge");
  f.mass = read_numbers<double>(j, "mass");
  read_rows(j, "bonds", f.bonds, [](const Json& r) { return BondTerm{u(r, 0), u(r, 1), d(r, 2), d(r, 3)}; });
  read_rows(j, "angles", f.angles, [](const Json& r) { return AngleTerm{u(r, 0), u(r, 1), u(r, 2), d(r, 3), d(r, 4)}; });
  auto tors = [](const Json& r) { return TorsionTerm{u(r, 0), u(r, 1), u(r, 2), u(r, 3), d(r, 4), n(r, 5), d(r, 6)}; };
  read_rows(j, "dihedrals", f.dihedrals, tors);
  read_rows(j, "impropers", f.impropers, tors);
  read_rows(j, "impropers_dlpoly", f.impropers_dlpoly, tors);
  read_rows(j, "impropers_harmonic", f.impropers_harmonic, [](const Json& r) { return HarmonicTorsion{u(r, 0), u(r, 1), u(r, 2), u(r, 3), d(r, 4), d(r, 5)}; });
  read_rows(j, "inversions", f.inversions, [](const Json& r) { return InversionTerm{u(r, 0), u(r, 1), u(r, 2), u(r, 3), d(r, 4), d(r, 5), n(r, 6)}; });
  read_rows(j, "bonds_x", f.bonds_x, [](const Json& r) { return BondX{u(r, 0), u(r, 1), n(r, 2), d(r, 3), d(r, 4), d(r, 5), d(r, 6)}; });
  read_rows(j, "angles_x", f.angles_x, [](const Json& r) { return AngleX{u(r, 0), u(r, 1), u(r, 2), n(r, 3), d(r, 4), d(r, 5)}; });
  read_rows(j, "urey_bradley", f.urey_bradley, [](const Json& r) { return UreyBradley{u(r, 0), u(r, 1), d(r, 2), d(r, 3)}; });
  read_rows(j, "cbt", f.cbt, [](const Json& r) {
    CbtTorsion t{u(r, 0), u(r, 1), u(r, 2), u(r, 3), {}};
    for (int q = 0; q < 5; ++q) t.a[q] = d(r, size_t(4 + q));
    return t;
  });
  read_rows(j, "lj_pairs", f.lj_pairs, [](const Json& r) { return PairLJ{u(r, 0), u(r, 1), d(r, 2), d(r, 3)}; });
  for (const auto& r : at(j, "pair_func").items()) f.pair_func[{n(r, 0), n(r, 1)}] = PairFunc{n(r, 2), d(r, 3), d(r, 4), d(r, 5)};
  read_rows(j, "lj14_types", f.lj14_types, [](const Json& r) { return PairType{d(r, 0), d(r, 1)}; });
  read_rows(j, "bonds2", f.bonds2, [](const Json& r) { return Class2Bond{u(r, 0), u(r, 1), d(r, 2), d(r, 3), d(r, 4), d(r, 5)}; });
  read_rows(j, "angles2", f.angles2, [](const Json& r) {
    return Class2Angle{u(r, 0), u(r, 1), u(r, 2), d(r, 3), d(r, 4), d(r, 5), d(r, 6), d(r, 7), d(r, 8), d(r, 9), d(r, 10), d(r, 11), d(r, 12), d(r, 13)};
  });
  read_rows(j, "dihedrals2", f.dihedrals2, [](const Json& r) {
    Class2Dihedral t{};
    t.i = u(r, 0), t.j = u(r, 1), t.k = u(r, 2), t.l = u(r, 3);
    t.k1 = d(r, 4), t.phi1 = d(r, 5), t.k2 = d(r, 6), t.phi2 = d(r, 7), t.k3 = d(r, 8), t.phi3 = d(r, 9);
    for (int q = 0; q < 3; ++q) t.mbt[q] = d(r, size_t(10 + q));
    t.mbt_r2 = d(r, 13);
    for (int q = 0; q < 3; ++q) t.ebt_b[q] = d(r, size_t(14 + q)), t.ebt_c[q] = d(r, size_t(17 + q));
    t.ebt_r1 = d(r, 20), t.ebt_r3 = d(r, 21);
    for (int q = 0; q < 3; ++q) t.at_d[q] = d(r, size_t(22 + q)), t.at_e[q] = d(r, size_t(25 + q));
    t.at_theta1 = d(r, 28), t.at_theta2 = d(r, 29);
    t.aat_m = d(r, 30), t.aat_theta1 = d(r, 31), t.aat_theta2 = d(r, 32);
    t.bb13_n = d(r, 33), t.bb13_r1 = d(r, 34), t.bb13_r3 = d(r, 35);
    return t;
  });
  read_rows(j, "impropers2", f.impropers2, [](const Json& r) {
    return Class2Improper{u(r, 0), u(r, 1), u(r, 2), u(r, 3), d(r, 4), d(r, 5), d(r, 6), d(r, 7), d(r, 8), d(r, 9), d(r, 10), d(r, 11)};
  });
  f.pair_form = j.text("pair_form", f.pair_form);
  f.mixing = j.text("mixing", f.mixing);
  f.native_pair = j.text("native_pair", ""), f.native_dihedral = j.text("native_dihedral", ""), f.native_improper = j.text("native_improper", "");
  f.improper_written = j.text("improper_written", "");
  f.impropers_dlpoly_ok = flag(j, "impropers_dlpoly_ok", true), f.impropers_dlpoly_set = flag(j, "impropers_dlpoly_set", false);
  f.native_cutoff = j.num("native_cutoff", 0);
  f.native_gromacs_lj = j.text("native_gromacs_lj", "");
  f.native_special = j.text("native_special", "");
  f.native_timestep = j.num("native_timestep", 0);
  for (const auto& r : at(j, "pair_override").items()) f.pair_override[{n(r, 0), n(r, 1)}] = PairType{d(r, 2), d(r, 3)};
  read_rows(j, "pairs14", f.pairs14, [](const Json& r) { return std::array<uint32_t, 2>{u(r, 0), u(r, 1)}; });
  f.pairs14_lj = read_numbers<double>(j, "pairs14_lj");
  f.pairs14_coul = read_numbers<double>(j, "pairs14_coul");
  f.type_part = read_numbers<int>(j, "type_part");
  read_rows(j, "part14", f.part14, [](const Json& r) { return std::array<double, 2>{d(r, 0), d(r, 1)}; });
  read_rows(j, "excluded", f.excluded, [](const Json& r) {
    std::vector<uint32_t> e;
    for (const auto& x : r.items()) e.push_back(uint32_t(x.number()));
    return e;
  });
  if (j.has("sw")) {
    const Json& s = j["sw"];
    f.sw.on = flag(s, "on", false);
    const Json& p = s["p"];
    f.sw.eps = d(p, 0), f.sw.sigma = d(p, 1), f.sw.a = d(p, 2), f.sw.lambda = d(p, 3), f.sw.gamma = d(p, 4), f.sw.cos0 = d(p, 5), f.sw.A = d(p, 6), f.sw.B = d(p, 7),
    f.sw.p = d(p, 8), f.sw.q = d(p, 9);
    f.sw.atom = read_numbers<char>(s, "atom");
  }
  if (j.has("manybody")) {
    const Json& m = j["manybody"];
    f.manybody.style = m.text("style", ""), f.manybody.file = m.text("file", ""), f.manybody.units = m.text("units", "");
    f.manybody.tagged = flag(m, "tagged", false);
    f.manybody.element = read_strings(m, "element"), f.manybody.citation = m.text("citation", ""), f.manybody.args = m.text("args", "");
    f.manybody.file2 = m.text("file2", ""), f.manybody.entry = read_strings(m, "entry"), f.manybody.extract = read_strings(m, "extract");
    f.manybody.metal_only = flag(m, "metal_only", false);
  }
  f.lj14 = j.num("lj14", f.lj14), f.coul14 = j.num("coul14", f.coul14), f.keep13 = flag(j, "keep13", false);
  for (const auto& r : at(j, "excluded_type_pairs").items()) f.excluded_type_pairs.insert({n(r, 0), n(r, 1)});
  f.coul_gromacs = flag(j, "coul_gromacs", false);
  f.coul_inner = j.num("coul_inner", 0), f.lj_inner = j.num("lj_inner", 0), f.dielectric = j.num("dielectric", 1), f.cutoff = j.num("cutoff", 0);
  f.coul_rf = flag(j, "coul_rf", false), f.eps_rf = j.num("eps_rf", 0), f.lj_shift = flag(j, "lj_shift", false), f.lj_fsw = flag(j, "lj_fsw", false);
  if (j.has("hbond")) {
    const Json& h = j["hbond"];
    f.hbond.inner = h.num("inner", f.hbond.inner), f.hbond.outer = h.num("outer", f.hbond.outer), f.hbond.cos_cut = h.num("cos_cut", 0);
    f.hbond.power = int(h.num("power", 4)), f.hbond.angle_deg = h.num("angle_deg", 90);
    read_rows(h, "hyd", f.hbond.hyd, [](const Json& r) {
      std::vector<uint32_t> e;
      for (const auto& x : r.items()) e.push_back(uint32_t(x.number()));
      return e;
    });
    f.hbond.acceptor = read_numbers<char>(h, "acceptor");
    for (const auto& r : at(h, "param").items()) f.hbond.param[{n(r, 0), n(r, 1)}] = {d(r, 2), d(r, 3), d(r, 4)};
    for (const auto& r : at(h, "htype").items()) f.hbond.htype[{n(r, 0), n(r, 1)}] = n(r, 2);
  }
  read_rows(j, "vsites", f.vsites, [](const Json& r) {
    VirtualSite v;
    v.site = uint32_t(r["site"].number());
    for (const auto& x : r["from"].items()) v.from.push_back(uint32_t(x.number()));
    for (const auto& x : r["w"].items()) v.w.push_back(x.number());
    v.c = r.num("c", 0.0);
    return v;
  });
  f.notes = read_strings(j, "notes");
  const size_t na = size_t(j.num("atoms", double(f.type_index.size())));
  if (f.type_index.size() != na || f.charge.size() != na || f.mass.size() != na) throw FieldError("force field file: per-atom lists of different lengths");
  return f;
}

}  // namespace caps
