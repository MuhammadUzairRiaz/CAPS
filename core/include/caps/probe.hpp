#pragma once
// Probes (design/boards/Probes): a point, plane, axis or ellipsoid made from a set of atoms (the molecule made whole first,
// each atom at its minimum image from the first), measured against another probe or atoms, frame by frame.
//   point      the mass-weighted centre
//   plane      through the centre, normal along the smallest principal direction of the mass-weighted covariance;
//              rms its atoms' distance from it (flatness)
//   axis       through the centre along the largest principal direction (a chain's long axis)
//   ellipsoid  the principal directions, semi-axes ∝ √λ scaled the least to enclose every atom
#include <array>
#include <string>
#include <vector>

#include "caps/system.hpp"

namespace caps {

enum class ProbeKind { Point, Plane, Axis, Ellipsoid };

struct Probe {
  ProbeKind kind = ProbeKind::Point;
  Vec3 centre{0, 0, 0};
  std::array<Vec3, 3> axes{};   // principal directions, largest variance first (unit)
  std::array<double, 3> semi{}; // ellipsoid semi-axes along them, Å
  double rms = 0;               // plane: the atoms' rms distance from it
  Vec3 direction() const { return kind == ProbeKind::Plane ? axes[2] : axes[0]; }   // the plane's normal, the axis
};

ProbeKind probe_kind(const std::string& name);   // "point" | "plane" | "axis" | "ellipsoid"; throws otherwise
Probe make_probe(const System& s, const std::vector<size_t>& atoms, ProbeKind kind);

// a measured against b: "distance" (Å; to a plane signed along its normal, else between centres, minimum image),
// "angle" (degrees, 0 – 90: between two lines, two planes' normals, a line and a plane), "rms" (a's flatness),
// "size" (a's largest semi-axis or, for a plane, its rms). b may be absent for rms and size.
double probe_measure(const System& s, const Probe& a, const Probe* b, const std::string& what);

}  // namespace caps
