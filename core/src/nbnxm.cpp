#include "nbnxm.hpp"

#include "parallel.hpp"
#include "nbnxm_fit.hpp"
#include "simd4.hpp"
#include <algorithm>
#include <atomic>
#include <array>
#include <cmath>
#include <cstdlib>
#include <string>
#include <type_traits>
#include <vector>

namespace caps::nbnxm {

void List::build(const std::vector<double>& x, const Cell& cell_in, double rl, const std::vector<std::vector<uint32_t>>& excluded,
                 const std::vector<int>& type_index, size_t ntypes, const std::vector<char>& skip_type, const std::vector<double>& q,
                 const std::vector<double>& f12, const std::vector<double>& f6, ThreadPool& pool, const float* c12, const float* c6,
                 const std::vector<std::vector<uint32_t>>* p14) {
  rlist = rl;
  const size_t n = x.size() / 3;
  Cell cell = cell_in;
  if (!cell.valid()) {   // no cell: a box around the atoms, not periodic
    Vec3 lo{1e300, 1e300, 1e300}, hi{-1e300, -1e300, -1e300};
    for (size_t i = 0; i < n; ++i)
      for (int k = 0; k < 3; ++k) { lo[k] = std::min(lo[k], x[3 * i + k]); hi[k] = std::max(hi[k], x[3 * i + k]); }
    if (n == 0) lo = hi = {0, 0, 0};
    cell.origin = lo - Vec3{1, 1, 1};
    cell.a = {hi[0] - lo[0] + 2, 0, 0};
    cell.b = {0, hi[1] - lo[1] + 2, 0};
    cell.c = {0, 0, hi[2] - lo[2] + 2};
    cell.periodic = {false, false, false};
  }
  const double v = cell.volume();
  const double w[3] = {v / norm(cross(cell.b, cell.c)), v / norm(cross(cell.c, cell.a)), v / norm(cross(cell.a, cell.b))};
  // columns along a and b about as wide as a cluster of four atoms at the mean density, atoms sorted along c in each
  const double side = std::cbrt(4.0 * v / std::max<double>(1.0, double(n)));
  const int nu = std::clamp(int(w[0] / side), 1, 1024), nv = std::clamp(int(w[1] / side), 1, 1024);
  std::vector<Vec3> fr(n);
  std::vector<uint32_t> col(n);
  std::vector<size_t> cstart(size_t(nu) * nv + 1, 0);
  for (size_t i = 0; i < n; ++i) {
    Vec3 f = cell.to_fractional({x[3 * i], x[3 * i + 1], x[3 * i + 2]});
    for (int k = 0; k < 3; ++k) {
      if (cell.periodic[k]) f[k] -= std::floor(f[k]);
      f[k] = std::clamp(f[k], 0.0, 1.0 - 1e-12);
    }
    fr[i] = f;
    const int cu = std::min(nu - 1, int(f[0] * nu)), cv = std::min(nv - 1, int(f[1] * nv));
    col[i] = uint32_t(cu * nv + cv);
    ++cstart[col[i] + 1];
  }
  const size_t ncols = size_t(nu) * nv;
  for (size_t c = 0; c < ncols; ++c) cstart[c + 1] += cstart[c];
  std::vector<uint32_t> order(n);
  {
    std::vector<size_t> fill(cstart.begin(), cstart.end() - 1);
    for (size_t i = 0; i < n; ++i) order[fill[col[i]]++] = uint32_t(i);
  }
  // clusters, column by column; each column's atoms by c, its last cluster padded
  std::vector<size_t> ccl(ncols + 1, 0);   // first cluster of each column
  for (size_t c = 0; c < ncols; ++c) {
    std::sort(order.begin() + std::ptrdiff_t(cstart[c]), order.begin() + std::ptrdiff_t(cstart[c + 1]),
              [&](uint32_t p, uint32_t r) { return fr[p][2] != fr[r][2] ? fr[p][2] < fr[r][2] : p < r; });
    ccl[c + 1] = ccl[c] + (cstart[c + 1] - cstart[c] + 3) / 4;
  }
  nclusters = ccl[ncols];
  const size_t ns = 4 * nclusters;
  perm.assign(ns, kPad);
  iperm.assign(n, 0);
  std::vector<uint32_t> ccol(nclusters);
  for (size_t c = 0; c < ncols; ++c)
    for (size_t k = cstart[c]; k < cstart[c + 1]; ++k) {
      const size_t s = 4 * ccl[c] + (k - cstart[c]);
      perm[s] = order[k];
      iperm[order[k]] = uint32_t(s);
    }
  for (size_t c = 0; c < ncols; ++c)
    for (size_t k = ccl[c]; k < ccl[c + 1]; ++k) ccol[k] = uint32_t(c);
  // wrapped Cartesian positions by slot (padding: those of the cluster's first atom), bounding boxes and c ranges
  std::vector<Vec3> cw(ns), frs(ns);
  for (size_t s = 0; s < ns; ++s) {
    const uint32_t a = perm[s] != kPad ? perm[s] : perm[s & ~size_t(3)];
    frs[s] = fr[a];
    cw[s] = cell.origin + cell.a * fr[a][0] + cell.b * fr[a][1] + cell.c * fr[a][2];
  }
  std::vector<Vec3> blo(nclusters), bhi(nclusters);
  std::vector<double> wlo(nclusters), whi(nclusters);
  for (size_t c = 0; c < nclusters; ++c) {
    blo[c] = bhi[c] = cw[4 * c];
    wlo[c] = whi[c] = frs[4 * c][2];
    for (size_t l = 1; l < 4; ++l)
      for (int k = 0; k < 3; ++k) {
        blo[c][k] = std::min(blo[c][k], cw[4 * c + l][k]);
        bhi[c][k] = std::max(bhi[c][k], cw[4 * c + l][k]);
      }
    for (size_t l = 1; l < 4; ++l) wlo[c] = std::min(wlo[c], frs[4 * c + l][2]), whi[c] = std::max(whi[c], frs[4 * c + l][2]);
  }
  // cluster origins (bounding-box centres), positions relative to them in single precision (the 16 atom pairs of a
  // candidate cluster pair are checked four at a time) and each cluster's real lanes
  org.resize(3 * nclusters);
  for (size_t c = 0; c < nclusters; ++c)
    for (int k = 0; k < 3; ++k) org[3 * c + k] = 0.5 * (blo[c][k] + bhi[c][k]);
  std::vector<float> bx(ns), by(ns), bz(ns);
  std::vector<uint8_t> real(nclusters, 0);
  for (size_t s = 0; s < ns; ++s) {
    const double* o = &org[3 * (s / 4)];
    bx[s] = float(cw[s][0] - o[0]), by[s] = float(cw[s][1] - o[1]), bz[s] = float(cw[s][2] - o[2]);
    if (perm[s] != kPad) real[s / 4] |= uint8_t(1u << (s % 4));
  }
  // excluded atom pairs by cluster pair: for each i cluster, its j clusters (sorted) with a mask of the excluded pairs,
  // so a candidate costs one look-up, not one per atom pair
  std::vector<size_t> xfirst(nclusters + 1, 0);
  std::vector<std::pair<uint32_t, uint16_t>> xm;
  std::vector<uint16_t> xm14;   // aligned with xm: the excluded pairs whose Coulomb keeps a fraction
  {
    struct X { uint32_t cj; uint16_t bit, is14; };
    std::vector<X> tmp;
    for (size_t c = 0; c < nclusters; ++c) {
      tmp.clear();
      for (int a = 0; a < 4; ++a) {
        const uint32_t ai = perm[4 * c + size_t(a)];
        if (ai == kPad) continue;
        for (uint32_t aj : excluded[ai]) {
          const uint32_t sj = iperm[aj];
          const uint16_t bit = uint16_t(1u << (4 * a + int(sj % 4)));
          const bool k14 = p14 && std::binary_search((*p14)[ai].begin(), (*p14)[ai].end(), aj);
          tmp.push_back({sj / 4, bit, uint16_t(k14 ? bit : 0)});
        }
      }
      std::sort(tmp.begin(), tmp.end(), [](const X& p, const X& r) { return p.cj != r.cj ? p.cj < r.cj : p.bit < r.bit; });
      for (const auto& e : tmp) {
        if (xm.size() > xfirst[c] && xm.back().first == e.cj) xm.back().second |= e.bit, xm14.back() |= e.is14;
        else xm.push_back({e.cj, e.bit}), xm14.push_back(e.is14);
      }
      xfirst[c + 1] = xm.size();
    }
  }
  // reach in columns and in images along c; image codes
  const int ru = cell.periodic[0] ? int(std::floor(rl * nu / w[0])) + 1 : nu;
  const int rv = cell.periodic[1] ? int(std::floor(rl * nv / w[1])) + 1 : nv;
  const int rq = cell.periodic[2] ? int(std::ceil(rl / w[2])) + 1 : 0;
  const double dw = rl / w[2];
  int qlo[3], span[3];
  qlo[0] = cell.periodic[0] ? int(std::floor(double(-ru) / nu)) : 0;
  qlo[1] = cell.periodic[1] ? int(std::floor(double(-rv) / nv)) : 0;
  qlo[2] = -rq;
  span[0] = (cell.periodic[0] ? int(std::floor(double(nu - 1 + ru) / nu)) : 0) - qlo[0] + 1;
  span[1] = (cell.periodic[1] ? int(std::floor(double(nv - 1 + rv) / nv)) : 0) - qlo[1] + 1;
  span[2] = 2 * rq + 1;
  if (size_t(span[0]) * span[1] * span[2] > 65536) throw std::runtime_error("the cell is too small for the cut-off: more than 65536 periodic images in range");
  qint.clear();
  for (int a = 0; a < span[0]; ++a)
    for (int b = 0; b < span[1]; ++b)
      for (int c = 0; c < span[2]; ++c) qint.insert(qint.end(), {qlo[0] + a, qlo[1] + b, qlo[2] + c});
  const double rl2 = rl * rl;
  const simd4::F4 rl2f = simd4::set1(float(rl2 * (1 + 1e-5)));
  const bool skipping = !skip_type.empty();
  // the cluster pairs of every i cluster, threads over contiguous ranges
  const int nth = pool.size();
  std::vector<std::vector<Entry>> ve(nth);
  std::vector<uint32_t> count(nclusters);
  pool.run(nclusters, [&](int t, size_t cb, size_t ce) {
    using namespace simd4;
    auto& E = ve[t];
    E.clear();
    for (size_t ci = cb; ci < ce; ++ci) {
      const size_t before = E.size();
      const int cu = int(ccol[ci]) / nv, cv = int(ccol[ci]) % nv;
      for (int du = -ru; du <= ru; ++du) {
        int U = cu + du, qa = 0;
        if (cell.periodic[0]) { qa = int(std::floor(double(U) / nu)); U -= qa * nu; }
        else if (U < 0 || U >= nu) continue;
        for (int dv = -rv; dv <= rv; ++dv) {
          int V = cv + dv, qb = 0;
          if (cell.periodic[1]) { qb = int(std::floor(double(V) / nv)); V -= qb * nv; }
          else if (V < 0 || V >= nv) continue;
          const size_t cc = size_t(U) * nv + V;
          // images along c whose range of c meets this cluster's, ± the list radius
          const int qc0 = cell.periodic[2] ? std::max(-rq, int(std::floor(wlo[ci] - dw))) : 0;
          const int qc1 = cell.periodic[2] ? std::min(rq, int(std::floor(whi[ci] + dw))) : 0;
          for (int qc = qc0; qc <= qc1; ++qc) {
            const Vec3 S = cell.a * double(qa) + cell.b * double(qb) + cell.c * double(qc);
            const bool zero = qa == 0 && qb == 0 && qc == 0;
            const bool positive = qa > 0 || (qa == 0 && (qb > 0 || (qb == 0 && qc > 0)));
            const uint16_t code = uint16_t(((qa - qlo[0]) * span[1] + (qb - qlo[1])) * span[2] + (qc - qlo[2]));
            // a column's clusters are in order of c (their low and high ends both rising): the ones in range by bisection
            const double lo = wlo[ci] - dw - qc, hi = whi[ci] + dw - qc;
            const size_t c0 = size_t(std::lower_bound(whi.begin() + std::ptrdiff_t(ccl[cc]), whi.begin() + std::ptrdiff_t(ccl[cc + 1]), lo) - whi.begin());
            for (size_t cj = c0; cj < ccl[cc + 1] && wlo[cj] <= hi; ++cj) {
              // half list: each cluster pair from one side (by the parity of ci + cj), a cluster with itself once per ± image
              if (cj != ci) {
                if ((cj > ci) != (((ci + cj) & 1) == 0)) continue;
              } else if (!zero && !positive) {
                continue;
              }
              double d2 = 0;
              for (int k = 0; k < 3; ++k) {
                const double g1 = blo[cj][k] + S[k] - bhi[ci][k], g2 = blo[ci][k] - (bhi[cj][k] + S[k]);
                const double g = g1 > 0 ? g1 : g2 > 0 ? g2 : 0.0;
                d2 += g * g;
              }
              if (d2 >= rl2) continue;
              // the atom pairs within the list radius (a hair wider in single precision, so rounding never drops one)
              const double* oi = &org[3 * ci];
              const double* oj = &org[3 * cj];
              const F4 xj = load(&bx[4 * cj]) + set1(float(S[0] + oj[0] - oi[0])), yj = load(&by[4 * cj]) + set1(float(S[1] + oj[1] - oi[1])),
                       zj = load(&bz[4 * cj]) + set1(float(S[2] + oj[2] - oi[2]));
              unsigned mask = 0;
              for (int a = 0; a < 4; ++a) {
                if (!((real[ci] >> a) & 1u)) continue;
                const size_t sa = 4 * ci + size_t(a);
                const F4 dx = xj - set1(bx[sa]), dy = yj - set1(by[sa]), dz = zj - set1(bz[sa]);
                unsigned row = bits(fma(dx, dx, fma(dy, dy, dz * dz)) < rl2f) & real[cj];
                if (cj == ci && zero) row &= (0xFu << (a + 1)) & 0xFu;   // a cluster with itself: each pair once
                if (!row) continue;
                if (skipping)
                  for (int l = 0; l < 4; ++l)
                    if (((row >> l) & 1u) && skip_type[size_t(type_index[perm[sa]]) * ntypes + size_t(type_index[perm[4 * cj + l]])]) row &= ~(1u << l);
                mask |= row << (4 * a);
              }
              // 1-2, 1-3 and 1-4 partners are left out at their nearest image only (as the atom pair list); with p14
              // they stay in the entry as excluded lanes, for the kernel's corrections
              unsigned exm = 0, ex14 = 0;
              if (mask) {
                const auto xb = xm.begin() + std::ptrdiff_t(xfirst[ci]), xe = xm.begin() + std::ptrdiff_t(xfirst[ci + 1]);
                const auto it = std::lower_bound(xb, xe, std::pair<uint32_t, uint16_t>{uint32_t(cj), 0});
                if (it != xe && it->first == cj && (it->second & mask)) {
                  const unsigned ex = it->second & mask, k14 = xm14[size_t(it - xm.begin())];
                  for (int b = 0; b < 16; ++b) {
                    if (!((ex >> b) & 1u)) continue;
                    const size_t sa = 4 * ci + size_t(b / 4), sb = 4 * cj + size_t(b % 4);
                    if (perm[sa] == perm[sb]) continue;   // an atom and its own image
                    const Vec3 d = cw[sb] + S - cw[sa];
                    Vec3 m = frs[sb] - frs[sa];
                    for (int k = 0; k < 3; ++k)
                      if (cell.periodic[k]) m[k] -= std::round(m[k]);
                    const Vec3 dm = cell.a * m[0] + cell.b * m[1] + cell.c * m[2];
                    if (dot(d, d) <= dot(dm, dm) + 1e-6) {
                      mask &= ~(1u << b);
                      if (p14) exm |= 1u << b, ex14 |= k14 & (1u << b);
                    }
                  }
                }
              }
              if (mask | exm) E.push_back({uint32_t(cj), code, uint16_t(mask), uint16_t(exm), uint16_t(ex14)});
            }
          }
        }
      }
      count[ci] = uint32_t(E.size() - before);
    }
  });
  first.assign(nclusters + 1, 0);
  for (size_t c = 0; c < nclusters; ++c) first[c + 1] = first[c] + count[c];
  entries.resize(first[nclusters]);
  pool.run(nclusters, [&](int t, size_t cb, size_t) {
    std::copy(ve[t].begin(), ve[t].end(), entries.begin() + std::ptrdiff_t(first[cb]));
  });
  // offsets to the binned positions, as lattice vectors kept by their integer coefficients
  off.assign(3 * ns, 0.0);
  ioff.assign(3 * ns, 0);
  for (size_t s = 0; s < ns; ++s) {
    if (perm[s] == kPad) continue;
    for (int k = 0; k < 3; ++k) off[3 * s + k] = cw[s][k] - x[3 * size_t(perm[s]) + k];
    if (cell_in.valid()) {
      const Vec3 o = cell_in.to_fractional(cell_in.origin + Vec3{off[3 * s], off[3 * s + 1], off[3 * s + 2]});
      for (int k = 0; k < 3; ++k) ioff[3 * s + k] = int32_t(std::lround(o[k]));
    }
  }
  qvec.resize(qint.size());
  for (size_t c = 0; c < qint.size(); c += 3) {
    const Vec3 sv = cell.a * qint[c] + cell.b * qint[c + 1] + cell.c * qint[c + 2];
    for (int k = 0; k < 3; ++k) qvec[c + k] = sv[k];
  }
  cells = cell_in;
  // per slot: type, charge, factors (padding: type 0, no charge, no factors)
  T.assign(ns, 0);
  Q.assign(ns, 0.0f);
  P12.assign(ns, 0.0f);
  P6.assign(ns, 0.0f);
  for (size_t s = 0; s < ns; ++s) {
    if (perm[s] == kPad) continue;
    T[s] = uint32_t(type_index[perm[s]]);
    Q[s] = float(q[perm[s]]);
    if (!f12.empty()) P12[s] = float(f12[T[s]]), P6[s] = float(f6[T[s]]);
  }
  X.assign(ns, 0.0f);
  Y.assign(ns, 0.0f);
  Z.assign(ns, 0.0f);
  CT12.clear();
  CT6.clear();
  if (f12.empty() && c12 && c6 && 8 * ntypes * ns * sizeof(float) <= (size_t(64) << 20)) {
    CT12.assign(ntypes * ns, 0.0f);
    CT6.assign(ntypes * ns, 0.0f);
    for (size_t t = 0; t < ntypes; ++t)
      for (size_t sl = 0; sl < ns; ++sl) {
        if (perm[sl] == kPad) continue;
        CT12[t * ns + sl] = c12[t * ntypes + T[sl]];
        CT6[t * ns + sl] = c6[t * ntypes + T[sl]];
      }
  }
  // chunks of i clusters, 4 per worker, taken by whichever worker is free: their entries shrink linearly from the first
  // to the last (chunk t gets a share ∝ nch − t), so the step does not end waiting on a big chunk on a slow core
  const size_t ne = entries.size(), nch = nth == 1 ? 1 : 4 * size_t(nth);
  tb.resize(nch);
  te.resize(nch);
  const double whole = double(nch) * double(nch + 1) / 2;
  auto bound = [&](size_t t) {   // the entry that starts chunk t
    const double before = double(t) * double(nch) - double(t) * double(t - (t > 0 ? 1 : 0)) / 2;   // Σ_{k<t} (nch − k)
    return size_t(double(ne) * before / whole);
  };
  for (size_t t = 0; t < nch; ++t) {
    tb[t] = size_t(std::lower_bound(first.begin(), first.end() - 1, t == 0 ? 0 : bound(t)) - first.begin());
    te[t] = t + 1 == nch ? nclusters : size_t(std::lower_bound(first.begin(), first.end() - 1, bound(t + 1)) - first.begin());
  }
  pe.assign(nch, {});
  pf.assign(nch, {});
}

void List::prune(double radius, ThreadPool& pool) {
  using namespace simd4;
  const F4 r2 = set1(float(radius * radius));
  std::atomic<size_t> next{0};
  pool.run(size_t(pool.size()), [&](int, size_t, size_t) {
    for (size_t t; (t = next.fetch_add(1, std::memory_order_relaxed)) < tb.size();) {
      auto& E = pe[t];
      auto& F = pf[t];
      E.clear();
      F.assign(te[t] - tb[t] + 1, 0);
      for (size_t ci = tb[t]; ci < te[t]; ++ci) {
        const size_t si = 4 * ci;
        for (size_t e = first[ci]; e < first[ci + 1]; ++e) {
          const Entry& en = entries[e];
          const size_t sj = 4 * size_t(en.cj);
          const double* sq = &qvec[3 * size_t(en.code)];
          const double* oi = &org[3 * ci];
          const double* oj = &org[3 * size_t(en.cj)];
          const float sx = float(sq[0] + oj[0] - oi[0]), sy = float(sq[1] + oj[1] - oi[1]), sz = float(sq[2] + oj[2] - oi[2]);
          const F4 xj = load(&X[sj]) + set1(sx), yj = load(&Y[sj]) + set1(sy), zj = load(&Z[sj]) + set1(sz);
          unsigned in = 0;   // the lanes within the radius
          const unsigned both = unsigned(en.mask) | unsigned(en.ex);
          for (int a = 0; a < 4; ++a) {
            const unsigned rowbits = (both >> (4 * a)) & 15u;
            if (!rowbits) continue;
            const F4 dx = xj - set1(X[si + a]), dy = yj - set1(Y[si + a]), dz = zj - set1(Z[si + a]);
            const F4 d2 = fma(dx, dx, fma(dy, dy, dz * dz));
            in |= (bits(d2 < r2) & rowbits) << (4 * a);
          }
          if (in) E.push_back({en.cj, en.code, uint16_t(en.mask & in), uint16_t(en.ex & in), uint16_t(en.ex14 & in)});
        }
        F[ci - tb[t] + 1] = uint32_t(E.size());
      }
    }
  });
}

void List::rescale(const Cell& cell) {
  for (size_t i = 0; i < off.size(); i += 3) {
    if (perm[i / 3] == kPad) continue;
    const Vec3 v = cell.a * ioff[i] + cell.b * ioff[i + 1] + cell.c * ioff[i + 2];
    off[i] = v[0]; off[i + 1] = v[1]; off[i + 2] = v[2];
  }
  for (size_t c = 0; c < qint.size(); c += 3) {
    const Vec3 v = cell.a * qint[c] + cell.b * qint[c + 1] + cell.c * qint[c + 2];
    qvec[c] = v[0]; qvec[c + 1] = v[1]; qvec[c + 2] = v[2];
  }
  cells = cell;
}

void List::positions(const std::vector<double>& x, ThreadPool& pool) {
  pool.run_slices(perm.size(), pool.slices(), [&](size_t, size_t sb, size_t se) {
    for (size_t s = sb; s < se; ++s) {
      // padding sits on its cluster's first atom (its pairs are masked)
      const size_t r = perm[s] != kPad ? s : (s & ~size_t(3));
      const size_t a = 3 * size_t(perm[r]);
      const double* o = &org[3 * (s / 4)];
      X[s] = float(x[a] + off[3 * r] - o[0]);
      Y[s] = float(x[a + 1] + off[3 * r + 1] - o[1]);
      Z[s] = float(x[a + 2] + off[3 * r + 2] - o[2]);
    }
  });
}

namespace {

using namespace simd4;

// the rational fits of nbnxm_fit.hpp, by Horner's rule
template <size_t N>
inline F4 horner(const float (&c)[N], F4 t) {
  F4 p = set1(c[N - 1]);
  for (size_t k = N - 1; k-- > 0;) p = fma(p, t, set1(c[k]));
  return p;
}
inline F4 ewald_gv(F4 t) { return horner(kGvP, t) / horner(kGvQ, t); }
inline F4 ewald_gf(F4 t) { return horner(kGfP, t) / horner(kGfQ, t); }

// kCoul: 0 none; 1 Ewald real space, 2 damped shifted force (rational fits, κ r_c ≤ 4); 3, 4 the same by erfc (any κ).
// kE: energies too (a step that reports them), else forces and virial only.
template <int kLj, bool kComb, int kCoul, bool kE>
void kernel_t(const List& L, const Params& p, int tw, float* fxo, float* fyo, float* fzo, Sums& out) {
  const F4 rc2 = set1(p.rc2), one = set1(1.0f), rcn = set1(p.rcn), rc6 = set1(p.rc6);
  const bool shifted = p.shift != 0;
  const F4 kap2 = set1(-p.kappa * p.kappa), epk = set1(kErfcP * p.kappa), kc = set1(p.kc), e0 = set1(p.e0), f0 = set1(p.f0),
           rcv = set1(p.rc);
  const F4 kap = set1(p.kappa), kap3 = set1(p.kappa * p.kappa * p.kappa), ku = set1(p.kappa * p.kappa / kFitUMax), f14 = set1(p.f14);
  const F4 twelve = set1(kLj == 2 ? 9.0f : 12.0f), six = set1(6.0f);
  // erfc by Abramowitz & Stegun 7.1.26 (|error| < 1.5e-7, the form LAMMPS uses without tables)
  const F4 a1 = set1(kErfcA[0]), a2 = set1(kErfcA[1]), a3 = set1(kErfcA[2]), a4 = set1(kErfcA[3]), a5 = set1(kErfcA[4]);
  const float* X = L.X.data();
  const float* Y = L.Y.data();
  const float* Z = L.Z.data();
  const float* Q = L.Q.data();
  const float* P12 = L.P12.data();
  const float* P6 = L.P6.data();
  const uint32_t* T = L.T.data();
  // The virial from the forces: Σ d ⊗ f over pairs (d = x_j + S − x_i, f on j) is Σ x ⊗ F over slots plus Σ S ⊗ (the
  // forces on j atoms taken through image S), the latter kept per image code
  const size_t ncodes = L.qint.size() / 3;
  std::vector<F4> fsh(3 * ncodes, zero());
  size_t zero_code = ncodes;
  for (size_t c = 0; c < ncodes; ++c)
    if (L.qint[3 * c] == 0 && L.qint[3 * c + 1] == 0 && L.qint[3 * c + 2] == 0) zero_code = c;
  const size_t cb = L.tb[size_t(tw)], ce = L.te[size_t(tw)];
  const Entry* ents = L.pe[size_t(tw)].data();
  const uint32_t* pfi = L.pf[size_t(tw)].data();
  const bool ct = !L.CT12.empty();
  const size_t nsl = L.perm.size();
  for (size_t ci = cb; ci < ce; ++ci) {
    const size_t si = 4 * ci;
    F4 fxa[4], fya[4], fza[4];
    for (int a = 0; a < 4; ++a) fxa[a] = fya[a] = fza[a] = zero();
    F4 ev = zero(), ec = zero();
    const F4 xi = load(X + si), yi = load(Y + si), zi = load(Z + si);
    const double* oi = &L.org[3 * ci];
    F4 qi = zero(), p12i = zero(), p6i = zero();
    if constexpr (kCoul != 0) qi = load(Q + si) * set1(p.kcoul);
    if constexpr (kComb) p12i = load(P12 + si), p6i = load(P6 + si);
    for (size_t e = pfi[ci - cb]; e < pfi[ci - cb + 1]; ++e) {
      const Entry& en = ents[e];
      const size_t sj = 4 * size_t(en.cj);
      // the shift: image plus the difference of the cluster origins, d = X_j − (X_i − shift)
      const double* sq = &L.qvec[3 * size_t(en.code)];
      const double* oj = &L.org[3 * size_t(en.cj)];
      const F4 xs = xi - set1(float(sq[0] + oj[0] - oi[0])), ys = yi - set1(float(sq[1] + oj[1] - oi[1])), zs = zi - set1(float(sq[2] + oj[2] - oi[2]));
      const F4 xj = load(X + sj), yj = load(Y + sj), zj = load(Z + sj);
      F4 qj = zero(), p12j = zero(), p6j = zero();
      if constexpr (kCoul != 0) qj = load(Q + sj);
      if constexpr (kComb) p12j = load(P12 + sj), p6j = load(P6 + sj);
      F4 fjx = zero(), fjy = zero(), fjz = zero();
      // kX: the entry has excluded lanes (a cluster pair along the diagonal, mostly); the others skip that work
      auto row = [&](auto ac, auto xc) {
        constexpr int A = decltype(ac)::value;
        constexpr bool kX = decltype(xc)::value;
        const unsigned bits = (en.mask >> (4 * A)) & 15u, xbits = kX ? (en.ex >> (4 * A)) & 15u : 0u;
        if (!(bits | xbits)) return;
        const F4 dx = xj - lane<A>(xs), dy = yj - lane<A>(ys), dz = zj - lane<A>(zs);
        F4 r2 = fma(dx, dx, fma(dy, dy, dz * dz));
        const M4 inr = r2 < rc2;
        // m: the pairs that interact (Lennard-Jones and Coulomb); mall: those and the excluded pairs (the Ewald / DSF
        // correction alone: Coulomb with the bare 1/r taken out, a 1-4 pair keeping f14 of it)
        const M4 m = mask_bits(bits) & inr, mall = kX ? mask_bits(bits | xbits) & inr : m;
        r2 = select(mall, r2, one);
        const F4 rinv = rsqrt(r2), r2i = rinv * rinv, r6i = r2i * r2i * r2i;
        F4 c12, c6;
        if constexpr (kComb) {
          c12 = lane<A>(p12i) * p12j;
          c6 = lane<A>(p6i) * p6j;
        } else if (ct) {
          const size_t at = size_t(T[si + A]) * nsl + sj;
          c12 = load(L.CT12.data() + at);
          c6 = load(L.CT6.data() + at);
        } else {
          const uint32_t tj[4] = {T[sj], T[sj + 1], T[sj + 2], T[sj + 3]};
          const size_t rw = size_t(T[si + A]) * p.nt;
          c12 = gather(p.c12 + rw, tj);
          c6 = gather(p.c6 + rw, tj);
        }
        F4 rn = r6i * r6i;   // r⁻¹² (12-6) or r⁻⁹ (9-6)
        if constexpr (kLj == 2) rn = r6i * r2i * rinv;
        const F4 Aa = c12 * rn, Bb = c6 * r6i;
        F4 fr = fnma(six, Bb, twelve * Aa) * r2i;
        if constexpr (kX) fr = masked(m, fr);
        if constexpr (kE) {
          F4 vlj = Aa - Bb;
          if (shifted) vlj = vlj - fnma(c6, rc6, c12 * rcn);
          ev = ev + masked(m, vlj);
        }
        // the weight of the bare 1/r: 1 for a pair that interacts, f14 for an excluded 1-4 pair, 0 for the others
        F4 wr = rinv;
        if constexpr (kX) wr = select(m, one, masked(mask_bits((en.ex14 >> (4 * A)) & 15u), f14)) * rinv;
        if constexpr (kCoul == 1 || kCoul == 2) {
          const F4 qq = lane<A>(qi) * qj;
          const F4 t = ku * r2;
          const F4 Fc = fma(kap3, ewald_gf(t), wr * r2i);   // −(dE/dr)/r
          if constexpr (kCoul == 1) {
            fr = fma(qq, Fc, fr);
            if constexpr (kE) ec = ec + masked(mall, qq * fnma(kap, ewald_gv(t), wr));
          } else {
            fr = fma(qq, fnma(f0, rinv, Fc), fr);
            if constexpr (kE) ec = ec + masked(mall, qq * fma(f0, r2 * rinv - rcv, fnma(kap, ewald_gv(t), wr) - e0));
          }
        }
        if constexpr (kCoul == 3 || kCoul == 4) {
          const F4 qq = lane<A>(qi) * qj;
          const F4 r = r2 * rinv;
          const F4 ex = exp_neg(kap2 * r2);
          const F4 t = one / fma(epk, r, one);
          const F4 poly = t * fma(t, fma(t, fma(t, fma(t, a5, a4), a3), a2), a1);
          const F4 E0 = poly * ex * rinv;                // erfc(κr)/r
          const F4 Fc0 = fma(kc, ex, E0) * r2i;          // −(dE/dr)/r
          const F4 E = E0 - (rinv - wr), Fc = Fc0 - (rinv - wr) * r2i;   // less the part of 1/r an excluded pair drops
          if constexpr (kCoul == 3) {
            if constexpr (kE) ec = ec + masked(mall, qq * E);
            fr = fma(qq, Fc, fr);
          } else {
            if constexpr (kE) ec = ec + masked(mall, qq * fma(f0, r - rcv, E - e0));
            fr = fma(qq, fnma(f0, rinv, Fc), fr);
          }
        }
        fr = masked(mall, fr);
        const F4 fx = fr * dx, fy = fr * dy, fz = fr * dz;
        fjx = fjx + fx; fjy = fjy + fy; fjz = fjz + fz;
        fxa[A] = fxa[A] - fx; fya[A] = fya[A] - fy; fza[A] = fza[A] - fz;
      };
      if (en.ex) {
        row(std::integral_constant<int, 0>{}, std::true_type{});
        row(std::integral_constant<int, 1>{}, std::true_type{});
        row(std::integral_constant<int, 2>{}, std::true_type{});
        row(std::integral_constant<int, 3>{}, std::true_type{});
      } else {
        row(std::integral_constant<int, 0>{}, std::false_type{});
        row(std::integral_constant<int, 1>{}, std::false_type{});
        row(std::integral_constant<int, 2>{}, std::false_type{});
        row(std::integral_constant<int, 3>{}, std::false_type{});
      }
      store(fxo + sj, load(fxo + sj) + fjx);
      store(fyo + sj, load(fyo + sj) + fjy);
      store(fzo + sj, load(fzo + sj) + fjz);
      if (en.code != zero_code) {
        F4* fs = &fsh[3 * size_t(en.code)];
        fs[0] = fs[0] + fjx; fs[1] = fs[1] + fjy; fs[2] = fs[2] + fjz;
      }
    }
    for (int a = 0; a < 4; ++a) {
      fxo[si + a] += hsum(fxa[a]);
      fyo[si + a] += hsum(fya[a]);
      fzo[si + a] += hsum(fza[a]);
    }
    out.evdw += hsum(ev);
    out.ecoul += hsum(ec);
  }
  // the forces taken through each image, for the virial
  if (out.fsh.size() < 3 * ncodes) out.fsh.resize(3 * ncodes, 0.0);
  for (size_t c = 0; c < ncodes; ++c) {
    if (c == zero_code) continue;
    for (int k = 0; k < 3; ++k) out.fsh[3 * c + size_t(k)] += hsum(fsh[3 * c + size_t(k)]);
  }
}

template <int kLj, bool kComb, bool kE>
void by_coul(const List& L, const Params& p, int tw, float* fx, float* fy, float* fz, Sums& out) {
  const bool fit = p.kappa * p.kappa * p.rc2 <= kFitUMax;   // u = (κr)² within the fits
  if (p.coul == 1) fit ? kernel_t<kLj, kComb, 1, kE>(L, p, tw, fx, fy, fz, out) : kernel_t<kLj, kComb, 3, kE>(L, p, tw, fx, fy, fz, out);
  else if (p.coul == 2) fit ? kernel_t<kLj, kComb, 2, kE>(L, p, tw, fx, fy, fz, out) : kernel_t<kLj, kComb, 4, kE>(L, p, tw, fx, fy, fz, out);
  else kernel_t<kLj, kComb, 0, kE>(L, p, tw, fx, fy, fz, out);
}

template <int kLj, bool kComb>
void by_energy(const List& L, const Params& p, int tw, float* fx, float* fy, float* fz, Sums& out) {
  if (p.energies) by_coul<kLj, kComb, true>(L, p, tw, fx, fy, fz, out);
  else by_coul<kLj, kComb, false>(L, p, tw, fx, fy, fz, out);
}

}  // namespace

void pair_virial(const List& L, const double* F, size_t sb, size_t se, double W[9]) {
  for (size_t s = sb; s < se; ++s) {
    const double* o = &L.org[3 * (s / 4)];
    const double r[3] = {double(L.X[s]) + o[0], double(L.Y[s]) + o[1], double(L.Z[s]) + o[2]};
    for (int u = 0; u < 3; ++u)
      for (int w = 0; w < 3; ++w) W[3 * u + w] += r[u] * F[3 * s + size_t(w)];
  }
}

void shift_virial(const List& L, const std::vector<double>& fsh, double W[9]) {
  for (size_t c = 0; 3 * c < fsh.size() && 3 * c < L.qvec.size(); ++c)
    for (int u = 0; u < 3; ++u)
      for (int w = 0; w < 3; ++w) W[3 * u + w] += L.qvec[3 * c + size_t(u)] * fsh[3 * c + size_t(w)];
}

namespace {
// AVX2 and FMA when built in and the processor reports them; CAPS_SIMD=sse2 keeps the four-lane kernel (to compare),
// CAPS_SIMD=avx2 takes the eight-lane one whatever the processor reports (an emulator that runs AVX2 without saying so)
bool use_avx2() {
  static const bool on = [] {
#if (defined(__x86_64__) || defined(_M_X64)) && defined(__GNUC__)
    const char* env = std::getenv("CAPS_SIMD");
    if (env && std::string(env) == "sse2") return false;
    if (env && std::string(env) == "avx2") return avx2_built();
    return avx2_built() && __builtin_cpu_supports("avx2") && __builtin_cpu_supports("fma");
#else
    return false;
#endif
  }();
  return on;
}
}  // namespace

const char* instruction_set() {
  if (use_avx2()) return "AVX2";
#if defined(CAPS_SIMD4_NEON)
  return "NEON";
#elif defined(CAPS_SIMD4_SSE)
  return "SSE2";
#else
  return "scalar";
#endif
}

void kernel(const List& L, const Params& p, int tw, float* fx, float* fy, float* fz, Sums& out) {
  if (use_avx2() && kernel_avx2(L, p, tw, fx, fy, fz, out)) return;
  if (p.lj == 2) {
    if (p.comb) by_energy<2, true>(L, p, tw, fx, fy, fz, out);
    else by_energy<2, false>(L, p, tw, fx, fy, fz, out);
  } else {
    if (p.comb) by_energy<1, true>(L, p, tw, fx, fy, fz, out);
    else by_energy<1, false>(L, p, tw, fx, fy, fz, out);
  }
}

}  // namespace caps::nbnxm
