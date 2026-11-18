// Vertex-clustering mesh decimation self-check / benchmark.
//
// Standalone script (not part of the library build). It:
//   1. builds four sample meshes (sphere, cube, plane, holed plate -- the
//      holed plate is triangulated with this repo's poly2tri CDT),
//   2. decimates each with vertex clustering at three grid cell sizes
//      (quantize vertices to a grid, merge each cell into one centroid
//      representative, rebuild triangles, drop degenerate/duplicate faces),
//   3. reports vertices / faces / surface area / volume before and after,
//   4. verifies: surface-area change < 5%, no edge shared by more than
//      two faces (non-manifold), and decimation of a ~100k-face mesh in
//      under 500 ms.
//
// Build (from repo root):
//   g++ -O2 -std=c++17 -I. tools/vertex_clustering_decimation.cpp \
//       poly2tri/common/shapes.cc poly2tri/sweep/*.cc -o /tmp/vcd
// Run:
//   /tmp/vcd

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <map>
#include <set>
#include <string>
#include <unordered_map>
#include <vector>

#include "poly2tri/poly2tri.h"

namespace {

struct Vec3 {
  double x, y, z;
};

struct Mesh {
  std::vector<Vec3> verts;
  std::vector<std::array<uint32_t, 3>> tris;
};

// ---------------------------------------------------------------------------
// Vertex clustering decimation
// ---------------------------------------------------------------------------

struct CellKey {
  int64_t x, y, z;
  bool operator==(const CellKey& o) const {
    return x == o.x && y == o.y && z == o.z;
  }
};

struct CellKeyHash {
  size_t operator()(const CellKey& k) const {
    uint64_t h = 1469598103934665603ull;
    auto mix = [&h](uint64_t v) {
      h ^= v;
      h *= 1099511628211ull;
    };
    mix(static_cast<uint64_t>(k.x));
    mix(static_cast<uint64_t>(k.y));
    mix(static_cast<uint64_t>(k.z));
    return static_cast<size_t>(h);
  }
};

// Quantize vertices to grid cells of side `cell`, merge each cell's vertices
// into one representative (the cell centroid), rebuild the triangles and drop
// degenerate (collapsed) and duplicate faces.
Mesh VertexCluster(const Mesh& in, double cell) {
  struct Accum {
    double sx = 0, sy = 0, sz = 0;
    uint32_t count = 0;
    uint32_t rep = 0;  // index into `reps`
  };
  std::unordered_map<CellKey, Accum, CellKeyHash> cells;
  cells.reserve(in.verts.size());

  const double inv = 1.0 / cell;
  std::vector<uint32_t> vert_cell(in.verts.size());
  for (size_t i = 0; i < in.verts.size(); ++i) {
    const Vec3& v = in.verts[i];
    CellKey key{static_cast<int64_t>(std::floor(v.x * inv)),
                static_cast<int64_t>(std::floor(v.y * inv)),
                static_cast<int64_t>(std::floor(v.z * inv))};
    auto it = cells.find(key);
    if (it == cells.end()) {
      Accum a;
      a.rep = static_cast<uint32_t>(cells.size());
      it = cells.emplace(key, a).first;
    }
    Accum& a = it->second;
    a.sx += v.x;
    a.sy += v.y;
    a.sz += v.z;
    ++a.count;
    vert_cell[i] = a.rep;
  }

  // Representatives ordered by first-seen cell id.
  std::vector<Vec3> reps(cells.size());
  for (const auto& kv : cells) {
    const Accum& a = kv.second;
    reps[a.rep] = Vec3{a.sx / a.count, a.sy / a.count, a.sz / a.count};
  }

  Mesh out;
  out.verts = std::move(reps);
  out.tris.reserve(in.tris.size());
  std::set<std::array<uint32_t, 3>> seen;  // duplicate-face filter
  for (const auto& t : in.tris) {
    uint32_t a = vert_cell[t[0]];
    uint32_t b = vert_cell[t[1]];
    uint32_t c = vert_cell[t[2]];
    if (a == b || b == c || a == c) continue;  // degenerate after merge
    std::array<uint32_t, 3> sorted{a, b, c};
    std::sort(sorted.begin(), sorted.end());
    if (!seen.insert(sorted).second) continue;  // identical face exists
    out.tris.push_back({a, b, c});
  }

  // Compact: keep only vertices referenced by a surviving face.
  std::vector<uint32_t> remap(out.verts.size(), UINT32_MAX);
  uint32_t next = 0;
  for (auto& t : out.tris)
    for (auto& i : t)
      if (remap[i] == UINT32_MAX) remap[i] = next++;
  std::vector<Vec3> used(next);
  for (size_t i = 0; i < out.verts.size(); ++i)
    if (remap[i] != UINT32_MAX) used[remap[i]] = out.verts[i];
  for (auto& t : out.tris)
    for (auto& i : t) i = remap[i];
  out.verts = std::move(used);
  return out;
}

// ---------------------------------------------------------------------------
// Mesh statistics
// ---------------------------------------------------------------------------

struct Stats {
  size_t verts = 0;
  size_t faces = 0;
  double area = 0;
  double volume = 0;      // signed volume (meaningful for closed meshes)
  int non_manifold = 0;   // edges shared by more than two faces
};

Stats Measure(const Mesh& m) {
  Stats s;
  s.verts = m.verts.size();
  s.faces = m.tris.size();
  std::unordered_map<uint64_t, int> edge_use;
  edge_use.reserve(m.tris.size() * 3);
  for (const auto& t : m.tris) {
    const Vec3& a = m.verts[t[0]];
    const Vec3& b = m.verts[t[1]];
    const Vec3& c = m.verts[t[2]];
    double ux = b.x - a.x, uy = b.y - a.y, uz = b.z - a.z;
    double vx = c.x - a.x, vy = c.y - a.y, vz = c.z - a.z;
    double cx = uy * vz - uz * vy;
    double cy = uz * vx - ux * vz;
    double cz = ux * vy - uy * vx;
    s.area += 0.5 * std::sqrt(cx * cx + cy * cy + cz * cz);
    s.volume += (a.x * (b.y * c.z - b.z * c.y) +
                 a.y * (b.z * c.x - b.x * c.z) +
                 a.z * (b.x * c.y - b.y * c.x)) /
                6.0;
    for (int e = 0; e < 3; ++e) {
      uint32_t p = t[e], q = t[(e + 1) % 3];
      if (p > q) std::swap(p, q);
      ++edge_use[(static_cast<uint64_t>(p) << 32) | q];
    }
  }
  for (const auto& kv : edge_use)
    if (kv.second > 2) ++s.non_manifold;
  return s;
}

double BboxDiag(const Mesh& m) {
  Vec3 lo{1e300, 1e300, 1e300}, hi{-1e300, -1e300, -1e300};
  for (const auto& v : m.verts) {
    lo.x = std::min(lo.x, v.x);
    lo.y = std::min(lo.y, v.y);
    lo.z = std::min(lo.z, v.z);
    hi.x = std::max(hi.x, v.x);
    hi.y = std::max(hi.y, v.y);
    hi.z = std::max(hi.z, v.z);
  }
  return std::sqrt((hi.x - lo.x) * (hi.x - lo.x) + (hi.y - lo.y) * (hi.y - lo.y) +
                   (hi.z - lo.z) * (hi.z - lo.z));
}

// ---------------------------------------------------------------------------
// Sample meshes
// ---------------------------------------------------------------------------

Mesh MakeSphere(double r, int nseg, int nring) {
  Mesh m;
  const double pi = std::acos(-1.0);
  m.verts.push_back({0, 0, r});  // north pole
  for (int j = 1; j < nring; ++j) {
    double th = pi * j / nring;
    for (int i = 0; i < nseg; ++i) {
      double ph = 2 * pi * i / nseg;
      m.verts.push_back(
          {r * std::sin(th) * std::cos(ph), r * std::sin(th) * std::sin(ph),
           r * std::cos(th)});
    }
  }
  uint32_t south = static_cast<uint32_t>(m.verts.size());
  m.verts.push_back({0, 0, -r});
  auto ring = [&](int j, int i) {  // j in [1, nring-1]
    return static_cast<uint32_t>(1 + (j - 1) * nseg + (i % nseg));
  };
  for (int i = 0; i < nseg; ++i)
    m.tris.push_back({0, ring(1, i + 1), ring(1, i)});
  for (int j = 1; j < nring - 1; ++j)
    for (int i = 0; i < nseg; ++i) {
      m.tris.push_back({ring(j, i), ring(j + 1, i), ring(j + 1, i + 1)});
      m.tris.push_back({ring(j, i), ring(j + 1, i + 1), ring(j, i + 1)});
    }
  for (int i = 0; i < nseg; ++i)
    m.tris.push_back({ring(nring - 1, i), ring(nring - 1, i + 1), south});
  return m;
}

// Weld helper: merge vertices whose coordinates agree to 1e-9.
struct Welder {
  Mesh* m;
  std::map<std::array<int64_t, 3>, uint32_t> index;
  uint32_t add(double x, double y, double z) {
    std::array<int64_t, 3> key{static_cast<int64_t>(std::llround(x * 1e9)),
                               static_cast<int64_t>(std::llround(y * 1e9)),
                               static_cast<int64_t>(std::llround(z * 1e9))};
    auto it = index.find(key);
    if (it != index.end()) return it->second;
    uint32_t id = static_cast<uint32_t>(m->verts.size());
    m->verts.push_back({x, y, z});
    index.emplace(key, id);
    return id;
  }
};

Mesh MakeCube(double half, int n) {
  Mesh m;
  Welder w{&m};
  // Six faces, each an n x n grid of quads. axes: (fixed axis, sign)
  for (int axis = 0; axis < 3; ++axis)
    for (int sign = -1; sign <= 1; sign += 2) {
      int u = (axis + 1) % 3, v = (axis + 2) % 3;
      auto point = [&](int i, int j) {
        double p[3] = {0, 0, 0};
        p[axis] = sign * half;
        p[u] = -half + 2.0 * half * i / n;
        p[v] = -half + 2.0 * half * j / n;
        return w.add(p[0], p[1], p[2]);
      };
      for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j) {
          uint32_t a = point(i, j), b = point(i + 1, j);
          uint32_t c = point(i + 1, j + 1), d = point(i, j + 1);
          if (sign > 0) {
            m.tris.push_back({a, b, c});
            m.tris.push_back({a, c, d});
          } else {
            m.tris.push_back({a, c, b});
            m.tris.push_back({a, d, c});
          }
        }
    }
  return m;
}

Mesh MakePlane(double half, int n) {
  Mesh m;
  for (int i = 0; i <= n; ++i)
    for (int j = 0; j <= n; ++j)
      m.verts.push_back(
          {-half + 2.0 * half * i / n, -half + 2.0 * half * j / n, 0});
  auto at = [&](int i, int j) {
    return static_cast<uint32_t>(i * (n + 1) + j);
  };
  for (int i = 0; i < n; ++i)
    for (int j = 0; j < n; ++j) {
      m.tris.push_back({at(i, j), at(i + 1, j), at(i + 1, j + 1)});
      m.tris.push_back({at(i, j), at(i + 1, j + 1), at(i, j + 1)});
    }
  return m;
}

// Square plate with a circular hole, triangulated by this repo's poly2tri.
Mesh MakeHoledPlate() {
  std::vector<p2t::Point*> keep;  // ownership
  auto P = [&keep](double x, double y) {
    keep.push_back(new p2t::Point(x, y));
    return keep.back();
  };
  std::vector<p2t::Point*> outer{P(-1, -1), P(1, -1), P(1, 1), P(-1, 1)};
  p2t::CDT cdt(outer);
  const double pi = std::acos(-1.0);
  std::vector<p2t::Point*> hole;
  for (int i = 0; i < 64; ++i) {
    double a = 2 * pi * i / 64;
    hole.push_back(P(0.4 * std::cos(a), 0.4 * std::sin(a)));
  }
  cdt.AddHole(hole);
  const int n = 97;  // Steiner grid
  for (int i = 0; i < n; ++i)
    for (int j = 0; j < n; ++j) {
      double x = -0.97 + 1.94 * i / (n - 1);
      double y = -0.97 + 1.94 * j / (n - 1);
      if (x * x + y * y < 0.45 * 0.45) continue;  // keep clear of the hole
      cdt.AddPoint(P(x, y));
    }
  cdt.Triangulate();

  Mesh m;
  std::unordered_map<const p2t::Point*, uint32_t> id;
  for (auto* t : cdt.GetTriangles()) {
    std::array<uint32_t, 3> tri;
    for (int k = 0; k < 3; ++k) {
      const p2t::Point* p = t->GetPoint(k);
      auto it = id.find(p);
      if (it == id.end()) {
        uint32_t nid = static_cast<uint32_t>(m.verts.size());
        m.verts.push_back({p->x, p->y, 0});
        it = id.emplace(p, nid).first;
      }
      tri[k] = it->second;
    }
    m.tris.push_back(tri);
  }
  for (auto* p : keep) delete p;
  return m;
}

// ---------------------------------------------------------------------------
// Driver
// ---------------------------------------------------------------------------

bool g_all_ok = true;

void Report(const std::string& name, const Mesh& m, bool closed) {
  const Stats before = Measure(m);
  const double diag = BboxDiag(m);
  std::printf(
      "\n### %s  (before: V=%zu F=%zu area=%.6f vol=%.6f%s)\n\n",
      name.c_str(), before.verts, before.faces, before.area,
      std::fabs(before.volume), closed ? "" : ", open mesh");
  std::printf(
      "| cell size | V after | F after | area after | dArea %% | vol after | "
      "dVol %% | non-manifold | time |\n");
  std::printf(
      "|-----------|---------|---------|------------|---------|-----------|"
      "---------|--------------|------|\n");
  for (double frac : {128.0, 64.0, 32.0}) {
    const double cell = diag / frac;
    const auto t0 = std::chrono::steady_clock::now();
    Mesh out = VertexCluster(m, cell);
    const auto t1 = std::chrono::steady_clock::now();
    const Stats after = Measure(out);
    const double ms =
        std::chrono::duration<double, std::milli>(t1 - t0).count();
    const double darea =
        100.0 * (after.area - before.area) / before.area;
    const double dvol =
        before.volume != 0
            ? 100.0 * (std::fabs(after.volume) - std::fabs(before.volume)) /
                  std::fabs(before.volume)
            : 0.0;
    const bool ok = std::fabs(darea) < 5.0 && after.non_manifold == 0;
    g_all_ok = g_all_ok && ok;
    std::printf(
        "| diag/%-6.0f | %7zu | %7zu | %10.6f | %+7.3f | %9.6f | %+7.3f | "
        "%12d | %5.1f ms | %s\n",
        frac, after.verts, after.faces, after.area, darea,
        std::fabs(after.volume), dvol, after.non_manifold, ms,
        ok ? "OK" : "FAIL <--");
  }
}

}  // namespace

int main() {
  Report("Sphere (r=1, 128x64)", MakeSphere(1.0, 128, 64), true);
  Report("Cube (side=2, 128x128/face)", MakeCube(1.0, 128), true);
  Report("Plane (2x2, 96x96)", MakePlane(1.0, 96), false);
  Report("Holed plate (poly2tri CDT)", MakeHoledPlate(), false);

  // Performance: ~100k-face sphere must decimate in < 500 ms.
  Mesh big = MakeSphere(1.0, 320, 160);
  const Stats bs = Measure(big);
  const double cell = BboxDiag(big) / 32.0;
  const auto t0 = std::chrono::steady_clock::now();
  Mesh out = VertexCluster(big, cell);
  const auto t1 = std::chrono::steady_clock::now();
  const Stats as = Measure(out);
  const double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
  const bool perf_ok = bs.faces >= 100000 && ms < 500.0 && as.non_manifold == 0;
  g_all_ok = g_all_ok && perf_ok;
  std::printf(
      "\n### Performance: %zu faces -> %zu faces in %.1f ms "
      "(limit 500 ms, non-manifold=%d) %s\n",
      bs.faces, as.faces, ms, as.non_manifold, perf_ok ? "OK" : "FAIL <--");

  std::printf("\nOverall: %s\n", g_all_ok ? "PASS" : "FAIL");
  return g_all_ok ? 0 : 1;
}
