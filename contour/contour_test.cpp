// Standalone contour-extraction (marching triangles) test for poly2tri.
//
// Pipeline (three stages):
//   1. Per triangle, classify vertex scalar values against the iso level to
//      find crossed edges. Degenerate (collinear, zero-area) triangles are
//      detected explicitly and skipped.
//   2. On each crossed edge, linearly interpolate the intersection point.
//   3. Chain the segments into polylines (closed rings or open lines) and
//      simplify them with a Douglas-Peucker threshold.
//
// Configurable: iso value, simplification threshold, closed/open handling.
//
// Build (no build-script changes needed):
//   g++ -std=c++14 -O2 -I. contour/contour_test.cpp
//       poly2tri/common/*.cc poly2tri/sweep/*.cc -o build/contour_test
// Run:
//   ./build/contour_test

#include <poly2tri/poly2tri.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <functional>
#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <utility>
#include <vector>

namespace {

constexpr double kPi = 3.14159265358979323846;

struct Vec2 {
  double x = 0.0;
  double y = 0.0;
};

// ---------------------------------------------------------------------------
// Configuration
// ---------------------------------------------------------------------------
struct ContourOptions {
  double iso = 0.0;           // iso value (等值)
  double simplify_eps = 0.0;  // Douglas-Peucker threshold (简化阈值)
  bool closed_only = false;   // keep only closed rings (是否闭合: 只保留闭合线)
  bool force_close = false;   // close open polylines by joining the ends
};

struct ContourLine {
  std::vector<Vec2> points;  // raw interpolated intersection points
  bool closed = false;
};

struct ContourResult {
  std::vector<ContourLine> lines;        // after simplification / filtering
  size_t raw_point_count = 0;            // total interpolated intersections
  size_t degenerate_triangles = 0;       // skipped collinear triangles
  double max_area_error = 0.0;           // max |A_simpl - A_raw| / A_raw (closed)
};

// ---------------------------------------------------------------------------
// Stage 1+2: per-triangle crossing test and edge interpolation
// ---------------------------------------------------------------------------
using HeightFn = std::function<double(double, double)>;

// Quantized key so intersection points computed from both sides of a shared
// edge land on the same bucket.
using Key = std::pair<int64_t, int64_t>;
Key MakeKey(const Vec2& p) {
  return {static_cast<int64_t>(std::llround(p.x * 1e9)),
          static_cast<int64_t>(std::llround(p.y * 1e9))};
}

struct Segment {
  Vec2 a, b;
};

// Extract one segment per non-degenerate, crossed triangle.
std::vector<Segment> ExtractSegments(const std::vector<p2t::Triangle*>& tris,
                                     const HeightFn& height, double iso,
                                     size_t* degenerate_count) {
  std::vector<Segment> segments;
  *degenerate_count = 0;
  for (p2t::Triangle* tri : tris) {
    const p2t::Point* p[3] = {tri->GetPoint(0), tri->GetPoint(1),
                              tri->GetPoint(2)};
    // Explicit degenerate-triangle handling: collinear vertices (zero area).
    const double d = (p[1]->x - p[0]->x) * (p[2]->y - p[0]->y) -
                     (p[2]->x - p[0]->x) * (p[1]->y - p[0]->y);
    if (std::abs(d) <= 1e-12) {
      ++(*degenerate_count);
      continue;
    }
    double v[3];
    bool below[3];
    for (int i = 0; i < 3; ++i) {
      v[i] = height(p[i]->x, p[i]->y);
      // A vertex exactly on the iso level is nudged above it so that no
      // vertex sits on the level; crossings then interpolate to (within
      // 1 ulp of) the vertex and merge into a single point when chained.
      if (v[i] == iso) v[i] = std::nextafter(iso, HUGE_VAL);
      below[i] = v[i] < iso;
    }
    Vec2 crossing[2];
    int n = 0;
    for (int e = 0; e < 3 && n < 2; ++e) {
      const int i = e, j = (e + 1) % 3;
      if (below[i] == below[j]) continue;
      // Linear interpolation of the crossing position on edge (i, j).
      const double t = (iso - v[i]) / (v[j] - v[i]);
      crossing[n++] = {p[i]->x + t * (p[j]->x - p[i]->x),
                       p[i]->y + t * (p[j]->y - p[i]->y)};
    }
    if (n == 2) segments.push_back({crossing[0], crossing[1]});
  }
  return segments;
}

// ---------------------------------------------------------------------------
// Stage 3a: chain segments into polylines
// ---------------------------------------------------------------------------
std::vector<ContourLine> ChainSegments(const std::vector<Segment>& segments) {
  std::map<Key, std::vector<size_t>> adjacency;
  std::vector<bool> used(segments.size(), false);
  for (size_t i = 0; i < segments.size(); ++i) {
    if (MakeKey(segments[i].a) == MakeKey(segments[i].b)) {
      used[i] = true;  // zero-length segment (e.g. at an on-iso vertex)
      continue;
    }
    adjacency[MakeKey(segments[i].a)].push_back(i);
    adjacency[MakeKey(segments[i].b)].push_back(i);
  }
  std::vector<ContourLine> lines;

  auto take_unused = [&](const Key& k, size_t* out) {
    auto it = adjacency.find(k);
    if (it == adjacency.end()) return false;
    for (size_t idx : it->second) {
      if (!used[idx]) {
        *out = idx;
        return true;
      }
    }
    return false;
  };

  for (size_t seed = 0; seed < segments.size(); ++seed) {
    if (used[seed]) continue;
    used[seed] = true;
    std::vector<Vec2> pts{segments[seed].a, segments[seed].b};
    // Extend the front, then the back, until no unused segment connects.
    for (int pass = 0; pass < 2; ++pass) {
      for (;;) {
        const Vec2 end = pass == 0 ? pts.back() : pts.front();
        size_t next;
        if (!take_unused(MakeKey(end), &next)) break;
        used[next] = true;
        const Vec2 a = segments[next].a, b = segments[next].b;
        const Vec2 other =
            (MakeKey(a) == MakeKey(end)) ? b : a;
        if (pass == 0)
          pts.push_back(other);
        else
          pts.insert(pts.begin(), other);
      }
    }
    ContourLine line;
    line.closed = MakeKey(pts.front()) == MakeKey(pts.back());
    if (line.closed && pts.size() > 1) pts.pop_back();  // drop duplicate join
    if (pts.size() < (line.closed ? 3 : 2)) continue;   // degenerate chain
    line.points = std::move(pts);
    lines.push_back(std::move(line));
  }
  return lines;
}

// ---------------------------------------------------------------------------
// Stage 3b: Douglas-Peucker simplification
// ---------------------------------------------------------------------------
double PointLineDistance(const Vec2& p, const Vec2& a, const Vec2& b) {
  const double dx = b.x - a.x, dy = b.y - a.y;
  const double len2 = dx * dx + dy * dy;
  if (len2 <= 0.0) return std::hypot(p.x - a.x, p.y - a.y);
  return std::abs(dy * p.x - dx * p.y + b.x * a.y - b.y * a.x) /
         std::sqrt(len2);
}

void DouglasPeucker(const std::vector<Vec2>& pts, size_t first, size_t last,
                    double eps, std::vector<bool>* keep) {
  if (last <= first + 1) return;
  double dmax = 0.0;
  size_t index = first;
  for (size_t i = first + 1; i < last; ++i) {
    const double d = PointLineDistance(pts[i], pts[first], pts[last]);
    if (d > dmax) {
      dmax = d;
      index = i;
    }
  }
  if (dmax > eps) {
    (*keep)[index] = true;
    DouglasPeucker(pts, first, index, eps, keep);
    DouglasPeucker(pts, index, last, eps, keep);
  }
}

std::vector<Vec2> Simplify(const std::vector<Vec2>& pts, bool closed,
                           double eps) {
  if (pts.size() < 3 || eps <= 0.0) return pts;
  // For closed rings, anchor the split at the point farthest from the
  // centroid so the simplification is not biased by the arbitrary start.
  std::vector<Vec2> work = pts;
  if (closed) {
    Vec2 c;
    for (const Vec2& p : pts) {
      c.x += p.x / pts.size();
      c.y += p.y / pts.size();
    }
    size_t far = 0;
    double best = -1.0;
    for (size_t i = 0; i < pts.size(); ++i) {
      const double d = std::hypot(pts[i].x - c.x, pts[i].y - c.y);
      if (d > best) {
        best = d;
        far = i;
      }
    }
    std::rotate(work.begin(), work.begin() + far, work.end());
    work.push_back(work.front());
  }
  std::vector<bool> keep(work.size(), false);
  keep.front() = keep.back() = true;
  DouglasPeucker(work, 0, work.size() - 1, eps, &keep);
  std::vector<Vec2> out;
  for (size_t i = 0; i < work.size(); ++i)
    if (keep[i]) out.push_back(work[i]);
  if (closed && out.size() > 1) out.pop_back();  // restore ring representation
  return out;
}

// ---------------------------------------------------------------------------
// Driver
// ---------------------------------------------------------------------------
double PolygonArea(const std::vector<Vec2>& pts) {
  double a = 0.0;
  for (size_t i = 0; i < pts.size(); ++i) {
    const Vec2& p = pts[i];
    const Vec2& q = pts[(i + 1) % pts.size()];
    a += p.x * q.y - q.x * p.y;
  }
  return 0.5 * a;
}

ContourResult ExtractContours(const std::vector<p2t::Triangle*>& tris,
                              const HeightFn& height,
                              const ContourOptions& opt) {
  ContourResult result;
  const std::vector<Segment> segments =
      ExtractSegments(tris, height, opt.iso, &result.degenerate_triangles);
  std::vector<ContourLine> lines = ChainSegments(segments);
  for (ContourLine& line : lines) {
    result.raw_point_count += line.points.size();
    if (opt.closed_only && !line.closed) continue;
    const double raw_area =
        (line.closed || opt.force_close) ? std::abs(PolygonArea(line.points))
                                         : 0.0;
    ContourLine out;
    out.closed = line.closed || opt.force_close;
    out.points = Simplify(line.points, line.closed, opt.simplify_eps);
    if (out.closed && out.points.size() >= 3 && raw_area > 1e-14) {
      const double err =
          std::abs(std::abs(PolygonArea(out.points)) - raw_area) / raw_area;
      result.max_area_error = std::max(result.max_area_error, err);
    }
    result.lines.push_back(std::move(out));
  }
  return result;
}

// ---------------------------------------------------------------------------
// Mesh construction via poly2tri (boundary square + interior Steiner grid)
// ---------------------------------------------------------------------------
struct Mesh {
  std::vector<p2t::Point*> owned;
  std::unique_ptr<p2t::CDT> cdt;  // owns the triangles; must outlive them
  std::vector<p2t::Triangle*> triangles;
  Mesh() = default;
  Mesh(Mesh&&) = default;
  Mesh& operator=(Mesh&&) = default;
  Mesh(const Mesh&) = delete;
  Mesh& operator=(const Mesh&) = delete;
  ~Mesh() {
    cdt.reset();  // delete triangles before their vertices
    for (p2t::Point* p : owned) delete p;
  }
};

Mesh BuildGridMesh(int n) {  // n x n interior Steiner points in [0,1]^2
  Mesh mesh;
  std::vector<p2t::Point*> boundary;
  for (const auto& c : {Vec2{0, 0}, Vec2{1, 0}, Vec2{1, 1}, Vec2{0, 1}}) {
    boundary.push_back(new p2t::Point(c.x, c.y));
  }
  for (p2t::Point* p : boundary) mesh.owned.push_back(p);
  mesh.cdt.reset(new p2t::CDT(boundary));
  for (int i = 1; i <= n; ++i) {
    for (int j = 1; j <= n; ++j) {
      auto* p = new p2t::Point(static_cast<double>(i) / (n + 1),
                               static_cast<double>(j) / (n + 1));
      mesh.owned.push_back(p);
      mesh.cdt->AddPoint(p);
    }
  }
  mesh.cdt->Triangulate();
  mesh.triangles = mesh.cdt->GetTriangles();
  return mesh;
}

// Snapshot of every vertex coordinate, sorted, for exact comparison.
std::vector<std::pair<double, double>> VertexSnapshot(
    const std::vector<p2t::Triangle*>& tris) {
  std::vector<std::pair<double, double>> out;
  for (p2t::Triangle* t : tris)
    for (int i = 0; i < 3; ++i)
      out.emplace_back(t->GetPoint(i)->x, t->GetPoint(i)->y);
  std::sort(out.begin(), out.end());
  return out;
}

// ---------------------------------------------------------------------------
// Test scaffolding
// ---------------------------------------------------------------------------
int g_failures = 0;

void Check(bool cond, const std::string& what) {
  if (!cond) {
    ++g_failures;
    std::cout << "  FAIL: " << what << '\n';
  }
}

struct TerrainReport {
  std::string name;
  size_t contours = 0;
  size_t closed = 0;
  size_t raw_points = 0;
  size_t simplified_points = 0;
  size_t degenerate = 0;
  double max_area_error = 0.0;
  double micros = 0.0;
};

TerrainReport RunTerrain(const std::string& name,
                         const std::vector<p2t::Triangle*>& tris,
                         const HeightFn& height, const ContourOptions& opt) {
  const auto before = VertexSnapshot(tris);
  const auto t0 = std::chrono::steady_clock::now();
  const ContourResult r = ExtractContours(tris, height, opt);
  const auto t1 = std::chrono::steady_clock::now();
  // The triangulation output must be untouched by extraction, vertex by
  // vertex (原三角化输出逐顶点不变).
  Check(VertexSnapshot(tris) == before,
        name + ": triangulation vertices changed by extraction");

  TerrainReport rep;
  rep.name = name;
  rep.contours = r.lines.size();
  rep.raw_points = r.raw_point_count;
  rep.degenerate = r.degenerate_triangles;
  rep.max_area_error = r.max_area_error;
  rep.micros = std::chrono::duration<double, std::micro>(t1 - t0).count();
  for (const ContourLine& l : r.lines) {
    if (l.closed) ++rep.closed;
    rep.simplified_points += l.points.size();
  }
  return rep;
}

void PrintReport(const std::vector<TerrainReport>& reports) {
  std::cout <<
      "\nterrain      contours  closed  raw_pts  simpl_pts  degen  max_area_err "
      " time_us\n";
  for (const TerrainReport& r : reports) {
    printf("%-12s %8zu  %6zu  %7zu  %9zu  %5zu  %12.3e  %8.1f\n",
           r.name.c_str(), r.contours, r.closed, r.raw_points,
           r.simplified_points, r.degenerate, r.max_area_error, r.micros);
  }
}

}  // namespace

int main() {
  std::vector<TerrainReport> reports;

  // -- Terrain 1: tilted plane h = 0.3x + 0.7y, iso 0.5 ---------------------
  // Expect exactly one open, perfectly straight contour.
  {
    Mesh mesh = BuildGridMesh(10);
    const HeightFn h = [](double x, double y) { return 0.3 * x + 0.7 * y; };
    ContourOptions opt;
    opt.iso = 0.5;
    opt.simplify_eps = 1e-3;
    const ContourResult r =
        ExtractContours(mesh.triangles, h, opt);  // checked via RunTerrain too
    TerrainReport rep = RunTerrain("plane", mesh.triangles, h, opt);
    Check(rep.contours == 1, "plane: expected exactly 1 contour");
    Check(rep.closed == 0, "plane: contour must be open");
    // Straightness: every raw intersection lies on 0.3x + 0.7y = 0.5.
    double max_dev = 0.0;
    for (const ContourLine& l : r.lines)
      for (const Vec2& p : l.points)
        max_dev = std::max(max_dev, std::abs(0.3 * p.x + 0.7 * p.y - 0.5));
    Check(max_dev < 1e-9, "plane: contour not straight");
    // closed_only must filter this open line out (config check).
    ContourOptions closed_opt = opt;
    closed_opt.closed_only = true;
    Check(ExtractContours(mesh.triangles, h, closed_opt).lines.empty(),
          "plane: closed_only should drop the open contour");
    reports.push_back(rep);
  }

  // -- Terrain 2: undulating h = 0.4 sin(2 pi x) sin(2 pi y), iso 0.15 ------
  // Two positive peaks fully inside the domain -> two closed rings.
  {
    Mesh mesh = BuildGridMesh(30);
    const HeightFn h = [](double x, double y) {
      return 0.4 * std::sin(2 * kPi * x) * std::sin(2 * kPi * y);
    };
    ContourOptions opt;
    opt.iso = 0.15;
    opt.simplify_eps = 1e-3;
    TerrainReport rep = RunTerrain("undulating", mesh.triangles, h, opt);
    Check(rep.contours == 2, "undulating: expected 2 contours");
    Check(rep.closed == 2, "undulating: both contours must be closed");
    Check(rep.max_area_error < 1e-3, "undulating: area error too large");
    reports.push_back(rep);
  }

  // -- Terrain 3: three Gaussian peaks, iso 0.3 ------------------------------
  {
    Mesh mesh = BuildGridMesh(30);
    const HeightFn h = [](double x, double y) {
      auto gauss = [](double x, double y, double cx, double cy, double a) {
        const double s = 0.08;
        const double dx = x - cx, dy = y - cy;
        return a * std::exp(-(dx * dx + dy * dy) / (2 * s * s));
      };
      return gauss(x, y, 0.30, 0.30, 1.0) + gauss(x, y, 0.70, 0.35, 0.8) +
             gauss(x, y, 0.50, 0.75, 0.9);
    };
    ContourOptions opt;
    opt.iso = 0.3;
    opt.simplify_eps = 1e-3;
    TerrainReport rep = RunTerrain("multipeak", mesh.triangles, h, opt);
    Check(rep.contours == 3, "multipeak: expected 3 contours");
    Check(rep.closed == 3, "multipeak: all contours must be closed");
    Check(rep.max_area_error < 5e-3,
          "multipeak: area error too large");  // ~eps/radius of smallest ring
    reports.push_back(rep);
  }

  // -- Terrain 4: degenerate mesh (collinear triangle vertices) --------------
  // Hand-built fan with two zero-area triangles mixed in; the extractor must
  // skip them explicitly and still contour the rest.
  {
    std::vector<p2t::Point> pts = {
        {0.0, 0.0}, {1.0, 0.0}, {2.0, 0.0},  // collinear triple
        {0.0, 1.0}, {1.0, 1.0}, {2.0, 1.0},
        {1.0, 0.5},                          // interior
    };
    std::vector<p2t::Triangle> storage;
    storage.emplace_back(pts[0], pts[1], pts[2]);  // degenerate (collinear)
    storage.emplace_back(pts[3], pts[4], pts[5]);  // degenerate (collinear)
    storage.emplace_back(pts[0], pts[1], pts[6]);
    storage.emplace_back(pts[1], pts[2], pts[6]);
    storage.emplace_back(pts[0], pts[6], pts[3]);
    storage.emplace_back(pts[3], pts[6], pts[4]);
    storage.emplace_back(pts[6], pts[5], pts[4]);
    storage.emplace_back(pts[2], pts[5], pts[6]);
    std::vector<p2t::Triangle*> tris;
    for (p2t::Triangle& t : storage) tris.push_back(&t);
    const HeightFn h = [](double x, double y) {
      const double dx = x - 1.0, dy = y - 0.5;
      return std::exp(-(dx * dx + dy * dy) / 0.2);
    };
    ContourOptions opt;
    opt.iso = 0.5;
    opt.simplify_eps = 1e-3;
    TerrainReport rep = RunTerrain("degenerate", tris, h, opt);
    Check(rep.degenerate == 2, "degenerate: expected 2 skipped triangles");
    Check(rep.contours == 1 && rep.closed == 1,
          "degenerate: expected 1 closed contour around the bump");
    reports.push_back(rep);
  }

  PrintReport(reports);
  if (g_failures == 0) {
    std::cout << "\nAll contour checks passed.\n";
    return 0;
  }
  std::cout << '\n' << g_failures << " check(s) FAILED.\n";
  return 1;
}
