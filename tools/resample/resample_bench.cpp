/*
 * Standalone checks and benchmark for p2t::ResamplePolyline.
 *
 * Uses only the C++ standard library and the poly2tri sources shipped in this
 * repository (no third-party test framework). Built via run.sh.
 *
 * Covers:
 *   - four sample contours (circle, rectangle, sharp-corner star,
 *     self-intersecting bowtie): spacing error, corner retention, area error,
 *     per-call runtime;
 *   - spacing == 0 (and other invalid arguments) must throw;
 *   - open polyline support, exact endpoint/corner preservation;
 *   - configurable spacing and corner retention angle;
 *   - determinism and absence of duplicate output points;
 *   - triangulation of the original inputs is unchanged vertex by vertex.
 */

#include <poly2tri/poly2tri.h>

#include <algorithm>
#include <cassert>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <list>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

int g_failures = 0;

#define CHECK(cond)                                                                            \
  do {                                                                                         \
    if (!(cond)) {                                                                             \
      std::cerr << "CHECK FAILED: " #cond " (line " << __LINE__ << ")\n";                      \
      ++g_failures;                                                                            \
    }                                                                                          \
  } while (0)

const double PI = 3.14159265358979323846;

double Distance(const p2t::Point& a, const p2t::Point& b)
{
  return std::hypot(a.x - b.x, a.y - b.y);
}

// Mirror of the turning-angle rule used by p2t::ResamplePolyline, so the
// expected corner counts are computed independently here.
double TurningAngle(const std::vector<p2t::Point>& p, std::size_t i, bool closed)
{
  const std::size_t n = p.size();
  if (!closed && (i == 0 || i == n - 1)) {
    return 0.0;
  }
  const std::size_t prv = (i == 0) ? n - 1 : i - 1;
  const std::size_t nxt = (i + 1 == n) ? (closed ? 0 : n - 1) : i + 1;
  if (prv == i || nxt == i) {
    return 0.0;
  }
  p2t::Point u = p[i] - p[prv];
  p2t::Point v = p[nxt] - p[i];
  if (!(u.Normalize() > 0.0) || !(v.Normalize() > 0.0)) {
    return 0.0;
  }
  double c = p2t::Dot(u, v);
  c = std::max(-1.0, std::min(1.0, c));
  return std::acos(c);
}

std::size_t CountCorners(const std::vector<p2t::Point>& p, double threshold, bool closed)
{
  std::size_t count = 0;
  for (std::size_t i = 0; i < p.size(); ++i) {
    if (TurningAngle(p, i, closed) >= threshold) {
      ++count;
    }
  }
  return count;
}

bool ContainsPoint(const std::vector<p2t::Point>& pts, const p2t::Point& q, double tol)
{
  for (const p2t::Point& r : pts) {
    if (Distance(r, q) <= tol) {
      return true;
    }
  }
  return false;
}

double SignedArea(const std::vector<p2t::Point>& p)
{
  double a = 0.0;
  for (std::size_t i = 0; i < p.size(); ++i) {
    const p2t::Point& q = p[(i + 1) % p.size()];
    a += p[i].x * q.y - q.x * p[i].y;
  }
  return 0.5 * a;
}

double ContourLength(const std::vector<p2t::Point>& p, bool closed)
{
  double len = 0.0;
  const std::size_t m = closed ? p.size() : p.size() - 1;
  for (std::size_t i = 0; i < m; ++i) {
    len += Distance(p[i], p[(i + 1) % p.size()]);
  }
  return len;
}

std::vector<p2t::Point> MakeCircle(std::size_t segments, double radius)
{
  std::vector<p2t::Point> p;
  for (std::size_t i = 0; i < segments; ++i) {
    const double a = 2.0 * PI * static_cast<double>(i) / static_cast<double>(segments);
    p.emplace_back(radius * std::cos(a), radius * std::sin(a));
  }
  return p;
}

// Simple concave decagon (non-self-intersecting five-tip star), CCW order.
std::vector<p2t::Point> MakeStar()
{
  std::vector<p2t::Point> p;
  const double outer = 2.0;
  const double inner = 0.9;
  for (std::size_t i = 0; i < 10; ++i) {
    const double a = -PI / 2.0 + PI * static_cast<double>(i) / 5.0;
    const double r = (i % 2 == 0) ? outer : inner;
    p.emplace_back(r * std::cos(a), r * std::sin(a));
  }
  return p;
}

// Bowtie: edges (0,0)->(4,4) and (4,0)->(0,4) cross, so it self-intersects.
std::vector<p2t::Point> MakeBowtie()
{
  return { p2t::Point(0, 0), p2t::Point(4, 4), p2t::Point(4, 0), p2t::Point(0, 4) };
}

std::vector<p2t::Point> MakeRectangle()
{
  return { p2t::Point(0, 0), p2t::Point(6, 0), p2t::Point(6, 3), p2t::Point(0, 3) };
}

struct Row {
  std::string name;
  std::vector<p2t::Point> in;
  double spacing;
  bool self_intersecting;
};

void PrintMetricRow(const Row& row, const p2t::ResampleConfig& cfg)
{
  const std::size_t corners_in = CountCorners(row.in, cfg.corner_angle_rad, cfg.closed);

  // Timing: many repeats, allocation included.
  const int repeats = 2000;
  std::vector<p2t::Point> out;
  const auto t0 = std::chrono::steady_clock::now();
  for (int k = 0; k < repeats; ++k) {
    out = p2t::ResamplePolyline(row.in, cfg);
  }
  const auto t1 = std::chrono::steady_clock::now();
  const double us =
    std::chrono::duration<double, std::micro>(t1 - t0).count() / repeats;

  // Spacing error over every output gap.
  const std::size_t gaps = cfg.closed ? out.size() : out.size() - 1;
  double max_rel = 0.0;
  double sum_rel = 0.0;
  double min_d = 1e300;
  double max_d = 0.0;
  for (std::size_t i = 0; i < gaps; ++i) {
    const double d = Distance(out[i], out[(i + 1) % out.size()]);
    min_d = std::min(min_d, d);
    max_d = std::max(max_d, d);
    const double rel = std::abs(d - cfg.spacing) / cfg.spacing;
    max_rel = std::max(max_rel, rel);
    sum_rel += rel;
  }
  const double mean_rel = sum_rel / static_cast<double>(gaps);

  // Corner retention: every input corner must appear verbatim in the output.
  std::size_t kept = 0;
  const double tol = 1e-9;
  for (std::size_t i = 0; i < row.in.size(); ++i) {
    if (TurningAngle(row.in, i, cfg.closed) >= cfg.corner_angle_rad &&
        ContainsPoint(out, row.in[i], tol)) {
      ++kept;
    }
  }

  const double area_in = SignedArea(row.in);
  const double area_out = SignedArea(out);
  const double area_abs = std::abs(area_out - area_in);
  const double area_rel =
    (std::abs(area_in) > 1e-12) ? area_abs / std::abs(area_in) : area_abs;

  std::cout << std::left << std::setw(16) << row.name << std::right << std::fixed
            << std::setprecision(3) << std::setw(8) << cfg.spacing << std::setw(9)
            << row.in.size() << " -> " << std::left << std::setw(5) << out.size()
            << std::right << std::setw(7) << corners_in << " -> " << std::left
            << std::setw(5) << kept << std::scientific << std::setprecision(2)
            << std::setw(10) << max_rel << std::setw(10) << mean_rel << std::setw(6)
            << min_d << ".." << std::left << std::setw(9) << max_d << std::right
            << std::setw(11) << area_rel << std::fixed << std::setprecision(2)
            << std::setw(9) << us << "\n";
}

// Triangulate and compare every output vertex against the expected coordinates.
void CheckTriangulationUnchanged()
{
  // Equivalent of the repository's BasicTest.
  std::list<p2t::Point> store;
  store.emplace_back(0, 0);
  store.emplace_back(1, 0);
  store.emplace_back(1, 1);
  std::vector<p2t::Point*> polyline{ &store.front(), &*(++store.begin()),
                                     &*(--store.end()) };
  p2t::CDT cdt{ polyline };
  cdt.Triangulate();
  const auto tris = cdt.GetTriangles();
  CHECK(tris.size() == 1);
  const p2t::Point* got[3] = { tris[0]->GetPoint(0), tris[0]->GetPoint(1),
                               tris[0]->GetPoint(2) };
  for (const p2t::Point* expected : polyline) {
    bool found = false;
    for (const p2t::Point* g : got) {
      if (g->x == expected->x && g->y == expected->y) {
        found = true;
      }
    }
    CHECK(found);
  }
  CHECK(p2t::IsDelaunay(tris));
}

// Structural triangulation checks: a simple n-gon yields n-2 triangles, every
// triangle vertex comes from the input contour and every input vertex is used.
// p2t::IsDelaunay itself is not used here because perfectly cocircular points
// (circle samples) trip its strict floating-point predicate.
void CheckResampledPolygonTriangulates(std::vector<p2t::Point> p)
{
  if (SignedArea(p) < 0.0) {
    std::reverse(p.begin(), p.end());
  }
  std::list<p2t::Point> store(p.begin(), p.end());
  std::vector<p2t::Point*> polyline;
  for (p2t::Point& q : store) {
    polyline.push_back(&q);
  }
  p2t::CDT cdt{ polyline };
  cdt.Triangulate();
  const auto tris = cdt.GetTriangles();
  CHECK(tris.size() == polyline.size() - 2);
  std::size_t used = 0;
  for (const p2t::Point* v : polyline) {
    bool referenced = false;
    for (p2t::Triangle* t : tris) {
      for (int k = 0; k < 3; ++k) {
        const p2t::Point* g = t->GetPoint(k);
        CHECK(ContainsPoint(p, *g, 0.0));
        if (g == v) {
          referenced = true;
        }
      }
    }
    if (referenced) {
      ++used;
    }
  }
  CHECK(used == polyline.size());
}

} // namespace

int main()
{
  const double corner_threshold = PI / 4.0;

  std::cout << "=== Arc-length resampling (corner threshold = 45 deg, closed) ===\n";
  std::cout << std::left << std::setw(16) << "contour" << std::right << std::setw(8)
            << "spacing" << std::setw(14) << "points in->out" << std::setw(13)
            << "corners in->kept" << std::setw(10) << "max ds/s" << std::setw(10)
            << "mean ds/s" << std::setw(17) << "gap range" << std::setw(11)
            << "area err" << std::setw(10) << "us/call"
            << "\n";

  const Row rows[] = {
    { "circle(32,r5)", MakeCircle(32, 5.0), 0.3, false },
    { "rectangle(6x3)", MakeRectangle(), 1.0, false },
    { "star(5 tips)", MakeStar(), 0.2, false },
    { "bowtie(self-x)", MakeBowtie(), 0.4, true },
  };
  for (const Row& row : rows) {
    p2t::ResampleConfig cfg;
    cfg.spacing = row.spacing;
    cfg.corner_angle_rad = corner_threshold;
    cfg.closed = true;
    PrintMetricRow(row, cfg);
  }

  std::cout << "\n=== Self checks ===\n";

  // 1. spacing == 0 must throw (also negative spacing).
  try {
    p2t::ResampleConfig cfg;
    cfg.spacing = 0.0;
    (void)p2t::ResamplePolyline(MakeRectangle(), cfg);
    CHECK(false);
  } catch (const std::invalid_argument& e) {
    std::cout << "[ok] spacing=0 rejected: " << e.what() << "\n";
  }
  try {
    p2t::ResampleConfig cfg;
    cfg.spacing = -1.0;
    (void)p2t::ResamplePolyline(MakeRectangle(), cfg);
    CHECK(false);
  } catch (const std::invalid_argument&) {
    std::cout << "[ok] negative spacing rejected\n";
  }

  // 2. Degenerate inputs must throw.
  try {
    p2t::ResampleConfig cfg;
    (void)p2t::ResamplePolyline({ p2t::Point(0, 0) }, cfg);
    CHECK(false);
  } catch (const std::invalid_argument&) {
    std::cout << "[ok] single-point closed contour rejected\n";
  }
  try {
    p2t::ResampleConfig cfg;
    cfg.closed = false;
    (void)p2t::ResamplePolyline({ p2t::Point(0, 0), p2t::Point(0, 0) }, cfg);
    CHECK(false);
  } catch (const std::invalid_argument&) {
    std::cout << "[ok] zero-length segment rejected\n";
  }
  {
    // Closed contour carrying a repeated closing vertex must be accepted.
    const std::vector<p2t::Point> sq{ p2t::Point(0, 0), p2t::Point(2, 0),
                                      p2t::Point(2, 2), p2t::Point(0, 2),
                                      p2t::Point(0, 0) };
    p2t::ResampleConfig cfg;
    cfg.spacing = 1.0;
    const auto out = p2t::ResamplePolyline(sq, cfg);
    CHECK(out.size() == 8);
    CHECK(Distance(out.front(), out.back()) > 1e-12);
    std::cout << "[ok] repeated closing vertex on closed contour accepted\n";
  }

  // 3. Open polyline: endpoints and sharp corners kept exactly, full length
  //    covered once with strictly positive gaps.
  {
    const std::vector<p2t::Point> line{ p2t::Point(0, 0), p2t::Point(2, 0),
                                        p2t::Point(2, 1), p2t::Point(4, 1) };
    p2t::ResampleConfig cfg;
    cfg.spacing = 0.5;
    cfg.corner_angle_rad = PI / 4.0;
    cfg.closed = false;
    const auto out = p2t::ResamplePolyline(line, cfg);
    CHECK(Distance(out.front(), line.front()) < 1e-12);
    CHECK(Distance(out.back(), line.back()) < 1e-12);
    for (const p2t::Point& q : line) {
      CHECK(ContainsPoint(out, q, 1e-9));
    }
    double covered = 0.0;
    for (std::size_t i = 1; i < out.size(); ++i) {
      const double d = Distance(out[i - 1], out[i]);
      CHECK(d > 0.0);
      covered += d;
    }
    CHECK(std::abs(covered - ContourLength(line, false)) < 1e-9);
    CHECK(CountCorners(out, cfg.corner_angle_rad, false) == 2);
    std::cout << "[ok] open polyline: endpoints/corners exact, length preserved\n";
  }

  // 4. Corner retention angle is configurable: a 100 deg threshold drops the
  //    rectangle's 90 deg corners; a 30 deg threshold keeps them. The checked
  //    corners avoid the half-perimeter midpoint, which the equidistant grid
  //    hits regardless of corner handling.
  {
    const auto rect = MakeRectangle();
    p2t::ResampleConfig loose;
    loose.spacing = 0.7;
    loose.corner_angle_rad = 100.0 * PI / 180.0;
    const auto dropped = p2t::ResamplePolyline(rect, loose);
    CHECK(!ContainsPoint(dropped, p2t::Point(6, 0), 1e-9));
    CHECK(!ContainsPoint(dropped, p2t::Point(0, 3), 1e-9));

    p2t::ResampleConfig tight;
    tight.spacing = 1.0;
    tight.corner_angle_rad = 30.0 * PI / 180.0;
    const auto kept = p2t::ResamplePolyline(rect, tight);
    for (const p2t::Point& q : rect) {
      CHECK(ContainsPoint(kept, q, 1e-9));
    }
    std::cout << "[ok] corner retention angle configurable (dropped at 100deg, "
                 "kept at 30deg)\n";
  }

  // 5. Closed contours: no duplicate points, all gaps positive, every input
  //    corner retained verbatim on all four sample shapes.
  for (const Row& row : rows) {
    p2t::ResampleConfig cfg;
    cfg.spacing = row.spacing;
    cfg.corner_angle_rad = corner_threshold;
    const auto out = p2t::ResamplePolyline(row.in, cfg);
    for (std::size_t i = 0; i < out.size(); ++i) {
      // Consecutive points (including the closing gap) must be distinct.
      const double gap = Distance(out[i], out[(i + 1) % out.size()]);
      CHECK(gap > 1e-12);
      // Non-self-intersecting contours must not visit the same spot twice.
      // Self-intersecting ones may legitimately pass the crossing point more
      // than once (at different arc positions), so the global check is skipped.
      if (!row.self_intersecting) {
        for (std::size_t j = i + 1; j < out.size(); ++j) {
          CHECK(Distance(out[i], out[j]) > 1e-12);
        }
      }
    }
    for (std::size_t i = 0; i < row.in.size(); ++i) {
      if (TurningAngle(row.in, i, true) >= corner_threshold) {
        CHECK(ContainsPoint(out, row.in[i], 1e-9));
      }
    }
    // Determinism.
    const auto out2 = p2t::ResamplePolyline(row.in, cfg);
    CHECK(out.size() == out2.size());
    for (std::size_t i = 0; i < out.size(); ++i) {
      CHECK(Distance(out[i], out2[i]) < 1e-15);
    }
  }
  std::cout << "[ok] closed contours: no duplicates, corners verbatim, deterministic\n";

  // 6. Spacing larger than the perimeter still yields a valid contour.
  {
    p2t::ResampleConfig cfg;
    cfg.spacing = 1000.0;
    const auto out = p2t::ResamplePolyline(MakeRectangle(), cfg);
    CHECK(out.size() == MakeRectangle().size());
    std::cout << "[ok] spacing > perimeter handled\n";
  }

  // 7. Original triangulation output is unchanged vertex by vertex, and
  //    resampled simple polygons still triangulate.
  CheckTriangulationUnchanged();
  std::cout << "[ok] original CDT output unchanged vertex by vertex\n";
  CheckResampledPolygonTriangulates(p2t::ResamplePolyline(MakeStar(),
                                                          p2t::ResampleConfig{}));
  p2t::ResampleConfig circle_cfg;
  circle_cfg.spacing = 0.5;
  CheckResampledPolygonTriangulates(
    p2t::ResamplePolyline(MakeCircle(32, 5.0), circle_cfg));
  std::cout << "[ok] resampled simple polygons triangulate (n-2 triangles, "
               "full vertex coverage)\n";

  if (g_failures == 0) {
    std::cout << "\nALL CHECKS PASSED\n";
    return 0;
  }
  std::cout << "\n" << g_failures << " CHECK(S) FAILED\n";
  return 1;
}
