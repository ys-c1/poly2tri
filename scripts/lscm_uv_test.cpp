// Standalone test/benchmark for p2t::LscmParameterize (approximate LSCM UV
// parameterization on top of poly2tri CDT output).
//
// Builds four small meshes with poly2tri, parameterizes them, and reports:
//   - number of flipped triangles in UV space
//   - area distortion (RMS log area ratio, 0 = area preserving)
//   - RMS error against the manual unwrap (normalized XY, after optimal
//     similarity alignment)
//   - solve time in milliseconds
// Also verifies that the triangulation itself is left untouched
// (per-vertex identical) and that singular systems raise errors.
//
// Exits with a non-zero status if any check fails.

#include <poly2tri/poly2tri.h>
#include <poly2tri/common/lscm.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <map>
#include <string>
#include <vector>

namespace {

struct Mesh {
  std::string name;
  std::vector<p2t::Point*> storage;       // owns every point
  std::vector<p2t::Point*> boundary;      // outer polyline (CCW)
  std::vector<std::vector<p2t::Point*>> holes;
  std::vector<p2t::Point*> steiner;
  double minx = 0.0, miny = 0.0, norm = 1.0; // manual-unwrap frame

  ~Mesh()
  {
    for (p2t::Point* p : storage) {
      delete p;
    }
  }

  p2t::Point* Add(double x, double y)
  {
    storage.push_back(new p2t::Point(x, y));
    return storage.back();
  }

  void ComputeFrame()
  {
    double maxx = -1e300, maxy = -1e300;
    minx = miny = 1e300;
    for (const p2t::Point* p : storage) {
      minx = std::min(minx, p->x);
      miny = std::min(miny, p->y);
      maxx = std::max(maxx, p->x);
      maxy = std::max(maxy, p->y);
    }
    norm = std::max(maxx - minx, maxy - miny);
  }

  // Manual unwrap: aspect-preserving normalized XY coordinates.
  double ManualU(const p2t::Point* p) const { return (p->x - minx) / norm; }
  double ManualV(const p2t::Point* p) const { return (p->y - miny) / norm; }
};

// Rectangle boundary polyline (CCW) subdivided with spacing ds.
void RectBoundary(Mesh& m, double w, double h, double ds,
                  std::vector<p2t::Point*>& out)
{
  const int nx = std::max(1, static_cast<int>(std::lround(w / ds)));
  const int ny = std::max(1, static_cast<int>(std::lround(h / ds)));
  for (int i = 0; i < nx; ++i) out.push_back(m.Add(w * i / nx, 0));
  for (int i = 0; i < ny; ++i) out.push_back(m.Add(w, h * i / ny));
  for (int i = 0; i < nx; ++i) out.push_back(m.Add(w - w * i / nx, h));
  for (int i = 0; i < ny; ++i) out.push_back(m.Add(0, h - h * i / ny));
}

Mesh MakeSquareGrid()
{
  Mesh m;
  m.name = "square_grid";
  RectBoundary(m, 1.0, 1.0, 0.2, m.boundary);
  for (int i = 1; i < 5; ++i) {
    for (int j = 1; j < 5; ++j) {
      m.steiner.push_back(m.Add(0.2 * i, 0.2 * j));
    }
  }
  return m;
}

Mesh MakeRectGrid()
{
  Mesh m;
  m.name = "rect_grid_2x1";
  RectBoundary(m, 2.0, 1.0, 0.25, m.boundary);
  for (int i = 1; i < 8; ++i) {
    for (int j = 1; j < 4; ++j) {
      m.steiner.push_back(m.Add(0.25 * i, 0.25 * j));
    }
  }
  return m;
}

Mesh MakeLShape()
{
  Mesh m;
  m.name = "l_shape";
  // L-shaped polygon (concave), CCW.
  const double pts[6][2] = {
    {0, 0}, {2, 0}, {2, 1}, {1, 1}, {1, 2}, {0, 2},
  };
  const double ds = 0.25;
  for (int e = 0; e < 6; ++e) {
    const double x0 = pts[e][0], y0 = pts[e][1];
    const double x1 = pts[(e + 1) % 6][0], y1 = pts[(e + 1) % 6][1];
    const double len = std::hypot(x1 - x0, y1 - y0);
    const int steps = std::max(1, static_cast<int>(std::lround(len / ds)));
    for (int i = 0; i < steps; ++i) {
      const double t = static_cast<double>(i) / steps;
      m.boundary.push_back(m.Add(x0 + t * (x1 - x0), y0 + t * (y1 - y0)));
    }
  }
  // Interior grid points strictly inside the L.
  const double eps = 1e-9;
  for (int i = 1; i < 8; ++i) {
    for (int j = 1; j < 8; ++j) {
      const double x = 0.25 * i, y = 0.25 * j;
      const bool inside = (x < 1.0 - eps && y < 2.0 - eps) ||
                          (x < 2.0 - eps && y < 1.0 - eps);
      if (inside) {
        m.steiner.push_back(m.Add(x, y));
      }
    }
  }
  return m;
}

Mesh MakeSquareHole()
{
  Mesh m;
  m.name = "square_hole";
  RectBoundary(m, 2.0, 2.0, 0.25, m.boundary);
  // Centered square hole [0.75, 1.25]^2.
  std::vector<p2t::Point*> hole;
  const double lo = 0.75, hi = 1.25;
  hole.push_back(m.Add(lo, lo));
  hole.push_back(m.Add(1.0, lo));
  hole.push_back(m.Add(hi, lo));
  hole.push_back(m.Add(hi, 1.0));
  hole.push_back(m.Add(hi, hi));
  hole.push_back(m.Add(1.0, hi));
  hole.push_back(m.Add(lo, hi));
  hole.push_back(m.Add(lo, 1.0));
  m.holes.push_back(hole);
  // Interior grid points between outer boundary and hole.
  const double eps = 1e-9;
  for (int i = 1; i < 8; ++i) {
    for (int j = 1; j < 8; ++j) {
      const double x = 0.25 * i, y = 0.25 * j;
      const bool in_hole = x > lo - eps && x < hi + eps && y > lo - eps &&
                           y < hi + eps;
      if (!in_hole) {
        m.steiner.push_back(m.Add(x, y));
      }
    }
  }
  return m;
}

// RMS error of computed UVs against the manual unwrap after the optimal
// similarity transform (translation + rotation + uniform scale).
double ProcrustesRms(const std::vector<p2t::Point*>& points,
                     const p2t::LscmResult& uv, const Mesh& mesh)
{
  const size_t n = points.size();
  double pmu = 0, pmv = 0, mmu = 0, mmv = 0;
  for (size_t i = 0; i < n; ++i) {
    pmu += uv.u[i];
    pmv += uv.v[i];
    mmu += mesh.ManualU(points[i]);
    mmv += mesh.ManualV(points[i]);
  }
  pmu /= n;
  pmv /= n;
  mmu /= n;
  mmv /= n;
  // Complex least squares: manual ~= a * computed + t.
  double num_re = 0, num_im = 0, den = 0;
  for (size_t i = 0; i < n; ++i) {
    const double pu = uv.u[i] - pmu, pv = uv.v[i] - pmv;
    const double mu = mesh.ManualU(points[i]) - mmu;
    const double mv = mesh.ManualV(points[i]) - mmv;
    num_re += mu * pu + mv * pv;
    num_im += mu * pv - mv * pu;
    den += pu * pu + pv * pv;
  }
  const double a_re = num_re / den, a_im = num_im / den;
  double accum = 0;
  for (size_t i = 0; i < n; ++i) {
    const double ru = a_re * uv.u[i] - a_im * uv.v[i] + (mmu - a_re * pmu + a_im * pmv);
    const double rv = a_im * uv.u[i] + a_re * uv.v[i] + (mmv - a_im * pmu - a_re * pmv);
    const double du = ru - mesh.ManualU(points[i]);
    const double dv = rv - mesh.ManualV(points[i]);
    accum += du * du + dv * dv;
  }
  return std::sqrt(accum / n);
}

int g_failures = 0;

void Check(bool ok, const std::string& what)
{
  if (!ok) {
    std::printf("  FAIL: %s\n", what.c_str());
    ++g_failures;
  }
}

void RunMesh(Mesh& mesh)
{
  mesh.ComputeFrame();
  p2t::CDT cdt(mesh.boundary);
  for (const auto& hole : mesh.holes) {
    cdt.AddHole(hole);
  }
  for (p2t::Point* p : mesh.steiner) {
    cdt.AddPoint(p);
  }
  cdt.Triangulate();
  std::vector<p2t::Triangle*> triangles = cdt.GetTriangles();

  // Snapshot the triangulation to prove it is left untouched.
  std::map<const p2t::Point*, std::pair<double, double>> before;
  for (p2t::Triangle* t : triangles) {
    for (int k = 0; k < 3; ++k) {
      const p2t::Point* p = t->GetPoint(k);
      before[p] = {p->x, p->y};
    }
  }
  std::vector<p2t::Point*> points(before.size());
  {
    size_t i = 0;
    for (const auto& kv : before) {
      points[i++] = const_cast<p2t::Point*>(kv.first);
    }
  }

  p2t::LscmOptions options; // defaults: 5 iterations, boundary weight 1e3
  p2t::LscmResult uv;
  const auto t0 = std::chrono::steady_clock::now();
  uv = p2t::LscmParameterize(points, triangles, options);
  const auto t1 = std::chrono::steady_clock::now();
  const double ms =
      std::chrono::duration<double, std::milli>(t1 - t0).count();

  // Triangulation must be per-vertex identical afterwards.
  std::vector<p2t::Triangle*> after = cdt.GetTriangles();
  Check(after.size() == triangles.size(), mesh.name + ": triangle count changed");
  bool coords_same = true, conn_same = true;
  for (size_t i = 0; i < after.size(); ++i) {
    for (int k = 0; k < 3; ++k) {
      const p2t::Point* p = after[i]->GetPoint(k);
      const auto it = before.find(p);
      if (it == before.end() || it->second.first != p->x ||
          it->second.second != p->y) {
        coords_same = false;
      }
      if (after[i]->GetPoint(k) != triangles[i]->GetPoint(k)) {
        conn_same = false;
      }
    }
  }
  Check(coords_same, mesh.name + ": vertex coordinates modified");
  Check(conn_same, mesh.name + ": triangle connectivity modified");

  const double err = ProcrustesRms(points, uv, mesh);
  std::printf("%-14s tris=%3zu  flipped=%d  area_distortion=%.4f  "
              "manual_err=%.3e  solve=%.2f ms\n",
              mesh.name.c_str(), triangles.size(), uv.flipped_triangles,
              uv.area_distortion, err, ms);

  Check(uv.flipped_triangles == 0, mesh.name + ": flipped triangles");
  // Sanity bound on distortion. Concave domains (l_shape) mapped to a
  // convex boundary inherently distort near reflex corners, so the bound
  // is loose on purpose; flipped triangles are the hard invariant.
  Check(uv.area_distortion < 1.0, mesh.name + ": excessive area distortion");
  Check(err < 0.2, mesh.name + ": excessive error vs manual unwrap");
}

void RunSingularChecks()
{
  // A valid small mesh plus one isolated vertex: the system is singular
  // and must throw instead of returning bad UVs.
  std::vector<p2t::Point*> pts;
  pts.push_back(new p2t::Point(0, 0));
  pts.push_back(new p2t::Point(1, 0));
  pts.push_back(new p2t::Point(0, 1));
  p2t::CDT cdt(std::vector<p2t::Point*>{pts[0], pts[1], pts[2]});
  cdt.Triangulate();
  std::vector<p2t::Point*> with_isolated = pts;
  p2t::Point* isolated = new p2t::Point(5, 5);
  with_isolated.push_back(isolated);
  bool threw = false;
  try {
    p2t::LscmParameterize(with_isolated, cdt.GetTriangles());
  } catch (const std::runtime_error&) {
    threw = true;
  }
  Check(threw, "isolated vertex: expected std::runtime_error");
  std::printf("%-14s isolated vertex -> %s\n", "singular_case",
              threw ? "std::runtime_error (OK)" : "NO ERROR (FAIL)");

  // Degenerate (zero-area) triangle must throw as well.
  std::vector<p2t::Point*> dpts{new p2t::Point(0, 0), new p2t::Point(1, 0),
                                new p2t::Point(2, 0)};
  p2t::Triangle* degenerate = new p2t::Triangle(*dpts[0], *dpts[1], *dpts[2]);
  threw = false;
  try {
    p2t::LscmParameterize(dpts, std::vector<p2t::Triangle*>{degenerate});
  } catch (const std::runtime_error&) {
    threw = true;
  }
  Check(threw, "degenerate triangle: expected std::runtime_error");
  std::printf("%-14s degenerate triangle -> %s\n", "singular_case",
              threw ? "std::runtime_error (OK)" : "NO ERROR (FAIL)");

  delete degenerate;
  delete isolated;
  for (p2t::Point* p : pts) {
    delete p;
  }
  for (p2t::Point* p : dpts) {
    delete p;
  }

  // Configurability smoke test: no iterations + hard boundary
  // (boundary_weight <= 0) must still produce a valid parameterization.
  {
    p2t::Point* steiner = new p2t::Point(0.5, 0.5);
    std::vector<p2t::Point*> qpts{new p2t::Point(0, 0), new p2t::Point(1, 0),
                                  new p2t::Point(1, 1), new p2t::Point(0, 1)};
    p2t::CDT qcdt(qpts);
    qcdt.AddPoint(steiner);
    qcdt.Triangulate();
    std::vector<p2t::Point*> tri_pts;
    for (p2t::Triangle* t : qcdt.GetTriangles()) {
      for (int k = 0; k < 3; ++k) {
        p2t::Point* p = t->GetPoint(k);
        if (std::find(tri_pts.begin(), tri_pts.end(), p) == tri_pts.end()) {
          tri_pts.push_back(p);
        }
      }
    }
    p2t::LscmOptions opt;
    opt.iterations = 0;
    opt.boundary_weight = 0.0; // hard boundary
    p2t::LscmResult uv;
    bool ok = true;
    try {
      uv = p2t::LscmParameterize(tri_pts, qcdt.GetTriangles(), opt);
    } catch (const std::exception&) {
      ok = false;
    }
    Check(ok && uv.flipped_triangles == 0,
          "options: iterations=0, hard boundary failed");
    std::printf("%-14s iterations=0, boundary_weight=0 -> %s\n", "options",
                ok ? "valid UVs (OK)" : "FAILURE");
    delete steiner;
    for (p2t::Point* p : qpts) {
      delete p;
    }
  }
}

} // namespace

int main()
{
  {
    Mesh m = MakeSquareGrid();
    RunMesh(m);
  }
  {
    Mesh m = MakeRectGrid();
    RunMesh(m);
  }
  {
    Mesh m = MakeLShape();
    RunMesh(m);
  }
  {
    Mesh m = MakeSquareHole();
    RunMesh(m);
  }
  RunSingularChecks();
  if (g_failures != 0) {
    std::printf("LSCM UV test: %d check(s) FAILED\n", g_failures);
    return 1;
  }
  std::printf("LSCM UV test: all checks passed\n");
  return 0;
}
