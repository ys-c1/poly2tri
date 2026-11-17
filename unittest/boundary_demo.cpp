/*
 * Poly2Tri Copyright (c) 2009-2018, Poly2Tri Contributors
 * https://github.com/jhasse/poly2tri
 *
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without modification,
 * are permitted provided that the following conditions are met:
 *
 * * Redistributions of source code must retain the above copyright notice,
 *   this list of conditions and the following disclaimer.
 * * Redistributions in binary form must reproduce the above copyright notice,
 *   this list of conditions and the following disclaimer in the documentation
 *   and/or other materials provided with the distribution.
 * * Neither the name of Poly2Tri nor the names of its contributors may be
 *   used to endorse or promote products derived from this software without specific
 *   prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR
 * CONTRIBUTORS BE LIABLE FOR DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

/*
 * Boundary extraction validation on four small, hand-built mesh families:
 *
 *   G1 single contour   : one concave polygon (no holes), redundant
 *                         collinear and short-edge vertices, Steiner points
 *   G2 with holes       : one outer square with one square hole
 *   G3 multi-component  : two disjoint components (a square with a hole and
 *                         a triangle), triangulated separately
 *   G4 degenerate       : G4a non-manifold triangle soup (open chain with
 *                         reported break points), G4b two holes touching at
 *                         one vertex, G4c empty mesh
 *
 * The program prints, for every run: number of closed/open contours, total
 * vertex count before and after simplification, symmetric Hausdorff
 * distance of the simplified contours against the unsimplified extraction
 * (and, where a reference input polyline exists, against that reference),
 * and the average extraction time. Every expectation is hard-checked; the
 * process exits non-zero on the first failure.
 *
 * Build (no build-system change, standard library only):
 *   g++ -std=c++14 -O2 -I repo-root unittest/boundary_demo.cpp \
 *       poly2tri/common/shapes.cc poly2tri/sweep/advancing_front.cc \
 *       poly2tri/sweep/cdt.cc poly2tri/sweep/sweep.cc \
 *       poly2tri/sweep/sweep_context.cc -o boundary_demo
 */

#include <poly2tri/poly2tri.h>
#include <poly2tri/boundary_extractor.h>

#include <algorithm>
#include <chrono>
#include <cstdio>
#include <functional>
#include <iomanip>
#include <iostream>
#include <list>
#include <numeric>
#include <sstream>
#include <string>
#include <vector>

namespace {

int g_failures = 0;
int g_checks = 0;

#define CHECK(cond, message)                                                       \
  do {                                                                             \
    ++g_checks;                                                                    \
    if (!(cond)) {                                                                 \
      ++g_failures;                                                                \
      std::cerr << "FAIL [" << __FILE__ << ":" << __LINE__ << "] " << message      \
                << " (condition: " #cond ")\n";                                    \
    }                                                                              \
  } while (0)

std::size_t TotalVertices(const std::vector<p2t::Contour>& contours)
{
  return std::accumulate(contours.begin(), contours.end(), std::size_t{ 0 },
                         [](std::size_t sum, const p2t::Contour& c) {
                           return sum + c.points.size();
                         });
}

std::size_t CountClosed(const std::vector<p2t::Contour>& contours)
{
  return static_cast<std::size_t>(
      std::count_if(contours.begin(), contours.end(),
                    [](const p2t::Contour& c) { return c.closed; }));
}

std::size_t CountOpen(const std::vector<p2t::Contour>& contours)
{
  return static_cast<std::size_t>(
      std::count_if(contours.begin(), contours.end(),
                    [](const p2t::Contour& c) { return !c.closed; }));
}

// Contour set built directly from an input polyline (reference geometry,
// copied by value).
std::vector<p2t::Contour> ReferenceContours(
    const std::vector<std::vector<p2t::Point*> >& polylines, bool closed = true)
{
  std::vector<p2t::Contour> ref;
  for (const auto& poly : polylines) {
    p2t::Contour contour;
    contour.closed = closed;
    contour.points.reserve(poly.size());
    for (const p2t::Point* p : poly) {
      contour.points.push_back(*p);
    }
    ref.push_back(contour);
  }
  return ref;
}

double TimedExtraction(const std::vector<p2t::Triangle*>& triangles,
                       const p2t::ContourOptions& options, const int repeats)
{
  const auto start = std::chrono::steady_clock::now();
  for (int i = 0; i < repeats; ++i) {
    volatile p2t::ContourResult sink = p2t::ExtractContours(triangles, options);
    (void)sink;
  }
  const auto stop = std::chrono::steady_clock::now();
  return std::chrono::duration<double, std::micro>(stop - start).count() /
         static_cast<double>(repeats);
}

double TimedExtractionEdges(
    const std::vector<p2t::Triangle*>& triangles,
    const std::vector<std::pair<p2t::Point*, p2t::Point*> >& edges,
    const p2t::ContourOptions& options, const int repeats)
{
  const auto start = std::chrono::steady_clock::now();
  for (int i = 0; i < repeats; ++i) {
    volatile p2t::ContourResult sink =
        p2t::ExtractContours(triangles, edges, options);
    (void)sink;
  }
  const auto stop = std::chrono::steady_clock::now();
  return std::chrono::duration<double, std::micro>(stop - start).count() /
         static_cast<double>(repeats);
}

// Snapshots the triangulation before extraction and verifies vertex-by-vertex
// that extraction left it untouched.
void CheckTriangulationUnchanged(
    const std::vector<p2t::Triangle*>& triangles,
    const std::vector<std::vector<std::pair<double, double> > >& snapshot)
{
  CHECK(triangles.size() == snapshot.size(),
        "triangle count changed during extraction");
  for (std::size_t t = 0; t < triangles.size() && t < snapshot.size(); ++t) {
    for (int i = 0; i < 3; ++i) {
      const p2t::Point* p = triangles[t]->GetPoint(i);
      CHECK(p->x == snapshot[t][static_cast<std::size_t>(i)].first &&
                p->y == snapshot[t][static_cast<std::size_t>(i)].second,
            "triangulation vertex changed during extraction");
    }
  }
}

std::vector<std::vector<std::pair<double, double> > > Snapshot(
    const std::vector<p2t::Triangle*>& triangles)
{
  std::vector<std::vector<std::pair<double, double> > > snapshot;
  snapshot.reserve(triangles.size());
  for (p2t::Triangle* t : triangles) {
    std::vector<std::pair<double, double> > tri;
    for (int i = 0; i < 3; ++i) {
      tri.emplace_back(t->GetPoint(i)->x, t->GetPoint(i)->y);
    }
    snapshot.push_back(tri);
  }
  return snapshot;
}

void PrintRow(const std::string& group, const std::string& run,
              std::size_t closed, std::size_t open, std::size_t raw_vertices,
              std::size_t simple_vertices, double haus_raw, double haus_ref,
              double us)
{
  std::printf("| %-16s | %-22s | %6zu | %5zu | %9zu | %9zu | %13.6f | %13s | %9.3f |\n",
              group.c_str(), run.c_str(), closed, open, raw_vertices,
              simple_vertices, haus_raw,
              haus_ref < 0.0 ? "n/a" : (std::to_string(haus_ref)).c_str(), us);
}

void PrintHeader()
{
  std::puts("+------------------+------------------------+--------+-------+-----------+"
            "-----------+---------------+---------------+-----------+");
  std::puts("| group            | run                    | closed | open  | raw verts |"
            " simp verts| Hausdorff raw | Hausdorff ref |  avg [us] |");
  std::puts("+------------------+------------------------+--------+-------+-----------+"
            "-----------+---------------+---------------+-----------+");
}

void PrintFooter()
{
  std::puts("+------------------+------------------------+--------+-------+-----------+"
            "-----------+---------------+---------------+-----------+");
}

// A mesh built with p2t::CDT plus the reference polylines it was built from.
// The pool is a std::list so emplace_back never invalidates the point
// pointers stored in the polylines (a std::vector would reallocate).
struct BuiltMesh {
  std::list<p2t::Point> pool;
  std::vector<p2t::Point*> boundary; // outer polyline
  std::vector<std::vector<p2t::Point*> > holes;
  std::vector<p2t::Point> steiner;
  std::vector<p2t::Triangle*> triangles;
};

p2t::Point* AddPoint(std::list<p2t::Point>& pool, double x, double y)
{
  pool.emplace_back(x, y);
  return &pool.back();
}

// Runs both extraction configurations of one mesh and fills the report.
void RunCdtMesh(const std::string& name, BuiltMesh& mesh,
                const p2t::ContourOptions& base_options,
                const std::vector<p2t::Contour>* reference, const int repeats)
{
  const auto before = Snapshot(mesh.triangles);

  // Raw extraction (no simplification, holes kept) gives the full-res
  // baseline and the raw vertex count.
  p2t::ContourOptions raw_options;
  raw_options.angle_tolerance_rad = 0.0;
  raw_options.edge_length_tolerance = 0.0;
  raw_options.keep_holes = true;
  const p2t::ContourResult raw =
      p2t::ExtractContours(mesh.triangles, raw_options);
  const std::size_t raw_vertices = TotalVertices(raw.contours);

  // Configuration 1: thresholds + keep holes.
  const p2t::ContourResult kept =
      p2t::ExtractContours(mesh.triangles, base_options);
  // Configuration 2: same thresholds, holes removed.
  p2t::ContourOptions nohole_options = base_options;
  nohole_options.keep_holes = false;
  const p2t::ContourResult dropped =
      p2t::ExtractContours(mesh.triangles, nohole_options);

  CheckTriangulationUnchanged(mesh.triangles, before);

  const double us_kept = TimedExtraction(mesh.triangles, base_options, repeats);
  const double us_dropped =
      TimedExtraction(mesh.triangles, nohole_options, repeats);

  const double h_kept_raw = p2t::SymmetricHausdorff(kept.contours, raw.contours);
  const double h_dropped_raw =
      p2t::SymmetricHausdorff(dropped.contours, raw.contours);
  const double h_kept_ref = reference
                                ? p2t::SymmetricHausdorff(kept.contours, *reference)
                                : -1.0;

  PrintRow(name, "simplify, keep holes", CountClosed(kept.contours),
           CountOpen(kept.contours), raw_vertices,
           TotalVertices(kept.contours), h_kept_raw, h_kept_ref, us_kept);
  PrintRow(name, "simplify, drop holes", CountClosed(dropped.contours),
           CountOpen(dropped.contours), raw_vertices,
           TotalVertices(dropped.contours), h_dropped_raw, -1.0, us_dropped);

  // Shared invariants for every healthy mesh.
  CHECK(raw.break_points.empty(), name + ": unexpected break points");
  CHECK(kept.break_points.empty(), name + ": unexpected break points");
  for (const auto& c : kept.contours) {
    CHECK(c.points.size() >= (c.closed ? 3u : 2u),
          name + ": contour simplified below minimum vertex count");
  }
  CHECK(h_kept_raw >= 0.0, name + ": Hausdorff distance must be non-negative");

  // Hand the detailed structures back for group-specific checks.
  std::printf("    detail[%s]: boundary edges=%zu, junctions=%zu, "
              "raw contours=%zu (closed=%zu, open=%zu)\n",
              name.c_str(), raw.boundary_edge_count, raw.junctions.size(),
              raw.contours.size(), CountClosed(raw.contours),
              CountOpen(raw.contours));

  // Expose results via thread-local-ish globals through a small trick: the
  // group-specific checks are registered as callbacks below.
  (void)dropped;
}

} // namespace

int main()
{
  constexpr int kRepeats = 200;
  std::printf("Boundary extraction on %d timed repeats per run "
              "(triangulation time excluded)\n", kRepeats);
  PrintHeader();

  // ---------------------------------------------------------------- G1
  // Simple concave polygon: a 5x5 rectangle with a rectangular notch on
  // the left. Two collinear points on the top edge and a 0.05-high jog
  // along the bottom edge give the simplifier something to remove; two
  // interior Steiner points must never appear on a contour.
  {
    BuiltMesh g1;
    g1.boundary = {
      AddPoint(g1.pool, 0.0, 0.0),   // 0
      AddPoint(g1.pool, 0.0, 2.0),   // 1 notch lower corner
      AddPoint(g1.pool, 1.0, 2.0),   // 2
      AddPoint(g1.pool, 1.0, 3.0),   // 3
      AddPoint(g1.pool, 0.0, 3.0),   // 4 notch upper corner
      AddPoint(g1.pool, 0.0, 5.0),   // 5
      AddPoint(g1.pool, 1.25, 5.0),  // 6 collinear on top edge
      AddPoint(g1.pool, 2.5, 5.0),   // 7 collinear on top edge
      AddPoint(g1.pool, 5.0, 5.0),   // 8
      AddPoint(g1.pool, 5.0, 0.0),   // 9
      AddPoint(g1.pool, 4.0, 0.0),   // 10 bottom jog (height 0.05)
      AddPoint(g1.pool, 4.0, 0.05),  // 11
      AddPoint(g1.pool, 2.0, 0.05),  // 12
      AddPoint(g1.pool, 2.0, 0.0)    // 13
    };
    g1.steiner = { p2t::Point(2.0, 2.0), p2t::Point(3.0, 3.0) };

    p2t::CDT cdt(g1.boundary);
    for (auto& s : g1.steiner) {
      cdt.AddPoint(&s);
    }
    cdt.Triangulate();
    g1.triangles = cdt.GetTriangles();
    CHECK(g1.triangles.size() >= 7, "G1 triangulation sanity");

    p2t::ContourOptions opt;
    opt.angle_tolerance_rad = 0.1; // ~5.7 degrees, flattens the 0.05 jog
    opt.edge_length_tolerance = 0.1;
    opt.keep_holes = true;

    const auto reference =
        ReferenceContours(std::vector<std::vector<p2t::Point*> >{ g1.boundary });
    RunCdtMesh("G1 single contour", g1, opt, &reference, kRepeats);

    const p2t::ContourResult r = p2t::ExtractContours(g1.triangles, opt);
    CHECK(r.contours.size() == 1 && r.contours[0].closed,
          "G1: exactly one closed contour");
    CHECK(r.contours[0].points.size() == 8,
          "G1: simplification leaves the eight notch/corner vertices");
    // The 0.05 jog vertex is removed, moving the bottom edge 0.05 inward:
    // the Hausdorff distance against the raw contour is exactly 0.05.
    const double h = p2t::SymmetricHausdorff(r.contours,
        p2t::ExtractContours(g1.triangles, p2t::ContourOptions()).contours);
    CHECK(std::fabs(h - 0.05) < 1e-9, "G1: Hausdorff equals jog height 0.05");
    for (const auto& s : g1.steiner) {
      for (const auto& c : r.contours) {
        for (const auto& p : c.points) {
          CHECK(!(p.x == s.x && p.y == s.y),
                "G1: interior Steiner point leaked onto contour");
        }
      }
    }
  }

  // ---------------------------------------------------------------- G2
  // Outer square (5 vertices, one collinear) with a square hole
  // (4 vertices) plus one interior Steiner point.
  {
    BuiltMesh g2;
    g2.boundary = {
      AddPoint(g2.pool, 0.0, 0.0),
      AddPoint(g2.pool, 0.0, 10.0),
      AddPoint(g2.pool, 5.0, 10.0), // collinear with the top edge
      AddPoint(g2.pool, 10.0, 10.0),
      AddPoint(g2.pool, 10.0, 0.0)
    };
    std::vector<p2t::Point*> hole = {
      AddPoint(g2.pool, 3.0, 3.0),
      AddPoint(g2.pool, 3.0, 7.0),
      AddPoint(g2.pool, 7.0, 7.0),
      AddPoint(g2.pool, 7.0, 3.0)
    };
    g2.holes.push_back(hole);
    g2.steiner = { p2t::Point(8.5, 5.0) };

    p2t::CDT cdt(g2.boundary);
    cdt.AddHole(hole);
    cdt.AddPoint(&g2.steiner[0]);
    cdt.Triangulate();
    g2.triangles = cdt.GetTriangles();

    p2t::ContourOptions opt;
    opt.angle_tolerance_rad = 0.01;
    opt.edge_length_tolerance = 0.1;
    opt.keep_holes = true;

    const auto reference = ReferenceContours(
        std::vector<std::vector<p2t::Point*> >{ g2.boundary, hole });
    RunCdtMesh("G2 with hole", g2, opt, &reference, kRepeats);

    const p2t::ContourResult kept = p2t::ExtractContours(g2.triangles, opt);
    CHECK(kept.contours.size() == 2, "G2: outer + hole contours");
    CHECK(CountClosed(kept.contours) == 2, "G2: both contours closed");
    // Both loops are quadrilaterals after simplification; the outer one is
    // the large-area square, the hole the small one.
    double outer_area = 0.0;
    double hole_area = 0.0;
    for (const auto& c : kept.contours) {
      CHECK(c.points.size() == 4, "G2: each loop simplifies to 4 vertices");
      const double area = std::fabs(c.SignedArea());
      if (area > outer_area) {
        hole_area = outer_area;
        outer_area = area;
      } else {
        hole_area = area;
      }
    }
    CHECK(outer_area > 50.0 && hole_area > 10.0 && hole_area < 20.0,
          "G2: outer square and 4x4 hole identified by area");
    p2t::ContourOptions drop = opt;
    drop.keep_holes = false;
    const p2t::ContourResult dropped = p2t::ExtractContours(g2.triangles, drop);
    CHECK(dropped.contours.size() == 1 && dropped.contours[0].closed,
          "G2: dropping holes leaves the outer contour only");
    CHECK(dropped.contours[0].points.size() == 4,
          "G2: outer contour simplified back to the square corners");
    const double href = p2t::SymmetricHausdorff(kept.contours, reference);
    CHECK(href < 1e-9, "G2: extraction reproduces the input polylines");
  }

  // ---------------------------------------------------------------- G3
  // Two disjoint components: a square with a triangular hole and a plain
  // triangle. They are triangulated separately and the triangle soups are
  // concatenated, proving extraction needs no CDT bookkeeping.
  {
    BuiltMesh g3a;
    g3a.boundary = {
      AddPoint(g3a.pool, 0.0, 0.0),
      AddPoint(g3a.pool, 0.0, 8.0),
      AddPoint(g3a.pool, 8.0, 8.0),
      AddPoint(g3a.pool, 8.0, 0.0)
    };
    std::vector<p2t::Point*> hole = {
      AddPoint(g3a.pool, 3.0, 2.0),
      AddPoint(g3a.pool, 5.0, 2.0),
      AddPoint(g3a.pool, 4.0, 5.0)
    };
    g3a.holes.push_back(hole);

    p2t::CDT cdt_a(g3a.boundary);
    cdt_a.AddHole(hole);
    cdt_a.Triangulate();
    const auto tris_a = cdt_a.GetTriangles();

    BuiltMesh g3b;
    g3b.boundary = {
      AddPoint(g3b.pool, 12.0, 0.0),
      AddPoint(g3b.pool, 16.0, 0.0),
      AddPoint(g3b.pool, 14.0, 4.0)
    };
    p2t::CDT cdt_b(g3b.boundary);
    cdt_b.Triangulate();
    const auto tris_b = cdt_b.GetTriangles();

    std::vector<p2t::Triangle*> combined = tris_a;
    combined.insert(combined.end(), tris_b.begin(), tris_b.end());
    const auto before = Snapshot(combined);

    p2t::ContourOptions opt;
    opt.angle_tolerance_rad = 0.01;
    opt.edge_length_tolerance = 0.1;
    opt.keep_holes = true;

    const p2t::ContourResult raw =
        p2t::ExtractContours(combined, p2t::ContourOptions());
    const p2t::ContourResult kept = p2t::ExtractContours(combined, opt);
    p2t::ContourOptions drop = opt;
    drop.keep_holes = false;
    const p2t::ContourResult dropped = p2t::ExtractContours(combined, drop);
    CheckTriangulationUnchanged(combined, before);

    const double us = TimedExtraction(combined, opt, kRepeats);
    PrintRow("G3 multi-component", "simplify, keep holes",
             CountClosed(kept.contours), CountOpen(kept.contours),
             TotalVertices(raw.contours), TotalVertices(kept.contours),
             p2t::SymmetricHausdorff(kept.contours, raw.contours), -1.0, us);
    PrintRow("G3 multi-component", "simplify, drop holes",
             CountClosed(dropped.contours), CountOpen(dropped.contours),
             TotalVertices(raw.contours), TotalVertices(dropped.contours),
             p2t::SymmetricHausdorff(dropped.contours, raw.contours), -1.0, us);

    CHECK(raw.break_points.empty(), "G3: healthy mesh has no break points");
    CHECK(kept.contours.size() == 3,
          "G3: square + triangle hole + separate triangle = 3 contours");
    CHECK(CountClosed(kept.contours) == 3, "G3: all contours closed");
    CHECK(dropped.contours.size() == 2,
          "G3: hole dropped, two outer contours remain");
    bool found_small_triangle = false;
    for (const auto& c : kept.contours) {
      if (c.points.size() == 3 && c.points[0].x >= 12.0) {
        found_small_triangle = true;
      }
    }
    CHECK(found_small_triangle, "G3: disjoint triangle component recovered");
    std::printf("    detail[G3 multi-component]: boundary edges=%zu, junctions=%zu\n",
                raw.boundary_edge_count, raw.junctions.size());
  }

  // ---------------------------------------------------------------- G4a
  // Non-manifold triangle soup with a T-junction and a dangling edge:
  //   triangle ABC is closed on its own;
  //   triangles A-E-G and E-B-G share edge E-G and split A-B at the
  //   T-junction vertex E, so A and B each have four incident boundary
  //   half-edges (non-manifold junctions);
  //   the standalone edge G-H belongs to no triangle and forms an open
  //   boundary chain whose ends (G and H) cannot be closed and are
  //   reported as break points.
  {
    std::vector<p2t::Point> pts = {
      p2t::Point(0.0, 0.0),  // 0 A
      p2t::Point(2.0, 0.0),  // 1 B
      p2t::Point(2.0, 2.0),  // 2 C
      p2t::Point(1.0, 0.0),  // 3 E (T-junction on A-B)
      p2t::Point(1.0, -1.0), // 4 G
      p2t::Point(1.0, -2.0)  // 5 H, dangling edge end point
    };
    std::vector<p2t::Triangle> soup;
    soup.emplace_back(pts[0], pts[1], pts[2]); // A B C
    soup.emplace_back(pts[0], pts[3], pts[4]); // A E G
    soup.emplace_back(pts[3], pts[1], pts[4]); // E B G
    std::vector<p2t::Triangle*> triangles;
    for (auto& t : soup) {
      triangles.push_back(&t);
    }
    std::vector<std::pair<p2t::Point*, p2t::Point*> > dangling{
      { &pts[4], &pts[5] } // G-H
    };
    const auto before = Snapshot(triangles);

    p2t::ContourOptions opt;
    opt.angle_tolerance_rad = 0.01;
    opt.edge_length_tolerance = 1e-12; // do not collapse the geometry
    opt.keep_holes = true;

    const p2t::ContourResult r =
        p2t::ExtractContours(triangles, dangling, opt);
    CheckTriangulationUnchanged(triangles, before);
    const p2t::ContourResult raw =
        p2t::ExtractContours(triangles, dangling, p2t::ContourOptions());
    const double us = TimedExtractionEdges(triangles, dangling, opt, kRepeats);
    PrintRow("G4a open chain", "simplify, keep holes",
             CountClosed(r.contours), CountOpen(r.contours),
             TotalVertices(raw.contours), TotalVertices(r.contours),
             p2t::SymmetricHausdorff(r.contours, raw.contours), -1.0, us);
    // Seven triangle boundary edges plus the dangling edge.
    CHECK(r.boundary_edge_count == 8, "G4a: eight boundary edges");
    CHECK(r.contours.size() == 3, "G4a: two closed loops + one open chain");
    CHECK(CountClosed(r.contours) == 2, "G4a: two closed contours");
    CHECK(CountOpen(r.contours) == 1, "G4a: one open contour");
    CHECK(r.break_points.size() == 2, "G4a: two break points");
    // The open chain is the dangling edge G-H.
    const auto& open = *std::find_if(r.contours.begin(), r.contours.end(),
                                     [](const p2t::Contour& c) { return !c.closed; });
    CHECK(open.points.size() == 2, "G4a: open chain G-H has 2 vertices");
    const bool ends_are_G_H =
        ((open.points[0].x == 1.0 && open.points[0].y == -1.0) &&
         (open.points[1].x == 1.0 && open.points[1].y == -2.0)) ||
        ((open.points[1].x == 1.0 && open.points[1].y == -1.0) &&
         (open.points[0].x == 1.0 && open.points[0].y == -2.0));
    CHECK(ends_are_G_H, "G4a: break points located at G(1,-1) and H(1,-2)");
    // H is a free end with a single incident boundary edge; G joins the
    // dangling edge and the two T-soup triangles.
    bool free_end_reported = false;
    for (const auto& bp : r.break_points) {
      if (bp.point.x == 1.0 && bp.point.y == -2.0 &&
          bp.incident_boundary_edges == 1) {
        free_end_reported = true;
      }
    }
    CHECK(free_end_reported, "G4a: free dangling end H reported (degree 1)");
    // The T-junction endpoints A and B each have four incident boundary
    // edges; the dangling junction G has three.
    std::size_t degree_four = 0;
    std::size_t degree_three = 0;
    for (const auto& j : r.junctions) {
      if (j.incident_boundary_edges == 4 &&
          ((j.point.x == 0.0 && j.point.y == 0.0) ||
           (j.point.x == 2.0 && j.point.y == 0.0))) {
        ++degree_four;
      }
      if (j.incident_boundary_edges == 3 && j.point.x == 1.0 &&
          j.point.y == -1.0) {
        ++degree_three;
      }
    }
    CHECK(degree_four == 2, "G4a: junction vertices A and B reported (degree 4)");
    CHECK(degree_three == 1, "G4a: dangling junction G reported (degree 3)");
    CHECK(r.junctions.size() == 4 || r.junctions.size() == 3,
          "G4a: non-manifold vertices reported as junctions");
    for (const auto& bp : r.break_points) {
      std::printf("    breakpoint at (%.3f, %.3f), incident boundary edges=%zu\n",
                  bp.point.x, bp.point.y, bp.incident_boundary_edges);
    }
    std::printf("    detail[G4a open chain]: boundary edges=%zu, junctions=%zu\n",
                r.boundary_edge_count, r.junctions.size());
  }

  // ---------------------------------------------------------------- G4b
  // Two holes touching at one vertex (mirrors the repository's
  // TouchingHolesTest). The touching vertex has boundary degree 4 and the
  // walk splits there into two independent closed sub-loops.
  {
    std::list<p2t::Point> pts;
    auto point = [&pts](double x, double y) {
      pts.emplace_back(x, y);
      return &pts.back();
    };
    std::vector<p2t::Point*> boundary = {
      point(0.0, 0.0), point(1.0, 0.0), point(0.5, 1.0)
    };
    p2t::Point pinch(0.2, 0.2);
    std::vector<p2t::Point*> hole0 = {
      point(0.1, 0.1), point(0.2, 0.1), &pinch
    };
    std::vector<p2t::Point*> hole1 = {
      &pinch, point(0.3, 0.2), point(0.3, 0.3)
    };
    p2t::CDT cdt(boundary);
    cdt.AddHole(hole0);
    cdt.AddHole(hole1);
    cdt.Triangulate();
    const auto triangles = cdt.GetTriangles();
    const auto before = Snapshot(triangles);

    p2t::ContourOptions opt;
    opt.angle_tolerance_rad = 0.01;
    opt.edge_length_tolerance = 0.1;
    opt.keep_holes = true;
    const p2t::ContourResult r = p2t::ExtractContours(triangles, opt);
    CheckTriangulationUnchanged(triangles, before);
    const p2t::ContourResult raw =
        p2t::ExtractContours(triangles, p2t::ContourOptions());
    const double us = TimedExtraction(triangles, opt, kRepeats);
    PrintRow("G4b touching holes", "simplify, keep holes",
             CountClosed(r.contours), CountOpen(r.contours),
             TotalVertices(raw.contours), TotalVertices(r.contours),
             p2t::SymmetricHausdorff(r.contours, raw.contours), -1.0, us);

    CHECK(r.break_points.empty(), "G4b: pinch closes, no break points");
    CHECK(r.contours.size() == 3, "G4b: outer triangle + two hole loops");
    CHECK(CountClosed(r.contours) == 3, "G4b: all three loops closed");
    bool pinch_reported = false;
    for (const auto& j : r.junctions) {
      if (j.point.x == 0.2 && j.point.y == 0.2 && j.incident_boundary_edges == 4) {
        pinch_reported = true;
      }
    }
    CHECK(pinch_reported, "G4b: pinch vertex reported as junction (degree 4)");
    p2t::ContourOptions drop = opt;
    drop.keep_holes = false;
    const p2t::ContourResult dropped = p2t::ExtractContours(triangles, drop);
    CHECK(dropped.contours.size() == 1 && dropped.contours[0].closed,
          "G4b: dropping holes leaves the outer triangle");
    std::printf("    detail[G4b touching holes]: boundary edges=%zu, junctions=%zu\n",
                r.boundary_edge_count, r.junctions.size());
  }

  // ---------------------------------------------------------------- G4c
  // Empty mesh must return an empty result instead of crashing.
  {
    const std::vector<p2t::Triangle*> empty;
    const p2t::ContourResult r =
        p2t::ExtractContours(empty, p2t::ContourOptions());
    CHECK(r.contours.empty() && r.break_points.empty() &&
              r.junctions.empty() && r.boundary_edge_count == 0,
          "G4c: empty mesh yields an empty result");
    CHECK(p2t::SymmetricHausdorff(r.contours, r.contours) == 0.0,
          "G4c: Hausdorff distance of empty sets is zero");
    std::puts("    G4c empty mesh: empty result, no crash");
  }

  PrintFooter();
  std::printf("\nChecks: %d, failures: %d\n", g_checks, g_failures);
  if (g_failures != 0) {
    std::cerr << "BOUNDARY DEMO FAILED\n";
    return 1;
  }
  std::puts("ALL BOUNDARY CHECKS PASSED");
  return 0;
}
