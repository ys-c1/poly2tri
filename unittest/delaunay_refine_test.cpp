/*
 * Acceptance tests for DelaunayRefine (Ruppert-style Delaunay refinement).
 *
 * Standalone: uses only the C++ standard library and poly2tri itself, no
 * third-party test framework. Build and run through run_refine_tests.sh.
 *
 * Four required inputs plus two robustness tests:
 *   1. long thin strip (acute angles)
 *   2. polygon with a hole
 *   3. polygon with many constrained edges (L-shape)
 *   4. already-good mesh (must produce zero insertions)
 *   5. point cap reached (must stop and report, never loop forever)
 *   6. inherently acute input (unimprovable corner; stop, report, intact)
 *
 * Ownership: each test keeps its original vertices in `orig*` vectors. The
 * DelaunayRefine allocates every edge-midpoint / interior Steiner point it
 * creates and deletes them in its destructor, so only the original vertices
 * are deleted by hand here.
 */

#include <poly2tri/sweep/delaunay_refine.h>

#include <cassert>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <string>
#include <vector>

using p2t::CDT;
using p2t::DelaunayRefine;
using p2t::Point;
using p2t::Triangle;

namespace {

/// Smallest interior angle (degrees) across all triangles.
double MinAngleDeg(std::vector<Triangle*> tris)
{
  double smallest = 180.0;
  for (Triangle* t : tris) {
    for (int i = 0; i < 3; ++i) {
      const Point& o = *t->GetPoint(i);
      const Point& p = *t->GetPoint((i + 1) % 3);
      const Point& q = *t->GetPoint((i + 2) % 3);
      const double v1x = p.x - o.x, v1y = p.y - o.y;
      const double v2x = q.x - o.x, v2y = q.y - o.y;
      const double n1 = std::sqrt(v1x * v1x + v1y * v1y);
      const double n2 = std::sqrt(v2x * v2x + v2y * v2y);
      double cosine = (v1x * v2x + v1y * v2y) / (n1 * n2);
      cosine = std::max(-1.0, std::min(1.0, cosine));
      smallest = std::min(smallest, std::acos(cosine) * 180.0 / M_PI);
    }
  }
  return smallest;
}

/// Unrefined CDT minimum angle, to show the "before" picture.
double RawMinAngleDeg(const std::vector<Point*>& boundary,
                      const std::vector<std::vector<Point*> >& holes)
{
  CDT cdt(boundary);
  for (const auto& h : holes) {
    cdt.AddHole(h);
  }
  cdt.Triangulate();
  return MinAngleDeg(cdt.GetTriangles());
}

bool OnSegment(const Point& p, const Point& a, const Point& b, double eps)
{
  const double cross = (b.x - a.x) * (p.y - a.y) - (b.y - a.y) * (p.x - a.x);
  if (std::fabs(cross) > eps * (1.0 + std::hypot(b.x - a.x, b.y - a.y))) {
    return false;
  }
  const double dot = (p.x - a.x) * (p.x - b.x) + (p.y - a.y) * (p.y - b.y);
  return dot <= eps;  // between a and b (inclusive)
}

/**
 * Constraint integrity, two independent checks:
 *  (a) every maintained constrained segment appears as an actual triangle
 *      edge in the final mesh (no segment was dropped);
 *  (b) no triangle edge crosses any constrained segment interior (no
 *      constrained edge was crossed).
 */
bool ConstraintsIntact(DelaunayRefine& refiner)
{
  auto& tris = refiner.GetTriangles();
  const auto segs = refiner.GetConstrainedSegments();
  const double eps = 1e-9;

  for (const auto& s : segs) {
    // (a) the segment endpoints must be joined by a triangle edge.
    bool found = false;
    for (Triangle* t : tris) {
      for (int i = 0; i < 3; ++i) {
        Point* p = t->GetPoint(i);
        Point* q = t->GetPoint((i + 1) % 3);
        if ((p == s.a && q == s.b) || (p == s.b && q == s.a)) {
          found = true;
          break;
        }
      }
      if (found) break;
    }
    if (!found) return false;

    // (b) no other triangle edge may properly cross the segment.
    for (Triangle* t : tris) {
      for (int i = 0; i < 3; ++i) {
        Point* e1 = t->GetPoint(i);
        Point* e2 = t->GetPoint((i + 1) % 3);
        // Edges sharing an endpoint with the constraint cannot cross it.
        if (e1 == s.a || e1 == s.b || e2 == s.a || e2 == s.b) continue;
        if (OnSegment(*e1, *s.a, *s.b, eps) || OnSegment(*e2, *s.a, *s.b, eps)) {
          return false;
        }
        const double d1 = (s.b->x - s.a->x) * (e1->y - s.a->y) -
                          (s.b->y - s.a->y) * (e1->x - s.a->x);
        const double d2 = (s.b->x - s.a->x) * (e2->y - s.a->y) -
                          (s.b->y - s.a->y) * (e2->x - s.a->x);
        const double d3 = (e2->x - e1->x) * (s.a->y - e1->y) -
                          (e2->y - e1->y) * (s.a->x - e1->x);
        const double d4 = (e2->x - e1->x) * (s.b->y - e1->y) -
                          (e2->y - e1->y) * (s.b->x - e1->x);
        if (((d1 > eps && d2 < -eps) || (d1 < -eps && d2 > eps)) &&
            ((d3 > eps && d4 < -eps) || (d3 < -eps && d4 > eps))) {
          return false;  // proper intersection
        }
      }
    }
  }
  return true;
}

void Report(const std::string& name, double raw_min,
            const DelaunayRefine::Result& r, bool intact)
{
  std::cout << "[" << name << "]\n"
            << "  raw min angle (before refine) : " << raw_min << " deg\n"
            << "  refined min angle              : " << r.min_angle_deg << " deg\n"
            << "  triangle count                 : " << r.triangle_count << "\n"
            << "  Steiner points inserted        : " << r.inserted_points << "\n"
            << "  target reached                 : "
            << (r.target_reached ? "yes" : "no") << "\n"
            << "  cap hit                        : "
            << (r.cap_hit ? "yes" : "no") << "\n"
            << "  constraint integrity           : "
            << (intact ? "OK (none crossed)" : "FAIL") << "\n"
            << std::endl;
}

void DeleteAll(const std::vector<Point*>& pts)
{
  for (Point* p : pts) delete p;
}

}  // namespace

int main()
{
  // ------------------------------------------------------------------
  // Test 1: long thin strip with acute angles. Target 20 deg; acceptance
  // requires refined min angle >= 18 deg and zero crossed constraints.
  // ------------------------------------------------------------------
  {
    std::vector<Point*> orig{
      new Point(0, 0), new Point(20, 0), new Point(20, 0.8), new Point(0, 0.8)
    };
    std::vector<Point*> boundary = orig;
    std::vector<std::vector<Point*> > holes;
    const double raw = RawMinAngleDeg(boundary, holes);

    DelaunayRefine::Result r;
    bool intact = false;
    {
      DelaunayRefine refiner(boundary, holes, 20.0, 1000);
      r = refiner.Refine();
      intact = ConstraintsIntact(refiner);
      Report("1 long-thin strip", raw, r, intact);
    }  // refiner frees the Steiner/midpoint points it allocated

    assert(r.inserted_points > 0);                 // it really refined
    assert(r.min_angle_deg >= 18.0);               // required acceptance bar
    assert(intact);                                // no constraint crossed
    assert(r.target_reached);
    DeleteAll(orig);
  }

  // ------------------------------------------------------------------
  // Test 2: polygon with a hole.
  // ------------------------------------------------------------------
  {
    std::vector<Point*> origB{
      new Point(0, 0), new Point(10, 0), new Point(10, 10), new Point(0, 10)
    };
    std::vector<Point*> origH{
      new Point(4, 4), new Point(6, 4), new Point(6, 6), new Point(4, 6)
    };
    std::vector<Point*> boundary = origB;
    std::vector<Point*> hole = origH;
    std::vector<std::vector<Point*> > holes{ hole };
    const double raw = RawMinAngleDeg(boundary, holes);

    DelaunayRefine::Result r;
    bool intact = false;
    {
      DelaunayRefine refiner(boundary, holes, 20.0, 1000);
      r = refiner.Refine();
      intact = ConstraintsIntact(refiner);
      Report("2 polygon with hole", raw, r, intact);
    }

    assert(r.target_reached);
    assert(r.min_angle_deg >= 18.0);
    assert(intact);
    DeleteAll(origB);
    DeleteAll(origH);
  }

  // ------------------------------------------------------------------
  // Test 3: many constrained edges (L-shaped concave polygon, 6 edges),
  // made skinny so it starts below target and genuinely needs refinement.
  // ------------------------------------------------------------------
  {
    std::vector<Point*> orig{
      new Point(0, 0), new Point(20, 0), new Point(20, 0.8),
      new Point(0.8, 0.8), new Point(0.8, 20), new Point(0, 20)
    };
    std::vector<Point*> boundary = orig;
    std::vector<std::vector<Point*> > holes;
    const double raw = RawMinAngleDeg(boundary, holes);

    DelaunayRefine::Result r;
    bool intact = false;
    std::size_t final_segments = 0;
    {
      DelaunayRefine refiner(boundary, holes, 20.0, 1000);
      r = refiner.Refine();
      intact = ConstraintsIntact(refiner);
      final_segments = refiner.GetConstrainedSegments().size();
      Report("3 L-shape multi constraints", raw, r, intact);
    }

    assert(r.target_reached);
    assert(r.min_angle_deg >= 18.0);
    assert(intact);
    assert(final_segments >= orig.size());  // splits only add segments
    DeleteAll(orig);
  }

  // ------------------------------------------------------------------
  // Test 4: already-good mesh -> zero insertions.
  // A rhombus of two equilateral triangles has a minimum angle of 60 deg.
  // ------------------------------------------------------------------
  {
    const double h = std::sqrt(3.0) / 2.0;
    std::vector<Point*> orig{
      new Point(0, 0), new Point(1, 0), new Point(1.5, h), new Point(0.5, h)
    };
    std::vector<Point*> boundary = orig;
    std::vector<std::vector<Point*> > holes;
    const double raw = RawMinAngleDeg(boundary, holes);

    DelaunayRefine::Result r;
    bool intact = false;
    {
      DelaunayRefine refiner(boundary, holes, 20.0, 1000);
      r = refiner.Refine();
      intact = ConstraintsIntact(refiner);
      Report("4 already-good mesh", raw, r, intact);
    }

    assert(r.inserted_points == 0);                // mandatory: zero insert
    assert(r.target_reached);
    assert(!r.cap_hit);
    assert(intact);
    DeleteAll(orig);
  }

  // ------------------------------------------------------------------
  // Test 5: point cap must stop the loop and be reported, never hang.
  // ------------------------------------------------------------------
  {
    std::vector<Point*> orig{
      new Point(0, 0), new Point(20, 0), new Point(20, 0.8), new Point(0, 0.8)
    };
    std::vector<Point*> boundary = orig;
    std::vector<std::vector<Point*> > holes;
    const double raw = RawMinAngleDeg(boundary, holes);

    const std::size_t cap = 5;
    DelaunayRefine::Result r;
    bool intact = false;
    {
      DelaunayRefine refiner(boundary, holes, 20.0, cap);
      r = refiner.Refine();
      intact = ConstraintsIntact(refiner);
      Report("5 point-cap termination", raw, r, intact);
    }

    assert(r.cap_hit);                             // stopped because of cap
    assert(r.inserted_points == cap);              // inserted exactly cap
    assert(!r.target_reached);
    assert(intact);
    DeleteAll(orig);
  }

  // ------------------------------------------------------------------
  // Test 6: inherently acute input (apex ~8 deg) with a 20 deg target.
  // Ruppert refinement cannot remove an input angle smaller than the
  // target, so the routine must stop without an infinite loop or a crash,
  // keep constraints intact, and report the target as not reached.
  // ------------------------------------------------------------------
  {
    std::vector<Point*> orig{
      new Point(0, 0), new Point(10, 0), new Point(5, 0.7)
    };
    std::vector<Point*> boundary = orig;
    std::vector<std::vector<Point*> > holes;
    const double raw = RawMinAngleDeg(boundary, holes);

    DelaunayRefine::Result r;
    bool intact = false;
    {
      DelaunayRefine refiner(boundary, holes, 20.0, 1000);
      r = refiner.Refine();
      intact = ConstraintsIntact(refiner);
      Report("6 acute-input (unimprovable corner)", raw, r, intact);
    }

    assert(!r.cap_hit);                            // stopped by skip, not cap
    assert(!r.target_reached);                     // honestly reported
    assert(r.min_angle_deg < 20.0);                // small input angle stays
    assert(intact);                                // constraints never crossed
    DeleteAll(orig);
  }

  std::cout << "ALL DELAUNAY REFINEMENT TESTS PASSED\n";
  return 0;
}
