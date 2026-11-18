/*
 * Poly2Tri Copyright (c) 2009-2024, Poly2Tri Contributors
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
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

/**
 * @file delaunay_refine.h
 *
 * @brief Ruppert-style Delaunay refinement on top of the poly2tri CDT.
 *
 * The sweep-line CDT in this library produces a *legal* constrained
 * triangulation but gives no quality guarantees: long thin input polygons
 * can produce triangles with very small minimum angles. This file adds a
 * quality refinement pass (Ruppert's algorithm, "Triangle: Engineering a 2D
 * Quality Mesh Generator and Delaunay Triangulator", Shewchuk 1996):
 *
 *   1. Any triangle whose minimum angle is below the configured target is
 *      "bad". Its circumcenter is the candidate Steiner point.
 *   2. A circumcenter inside the diametral lens of a constrained segment
 *      would encroach on (cross) that segment. The segment is split at its
 *      midpoint instead and the midpoint becomes a Steiner vertex on the
 *      constraint, so constrained edges are never crossed.
 *   3. A circumcenter outside the domain (possible next to an acute input
 *      corner) cannot be inserted; the nearest constrained segment is split
 *      so the loop keeps making progress.
 *   4. Otherwise the circumcenter is inserted strictly inside the domain.
 *   5. The constrained Delaunay triangulation is rebuilt after every
 *      insertion, so every intermediate mesh is again a legal CDT.
 *
 * The constraint polylines are the single source of truth: edge splits are
 * recorded by inserting the midpoint into the corresponding boundary/hole
 * polyline, and constrained segments are derived from those polylines on
 * demand rather than maintained in a parallel container.
 *
 * Small input angles: plain Ruppert refinement does not converge when the
 * input itself contains angles below the target. To keep the routine
 * robust (and to avoid ever producing coincident vertices that the
 * sweep-line triangulator cannot handle) two safeguards are added, matching
 * the practice in Shewchuk's Triangle:
 *   - a minimum edge-length floor: a constrained segment shorter than the
 *     floor is never split; and
 *   - a bad triangle whose only remedy is forbidden (an out-of-domain
 *     center or a sub-floor segment split) is skipped rather than retried.
 * Such triangles stay in the mesh and are reported as unimproved instead of
 * causing an infinite loop or a crash. The point cap is a second, hard
 * termination guarantee.
 *
 * Termination: the loop stops as soon as no remaining triangle can be
 * improved, or when the number of inserted Steiner points reaches the
 * configured upper bound. DelaunayRefine::Result distinguishes "all
 * triangles meet the target" from "stopped because of the cap /
 * unimprovable input angles".
 *
 * Ownership: the caller owns the points it passes in through the polylines.
 * Every midpoint/interior Steiner point created here is allocated by this
 * class and deleted in its destructor; because edge-split midpoints are
 * inserted into the polylines by pointer, the caller must not delete them.
 *
 * Header-only and built entirely on the existing public CDT interface, so
 * no build scripts or existing public signatures change.
 */

#pragma once

#include "cdt.h"
#include "../common/shapes.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <set>
#include <vector>

// C99 removes M_PI from math.h
#ifndef M_PI
#define M_PI 3.14159265358979323846264338327
#endif

namespace p2t {

class DelaunayRefine {
public:

  /// A straight constrained segment between two mesh vertices.
  struct Segment {
    Point* a;
    Point* b;
  };

  /// Outcome / statistics of a refinement run.
  struct Result {
    /// True when every interior triangle meets the minimum-angle target.
    /// False if the cap was hit or an input angle made some triangle
    /// unimprovable (see min_angle_deg for the achieved value).
    bool target_reached;
    /// True when refinement stopped because the Steiner-point cap was hit.
    bool cap_hit;
    /// Number of Steiner points inserted by the refinement pass.
    std::size_t inserted_points;
    /// Minimum interior angle (degrees) found in the final mesh.
    double min_angle_deg;
    /// Number of triangles in the final mesh.
    std::size_t triangle_count;
  };

  /**
   * Constructor.
   *
   * @param boundary outer boundary polygon (closed without repeating the
   *                 first point). Edge-split midpoints are inserted here.
   * @param holes    interior hole polylines, may be empty. Edge-split
   *                 midpoints are inserted into the affected hole.
   * @param min_angle_deg desired minimum angle in degrees.
   * @param max_steiner_points hard upper bound on inserted Steiner points.
   */
  DelaunayRefine(std::vector<Point*>& boundary,
                 std::vector<std::vector<Point*> >& holes,
                 double min_angle_deg = 20.0,
                 std::size_t max_steiner_points = 1000)
    : boundary_(boundary),
      holes_(holes),
      min_angle_rad_(min_angle_deg * M_PI / 180.0),
      max_steiner_points_(max_steiner_points),
      inserted_(0),
      cap_hit_(false),
      min_edge_length_(0.0)
  {
  }

  /// Run the refinement. Call once per instance.
  Result Refine()
  {
    // Floor is the shortest input edge. No constrained segment is ever
    // split below this length, which prevents unbounded refinement toward
    // acute input corners.
    min_edge_length_ = ShortestInputEdgeLength();

    // First legal CDT of the un-refined input.
    Rebuild();

    Result result;
    result.inserted_points = 0;
    result.min_angle_deg = MinAngleDeg();
    result.triangle_count = triangles_.size();

    // Already good: zero insertions, mesh returned untouched.
    if (FindBadTriangle() == nullptr) {
      result.target_reached = true;
      result.cap_hit = false;
      return result;
    }

    while (true) {
      Triangle* bad = FindBadTriangle();
      if (bad == nullptr) {
        break;  // every remaining triangle meets the angle target
      }
      if (inserted_ >= max_steiner_points_) {
        cap_hit_ = true;
        break;
      }

      Point center = Circumcenter(*bad);

      // Constrained segments may not be crossed: if the circumcenter
      // encroaches on one, split that segment instead of inserting it.
      Segment encroached;
      bool progressed = false;
      if (FindEncroachedSegment(center, encroached)) {
        if (SegmentLength(encroached) >= min_edge_length_) {
          SplitConstrainedSegment(encroached);
          progressed = true;
        }
      } else if (!PointInDomain(center)) {
        // Center outside the polygon (acute input corner): split the
        // nearest constrained segment if it is still long enough.
        Segment nearest;
        if (FindNearestSegment(center, nearest) &&
            SegmentLength(nearest) >= min_edge_length_) {
          SplitConstrainedSegment(nearest);
          progressed = true;
        }
      } else {
        InsertInteriorPoint(new Point(center.x, center.y));
        progressed = true;
      }

      if (!progressed) {
        // This bad triangle cannot be improved without creating a
        // sub-floor edge (small input angle). Skip it so the loop cannot
        // spin forever; it remains in the mesh and in the angle report.
        skipped_.insert(TriangleKey(*bad));
        continue;
      }

      Rebuild();
    }

    // Ground truth: did any final triangle end up below the target? The
    // skip set is ignored here on purpose.
    result.target_reached = !cap_hit_ && !AnyTriangleBelowTarget();
    result.cap_hit = cap_hit_;
    result.inserted_points = inserted_;
    result.min_angle_deg = MinAngleDeg();
    result.triangle_count = triangles_.size();
    return result;
  }

  /// Final interior triangles, valid for the lifetime of this object.
  std::vector<Triangle*>& GetTriangles()
  {
    return triangles_;
  }

  /**
   * The constrained segments of the current refinement state, derived from
   * the (possibly edge-split) boundary and hole polylines. Useful for
   * constraint-integrity checks.
   */
  std::vector<Segment> GetConstrainedSegments() const
  {
    std::vector<Segment> out;
    CollectLoopSegments(boundary_, out);
    for (const auto& hole : holes_) {
      CollectLoopSegments(hole, out);
    }
    return out;
  }

  /// Steiner points created by the refiner (owned and deleted here).
  const std::vector<Point*>& GetSteinerPoints() const
  {
    return owned_points_;
  }

  ~DelaunayRefine()
  {
    delete cdt_;
    for (Point* p : owned_points_) {
      delete p;
    }
  }

private:

  // Small angular slack (radians) so floating-point noise cannot make a
  // triangle that exactly meets the target look bad.
  static double AngleTolRad() { return 1e-9 * M_PI / 180.0; }

  std::vector<Point*>& boundary_;
  std::vector<std::vector<Point*> >& holes_;
  double min_angle_rad_;
  std::size_t max_steiner_points_;
  std::size_t inserted_;
  bool cap_hit_;
  /// Shortest input edge; segments shorter than this are never split.
  double min_edge_length_;

  CDT* cdt_ = nullptr;
  std::vector<Triangle*> triangles_;
  /// Steiner points strictly inside the domain (edge midpoints are conveyed
  /// through the polylines instead).
  std::vector<Point*> interior_steiner_;
  /// Every point allocated by this class (deleted in destructor).
  std::vector<Point*> owned_points_;
  /// Vertex-pointer triples of triangles given up as unimprovable.
  std::set<std::set<const Point*> > skipped_;

  static void CollectLoopSegments(const std::vector<Point*>& loop,
                                  std::vector<Segment>& out)
  {
    for (std::size_t i = 0; i < loop.size(); ++i) {
      out.push_back(Segment{ loop[i], loop[(i + 1) % loop.size()] });
    }
  }

  double ShortestInputEdgeLength() const
  {
    std::vector<Segment> all;
    CollectLoopSegments(boundary_, all);
    for (const auto& hole : holes_) {
      CollectLoopSegments(hole, all);
    }
    double shortest = 0.0;
    for (std::size_t i = 0; i < all.size(); ++i) {
      const double len = SegmentLength(all[i]);
      if (i == 0 || len < shortest) {
        shortest = len;
      }
    }
    return shortest;
  }

  static double SegmentLength(const Segment& s)
  {
    return std::hypot(s.b->x - s.a->x, s.b->y - s.a->y);
  }

  /// Build a fresh CDT from the current polylines and interior Steiner
  /// points. The previous CDT (and all its Triangle objects) is destroyed.
  void Rebuild()
  {
    delete cdt_;
    // poly2tri's Edge objects register themselves in each Point's
    // edge_list, but SweepContext deletes its Edges without clearing those
    // lists. Clear the dangling Edge pointers here; the new CDT rebuilds
    // them from the constraint polylines.
    for (Point* p : boundary_) {
      p->edge_list.clear();
    }
    for (const auto& hole : holes_) {
      for (Point* p : hole) {
        p->edge_list.clear();
      }
    }
    for (Point* p : interior_steiner_) {
      p->edge_list.clear();
    }

    cdt_ = new CDT(boundary_);
    for (const auto& hole : holes_) {
      cdt_->AddHole(hole);
    }
    for (Point* p : interior_steiner_) {
      cdt_->AddPoint(p);
    }
    cdt_->Triangulate();
    triangles_ = cdt_->GetTriangles();
    // Topology changed: re-evaluate every triangle rather than honouring
    // skip decisions made for the previous mesh. Termination still holds
    // because every rebuild consumed one insertion, bounded by the point cap.
    skipped_.clear();
  }

  void InsertInteriorPoint(Point* p)
  {
    owned_points_.push_back(p);
    interior_steiner_.push_back(p);
    ++inserted_;
  }

  /**
   * Insert @p mid between the consecutive vertices @p a and @p b (unordered)
   * of @p loop. Returns false if the pair is not part of this loop.
   */
  static bool InsertVertexOnLoop(std::vector<Point*>& loop, Point* a,
                                 Point* b, Point* mid)
  {
    for (std::size_t i = 0; i < loop.size(); ++i) {
      Point* p = loop[i];
      Point* q = loop[(i + 1) % loop.size()];
      if ((p == a && q == b) || (p == b && q == a)) {
        loop.insert(loop.begin() + static_cast<std::ptrdiff_t>(i) + 1, mid);
        return true;
      }
    }
    return false;
  }

  /// Split a constrained segment at its midpoint (edge split).
  void SplitConstrainedSegment(const Segment& seg)
  {
    Point* mid = new Point(0.5 * (seg.a->x + seg.b->x),
                           0.5 * (seg.a->y + seg.b->y));
    owned_points_.push_back(mid);
    ++inserted_;

    if (!InsertVertexOnLoop(boundary_, seg.a, seg.b, mid)) {
      for (auto& hole : holes_) {
        if (InsertVertexOnLoop(hole, seg.a, seg.b, mid)) {
          break;
        }
      }
    }
  }

  /**
   * Diametral-lens (encroachment) test: @p p is inside the open diametral
   * circle of segment s when it is closer to the segment midpoint than half
   * the segment length.
   */
  static bool InDiametralLens(const Point& p, const Segment& s)
  {
    const double mx = 0.5 * (s.a->x + s.b->x);
    const double my = 0.5 * (s.a->y + s.b->y);
    const double dx = s.b->x - s.a->x;
    const double dy = s.b->y - s.a->y;
    const double lens_r2 = 0.25 * (dx * dx + dy * dy);
    const double px = p.x - mx;
    const double py = p.y - my;
    return px * px + py * py < lens_r2;
  }

  /// All current constrained segments (boundary loops plus hole loops).
  std::vector<Segment> AllSegments() const
  {
    std::vector<Segment> all;
    CollectLoopSegments(boundary_, all);
    for (const auto& hole : holes_) {
      CollectLoopSegments(hole, all);
    }
    return all;
  }

  /// Find a constrained segment encroached by @p p, if any.
  bool FindEncroachedSegment(const Point& p, Segment& out) const
  {
    const std::vector<Segment> all = AllSegments();
    for (const Segment& s : all) {
      if (InDiametralLens(p, s)) {
        out = s;
        return true;
      }
    }
    return false;
  }

  /// Find the constrained segment whose midpoint is closest to @p p.
  bool FindNearestSegment(const Point& p, Segment& out) const
  {
    const std::vector<Segment> all = AllSegments();
    if (all.empty()) {
      return false;
    }
    std::size_t best = 0;
    double best_d2 = 0.0;
    for (std::size_t i = 0; i < all.size(); ++i) {
      const Segment& s = all[i];
      const double mx = 0.5 * (s.a->x + s.b->x);
      const double my = 0.5 * (s.a->y + s.b->y);
      const double d2 = (p.x - mx) * (p.x - mx) + (p.y - my) * (p.y - my);
      if (i == 0 || d2 < best_d2) {
        best_d2 = d2;
        best = i;
      }
    }
    out = all[best];
    return true;
  }

  /// Circumcenter of a triangle.
  static Point Circumcenter(Triangle& t)
  {
    const Point& a = *t.GetPoint(0);
    const Point& b = *t.GetPoint(1);
    const Point& c = *t.GetPoint(2);
    const double d = 2.0 * (a.x * (b.y - c.y) + b.x * (c.y - a.y) +
                            c.x * (a.y - b.y));
    const double a2 = a.x * a.x + a.y * a.y;
    const double b2 = b.x * b.x + b.y * b.y;
    const double c2 = c.x * c.x + c.y * c.y;
    const double ux = (a2 * (b.y - c.y) + b2 * (c.y - a.y) +
                       c2 * (a.y - b.y)) / d;
    const double uy = (a2 * (c.x - b.x) + b2 * (a.x - c.x) +
                       c2 * (b.x - a.x)) / d;
    return Point(ux, uy);
  }

  /// Minimum interior angle (radians) of a triangle.
  static double MinAngleRad(Triangle& t)
  {
    double smallest = M_PI;
    for (int i = 0; i < 3; ++i) {
      const Point& o = *t.GetPoint(i);
      const Point& p = *t.GetPoint((i + 1) % 3);
      const Point& q = *t.GetPoint((i + 2) % 3);
      const double v1x = p.x - o.x, v1y = p.y - o.y;
      const double v2x = q.x - o.x, v2y = q.y - o.y;
      const double n1 = std::sqrt(v1x * v1x + v1y * v1y);
      const double n2 = std::sqrt(v2x * v2x + v2y * v2y);
      double cosine = (v1x * v2x + v1y * v2y) / (n1 * n2);
      cosine = std::max(-1.0, std::min(1.0, cosine));
      smallest = std::min(smallest, std::acos(cosine));
    }
    return smallest;
  }

  double MinAngleDeg()
  {
    double smallest = M_PI;
    for (Triangle* t : triangles_) {
      smallest = std::min(smallest, MinAngleRad(*t));
    }
    return smallest * 180.0 / M_PI;
  }

  /// Stable identity of a triangle from its (unordered) vertex set.
  static std::set<const Point*> TriangleKey(Triangle& t)
  {
    return std::set<const Point*>{ t.GetPoint(0), t.GetPoint(1),
                                   t.GetPoint(2) };
  }

  /// Ground-truth check: any triangle (including skipped ones) below target?
  bool AnyTriangleBelowTarget()
  {
    for (Triangle* t : triangles_) {
      if (MinAngleRad(*t) + AngleTolRad() < min_angle_rad_) {
        return true;
      }
    }
    return false;
  }

  /// Worst triangle (smallest minimum angle) that fails the angle target
  /// and has not been given up as unimprovable.
  Triangle* FindBadTriangle()
  {
    Triangle* worst = nullptr;
    double worst_angle = min_angle_rad_;
    for (Triangle* t : triangles_) {
      if (skipped_.count(TriangleKey(*t)) != 0) {
        continue;
      }
      const double a = MinAngleRad(*t);
      if (a + AngleTolRad() < min_angle_rad_ &&
          (worst == nullptr || a < worst_angle)) {
        worst = t;
        worst_angle = a;
      }
    }
    return worst;
  }

  /**
   * Point-in-domain test for a (possibly holey) region defined by the
   * boundary and hole loops. Standard ray casting; boundary points count as
   * inside.
   */
  bool PointInDomain(const Point& p) const
  {
    if (!PointInLoop(p, boundary_)) {
      return false;
    }
    for (const auto& hole : holes_) {
      if (PointInLoop(p, hole)) {
        return false;
      }
    }
    return true;
  }

  static bool PointInLoop(const Point& p, const std::vector<Point*>& loop)
  {
    bool inside = false;
    for (std::size_t i = 0, j = loop.size() - 1; i < loop.size(); j = i++) {
      const Point& vi = *loop[i];
      const Point& vj = *loop[j];
      if (((vi.y > p.y) != (vj.y > p.y)) &&
          (p.x < (vj.x - vi.x) * (p.y - vi.y) / (vj.y - vi.y) + vi.x)) {
        inside = !inside;
      }
    }
    return inside;
  }
};

} // namespace p2t
