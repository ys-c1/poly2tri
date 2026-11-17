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

/**
 * @file boundary_extractor.h
 *
 * Reverse-engineer boundary contours from a triangle mesh produced by
 * p2t::CDT (or any soup of p2t::Triangle objects). The extraction is
 * read-only: triangulation data is never modified, and every point of the
 * result is copied by value so the result outlives the mesh.
 *
 * Three stages:
 *
 * 1. Mark boundary edges. Every directed triangle edge is canonicalised and
 *    counted. An edge used by exactly one triangle is a boundary edge. This
 *    is purely combinatorial (pointer identity), so it also works on a raw
 *    triangle soup, including non-manifold input.
 *
 * 2. Chain the boundary edges into contours. Each boundary edge is turned
 *    into a directed half-edge using the orientation of its incident
 *    triangle (the triangle interior lies to the left of the half-edge).
 *    Walking the half-edges, at every vertex the outgoing half-edge with
 *    the smallest left turn relative to the incoming one is selected,
 *    which keeps the filled domain on the left. This is purely geometric
 *    and gives correct results at non-manifold junctions (e.g. two holes
 *    touching at one vertex), where the walk splits into independent
 *    closed loops. Standalone edges that belong to no triangle start open
 *    chains; if a walk cannot be closed the open chain is kept and both
 *    of its ends are reported as break points, together with the number
 *    of incident boundary edges. Vertices whose boundary degree differs
 *    from two are reported as junctions.
 *
 * 3. Simplify each contour with two configurable thresholds:
 *    - angle_tolerance_rad: a vertex whose turning angle is smaller than
 *      this value is removable (collinear / nearly-collinear points);
 *    - edge_length_tolerance: a vertex adjacent to an edge shorter than
 *      this value is removable.
 *    Removal iterates until no vertex is removable. Closed contours keep
 *    at least three vertices, open contours at least two, and open contour
 *    end points are never removed.
 *
 * Hole classification is orientation independent (winding-number parity).
 * With ContourOptions::keep_holes == false every closed contour nested
 * inside another one is discarded.
 *
 * The file is header-only and depends solely on the C++14 standard library.
 */

#pragma once

#include "common/shapes.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <map>
#include <utility>
#include <vector>

namespace p2t {

/// Configuration of contour extraction and simplification.
struct ContourOptions {
  /// Vertices with a turning angle below this value (radians) may be removed.
  double angle_tolerance_rad = 0.0;
  /// Vertices adjacent to an edge shorter than this value may be removed.
  double edge_length_tolerance = 0.0;
  /// Keep closed contours that are holes (nested inside another contour).
  bool keep_holes = true;
};

/// A polygonal contour. Closed contours store each vertex once (no repeated
/// closing vertex).
struct Contour {
  bool closed = true;
  std::vector<Point> points;

  /// Signed area (shoelace formula). Positive for counter-clockwise loops.
  /// Zero for open contours treated as a closed fan.
  double SignedArea() const
  {
    double area = 0.0;
    for (std::size_t i = 0; i < points.size(); ++i) {
      const Point& a = points[i];
      const Point& b = points[(i + 1) % points.size()];
      area += a.x * b.y - b.x * a.y;
    }
    return 0.5 * area;
  }
};

/// End point of a boundary chain that could not be closed.
struct BreakPoint {
  Point point;
  /// Number of boundary edges incident to the vertex in the full mesh.
  std::size_t incident_boundary_edges = 0;
};

/// Vertex with more than two incident boundary edges (non-manifold junction,
/// e.g. two holes touching at a point).
struct Junction {
  Point point;
  std::size_t incident_boundary_edges = 0;
};

struct ContourResult {
  /// Number of distinct edges used by exactly one triangle.
  std::size_t boundary_edge_count = 0;
  /// All extracted contours (holes filtered according to the options).
  std::vector<Contour> contours;
  /// Ends of chains that could not be closed. Empty on a healthy mesh.
  std::vector<BreakPoint> break_points;
  /// Vertices with a boundary degree different from two.
  std::vector<Junction> junctions;
};

namespace boundary_detail {

inline std::pair<Point*, Point*> EdgeKey(Point* a, Point* b)
{
  return (a < b) ? std::make_pair(a, b) : std::make_pair(b, a);
}

inline double DistancePointSegment(const Point& p, const Point& a, const Point& b)
{
  const double vx = b.x - a.x;
  const double vy = b.y - a.y;
  const double wx = p.x - a.x;
  const double wy = p.y - a.y;
  const double len2 = vx * vx + vy * vy;
  double t = 0.0;
  if (len2 > 0.0) {
    t = (wx * vx + wy * vy) / len2;
    if (t > 1.0) {
      t = 1.0;
    } else if (t < 0.0) {
      t = 0.0;
    }
  }
  const double dx = wx - t * vx;
  const double dy = wy - t * vy;
  return std::sqrt(dx * dx + dy * dy);
}

/// Winding number of a closed contour around a point. Non-zero means the
/// point is inside (orientation independent).
inline int WindingNumber(const std::vector<Point>& poly, double qx, double qy)
{
  int wn = 0;
  const std::size_t n = poly.size();
  for (std::size_t i = 0; i < n; ++i) {
    const Point& a = poly[i];
    const Point& b = poly[(i + 1) % n];
    if (a.y <= qy) {
      if (b.y > qy && (b.x - a.x) * (qy - a.y) - (b.y - a.y) * (qx - a.x) > 0.0) {
        ++wn;
      }
    } else {
      if (b.y <= qy && (b.x - a.x) * (qy - a.y) - (b.y - a.y) * (qx - a.x) < 0.0) {
        --wn;
      }
    }
  }
  return wn;
}

/// Distance based "on boundary" test with a tolerance scaled to the
/// contour's coordinate magnitude.
inline bool OnBoundary(const std::vector<Point>& poly, double qx, double qy, double eps)
{
  const std::size_t n = poly.size();
  const Point q(qx, qy);
  for (std::size_t i = 0; i < n; ++i) {
    if (DistancePointSegment(q, poly[i], poly[(i + 1) % n]) <= eps) {
      return true;
    }
  }
  return false;
}

inline bool ContourContains(const std::vector<Point>& poly, const Point& q)
{
  double scale = 1.0;
  for (const Point& p : poly) {
    scale = std::max(scale, std::max(std::fabs(p.x), std::fabs(p.y)));
  }
  const double eps = 1e-10 * scale;
  return OnBoundary(poly, q.x, q.y, eps) || WindingNumber(poly, q.x, q.y) != 0;
}

/// Turning angle at cur for the walk prev -> cur -> next, in [0, pi].
inline double TurnAngle(const Point& prev, const Point& cur, const Point& next)
{
  const double ux = cur.x - prev.x;
  const double uy = cur.y - prev.y;
  const double vx = next.x - cur.x;
  const double vy = next.y - cur.y;
  const double lu = std::sqrt(ux * ux + uy * uy);
  const double lv = std::sqrt(vx * vx + vy * vy);
  double cos_a = 0.0;
  if (lu > 0.0 && lv > 0.0) {
    cos_a = (ux * vx + uy * vy) / (lu * lv);
  }
  if (cos_a < -1.0) {
    cos_a = -1.0;
  } else if (cos_a > 1.0) {
    cos_a = 1.0;
  }
  return std::acos(cos_a);
}

/// Signed area of a triangle; positive when (a, b, c) is counter-clockwise.
inline double SignedArea3(const Point& a, const Point& b, const Point& c)
{
  return 0.5 * ((b.x - a.x) * (c.y - a.y) - (c.x - a.x) * (b.y - a.y));
}

/// Polar angle of the vector (from -> to) in (-pi, pi].
inline double PolarAngle(const Point& from, const Point& to)
{
  return std::atan2(to.y - from.y, to.x - from.x);
}

/// Left-turn angle when going from direction a_in->cur to cur->cand,
/// measured in [0, 2*pi). The smallest value selects the sharpest left
/// turn, which is the next boundary half-edge keeping the triangle on the
/// left side of the walk.
inline double LeftTurnAngle(const Point& a_in, const Point& cur, const Point& cand)
{
  static const double kTwoPi = 2.0 * std::acos(-1.0);
  double angle = PolarAngle(cur, cand) - PolarAngle(a_in, cur);
  if (angle < 0.0) {
    angle += kTwoPi;
  }
  return angle;
}

inline Contour MakeContour(const std::vector<Point*>& chain, bool closed)
{
  Contour contour;
  contour.closed = closed;
  contour.points.reserve(chain.size());
  for (const Point* p : chain) {
    contour.points.push_back(*p);
  }
  return contour;
}

inline void SimplifyContour(Contour& contour, double angle_tol, double length_tol)
{
  const std::size_t min_vertices = contour.closed ? 3 : 2;
  if (contour.points.size() < min_vertices) {
    return;
  }
  bool removed_any = true;
  while (removed_any && contour.points.size() > min_vertices) {
    removed_any = false;
    const std::size_t n = contour.points.size();
    for (std::size_t i = 0; i < n && contour.points.size() > min_vertices;) {
      if (!contour.closed && (i == 0 || i == n - 1)) {
        ++i;
        continue;
      }
      const Point& prev = contour.points[(i + n - 1) % n];
      const Point& cur = contour.points[i];
      const Point& next = contour.points[(i + 1) % n];
      const double l1 = std::sqrt((cur.x - prev.x) * (cur.x - prev.x) +
                                  (cur.y - prev.y) * (cur.y - prev.y));
      const double l2 = std::sqrt((next.x - cur.x) * (next.x - cur.x) +
                                  (next.y - cur.y) * (next.y - cur.y));
      const bool short_edge = (std::min(l1, l2) < length_tol);
      const bool flat = (TurnAngle(prev, cur, next) < angle_tol);
      if (flat || short_edge) {
        contour.points.erase(contour.points.begin() +
                             static_cast<std::ptrdiff_t>(i));
        removed_any = true;
        break; // restart: neighbourhood of the neighbours changed
      }
      ++i;
    }
  }
}

} // namespace boundary_detail

/**
 * Extract boundary contours from a triangle mesh plus a set of standalone
 * edges. A standalone edge is not part of any triangle (for example a
 * dangling constraint edge); it is always treated as a boundary edge and
 * allows open chains whose ends are reported as break points.
 *
 * @param triangles       Mesh triangles; read-only, never modified.
 * @param standalone_edges Extra edges given as unordered endpoint pairs.
 * @param options         Simplification thresholds and hole filtering.
 */
inline ContourResult
ExtractContours(const std::vector<Triangle*>& triangles,
                const std::vector<std::pair<Point*, Point*> >& standalone_edges,
                const ContourOptions& options = ContourOptions())
{
  using boundary_detail::EdgeKey;
  using boundary_detail::LeftTurnAngle;
  using boundary_detail::MakeContour;
  using boundary_detail::PolarAngle;
  using boundary_detail::SignedArea3;
  using boundary_detail::SimplifyContour;

  ContourResult result;

  // Stage 1: mark boundary edges (used by exactly one triangle) and orient
  // them as half-edges. For a counter-clockwise triangle the triangle lies
  // to the left of each of its directed edges; for a clockwise triangle the
  // direction is reversed. Degenerate (zero-area) triangles are skipped.
  std::map<std::pair<Point*, Point*>, int> edge_use;
  std::map<std::pair<Point*, Point*>, Point*> oriented; // key -> directed tail->head
  for (Triangle* triangle : triangles) {
    Point* p[3] = { triangle->GetPoint(0), triangle->GetPoint(1),
                    triangle->GetPoint(2) };
    if (SignedArea3(*p[0], *p[1], *p[2]) == 0.0) {
      continue;
    }
    const bool ccw = SignedArea3(*p[0], *p[1], *p[2]) > 0.0;
    for (int i = 0; i < 3; ++i) {
      Point* tail = p[i];
      Point* head = p[(i + 1) % 3];
      if (!ccw) {
        std::swap(tail, head);
      }
      const auto key = EdgeKey(tail, head);
      ++edge_use[key];
      oriented[key] = tail; // direction tail -> head (the key itself is sorted)
    }
  }

  // Standalone edges are boundary edges by definition. Their orientation is
  // chosen arbitrarily but deterministically (from the sorted key).
  for (const auto& edge : standalone_edges) {
    const auto key = EdgeKey(edge.first, edge.second);
    edge_use[key] += 1;
    if (oriented.find(key) == oriented.end()) {
      oriented[key] = key.first;
    }
  }

  // Directed boundary half-edges: tail -> head, with the domain on the left.
  // outgoing[vertex] lists every boundary half-edge starting at the vertex.
  std::map<Point*, std::vector<Point*> > outgoing;
  std::size_t boundary_edges = 0;
  for (const auto& entry : edge_use) {
    if (entry.second != 1) {
      continue;
    }
    ++boundary_edges;
    Point* a = entry.first.first;
    Point* b = entry.first.second;
    if (oriented[entry.first] == a) {
      outgoing[a].push_back(b);
    } else {
      outgoing[b].push_back(a);
    }
  }
  result.boundary_edge_count = boundary_edges;

  // Undirected adjacency, needed for junction/break-point bookkeeping.
  std::map<Point*, std::vector<Point*> > adj;
  for (const auto& entry : edge_use) {
    if (entry.second == 1) {
      adj[entry.first.first].push_back(entry.first.second);
      adj[entry.first.second].push_back(entry.first.first);
    }
  }

  // Stage 2: walk the directed half-edges. At every vertex the outgoing
  // half-edge with the smallest left turn relative to the incoming one is
  // taken. This keeps the triangle domain on the left and gives correct
  // results also at non-manifold junctions. Half-edges are removed when
  // consumed; vertices with no outgoing half-edge terminate an open chain.
  std::vector<Contour> closed_loops;
  std::vector<Contour> open_chains;
  while (true) {
    Point* start = nullptr;
    for (auto& entry : outgoing) {
      if (!entry.second.empty()) {
        start = entry.first;
        break;
      }
    }
    if (start == nullptr) {
      break;
    }

    std::vector<Point*> chain;
    chain.push_back(start);
    Point* cur = start;
    Point* prev = nullptr;
    bool closed = false;

    while (true) {
      auto& candidates = outgoing[cur];
      if (candidates.empty()) {
        break; // open chain end: no boundary half-edge leaves this vertex
      }
      Point* next;
      if (prev == nullptr) {
        // First step: deterministic choice by polar angle.
        next = *std::min_element(candidates.begin(), candidates.end(),
                                 [&](Point* x, Point* y) {
                                   return PolarAngle(*cur, *x) <
                                          PolarAngle(*cur, *y);
                                 });
      } else {
        // Sharpest left turn keeps the filled domain on the left.
        next = *std::min_element(candidates.begin(), candidates.end(),
                                 [&](Point* x, Point* y) {
                                   return LeftTurnAngle(*prev, *cur, *x) <
                                          LeftTurnAngle(*prev, *cur, *y);
                                 });
      }
      candidates.erase(std::find(candidates.begin(), candidates.end(), next));

      if (next == start) {
        closed = true; // returned to the origin: closed loop
        break;
      }
      const auto hit = std::find(chain.begin(), chain.end(), next);
      if (hit != chain.end()) {
        // Pinch point: a sub-loop closes at an already visited vertex.
        const std::size_t k = static_cast<std::size_t>(hit - chain.begin());
        std::vector<Point*> sub(chain.begin() + static_cast<std::ptrdiff_t>(k),
                                chain.end());
        closed_loops.push_back(MakeContour(sub, true));
        chain.resize(k + 1);
        prev = (chain.size() >= 2) ? chain[chain.size() - 2] : nullptr;
        cur = chain.back();
        continue;
      }
      chain.push_back(next);
      prev = cur;
      cur = next;
    }

    if (closed) {
      if (chain.size() >= 3) {
        closed_loops.push_back(MakeContour(chain, true));
      }
    } else {
      open_chains.push_back(MakeContour(chain, false));
    }
  }

  // Non-manifold vertices: boundary degree different from two.
  for (const auto& entry : adj) {
    if (entry.second.size() != 2) {
      Junction junction;
      junction.point = *entry.first;
      junction.incident_boundary_edges = entry.second.size();
      result.junctions.push_back(junction);
    }
  }
  std::sort(result.junctions.begin(), result.junctions.end(),
            [](const Junction& a, const Junction& b) {
              return (a.point.x != b.point.x) ? (a.point.x < b.point.x)
                                              : (a.point.y < b.point.y);
            });

  // Break points: the two ends of every open chain.
  for (const Contour& chain : open_chains) {
    for (const Point* end : { &chain.points.front(), &chain.points.back() }) {
      BreakPoint bp;
      bp.point = *end;
      bp.incident_boundary_edges = 0;
      for (const auto& entry : adj) {
        if (entry.first->x == end->x && entry.first->y == end->y) {
          bp.incident_boundary_edges = entry.second.size();
          break;
        }
      }
      result.break_points.push_back(bp);
    }
  }

  // Hole classification on the closed loops (winding-number parity, so the
  // input orientation of holes does not matter). Probe points are edge
  // midpoints rather than vertices: at a non-manifold pinch point several
  // loops share a vertex, which would make a vertex-based containment test
  // ambiguous.
  std::vector<bool> is_hole(closed_loops.size(), false);
  for (std::size_t i = 0; i < closed_loops.size(); ++i) {
    const std::vector<Point>& poly_i = closed_loops[i].points;
    const std::size_t n = poly_i.size();
    std::size_t enclosing = 0;
    for (std::size_t k = 0; k < n; ++k) {
      const Point& a = poly_i[k];
      const Point& b = poly_i[(k + 1) % n];
      const Point mid(0.5 * (a.x + b.x), 0.5 * (a.y + b.y));
      std::size_t containing = 0;
      for (std::size_t j = 0; j < closed_loops.size(); ++j) {
        if (i != j && boundary_detail::ContourContains(closed_loops[j].points, mid)) {
          ++containing;
        }
      }
      enclosing = std::max(enclosing, containing);
    }
    is_hole[i] = (enclosing % 2u == 1u);
  }

  result.contours = open_chains; // open chains are always kept
  for (std::size_t i = 0; i < closed_loops.size(); ++i) {
    if (options.keep_holes || !is_hole[i]) {
      result.contours.push_back(closed_loops[i]);
    }
  }

  // Stage 3: per-contour simplification.
  for (Contour& contour : result.contours) {
    SimplifyContour(contour, options.angle_tolerance_rad,
                    options.edge_length_tolerance);
  }

  return result;
}

/**
 * Extract boundary contours from a triangle mesh.
 *
 * @param triangles Mesh triangles; read-only, never modified.
 * @param options   Simplification thresholds and hole filtering.
 */
inline ContourResult ExtractContours(const std::vector<Triangle*>& triangles,
                                     const ContourOptions& options = ContourOptions())
{
  return ExtractContours(triangles, std::vector<std::pair<Point*, Point*> >(),
                         options);
}

/**
 * Symmetric Hausdorff distance between two contour sets. Every vertex of one
 * set is compared against every segment of the other set and vice versa;
 * the maximum of both directed distances is returned.
 */
inline double SymmetricHausdorff(const std::vector<Contour>& a,
                                 const std::vector<Contour>& b)
{
  using boundary_detail::DistancePointSegment;

  auto segments_of = [](const std::vector<Contour>& set,
                        std::vector<std::pair<Point, Point> >& segments) {
    for (const Contour& contour : set) {
      const std::size_t n = contour.points.size();
      if (n < 2) {
        continue;
      }
      const std::size_t last = contour.closed ? n : n - 1;
      for (std::size_t i = 0; i < last; ++i) {
        segments.emplace_back(contour.points[i], contour.points[(i + 1) % n]);
      }
    }
  };

  std::vector<std::pair<Point, Point> > seg_a;
  std::vector<std::pair<Point, Point> > seg_b;
  segments_of(a, seg_a);
  segments_of(b, seg_b);

  double distance = 0.0;
  for (int pass = 0; pass < 2; ++pass) {
    const std::vector<Contour>& vertices_set = (pass == 0) ? a : b;
    const std::vector<std::pair<Point, Point> >& target = (pass == 0) ? seg_b : seg_a;
    if (target.empty()) {
      continue;
    }
    for (const Contour& contour : vertices_set) {
      for (const Point& p : contour.points) {
        double nearest = 0.0;
        for (std::size_t i = 0; i < target.size(); ++i) {
          const double d = DistancePointSegment(p, target[i].first, target[i].second);
          if (i == 0 || d < nearest) {
            nearest = d;
          }
        }
        distance = std::max(distance, nearest);
      }
    }
  }
  return distance;
}

} // namespace p2t
