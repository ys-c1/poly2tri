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
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#pragma once

#include "dll_symbol.h"
#include "shapes.h"

#include <cstddef>
#include <vector>

namespace p2t {

/// One self-intersection of a polyline / polygon outline.
struct P2T_DLL_SYMBOL SelfIntersection {
  double x, y;        ///< Intersection point. For `overlap`, the lexicographically
                      ///< smallest point of the shared interval.
  std::size_t edge_a; ///< Index of the first edge (edge i connects vertex i to i+1,
                      ///< with wrap-around for closed outlines). edge_a < edge_b.
  std::size_t edge_b;
  bool overlap;       ///< true: the two edges are collinear and share an interval
                      ///< of positive length; false: they meet at a single point
                      ///< (proper crossing or endpoint/interior touching).
};

/// Sweep-line (Bentley-Ottmann) self-intersection detection.
///
/// Semantics (identical for both functions):
///  - Edge i connects vertex i to vertex (i+1) mod n (closed) or i+1 (open).
///  - Zero-length edges are ignored.
///  - Edges adjacent in the outline (indices differing by 1, or {0, m-1} for a
///    closed outline) never count as intersecting each other.
///  - Any other pair of edges that shares at least one point is reported once:
///    proper crossings, endpoint-to-endpoint touches, endpoint-on-interior
///    touches (T-junctions) and collinear overlaps.
///
/// `out` receives up to `max_report` records, sorted by (x, y, edge_a, edge_b)
/// before return. Returns true iff at least one self-intersection exists.
///
/// Truncation rule: the sweep collects intersections in discovery order
/// (increasing sweep position) and stops as soon as `max_report` records have
/// been found; `truncated` (if non-null) is then set to true. Because discovery
/// order is not exactly the final sort order, a truncated result is guaranteed
/// to be a subset of the true intersections, but not necessarily the 100
/// lexicographically smallest ones. `truncated` is false iff the result is
/// complete.
///
/// Complexity: O((n + k) log n) with k capped by max_report (the sweep aborts
/// early), versus O(n^2) for the brute-force reference below.
bool P2T_DLL_SYMBOL FindSelfIntersections(const std::vector<Point>& polyline,
                                          bool closed,
                                          std::size_t max_report,
                                          std::vector<SelfIntersection>& out,
                                          bool* truncated = nullptr);

/// O(n^2) brute-force reference implementation with identical semantics.
/// Kept as the validation oracle and as the performance baseline.
bool P2T_DLL_SYMBOL FindSelfIntersectionsBruteForce(const std::vector<Point>& polyline,
                                                    bool closed,
                                                    std::size_t max_report,
                                                    std::vector<SelfIntersection>& out,
                                                    bool* truncated = nullptr);

} // namespace p2t
