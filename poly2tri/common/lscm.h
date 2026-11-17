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
 *   this list of conditions and the following disclaimer.
 * * Neither the name of Poly2Tri nor the names of its contributors may be
 *   used to endorse or promote products derived from this software without specific
 *   prior written permission.
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

// Approximate LSCM (Least Squares Conformal Maps) UV parameterization for
// planar triangle meshes, e.g. the output of p2t::CDT.
//
// Method:
//   1. Build a sparse linear system from cotangent edge weights (discrete
//      conformal / harmonic energy, the intrinsic part of LSCM).
//   2. Solve it twice (once for U, once for V) with boundary vertices
//      pinned to a convex chord-length parameterization (bounding rectangle
//      for the outer loop, circles for holes).
//   3. Optionally iterate: re-weight each triangle's contribution to reduce
//      area distortion, keeping the boundary soft-constrained towards its
//      convex target with a configurable weight.
//
// The input mesh is never modified; UVs are returned as separate arrays
// indexed like the input point array.
//
// Error handling: any singular / non-positive-definite / non-converging
// system (isolated vertices, degenerate triangles, disconnected components
// without boundary, ...) throws std::runtime_error instead of returning
// invalid UVs.

#include "shapes.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <map>
#include <stdexcept>
#include <unordered_map>
#include <utility>
#include <vector>

namespace p2t {

struct LscmOptions {
  /// Number of area-distortion minimization iterations after the initial
  /// conformal solve. 0 disables the iteration entirely.
  int iterations = 5;
  /// Weight of the soft boundary constraint during the iterations
  /// (boundary vertices are pulled towards their convex target UVs with
  /// this weight). Values <= 0 make the boundary hard-fixed instead.
  double boundary_weight = 1.0e3;
  /// Relative residual tolerance of the linear solver.
  double solver_tolerance = 1.0e-10;
  /// Iteration cap of the linear solver; exceeding it throws.
  int max_solver_iterations = 100000;
};

struct LscmResult {
  /// UV coordinates, same order and size as the input point array.
  std::vector<double> u;
  std::vector<double> v;
  /// Number of triangles with negative orientation in UV space.
  int flipped_triangles = 0;
  /// Area distortion: RMS over triangles of log(A_uv / A_3d) after
  /// normalizing the total UV area to the total 3D area (0 = perfectly
  /// area preserving).
  double area_distortion = 0.0;
};

namespace detail {

struct LscmCsrMatrix {
  int n = 0;
  std::vector<int> row; // size n + 1
  std::vector<int> col;
  std::vector<double> val;
  std::vector<double> diag; // size n
};

inline LscmCsrMatrix LscmBuildCsr(const std::vector<std::map<int, double>>& rows)
{
  LscmCsrMatrix m;
  m.n = static_cast<int>(rows.size());
  m.row.assign(static_cast<size_t>(m.n) + 1, 0);
  m.diag.assign(m.n, 0.0);
  for (int i = 0; i < m.n; ++i) {
    m.row[static_cast<size_t>(i) + 1] = m.row[i] + static_cast<int>(rows[i].size());
  }
  m.col.resize(m.row.back());
  m.val.resize(m.row.back());
  for (int i = 0; i < m.n; ++i) {
    int k = m.row[i];
    for (const auto& kv : rows[i]) {
      m.col[k] = kv.first;
      m.val[k] = kv.second;
      if (kv.first == i) {
        m.diag[i] = kv.second;
      }
      ++k;
    }
  }
  return m;
}

inline void LscmMatVec(const LscmCsrMatrix& A, const std::vector<double>& x,
                       std::vector<double>& y)
{
  std::fill(y.begin(), y.end(), 0.0);
  for (int i = 0; i < A.n; ++i) {
    double s = 0.0;
    for (int k = A.row[i]; k < A.row[i + 1]; ++k) {
      s += A.val[k] * x[A.col[k]];
    }
    y[i] = s;
  }
}

/// Preconditioned conjugate gradient for symmetric positive definite systems.
/// Throws std::runtime_error on breakdown (singular / indefinite matrix) or
/// when the residual target is not met within max_iter iterations.
inline void LscmSolveSpd(const LscmCsrMatrix& A, const std::vector<double>& b,
                         std::vector<double>& x, double tol, int max_iter)
{
  const int n = A.n;
  x.assign(n, 0.0);
  if (n == 0) {
    return;
  }
  for (int i = 0; i < n; ++i) {
    if (!(A.diag[i] > 0.0) || !std::isfinite(A.diag[i])) {
      throw std::runtime_error(
          "LscmParameterize: singular system (non-positive diagonal entry)");
    }
  }
  double bnorm = 0.0;
  for (double bi : b) {
    bnorm = std::max(bnorm, std::fabs(bi));
  }
  const double target = tol * std::max(1.0, bnorm);

  std::vector<double> r(b), z(n), p(n), Ap(n), minv(n);
  for (int i = 0; i < n; ++i) {
    minv[i] = 1.0 / A.diag[i];
    z[i] = r[i] * minv[i];
    p[i] = z[i];
  }
  double rz = 0.0;
  for (int i = 0; i < n; ++i) {
    rz += r[i] * z[i];
  }
  double rnorm = 0.0;
  for (double ri : r) {
    rnorm = std::max(rnorm, std::fabs(ri));
  }
  int iter = 0;
  while (rnorm > target && iter < max_iter) {
    LscmMatVec(A, p, Ap);
    double pAp = 0.0;
    for (int i = 0; i < n; ++i) {
      pAp += p[i] * Ap[i];
    }
    if (!(pAp > 0.0) || !std::isfinite(pAp)) {
      throw std::runtime_error(
          "LscmParameterize: singular system (CG breakdown, matrix not "
          "positive definite)");
    }
    const double alpha = rz / pAp;
    double next_rz = 0.0;
    rnorm = 0.0;
    for (int i = 0; i < n; ++i) {
      x[i] += alpha * p[i];
      r[i] -= alpha * Ap[i];
      z[i] = r[i] * minv[i];
      next_rz += r[i] * z[i];
      rnorm = std::max(rnorm, std::fabs(r[i]));
    }
    const double beta = next_rz / rz;
    for (int i = 0; i < n; ++i) {
      p[i] = z[i] + beta * p[i];
    }
    rz = next_rz;
    ++iter;
  }
  // Verify the residual explicitly instead of trusting the recurrence.
  std::vector<double> check(n);
  LscmMatVec(A, x, check);
  double res = 0.0;
  for (int i = 0; i < n; ++i) {
    res = std::max(res, std::fabs(check[i] - b[i]));
  }
  if (!(res <= 100.0 * target) || !std::isfinite(res)) {
    throw std::runtime_error(
        "LscmParameterize: singular system (linear solver did not converge)");
  }
}

} // namespace detail

/// Parameterize a planar triangle mesh into the UV plane.
///
/// @param points    Mesh vertices. Every triangle vertex must be present in
///                  this array and every entry must be referenced by at
///                  least one triangle (otherwise the system is singular
///                  and an exception is thrown).
/// @param triangles Mesh triangles (e.g. from p2t::CDT::GetTriangles()).
/// @param options   Iteration count, boundary weight and solver controls.
/// @return          UVs plus quality metrics. The input mesh is not modified.
/// @throws std::invalid_argument on structurally invalid input and
///         std::runtime_error when the linear system is singular.
inline LscmResult LscmParameterize(const std::vector<Point*>& points,
                                   const std::vector<Triangle*>& triangles,
                                   const LscmOptions& options = LscmOptions())
{
  constexpr double kPi = 3.14159265358979323846;
  const int n = static_cast<int>(points.size());
  if (n < 3) {
    throw std::invalid_argument("LscmParameterize: need at least 3 points");
  }
  if (triangles.empty()) {
    throw std::invalid_argument("LscmParameterize: no triangles");
  }

  std::unordered_map<const Point*, int> index;
  index.reserve(points.size());
  for (int i = 0; i < n; ++i) {
    if (points[i] == nullptr || !index.emplace(points[i], i).second) {
      throw std::invalid_argument(
          "LscmParameterize: null or duplicate point in input");
    }
  }

  // Index triangles and validate against degenerate (zero-area) input.
  std::vector<std::array<int, 3>> tris(triangles.size());
  std::vector<double> area3d(triangles.size());
  std::vector<char> used(n, 0);
  double total_area3d = 0.0;
  for (size_t t = 0; t < triangles.size(); ++t) {
    if (triangles[t] == nullptr) {
      throw std::invalid_argument("LscmParameterize: null triangle in input");
    }
    for (int k = 0; k < 3; ++k) {
      const auto it = index.find(triangles[t]->GetPoint(k));
      if (it == index.end()) {
        throw std::invalid_argument(
            "LscmParameterize: triangle references a point not in the point "
            "array");
      }
      tris[t][k] = it->second;
      used[it->second] = 1;
    }
    const Point& a = *points[tris[t][0]];
    const Point& b = *points[tris[t][1]];
    const Point& c = *points[tris[t][2]];
    const double cross = (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x);
    const double scale = (b - a).Length() * (c - a).Length();
    if (!(std::fabs(cross) > 1e-12 * scale) || !std::isfinite(cross)) {
      throw std::runtime_error(
          "LscmParameterize: singular system (degenerate triangle)");
    }
    // Orient CCW so that 3D areas are positive.
    if (cross < 0.0) {
      std::swap(tris[t][1], tris[t][2]);
    }
    area3d[t] = 0.5 * std::fabs(cross);
    total_area3d += area3d[t];
  }
  for (int i = 0; i < n; ++i) {
    if (!used[i]) {
      throw std::runtime_error(
          "LscmParameterize: singular system (isolated vertex)");
    }
  }

  // Cotangent edge weights and edge incidence counts.
  std::map<std::pair<int, int>, double> edge_weight;
  std::map<std::pair<int, int>, int> edge_count;
  std::vector<std::vector<int>> adjacency(n);
  for (size_t t = 0; t < tris.size(); ++t) {
    for (int k = 0; k < 3; ++k) {
      const int i = tris[t][k];
      const int j = tris[t][(k + 1) % 3];
      const int l = tris[t][(k + 2) % 3];
      // Angle at vertex l, opposite edge (i, j).
      const Point a = *points[i] - *points[l];
      const Point b = *points[j] - *points[l];
      const double cross = a.x * b.y - a.y * b.x;
      const double cot = Dot(a, b) / std::fabs(cross);
      const auto key = std::minmax(i, j);
      edge_weight[key] += 0.5 * cot; // base weight, re-weighted later
      if (edge_count[key]++ == 0) {
        adjacency[i].push_back(j);
        adjacency[j].push_back(i);
      }
    }
  }

  // Boundary edges belong to exactly one triangle.
  std::map<int, std::vector<int>> boundary_nbr;
  for (const auto& kv : edge_count) {
    if (kv.second == 1) {
      boundary_nbr[kv.first.first].push_back(kv.first.second);
      boundary_nbr[kv.first.second].push_back(kv.first.first);
    }
  }
  if (boundary_nbr.size() < 3) {
    throw std::runtime_error(
        "LscmParameterize: singular system (no boundary loop)");
  }
  std::vector<char> is_boundary(n, 0);
  for (const auto& kv : boundary_nbr) {
    if (kv.second.size() != 2) {
      throw std::runtime_error(
          "LscmParameterize: singular system (non-manifold boundary)");
    }
    is_boundary[kv.first] = 1;
  }

  // The graph must be connected and every component must touch the
  // boundary, otherwise the Dirichlet system has a zero eigenvalue.
  {
    std::vector<char> seen(n, 0);
    std::vector<int> stack;
    for (const auto& kv : boundary_nbr) {
      if (!seen[kv.first]) {
        seen[kv.first] = 1;
        stack.push_back(kv.first);
      }
    }
    while (!stack.empty()) {
      const int v = stack.back();
      stack.pop_back();
      for (const int w : adjacency[v]) {
        if (!seen[w]) {
          seen[w] = 1;
          stack.push_back(w);
        }
      }
    }
    for (int i = 0; i < n; ++i) {
      if (!seen[i]) {
        throw std::runtime_error(
            "LscmParameterize: singular system (mesh component without "
            "boundary)");
      }
    }
  }

  // Boundary loops, ordered.
  std::vector<std::vector<int>> loops;
  {
    std::map<std::pair<int, int>, bool> visited_edge;
    for (const auto& kv : boundary_nbr) {
      const int start = kv.first;
      for (const int next : kv.second) {
        if (visited_edge.count(std::minmax(start, next))) {
          continue;
        }
        std::vector<int> loop;
        int prev = start;
        int cur = next;
        loop.push_back(start);
        while (cur != start) {
          visited_edge[std::minmax(prev, cur)] = true;
          loop.push_back(cur);
          const auto& nbrs = boundary_nbr[cur];
          const int nxt = (nbrs[0] != prev) ? nbrs[0] : nbrs[1];
          prev = cur;
          cur = nxt;
        }
        visited_edge[std::minmax(prev, cur)] = true;
        loops.push_back(std::move(loop));
      }
    }
  }

  // Normalization frame of the 3D mesh.
  double minx = std::numeric_limits<double>::max();
  double miny = std::numeric_limits<double>::max();
  double maxx = std::numeric_limits<double>::lowest();
  double maxy = std::numeric_limits<double>::lowest();
  for (const Point* p : points) {
    minx = std::min(minx, p->x);
    miny = std::min(miny, p->y);
    maxx = std::max(maxx, p->x);
    maxy = std::max(maxy, p->y);
  }
  const double norm = std::max(maxx - minx, maxy - miny);
  if (!(norm > 0.0) || !std::isfinite(norm)) {
    throw std::runtime_error("LscmParameterize: singular system (zero extent)");
  }
  const double rect_w = (maxx - minx) / norm;
  const double rect_h = (maxy - miny) / norm;

  // Signed 3D area of a loop.
  const auto loop_area = [&points](const std::vector<int>& loop) {
    double a = 0.0;
    for (size_t k = 0; k < loop.size(); ++k) {
      const Point& p = *points[loop[k]];
      const Point& q = *points[loop[(k + 1) % loop.size()]];
      a += p.x * q.y - p.y * q.x;
    }
    return 0.5 * a;
  };

  // Outer loop = largest absolute area; it is oriented CCW, holes CW.
  size_t outer = 0;
  double outer_area = 0.0;
  for (size_t l = 0; l < loops.size(); ++l) {
    if (std::fabs(loop_area(loops[l])) > std::fabs(outer_area)) {
      outer = l;
      outer_area = loop_area(loops[l]);
    }
  }

  // Convex boundary target UVs.
  LscmResult result;
  result.u.assign(n, 0.0);
  result.v.assign(n, 0.0);
  for (size_t l = 0; l < loops.size(); ++l) {
    std::vector<int>& loop = loops[l];
    double signed_area = loop_area(loop);
    const bool is_outer = (l == outer);
    // Outer loop must be CCW, holes CW, to keep the map orientation
    // preserving.
    if ((is_outer && signed_area < 0.0) || (!is_outer && signed_area > 0.0)) {
      std::reverse(loop.begin(), loop.end());
      signed_area = -signed_area;
    }
    // Chord-length parameter along the loop.
    std::vector<double> cum(loop.size(), 0.0);
    for (size_t k = 1; k < loop.size(); ++k) {
      cum[k] = cum[k - 1] +
               (*points[loop[k]] - *points[loop[k - 1]]).Length();
    }
    const double total =
        cum.back() + (*points[loop[0]] - *points[loop.back()]).Length();
    if (!(total > 0.0)) {
      throw std::runtime_error(
          "LscmParameterize: singular system (degenerate boundary loop)");
    }
    const auto map_to_circle = [&](double orientation) {
      // Chord-length map onto a circle around the loop's normalized 3D
      // centroid. The radius follows the loop's chord length but is capped
      // to fit (with margin) inside the loop's own bounding box: an
      // inflated hole circle would pull interior vertices across the
      // circle and fold the mesh.
      double cx = 0.0, cy = 0.0;
      double bx0 = std::numeric_limits<double>::max();
      double by0 = std::numeric_limits<double>::max();
      double bx1 = std::numeric_limits<double>::lowest();
      double by1 = std::numeric_limits<double>::lowest();
      for (const int idx : loop) {
        const double px = (points[idx]->x - minx) / norm;
        const double py = (points[idx]->y - miny) / norm;
        cx += px;
        cy += py;
        bx0 = std::min(bx0, px);
        by0 = std::min(by0, py);
        bx1 = std::max(bx1, px);
        by1 = std::max(by1, py);
      }
      cx /= loop.size();
      cy /= loop.size();
      const double chord_radius = total / (2.0 * kPi) / norm;
      const double fit_radius =
          0.45 * std::min(bx1 - bx0, by1 - by0);
      const double radius = std::min(chord_radius, fit_radius);
      // Align the rotational phase of the circle with the 3D geometry
      // (optimal rotation in the least-squares sense). Without this the
      // circle is twisted against the surrounding mesh, which folds the
      // triangles between the loop and the rest of the mesh.
      double rot_re = 0.0, rot_im = 0.0;
      for (size_t k = 0; k < loop.size(); ++k) {
        const double theta = 2.0 * kPi * cum[k] / total;
        const double zx = std::cos(theta);
        const double zy = orientation * std::sin(theta);
        const double wx = (points[loop[k]]->x - minx) / norm - cx;
        const double wy = (points[loop[k]]->y - miny) / norm - cy;
        rot_re += wx * zx + wy * zy;
        rot_im += wy * zx - wx * zy;
      }
      const double phase = std::atan2(rot_im, rot_re);
      const double cp = std::cos(phase), sp = std::sin(phase);
      for (size_t k = 0; k < loop.size(); ++k) {
        const double theta = 2.0 * kPi * cum[k] / total;
        const double zx = std::cos(theta);
        const double zy = orientation * std::sin(theta);
        result.u[loop[k]] = cx + radius * (zx * cp - zy * sp);
        result.v[loop[k]] = cy + radius * (zx * sp + zy * cp);
      }
    };
    if (!is_outer) {
      map_to_circle(-1.0); // holes: CW
      continue;
    }
    // Outer loop. Detect sharp convex corners (turning angle above
    // threshold). With exactly four of them the loop is mapped onto the
    // bounding rectangle with the corners aligned to the rectangle's
    // corners (a plain chord-length map would place 3D corners in the
    // middle of rectangle edges, collapsing the adjacent triangles).
    // Otherwise fall back to a corner-free circle map.
    std::vector<size_t> corners;
    for (size_t k = 0; k < loop.size(); ++k) {
      const Point p0 = *points[loop[(k + loop.size() - 1) % loop.size()]];
      const Point p1 = *points[loop[k]];
      const Point p2 = *points[loop[(k + 1) % loop.size()]];
      const Point e0 = p1 - p0;
      const Point e1 = p2 - p1;
      const double turn =
          std::atan2(Cross(e0, e1), Dot(e0, e1)); // CCW loop: convex > 0
      if (turn > 0.3) {
        corners.push_back(k);
      }
    }
    if (corners.size() != 4) {
      map_to_circle(1.0); // CCW
      continue;
    }
    // Rectangle corners in CCW order. Try all four rotations of the
    // detected corners and keep the one whose side lengths best match the
    // rectangle's edge lengths.
    const double rect_x[4] = {0.0, rect_w, rect_w, 0.0};
    const double rect_y[4] = {0.0, 0.0, rect_h, rect_h};
    const double rect_edge[4] = {rect_w, rect_h, rect_w, rect_h};
    double side_len[4];
    for (int m = 0; m < 4; ++m) {
      const size_t k0 = corners[m];
      const size_t k1 = corners[(m + 1) % 4];
      const double len = (k1 > k0 ? cum[k1] - cum[k0]
                                  : total - cum[k0] + cum[k1]);
      side_len[m] = len / total * 2.0 * (rect_w + rect_h);
    }
    int best_rot = 0;
    double best_cost = std::numeric_limits<double>::max();
    for (int rot = 0; rot < 4; ++rot) {
      double cost = 0.0;
      for (int m = 0; m < 4; ++m) {
        const double d = side_len[(m + rot) % 4] - rect_edge[m];
        cost += d * d;
      }
      if (cost < best_cost) {
        best_cost = cost;
        best_rot = rot;
      }
    }
    // Piecewise chord-length map: corner m (matched to rectangle corner
    // (m - best_rot) mod 4) and the loop segment after it onto the
    // corresponding rectangle edge.
    for (int m = 0; m < 4; ++m) {
      const int edge = ((m - best_rot) % 4 + 4) % 4;
      const size_t k0 = corners[m];
      const size_t k1 = corners[(m + 1) % 4];
      const double seg_len = (k1 > k0 ? cum[k1] - cum[k0]
                                      : total - cum[k0] + cum[k1]);
      if (!(seg_len > 0.0)) {
        throw std::runtime_error(
            "LscmParameterize: singular system (degenerate boundary side)");
      }
      const double x0 = rect_x[edge], y0 = rect_y[edge];
      const double dx = rect_x[(edge + 1) % 4] - x0;
      const double dy = rect_y[(edge + 1) % 4] - y0;
      for (size_t k = k0;; k = (k + 1) % loop.size()) {
        const double s = (k >= k0 ? cum[k] - cum[k0]
                                  : total - cum[k0] + cum[k]) /
                         seg_len;
        result.u[loop[k]] = x0 + s * dx;
        result.v[loop[k]] = y0 + s * dy;
        if (k == k1) {
          break;
        }
      }
    }
  }

  // Assemble the weighted Laplacian. triangle_factor optionally re-weights
  // each triangle's contribution (area-distortion minimization).
  const auto build_laplacian = [&](const std::vector<double>* tri_factor) {
    // Accumulate cotangent contributions per undirected edge.
    std::map<std::pair<int, int>, double> edge_w;
    for (size_t t = 0; t < tris.size(); ++t) {
      const double f = tri_factor ? (*tri_factor)[t] : 1.0;
      for (int k = 0; k < 3; ++k) {
        const int i = tris[t][k];
        const int j = tris[t][(k + 1) % 3];
        const int l = tris[t][(k + 2) % 3];
        const Point a = *points[i] - *points[l];
        const Point b = *points[j] - *points[l];
        const double cross = a.x * b.y - a.y * b.x;
        const double w = 0.5 * Dot(a, b) / std::fabs(cross) * f;
        edge_w[std::minmax(i, j)] += w;
      }
    }
    // Symmetric Laplacian: off-diagonal -w (clamped to a small positive
    // weight, keeping the Tutte embedding valid and the matrix positive
    // definite), diagonal = sum of incident weights.
    std::vector<std::map<int, double>> L(n);
    for (const auto& kv : edge_w) {
      const double w = std::max(kv.second, 1e-12);
      const int i = kv.first.first;
      const int j = kv.first.second;
      L[i][j] -= w;
      L[j][i] -= w;
      L[i][i] += w;
      L[j][j] += w;
    }
    return L;
  };

  std::vector<int> free_vertices;
  for (int i = 0; i < n; ++i) {
    if (!is_boundary[i]) {
      free_vertices.push_back(i);
    }
  }

  // Solve the two scalar systems (U and V) with hard Dirichlet boundary.
  const auto solve_hard = [&](const std::vector<std::map<int, double>>& L,
                              std::vector<double>& u, std::vector<double>& v) {
    const int nf = static_cast<int>(free_vertices.size());
    if (nf == 0) {
      return;
    }
    std::vector<int> free_of(n, -1);
    for (int i = 0; i < nf; ++i) {
      free_of[free_vertices[i]] = i;
    }
    std::vector<std::map<int, double>> sub(nf);
    std::vector<double> bu(nf, 0.0), bv(nf, 0.0);
    for (int i = 0; i < nf; ++i) {
      const int gv = free_vertices[i];
      for (const auto& kv : L[gv]) {
        if (kv.first == gv) {
          sub[i][i] += kv.second;
        } else if (free_of[kv.first] >= 0) {
          sub[i][free_of[kv.first]] += kv.second;
        } else {
          // L[gv][fixed] == -w, equation: ... - w*x_fixed moved to RHS.
          bu[i] -= kv.second * u[kv.first];
          bv[i] -= kv.second * v[kv.first];
        }
      }
    }
    const detail::LscmCsrMatrix A = detail::LscmBuildCsr(sub);
    std::vector<double> fu, fv;
    detail::LscmSolveSpd(A, bu, fu, options.solver_tolerance,
                         options.max_solver_iterations);
    detail::LscmSolveSpd(A, bv, fv, options.solver_tolerance,
                         options.max_solver_iterations);
    for (int i = 0; i < nf; ++i) {
      u[free_vertices[i]] = fu[i];
      v[free_vertices[i]] = fv[i];
    }
  };

  // Soft-boundary solve used by the iterations: (L + w_b B) x = w_b B x0
  // where B is the boundary indicator diagonal and x0 the convex target.
  const auto solve_soft = [&](const std::vector<std::map<int, double>>& L,
                              double boundary_weight,
                              const std::vector<double>& tu,
                              const std::vector<double>& tv,
                              std::vector<double>& u, std::vector<double>& v) {
    std::vector<std::map<int, double>> M = L;
    std::vector<double> bu(n, 0.0), bv(n, 0.0);
    for (int i = 0; i < n; ++i) {
      if (is_boundary[i]) {
        M[i][i] += boundary_weight;
        bu[i] = boundary_weight * tu[i];
        bv[i] = boundary_weight * tv[i];
      }
    }
    const detail::LscmCsrMatrix A = detail::LscmBuildCsr(M);
    detail::LscmSolveSpd(A, bu, u, options.solver_tolerance,
                         options.max_solver_iterations);
    detail::LscmSolveSpd(A, bv, v, options.solver_tolerance,
                         options.max_solver_iterations);
  };

  const auto triangle_uv_areas = [&](const std::vector<double>& u,
                                     const std::vector<double>& v) {
    std::vector<double> areas(tris.size());
    for (size_t t = 0; t < tris.size(); ++t) {
      const int i = tris[t][0], j = tris[t][1], k = tris[t][2];
      areas[t] = 0.5 * ((u[j] - u[i]) * (v[k] - v[i]) -
                        (v[j] - v[i]) * (u[k] - u[i]));
    }
    return areas;
  };

  // Initial conformal solve (solve twice: once for U, once for V).
  {
    const auto L = build_laplacian(nullptr);
    solve_hard(L, result.u, result.v);
  }
  const std::vector<double> target_u = result.u;
  const std::vector<double> target_v = result.v;

  // Area-distortion minimization iterations with (soft) boundary fixing.
  for (int it = 0; it < options.iterations; ++it) {
    std::vector<double> auv = triangle_uv_areas(result.u, result.v);
    double total_uv = 0.0;
    for (double a : auv) {
      total_uv += a;
    }
    if (!(total_uv > 0.0) || !std::isfinite(total_uv)) {
      throw std::runtime_error(
          "LscmParameterize: singular system (parameterization collapsed)");
    }
    const double area_scale = total_area3d / total_uv;
    std::vector<double> factor(tris.size());
    double fmean = 0.0;
    for (size_t t = 0; t < tris.size(); ++t) {
      const double scaled_uv =
          std::max(auv[t] * area_scale, 1e-30 * area3d[t]);
      // Boost triangles that are compressed in UV space.
      factor[t] = std::sqrt(area3d[t] / scaled_uv);
      factor[t] = std::min(5.0, std::max(0.2, factor[t]));
      fmean += area3d[t] * factor[t];
    }
    fmean /= total_area3d;
    for (double& f : factor) {
      f /= fmean;
    }
    const auto L = build_laplacian(&factor);
    if (options.boundary_weight > 0.0) {
      solve_soft(L, options.boundary_weight, target_u, target_v, result.u,
                 result.v);
    } else {
      for (int i = 0; i < n; ++i) {
        if (is_boundary[i]) {
          result.u[i] = target_u[i];
          result.v[i] = target_v[i];
        }
      }
      solve_hard(L, result.u, result.v);
    }
  }

  // Quality metrics on the final UVs.
  {
    const std::vector<double> auv = triangle_uv_areas(result.u, result.v);
    double total_uv = 0.0;
    for (double a : auv) {
      total_uv += a;
    }
    if (!(total_uv > 0.0) || !std::isfinite(total_uv)) {
      throw std::runtime_error(
          "LscmParameterize: singular system (parameterization collapsed)");
    }
    const double area_scale = total_area3d / total_uv;
    double accum = 0.0;
    for (size_t t = 0; t < tris.size(); ++t) {
      const double scaled = auv[t] * area_scale;
      if (auv[t] < -1e-12 * area3d[t]) {
        ++result.flipped_triangles;
      }
      const double ratio =
          std::max(scaled, 1e-30 * area3d[t]) / area3d[t];
      const double logr = std::log(ratio);
      accum += area3d[t] * logr * logr;
    }
    result.area_distortion = std::sqrt(accum / total_area3d);
  }
  for (int i = 0; i < n; ++i) {
    if (!std::isfinite(result.u[i]) || !std::isfinite(result.v[i])) {
      throw std::runtime_error(
          "LscmParameterize: singular system (non-finite UVs)");
    }
  }
  return result;
}

} // namespace p2t
