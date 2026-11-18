/*
 * Poly2Tri Copyright (c) 2009-2022, Poly2Tri Contributors
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
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#pragma once

// Header-only mesh topology inspection for indexed triangle meshes.
//
// Given a triangle soup (each triangle is three vertex indices) this
// computes:
//   * V / E / F counts and the Euler characteristic chi = V - E + F
//   * boundary edges (used by exactly 1 triangle)
//   * non-manifold edges (used by 3 or more triangles)
//   * connected components (two triangles are connected when they share
//     an edge), with per-component V/E/F, chi, boundary loop count and
//     genus derived from the Euler characteristic
//
// Complexity: expected O(F) time and O(V + E) memory.
//
// Scope note: "manifold" here means edge-manifold, i.e. every edge is
// incident to at most two triangles, matching the task definition. Two
// sheets touching at a single vertex (a pinched vertex) keep all edges
// manifold but form separate edge-connected components; this is visible
// in the component breakdown and is not classified as a non-manifold
// edge. Orientation consistency (orientability) is not checked.
//
// The implementation only uses the C++ standard library.

#include <array>
#include <cstddef>
#include <cstdint>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace p2t {
namespace topology {

// One undirected edge (canonicalized so that a < b) and how many
// triangles reference it.
struct EdgeInfo {
  std::size_t a;
  std::size_t b;
  int face_count;
};

// Topology verdict for a single edge-connected component.
struct ComponentReport {
  std::size_t index = 0;
  std::size_t num_vertices = 0;
  std::size_t num_edges = 0;
  std::size_t num_faces = 0;
  long long euler_characteristic = 0;  // V - E + F
  std::size_t num_boundary_edges = 0;
  std::size_t num_non_manifold_edges = 0;
  long long num_boundary_loops = -1;   // -1 = not applicable / irregular
  long long genus = -1;                // -1 = undefined (non-manifold)
  bool is_manifold = false;            // every edge used by at most 2 faces
  bool is_closed = false;              // manifold and no boundary edges
  std::string conclusion;
};

// Whole-mesh report. Global counts describe the mesh as given; per
// component counts describe each edge-connected component separately.
struct TopologyReport {
  std::size_t num_vertices = 0;
  std::size_t num_edges = 0;
  std::size_t num_faces = 0;
  long long euler_characteristic = 0;  // V - E + F
  std::size_t num_boundary_edges = 0;
  std::size_t num_non_manifold_edges = 0;
  std::size_t num_components = 0;
  bool is_manifold = false;
  std::string conclusion;
  std::vector<EdgeInfo> boundary_edges;
  std::vector<EdgeInfo> non_manifold_edges;
  std::vector<ComponentReport> components;
};

namespace detail {

static const std::uint32_t kNoComponent = UINT32_MAX;

inline std::uint64_t PackEdge(std::uint64_t a, std::uint64_t b)
{
  if (a > b) {
    const std::uint64_t tmp = a;
    a = b;
    b = tmp;
  }
  return (a << 32) | b;
}

struct DisjointSet {
  std::vector<std::uint32_t> parent;
  std::vector<std::uint32_t> rank_;

  explicit DisjointSet(std::uint32_t n) : parent(n), rank_(n, 0)
  {
    std::iota(parent.begin(), parent.end(), 0);
  }

  std::uint32_t Find(std::uint32_t x)
  {
    while (parent[x] != x) {
      parent[x] = parent[parent[x]];
      x = parent[x];
    }
    return x;
  }

  void Unite(std::uint32_t x, std::uint32_t y)
  {
    x = Find(x);
    y = Find(y);
    if (x == y) {
      return;
    }
    if (rank_[x] < rank_[y]) {
      parent[x] = y;
    } else if (rank_[x] > rank_[y]) {
      parent[y] = x;
    } else {
      parent[y] = x;
      ++rank_[x];
    }
  }
};

struct EdgeRecord {
  std::uint32_t count;
  std::uint32_t first_face;
};

// Maps a vertex index to a component id. Uses a dense vector when the
// index range is compact (the common case) and falls back to a hash map
// for sparse index sets so non-compact inputs cannot blow up memory.
template <typename Index>
class VertexComponentMap {
public:
  VertexComponentMap(Index max_index, std::size_t num_faces)
  {
    const std::uint64_t span = static_cast<std::uint64_t>(max_index) + 1;
    dense_mode_ = span <= 16 * num_faces + 4096;
    if (dense_mode_) {
      dense_.assign(static_cast<std::size_t>(span), kNoComponent);
    } else {
      sparse_.reserve(num_faces * 2);
    }
  }

  std::uint32_t Get(Index vertex) const
  {
    if (dense_mode_) {
      return dense_[static_cast<std::size_t>(vertex)];
    }
    typename std::unordered_map<Index, std::uint32_t>::const_iterator it =
        sparse_.find(vertex);
    return it == sparse_.end() ? kNoComponent : it->second;
  }

  void Set(Index vertex, std::uint32_t component)
  {
    if (dense_mode_) {
      dense_[static_cast<std::size_t>(vertex)] = component;
    } else {
      sparse_[vertex] = component;
    }
  }

  std::size_t CountNonEmpty() const
  {
    if (dense_mode_) {
      std::size_t count = 0;
      for (std::uint32_t value : dense_) {
        if (value != kNoComponent) {
          ++count;
        }
      }
      return count;
    }
    return sparse_.size();
  }

private:
  bool dense_mode_;
  std::vector<std::uint32_t> dense_;
  std::unordered_map<Index, std::uint32_t> sparse_;
};

struct BoundaryVertex {
  std::size_t neighbor[2];
  std::uint32_t degree;
};

// Counts the boundary loops of one edge-manifold component by walking
// its boundary edges. In a manifold mesh every boundary vertex has
// exactly two boundary edges, so the edges close into disjoint loops.
// Returns -1 if an irregular boundary is encountered.
inline long long CountBoundaryLoops(
    const std::vector<std::uint64_t>& boundary_keys)
{
  if (boundary_keys.empty()) {
    return 0;
  }

  std::unordered_map<std::size_t, BoundaryVertex> adjacency;
  adjacency.reserve(boundary_keys.size() * 2);
  for (std::uint64_t key : boundary_keys) {
    const std::size_t a = static_cast<std::size_t>(key >> 32);
    const std::size_t b = static_cast<std::size_t>(key & 0xFFFFFFFFULL);
    BoundaryVertex& va = adjacency[a];
    BoundaryVertex& vb = adjacency[b];
    if (va.degree >= 2 || vb.degree >= 2) {
      return -1;
    }
    va.neighbor[va.degree++] = b;
    vb.neighbor[vb.degree++] = a;
  }

  std::unordered_set<std::uint64_t> visited;
  visited.reserve(boundary_keys.size() * 2);
  long long loops = 0;

  for (std::uint64_t key : boundary_keys) {
    if (visited.count(key) != 0) {
      continue;
    }
    ++loops;

    std::size_t start = static_cast<std::size_t>(key >> 32);
    std::size_t previous = start;
    std::size_t current = static_cast<std::size_t>(key & 0xFFFFFFFFULL);
    visited.insert(key);

    while (current != start) {
      typename std::unordered_map<std::size_t, BoundaryVertex>::iterator it =
          adjacency.find(current);
      if (it == adjacency.end() || it->second.degree != 2) {
        return -1;
      }
      const BoundaryVertex& bv = it->second;
      const std::size_t next =
          bv.neighbor[0] == previous ? bv.neighbor[1] : bv.neighbor[0];
      if (visited.insert(PackEdge(current, next)).second == false) {
        return -1;
      }
      previous = current;
      current = next;
    }
  }
  return loops;
}

inline std::string MakeConclusion(const ComponentReport& component)
{
  std::ostringstream out;
  if (!component.is_manifold) {
    out << "non-manifold: " << component.num_non_manifold_edges
        << " edge(s) shared by 3+ triangles; genus is undefined";
    return out.str();
  }
  if (component.is_closed) {
    const long long twice_genus = 2 - component.euler_characteristic;
    if (twice_genus < 0 || (twice_genus & 1LL) != 0) {
      out << "manifold, closed, but chi=" << component.euler_characteristic
          << " is inconsistent with an orientable closed surface";
    } else {
      out << "closed manifold surface, genus=" << component.genus
          << " (chi=" << component.euler_characteristic
          << ", expected 2-2g)";
    }
    return out.str();
  }

  const long long twice_genus =
      2 - component.num_boundary_loops - component.euler_characteristic;
  if (component.num_boundary_loops < 0 || twice_genus < 0 ||
      (twice_genus & 1LL) != 0) {
    out << "manifold with boundary, but chi="
        << component.euler_characteristic
        << " with " << component.num_boundary_edges
        << " boundary edges is inconsistent with an orientable surface";
  } else {
    out << "manifold surface with boundary: "
        << component.num_boundary_loops << " boundary loop(s), genus="
        << component.genus << " (chi=" << component.euler_characteristic
        << ", expected 2-2g-b)";
  }
  return out.str();
}

}  // namespace detail

// Analyzes an indexed triangle mesh. Each element of `triangles` holds
// three vertex indices; vertex indices must fit in 32 bits and the three
// indices of one triangle must be distinct.
//
// Expected O(F) complexity. Throws std::invalid_argument on malformed
// input (empty mesh is allowed and yields an empty report).
template <typename Index>
inline TopologyReport AnalyzeMesh(
    const std::vector<std::array<Index, 3> >& triangles)
{
  TopologyReport report;
  const std::size_t num_faces = triangles.size();
  report.num_faces = num_faces;
  if (num_faces == 0) {
    report.num_components = 0;
    report.is_manifold = true;
    report.conclusion = "empty mesh";
    return report;
  }
  if (num_faces >= UINT32_MAX) {
    throw std::invalid_argument(
        "AnalyzeMesh: too many triangles (32-bit face id exhausted)");
  }

  Index max_index = 0;
  for (std::size_t f = 0; f < num_faces; ++f) {
    for (int k = 0; k < 3; ++k) {
      const Index vertex = triangles[f][k];
      if (static_cast<std::uint64_t>(vertex) >= 0xFFFFFFFFULL) {
        throw std::invalid_argument(
            "AnalyzeMesh: vertex index does not fit in 32 bits");
      }
      if (vertex > max_index) {
        max_index = vertex;
      }
    }
    if (triangles[f][0] == triangles[f][1] ||
        triangles[f][1] == triangles[f][2] ||
        triangles[f][0] == triangles[f][2]) {
      throw std::invalid_argument(
          "AnalyzeMesh: triangle with repeated vertex indices");
    }
  }

  // Single hash pass: count faces per undirected edge and union the
  // faces sharing an edge (builds the edge-connected components).
  std::unordered_map<std::uint64_t, detail::EdgeRecord> edge_records;
  edge_records.reserve(num_faces * 4);
  edge_records.max_load_factor(0.7f);

  detail::DisjointSet dsu(static_cast<std::uint32_t>(num_faces));

  for (std::uint32_t f = 0; f < num_faces; ++f) {
    const std::uint64_t v0 = static_cast<std::uint64_t>(triangles[f][0]);
    const std::uint64_t v1 = static_cast<std::uint64_t>(triangles[f][1]);
    const std::uint64_t v2 = static_cast<std::uint64_t>(triangles[f][2]);
    const std::uint64_t keys[3] = {
      detail::PackEdge(v0, v1),
      detail::PackEdge(v1, v2),
      detail::PackEdge(v2, v0)
    };
    for (int k = 0; k < 3; ++k) {
      std::pair<typename std::unordered_map<
          std::uint64_t, detail::EdgeRecord>::iterator, bool> inserted =
          edge_records.insert(std::make_pair(
              keys[k], detail::EdgeRecord{ 1, f }));
      if (!inserted.second) {
        detail::EdgeRecord& record = inserted.first->second;
        ++record.count;
        dsu.Unite(f, record.first_face);
      }
    }
  }

  // Compact component ids from DSU roots.
  std::vector<std::uint32_t> component_of_face(num_faces);
  std::unordered_map<std::uint32_t, std::uint32_t> root_to_component;
  root_to_component.reserve(num_faces);
  for (std::uint32_t f = 0; f < num_faces; ++f) {
    const std::uint32_t root = dsu.Find(f);
    std::pair<typename std::unordered_map<std::uint32_t,
          std::uint32_t>::iterator, bool> inserted =
        root_to_component.insert(
            std::make_pair(root,
              static_cast<std::uint32_t>(root_to_component.size())));
    component_of_face[f] = inserted.first->second;
  }

  const std::uint32_t num_components =
      static_cast<std::uint32_t>(root_to_component.size());
  report.num_components = num_components;
  report.components.resize(num_components);
  for (std::uint32_t c = 0; c < num_components; ++c) {
    report.components[c].index = c;
  }

  std::vector<std::size_t> component_vertices(num_components, 0);
  std::vector<std::size_t> component_edges(num_components, 0);
  std::vector<std::size_t> component_boundary(num_components, 0);
  std::vector<std::size_t> component_non_manifold(num_components, 0);
  std::vector<std::vector<std::uint64_t> > component_boundary_keys(
      num_components);

  // Face and vertex counts per component. A vertex belongs to one
  // component except for a pinched vertex shared by edge-disconnected
  // sheets; those rare conflicts are deduplicated explicitly.
  detail::VertexComponentMap<Index> vertex_component(max_index, num_faces);
  std::unordered_set<std::uint64_t> shared_vertex_pairs;
  for (std::uint32_t f = 0; f < num_faces; ++f) {
    const std::uint32_t component = component_of_face[f];
    ++report.components[component].num_faces;
    for (int k = 0; k < 3; ++k) {
      const Index vertex = triangles[f][k];
      const std::uint32_t owner = vertex_component.Get(vertex);
      if (owner == detail::kNoComponent) {
        vertex_component.Set(vertex, component);
        ++component_vertices[component];
      } else if (owner != component) {
        const std::uint64_t pair_key =
            (static_cast<std::uint64_t>(component) << 32) |
            static_cast<std::uint64_t>(vertex);
        if (shared_vertex_pairs.insert(pair_key).second) {
          ++component_vertices[component];
        }
      }
    }
  }
  report.num_vertices = vertex_component.CountNonEmpty();

  // Edge counts and edge classification per component and globally.
  for (typename std::unordered_map<std::uint64_t,
        detail::EdgeRecord>::const_iterator it = edge_records.begin();
      it != edge_records.end(); ++it) {
    const std::uint64_t key = it->first;
    const detail::EdgeRecord& record = it->second;
    const std::uint32_t component =
        component_of_face[record.first_face];
    ++component_edges[component];

    const EdgeInfo info{
      static_cast<std::size_t>(key >> 32),
      static_cast<std::size_t>(key & 0xFFFFFFFFULL),
      static_cast<int>(record.count)
    };
    if (record.count == 1) {
      ++report.num_boundary_edges;
      ++component_boundary[component];
      report.boundary_edges.push_back(info);
      component_boundary_keys[component].push_back(key);
    } else if (record.count >= 3) {
      ++report.num_non_manifold_edges;
      ++component_non_manifold[component];
      report.non_manifold_edges.push_back(info);
    }
  }

  report.num_edges = edge_records.size();
  report.euler_characteristic =
      static_cast<long long>(report.num_vertices) -
      static_cast<long long>(report.num_edges) +
      static_cast<long long>(report.num_faces);
  report.is_manifold = report.num_non_manifold_edges == 0;

  // Assemble per-component reports and derive loops / genus.
  for (std::uint32_t c = 0; c < num_components; ++c) {
    ComponentReport& component = report.components[c];
    component.num_vertices = component_vertices[c];
    component.num_edges = component_edges[c];
    component.num_boundary_edges = component_boundary[c];
    component.num_non_manifold_edges = component_non_manifold[c];
    component.euler_characteristic =
        static_cast<long long>(component.num_vertices) -
        static_cast<long long>(component.num_edges) +
        static_cast<long long>(component.num_faces);
    component.is_manifold = component.num_non_manifold_edges == 0;
    component.is_closed =
        component.is_manifold && component.num_boundary_edges == 0;

    if (component.is_manifold) {
      component.num_boundary_loops =
          detail::CountBoundaryLoops(component_boundary_keys[c]);
      if (component.is_closed) {
        // chi = 2 - 2g
        component.genus =
            (2 - component.euler_characteristic) / 2;
      } else if (component.num_boundary_loops >= 0) {
        // chi = 2 - 2g - b, where b is the number of boundary loops
        component.genus =
            (2 - component.num_boundary_loops -
              component.euler_characteristic) / 2;
      }
    }
    component.conclusion = detail::MakeConclusion(component);
  }

  std::ostringstream out;
  if (report.is_manifold) {
    out << "edge-manifold mesh, " << num_components << " component(s), chi="
        << report.euler_characteristic << ", "
        << report.num_boundary_edges << " boundary edge(s)";
  } else {
    out << "NON-MANIFOLD mesh: " << report.num_non_manifold_edges
        << " edge(s) shared by 3+ triangles across " << num_components
        << " component(s)";
  }
  report.conclusion = out.str();
  return report;
}

// Multi-line human-readable rendering of a TopologyReport, suitable for
// mesh quality reports.
inline std::string Summary(const TopologyReport& report)
{
  std::ostringstream out;
  out << "V=" << report.num_vertices << " E=" << report.num_edges
      << " F=" << report.num_faces
      << " chi=" << report.euler_characteristic
      << " boundary_edges=" << report.num_boundary_edges
      << " non_manifold_edges=" << report.num_non_manifold_edges
      << " components=" << report.num_components << "\n";
  for (std::size_t c = 0; c < report.components.size(); ++c) {
    const ComponentReport& component = report.components[c];
    out << "  [#" << component.index << "] V=" << component.num_vertices
        << " E=" << component.num_edges
        << " F=" << component.num_faces
        << " chi=" << component.euler_characteristic
        << " boundary=" << component.num_boundary_edges
        << " loops=" << component.num_boundary_loops
        << " non_manifold=" << component.num_non_manifold_edges
        << " genus=" << component.genus
        << " -> " << component.conclusion << "\n";
  }
  out << "  verdict: " << report.conclusion;
  return out.str();
}

}  // namespace topology
}  // namespace p2t
