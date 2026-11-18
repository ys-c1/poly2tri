/*
 * Standalone tests for the header-only mesh topology inspector
 * (poly2tri/common/topology.h).
 *
 * This file deliberately avoids the Boost test framework so that it can
 * be built and run with a bare C++ toolchain (see run_topology_test.sh).
 * It contains:
 *   1. Four small hand-built meshes: tetrahedron (closed sphere), square
 *      annulus (one hole / two boundary loops), two disconnected
 *      tetrahedra and a 3-page "book" (non-manifold edge).
 *   2. An independent brute-force O(F^2) statistics implementation that
 *      every small mesh is checked against, value for value.
 *   3. A 100,000-triangle planar grid performance test.
 *
 * Only the C++ standard library is used.
 */

#include <poly2tri/common/topology.h>

#include <algorithm>
#include <array>
#include <cassert>
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {

using Triangle = std::array<int, 3>;
using Mesh = std::vector<Triangle>;

// ---------------------------------------------------------------------
// Independent brute-force reference implementation. Deliberately uses
// std::set / pair scans so nothing is shared with the hash-based code
// under test. O(F^2) component merge, O(F log F) everything else.
// ---------------------------------------------------------------------
struct BruteForceStats {
  std::size_t vertices;
  std::size_t edges;
  std::size_t faces;
  long long chi;
  std::size_t boundary_edges;
  std::size_t non_manifold_edges;
  std::size_t components;
};

BruteForceStats BruteForce(const Mesh& mesh)
{
  BruteForceStats stats{ 0, 0, mesh.size(), 0, 0, 0, 0 };

  std::set<int> used_vertices;
  std::set<std::pair<int, int> > unique_edges;
  for (const Triangle& triangle : mesh) {
    for (int k = 0; k < 3; ++k) {
      used_vertices.insert(triangle[k]);
      int a = triangle[k];
      int b = triangle[(k + 1) % 3];
      if (a > b) {
        std::swap(a, b);
      }
      unique_edges.insert(std::make_pair(a, b));
    }
  }
  stats.vertices = used_vertices.size();
  stats.edges = unique_edges.size();

  for (const std::pair<int, int>& edge : unique_edges) {
    std::size_t count = 0;
    for (const Triangle& triangle : mesh) {
      for (int k = 0; k < 3; ++k) {
        int a = triangle[k];
        int b = triangle[(k + 1) % 3];
        if (a > b) {
          std::swap(a, b);
        }
        if (a == edge.first && b == edge.second) {
          ++count;
          break;
        }
      }
    }
    if (count == 1) {
      ++stats.boundary_edges;
    } else if (count >= 3) {
      ++stats.non_manifold_edges;
    }
  }

  // O(F^2) union of faces that share an edge.
  std::vector<int> component(mesh.size());
  for (std::size_t i = 0; i < mesh.size(); ++i) {
    component[i] = static_cast<int>(i);
  }
  auto CanonicalEdges = [](const Triangle& triangle) {
    std::set<std::pair<int, int> > edges;
    for (int k = 0; k < 3; ++k) {
      int a = triangle[k];
      int b = triangle[(k + 1) % 3];
      if (a > b) {
        std::swap(a, b);
      }
      edges.insert(std::make_pair(a, b));
    }
    return edges;
  };
  std::vector<std::set<std::pair<int, int> > > face_edges(mesh.size());
  for (std::size_t i = 0; i < mesh.size(); ++i) {
    face_edges[i] = CanonicalEdges(mesh[i]);
  }
  for (std::size_t i = 0; i < mesh.size(); ++i) {
    for (std::size_t j = i + 1; j < mesh.size(); ++j) {
      bool shares = false;
      for (const std::pair<int, int>& edge : face_edges[i]) {
        if (face_edges[j].count(edge) != 0) {
          shares = true;
          break;
        }
      }
      if (shares) {
        const int root_j = component[j];
        for (int& root : component) {
          if (root == root_j) {
            root = component[i];
          }
        }
      }
    }
  }
  std::set<int> roots(component.begin(), component.end());
  stats.components = roots.size();

  stats.chi = static_cast<long long>(stats.vertices) -
              static_cast<long long>(stats.edges) +
              static_cast<long long>(stats.faces);
  return stats;
}

// ---------------------------------------------------------------------
// The four acceptance meshes plus a planar disc used by the grid test.
// ---------------------------------------------------------------------

// Tetrahedron surface: 4 vertices, 4 triangles, closed sphere.
Mesh MakeTetrahedron()
{
  return Mesh{
    { 0, 1, 2 }, { 0, 3, 1 }, { 0, 2, 3 }, { 1, 3, 2 }
  };
}

// Square annulus in the plane: outer square (0..3) and inner hole
// (4..7), 8 triangles forming a ring: one hole, two boundary loops.
Mesh MakeAnnulus()
{
  return Mesh{
    { 0, 1, 5 }, { 0, 5, 4 },
    { 1, 2, 6 }, { 1, 6, 5 },
    { 2, 3, 7 }, { 2, 7, 6 },
    { 3, 0, 4 }, { 3, 4, 7 }
  };
}

// Two edge-disconnected tetrahedra using disjoint vertex index ranges.
Mesh MakeTwoTetrahedra()
{
  const Mesh first = MakeTetrahedron();
  Mesh mesh = first;
  for (const Triangle& triangle : first) {
    mesh.push_back(Triangle{ triangle[0] + 4, triangle[1] + 4,
                             triangle[2] + 4 });
  }
  return mesh;
}

// Three triangles fanned around the shared edge (0,1): that edge has
// three incident faces, so it is non-manifold. A "3-page book".
Mesh MakeBook()
{
  return Mesh{
    { 0, 1, 2 }, { 0, 1, 3 }, { 0, 1, 4 }
  };
}

// Planar rectangular triangulation with `cols` x `rows` quad cells,
// two triangles per cell -> 2*cols*rows faces.
Mesh MakeGrid(int cols, int rows)
{
  Mesh mesh;
  mesh.reserve(static_cast<std::size_t>(2) * cols * rows);
  for (int y = 0; y < rows; ++y) {
    for (int x = 0; x < cols; ++x) {
      const int v00 = y * (cols + 1) + x;
      const int v10 = v00 + 1;
      const int v01 = v00 + cols + 1;
      const int v11 = v01 + 1;
      mesh.push_back(Triangle{ v00, v10, v11 });
      mesh.push_back(Triangle{ v00, v11, v01 });
    }
  }
  return mesh;
}

// ---------------------------------------------------------------------
// Checking helpers
// ---------------------------------------------------------------------
int g_failures = 0;
int g_checks = 0;

void Check(bool condition, const std::string& message)
{
  ++g_checks;
  if (!condition) {
    ++g_failures;
    std::cout << "  [FAIL] " << message << "\n";
  }
}

void CheckAgainstBruteForce(const std::string& name, const Mesh& mesh)
{
  const p2t::topology::TopologyReport report =
      p2t::topology::AnalyzeMesh(Mesh(mesh));
  const BruteForceStats expected = BruteForce(mesh);

  std::cout << "== " << name << " ==\n";
  std::cout << p2t::topology::Summary(report) << "\n";

  Check(report.num_faces == expected.faces,
      name + ": F mismatch");
  Check(report.num_vertices == expected.vertices,
      name + ": V mismatch");
  Check(report.num_edges == expected.edges,
      name + ": E mismatch");
  Check(report.euler_characteristic == expected.chi,
      name + ": chi mismatch");
  Check(report.num_boundary_edges == expected.boundary_edges,
      name + ": boundary edge count mismatch");
  Check(report.num_non_manifold_edges == expected.non_manifold_edges,
      name + ": non-manifold edge count mismatch");
  Check(report.num_components == expected.components,
      name + ": component count mismatch");

  // The boundary/non-manifold edge lists must agree with the counts and
  // contain each classified edge exactly once.
  Check(report.boundary_edges.size() == report.num_boundary_edges,
      name + ": boundary edge list size");
  Check(report.non_manifold_edges.size() ==
      report.num_non_manifold_edges,
      name + ": non-manifold edge list size");
  for (const p2t::topology::EdgeInfo& edge : report.boundary_edges) {
    Check(edge.face_count == 1, name + ": boundary edge face_count");
  }
  for (const p2t::topology::EdgeInfo& edge :
      report.non_manifold_edges) {
    Check(edge.face_count >= 3,
        name + ": non-manifold edge face_count");
  }

  // Per-component V/E/F must add up to the global counts (each edge and
  // each face belongs to exactly one component; pinched vertices may be
  // counted in more than one component, hence the >= relation).
  std::size_t sum_v = 0;
  std::size_t sum_e = 0;
  std::size_t sum_f = 0;
  long long sum_chi = 0;
  for (const p2t::topology::ComponentReport& component :
      report.components) {
    sum_v += component.num_vertices;
    sum_e += component.num_edges;
    sum_f += component.num_faces;
    sum_chi += component.euler_characteristic;
    Check(component.num_edges ==
            component.num_vertices + component.num_faces -
            static_cast<std::size_t>(component.euler_characteristic),
        name + ": component chi identity");
  }
  Check(sum_e == report.num_edges, name + ": component edges sum");
  Check(sum_f == report.num_faces, name + ": component faces sum");
  Check(sum_v >= report.num_vertices, name + ": component vertices sum");
  Check(sum_chi >= report.euler_characteristic,
      name + ": component chi sum");
}

// ---------------------------------------------------------------------
// Entry point
// ---------------------------------------------------------------------
}  // namespace

int main()
{
  std::cout << "==== Small meshes: checked against brute force ====\n\n";

  // 1. Tetrahedron: closed sphere, chi=2, genus 0.
  CheckAgainstBruteForce("tetrahedron (closed sphere)",
      MakeTetrahedron());
  {
    const auto report = p2t::topology::AnalyzeMesh(MakeTetrahedron());
    Check(report.is_manifold, "tetrahedron: must be manifold");
    Check(report.num_components == 1, "tetrahedron: one component");
    Check(report.components[0].is_closed,
        "tetrahedron: component is closed");
    Check(report.components[0].genus == 0,
        "tetrahedron: genus 0");
    Check(report.euler_characteristic == 2,
        "tetrahedron: chi == 2 (sphere)");
  }

  // 2. Square annulus: one hole, chi=0, two boundary loops, genus 0.
  CheckAgainstBruteForce("square annulus (one hole)",
      MakeAnnulus());
  {
    const auto report = p2t::topology::AnalyzeMesh(MakeAnnulus());
    Check(report.is_manifold, "annulus: must be manifold");
    Check(report.num_boundary_edges == 8,
        "annulus: 8 boundary edges");
    Check(report.components[0].num_boundary_loops == 2,
        "annulus: two boundary loops (outer + hole)");
    Check(report.components[0].genus == 0,
        "annulus: genus 0");
    Check(report.euler_characteristic == 0,
        "annulus: chi == 0");
  }

  // 3. Two disconnected tetrahedra: two components, chi=4.
  CheckAgainstBruteForce("two disconnected tetrahedra",
      MakeTwoTetrahedra());
  {
    const auto report =
        p2t::topology::AnalyzeMesh(MakeTwoTetrahedra());
    Check(report.num_components == 2,
        "two tetrahedra: two components");
    Check(report.components[0].num_faces == 4 &&
          report.components[1].num_faces == 4,
        "two tetrahedra: 4 faces per component");
    Check(report.components[0].euler_characteristic == 2 &&
          report.components[1].euler_characteristic == 2,
        "two tetrahedra: each component chi=2");
    Check(report.euler_characteristic == 4,
        "two tetrahedra: global chi=4");
  }

  // 4. Three-page book: one edge shared by 3 triangles -> non-manifold.
  CheckAgainstBruteForce("non-manifold book (edge shared by 3 faces)",
      MakeBook());
  {
    const auto report = p2t::topology::AnalyzeMesh(MakeBook());
    Check(!report.is_manifold, "book: must be non-manifold");
    Check(report.num_non_manifold_edges == 1,
        "book: exactly one non-manifold edge");
    Check(report.non_manifold_edges[0].a == 0 &&
          report.non_manifold_edges[0].b == 1 &&
          report.non_manifold_edges[0].face_count == 3,
        "book: edge (0,1) used by 3 faces");
    Check(report.num_boundary_edges == 6,
        "book: 6 boundary edges");
    Check(report.components[0].genus == -1,
        "book: genus undefined for non-manifold component");
  }

  // Extra: empty mesh must be accepted and reported as empty.
  {
    const auto report = p2t::topology::AnalyzeMesh(Mesh{});
    Check(report.num_faces == 0 && report.num_edges == 0 &&
          report.num_vertices == 0 && report.num_components == 0 &&
          report.is_manifold,
        "empty mesh: zeroed manifold report");
  }

  // Extra: two triangles touching at a single vertex. Every edge is
  // manifold but they form two edge-connected components sharing a
  // pinched vertex.
  {
    const Mesh pinched{ { 0, 1, 2 }, { 0, 3, 4 } };
    CheckAgainstBruteForce("pinched vertex (two triangles, one shared vertex)",
        pinched);
    const auto report = p2t::topology::AnalyzeMesh(pinched);
    Check(report.num_components == 2, "pinched: two components");
    Check(report.is_manifold, "pinched: all edges still manifold");
    Check(report.num_vertices == 5, "pinched: 5 distinct vertices");
    Check(report.components[0].num_vertices +
          report.components[1].num_vertices == 6,
        "pinched: shared vertex counted once per component (6 vs V=5)");
  }

  // Extra: torus mesh (closed, genus 1 -> chi = 2 - 2g = 0).
  // Built as a parametric grid with both dimensions wrapped:
  // cols*rows vertices, cols*rows*2 triangles, no boundary edges.
  {
    const int cols = 8;
    const int rows = 6;
    Mesh torus;
    auto torus_id = [](int x, int y) {
      return (y % rows) * cols + (x % cols);
    };
    for (int y = 0; y < rows; ++y) {
      for (int x = 0; x < cols; ++x) {
        const int v00 = torus_id(x, y);
        const int v10 = torus_id(x + 1, y);
        const int v01 = torus_id(x, y + 1);
        const int v11 = torus_id(x + 1, y + 1);
        torus.push_back(Triangle{ v00, v10, v11 });
        torus.push_back(Triangle{ v00, v11, v01 });
      }
    }
    CheckAgainstBruteForce("torus (closed genus 1)", torus);
    const auto report = p2t::topology::AnalyzeMesh(torus);
    Check(report.is_manifold && report.components[0].is_closed,
        "torus: closed manifold");
    Check(report.num_boundary_edges == 0, "torus: no boundary edges");
    Check(report.euler_characteristic == 0, "torus: chi == 0");
    Check(report.components[0].genus == 1, "torus: genus 1");
  }

  // Extra: sparse vertex indices force the hash-map fallback inside
  // VertexComponentMap; results must stay identical.
  {
    Mesh sparse = MakeTetrahedron();
    for (int& vertex : sparse[0]) {
      vertex += 1000000;
    }
    for (int& vertex : sparse[1]) {
      vertex += 1000000;
    }
    for (int& vertex : sparse[2]) {
      vertex += 1000000;
    }
    for (int& vertex : sparse[3]) {
      vertex += 1000000;
    }
    const auto report = p2t::topology::AnalyzeMesh(sparse);
    Check(report.num_vertices == 4 && report.num_edges == 6 &&
          report.num_faces == 4 && report.euler_characteristic == 2 &&
          report.num_components == 1 && report.is_manifold,
        "sparse indices: same topology as tetrahedron");
  }

  // Extra: malformed input must throw instead of producing garbage.
  {
    const Mesh degenerate{ { 0, 0, 1 }, { 0, 1, 2 } };
    bool threw = false;
    try {
      p2t::topology::AnalyzeMesh(degenerate);
    } catch (const std::invalid_argument&) {
      threw = true;
    }
    Check(threw, "degenerate triangle: invalid_argument thrown");
  }

  std::cout << "==== Performance: 100,000-face planar grid ====\n\n";

  // A medium grid is still fully cross-checked with the O(F^2) brute
  // force reference; the 100k grid is checked via analytic invariants
  // (its counts are known exactly) and used purely for timing.
  const Mesh medium_grid = MakeGrid(50, 40);  // 4000 faces
  CheckAgainstBruteForce("medium planar grid (4000 faces, full brute force)",
      medium_grid);

  const Mesh grid = MakeGrid(250, 200);  // 250*200*2 = 100000 faces
  assert(grid.size() == 100000);

  const auto warm = p2t::topology::AnalyzeMesh(grid);
  Check(warm.num_faces == 100000, "grid: face count");
  Check(warm.is_manifold, "grid: planar mesh is manifold");
  Check(warm.num_components == 1, "grid: one component");
  Check(warm.num_non_manifold_edges == 0,
      "grid: no non-manifold edges");
  Check(warm.num_boundary_edges == 2 * (250 + 200),
      "grid: perimeter boundary edges");
  Check(warm.num_vertices == 251 * 201, "grid: vertex count");
  Check(warm.euler_characteristic == 1,
      "grid: disc chi == 1");
  Check(warm.components[0].num_boundary_loops == 1,
      "grid: one boundary loop");
  Check(warm.components[0].genus == 0, "grid: genus 0");

  const int iterations = 9;
  double min_ms = 1e300;
  double max_ms = 0.0;
  double sum_ms = 0.0;
  for (int i = 0; i < iterations; ++i) {
    const auto t0 = std::chrono::high_resolution_clock::now();
    const auto report = p2t::topology::AnalyzeMesh(grid);
    const auto t1 = std::chrono::high_resolution_clock::now();
    if (report.num_faces != grid.size()) {
      Check(false, "grid: timed run returned inconsistent result");
    }
    const double ms =
        std::chrono::duration<double, std::milli>(t1 - t0).count();
    min_ms = std::min(min_ms, ms);
    max_ms = std::max(max_ms, ms);
    sum_ms += ms;
  }
  const double avg_ms = sum_ms / iterations;
  std::cout << "grid V=" << warm.num_vertices
            << " E=" << warm.num_edges
            << " F=" << warm.num_faces
            << " chi=" << warm.euler_characteristic
            << " boundary=" << warm.num_boundary_edges
            << " components=" << warm.num_components << "\n";
  std::cout << "AnalyzeMesh x" << iterations
            << ": min=" << min_ms << " ms"
            << " avg=" << avg_ms << " ms"
            << " max=" << max_ms << " ms"
            << " (budget: 100 ms)\n";
  Check(min_ms <= 100.0, "grid: min runtime must stay under 100 ms");

  std::cout << "\n==== Result: " << g_checks << " checks, "
            << g_failures << " failures ====\n";
  return g_failures == 0 ? 0 : 1;
}
