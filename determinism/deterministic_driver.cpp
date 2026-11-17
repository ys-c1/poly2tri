/*
 * Deterministic-output verification driver for poly2tri.
 *
 * Self-contained: uses only the C++ standard library and the poly2tri
 * headers/sources shipped in this repository. No third-party dependency is
 * added. Hashes are produced externally by sha256sum over the serialized
 * byte stream emitted here.
 *
 * Usage:
 *   deterministic_driver <fixture-name>
 *
 * Fixtures (all tiny, generated in-code, no data files):
 *   regular_grid   - square boundary + interior Steiner grid
 *   shared_point   - two holes touching at one shared vertex
 *   near_coincident- vertices separated by sub-eps / near-equal coordinates
 *   hole           - polygon with one interior hole
 *
 * Output is a fully canonical, locale-independent text serialization:
 *   - every vertex referenced by a canonical index assigned by explicit rules
 *     (quantized (qx,qy) ascending, then insertion order for ties);
 *   - coordinates printed with %a (hexadecimal float, exact IEEE-754 round
 *     trip, identical bytes on every conforming implementation);
 *   - triangles sorted by their canonical index triples;
 *   - constrained edges sorted, with per-edge constrained-flag counts;
 *   - a CONSTRAINTS section whose PASS/FAIL line reports constraint integrity.
 */

#include <poly2tri/poly2tri.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

// ---------------------------------------------------------------------------
// Quantization (mirrors the rule compiled into the library; default 0 = exact)
// Both the library and this driver read the same P2T_QUANTIZATION_PRECISION
// environment variable, so canonical indices always follow the same rule as
// the library's vertex merging.
// ---------------------------------------------------------------------------
int DriverPrecision()
{
  const char* raw = std::getenv("P2T_QUANTIZATION_PRECISION");
  if (raw == nullptr) return 0;
  return std::stoi(raw); // validated by the library before triangulation
}

unsigned long long Quantize(double v)
{
  const int precision = DriverPrecision();
  if (precision == 0) {
    // Exact mode: the raw ordered bit pattern of the double; -0.0 == +0.0.
    unsigned long long bits;
    std::memcpy(&bits, &v, sizeof(bits));
    if (bits == 0x8000000000000000ULL) bits = 0ULL;
    const unsigned long long mask = 0ULL - (bits >> 63);
    return bits ^ (mask | (1ULL << 63));
  }
  return static_cast<unsigned long long>(
             std::llround(v * std::pow(10.0, precision)))
         + (1ULL << 63);
}

// Canonical point key: quantized coordinates first, insertion order second.
struct PointKey {
  unsigned long long qx, qy;
  std::size_t order;
};

struct PointKeyLess {
  bool operator()(const PointKey& a, const PointKey& b) const
  {
    if (a.qy != b.qy) return a.qy < b.qy;
    if (a.qx != b.qx) return a.qx < b.qx;
    return a.order < b.order;
  }
};

// Exact, locale-independent coordinate rendering.
std::string HexDouble(double v)
{
  char buf[64];
  std::snprintf(buf, sizeof(buf), "%a", v);
  return std::string(buf);
}

// ---------------------------------------------------------------------------
// Fixtures
// ---------------------------------------------------------------------------
struct Fixture {
  std::vector<p2t::Point*> boundary;
  std::vector<std::vector<p2t::Point*>> holes;
  std::vector<p2t::Point*> steiner;
  // Expected set of constrained edges, as pairs of vertex coordinates
  // (raw doubles, exactly as given to the CDT).
  std::vector<std::pair<p2t::Point*, p2t::Point*>> expected_constraints;
};

void CloseConstraints(Fixture& f, const std::vector<p2t::Point*>& loop)
{
  for (std::size_t i = 0; i < loop.size(); ++i) {
    f.expected_constraints.emplace_back(loop[i], loop[(i + 1) % loop.size()]);
  }
}

Fixture MakeFixture(const std::string& name)
{
  Fixture f;
  if (name == "regular_grid") {
    // Square boundary (CCW), 2x2 interior Steiner points => 9 triangles.
    p2t::Point* a = new p2t::Point(0.0, 0.0);
    p2t::Point* b = new p2t::Point(2.0, 0.0);
    p2t::Point* c = new p2t::Point(2.0, 2.0);
    p2t::Point* d = new p2t::Point(0.0, 2.0);
    f.boundary = { a, b, c, d };
    f.steiner = {
      new p2t::Point(1.0, 1.0), new p2t::Point(1.0, 0.5),
      new p2t::Point(0.5, 1.0), new p2t::Point(1.5, 1.0)
    };
    CloseConstraints(f, f.boundary);
  } else if (name == "shared_point") {
    // Triangle boundary; two holes that share one vertex (supported input).
    f.boundary = {
      new p2t::Point(0.0, 0.0), new p2t::Point(4.0, 0.0),
      new p2t::Point(2.0, 4.0)
    };
    p2t::Point* shared = new p2t::Point(1.5, 1.2);
    f.holes.push_back({
      new p2t::Point(1.2, 1.0), new p2t::Point(1.8, 1.0), shared
    });
    f.holes.push_back({
      shared, new p2t::Point(1.8, 1.6), new p2t::Point(1.2, 1.6)
    });
    CloseConstraints(f, f.boundary);
    CloseConstraints(f, f.holes[0]);
    CloseConstraints(f, f.holes[1]);
  } else if (name == "near_coincident") {
    // Near-coincident consecutive boundary vertices: two points separated
    // by ~2e-9 (well above EPSILON 1e-12) with tiny perpendicular
    // perturbations so no three boundary points are collinear (poly2tri
    // skips collinear constraint-edge vertices by design). Plus a Steiner
    // point whose x nearly matches the near-coincident pair's x.
    f.boundary = {
      new p2t::Point(0.0, 0.0),
      new p2t::Point(3.0, 0.0),
      new p2t::Point(3.0, 2.0),
      new p2t::Point(1.000000001, 2.000000001),
      new p2t::Point(1.0, 1.999999999),
      new p2t::Point(0.0, 2.0)
    };
    f.steiner = { new p2t::Point(0.999999999, 0.5) };
    CloseConstraints(f, f.boundary);
  } else if (name == "hole") {
    // Square boundary with a centered square hole.
    f.boundary = {
      new p2t::Point(0.0, 0.0), new p2t::Point(3.0, 0.0),
      new p2t::Point(3.0, 3.0), new p2t::Point(0.0, 3.0)
    };
    f.holes.push_back({
      new p2t::Point(1.0, 1.0), new p2t::Point(2.0, 1.0),
      new p2t::Point(2.0, 2.0), new p2t::Point(1.0, 2.0)
    });
    CloseConstraints(f, f.boundary);
    CloseConstraints(f, f.holes[0]);
  } else {
    throw std::runtime_error("unknown fixture: " + name);
  }
  return f;
}

void DeleteFixture(Fixture& f)
{
  // A point may be shared across loops (touching holes): delete each once.
  std::set<p2t::Point*> owned;
  for (p2t::Point* p : f.boundary) owned.insert(p);
  for (auto& h : f.holes)
    for (p2t::Point* p : h) owned.insert(p);
  for (p2t::Point* p : f.steiner) owned.insert(p);
  for (p2t::Point* p : owned) delete p;
}

// ---------------------------------------------------------------------------
// Canonical serialization
// ---------------------------------------------------------------------------
struct Canonical {
  // point identity (pointer) -> canonical index
  std::map<const p2t::Point*, std::size_t> index;
};

std::string CanonicalEdge(std::size_t i, std::size_t j)
{
  if (i > j) std::swap(i, j);
  std::ostringstream os;
  os << i << "," << j;
  return os.str();
}

int Run(const std::string& fixture_name)
{
  Fixture f = MakeFixture(fixture_name);

  // Insertion order: boundary first (in given order), then holes (each in
  // given order), then Steiner points. This is the order in which the CDT
  // receives them.
  std::vector<p2t::Point*> insertion;
  for (p2t::Point* p : f.boundary) insertion.push_back(p);
  for (auto& h : f.holes)
    for (p2t::Point* p : h) insertion.push_back(p);
  for (p2t::Point* p : f.steiner) insertion.push_back(p);

  p2t::CDT cdt(f.boundary);
  for (auto& h : f.holes) cdt.AddHole(h);
  for (p2t::Point* p : f.steiner) cdt.AddPoint(p);
  cdt.Triangulate();

  std::vector<p2t::Triangle*> triangles = cdt.GetTriangles();

  // Assign canonical indices with the explicit rule:
  // quantized (qx,qy) ascending, insertion order as the tie breaker.
  std::map<const p2t::Point*, std::size_t> insertion_order;
  for (std::size_t i = 0; i < insertion.size(); ++i) {
    insertion_order[insertion[i]] = i;
  }
  struct Rec {
    const p2t::Point* p;
    PointKey key;
  };
  std::vector<Rec> recs;
  std::set<const p2t::Point*> seen;
  for (p2t::Triangle* t : triangles) {
    for (int k = 0; k < 3; ++k) {
      const p2t::Point* p = t->GetPoint(k);
      if (seen.insert(p).second) {
        Rec r;
        r.p = p;
        r.key = { Quantize(p->x), Quantize(p->y), insertion_order[p] };
        recs.push_back(r);
      }
    }
  }
  std::sort(recs.begin(), recs.end(), [](const Rec& a, const Rec& b) {
    return PointKeyLess()(a.key, b.key);
  });
  Canonical canon;
  for (std::size_t i = 0; i < recs.size(); ++i) {
    canon.index[recs[i].p] = i;
  }

  // Canonical triangles: sorted index triples.
  std::vector<std::array<std::size_t, 3>> ct;
  for (p2t::Triangle* t : triangles) {
    std::array<std::size_t, 3> tri = {
      canon.index[t->GetPoint(0)], canon.index[t->GetPoint(1)],
      canon.index[t->GetPoint(2)]
    };
    std::sort(tri.begin(), tri.end());
    ct.push_back(tri);
  }
  std::sort(ct.begin(), ct.end(),
            [](const std::array<std::size_t, 3>& a,
               const std::array<std::size_t, 3>& b) {
              for (int k = 0; k < 3; ++k)
                if (a[k] != b[k]) return a[k] < b[k];
              return false;
            });

  // Collect constrained edges actually present in the output.
  // constrained_edge[k] marks the edge OPPOSITE point k, i.e. between
  // points (k+1)%3 and (k+2)%3 (see Triangle::MarkConstrainedEdge).
  std::map<std::string, int> constrained_edges;
  for (p2t::Triangle* t : triangles) {
    for (int k = 0; k < 3; ++k) {
      if (t->constrained_edge[k]) {
        const p2t::Point* p1 = t->GetPoint((k + 1) % 3);
        const p2t::Point* p2 = t->GetPoint((k + 2) % 3);
        constrained_edges[CanonicalEdge(canon.index[p1], canon.index[p2])]++;
      }
    }
  }

  // Verify constraint integrity in both directions:
  //  (a) every requested constraint edge appears, marked constrained on at
  //      least one triangle;
  //  (b) every constrained edge in the mesh was actually requested (no
  //      spurious constrained edges introduced by the sweep).
  bool constraints_pass = true;
  std::vector<std::string> constraint_failures;
  std::set<std::string> expected_keys;
  for (auto& e : f.expected_constraints) {
    auto it1 = canon.index.find(e.first);
    auto it2 = canon.index.find(e.second);
    if (it1 == canon.index.end() || it2 == canon.index.end()) {
      constraints_pass = false;
      constraint_failures.push_back("missing-vertex");
      continue;
    }
    const std::string key = CanonicalEdge(it1->second, it2->second);
    expected_keys.insert(key);
    auto ce = constrained_edges.find(key);
    if (ce == constrained_edges.end() || ce->second < 1) {
      constraints_pass = false;
      constraint_failures.push_back("missing:" + key);
    }
  }
  for (const auto& kv : constrained_edges) {
    if (expected_keys.find(kv.first) == expected_keys.end()) {
      constraints_pass = false;
      constraint_failures.push_back("spurious:" + kv.first);
    }
  }

  // ---- Deterministic byte stream -------------------------------------------
  std::ostringstream out;
  out << "FIXTURE " << fixture_name << "\n";
  out << "QUANTIZATION_PRECISION " << DriverPrecision() << "\n";
  out << "VERTICES " << recs.size() << "\n";
  for (const Rec& r : recs) {
    out << canon.index[r.p] << " " << HexDouble(r.p->x) << " "
        << HexDouble(r.p->y) << "\n";
  }
  out << "TRIANGLES " << ct.size() << "\n";
  for (const auto& tri : ct) {
    out << tri[0] << " " << tri[1] << " " << tri[2] << "\n";
  }
  out << "CONSTRAINED_EDGES " << constrained_edges.size() << "\n";
  for (const auto& kv : constrained_edges) {
    out << kv.first << " flags=" << kv.second << "\n";
  }
  out << "EXPECTED_CONSTRAINTS " << f.expected_constraints.size() << "\n";
  out << "CONSTRAINTS " << (constraints_pass ? "PASS" : "FAIL");
  if (!constraint_failures.empty()) {
    out << " missing:";
    for (const auto& s : constraint_failures) out << " " << s;
  }
  out << "\n";
  out << "END\n";

  const std::string text = out.str();
  std::fwrite(text.data(), 1, text.size(), stdout);

  DeleteFixture(f);
  return constraints_pass ? 0 : 2;
}

// ---------------------------------------------------------------------------
// Raw serialization: emits GetTriangles() exactly as the library returns it
// - triangle order and per-triangle vertex slot order are NOT reordered;
// - each vertex is tagged with its driver-side insertion id, so even when two
//   distinct Point objects share coordinates, the chosen object identity is
//   visible in the byte stream.
// Used to prove the default library output is unchanged vertex-by-vertex.
// ---------------------------------------------------------------------------
int RunRaw(const std::string& fixture_name)
{
  Fixture f = MakeFixture(fixture_name);

  std::vector<p2t::Point*> insertion;
  for (p2t::Point* p : f.boundary) insertion.push_back(p);
  for (auto& h : f.holes)
    for (p2t::Point* p : h) insertion.push_back(p);
  for (p2t::Point* p : f.steiner) insertion.push_back(p);
  std::map<const p2t::Point*, std::size_t> insertion_id;
  for (std::size_t i = 0; i < insertion.size(); ++i) {
    insertion_id[insertion[i]] = i;
  }

  p2t::CDT cdt(f.boundary);
  for (auto& h : f.holes) cdt.AddHole(h);
  for (p2t::Point* p : f.steiner) cdt.AddPoint(p);
  cdt.Triangulate();
  std::vector<p2t::Triangle*> triangles = cdt.GetTriangles();

  std::ostringstream out;
  out << "RAW " << fixture_name << "\n";
  out << "TRIANGLES " << triangles.size() << "\n";
  for (std::size_t ti = 0; ti < triangles.size(); ++ti) {
    p2t::Triangle* t = triangles[ti];
    out << ti;
    for (int k = 0; k < 3; ++k) {
      const p2t::Point* p = t->GetPoint(k);
      out << " | id=" << insertion_id[p] << " " << HexDouble(p->x) << " "
          << HexDouble(p->y);
    }
    out << "\n";
  }
  out << "END\n";
  const std::string text = out.str();
  std::fwrite(text.data(), 1, text.size(), stdout);
  DeleteFixture(f);
  return 0;
}

} // namespace

int main(int argc, char** argv)
{
  if (argc != 3) {
    std::fprintf(stderr,
                 "usage: %s <canonical|raw> <fixture-name>\n", argv[0]);
    return 64;
  }
  try {
    const std::string mode = argv[1];
    if (mode == "canonical") return Run(argv[2]);
    if (mode == "raw") return RunRaw(argv[2]);
    std::fprintf(stderr, "unknown mode: %s\n", mode.c_str());
    return 64;
  } catch (const std::exception& e) {
    std::fprintf(stderr, "error: %s\n", e.what());
    return 1;
  }
}
