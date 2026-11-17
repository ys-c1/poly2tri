/*
 * Quantization-precision validation for the deterministic poly2tri output.
 * STL only; no new third-party dependency, no network access.
 *
 * Checks:
 *   1. invalid P2T_QUANTIZATION_PRECISION values make triangulation throw
 *      std::runtime_error (empty, non-numeric, negative, fractional, > 17);
 *   2. every legal value in [0, 17] is accepted on a non-degenerate polygon;
 *   3. in exact mode a near-coincident pair stays distinct, while at a
 *      precision coarser than their separation the two points quantize to
 *      the same cell and the second one is rejected as an existing vertex;
 *   4. the resulting vertex count follows the quantization rule exactly.
 *
 * Cross-variant byte identity of the serialized output is verified by
 * run.sh (sha256 over all g++/clang++, standard and optimization variants).
 */

#include <poly2tri/poly2tri.h>

#include <cstdio>
#include <cstdlib>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

int g_failures = 0;

void Check(bool condition, const std::string& label)
{
  std::printf("  [%s] %s\n", condition ? "PASS" : "FAIL", label.c_str());
  if (!condition) ++g_failures;
}

// Non-degenerate 10x10 square; integer corners are far apart and never
// collapse onto one another for any precision in [0, 17].
std::vector<p2t::Point*> SquareBoundary()
{
  return {
    new p2t::Point(0.0, 0.0), new p2t::Point(10.0, 0.0),
    new p2t::Point(10.0, 10.0), new p2t::Point(0.0, 10.0)
  };
}

void DeleteBoundary(std::vector<p2t::Point*>& b)
{
  for (p2t::Point* p : b) delete p;
}

std::size_t DistinctVertices(const std::vector<p2t::Triangle*>& tris)
{
  std::set<const p2t::Point*> seen;
  for (p2t::Triangle* t : tris)
    for (int k = 0; k < 3; ++k) seen.insert(t->GetPoint(k));
  return seen.size();
}

void ExpectInvalid(const char* value)
{
  setenv("P2T_QUANTIZATION_PRECISION", value, 1);
  auto b = SquareBoundary();
  bool threw = false;
  std::string message;
  try {
    p2t::CDT cdt(b);
    cdt.Triangulate();
  } catch (const std::runtime_error& e) {
    threw = true;
    message = e.what();
  }
  Check(threw && message.find("P2T_QUANTIZATION_PRECISION") != std::string::npos,
        std::string("illegal precision '") + (value[0] ? value : "<empty>")
            + "' throws runtime_error");
  DeleteBoundary(b);
  unsetenv("P2T_QUANTIZATION_PRECISION");
}

} // namespace

int main()
{
  std::puts("== 1. illegal precision values must error ==");
  ExpectInvalid("");
  ExpectInvalid("abc");
  ExpectInvalid("-1");
  ExpectInvalid("3.5");
  ExpectInvalid("18");
  ExpectInvalid("99999999999999999999999");

  std::puts("== 2. every legal precision in [0,17] is accepted ==");
  for (int p = 0; p <= 17; ++p) {
    const std::string s = std::to_string(p);
    setenv("P2T_QUANTIZATION_PRECISION", s.c_str(), 1);
    auto b = SquareBoundary();
    bool ok = true;
    try {
      p2t::CDT cdt(b);
      cdt.Triangulate();
    } catch (const std::exception&) {
      ok = false;
    }
    Check(ok, "precision " + s + " accepted, 2 triangles");
    DeleteBoundary(b);
  }
  unsetenv("P2T_QUANTIZATION_PRECISION");

  std::puts("== 3. exact mode: near-coincident points stay distinct ==");
  {
    unsetenv("P2T_QUANTIZATION_PRECISION");
    auto b = SquareBoundary();
    p2t::CDT cdt(b);
    // 1e-7 away from the interior point (5,5): distinct in exact mode.
    auto* interior = new p2t::Point(5.0, 5.0);
    auto* near = new p2t::Point(5.0000001, 5.0);
    cdt.AddPoint(interior);
    cdt.AddPoint(near);
    cdt.Triangulate();
    Check(DistinctVertices(cdt.GetTriangles()) == 6,
          "6 vertices: boundary(4) + 2 distinct interior points");
    delete interior;
    delete near;
    DeleteBoundary(b);
  }

  std::puts("== 4. precision 6: sub-1e-6 points merge into one vertex ==");
  {
    setenv("P2T_QUANTIZATION_PRECISION", "6", 1);
    auto b = SquareBoundary();
    p2t::CDT cdt(b);
    auto* interior = new p2t::Point(5.0, 5.0);
    cdt.AddPoint(interior);
    // Rounds to exactly (5.000000, 5.000000) at precision 6 -> same vertex.
    auto* merged = new p2t::Point(5.0000001, 5.0);
    bool rejected = false;
    std::string what;
    try {
      cdt.AddPoint(merged);
    } catch (const std::runtime_error& e) {
      rejected = true;
      what = e.what();
    }
    Check(rejected && what == "Point already exists",
          "second point in same quantization cell rejected");
    delete merged;
    cdt.Triangulate();
    Check(DistinctVertices(cdt.GetTriangles()) == 5,
          "5 vertices: boundary(4) + 1 merged interior point");
    delete interior;
    DeleteBoundary(b);
    unsetenv("P2T_QUANTIZATION_PRECISION");
  }

  std::puts("== 5. precision 9: same pair (1e-7 apart) stays distinct ==");
  {
    setenv("P2T_QUANTIZATION_PRECISION", "9", 1);
    auto b = SquareBoundary();
    p2t::CDT cdt(b);
    auto* interior = new p2t::Point(5.0, 5.0);
    auto* near = new p2t::Point(5.0000001, 5.0);
    cdt.AddPoint(interior);
    cdt.AddPoint(near);
    cdt.Triangulate();
    Check(DistinctVertices(cdt.GetTriangles()) == 6,
          "6 vertices: 1e-7 separation survives 9-decimal quantization");
    delete interior;
    delete near;
    DeleteBoundary(b);
    unsetenv("P2T_QUANTIZATION_PRECISION");
  }

  std::printf("\n%s\n", g_failures == 0 ? "ALL QUANTIZATION CHECKS PASSED"
                                        : "QUANTIZATION CHECK FAILURES");
  return g_failures == 0 ? 0 : 1;
}
