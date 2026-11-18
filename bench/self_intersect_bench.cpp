// Differential test + degenerate-case test + benchmark for the sweep-line
// self-intersection detection vs. the O(n^2) brute-force oracle.
// Standalone: no third-party dependencies. Build:
//   g++ -O2 -std=c++14 -I. bench/self_intersect_bench.cpp
//       poly2tri/common/self_intersect.cc poly2tri/common/shapes.cc -o si_bench
// Run with argument "quick" to skip the 50k-edge benchmark.

#include "poly2tri/common/self_intersect.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <random>
#include <string>
#include <vector>

using p2t::Point;
using p2t::SelfIntersection;

namespace {

const double kEps = 1e-9;

std::vector<Point> P(std::initializer_list<std::pair<double, double>> v)
{
  std::vector<Point> r;
  for (auto& p : v) r.push_back(Point(p.first, p.second));
  return r;
}

bool RecMatch(const SelfIntersection& a, const SelfIntersection& b)
{
  return a.edge_a == b.edge_a && a.edge_b == b.edge_b &&
         a.overlap == b.overlap && std::fabs(a.x - b.x) <= kEps &&
         std::fabs(a.y - b.y) <= kEps;
}

// Multiset containment: every record in `sub` has a distinct match in `sup`.
bool IsSubset(const std::vector<SelfIntersection>& sub,
              const std::vector<SelfIntersection>& sup)
{
  std::vector<bool> used(sup.size(), false);
  for (const auto& a : sub) {
    bool found = false;
    for (size_t i = 0; i < sup.size(); ++i) {
      if (!used[i] && RecMatch(a, sup[i])) {
        used[i] = true;
        found = true;
        break;
      }
    }
    if (!found) return false;
  }
  return true;
}

bool SameSet(const std::vector<SelfIntersection>& a,
             const std::vector<SelfIntersection>& b)
{
  return a.size() == b.size() && IsSubset(a, b) && IsSubset(b, a);
}

void Dump(const char* tag, const std::vector<SelfIntersection>& v)
{
  std::printf("  %s (%zu):", tag, v.size());
  for (const auto& r : v)
    std::printf(" [(%g,%g) e%zu-e%zu%s]", r.x, r.y, r.edge_a, r.edge_b,
                r.overlap ? " ov" : "");
  std::printf("\n");
}

// Records in a that have no match in b.
std::vector<SelfIntersection> SetDiff(const std::vector<SelfIntersection>& a,
                                      const std::vector<SelfIntersection>& b)
{
  std::vector<bool> used(b.size(), false);
  std::vector<SelfIntersection> r;
  for (const auto& x : a) {
    bool found = false;
    for (size_t i = 0; i < b.size(); ++i) {
      if (!used[i] && RecMatch(x, b[i])) {
        used[i] = true;
        found = true;
        break;
      }
    }
    if (!found) r.push_back(x);
  }
  return r;
}

int g_failures = 0;

// Run both implementations on one case and verify the acceptance rule:
// identical boolean answer; if the oracle finds <= 100 records the sweep
// must return exactly the same set, otherwise exactly 100 records that are
// a subset of the oracle's set (plus the truncated flag).
void DiffOne(const char* name, const std::vector<Point>& pl, bool closed,
             bool verbose_on_fail = true)
{
  std::vector<SelfIntersection> brute, sweep;
  bool tr_brute = false, tr_sweep = false;
  const bool found_brute = p2t::FindSelfIntersectionsBruteForce(
      pl, closed, SIZE_MAX, brute, &tr_brute);
  const bool found_sweep =
      p2t::FindSelfIntersections(pl, closed, 100, sweep, &tr_sweep);

  bool ok = true;
  std::string why;
  if (found_brute != found_sweep) {
    ok = false;
    why = "boolean mismatch";
  } else if (brute.size() <= 100) {
    if (tr_sweep) {
      ok = false;
      why = "sweep truncated although <=100 exist";
    } else if (!SameSet(brute, sweep)) {
      ok = false;
      why = "record sets differ";
    }
  } else {
    if (!tr_sweep || sweep.size() != 100) {
      ok = false;
      why = "sweep did not truncate at 100";
    } else if (!IsSubset(sweep, brute)) {
      ok = false;
      why = "truncated sweep result not a subset of oracle";
    }
  }
  if (!ok) {
    ++g_failures;
    if (verbose_on_fail) {
      std::printf("FAIL %s: %s (brute=%zu sweep=%zu)\n", name, why.c_str(),
                  brute.size(), sweep.size());
      Dump("brute-only", SetDiff(brute, sweep));
      Dump("sweep-only", SetDiff(sweep, brute));
    }
  }
}

// Crafted degenerate case with a hand-computed expected result.
void Expect(const char* name, const std::vector<Point>& pl, bool closed,
            const std::vector<SelfIntersection>& expected)
{
  std::vector<SelfIntersection> sweep, brute;
  p2t::FindSelfIntersections(pl, closed, 1000, sweep);
  p2t::FindSelfIntersectionsBruteForce(pl, closed, 1000, brute);
  const bool ok_sweep = SameSet(sweep, expected);
  const bool ok_brute = SameSet(brute, expected);
  if (!ok_sweep || !ok_brute) {
    ++g_failures;
    std::printf("FAIL %s: %s\n", name,
                !ok_sweep ? "sweep != expected" : "brute != expected");
    Dump("expected", expected);
    Dump("sweep   ", sweep);
    Dump("brute   ", brute);
  } else {
    std::printf("ok   %-28s -> %zu intersection(s)\n", name, sweep.size());
  }
}

SelfIntersection R(double x, double y, size_t a, size_t b, bool ov = false)
{
  SelfIntersection r;
  r.x = x;
  r.y = y;
  r.edge_a = a;
  r.edge_b = b;
  r.overlap = ov;
  return r;
}

// ---------------------------------------------------------------------------
// Degenerate cases
// ---------------------------------------------------------------------------

void DegenerateTests()
{
  std::printf("--- degenerate cases ---\n");

  // Horizontal edges crossed by verticals (open polyline).
  Expect("horizontal+vertical crossings",
         P({{0, 0}, {4, 0}, {4, 1}, {1, 1}, {1, -1}, {2, -1}, {2, 2}}),
         false,
         {R(1, 0, 0, 3), R(2, 0, 0, 5), R(2, 1, 2, 5)});

  // Collinear overlap: the polyline walks back over an earlier horizontal.
  Expect("collinear overlap (backtrack)",
         P({{0, 0}, {5, 0}, {5, 1}, {2, 1}, {2, 0}, {7, 0}}), false,
         {R(2, 0, 0, 3), R(2, 0, 0, 4, true), R(5, 0, 1, 4)});

  // Adjacent collinear continuation must NOT count (open and closed).
  Expect("adjacent collinear (open)", P({{0, 0}, {2, 0}, {4, 0}}), false, {});
  Expect("adjacent collinear (closed)", P({{0, 0}, {2, 0}, {4, 0}}), true,
         {});

  // Collinear edges sharing exactly one endpoint (non-adjacent).
  Expect("collinear single-point touch",
         P({{0, 0}, {4, 0}, {4, 2}, {6, 2}, {6, 0}, {4, 0}}), false,
         {R(4, 0, 0, 4), R(4, 0, 1, 4)});

  // Same-line horizontals, disjoint: no intersection.
  Expect("collinear disjoint horizontals",
         P({{0, 0}, {2, 0}, {2, 1}, {5, 1}, {5, 0}, {7, 0}}), false, {});

  // Endpoint touching: closed figure-eight revisiting vertex (0,0).
  Expect("endpoint touch (figure-8, closed)",
         P({{0, 0}, {2, 0}, {2, 2}, {0, 0}, {-2, 2}, {-2, 0}}), true,
         {R(0, 0, 0, 2), R(0, 0, 0, 3), R(0, 0, 2, 5), R(0, 0, 3, 5)});

  // T-junction: endpoint of one edge in the interior of another.
  Expect("T-junction (endpoint on interior)",
         P({{0, 0}, {4, 0}, {4, 2}, {2, 2}, {2, 0}}), false,
         {R(2, 0, 0, 3)});

  // Convex square: clean, no intersections.
  Expect("convex square (closed)",
         P({{0, 0}, {4, 0}, {4, 4}, {0, 4}}), true, {});

  // Three edges concurrent at (2,2) plus an irrational crossing at (4/3,4/3).
  Expect("concurrent crossings",
         P({{0, 0}, {4, 4}, {4, 0}, {0, 4}, {2, 0}, {2, 4}}), false,
         {R(2, 2, 0, 2), R(4.0 / 3.0, 4.0 / 3.0, 0, 3), R(2, 2, 0, 4),
          R(2, 2, 2, 4)});

  // Closed 5-point star: 5 proper crossings.
  {
    std::vector<Point> star;
    for (int i = 0; i < 5; ++i) {
      const double a = M_PI / 2 + i * 4 * M_PI / 5;
      star.push_back(Point(std::cos(a), std::sin(a)));
    }
    std::vector<SelfIntersection> sw, bf;
    p2t::FindSelfIntersections(star, true, 100, sw);
    p2t::FindSelfIntersectionsBruteForce(star, true, 100, bf);
    if (sw.size() == 5 && SameSet(sw, bf)) {
      std::printf("ok   %-28s -> 5 intersections\n", "5-point star (closed)");
    } else {
      ++g_failures;
      std::printf("FAIL 5-point star: sweep=%zu brute=%zu\n", sw.size(),
                  bf.size());
      Dump("sweep", sw);
      Dump("brute", bf);
    }
  }
}

// ---------------------------------------------------------------------------
// Differential test: 2000 random polylines
// ---------------------------------------------------------------------------

void DifferentialTests()
{
  std::printf("--- differential test: 2000 random polylines ---\n");
  int counts[4] = {0, 0, 0, 0};
  for (int i = 0; i < 2000; ++i) {
    std::mt19937_64 rng(12345 + i);
    const int cat = i % 4;
    const bool closed = (cat % 2 == 1);
    const bool grid = (cat >= 2);
    const int n = 5 + static_cast<int>(rng() % 41); // 5..45 vertices
    std::vector<Point> pl;
    double x = 0, y = 0;
    if (grid) {
      x = static_cast<double>(rng() % 9);
      y = static_cast<double>(rng() % 9);
    }
    pl.push_back(Point(x, y));
    for (int j = 1; j < n; ++j) {
      if (grid) {
        // Integer grid 0..8: dense in collinearities, touches, revisits.
        x = static_cast<double>(rng() % 9);
        y = static_cast<double>(rng() % 9);
      } else {
        const double ang =
            (rng() % 1000000) / 1000000.0 * 2 * M_PI;
        const double len = 0.2 + (rng() % 1000) / 1000.0;
        x += len * std::cos(ang);
        y += len * std::sin(ang);
      }
      pl.push_back(Point(x, y));
    }
    char name[64];
    std::snprintf(name, sizeof(name), "case#%d cat%d %s", i, cat,
                  closed ? "closed" : "open");
    DiffOne(name, pl, closed);
    counts[cat]++;
  }
  std::printf("ran %d cases (float-open %d, float-closed %d, grid-open %d, "
              "grid-closed %d)\n",
              counts[0] + counts[1] + counts[2] + counts[3], counts[0],
              counts[1], counts[2], counts[3]);
}

// ---------------------------------------------------------------------------
// Benchmark
// ---------------------------------------------------------------------------

std::vector<Point> MakeSpiral(int n)
{
  // Archimedean spiral with angular step 0.05 rad: the per-turn radial gap
  // (~0.25) stays well above the chord sagitta (~r*3.1e-4, r <= 101), so the
  // polyline is genuinely self-intersection free (verified by the oracle).
  std::vector<Point> pl;
  for (int i = 0; i < n; ++i) {
    const double a = i * 0.05;
    const double r = 1.0 + 0.002 * i;
    pl.push_back(Point(r * std::cos(a), r * std::sin(a)));
  }
  return pl;
}

std::vector<Point> MakeRandomWalk(int n, uint64_t seed)
{
  std::mt19937_64 rng(seed);
  std::vector<Point> pl;
  double x = 0, y = 0;
  pl.push_back(Point(x, y));
  for (int i = 1; i < n; ++i) {
    const double a = (rng() % 1000000) / 1000000.0 * 2 * M_PI;
    x += std::cos(a);
    y += std::sin(a);
    pl.push_back(Point(x, y));
  }
  return pl;
}

template <typename F>
double TimeOnce(F&& f)
{
  const auto t0 = std::chrono::steady_clock::now();
  f();
  const auto t1 = std::chrono::steady_clock::now();
  return std::chrono::duration<double>(t1 - t0).count();
}

void BenchOne(const char* name, const std::vector<Point>& pl)
{
  std::printf("--- benchmark: %s (%zu edges) ---\n", name, pl.size() - 1);

  std::vector<SelfIntersection> brute;
  bool tr_brute = false;
  double t_brute = 0;
  bool found_brute = false;
  t_brute = TimeOnce([&] {
    found_brute = p2t::FindSelfIntersectionsBruteForce(pl, false, SIZE_MAX,
                                                       brute, &tr_brute);
  });

  std::vector<SelfIntersection> sweep;
  bool tr_sweep = false;
  bool found_sweep = false;
  double t_sweep = 1e300;
  for (int rep = 0; rep < 5; ++rep) {
    sweep.clear();
    t_sweep = std::min(t_sweep, TimeOnce([&] {
      found_sweep = p2t::FindSelfIntersections(pl, false, 100, sweep,
                                               &tr_sweep);
    }));
  }

  // Correctness at scale.
  bool ok = (found_brute == found_sweep);
  if (brute.size() <= 100) {
    ok = ok && SameSet(brute, sweep);
  } else {
    ok = ok && tr_sweep && sweep.size() == 100 && IsSubset(sweep, brute);
  }
  if (!ok) {
    ++g_failures;
    std::printf("FAIL benchmark consistency for %s\n", name);
  }

  std::printf("brute force (full O(n^2) scan): %8.3f s, %zu intersections\n",
              t_brute, brute.size());
  std::printf("sweep line  (first 100):        %8.3f s, %zu reported%s\n",
              t_sweep, sweep.size(), tr_sweep ? " (truncated)" : "");
  std::printf("speedup: %.1fx  [%s]\n", t_brute / t_sweep,
              ok ? "consistent" : "INCONSISTENT");
}

} // namespace

int main(int argc, char** argv)
{
  const bool quick = (argc > 1 && std::string(argv[1]) == "quick");
  DegenerateTests();
  DifferentialTests();
  if (!quick) {
    BenchOne("spiral, no self-intersections", MakeSpiral(50001));
    BenchOne("random walk, many self-intersections", MakeRandomWalk(50001, 7));
  }

  if (g_failures == 0) {
    std::printf("ALL TESTS PASSED\n");
    return 0;
  }
  std::printf("%d FAILURES\n", g_failures);
  return 1;
}
