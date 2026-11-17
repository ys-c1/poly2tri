/*
 * Demonstrates the non-determinism hazard of the original point comparator
 * p2t::cmp: it defines a strict weak ordering where points with equal
 * coordinates are EQUIVALENT, so std::sort is permitted to place them in any
 * relative order. Different STL implementations (or the same implementation
 * on different input permutations / sizes) can therefore sweep equivalent
 * points in different orders, which changes the resulting triangle sequence.
 *
 * This file is test scaffolding only; it does not touch the library.
 */
#include <poly2tri/common/shapes.h>
#include <algorithm>
#include <cstdio>
#include <vector>

static void Show(const char* label, const std::vector<p2t::Point*>& v)
{
  std::printf("%s:", label);
  for (auto* p : v) std::printf(" (%g,%g)", p->x, p->y);
  std::printf("\n");
}

int main()
{
  // Five distinct Point objects; P2..P4 share the same y (and x) -> they are
  // equivalent under cmp, though distinguishable by insertion identity.
  p2t::Point p0(0.0, 0.0), p1(1.0, 0.0), p2(2.0, 1.0), p3(2.0, 1.0),
      p4(2.0, 1.0), p5(3.0, 2.0);
  std::vector<p2t::Point*> base{ &p0, &p1, &p2, &p3, &p4, &p5 };

  std::vector<p2t::Point*> a = base;
  std::sort(a.begin(), a.end(), p2t::cmp);
  Show("std::sort(cmp) input order A      ", a);

  std::vector<p2t::Point*> b = base;
  std::reverse(b.begin(), b.end());
  std::sort(b.begin(), b.end(), p2t::cmp);
  Show("std::sort(cmp) input order B (rev) ", b);

  // Both outputs satisfy cmp. Compare object identity order at the tie group.
  const bool same = (a[2] == b[2] && a[3] == b[3] && a[4] == b[4]);
  std::printf("equivalent-point object order stable across permutations: %s\n",
              same ? "yes" : "NO -> sweep order depends on sort/input");
  return same ? 0 : 1;
}
