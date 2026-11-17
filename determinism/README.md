# Deterministic poly2tri output

The same triangulation input could previously produce different mesh byte
streams on different machines / standard-library implementations, so exported
meshes could not be compared with a binary diff. This directory contains the
reproduction and the verification harness for the fix.

## Root causes (both inside the library)

1. `SweepContext::InitTriangulation` sorted points with `std::sort` and the
   equivalence-only comparator `p2t::cmp`. Points with equal coordinates are
   *equivalent* under that comparator, so the standard allows their relative
   order to differ between implementations. The sweep order then changes the
   produced triangle sequence. `sort_instability_demo.cpp` demonstrates the
   reordering directly (it exits non-zero when equivalent points move).
2. De-duplication used `std::unordered_set<Point*>`, whose layout depends on
   the hash implementation and allocator addresses (hence on the platform and
   on ASLR).

## The explicit rules (implemented in `poly2tri/sweep/sweep_context.cc`)

* **Ordering** — quantized coordinates first: `(y, x)` ascending.
* **Equal-coordinate tie break** — insertion order, made explicit with
  `std::stable_sort` instead of `std::sort`.
* **Equality** — in the default exact mode the historical pointer-identity
  semantics are kept; with a quantization grid configured, two points whose
  coordinates round to the same cell are the same vertex.
* The hash set is replaced by an insertion-order linear scan, so no data
  structure's internal traversal order can influence the result.

Exact-mode ordering uses an unsigned monotone key of the IEEE-754 bit pattern
(`OrderPreservingKey`), so floating-point order is preserved exactly without
any floating-point comparison in the comparator.

## Configurable quantization precision

Environment variable `P2T_QUANTIZATION_PRECISION` (no public signature
changes):

| value | meaning |
|-------|---------|
| unset / `0` | exact mode — default; output is byte-identical to the old library |
| `1`..`17` | round each coordinate to that many decimals before compare/merge |

Illegal values (empty, non-numeric, negative, fractional, or `> 17`) make the
triangulation throw `std::runtime_error`; `quantization_test.cpp` covers this.
Coarse precision can legitimately collapse distinct boundary vertices onto one
cell (a zero-length constraint edge); that geometry error is then reported
deterministically on every platform.

## Files

| file | purpose |
|------|---------|
| `deterministic_driver.cpp` | canonical + raw serialization of the four fixtures |
| `quantization_test.cpp` | illegal-precision errors and merge semantics |
| `sort_instability_demo.cpp` | reproduces the old comparator hazard |
| `run.sh` | one-shot build + full verification, STL/coreutils only |

Fixtures are generated in code (no data files): `regular_grid`,
`shared_point` (two holes touching at one vertex), `near_coincident`
(sub-epsilon-separated vertices), `hole` (polygon with an interior hole).

## Run

```sh
./run.sh
```

It builds 8 toolchain variants (g++/clang++ × C++14/17 × -O0/-O2), runs the
quantization tests on all of them, hashes every fixture at precisions
0/3/6/9, checks bidirectional constraint integrity, and stress-tests raw
output. The summary lands in `out/report.txt`, per-run hashes in
`out/*.sha256`.
