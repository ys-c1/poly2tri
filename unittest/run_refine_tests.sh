#!/usr/bin/env bash
#
# Build and run the Delaunay refinement acceptance tests.
#
# Uses only the repository's own poly2tri sources and the C++ standard
# library (no third-party test framework, no network access). The existing
# CMake/meson build files are untouched.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BUILD_DIR="${ROOT}/build"
mkdir -p "${BUILD_DIR}"

SOURCES=(
  "${ROOT}/poly2tri/common/shapes.cc"
  "${ROOT}/poly2tri/sweep/advancing_front.cc"
  "${ROOT}/poly2tri/sweep/cdt.cc"
  "${ROOT}/poly2tri/sweep/sweep.cc"
  "${ROOT}/poly2tri/sweep/sweep_context.cc"
  "${ROOT}/unittest/delaunay_refine_test.cpp"
)

g++ -std=c++14 -O2 -DP2T_STATIC_EXPORTS -I"${ROOT}" \
    "${SOURCES[@]}" -o "${BUILD_DIR}/delaunay_refine_test"

"${BUILD_DIR}/delaunay_refine_test"
