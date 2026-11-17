#!/bin/sh
# Build and run the standalone LSCM UV parameterization test.
# Uses only the repository sources and the system compiler (no build-system
# changes, no third-party dependencies, no network access).
set -e

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
OUT_DIR="$ROOT/artifacts"
mkdir -p "$OUT_DIR"

g++ -std=c++14 -O2 -Wall -Wextra -DP2T_STATIC_EXPORTS -I"$ROOT" \
    "$ROOT"/poly2tri/common/shapes.cc \
    "$ROOT"/poly2tri/sweep/advancing_front.cc \
    "$ROOT"/poly2tri/sweep/cdt.cc \
    "$ROOT"/poly2tri/sweep/sweep.cc \
    "$ROOT"/poly2tri/sweep/sweep_context.cc \
    "$ROOT"/scripts/lscm_uv_test.cpp \
    -o "$OUT_DIR/lscm_uv_test"

"$OUT_DIR/lscm_uv_test"
