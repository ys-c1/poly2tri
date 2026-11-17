#!/usr/bin/env bash
# Build and run the standalone resampling checks/benchmark.
# Uses only the C++ standard library and poly2tri sources in this repository.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
BUILD_DIR="$(mktemp -d)"
trap 'rm -rf "${BUILD_DIR}"' EXIT

g++ -std=c++14 -O2 -Wall -Wextra -Werror -DP2T_STATIC_EXPORTS \
  -I"${REPO_ROOT}" \
  "${SCRIPT_DIR}/resample_bench.cpp" \
  "${REPO_ROOT}/poly2tri/common/shapes.cc" \
  "${REPO_ROOT}/poly2tri/sweep/advancing_front.cc" \
  "${REPO_ROOT}/poly2tri/sweep/cdt.cc" \
  "${REPO_ROOT}/poly2tri/sweep/sweep.cc" \
  "${REPO_ROOT}/poly2tri/sweep/sweep_context.cc" \
  -o "${BUILD_DIR}/resample_bench"

"${BUILD_DIR}/resample_bench"
