#!/usr/bin/env bash
# Build and run the standalone topology inspection tests.
# Uses only the system C++ toolchain and the C++ standard library;
# the project's CMake/meson build files are untouched.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
BUILD_DIR="${REPO_ROOT}/build/topology-test"

mkdir -p "${BUILD_DIR}"

CXX="${CXX:-g++}"
"${CXX}" -std=c++14 -O2 -Wall -Wextra -pedantic \
    -I "${REPO_ROOT}" \
    "${SCRIPT_DIR}/topology_test.cpp" \
    -o "${BUILD_DIR}/topology_test"

exec "${BUILD_DIR}/topology_test"
