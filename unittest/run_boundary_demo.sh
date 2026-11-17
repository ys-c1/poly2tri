#!/bin/sh
# Build and run the boundary extraction validation program without touching
# the CMake/meson build files. Standard library only; no new dependencies.
set -e

repo_root=$(cd "$(dirname "$0")/.." && pwd)
build_dir="$repo_root/artifacts/boundary-demo"
mkdir -p "$build_dir"

g++ -std=c++14 -O2 -Wall -Wextra -Wpedantic \
    -I"$repo_root" \
    "$repo_root/unittest/boundary_demo.cpp" \
    "$repo_root/poly2tri/common/shapes.cc" \
    "$repo_root/poly2tri/sweep/advancing_front.cc" \
    "$repo_root/poly2tri/sweep/cdt.cc" \
    "$repo_root/poly2tri/sweep/sweep.cc" \
    "$repo_root/poly2tri/sweep/sweep_context.cc" \
    -o "$build_dir/boundary_demo"

"$build_dir/boundary_demo"
