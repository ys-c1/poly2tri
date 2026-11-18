#!/bin/sh
# 编译并运行自交修复验证工具。不改动任何构建脚本,直接用 g++ 编译仓库现有源码。
set -e
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
BIN="$(mktemp -d)/self_intersect_fix"
g++ -std=c++14 -O2 -Wall -Wextra -I"$ROOT" \
  "$ROOT/tools/self_intersect_fix.cc" \
  "$ROOT/poly2tri/common/shapes.cc" \
  "$ROOT/poly2tri/sweep/cdt.cc" \
  "$ROOT/poly2tri/sweep/sweep.cc" \
  "$ROOT/poly2tri/sweep/sweep_context.cc" \
  "$ROOT/poly2tri/sweep/advancing_front.cc" \
  -o "$BIN"
"$BIN" "$@"
