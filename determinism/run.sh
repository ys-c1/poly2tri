#!/usr/bin/env bash
#
# Deterministic-output verification for poly2tri.
# Runs everything in one shot; uses only the repository's own sources and the
# system C++ toolchain (g++, clang++) plus coreutils (sha256sum). No network
# access and no third-party library are required.
#
# Outputs:
#   build/<variant>/driver            four-fixture canonical/raw driver
#   build/<variant>/quantization_test quantization validation
#   build/sortdemo                     old-comparator instability demo
#   out/<fixture>__<variant>__p<precision>.sha256
#   out/report.txt                     hash / triangle / constraint table
#
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "${HERE}/.." && pwd)"
cd "${HERE}"

SRCS="${ROOT}/poly2tri/common/shapes.cc \
      ${ROOT}/poly2tri/sweep/advancing_front.cc \
      ${ROOT}/poly2tri/sweep/cdt.cc \
      ${ROOT}/poly2tri/sweep/sweep.cc \
      ${ROOT}/poly2tri/sweep/sweep_context.cc"
FIXTURES=(regular_grid shared_point near_coincident hole)
VARIANTS=(g++-c++14-O0 g++-c++14-O2 g++-c++17-O0 g++-c++17-O2
          clang++-c++14-O0 clang++-c++14-O2 clang++-c++17-O0 clang++-c++17-O2)
PRECISIONS=(0 3 6 9)

mkdir -p build out
rm -f out/*.sha256 out/report.txt

echo "=== 1/5  comparator instability demo (old p2t::cmp hazard) ==="
g++ -std=c++17 -I"${ROOT}" sort_instability_demo.cpp \
    "${ROOT}/poly2tri/common/shapes.cc" -o build/sortdemo
set +e
build/sortdemo
demo_rc=$?
set -e
if [ "${demo_rc}" -eq 1 ]; then
  echo "    confirmed: equivalent points reorder under std::sort(cmp)"
else
  echo "    WARNING: demo did not trigger on this platform"
fi

echo
echo "=== 2/5  build drivers across ${#VARIANTS[@]} toolchain variants ==="
for variant in "${VARIANTS[@]}"; do
  compiler="${variant%%-*}"
  rest="${variant#*-}"
  standard="${rest%%-*}"
  opt="${variant##*-}"
  mkdir -p "build/${variant}"
  "${compiler}" "-std=${standard}" "-${opt}" -I"${ROOT}" \
      deterministic_driver.cpp ${SRCS} -o "build/${variant}/driver"
  "${compiler}" "-std=${standard}" "-${opt}" -I"${ROOT}" \
      quantization_test.cpp ${SRCS} -o "build/${variant}/quantization_test"
  echo "    built ${variant}"
done

echo
echo "=== 3/5  quantization validation (illegal values must error) ==="
for variant in "${VARIANTS[@]}"; do
  "build/${variant}/quantization_test" | sed "s/^/    [${variant}] /"
done

echo
echo "=== 4/5  fixture hashes, triangle counts, constraint integrity ==="
: > out/report.txt
for precision in "${PRECISIONS[@]}"; do
  export P2T_QUANTIZATION_PRECISION="${precision}"
  for fixture in "${FIXTURES[@]}"; do
    declare -A hashes=()
    declare -A statuses=()
    triangles=""
    constraints=""
    for variant in "${VARIANTS[@]}"; do
      set +e
      text="$("build/${variant}/driver" canonical "${fixture}" 2>out/last.err)"
      rc=$?
      set -e
      if [ "${rc}" -eq 2 ]; then
        statuses["${variant}"]="CONSTRAINT-FAIL"
      elif [ "${rc}" -ne 0 ]; then
        # Coarse quantization can collapse distinct boundary vertices onto
        # one cell (zero-length constraint edge): a geometry-level error that
        # is itself deterministic. Record it rather than mislabel it.
        statuses["${variant}"]="DEGENERATE:$(cat out/last.err)"
        hashes["${variant}"]="DEGENERATE"
        continue
      else
        statuses["${variant}"]="OK"
      fi
      hash="$(printf '%s' "${text}" | sha256sum | cut -d' ' -f1)"
      hashes["${variant}"]="${hash}"
      triangles="$(printf '%s' "${text}" | sed -n 's/^TRIANGLES //p')"
      constraints="$(printf '%s' "${text}" | sed -n 's/^CONSTRAINTS //p' | awk '{print $1}')"
      printf '%s  %s  %s  p%s\n' "${hash}" "${fixture}" "${variant}" "${precision}" \
          > "out/${fixture}__${variant}__p${precision}.sha256"
    done
    unique_hashes="$(printf '%s\n' "${hashes[@]}" | sort -u | wc -l)"
    unique_status="$(printf '%s\n' "${statuses[@]}" | sort -u | wc -l)"
    if [ "${unique_hashes}" -eq 1 ] && [ "${unique_status}" -eq 1 ]; then
      verdict="IDENTICAL"
    else
      verdict="MISMATCH"
    fi
    representative="$(printf '%s\n' "${statuses[@]}" | sort -u | head -1)"
    if [ "${representative}" = "OK" ]; then
      printf '%-16s p=%d  triangles=%-3s constraints=%-4s hashes=%s (%s)\n' \
          "${fixture}" "${precision}" "${triangles}" "${constraints}" \
          "${unique_hashes}" "${verdict}" | tee -a out/report.txt
    else
      printf '%-16s p=%d  %s across all %d variants (%s)\n' \
          "${fixture}" "${precision}" "${representative}" \
          "${#VARIANTS[@]}" "${verdict}" | tee -a out/report.txt
    fi
    [ "${verdict}" = "IDENTICAL" ] || { echo "DETERMINISM FAILURE"; exit 1; }
  done
done
unset P2T_QUANTIZATION_PRECISION
rm -f out/last.err

echo
echo "=== 5/5  raw (un-canonicalized) library output, exact mode ==="
for fixture in "${FIXTURES[@]}"; do
  declare -A rawhashes=()
  for variant in "${VARIANTS[@]}"; do
    rawhashes["${variant}"]="$("build/${variant}/driver" raw "${fixture}" \
        | sha256sum | cut -d' ' -f1)"
  done
  unique_hashes="$(printf '%s\n' "${rawhashes[@]}" | sort -u | wc -l)"
  [ "${unique_hashes}" -eq 1 ] && verdict="IDENTICAL" || verdict="MISMATCH"
  printf '%-16s raw hashes=%s (%s)\n' "${fixture}" "${unique_hashes}" "${verdict}" \
      | tee -a out/report.txt
done

echo
echo "All artifacts written to ${HERE}/out (report.txt + .sha256 files)."
