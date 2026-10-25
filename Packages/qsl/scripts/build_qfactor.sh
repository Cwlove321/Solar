#!/bin/sh
set -eu

package_dir=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
compiler=${FC:-gfortran}
build_dir=$(mktemp -d "${TMPDIR:-/tmp}/solar-qsl-build.XXXXXX")

cleanup() {
    find "$build_dir" -type f -delete
    rmdir "$build_dir"
}
trap cleanup EXIT HUP INT TERM

mkdir -p "$package_dir/src/qsl/bin"
cd "$build_dir"
"$compiler" -O3 -fopenmp -march=native \
    -J "$build_dir" \
    -I "$package_dir/fortran" \
    -o "$package_dir/src/qsl/bin/qfactor.x" \
    "$package_dir/fortran/qfactor.f90"

echo "Built $package_dir/src/qsl/bin/qfactor.x"
