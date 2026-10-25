#!/bin/sh
set -eu

package_dir=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
compiler=${FC:-gfortran}
build_dir=$(mktemp -d "${TMPDIR:-/tmp}/fastqsl2-build.XXXXXX")

cleanup() {
    find "$build_dir" -type f -delete
    rmdir "$build_dir"
}
trap cleanup EXIT HUP INT TERM

mkdir -p "$package_dir/bin"
cd "$build_dir"
"$compiler" -O3 -fopenmp -march=native \
    -J "$build_dir" \
    -I "$package_dir" \
    -o "$package_dir/bin/fastqsl.x" \
    "$package_dir/fastqsl.f90"

echo "Built $package_dir/bin/fastqsl.x"
