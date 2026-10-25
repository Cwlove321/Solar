#!/bin/bash
set -euo pipefail

ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
: "${AMRVAC_DIR:?Set AMRVAC_DIR to the MPI-AMRVAC source directory}"
export GLOBAL_FMHD_LIBDIR=${GLOBAL_FMHD_LIBDIR:-/share/home/guoyang1/apps/global_fmhd_deps/chenpf_libs}
PFSS_STAGE=${PFSS_STAGE:-01_pfss}

for stage in "${PFSS_STAGE}" 02_relax 03_insert 04_fullmhd; do
  echo "==> building ${stage}"
  cd "${ROOT}/${stage}"
  "${AMRVAC_DIR}/setup.pl" -d=3
  make ARCH=intel -j "${MAKE_JOBS:-8}"
done
