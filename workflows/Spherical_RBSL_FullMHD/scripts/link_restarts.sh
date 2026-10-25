#!/bin/bash
set -euo pipefail

ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
PFSS_STAGE=${PFSS_STAGE:-01_pfss}
if [[ "${PFSS_STAGE}" == "01_pfss_global_ct" ]]; then
  DEFAULT_PFSS_DAT=${ROOT}/01_pfss_global_ct/data/pfss_global_ct_0000.dat
else
  DEFAULT_PFSS_DAT=${ROOT}/01_pfss/data/pfss_0000.dat
fi
PFSS_DAT=${1:-${DEFAULT_PFSS_DAT}}
RELAX_DAT=${2:-${ROOT}/02_relax/data/relax_0020.dat}
INSERT_DAT=${3:-${ROOT}/03_insert/data/rbsl_insert_load_0011.dat}
FULLMHD_SOURCE=${FULLMHD_SOURCE:-insert}

ln -sfn "${PFSS_DAT}" "${ROOT}/02_relax/restart.dat"
ln -sfn "${RELAX_DAT}" "${ROOT}/03_insert/restart.dat"

case "${FULLMHD_SOURCE}" in
  insert) FULLMHD_DAT=${INSERT_DAT} ;;
  relax)
    FULLMHD_DAT=${RELAX_DAT}
    echo "warning: direct relax -> fullmhd requires PFSS without RBSL Br subtraction" >&2
    ;;
  *)
    echo "FULLMHD_SOURCE must be insert or relax: ${FULLMHD_SOURCE}" >&2
    exit 2
    ;;
esac
ln -sfn "${FULLMHD_DAT}" "${ROOT}/04_fullmhd/restart.dat"

echo "02_relax/restart.dat -> ${PFSS_DAT}"
echo "03_insert/restart.dat -> ${RELAX_DAT}"
echo "04_fullmhd/restart.dat -> ${FULLMHD_DAT} (${FULLMHD_SOURCE})"
