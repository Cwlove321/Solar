#!/bin/bash
set -euo pipefail

ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
PFSS_STAGE=${PFSS_STAGE:-01_pfss}

FULLMHD_SOURCE=${FULLMHD_SOURCE:-insert}

for stage in "${PFSS_STAGE}" 02_relax 04_fullmhd; do
  [[ -x "${ROOT}/${stage}/amrvac" ]] || {
    echo "Missing executable: ${ROOT}/${stage}/amrvac" >&2
    exit 2
  }
done
if [[ "${FULLMHD_SOURCE}" == "insert" ]]; then
  [[ -x "${ROOT}/03_insert/amrvac" ]] || {
    echo "Missing executable: ${ROOT}/03_insert/amrvac" >&2
    exit 2
  }
elif [[ "${FULLMHD_SOURCE}" != "relax" ]]; then
  echo "FULLMHD_SOURCE must be insert or relax: ${FULLMHD_SOURCE}" >&2
  exit 2
fi
if [[ "${FULLMHD_SOURCE}" == "relax" ]]; then
  echo "warning: direct relax -> fullmhd requires PFSS without RBSL Br subtraction" >&2
fi

extract_job_id() {
  sed -n 's/.*Job <\([0-9][0-9]*\)>.*/\1/p'
}

PFSS_JOB=$(cd "${ROOT}/${PFSS_STAGE}" && bsub < job | extract_job_id)
RELAX_JOB=$(cd "${ROOT}/02_relax" && bsub -w "done(${PFSS_JOB})" < job | extract_job_id)
if [[ "${FULLMHD_SOURCE}" == "insert" ]]; then
  INSERT_JOB=$(cd "${ROOT}/03_insert" && bsub -w "done(${RELAX_JOB})" < job | extract_job_id)
  FULL_JOB=$(cd "${ROOT}/04_fullmhd" && bsub -w "done(${INSERT_JOB})" < job | extract_job_id)
else
  FULL_JOB=$(cd "${ROOT}/04_fullmhd" && bsub -w "done(${RELAX_JOB})" < job | extract_job_id)
fi

echo "PFSS:    ${PFSS_JOB} (${PFSS_STAGE})"
echo "Relax:   ${RELAX_JOB} depends on ${PFSS_JOB}"
if [[ "${FULLMHD_SOURCE}" == "insert" ]]; then
  echo "Insert:  ${INSERT_JOB} depends on ${RELAX_JOB}"
  echo "FullMHD: ${FULL_JOB} depends on ${INSERT_JOB}"
else
  echo "FullMHD: ${FULL_JOB} depends on ${RELAX_JOB} (no RBSL insertion)"
fi
