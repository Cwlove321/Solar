#!/usr/bin/env bash

# Run MPI-AMRVAC's native optically thin AIA synthesis on a spherical,
# stretched-grid snapshot. The two built-in views match the rbsl14 DDA views.

set -euo pipefail

usage() {
  cat <<'EOF'
Usage:
  synthesize_aia_amrvac_native_spherical_stretched.sh \
    BASE_PAR INPUT_DAT AMRVAC_EXE OUTPUT_DIR [WAVELENGTH]

Environment variables:
  MPI_RANKS                    MPI ranks used by each view (default: 64)
  MPIRUN                       MPI launcher (default: mpirun)
  PYTHON                       Python used to render VTI files
  INSTRUMENT_RESOLUTION_FACTOR AMRVAC resolution factor (default: 0.25)
  R_OPT_THICK                  opaque solar radius in Rsun (default: 1.0)
  REUSE_EXISTING_VTI           reuse matching non-empty VTI files (default: false)

The script runs two views sequentially:
  top   radial inward above the rbsl14 target centre
  side  south-to-north tangent through the rbsl14 flux rope
EOF
}

if [[ $# -lt 4 || $# -gt 5 ]]; then
  usage >&2
  exit 2
fi

BASE_PAR=$1
INPUT_DAT=$2
AMRVAC_EXE=$3
OUTPUT_DIR=$4
WAVELENGTH=${5:-193}

MPI_RANKS=${MPI_RANKS:-64}
MPIRUN=${MPIRUN:-mpirun}
PYTHON=${PYTHON:-/share/home/guoyang1/apps/miniconda3/envs/zcc/bin/python}
INSTRUMENT_RESOLUTION_FACTOR=${INSTRUMENT_RESOLUTION_FACTOR:-0.25}
R_OPT_THICK=${R_OPT_THICK:-1.0}
REUSE_EXISTING_VTI=${REUSE_EXISTING_VTI:-false}

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
RENDER_SCRIPT="$SCRIPT_DIR/render_amrvac_native_aia.py"

for path in "$BASE_PAR" "$INPUT_DAT" "$AMRVAC_EXE" "$RENDER_SCRIPT"; do
  if [[ ! -f "$path" ]]; then
    echo "Required file not found: $path" >&2
    exit 1
  fi
done

case "$WAVELENGTH" in
  94|131|171|193|211|304|335) ;;
  *)
    echo "Unsupported native AIA wavelength: $WAVELENGTH" >&2
    exit 2
    ;;
esac

mkdir -p "$OUTPUT_DIR"/{data,params,vti,png,logs}
OUTPUT_DIR=$(cd "$OUTPUT_DIR" && pwd)
BASE_PAR=$(cd "$(dirname "$BASE_PAR")" && pwd)/$(basename "$BASE_PAR")
INPUT_DAT=$(cd "$(dirname "$INPUT_DAT")" && pwd)/$(basename "$INPUT_DAT")
AMRVAC_EXE=$(cd "$(dirname "$AMRVAC_EXE")" && pwd)/$(basename "$AMRVAC_EXE")

FRAME=$(basename "$INPUT_DAT" .dat)

# AMRVAC stores command-line arguments in a short fixed-length Fortran string.
# Short links prevent long cluster paths from being silently truncated.
ln -sfn "$BASE_PAR" "$OUTPUT_DIR/base.par"
ln -sfn "$INPUT_DAT" "$OUTPUT_DIR/${FRAME}.dat"
BASE_DIR=$(dirname "$BASE_PAR")
EXE_DIR=$(dirname "$AMRVAC_EXE")
for runtime_input in path.dat parameters.txt; do
  for runtime_dir in "$BASE_DIR" "$EXE_DIR"; do
    if [[ -f "$runtime_dir/$runtime_input" ]]; then
      ln -sfn "$runtime_dir/$runtime_input" "$OUTPUT_DIR/$runtime_input"
      break
    fi
  done
done

write_override() {
  local view=$1
  local theta=$2
  local phi=$3
  local par="$OUTPUT_DIR/params/aia${WAVELENGTH}_${FRAME}_${view}.par"
  local prefix="vti/aia${WAVELENGTH}_${FRAME}_${view}_"

  cat > "$par" <<EOF
&filelist
  autoconvert = .false.
  convert_type = 'EIvtiCCmpi'
/

&Mag_list
  do_rbsl = .false.
  liakh_mass_loading = .false.
/

&emissionlist
  filename_euv = '${prefix}'
  wavelength = ${WAVELENGTH}
  LOS_theta = ${theta}d0
  LOS_phi = ${phi}d0
  image_rotate = 0.0d0
  dat_resolution = .false.
  instrument_resolution_factor = ${INSTRUMENT_RESOLUTION_FACTOR}d0
  activate_unit_arcsec = .true.
  R_opt_thick = ${R_OPT_THICK}d0
/
EOF

  printf '%s\n' "params/$(basename "$par")"
}

run_view() {
  local view=$1
  local theta=$2
  local phi=$3
  local par
  local prefix="$OUTPUT_DIR/vti/aia${WAVELENGTH}_${FRAME}_${view}_"
  local log="$OUTPUT_DIR/logs/aia${WAVELENGTH}_${FRAME}_${view}.out"
  local err="$OUTPUT_DIR/logs/aia${WAVELENGTH}_${FRAME}_${view}.err"
  local vti

  par=$(write_override "$view" "$theta" "$phi")
  if [[ "$REUSE_EXISTING_VTI" == "true" ]]; then
    vti=$(find "$OUTPUT_DIR/vti" -maxdepth 1 -type f \
      -name "$(basename "$prefix")*.vti" -print | sort | tail -1)
    if [[ -n "$vti" && -s "$vti" ]]; then
      printf '%s\n' "$vti"
      return
    fi
  fi
  echo "Running native AIA ${WAVELENGTH} view=${view}, theta=${theta}, phi=${phi}"
  (
    cd "$OUTPUT_DIR"
    "$MPIRUN" -np "$MPI_RANKS" "$AMRVAC_EXE" \
      -i base.par "$par" -if "${FRAME}.dat" -convert
  ) >"$log" 2>"$err"

  vti=$(find "$OUTPUT_DIR/vti" -maxdepth 1 -type f \
    -name "$(basename "$prefix")*.vti" -print | sort | tail -1)
  if [[ -z "$vti" || ! -s "$vti" ]]; then
    echo "AMRVAC did not create a non-empty VTI for view=$view" >&2
    exit 1
  fi
  if ! grep -q "Systhesizing EUV image" "$log"; then
    echo "AMRVAC log does not contain the EUV synthesis marker: $log" >&2
    exit 1
  fi
  printf '%s\n' "$vti"
}

# AMRVAC defines vec_LOS as minus the spherical unit vector specified here.
# These values reproduce the directions stored by the rbsl14 DDA workflow.
TOP_VTI=$(run_view top 72.0 195.0 | tail -1)
SIDE_VTI=$(run_view side 159.92562506233102 204.95737169575634 | tail -1)

"$PYTHON" "$RENDER_SCRIPT" \
  --top "$TOP_VTI" \
  --side "$SIDE_VTI" \
  --output-dir "$OUTPUT_DIR/png" \
  --wave "$WAVELENGTH" \
  --frame "$FRAME"

printf 'Native AMRVAC AIA synthesis completed: %s\n' "$OUTPUT_DIR"
