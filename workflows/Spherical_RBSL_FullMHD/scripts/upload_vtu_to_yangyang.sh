#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Usage:
  upload_vtu_to_yangyang.sh '<guoguo_source_data_glob>' [yangyang_target_dir]
  upload_vtu_to_yangyang.sh <guoguo_source_data_dir> <yangyang_target_dir> [file_glob]

Examples:
  upload_vtu_to_yangyang.sh \
    '/share/home/guoyang1/zhouchang/20241009/rbsl12/internal_e_16/data/*vtu'

  upload_vtu_to_yangyang.sh \
    '/share/home/guoyang1/zhouchang/20241009/rbsl12/internal_e_16/data/*vtu' \
    /home/guoyang/custom/internal_e_16

  upload_vtu_to_yangyang.sh \
    /share/home/guoyang1/zhouchang/20241009/theoretical_simulation/0beta2/data \
    /home/guoyang/zhouchang/20241009/theoretical_simulation/0beta2 \
    'eru*.vtu'

Environment:
  SOURCE_SSH=guoguo     SSH alias for the source server
  TARGET_SSH=yangyang   SSH alias for the ParaView server
  STREAM_COMPRESSION=none|gzip
                        Optional payload compression (default: none). Use gzip
                        for large ASCII VTU files; the target file is unchanged.

Behavior:
  - With a source glob only, maps /share/home/guoyang1/zhouchang/.../<run>/data
    to /home/guoyang/zhouchang/.../<run> automatically.
  - With a source glob and target directory, uses the supplied target directory.
  - Creates the target directory.
  - Streams files through the local machine without saving VTU files locally.
  - Writes target files as <name>.part first.
  - Compares SHA256 on source and target .part.
  - Moves .part to the final filename only after hashes match.
  - Skips existing target files when their SHA256 already matches.
EOF
}

if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
  usage
  exit 0
fi

if [[ $# -lt 1 || $# -gt 3 ]]; then
  usage
  exit 2
fi

source_ssh="${SOURCE_SSH:-guoguo}"
target_ssh="${TARGET_SSH:-yangyang}"
stream_compression="${STREAM_COMPRESSION:-none}"

case "$stream_compression" in
  none|gzip) ;;
  *)
    echo "Unsupported STREAM_COMPRESSION: ${stream_compression}" >&2
    echo "Expected: none or gzip" >&2
    exit 2
    ;;
esac

if [[ $# -eq 1 || ( $# -eq 2 && ( "$1" == *"*"* || "$1" == *"?"* || "$1" == *"["* || "$1" == *.vtu ) ) ]]; then
  source_spec="$1"
  source_dir="$(dirname "$source_spec")"
  file_glob="$(basename "$source_spec")"

  if [[ $# -eq 2 ]]; then
    target_dir="$2"
  else
    case "$source_dir" in
      /share/home/guoyang1/zhouchang/*/data)
        relative_dir="${source_dir#/share/home/guoyang1/zhouchang/}"
        relative_dir="${relative_dir%/data}"
        target_dir="/home/guoyang/zhouchang/${relative_dir}"
        ;;
      *)
        echo "Cannot infer target directory from: ${source_dir}" >&2
        echo "Expected: /share/home/guoyang1/zhouchang/.../<run>/data" >&2
        exit 2
        ;;
    esac
  fi
else
  source_dir="$1"
  target_dir="$2"
  file_glob="${3:-*.vtu}"
fi

list_file="$(mktemp -t vtu-upload-list.XXXXXX)"
cleanup() {
  rm -f "$list_file"
}
trap cleanup EXIT

echo "Source: ${source_ssh}:${source_dir}"
echo "Target: ${target_ssh}:${target_dir}"
echo "Pattern: ${file_glob}"
echo "Stream compression: ${stream_compression}"

ssh -n "$target_ssh" "mkdir -p '$target_dir' && ls -ld '$target_dir'"
ssh -n "$source_ssh" "cd '$source_dir' && ls $file_glob 2>/dev/null | sort" > "$list_file"

count="$(wc -l < "$list_file" | tr -d ' ')"
if [[ "$count" == "0" ]]; then
  echo "No matching VTU files found."
  exit 0
fi

echo "Found ${count} file(s)."

idx=0
while IFS= read -r name; do
  [[ -n "$name" ]] || continue
  idx=$((idx + 1))

  src="${source_dir}/${name}"
  dst="${target_dir}/${name}"

  source_hash="$(ssh -n "$source_ssh" "sha256sum '$src'" | awk '{print $1}')"

  if ssh -n "$target_ssh" "test -f '$dst'"; then
    target_hash="$(ssh -n "$target_ssh" "sha256sum '$dst'" | awk '{print $1}')"
    if [[ "$source_hash" == "$target_hash" ]]; then
      echo "[$idx/$count] skip ${name} hash=${source_hash}"
      continue
    fi
  fi

  echo "[$idx/$count] upload ${name}"
  if [[ "$stream_compression" == "gzip" ]]; then
    ssh -n "$source_ssh" "gzip -1 -c '$src'" | ssh "$target_ssh" "gzip -d > '$dst.part'"
  else
    ssh -n "$source_ssh" "cat '$src'" | ssh "$target_ssh" "cat > '$dst.part'"
  fi

  part_hash="$(ssh -n "$target_ssh" "sha256sum '$dst.part'" | awk '{print $1}')"
  if [[ "$source_hash" != "$part_hash" ]]; then
    echo "HASH_MISMATCH ${name} source=${source_hash} target_part=${part_hash}" >&2
    exit 1
  fi

  ssh -n "$target_ssh" "mv '$dst.part' '$dst'"
  echo "[$idx/$count] ok ${name} hash=${source_hash}"
done < "$list_file"

echo "Checking for .part residuals..."
ssh -n "$target_ssh" "find '$target_dir' -maxdepth 1 -name '*.part' -print"
echo "Done."
