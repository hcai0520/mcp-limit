#!/usr/bin/env bash
set -euo pipefail

RESULTS_DIR=${1:-outputs/raw}
OUTPUT_FILE=${2:-outputs/dy_mcp_merged.root}

mkdir -p "$(dirname "$OUTPUT_FILE")"
mapfile -d '' files < <(find "$RESULTS_DIR" -type f -name '*.root' -print0)
if [[ ${#files[@]} -eq 0 ]]; then
    echo "No ROOT files found under $RESULTS_DIR" >&2
    exit 1
fi

hadd -f "$OUTPUT_FILE" "${files[@]}"
