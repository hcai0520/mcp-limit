#!/usr/bin/env bash
set -euo pipefail

JOB_LIST=${JOB_LIST:-job_list.csv}
JOB_ID=${1:-${SLURM_ARRAY_TASK_ID:-}}
EXEC=${EXEC:-./build-gcc14/dy_mcp_production}
OVERWRITE=${OVERWRITE:-0}

[[ -n "$JOB_ID" ]] || { echo "Usage: $0 JOB_ID"; exit 2; }

row=$(python3 - "$JOB_LIST" "$JOB_ID" <<'PY'
import csv
import sys

with open(sys.argv[1]) as source:
    for row in csv.DictReader(source):
        if row["job_id"] == sys.argv[2]:
            fields = [
                "seed",
                "n_events",
                "mcp_mass_GeV",
                "geometry",
                "emitter_name",
                "production_mode",
                "output_file",
                "job_id",
            ]
            print("|".join(row[field] for field in fields))
            break
    else:
        raise SystemExit(f"Job ID {sys.argv[2]} not found")
PY
)

IFS='|' read -r seed events mass geometry emitter mode output job_id <<< "$row"
mkdir -p "$(dirname "$output")" logs/jobs

if [[ -s "$output" && "$OVERWRITE" != 1 ]]; then
    echo "Output exists, skipping: $output"
    exit 0
fi

"$EXEC" \
    "$seed" 1 "$events" "$mass" "$geometry" "$emitter" "$mode" "$output" \
    --mode fixed-events \
    --n-events "$events" \
    --epsilon "${EPSILON:-0.01}" \
    --write-spectra "${WRITE_SPECTRA:-accepted}" \
    --spectra-prescale "${SPECTRA_PRESCALE:-1}" \
    --job-id "$job_id" \
    --batch \
    > "logs/jobs/job_${job_id}.log" 2>&1
