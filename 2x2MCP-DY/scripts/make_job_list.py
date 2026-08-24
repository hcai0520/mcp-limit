#!/usr/bin/env python3

import argparse
import csv
from pathlib import Path


def mass_tag(value):
    return f"{value:.6f}".replace(".", "p")


parser = argparse.ArgumentParser()
parser.add_argument("--masses", default="masses.txt")
parser.add_argument("--output", default="job_list.csv")
parser.add_argument("--output-dir", default="outputs/raw")
parser.add_argument("--n-events", type=int, default=10000)
parser.add_argument("--emitters", default="dy")
parser.add_argument("--geometry", default="2x2")
parser.add_argument("--seed-base", type=int, default=1000)
parser.add_argument("--shards", type=int, default=1)
parser.add_argument("--include-closed", action="store_true")
args = parser.parse_args()

if args.emitters.strip() != "dy":
    raise SystemExit("This independent generator supports --emitters dy only")

masses = []
for line in Path(args.masses).read_text().splitlines():
    line = line.strip()
    if line and not line.startswith("#"):
        masses.append(float(line))

rows = []
job_id = 0
for mass in masses:
    for shard in range(args.shards):
        output = (
            Path(args.output_dir)
            / "emitter_dy"
            / f"mass_{mass_tag(mass)}"
            / f"job_{job_id:06d}.root"
        )
        rows.append({
            "job_id": job_id,
            "emitter_name": "dy",
            "emitter_pdg": 0,
            "production_mode": "drell_yan",
            "seed": args.seed_base + job_id,
            "n_events": args.n_events,
            "mcp_mass_GeV": mass,
            "geometry": args.geometry,
            "output_file": output,
        })
        job_id += 1

fields = [
    "job_id",
    "emitter_name",
    "emitter_pdg",
    "production_mode",
    "mcp_mass_GeV",
    "seed",
    "n_events",
    "geometry",
    "output_file",
]
Path(args.output).parent.mkdir(parents=True, exist_ok=True)
with open(args.output, "w", newline="") as output:
    writer = csv.DictWriter(output, fieldnames=fields)
    writer.writeheader()
    writer.writerows(rows)

print(f"Wrote {len(rows)} jobs to {args.output}")
