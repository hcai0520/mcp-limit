#!/usr/bin/env python3

import argparse
import csv
from pathlib import Path


def mass_tag(value):
    return f"{value:.6f}".replace(".", "p")


parser = argparse.ArgumentParser()
parser.add_argument("--masses", default="masses.txt")
parser.add_argument("--output", default="job_list_production.csv")
parser.add_argument("--output-dir", default="outputs/raw")
parser.add_argument("--n-trials", type=int, default=1_000_000)
parser.add_argument("--geometry", default="2x2")
parser.add_argument("--seed-base", type=int, default=1000)
parser.add_argument("--shards", type=int, default=10)
args = parser.parse_args()

masses = [
    float(line.strip())
    for line in Path(args.masses).read_text().splitlines()
    if line.strip() and not line.lstrip().startswith("#")
]

rows = []
for mass in masses:
    for shard in range(args.shards):
        job_id = len(rows)
        output = (
            Path(args.output_dir)
            / "emitter_pb"
            / f"mass_{mass_tag(mass)}"
            / f"job_{job_id:06d}.root"
        )
        rows.append({
            "job_id": job_id,
            "emitter_name": "pb",
            "emitter_pdg": 22,
            "production_mode": "proton_bremsstrahlung",
            "mcp_mass_GeV": mass,
            "seed": args.seed_base + job_id,
            "n_trials": args.n_trials,
            "geometry": args.geometry,
            "output_file": output,
        })

fields = [
    "job_id",
    "emitter_name",
    "emitter_pdg",
    "production_mode",
    "mcp_mass_GeV",
    "seed",
    "n_trials",
    "geometry",
    "output_file",
]
Path(args.output).parent.mkdir(parents=True, exist_ok=True)
with open(args.output, "w", newline="") as output:
    writer = csv.DictWriter(output, fieldnames=fields)
    writer.writeheader()
    writer.writerows(rows)

print(f"Wrote {len(rows)} jobs to {args.output}")

