#!/usr/bin/env python3

import argparse
import csv
import math
from pathlib import Path

try:
    import ROOT
except Exception as error:
    raise SystemExit(f"PyROOT is required: {error}")


parser = argparse.ArgumentParser()
parser.add_argument("inputs", nargs="*", default=["outputs/raw"])
parser.add_argument("--csv", default="outputs/aggregate_summary.csv")
parser.add_argument("--root", default="outputs/aggregate_summary.root")
args = parser.parse_args()

files = []
for item in args.inputs:
    path = Path(item)
    files += [path] if path.is_file() else list(path.rglob("*.root"))

count_fields = [
    "n_events_generated",
    "n_emitter_record_entries",
    "n_emitter_decayed_to_mcp",
    "n_mcp_pairs",
    "n_mcp_pairing_anomalies",
    "n_emitter_total",
    "n_mcp_total",
    "n_mcp_accepted",
    "n_mcp_wrong_mother",
    "n_mcp_no_mother",
]

aggregate = {}
for filename in files:
    source = ROOT.TFile.Open(str(filename))
    if not source or source.IsZombie():
        continue
    tree = source.Get("mcp_summary")
    if not tree:
        source.Close()
        continue

    for row in tree:
        key = (
            round(float(row.mcp_mass_GeV), 12),
            int(row.emitter_pdg),
            int(row.production_mode),
            int(row.geometry_id),
        )
        values = aggregate.setdefault(
            key,
            {field: 0 for field in count_fields}
            | {
                "sigma_gen_mb_sum": 0.0,
                "sigma_err_mb_squared_sum": 0.0,
                "sigma_count": 0,
                "epsilon_gen": float(row.epsilon_gen),
                "n_input_files": 0,
            },
        )
        for field in count_fields:
            if hasattr(row, field):
                values[field] += int(getattr(row, field))
        values["sigma_gen_mb_sum"] += float(row.sigma_gen_mb)
        values["sigma_err_mb_squared_sum"] += float(row.sigma_err_mb) ** 2
        values["sigma_count"] += 1
        values["n_input_files"] += 1
    source.Close()

Path(args.csv).parent.mkdir(parents=True, exist_ok=True)
fields = [
    "mcp_mass_GeV",
    "emitter_pdg",
    "production_mode",
    "geometry_id",
    "epsilon_gen",
] + count_fields + [
    "acceptance_fraction",
    "acceptance_uncertainty_binomial",
    "accepted_mcp_per_event",
    "n_input_files",
    "sigma_gen_mb_mean",
    "sigma_err_mb_mean",
    "sigma_over_epsilon2_mb",
    "accepted_mcp_sigma_over_epsilon2_mb",
]

rows = []
for key, values in sorted(aggregate.items()):
    total = values["n_mcp_total"]
    events = values["n_events_generated"]
    accepted = values["n_mcp_accepted"]
    acceptance = accepted / total if total else 0.0
    acceptance_error = (
        math.sqrt(max(0.0, acceptance * (1.0 - acceptance) / total))
        if total else 0.0
    )
    accepted_per_event = accepted / events if events else 0.0
    sigma_count = values["sigma_count"]
    sigma_mean = values["sigma_gen_mb_sum"] / sigma_count
    sigma_error = math.sqrt(values["sigma_err_mb_squared_sum"]) / sigma_count
    epsilon = values["epsilon_gen"]
    sigma_over_epsilon2 = sigma_mean / (epsilon * epsilon)

    row = {
        "mcp_mass_GeV": key[0],
        "emitter_pdg": key[1],
        "production_mode": key[2],
        "geometry_id": key[3],
        "epsilon_gen": epsilon,
    }
    for field in count_fields:
        row[field] = values[field]
    row.update({
        "acceptance_fraction": acceptance,
        "acceptance_uncertainty_binomial": acceptance_error,
        "accepted_mcp_per_event": accepted_per_event,
        "n_input_files": values["n_input_files"],
        "sigma_gen_mb_mean": sigma_mean,
        "sigma_err_mb_mean": sigma_error,
        "sigma_over_epsilon2_mb": sigma_over_epsilon2,
        "accepted_mcp_sigma_over_epsilon2_mb": (
            sigma_over_epsilon2 * accepted_per_event
        ),
    })
    rows.append(row)

with open(args.csv, "w", newline="") as output:
    writer = csv.DictWriter(output, fieldnames=fields)
    writer.writeheader()
    writer.writerows(rows)

root_output = ROOT.TFile(args.root, "RECREATE")
tree_output = ROOT.TTree("mcp_summary", "Aggregated DY MCP summary")
from array import array

buffers = {}
integer_fields = {
    "emitter_pdg", "production_mode", "geometry_id", "n_input_files"
}
count_field_set = set(count_fields)
for field in fields:
    if field in integer_fields:
        buffers[field] = array("i", [0])
        leaf = f"{field}/I"
    elif field in count_field_set:
        buffers[field] = array("l", [0])
        leaf = f"{field}/L"
    else:
        buffers[field] = array("d", [0.0])
        leaf = f"{field}/D"
    tree_output.Branch(field, buffers[field], leaf)

for row in rows:
    for field, value in row.items():
        buffers[field][0] = value
    tree_output.Fill()

tree_output.Write()
root_output.Close()
print(f"Wrote {len(rows)} mass points to {args.csv} and {args.root}")
