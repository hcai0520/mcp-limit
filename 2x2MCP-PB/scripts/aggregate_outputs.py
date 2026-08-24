#!/usr/bin/env python3

import argparse
import csv
import math
from array import array
from pathlib import Path

try:
    import ROOT
except Exception as error:
    raise SystemExit(f"PyROOT is required: {error}")


def standard_error(total, total2, count):
    if count < 2:
        return 0.0
    numerator = max(0.0, total2 - total * total / count)
    return math.sqrt(numerator / (count * (count - 1)))


def acceptance_error(values):
    count = values["n_events_generated"]
    sum_x = values["sum_pair_weights"]
    sum_y = values["sum_accepted_weights"]
    if count < 2 or sum_x <= 0.0:
        return 0.0

    mean_x = sum_x / count
    mean_y = sum_y / count
    var_mean_x = max(
        0.0, values["sum_pair_weights2"] - count * mean_x * mean_x
    ) / (count * (count - 1))
    var_mean_y = max(
        0.0, values["sum_accepted_weights2"] - count * mean_y * mean_y
    ) / (count * (count - 1))
    cov_mean = (
        values["sum_pair_accepted_weights"] - count * mean_x * mean_y
    ) / (count * (count - 1))
    grad_x = -mean_y / (2.0 * mean_x * mean_x)
    grad_y = 1.0 / (2.0 * mean_x)
    variance = (
        grad_x * grad_x * var_mean_x
        + grad_y * grad_y * var_mean_y
        + 2.0 * grad_x * grad_y * cov_mean
    )
    return math.sqrt(max(0.0, variance))


parser = argparse.ArgumentParser()
parser.add_argument("inputs", nargs="*", default=["outputs/raw"])
parser.add_argument("--csv", default="outputs/aggregate_summary.csv")
parser.add_argument("--root", default="outputs/aggregate_summary.root")
args = parser.parse_args()

files = []
for item in args.inputs:
    path = Path(item)
    files += [path] if path.is_file() else list(path.rglob("*.root"))

sum_fields = [
    "n_events_generated",
    "n_phase_space_valid",
    "n_mcp_pairs",
    "n_mcp_total",
    "n_mcp_accepted",
    "n_spectra_written",
    "cut_pt_relative",
    "cut_pt_absolute",
    "cut_qmin",
    "cut_energy",
    "cut_sprime",
    "cut_kernel",
    "sum_pair_weights",
    "sum_pair_weights2",
    "sum_accepted_weights",
    "sum_accepted_weights2",
    "sum_pair_accepted_weights",
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

    for entry in tree:
        key = (
            round(float(entry.mcp_mass_GeV), 12),
            int(entry.emitter_pdg),
            int(entry.production_mode),
            int(entry.geometry_id),
        )
        values = aggregate.setdefault(
            key,
            {field: 0.0 for field in sum_fields}
            | {
                "epsilon_gen": float(entry.epsilon_gen),
                "beam_energy_GeV": float(entry.beam_energy_GeV),
                "sigma_inelastic_mb": float(entry.sigma_inelastic_mb),
                "interaction_cross_section_ratio": float(
                    entry.interaction_cross_section_ratio
                ),
                "n_input_files": 0,
            },
        )
        for field in sum_fields:
            values[field] += float(getattr(entry, field))
        values["n_input_files"] += 1
    source.Close()

fields = [
    "mcp_mass_GeV",
    "emitter_pdg",
    "production_mode",
    "geometry_id",
    "epsilon_gen",
    "beam_energy_GeV",
    "sigma_inelastic_mb",
    "interaction_cross_section_ratio",
    "n_input_files",
] + sum_fields + [
    "probability_pair_over_epsilon2",
    "probability_pair_error_over_epsilon2",
    "accepted_mcp_probability_over_epsilon2",
    "accepted_mcp_probability_error_over_epsilon2",
    "acceptance_fraction",
    "acceptance_uncertainty_mc",
    "effective_sample_size",
    "sigma_gen_mb",
    "sigma_err_mb",
    "sigma_over_epsilon2_mb",
    "sigma_error_over_epsilon2_mb",
    "accepted_mcp_sigma_over_epsilon2_mb",
    "accepted_mcp_sigma_error_over_epsilon2_mb",
]

rows = []
for key, values in sorted(aggregate.items()):
    trials = int(values["n_events_generated"])
    sum_pair = values["sum_pair_weights"]
    sum_accepted = values["sum_accepted_weights"]
    probability = sum_pair / trials if trials else 0.0
    probability_error = standard_error(
        sum_pair, values["sum_pair_weights2"], trials
    )
    accepted_probability = sum_accepted / trials if trials else 0.0
    accepted_probability_error = standard_error(
        sum_accepted, values["sum_accepted_weights2"], trials
    )
    acceptance = sum_accepted / (2.0 * sum_pair) if sum_pair else 0.0
    sigma_reference = values["sigma_inelastic_mb"]
    epsilon = values["epsilon_gen"]
    sigma_over_epsilon2 = sigma_reference * probability
    sigma_error_over_epsilon2 = sigma_reference * probability_error
    accepted_sigma_over_epsilon2 = sigma_reference * accepted_probability
    accepted_sigma_error_over_epsilon2 = (
        sigma_reference * accepted_probability_error
    )

    row = {
        "mcp_mass_GeV": key[0],
        "emitter_pdg": key[1],
        "production_mode": key[2],
        "geometry_id": key[3],
        "epsilon_gen": epsilon,
        "beam_energy_GeV": values["beam_energy_GeV"],
        "sigma_inelastic_mb": sigma_reference,
        "interaction_cross_section_ratio": values[
            "interaction_cross_section_ratio"
        ],
        "n_input_files": values["n_input_files"],
    }
    for field in sum_fields:
        row[field] = int(values[field]) if field.startswith(("n_", "cut_")) else values[field]
    row.update({
        "probability_pair_over_epsilon2": probability,
        "probability_pair_error_over_epsilon2": probability_error,
        "accepted_mcp_probability_over_epsilon2": accepted_probability,
        "accepted_mcp_probability_error_over_epsilon2": accepted_probability_error,
        "acceptance_fraction": acceptance,
        "acceptance_uncertainty_mc": acceptance_error(values),
        "effective_sample_size": (
            sum_pair * sum_pair / values["sum_pair_weights2"]
            if values["sum_pair_weights2"] else 0.0
        ),
        "sigma_gen_mb": sigma_over_epsilon2 * epsilon * epsilon,
        "sigma_err_mb": sigma_error_over_epsilon2 * epsilon * epsilon,
        "sigma_over_epsilon2_mb": sigma_over_epsilon2,
        "sigma_error_over_epsilon2_mb": sigma_error_over_epsilon2,
        "accepted_mcp_sigma_over_epsilon2_mb": accepted_sigma_over_epsilon2,
        "accepted_mcp_sigma_error_over_epsilon2_mb": (
            accepted_sigma_error_over_epsilon2
        ),
    })
    rows.append(row)

Path(args.csv).parent.mkdir(parents=True, exist_ok=True)
with open(args.csv, "w", newline="") as output:
    writer = csv.DictWriter(output, fieldnames=fields)
    writer.writeheader()
    writer.writerows(rows)

Path(args.root).parent.mkdir(parents=True, exist_ok=True)
root_output = ROOT.TFile(args.root, "RECREATE")
tree_output = ROOT.TTree("mcp_summary", "Aggregated FWW PB summary")
integer_fields = {
    "emitter_pdg",
    "production_mode",
    "geometry_id",
    "n_input_files",
} | {field for field in sum_fields if field.startswith(("n_", "cut_"))}
buffers = {}
for field in fields:
    if field in integer_fields:
        buffers[field] = array("q", [0])
        leaf = f"{field}/L"
    else:
        buffers[field] = array("d", [0.0])
        leaf = f"{field}/D"
    tree_output.Branch(field, buffers[field], leaf)

for row in rows:
    for field in fields:
        buffers[field][0] = row[field]
    tree_output.Fill()

tree_output.Write()
root_output.Close()
print(f"Wrote {len(rows)} mass points to {args.csv} and {args.root}")

