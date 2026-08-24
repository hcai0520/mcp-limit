# 2x2 MCP proton bremsstrahlung production

This directory is independent of `2x2MCP-PythiaGen` and `2x2MCP-DY`. It generates

```text
p -> p gamma*, gamma* -> chi chibar
```

with the conventional Fermi-Weizsaecker-Williams splitting kernel and the full time-like proton VMD form factor from Appendix B of arXiv:2406.01668.

The baseline is a 120 GeV fixed-target proton beam. The MCP is a Dirac fermion. Production weights are calculated per epsilon squared, and the reference `epsilon` is used only to provide the corresponding physical `sigma_gen_mb` field.

The implementation applies the FWW validity requirements event by event:

```text
pT < 0.1 E(gamma*)
pT < 1 GeV
|qmin^2| < LambdaQCD^2, LambdaQCD = 0.25 GeV
E(gamma*) and Ep-E(gamma*) > 5 mp
E(gamma*) and Ep-E(gamma*) > 5 m(gamma*)
```

It uses the first-version approximation `sigma(s') / sigma(s) = 1`. The reference inelastic cross section defaults to `38.4538 mb`, matching the current Bruno normalization used in this reproduction.

## NERSC environment and build

```bash
module load PrgEnv-gnu
module load python
conda activate "$SCRATCH/mcp-env"

cmake -S . -B build-gcc14 \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_CXX_COMPILER="$(which g++)" \
  -DCMAKE_PREFIX_PATH="$CONDA_PREFIX"

cmake --build build-gcc14 -j 4
```

## One run

Run this inside a CPU allocation:

```bash
srun -n 1 -c 1 ./build-gcc14/pb_fww_production \
  12345 1 1000000 0.10 2x2 pb proton_bremsstrahlung \
  outputs/pb_m0p1.root \
  --write-spectra accepted
```

## Full 50-point production

Create 10 independent one-million-trial shards per mass:

```bash
python scripts/make_job_list.py \
  --masses masses.txt \
  --n-trials 1000000 \
  --shards 10 \
  --geometry 2x2 \
  --output job_list_production.csv \
  --output-dir "$SCRATCH/mcp-pb-production/raw"
```

This creates 500 jobs. Submit at most 64 simultaneous one-core tasks:

```bash
mkdir -p logs/slurm logs/jobs
N=$(($(wc -l < job_list_production.csv) - 2))

JOB_LIST="$PWD/job_list_production.csv" \
EXEC="$PWD/build-gcc14/pb_fww_production" \
sbatch --array=0-${N}%64 scripts/submit_nersc_array.slurm
```

Existing nonempty ROOT outputs are skipped when the array is resubmitted.

## Aggregate the completed shards

```bash
mkdir -p outputs
python scripts/aggregate_outputs.py \
  "$SCRATCH/mcp-pb-production/raw" \
  --csv outputs/aggregate_summary.csv \
  --root outputs/aggregate_summary.root
```

The main columns are:

```text
sigma_over_epsilon2_mb
accepted_mcp_sigma_over_epsilon2_mb
acceptance_fraction
probability_pair_over_epsilon2
accepted_mcp_probability_over_epsilon2
```

`acceptance_fraction` is the weighted single-MCP geometric acceptance. It is not the raw fraction of retained Monte Carlo points.

## Channel separation

The full VMD form factor contains rho and omega resonance contributions. Keep this PB result separate from explicit vector-meson decay results during validation. Do not directly add full-VMD PB to rho/omega meson production until a double-counting prescription has been selected.

