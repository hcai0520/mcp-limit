# 2x2 MCP Drell-Yan production

This directory is independent of `2x2MCP-PythiaGen`. It generates

```text
q qbar -> gamma* -> chi chibar
```

for a Dirac millicharged fermion. The ROOT tree and geometry field names are compatible with Bruno's generator, but the production and normalization remain separate.

## NERSC environment

```bash
module load PrgEnv-gnu
module load python
conda activate "$SCRATCH/mcp-env"
```

## Build

```bash
cmake -S . -B build-gcc14 \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_CXX_COMPILER="$(which g++)" \
  -DCMAKE_PREFIX_PATH="$CONDA_PREFIX"

cmake --build build-gcc14 -j 4
```

## One run

```bash
srun -n 1 -c 1 ./build-gcc14/dy_mcp_production \
  12345 1 50000 0.10 2x2 dy drell_yan \
  outputs/dy_m0p1.root --epsilon 0.01 \
  --write-spectra all
```

Arguments are

```text
seed nThreads nEvents mcpMassGeV geometry emitterName productionMode outputFile
```

The default hard process includes PYTHIA ISR, beam remnants and hadronization. Use `--hard-only` only for a hard-process validation comparison.

The default PDF setting is `PDF:pSet = 5`, matching the validated DY setup. The MCP charge is set to zero in PYTHIA's particle table because the physical millicharge `epsilon` is already included in the custom hard matrix element; this prevents PYTHIA from applying unit-charge QED radiation.

The `2x2` geometry is kept identical to Bruno's implementation: the detector plane is 1.04 km downstream, with no detector changes made in this standalone generator.

## Production mass scan

The included `masses.txt` is identical to Bruno's mass grid. Generate one million events per shard and eight independent shards per mass with:

```bash
python scripts/make_job_list.py \
  --masses masses.txt \
  --n-events 1000000 \
  --shards 8 \
  --geometry 2x2 \
  --output job_list.csv \
  --output-dir outputs/raw
```

This produces 168 jobs: 21 masses times 8 shards. Submit the complete scan with at most 64 simultaneous one-core jobs:

```bash
mkdir -p logs/slurm
N=$(($(wc -l < job_list.csv) - 2))
sbatch --array=0-${N}%64 scripts/submit_nersc_array.slurm
```

Each task uses `nThreads = 1`, a distinct seed, and writes only detector-accepted MCP spectra. Existing nonempty ROOT files are skipped when a task is resubmitted.

Aggregate the finished shards with the same summary-tree workflow used by Bruno:

```bash
python scripts/aggregate_outputs.py outputs/raw \
  --csv outputs/aggregate_summary.csv \
  --root outputs/aggregate_summary.root
```

Do not add the raw DY MCP count to the meson MCP count. The DY contribution must first be normalized by its generated cross section relative to the common proton-interaction cross section. The physical merge is deliberately left for the later analysis stage.
