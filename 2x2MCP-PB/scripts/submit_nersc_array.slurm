#!/bin/bash
#SBATCH -C cpu
#SBATCH -q shared
#SBATCH -t 02:00:00
#SBATCH -N 1
#SBATCH -c 1
#SBATCH -J mcp_pb
#SBATCH -o logs/slurm/%x_%A_%a.out
#SBATCH -e logs/slurm/%x_%A_%a.err

set -eo pipefail

module load PrgEnv-gnu
module load python
conda activate "$SCRATCH/mcp-env"

set -u
export JOB_LIST=${JOB_LIST:-job_list_production.csv}
export EXEC=${EXEC:-./build-gcc14/pb_fww_production}

bash ./scripts/run_one_job.sh "${SLURM_ARRAY_TASK_ID}"

