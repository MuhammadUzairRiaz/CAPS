#!/bin/bash
# The GPU half of the pilot on Alpha Centauri: the pilot's xyz files split over ARRAY tasks, one A100 each.
# Everything (inputs, results, logs, caches, scratch) is in the horse workspace; nothing is written to $HOME.
#   cd $WS/pilot && sbatch alpha_pilot.sbatch
#SBATCH --job-name=caps-ml-pilot
#SBATCH --account=p_yourproject
#SBATCH --partition=alpha
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --gres=gpu:1
#SBATCH --mem=60G
#SBATCH --time=03:00:00
#SBATCH --array=0-3
#SBATCH --output=/data/horse/ws/USER-caps-ml/logs/pilot-%A_%a.out

source /data/horse/ws/USER-caps-ml/env.sh
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
cd $WS/pilot
mapfile -t ALL < <(ls xyz/*.xyz | sort)
N=$SLURM_ARRAY_TASK_COUNT
MINE=()
for i in "${!ALL[@]}"; do
  if (( i % N == SLURM_ARRAY_TASK_ID )); then MINE+=("${ALL[$i]}"); fi
done
nvidia-smi --query-gpu=name,memory.total --format=csv,noheader
python3 run_gpu4pyscf.py gpu "${MINE[@]}"
