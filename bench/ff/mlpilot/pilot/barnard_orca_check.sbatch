#!/bin/bash
# The ORCA half of the pilot on Barnard; inputs, scratch, results and logs all in the horse workspace.
#   cd /data/horse/ws/USER-caps-ml/pilot && sbatch barnard_orca.sbatch
#SBATCH --job-name=caps-ml-orcachk
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G
#SBATCH --time=06:00:00
#SBATCH --array=1-6
#SBATCH --output=/data/horse/ws/USER-caps-ml/logs/orcachk-%A_%a.out

module purge >/dev/null 2>&1
module load release/24.04 GCC/12.3.0 OpenMPI/4.1.5 ORCA/6.0.1-avx2 >/dev/null 2>&1 || module load ORCA/6.0.1-avx2
d=/data/horse/ws/USER-caps-ml/pilot/orca_check/$(sed -n "${SLURM_ARRAY_TASK_ID}p" /data/horse/ws/USER-caps-ml/pilot/orca_check/list.txt)
grep -q "ORCA TERMINATED NORMALLY" "$d/point.out" 2>/dev/null && exit 0
cd "$d" && /usr/bin/time -f "wall %e s" $(which orca) point.inp > point.out 2>&1
