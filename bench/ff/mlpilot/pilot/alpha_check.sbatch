#!/bin/bash
# The GPU side of the outlier check: the three outliers with a denser grid and without density fitting
#SBATCH --job-name=caps-ml-gpuchk
#SBATCH --account=p_thermo_frc
#SBATCH --partition=alpha
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --gres=gpu:1
#SBATCH --mem=60G
#SBATCH --time=01:00:00
#SBATCH --output=/data/horse/ws/muri912g-caps-ml/logs/gpuchk-%j.out
source /data/horse/ws/muri912g-caps-ml/env.sh
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
cd $WS/pilot
F="xyz/12f9bb3a1aa1_0.xyz xyz/2d9204ae7408_0.xyz xyz/bafc1feb3458_0.xyz"
python3 run_gpu4pyscf.py gpu_grid150 $F --grid 150,974
python3 run_gpu4pyscf.py gpu_nodf $F --no-df
