#!/bin/bash
# The training set of CAPS's own ML potential on Capella: one task per tasks/NN.txt, one H100 each.
# Conformers (AIMNet2 relaxed scans + distorted copies), then GPU4PySCF energies and forces for all of them.
# Everything (inputs, results, logs, caches, scratch) is in the horse workspace; nothing is written to $HOME.
# Resubmitting continues where a task stopped.
#   cd $WS/trainset && sbatch alpha_trainset.sbatch
#SBATCH --job-name=caps-ml-trainset
#SBATCH --account=p_yourproject
#SBATCH --partition=capella
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:1
#SBATCH --mem=110G
#SBATCH --time=24:00:00
#SBATCH --array=0-39
#SBATCH --output=/data/horse/ws/USER-caps-ml/logs/trainset-capella-%A_%a.out

source /data/horse/ws/USER-caps-ml/env-capella.sh
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
cd $WS/trainset
mapfile -t ENVS < tasks/$(printf %02d $SLURM_ARRAY_TASK_ID).txt
nvidia-smi --query-gpu=name,memory.total --format=csv,noheader
echo "== conformers $(date +%T)"
python3 gen_conformers.py . "${ENVS[@]}"
echo "== DFT $(date +%T)"
python3 dft_confs.py . "${ENVS[@]}"
echo "== end $(date +%T)"
