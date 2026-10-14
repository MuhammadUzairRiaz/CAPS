#!/bin/bash
#SBATCH --account=p_thermo_frc
#SBATCH --partition=alpha
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --gres=gpu:1
#SBATCH --mem=110G
#SBATCH --job-name=caps-mace-v1
#SBATCH --time=24:00:00
#SBATCH --output=/data/horse/ws/muri912g-caps-ml/logs/mace-v1-%j.out
NAME=caps-mace-v1 EPOCHS=150 START_SWA=120 ./train.sh
