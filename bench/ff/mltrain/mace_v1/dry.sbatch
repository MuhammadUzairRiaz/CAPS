#!/bin/bash
#SBATCH --account=p_thermo_frc
#SBATCH --partition=alpha
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --gres=gpu:1
#SBATCH --mem=110G
#SBATCH --job-name=caps-mace-v1-dry
#SBATCH --time=01:00:00
#SBATCH --output=/data/horse/ws/muri912g-caps-ml/logs/mace-v1-dry-%j.out
NAME=caps-mace-v1-dry EPOCHS=1 START_SWA=1 ./train.sh
