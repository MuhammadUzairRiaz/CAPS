#!/bin/bash
#SBATCH --job-name=caps-mace-v1-eval
#SBATCH --account=p_thermo_frc
#SBATCH --partition=alpha
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --gres=gpu:1
#SBATCH --mem=40G
#SBATCH --time=01:00:00
#SBATCH --output=/data/horse/ws/muri912g-caps-ml/logs/mace-v1-eval-%j.out
source /data/horse/ws/muri912g-caps-ml/env.sh
cd $WS/mace_v1
python3 -u eval_torsions.py test.xyz models/caps-mace-v1_stagetwo.model --aimnet
