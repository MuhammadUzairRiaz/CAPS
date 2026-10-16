#!/bin/bash
#SBATCH --job-name=caps-mace-v2-eval
#SBATCH --account=p_thermo_frc
#SBATCH --partition=alpha
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --gres=gpu:1
#SBATCH --mem=40G
#SBATCH --time=01:00:00
#SBATCH --output=/data/horse/ws/muri912g-caps-ml/logs/mace-v2-eval-%j.out
# v1 and v2 torsion profiles against DFT on the same held-out fragments: v1's test set (the conformers v1 was judged on)
# and v2's (the same 131 fragments with round 2's conformers added), each model on each
source /data/horse/ws/muri912g-caps-ml/env.sh
cd $WS/mace_v2
V1=$WS/mace_v1/models/caps-mace-v1_stagetwo.model
V2=$WS/mace_v2/models/caps-mace-v2_stagetwo.model
for m in v1 v2; do
  for t in v1 v2; do
    model=V${m#v}; mkdir -p eval/${m}_on_${t}test
    python3 -u $WS/mace_v1/eval_torsions.py $WS/mace_$t/test.xyz ${!model} --out eval/${m}_on_${t}test
  done
done
