#!/bin/bash
# CAPS MACE v1 on Alpha (one A100): the first-round training set (whole environments held out in test.xyz).
# EPOCHS and NAME come from the job script, so the dry run and the training run the same command.
source /data/horse/ws/USER-caps-ml/env.sh
cd $WS/mace_v1
date +%T
mace_run_train --name $NAME --model MACE --train_file train.xyz --valid_fraction 0.05 --test_file test.xyz \
  --energy_key REF_energy --forces_key REF_forces --E0s average \
  --hidden_irreps "128x0e + 128x1o" --r_max 5.0 --batch_size 32 --max_num_epochs $EPOCHS \
  --swa --start_swa $START_SWA --ema --ema_decay 0.99 --amsgrad --lr 0.01 --scheduler_patience 5 --patience 30 \
  --default_dtype float64 --device cuda --seed 1 --num_workers 4 --error_table PerAtomRMSE --save_cpu --restart_latest \
  --model_dir models --checkpoints_dir checkpoints --results_dir results --log_dir logs
date +%T
