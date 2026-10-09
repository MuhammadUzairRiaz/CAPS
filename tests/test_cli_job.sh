#!/bin/bash
# caps job on a fake SLURM cluster (tests/fake_cluster/bin: sbatch, squeue, sacct, scancel, ws_allocate, ws_release
# run the job scripts locally): a job made and submitted, followed while it runs, copied back from its workspace and the
# workspace released; an array of two tasks; a failing command; a cancelled job. Short runs only.
# usage: test_cli_job.sh CAPS SOURCE_DIR
set -u
CAPS_BIN=$1; SRC=$2
T=$(mktemp -d); trap 'rm -rf "$T"' EXIT
export HOME=$T/home FAKE_STATE=$T/state FAKE_WS=$T/ws CAPS_HOME=$SRC PATH=$SRC/tests/fake_cluster/bin:$PATH
mkdir -p "$HOME/CAPS/bin" "$FAKE_STATE" "$FAKE_WS"; ln -s "$CAPS_BIN" "$HOME/CAPS/bin/caps"
C=$HOME/CAPS/bin/caps
fail() { echo "FAIL: $*"; exit 1; }
wait_state() {   # wait_state DIR STATE
  for _ in $(seq 1 120); do grep -q "\"state\": *\"$2\"" "$1/caps-job.json" 2>/dev/null && return 0; sleep 0.5; done
  cat "$1/caps-job.json"; fail "$1 never reached $2"
}
"$C" job profile --preset slurm-workspace > /dev/null || fail profile
"$C" job profile --set "cpus=2,time=00:10:00,ws_filesystem=" > /dev/null || fail "profile --set"
cp "$SRC/samples/ps_melt.data" "$T/cell.data" 2> /dev/null || "$C" grow -o "$T/cell.data" --chains 2 --dp 4 --density 0.5 --seed 1 > /dev/null
# 1 a job: made, submitted, finished, copied back, the workspace released
"$C" job new --title "PS cell" --kind md --input "$T/cell.data" --submit -- caps md cell.data -o md.data --steps 200 --thermo 50 \
  --dump traj.lammpstrj --every 100 --log thermo.csv --quiet > "$T/new.txt" || { cat "$T/new.txt"; fail "job new"; }
J=$HOME/CAPS/PS_cell/md-1
[ -f "$J/job.sh" ] || fail "no job.sh"
wait_state "$J" finished
for f in md.data thermo.csv progress.jsonl traj.lammpstrj run.log; do [ -s "$J/out/$f" ] || fail "out/$f missing"; done
[ -z "$(ls "$FAKE_WS")" ] || fail "the workspace was not released"
"$C" job status "$J" | grep -q "reached  md · step 200" || { "$C" job status "$J"; fail status; }
"$C" job poll "$J" | grep -q '"state":"finished"' || fail poll
"$C" job list | grep -q "PS_cell/md-1" || fail list
# 2 an array of two tasks
printf 'A caps md cell.data -o md.data --steps 100 --seed 1 --quiet\nB caps md cell.data -o md.data --steps 100 --seed 2 --quiet\n' > "$T/list.txt"
"$C" job new --kind md --array "$T/list.txt" --input "$T/cell.data" --submit > /dev/null || fail "array"
wait_state "$HOME/CAPS/A/md-1" finished
wait_state "$HOME/CAPS/B/md-1" finished
[ -s "$HOME/CAPS/B/md-1/out/md.data" ] || fail "array task B: no result"
# 3 a failing command: failed, its exit code kept
"$C" job new --title bad --kind md --submit -- caps md missing.data -o x.data > /dev/null
wait_state "$HOME/CAPS/bad/md-1" failed
grep -q '"exit": *[1-9]' "$HOME/CAPS/bad/md-1/caps-job.json" || fail "the exit code"
# 4 a cancelled job
"$C" job new --title slow --kind test --submit -- bash -c 'sleep 30' > /dev/null
sleep 1
"$C" job cancel "$HOME/CAPS/slow/test-1" > /dev/null || fail cancel
grep -q '"state": *"cancelled"' "$HOME/CAPS/slow/test-1/caps-job.json" || fail "cancelled state"
echo "caps job on the fake cluster: ok"
