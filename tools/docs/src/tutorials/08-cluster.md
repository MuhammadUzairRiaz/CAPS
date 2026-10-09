---
title: Running CAPS on your own cluster
lede: Build the caps command line on a SLURM or PBS cluster, send runs there from the command line or from the Studio, follow them live, and continue a run that hit its time limit.
time: 20 minutes, plus queue time
level: intermediate
uses: caps job · Settings › Compute & remote · Jobs
---

## What you need

An account on a cluster you can reach with `ssh` using a key (an SSH agent holds it; CAPS never asks for or stores a password), a C++20 compiler with gfortran, CMake 3.24 or newer, and git. On most clusters these come from modules. A workstation you can ssh into works too: CAPS then runs jobs in the background without a queue.

CAPS runs threads, not MPI ranks. One run uses one node, so a job asks for one node, one task and N cores.

## 1. Build caps on the cluster

Log in and build the command line once, in your home folder. Load your site's compiler and CMake modules first.

```bash
module load GCC CMake            # your site's names
git clone https://github.com/MuhammadUzairRiaz/CAPS.git ~/CAPS/src/CAPS
cd ~/CAPS/src/CAPS
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DCAPS_BUILD_TESTS=OFF
cmake --build build -j 8 --target caps_cli
mkdir -p ~/CAPS/bin && ln -sf ~/CAPS/src/CAPS/build/cli/caps ~/CAPS/bin/caps
~/CAPS/bin/caps --version
```

Build inside a short interactive or batch job if your site does not allow compiling on login nodes. The Studio's Install / update caps (step 5) does exactly this as a batch job.

## 2. Describe your cluster once

`caps job profile` writes `~/CAPS/host.json`. Start from a generic preset and fill in your own values:

```bash
~/CAPS/bin/caps job profile --preset slurm
~/CAPS/bin/caps job profile --set "account=myproject,partition=compute,cpus=32,mem_per_cpu=2G,time=24:00:00"
~/CAPS/bin/caps job profile --set "modules=module purge; module load GCC"
```

The presets are `slurm`, `slurm-workspace` (sites with `ws_allocate`), `pbs` and `workstation`. `scratch` is `workspace`, `env` (a folder under `scratch_var`, such as `$SCRATCH` or `$TMPDIR`) or `none`. The profile stays on the cluster; nothing in it goes anywhere else.

## 3. Send a run

```bash
~/CAPS/bin/caps job new --title PS_melt --kind md --input cell.data --submit -- \
  caps md cell.data -o md.data --steps 1000000 --dump traj.lammpstrj --every 5000 --log thermo.csv
```

```text
/home/you/CAPS/PS_melt/md-1
submitted: 4815162
```

The job folder holds the inputs, `cmd.txt`, `job.sh` (open it to see exactly what runs) and `caps-job.json`. The script copies the inputs to scratch, runs the command there, copies everything back to `out/` and releases a workspace only after that copy succeeded.

Many runs at once go as one array job. Write one line per task, a title and a command:

```text
PS_seed1 caps md cell.data -o md.data --steps 1000000 --seed 1
PS_seed2 caps md cell.data -o md.data --steps 1000000 --seed 2
```

```bash
~/CAPS/bin/caps job new --kind md --array seeds.txt --input cell.data --submit
```

## 4. Follow it

```bash
~/CAPS/bin/caps job status ~/CAPS/PS_melt/md-1
~/CAPS/bin/caps job tail ~/CAPS/PS_melt/md-1 -f
~/CAPS/bin/caps job list
```

`status` shows the queue state, the node, the step reached with its temperature, pressure and density, and an ETA:

$$
t_\mathrm{left} = t_\mathrm{elapsed}\,\frac{1 - f}{f}
$$

where f is the fraction done. `tail -f` follows `run.log` in the scratch folder while the job runs.

## 5. The same from the Studio

In Settings › Compute & remote, add a host: its name, the login node, your user name and the port. Choose a preset under "Jobs on this host" and fill in your account, cores, memory, time limit, scratch kind and modules. Then:

1. Test connection reads caps's version, the node size, the queue and whether workspaces exist.
2. Install / update caps builds the command line on the host at the Studio's own version, as a short job.
3. On any simulation page (Dynamics, Equilibrate, Relax, React, Grow, Glass, Mechanics, the pull-out, Pack, Coarse-grain), choose the host in Run where. A panel shows the job's title, cores, memory, time limit, how often frames and thermo rows are written, and the exact caps command that will run.

Jobs then shows the queue, the run's own log (follow, search, errors in red), temperature and density as they come in, and the latest frame on request. Copy terminal commands gives you the `ssh` and `caps job tail -f` lines for the same job. Find my jobs adds jobs you started by hand on the cluster.

## 6. When a run hits its time limit

Ten minutes before the limit SLURM sends USR1 (the job script asks for it). caps writes a checkpoint, the results come back, and the job is marked `timeout`. Continue it as a new job:

```bash
~/CAPS/bin/caps job resume ~/CAPS/PS_melt/md-1 --submit
```

or press Resume from checkpoint in Jobs. MD continues from its last step, an equilibration from its stage and step, a Tg scan from the temperatures done, a tensile pull from the strain reached, and a recipe from the stages done. Logs and trajectories keep growing in the new job.

## 7. How many cores to ask for

```bash
~/CAPS/bin/caps job scaling --input cell.data --submit
~/CAPS/bin/caps job collect-bench ~/CAPS/_arrays/bench-1
```

```text
threads   ns/day   speedup  efficiency
      1   74.234      1.00        100 %
      2  128.421      1.73         86 %
      4  203.964      2.75         69 %
      8  212.021      2.86         36 %
suggested: 2 threads (the fastest at 70 % parallel efficiency or more)
```

That table is from a laptop with four performance cores, running a 648-atom polystyrene cell with the `workstation` preset (the tasks run one after another). Past four threads the run spills onto the slower cores and barely gains. A cluster node and a larger cell give a different curve. Run the check with a cell of the size you will simulate.

## Next

- [Theory: Cluster jobs](../theory/cluster-jobs.html) covers the job script, checkpoints and the scaling table.
- `caps job help` lists every subcommand.
