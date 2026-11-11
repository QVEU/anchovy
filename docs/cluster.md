# Cluster execution

[Back to the README](../README.md)

See [performance tuning](performance.md) for measured memory requirements.

## Allocate a node

**Get a node, then run the whole pipeline inside it.** One queue wait, then
full speed. Measured: 2 samples, 11,601 cells, 23,220 steps in **25 min 32 s**
on one 64-core node.

Interactively:

```bash
salloc --cpus-per-task=64 --mem=200G --time=24:00:00
snakemake -s workflow/Snakefile --configfile <config> --cores 64 --keep-going
```

`--cores` is required — snakemake refuses a local run without it (*cores have
to be specified for local execution*) — and it is what the thread counts
follow, so match it to `--cpus-per-task` or the node sits half idle.

Or as a batch job, taking the core count from SLURM so it cannot drift from
the allocation:

```bash
#!/bin/bash
#SBATCH --cpus-per-task=64
#SBATCH --mem=200G
#SBATCH --time=24:00:00
#SBATCH --job-name=anchovy

snakemake -s workflow/Snakefile --configfile <config> \
    --cores "$SLURM_CPUS_PER_TASK" --keep-going
```

**When you also need `--resources mem_mb=`.** Not above: `extract` takes every
core, so only one runs at a time whatever its memory reservation says, and the
cores alone keep the allocation safe. It matters once `extract_threads` is set
*below* `--cores`, because then several extract jobs start at once — four
16-thread extracts against a v3 whitelist is ~104 GB. Add
`--resources mem_mb=<your --mem in MB>` there and snakemake serialises what
will not fit; without a global budget the per-rule reservations are inert.

**Why not `--executor slurm`.** It submits each job or job group separately, so
every one takes its own turn in the shared queue — on a busy cluster, hours of
waiting against the 25 minutes above. An allocation costs one queue wait; the
executor costs many. It stays an option for when a single node genuinely is not
enough, and needs `snakemake-executor-plugin-slurm`, which `environment.yml`
declares (an older environment needs `conda install -n anchovy -c bioconda
-c conda-forge snakemake-executor-plugin-slurm`, or the run stops at argument
parsing with `invalid choice: 'slurm'`).

**Every rule declares `mem_mb` and `runtime`.** Inside an allocation those are
what `--resources mem_mb=` schedules against, so a job too big to fit beside
another waits rather than overrunning the node; under `--executor slurm` they
are what each job reserves, and without them a job goes in at the partition
default and `extract` is killed on a node with plenty free. The defaults are
computed in the Snakefile from measured figures — `extract`'s follows
`chemistry`, `extract_threads`, `extract_chunk_size` and the whitelist's own
line count automatically. Override any of them from the config:

```yaml
reads_estimate: 10000000   # reads expected after mapping; sizes the fasta stage
mem_mb:
  extract: 200000
runtime:                   # minutes
  extract: 480
```

`reads_estimate` matters because `fasta` loads the whole anchovy table — every
read's mapped sequence — at a measured 6.3 KB per read. At 10 million reads that
is ~84 GB, the second largest reservation in the pipeline. Nothing else scales
with it: `extract` is chunked, and the per-cell rules see one cell each.

**Threads follow `--cores`, and so does `extract`'s memory.** `map_threads`
and `extract_threads` default to the run's core count rather than a fixed 8 and
16, so `--cores 64` uses 64. Because the reservation is computed from
`extract_threads`, it follows too: against a v3 whitelist, 16 workers reserves
~26 GB and 64 reserves ~102 GB. Under `--executor slurm` there is no `--cores`
to follow, so set both in the config — `config_cluster.yaml` does.

**Only the per-cell stages scale out.** `map_reads`, `extract` and everything
after `merge_consensus` are single jobs on one node, so `extract` is capped at
that node's core count however many nodes you have.
