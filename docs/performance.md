# Performance and troubleshooting

[Back to the README](../README.md)

For allocation and submission commands, see [cluster execution](cluster.md).

## If `extract` is killed

A run that stops with `died with <Signals.SIGKILL: 9>` and nothing else was
killed by the kernel for using too much memory. SIGKILL cannot be caught, so the
stage gets no chance to say so itself.

`extract` processes the SAM a chunk at a time, so **its footprint does not
depend on how big the run is** — only on the chunk and the worker count:

```
GB  =  W  +  extract_threads × ( W  +  4.6 KB × extract_chunk_size )

       W = 0.4 GB   v2 whitelist    (737,280 barcodes)
       W = 2.8 GB   v3 whitelist  (6,794,880 barcodes)
       W ≈ 0        a small run-specific list
```

**`W` is the whitelist's lookup tables, and every worker holds its own copy.**
They are read-only and identical, but fork's copy-on-write does not save them —
CPython's refcounts touch every object. Measured: 0.42 GB per worker at v2,
2.75 GB at v3, linear in the worker count.

On a full 10x whitelist `W` dwarfs the chunk term (2.75 GB against 0.46 GB at a
100,000-read chunk), so **`extract_threads` is effectively the only memory knob
that matters** and chunking buys much less than it does against a small list:

`W` also depends on `max_barcode_errors`: the per-barcode template array is
only built when something will scan it, which is only when that setting is
unset. A bounded search resolves to one index and rebuilds one template, so
setting it cuts memory as well as runtime.

| | `W` bounded | `W` unset |
|---|---:|---:|
| v2 | 0.2 GB | 0.4 GB |
| v3 | **1.1 GB** | 2.8 GB |

| v3, 100,000-read chunks | `max_barcode_errors: 2` | unset |
|---|---:|---:|
| 4 threads | 7 GB | 16 GB |
| 8 threads | 14 GB | 28 GB |
| 16 threads | 26 GB | 54 GB |
| 64 threads | 102 GB | 208 GB |

The v3 tables also build in 10 s rather than 31 s, once, before the first read.


It does not depend on how big the run is — that is what chunking bought. It does
depend heavily on the whitelist: at the defaults (16 workers, 100,000-read
chunks) a small custom list costs about 7 GB, a v2 whitelist about 14 GB and a
v3 whitelist about 54 GB.

Lower `extract_threads` first on a tight node — against a full whitelist it is
the term that counts, and it costs wall time in proportion.
`extract_chunk_size` is worth trimming only when the whitelist is small.

The stage prints what it is about to need before it starts, so a subsequent kill
is at least diagnosable.

## If most barcodes come back 2 errors off

Check the line `extract` prints after the match-quality histogram. If it warns
that matched blocks do not begin with the signature's 5' handle, the search
window is misplaced and the barcodes are fine — the block is shifted, so the
barcode was sliced from the wrong 16 bases and every affected read scores the
same fixed penalty. It looks exactly like two sequencing errors, which is the
point of the warning.

The usual cause is the aligner starting a base or two inside the construct's
3' handle: `TTTCTTATAT` is T-rich and the cDNA after it often begins with
polyT, so minimap2 extends back into it. The window reaches `downstream_window`
bases past the alignment offset to absorb that; raise it if the warning
persists. If the warning covers most of the run, check `chemistry` and
`signature` instead — a v2 signature against v3 data shifts every read.

## If `extract` is slow

Barcode assignment is about **90% of the stage's time** on real data — the
signature search is noise beside it — and nearly all of that is spent on reads
whose barcode is *not* an exact whitelist hit. So `max_barcode_errors` is the
lever, not the hardware. Measured end to end on 40,000 reads with 44% of them
inexact, against a 100,000-barcode whitelist, 4 workers:

| `max_barcode_errors` | whole stage | reads assigned |
|---|---:|---:|
| `0` | 1.8 s | 56% |
| `1` | 2.1 s | 85% |
| **`2`** | **5.4 s** | **100%** |

**`2` now costs 2.5× what `1` does, not the 5.2× it used to.** The search stops
at the first shell that hits, so a read whose barcode is one error out is
settled in 48 candidates instead of 1,128 — the 2-error neighbourhood is only
enumerated for reads that actually need it. That halved the `2` row (11.6 s →
5.4 s) and left `1` alone.

That changes the recommendation. `1` used to be the pragmatic default because
`2` was punitive; now `2` keeps every read for a cost much closer to `1`'s.
Prefer `2` unless the stage is genuinely your bottleneck. The 15% that only `2`
recovers are the reads whose barcode needed *two* corrections — the least
trustworthy assignments in the run — so dropping them is a judgement about your
experiment, and it is no longer mostly a performance question.

Note also that the bounded search is **not** unconditionally faster than the
scan it replaces. Its cost is fixed whatever the whitelist holds, while the scan
is linear in the list — so against 10x's 737,280 entries the bound wins
enormously, but against a run-specific Cell Ranger list of a couple of thousand
barcodes, leaving it unset was *faster* than setting it to `2` (8.4 s against
10.8 s for that step), and loses no reads. Check `wc -l` on your whitelist first.

**More cores help; more nodes do not.** `extract` is a single Snakemake job
using one machine's process pool — `--executor slurm` distributes the per-cell
jobs, not this one — so `extract_threads` is the hardware lever, and it is
capped at one node.

Memory is bounded by the chunk rather than the run, but it is *not* free of the
thread count: each worker holds its own copy of the whitelist's lookup tables
(2.75 GB against a v3 list), so raising `extract_threads` raises the
reservation roughly in proportion. Lowering `extract_chunk_size` does not pay
for it on a full whitelist — that term is 0.46 GB per worker against the
tables' 2.75. Budget for `threads × W` and check what the stage prints on the
way in.
