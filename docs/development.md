# Development and reproducibility

[Back to the README](../README.md)

## For developers

anchovy comes with an automated test suite that checks each step still produces
the expected results. To run it:

```bash
pytest                              # run all tests
pytest tests/test_annotate.py -v   # run the tests for one step
```

Test data lives in `tests/data/`. The `tests/make_*_fixtures.py` scripts
regenerate it if ever needed.

`anchovy_v3/` is the archived original this version was migrated from. It is
kept because `tests/data/golden/` is frozen from its output, so the port can be
checked against it. Nothing in the pipeline runs it.

## Dependency environments

`environment.yml` declares supported dependency ranges for development and
platforms beyond Linux x86_64. The checked-in
[`locks/linux-64/explicit.txt`](../locks/linux-64/explicit.txt) records exact Conda
package URLs, builds and SHA-256 hashes from a tested Linux x86_64 environment.
It includes Python, mapping tools, Snakemake, R reporting, and test tools. The
Anchovy source itself is installed separately from this checkout.

For a new Linux x86_64 environment, run from the repository root:

```bash
conda create -n anchovy-locked --file locks/linux-64/explicit.txt
conda activate anchovy-locked
python -m pip install -e . --no-deps --no-build-isolation
pytest -v
```

`--no-build-isolation` uses the build tools pinned in the lock instead of fetching
new setuptools/wheel versions. Record the Anchovy Git commit alongside the lock
when archiving an analysis. The lock does not pin the OS, kernel, input data or
configuration, and does not guarantee byte-identical figures across machines.
Other architectures should solve `environment.yml` in their native environment;
do not use the Linux x86_64 lock on macOS, ARM or Windows.

To refresh the lock after changing dependency requirements, create a **new**
environment from `environment.yml`, install this checkout, and validate it before
exporting. Keep the previous environment until the replacement passes:

```bash
conda env create -n anchovy-lock-refresh -f environment.yml
conda activate anchovy-lock-refresh
python -m pip install -e . --no-deps --no-build-isolation
pytest -v
conda list --explicit --sha256 > locks/linux-64/explicit.txt
```

Then create a separate environment from the exported lock and rerun the suite.
The suite includes the report-enabled bundled workflow; check that it executed
rather than skipped. Commit the dependency-spec changes and refreshed lock
together. Use a Conda version supporting SHA-256 explicit exports; never disable
artifact or TLS verification to refresh a lock.

## CI and integration checks

CI installs `environment.yml` so new compatible dependency releases are tested.
It requires minimap2, samtools, Snakemake, R and Pandoc, and loads the reporting
packages before pytest. The suite runs the bundled workflow both with reports
disabled and enabled. The report smoke test checks visible sections, genotype
labels, and embedded plots. Structural R Markdown tests separately guard colour
and network conventions.

Without the workflow tools, local integration tests skip; without R, the report
smoke test skips. Pure Python tests remain usable. A skip is not evidence that
reporting works.

## Restricted home directories and managed networks

If Snakemake cannot create its user cache, point it at a writable directory:

```bash
export XDG_CACHE_HOME="$PWD/.snakemake/user-cache"
mkdir -p "$XDG_CACHE_HOME"
```

For Conda behind an HTTPS proxy, use the administrator-provided CA bundle through
`ssl_verify` and retain the inherited proxy settings. Do not disable TLS
verification. Home-directory environment-registration permissions may still be
needed when creating Conda environments.
