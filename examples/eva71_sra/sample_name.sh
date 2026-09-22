# sample_name.sh -- how a reads path becomes a sample name.
#
# SOURCED, NOT RUN. fetch.sh uses this to name the SAM it writes; run_cluster.sh
# uses it to tell the workflow which SAM to read. They have to agree, and when
# they did not the failure was silent and expensive: fetch.sh mapped the new
# FASTQ to its own name while the workflow read `sample` from the config, found
# the PREVIOUS run's SAM already sitting there, and reported everything up to
# date -- a full mapping run discarded, and results that look like the new
# sample but are not.
#
# One definition, sourced by both, so they cannot drift apart again.

sample_name_from_fastq() {
    local base
    base=$(basename "$1")
    base="${base%.gz}"
    base="${base%.fastq}"
    base="${base%.fq}"
    printf '%s' "$base"
}
