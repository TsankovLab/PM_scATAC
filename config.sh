# config.sh: root folders read and written by the shell / job scripts; the same roots as config.R (see its header for what each one holds).
# Sourced by every shell script: found by searching upward from the working directory (run or submit scripts from inside the repository),
# or set MESO_CONFIG=/path/to/config.sh. Set the roots by editing the placeholders, by environment variables (MESO_DATA_DIR, ...), or in
# config.local.sh next to this file (git-ignored, sourced last). The roots are exported, so Python steps see them too.
REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DATA_DIR="${MESO_DATA_DIR:-/path/to/data}"
RESULTS_DIR="${MESO_RESULTS_DIR:-/path/to/results}"
REF_DIR="${MESO_REF_DIR:-/path/to/references}"
GENOME_DIR="${MESO_GENOME_DIR:-/path/to/genomes}"
SCRATCH_DIR="${MESO_SCRATCH_DIR:-/path/to/scratch}"
CONDA_ENVS="${MESO_CONDA_ENVS:-/path/to/conda/envs}"
CONDA_INIT="${MESO_CONDA_INIT:-/path/to/conda/etc/profile.d/conda.sh}"
LSF_ACCOUNT="${MESO_LSF_ACCOUNT:-<LSF_ACCOUNT>}"
if [ -f "$REPO_DIR/config.local.sh" ]; then source "$REPO_DIR/config.local.sh"; fi
export REPO_DIR DATA_DIR RESULTS_DIR REF_DIR GENOME_DIR SCRATCH_DIR CONDA_ENVS CONDA_INIT LSF_ACCOUNT
for v in DATA_DIR RESULTS_DIR REF_DIR GENOME_DIR SCRATCH_DIR; do case "${!v}" in /path/to/*) echo "config.sh: $v not set yet (see config.sh / config.R)" >&2;; esac; done
