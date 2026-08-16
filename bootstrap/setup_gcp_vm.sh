#!/usr/bin/env bash
set -euo pipefail

environment_name="setd1b-rnaseq-workflow"
miniforge_dir="${HOME}/miniforge3"
installer_url="https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh"

if [[ ! -f bootstrap/environment.yaml || ! -f workflow/Snakefile ]]; then
    echo "Run this script from the repository root." >&2
    exit 1
fi

if [[ "$(uname -s)" != "Linux" || "$(uname -m)" != "x86_64" ]]; then
    echo "This setup script requires 64-bit Linux (x86_64)." >&2
    exit 1
fi

for required_command in curl bzip2; do
    if ! command -v "${required_command}" >/dev/null 2>&1; then
        echo "Missing ${required_command}. Install prerequisites first:" >&2
        echo "  sudo apt-get update && sudo apt-get install -y git curl bzip2 tmux make" >&2
        exit 1
    fi
done

memory_kib="$(awk '/^MemTotal:/ {print $2}' /proc/meminfo)"
disk_kib="$(df -Pk . | awk 'NR == 2 {print $4}')"
if (( memory_kib < 58000000 )); then
    echo "Warning: less than about 58 GB RAM is available; STAR indexing may fail." >&2
fi
if (( disk_kib < 400000000 )); then
    echo "Warning: less than about 400 GB disk is available; monitor disk usage closely." >&2
fi

if [[ ! -x "${miniforge_dir}/bin/conda" ]]; then
    installer_path="$(mktemp --suffix=.sh)"
    trap 'rm -f "${installer_path:-}"' EXIT
    echo "Installing Miniforge in ${miniforge_dir}..."
    curl --fail --location --retry 5 "${installer_url}" --output "${installer_path}"
    bash "${installer_path}" -b -p "${miniforge_dir}"
fi

source "${miniforge_dir}/etc/profile.d/conda.sh"

if conda env list | awk '{print $1}' | grep -Fxq "${environment_name}"; then
    echo "Updating existing ${environment_name} environment..."
    CONDA_CHANNEL_PRIORITY=strict conda env update \
        --name "${environment_name}" \
        --file bootstrap/environment.yaml \
        --prune
else
    echo "Creating ${environment_name} environment..."
    CONDA_CHANNEL_PRIORITY=strict conda env create --file bootstrap/environment.yaml
fi

conda run --name "${environment_name}" snakemake --version

cat <<'EOF'

Setup completed.

For this terminal, activate the workflow with:
  source "$HOME/miniforge3/etc/profile.d/conda.sh"
  conda activate setd1b-rnaseq-workflow

Then validate the workflow with:
  make gcp-dryrun
EOF
