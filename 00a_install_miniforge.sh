#!/bin/bash
# ============================================================
# Step 00a (once, on the login node, ~2 min):
# install your own up-to-date conda + mamba (Miniforge) in $HOME.
# MSI's "module load conda" ships mamba 0.5.1 (2020), which can no
# longer download today's conda-forge/bioconda package lists.
#
#   bash 00a_install_miniforge.sh
#
# After this, every pipeline script uses ~/miniforge3 automatically.
# Your old envs (e.g. pu_RNAseq) are untouched.
# ============================================================
set -euo pipefail
mkdir -p log
LOG="log/00a_install_miniforge_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "${LOG}") 2>&1

PREFIX="$HOME/miniforge3"
if [[ -x "${PREFIX}/bin/conda" ]]; then
    echo "Miniforge already installed at ${PREFIX}"
else
    cd "${TMPDIR:-/tmp}"
    curl -fsSL -o Miniforge3.sh \
        "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh"
    bash Miniforge3.sh -b -p "${PREFIX}"     # -b: no prompts, does NOT edit ~/.bashrc
    rm -f Miniforge3.sh
fi

source "${PREFIX}/etc/profile.d/conda.sh"
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority strict
echo "conda : $(conda --version)"
echo "mamba : $(mamba --version | head -1)"
echo "Done. Next: sbatch 00_setup_env.sh"
