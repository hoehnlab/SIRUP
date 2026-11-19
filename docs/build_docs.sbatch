#!/bin/bash
# Build the generated docs pages on Discovery. Run it from the package root:
#
#   sbatch docs/build_docs.sbatch
#
# To build only some vignettes (faster), set two variables:
#   sbatch --export=ALL,DOCS_STEPS=vignettes,DOCS_VIGNETTES=regions:multisample docs/build_docs.sbatch
# A partial build skips the nav, the extras and the final check.
#
# The job knits the 6 vignettes and converts man/*.Rd to docs/topics/*.md.
# Then it makes the nav, index and llms files. When it ends, check the pages and
# commit docs/, mkdocs.yml, llms.txt and llms-full.txt.
#
# Resources: measured on job 10924415 (2026-11-08). It took 2 min 33 s, used about
# 0.8 CPU cores, and peaked at 1.8 GB. The request below adds headroom. To check
# a later run, use:
#   sacct -j <jobid> --format=Elapsed,TotalCPU,MaxRSS
#SBATCH --job-name=sirup-docs
#SBATCH --account=hoehnlab-share
#SBATCH --partition=standard
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:30:00
#SBATCH --output=docs/build_docs_%j.log

set -euo pipefail

# sbatch can copy a stale module state from the login shell. That makes `module load`
# fail with "Loaded environment state is inconsistent". Clear the state first.
unset LOADEDMODULES _LMFILES_

# R 4.5.1 loads the personal library with Seurat. R 4.4.3 and 4.6.1 do not.
module load R/4.5.1

# The system pandoc is 2.0.6. The topic converter needs 2.11.2 or later.
# Use the pandoc 3.8.3 from the VisiumHD_Analysis conda environment.
export RSTUDIO_PANDOC=/dartfs-hpc/rc/home/j/f007p0j/.conda/envs/VisiumHD_Analysis/bin

# The tree chunks in the vignettes need the IMGT germline fasta files.
export SIRUP_IMGT_DIR=/dartfs/rc/lab/H/HoehnK/Hunter/germlines/imgt/human/vdj

# pRESTO and Change-O tools, and IgBLAST, in case a chunk calls them.
export PATH=/dartfs-hpc/rc/home/j/f007p0j/.conda/envs/DiscoveryPythonEnv/bin:$HOME/software/igblast/bin:$PATH

cd "${SLURM_SUBMIT_DIR:-.}"
test -f DESCRIPTION || { echo "Run this from the SIRUP package root." >&2; exit 1; }

# DOCS_STEPS is a space-separated list of steps (topics vignettes nav extras). Empty means all.
# shellcheck disable=SC2086
Rscript docs/build.R ${DOCS_STEPS:-}

# After a full build, also make a static copy of the site for the lab, until the package
# is public. The copy goes to ../docs_preview. A failure here does not fail the job.
if [ -z "${DOCS_STEPS:-}" ] && [ -z "${DOCS_VIGNETTES:-}" ]; then
  MKDOCS_PYTHON=/dartfs/rc/lab/H/HoehnK/Hunter/spatial_bcr_tcr_package/tmp/mkdocs-venv/bin/python \
    bash docs/build_preview.sh || echo "The preview build failed. The docs pages are built." >&2
fi
if [ -z "${DOCS_STEPS:-}" ]; then
  python3 docs/docs_tools.py check --no-build
fi
