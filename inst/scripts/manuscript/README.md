# Manuscript Analysis Scripts

This directory holds figure-oriented scripts for the SIRAP Genome Research
spatial omics manuscript. Keep scripts reproducible, parameterized by sample
ID where possible, and write tabular outputs alongside plots so manuscript
claims can be checked without rerunning every figure.

Manuscript framing:

- SIRAP has two arms: processing spAIRR data and visualizing/analyzing spAIRR
  data.
- RAD is the companion demultiplexing/read-structure package that SIRAP wraps
  from R; avoid describing RAD as an internal SIRAP method.
- SPaCeD is its own spatial pairing method. SIRAP includes a
  re-implementation/integration of SPaCeD-style heavy/light pairing for use
  alongside the rest of the analysis workflow.

Suggested output root:

```r
manuscript_out_dir <- file.path(thymusDir, "SIRAP_manuscript")
```

Current script targets:

- `figure_01_workflow_and_inventory.R`: dataset inventory and representative receptor overlays.
- `figure_02_thymus_gc_clonality.R`: thymic GC polygons, clone dominance, clone localization, lineage examples.
- `figure_05_thymus_pbmc_connections.R`: thymus versus PBMC shared-clone summaries and representative lineage trees.
