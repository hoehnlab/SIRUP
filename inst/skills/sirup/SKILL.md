---
name: sirup
description: Run SIRUP (Spatial Immune Receptor Utilities Package) in R. SIRUP analyzes spatial B cell receptor (BCR) and T cell receptor (TCR) sequencing together with spatial gene expression from the same tissue section. Use this skill when the user has a Seurat spatial object and an AIRR receptor table, works with a SIRUP object (createSIRUP, updateReceptors, readSIRUP), or asks about clones on tissue, clonal dominance or diversity in regions, lineage trees on tissue, spatial statistics of genes, cell types, or clones, or processing spatial long-read receptor data.
---

# Use SIRUP

SIRUP joins three kinds of data in one R object. These are spatial gene expression (a Seurat object), receptor sequences (an AIRR table), and tissue regions that the user draws as polygons. SIRUP then measures clones, regions, trees, and spatial patterns. SIRUP does not call clones or build trees itself. It runs Immcantation tools (scoper, dowser, shazam) on the data inside the object and adds the spatial part.

## Start here

1. **Find out what the user has.** There are three cases.
   - Demultiplexed long reads with spatial barcodes: start with the processing workflow (`references/workflows.md`, recipe F).
   - A Seurat spatial object and an AIRR table: build a SIRUP object (recipe A).
   - A saved SIRUP object: load it with `readSIRUP()` and go to the analysis the user wants.
2. **Check the setup.** SIRUP needs R 4.1 or later. Run `library(SIRUP)` and `packageVersion("SIRUP")`. If SIRUP is missing, install it with BiocManager: `BiocManager::install("hoehnlab/SIRUP")`. This installs only the required packages. The standard workflow also needs the suggested packages `airr`, `Seurat`, and `scoper`. Check each package with `requireNamespace("name", quietly = TRUE)` before you use it, and see "Packages by task" below.
3. **Read the reference that matches the task.** Do not guess function names or arguments. Read the usage in `references/functions-*.md`. Start from `references/function-index.md`.

## The SIRUP object

A SIRUP object is a list with five parts, in this order:

| Part | Holds | Read it with |
|---|---|---|
| `spatial` | The Seurat object | `getSpatial(x)` |
| `receptors` | The AIRR table, as a tibble | `getReceptors(x)` |
| `metadata` | Column names and image settings | `x$metadata` |
| `polygons` | Tissue regions, one row for each polygon | `getPolygons(x)` |
| `analyses` | Results that you store | `getAnalysis(x, name)` |

Follow these rules:

- **Change the object only through SIRUP functions.** Use `updateReceptors()`, `updateSpatial()`, `updateAnalysis()`, `setAnalysis()`, and `addPolygons()`. Each one validates the object before it returns. Assignment with `$<-` or `[[<-` skips validation.
- **Set the column roles once.** `createSIRUP()` stores which receptor column holds the spatial ID, the clone, and the locus. Later functions find those columns without arguments.
- **Set the image settings once.** Pass `image_scale` (`"lowres"` or `"hires"`) and `swap_coords` to `createSIRUP()`. Plots and spatial statistics read them from the object.
- **Store results in the object.** Use `setAnalysis(x, "name", result)`. Save the whole object with `saveSIRUP()`.

`references/object.md` has the details: field names, validation rules, polygon columns, and multi-sample objects.

## The standard workflow

This code builds an object and calls clones. It is the core of every analysis. Replace the file names with the user's files.

```r
library(SIRUP)

spatial <- readRDS("spatial.rds")                          # Seurat object
bcrs <- airr::read_rearrangement("bcr.tsv.gz")             # AIRR table

# Copy the spatial barcode out of sequence_id. Match the Seurat spot names.
bcrs <- extractSpatialIds(
  bcrs, id = "sequence_id", spatial_id = "barcode_10x",
  regex = ",([ACGT]{16})$", suffix = "-1"
)
bcrs <- filterAirr(
  bcrs, locus_val = "IGH", min_v_length = 250, min_j_length = 30,
  min_junction_length = 10, productive_only = TRUE,
  min_consensus_count = NULL, spatial_id = "barcode_10x"
)

sirup <- createSIRUP(
  spatial = spatial, receptors = bcrs, spatial_id = "barcode_10x",
  image_scale = "lowres", swap_coords = FALSE
)
sirup <- updateReceptors(
  sirup, scoper::hierarchicalClones,
  threshold = 0.1, summarize_clones = FALSE, nproc = 1
)

validateSIRUP(sirup, verbose = TRUE)
summary(sirup)
```

The barcode pattern and the filter thresholds above come from the SIRUP paper for 10x Visium data. Check the user's `sequence_id` format before you reuse the pattern. Use `scoper::identicalClones` for TCRs.

Next, choose the analysis by the user's question. `references/choose-function.md` maps questions to functions. `references/workflows.md` has a full recipe for each analysis.

## Rules for agents

- **Check that spatial IDs match.** Run `mean(bcrs$barcode_10x %in% colnames(spatial))`. A value near 0 means the barcode pattern or the `suffix` is wrong. Fix it before you continue.
- **Do not run `selectSpatialPolygons()`.** It opens an interactive app that needs mouse clicks. Ask the user to draw the polygons and save a CSV. Then load the file with `addPolygons(sirup, polygons_path = "polygons.csv")`. The CSV needs the columns `x`, `y`, and `polygon` (or `polygon_id`).
- **Use the right clone caller.** scoper stops with an error if a junction has a character other than A, C, G, T, N, or `?`. It also stops if the table has a `clone_id` column. Remove that column before you call clones again. scoper drops sequences with an `N` and prints a warning. Count the receptors before and after, and report the change to the user.
- **Set a seed before tree building.** Tree building is random. Run `set.seed(1)` before `getSpatialTrees()`, so the user gets the same trees on each run.
- **Tell the user when a step needs outside data or tools.** Germline reconstruction needs the IMGT reference folder (the user sets `SIRUP_IMGT_DIR`). Read processing needs IgBLAST, pRESTO, vsearch, and MUSCLE. Do not install or download these without the user's approval.
- **Ask before long jobs.** Read processing and tree building for large clones can run for hours. Offer to write a script, and let the user run it on a compute node.
- **Use imputed expression for plots only.** The `alra` assay in the example data gives smoother maps. Use the `Spatial` assay for statistics.
- **Use `collapseIdenticalNeighbors()` on Visium HD data only.** On standard Visium, neighboring spots with the same sequence are more likely two related cells than spatial bleed.
- **Report what you ran.** Print the receptor count and the clone count after each change to the object. Look at each plot before you describe it. State which checks you did not run.

## Check your work

After each step, check these points:

1. `validateSIRUP(x, verbose = TRUE)` returns without an error.
2. `summary(x)` shows the numbers that you expect: spatial IDs, receptor records, matched IDs, stored analyses, and polygons.
3. The receptor table has a `clone_id` column after clone calling. Run `dplyr::n_distinct(getReceptors(x)$clone_id)`.
4. A map shows spots on the tissue image. If the spots sit rotated or mirrored over the image, set `swap_coords = TRUE` when you create the object.
5. A region summary has one row for each polygon that you expect.

## Packages by task

`BiocManager::install()` does not install suggested packages. Install the ones that the task needs.

| Task | Package | Installed with SIRUP? |
|---|---|---|
| Read an AIRR table with `airr::read_rearrangement()` | airr | No (suggested) |
| Read, edit, or plot Seurat data (`Seurat::Load10X_Spatial()`, `Seurat::AddMetaData()`) | Seurat | No (suggested). SeuratObject is installed. |
| Call clones | scoper | No (suggested) |
| Reconstruct germlines, build trees | dowser | Yes. The IMGT reference folder is also needed. |
| Count mutations, find chimeras | shazam | No (suggested) |
| Moran's I, Lee's L, join counts | spdep | Yes |
| Ripley's cross K | spatstat.explore, spatstat.geom, spatstat.random | No (suggested) |
| 3D clone density | rgl | No (suggested) |
| Read a Visium HD barcode mapping from a parquet file | arrow | No (suggested) |

## Where to read more

- `references/object.md`: the parts of the object, field names, validation rules, and the update functions.
- `references/workflows.md`: step-by-step recipes with code.
- `references/choose-function.md`: which function answers which question.
- `references/function-index.md` and `references/functions-*.md`: every function, with its usage and arguments. The build generates these files.
- `references/tutorials.md`: the tutorials and guides on the documentation site.
- The documentation site is at https://sirup.readthedocs.io. The files `llms.json` and `llms-full.json` in the repository hold the whole documentation in machine-readable form.
