# SIRUP

**S**patial **I**mmune **R**eceptor **U**tilities **P**ackage

SIRUP is an R package for spatial adaptive immune receptor repertoire
sequencing (spAIRR-seq). It works with B cell receptor (BCR) and T cell
receptor (TCR) sequences that carry a spatial barcode, together with spatial
gene expression from the same tissue section.

SIRUP covers three stages:

1. **Process.** Turn demultiplexed long reads into consensus receptor
   sequences. IgBLAST and pRESTO build an AIRR rearrangement table with one
   spatial ID per sequence.
2. **Integrate.** Put the receptor table and the spatial Seurat object into
   one validated SIRUP object. One object can hold several tissue sections,
   bulk AIRR-seq, and blood.
3. **Analyze and visualize.** Draw tissue regions, measure clonal dominance
   and diversity, test spatial patterns, and show lineage trees on the tissue.

SIRUP is in active development. Function names and arguments can change
before the first release.

## Installation

Install SIRUP from GitHub:

```r
if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
remotes::install_github("hjm2217/SIRUP")
```

Most analyses also use Seurat, packages from the
[Immcantation](https://immcantation.readthedocs.io) framework, and spatial
statistics packages. These are suggested, not required. Install the ones that
you need:

```r
install.packages(c(
  "Seurat", "airr", "alakazam", "scoper", "shazam", "dowser",
  "spdep", "spatstat.explore", "spatstat.geom", "spatstat.random",
  "cowplot", "patchwork"
))
BiocManager::install("ggtree")
```

### Tools for read processing

You need these command-line tools only to process raw reads. You do not need
them to analyze an AIRR table that you already have.

- [IgBLAST](https://ncbi.github.io/igblast/) (`igblastn`) and IMGT germline
  references
- [pRESTO and Change-O](https://immcantation.readthedocs.io)
  (`pip install presto changeo`)
- [vsearch](https://github.com/torognes/vsearch) and
  [MUSCLE](https://drive5.com/muscle/)

The [Immcantation Docker image](https://immcantation.readthedocs.io/en/stable/docker/intro.html)
contains all of them.

Demultiplexing happens before SIRUP, in the separate
[rrad](https://github.com/indianewok/rrad) package. rrad finds and corrects
the spatial barcode and UMI of each read, and writes a FASTA file that SIRUP
reads.

## Quick start

SIRUP includes a small spAIRR-seq data set: two Visium sections of a
myasthenia gravis thymus that contains ectopic germinal centers. This example
loads one section, calls clones, adds the germinal-center outlines, and
measures clonal dominance in each one.

```r
library(SIRUP)
extdata <- system.file("extdata", package = "SIRUP")

# Spatial gene expression (Seurat) and BCR consensus sequences (AIRR)
thymus <- readRDS(file.path(extdata, "mg_thymus_slice2_spatial.rds"))
bcrs <- airr::read_rearrangement(file.path(extdata, "mg_thymus_slice2_bcr.tsv.gz"))

# Move the spot barcode out of sequence_id, and keep high-quality sequences
bcrs <- extractSpatialIds(
  bcrs,
  id = "sequence_id",
  spatial_id = "barcode_10x",
  regex = ",([ACGT]{16})$",
  suffix = "-1"
)
bcrs <- filterAirr(
  bcrs,
  locus_val = "IGH",
  min_v_length = 250,
  min_j_length = 30,
  min_junction_length = 10,
  productive_only = TRUE,
  min_consensus_count = NULL,
  spatial_id = "barcode_10x"
)
bcrs <- bcrs[!grepl("[^ACGT]", bcrs$junction), ]

# Pair the receptors with the tissue
sirup <- createSIRUP(
  spatial = thymus,
  receptors = bcrs,
  spatial_id = "barcode_10x",
  image_scale = "lowres",
  swap_coords = FALSE
)

# Run any function on the receptor table, here clone calling with scoper
sirup <- updateReceptors(
  sirup,
  scoper::hierarchicalClones,
  threshold = 0.1,
  summarize_clones = FALSE
)

# Add the germinal-center outlines and plot them on the H&E image
sirup <- addPolygons(
  sirup,
  polygons_path = file.path(extdata, "mg_thymus_slice2_gc_polygons.csv")
)
plotSpatialPolygons(sirup)

# Fraction of each germinal center's BCRs that belongs to each clone
dominance <- spatialCloneDominance(object = sirup)
plotSpatialPolygonPies(dominance = dominance, object = sirup, top_n = 5)$page_1

# Keep the result with the data, and save everything to one file
sirup <- setAnalysis(sirup, "gc_dominance", dominance)
saveSIRUP(sirup, "mg_thymus_slice2.sirup.rds")
```

## Process your own reads

`processSpatialSample()` runs the full pipeline on one sample. It starts from
the rrad FASTA, a read-level AIRR table, or a FASTA that is already filtered:

| Step | Tool | What it does |
|---|---|---|
| 1 | IgBLAST | Annotates each read with V(D)J genes and a constant region |
| 2 | `filterAirr()` | Removes reads that are not receptors or are too short |
| 3 | pRESTO `ClusterSets.py` | Groups the reads of each spot by sequence identity |
| 4 | pRESTO `AlignSets.py` | Aligns the reads in each group |
| 5 | pRESTO `BuildConsensus.py` | Builds one consensus sequence per group |
| 6 | IgBLAST and Change-O `MakeDb.py` | Annotates the consensus sequences and writes the AIRR table |

```r
result <- processSpatialSample(
  input = "prepared/sample-igblast.fasta",
  output_dir = "processed/bcr",
  receptor = "ig",            # "tr" for TCRs
  igblast_db_dir = "/path/to/igblast/database",
  igblast_base = "/path/to/igblast",
  germline_glob = "/path/to/germlines/imgt/human/vdj/imgt_human_IG*.fasta",
  identity = 1.0,             # 0.99 for TCRs
  consensus_freq = 0.51,
  nproc = 8
)
```

Use `umi_support = TRUE` to build one consensus per molecule first. Use
`processSpatialSamples()` to run several samples with the same settings. For
Visium HD, `readSpatialIdMapping()` maps 2 µm barcodes to 8 µm and 16 µm bins.
The processing vignette explains each option.

## Functions by task

| Task | Functions |
|---|---|
| Process reads | `processSpatialSample()`, `processSpatialSamples()`, `runIgblast()`, `parseRadSequenceIds()` |
| Clean AIRR tables | `extractSpatialIds()`, `filterAirr()`, `plotLengthQc()`, `addSequenceQualityMetrics()`, `collapseTopConsensusPerSpatialId()`, `collapseIdenticalNeighbors()`, `readSpatialIdMapping()` |
| Build and manage objects | `createSIRUP()`, `mergeSIRUP()`, `subsetSIRUP()`, `validateSIRUP()`, `spatialSamples()`, `saveSIRUP()`, `readSIRUP()` |
| Access and change data | `getSpatial()`, `getReceptors()`, `updateSpatial()`, `updateReceptors()`, `setAnalysis()`, `getAnalysis()`, `listAnalyses()` |
| Plot on the tissue | `spatialExpressionPlot()`, `spatialExpressionBlendPlot()`, `spatialDiscretePlot()`, `plotReceptors()`, `plotReceptorsExpression()` |
| Tissue regions | `selectSpatialPolygons()`, `addPolygons()`, `getPolygons()`, `plotSpatialPolygons()`, `withPolygons()`, `expandPolygons()` |
| Regional repertoires | `polygonReceptorSummary()`, `spatialCloneDominance()`, `plotSpatialPolygonPies()`, `polygonDiversity()`, `polygonAbundance()`, `plotSpatialCloneDensity3D()` |
| Spatial statistics | `spatialMoran()`, `spatialLocalMoran()`, `spatialLee()`, `spatialRipleyCrossK()`, `spatialCloneJoinCount()`, `spatialCloneSpread()` |
| Lineage trees | `formatSpatialClones()`, `getSpatialTrees()`, `plotSpatialCloneTreePairs()`, `spatialCloneTreePairsToPDF()`, `tipPatristicMatrix()` |
| Heavy-light chain pairing | `pairHeavyLightSpaced()`, `buildSpacedCloneTable()`, `spacedSpatialAffinity()` |

SIRUP does not wrap every AIRR tool. Use `updateReceptors()` to run functions
from scoper, dowser, shazam, or your own code on the stored receptor table.
SIRUP validates the object again after each change.

## Documentation

The vignettes use the example data and run from start to finish. Install with
`build_vignettes = TRUE` to read them in R:

```r
remotes::install_github("hjm2217/SIRUP", build_vignettes = TRUE)
browseVignettes("SIRUP")
```

| Vignette | Topic |
|---|---|
| `vignette("SIRUP")` | Get started: the SIRUP object, plots, and a first regional analysis |
| `vignette("processing")` | Reads to AIRR tables with IgBLAST and pRESTO |
| `vignette("regions")` | Tissue regions, clonal dominance, and diversity |
| `vignette("spatial-statistics")` | Moran's I, Lee's L, Ripley's cross K, and join counts |
| `vignette("clones-and-trees")` | Clones, germlines, and lineage trees on the tissue |
| `vignette("multisample")` | Several sections, bulk AIRR-seq, and blood |

More guides are in [`docs/readthedocs`](docs/readthedocs):

- [Preparing FASTA with rrad](docs/readthedocs/markdown/demultiplexing-with-rad.md)
- [Custom read layouts and whitelists](docs/readthedocs/markdown/custom-read-layouts-with-rad.md)
- [SPaCeD heavy-light chain pairing](docs/readthedocs/markdown/spaced-pairing.md)

## Citation

The SIRUP manuscript is in preparation. Until it is published, please cite
this repository.

## Authors

- Hunter J. Melton (maintainer)
- Chinmay Vaidya

## License

MIT. See [LICENSE](LICENSE).
