Process reads into spatial AIRR tables
================

This vignette takes you from demultiplexed long reads to an AIRR table
in which each consensus receptor sequence has a spatial ID. The next
step, `createSIRUP()`, is in `vignette("SIRUP", package = "SIRUP")`.

The pipeline calls external programs (IgBLAST, pRESTO, vsearch, and
MUSCLE). A vignette cannot run these programs on every computer, so the
chunks that call them show code only. The other chunks run on small
samples of the MG thymus data that come with SIRUP.

## The pipeline

Spatial AIRR-seq libraries are long-read sequencing libraries of
receptor transcripts. Each read carries a spatial barcode, which
identifies its spot. In most protocols it also carries a UMI, which
identifies its original molecule. The processing has two parts:

1.  **Demultiplexing**, in the separate
    [rrad](https://github.com/indianewok/rrad) package. rrad finds the
    barcode and UMI in each read, corrects them against the whitelist,
    and writes them into the read header.
2.  **Consensus building**, in SIRUP. `processSpatialSample()` runs
    these steps:

| Step | Tool                             | What it does                                                 |
| ---- | -------------------------------- | ------------------------------------------------------------ |
| 1    | IgBLAST                          | Annotates each read with V(D)J genes and a constant region   |
| 2    | `filterAirr()`                   | Removes reads that are not receptors or are too short to use |
| 3    | pRESTO `ClusterSets.py`          | Groups the reads of each spot by sequence identity (vsearch) |
| 4    | pRESTO `AlignSets.py`            | Aligns the reads in each group (MUSCLE)                      |
| 5    | pRESTO `BuildConsensus.py`       | Builds one consensus sequence per group                      |
| 6    | IgBLAST and Change-O `MakeDb.py` | Annotates the consensus sequences and writes the AIRR table  |

Long-read sequencing has a high error rate. A consensus of several reads
removes most random errors. The spot barcode stays in the sequence
header through every step, so each consensus sequence keeps its spot.

## Install the external tools

You need these programs on your computer or cluster:

  - **IgBLAST** (`igblastn`) from NCBI, with its `internal_data` and
    `optional_file` folders.
  - **IMGT germline references**, both as IgBLAST databases and as FASTA
    files for Change-O. The Immcantation scripts `fetch_imgtdb.sh` and
    `imgt2igblast.sh` download and convert them.
  - **pRESTO** and **Change-O** (`pip install presto changeo`).
  - **vsearch**, which `ClusterSets.py` calls.
  - **MUSCLE**, which `AlignSets.py` calls. The SIRUP paper used MUSCLE
    3.8.31.

The [Immcantation Docker image](https://immcantation.readthedocs.io)
contains all of them. If you use it, run R inside the container.

Make sure R can find the programs:

``` r
Sys.which(c(
  "igblastn", "ClusterSets.py", "ParseHeaders.py", "AlignSets.py",
  "BuildConsensus.py", "AssignGenes.py", "MakeDb.py", "vsearch"
))
```

An empty string means R cannot find that program. Add its folder to
`PATH`, or give its full path to `processSpatialSample()`
(`igblast_base`, `igblast_exec`, `muscle_exec`).

## Demultiplex with rrad

rrad is a separate package. Install it from GitHub:

``` r
remotes::install_github("indianewok/rrad")
```

`rrad::rad_demux()` finds and corrects the barcodes and UMIs.
`rrad::rad_reformat()` then writes a FASTA file that IgBLAST can read.
Choose the `layout` that matches your library:

``` r
demux <- rrad::rad_demux(
  layout = "five_prime",
  fastq = "reads.fq.gz",
  out_dir = "prepared",
  output = "sample"
)

rrad::rad_reformat(
  demux,
  output_fastq = "prepared/sample-igblast.fasta",
  to_fasta = TRUE,
  header_format = "presto"
)
```

Use `header_format = "presto"`. It puts the barcode and UMI into the
read identifier, which IgBLAST keeps:

``` text
>read1|BARCODE=AAACAAGTATCTCCCA|UMI=TCACATTCTGTA
```

`parseRadSequenceIds()` reads these tags into columns.
`processSpatialSample()` calls it for you, but you can use it to check a
file:

``` r
library(SIRUP)

headers <- data.frame(sequence_id = c(
  "read1|BARCODE=AAACAAGTATCTCCCA|UMI=TCACATTCTGTA",
  "read2|BARCODE=AAACACCAATAACTGC|UMI=AAAAAGGCTGGT"
))
parseRadSequenceIds(headers)
#>   sequence_id          barcode          umi
#> 1       read1 AAACAAGTATCTCCCA TCACATTCTGTA
#> 2       read2 AAACACCAATAACTGC AAAAAGGCTGGT
```

For Visium HD, rrad writes a coordinate barcode such as
`s_002um_00001_00002-1` in place of the 16-base barcode. The [rrad
documentation](https://github.com/indianewok/rrad) explains the Visium
HD layout.

## Look at the read-level data

SIRUP includes 2,000 read-level IgBLAST rows from MG thymus section 2.
The RAD run for these reads wrote the barcode and UMI as the columns
`cell_barcode` and `umi`, not in the header. We enriched this sample for
usable reads: 750 of the rows pass the IGH filter below. In the full
section, only 62,738 of 15.6 million reads (0.4%) pass.

``` r
library(dplyr)

extdata <- system.file("extdata", package = "SIRUP")
reads <- airr::read_rearrangement(
  file.path(extdata, "mg_thymus_slice2_reads.tsv.gz")
)

dim(reads)
#> [1] 2000  130
count(reads, locus, sort = TRUE)
#> # A tibble: 4 × 2
#>   locus     n
#>   <chr> <int>
#> 1 IGH    1259
#> 2 IGL     480
#> 3 IGK     205
#> 4 <NA>     56
```

Most reads are fragments, light chains, or other transcripts.
`filterAirr()` keeps the reads worth building consensus sequences from.
These settings match the SIRUP paper for IGH reads:

``` r
reads$barcode <- reads$cell_barcode

igh_reads <- filterAirr(
  reads,
  locus_val = "IGH",
  min_v_length = 100,
  min_j_length = 10,
  min_junction_length = 1,
  productive_only = FALSE,
  min_consensus_count = NULL,
  spatial_id = "barcode"
)
#> filterAirr kept 750 of 2000 rows (dropped 1250).

spatialReceptorSummary(igh_reads, spatial_id = "barcode")
#> # A tibble: 1 × 4
#>   n_sequences n_spatial_ids n_loci n_productive
#>         <int>         <int>  <int>        <int>
#> 1         750           451      1          307
```

Keep `productive_only = FALSE` at the read level. A single read can look
unproductive because of a sequencing error that the consensus step
corrects. Filter on productivity after the consensus step. (With
`umi_support = TRUE`, SIRUP builds each molecule consensus from the
productive reads of that molecule.)

`plotLengthQc()` shows the V, J, and junction lengths. Use it to choose
the thresholds:

``` r
plotLengthQc(igh_reads, bin_width = 10, title = "IGH reads, MG thymus section 2")
```

<img src="/dartfs/rc/lab/H/HoehnK/Hunter/spatial_bcr_tcr_package/SIRUP/docs/readthedocs/markdown/processing_files/figure-gfm/length-qc-1.jpeg" alt="" width="100%" />

## Run the consensus pipeline

`processSpatialSample()` runs every step from the table above. It prints
each command before it runs it. If a command fails, it stops and shows
the command.

You can start from three types of input:

| `input`                                                     | `input_type`              | First step that runs |
| ----------------------------------------------------------- | ------------------------- | -------------------- |
| FASTA or FASTQ from rrad                                    | `"auto"` or `"raw_reads"` | IgBLAST on each read |
| AIRR table (path or data frame) from a first IgBLAST pass   | `"auto"` or `"airr"`      | `filterAirr()`       |
| FASTA that is already filtered, with the spot in the header | `"spatial_fasta"`         | pRESTO clustering    |

### From rrad FASTA

``` r
result <- processSpatialSample(
  input = "prepared/sample-igblast.fasta",
  output_dir = "processed/bcr",
  receptor = "ig",
  organism = "human",
  igblast_db_dir = "~/share/igblast/database",
  igblast_base = "~/share/igblast",
  germline_glob = "~/share/germlines/imgt/human/vdj/imgt_human_IG*.fasta",
  muscle_exec = "~/bin/muscle3.8.31_i86linux64",
  identity = 1.0,
  consensus_freq = 0.51,
  consensus_min_count = 2,
  nproc = 8
)

result$db_pass
```

Use `receptor = "tr"` and the TR germline files for TCRs. Run BCRs and
TCRs as separate calls.

### With UMI support

Reads of one molecule can differ by sequencing errors. At `identity
= 1.0`, those reads fall into different clusters, and each cluster can
then fail the two-read minimum. With `umi_support = TRUE`, SIRUP first
groups the reads of each spot and UMI at a looser identity
(`umi_identity`, default 0.95). It builds one consensus per molecule
from the productive reads only. The spot-level step then combines the
molecules and sums their reads.

`umi_support` needs an AIRR table with `umi` and `productive` columns,
so give it the output of the first IgBLAST pass. This is the call we
used for the MG thymus BCRs in the SIRUP paper:

``` r
result <- processSpatialSample(
  input = igblast_reads,           # read-level AIRR data frame
  input_type = "airr",
  output_dir = "processed/kh02_bcr",
  receptor = "ig",
  identity = 1.0,
  consensus_freq = 0.51,
  cluster_field = "barcode",
  germline_glob = "~/share/germlines/imgt/human/vdj/imgt_human_IG*.fasta",
  igblast_base = "~/share/igblast",
  muscle_exec = "~/bin/muscle3.8.31_i86linux64",
  nproc = 8,
  umi_support = TRUE,
  umi_field = "umi",
  umi_identity = 0.95,
  min_umi_support = 2,
  consensus_copy_fields = "C_CALL",
  consensus_copy_action = "majority",
  filter_args = list(
    locus_val = "IGH",
    min_v_length = 100,
    min_j_length = 10,
    productive_only = FALSE,
    complete_vdj = FALSE,
    min_consensus_count = NULL,
    spatial_id = "barcode",
    add_spatial_sequence_id = FALSE
  )
)

result$db_pass_umi
```

With `umi_support = TRUE`, the result has two AIRR tables. `db_pass`
holds every consensus sequence. `db_pass_umi` holds the sequences with
at least `min_umi_support` reads. You can change the threshold later
without running the pipeline again.

### Keep the constant region

The consensus covers the V(D)J region only, so the second IgBLAST pass
cannot call the constant region. To keep the isotype, put the read-level
`c_call` into the FASTA header. Then give `consensus_copy_fields =
"C_CALL"` and `consensus_copy_action = "majority"`. Each consensus
sequence gets the constant region of most of its reads, in `c_call`, and
the fraction of reads that agree, in `c_call_freq`.

### Many samples

`processSpatialSamples()` runs `processSpatialSample()` on each input,
with the same settings:

``` r
processSpatialSamples(
  inputs = c("prepared/slice1.fasta", "prepared/slice2.fasta"),
  output_dirs = c("processed/slice1", "processed/slice2"),
  receptor = "ig",
  igblast_db_dir = "~/share/igblast/database",
  germline_glob = "~/share/germlines/imgt/human/vdj/imgt_human_IG*.fasta",
  nproc = 8
)
```

On a cluster, submit one job per sample and chain instead. Each job then
calls `processSpatialSample()` once.

## Choose the parameters

Two parameters have the largest effect on the result:

  - `identity` is the similarity that two reads need to go into the same
    cluster in a spot.
  - `consensus_freq` is the fraction of reads that must agree at a
    position. If fewer reads agree, the consensus gets an `N` at that
    position.

A standard Visium spot is 55 µm wide and holds many cells. Two BCRs in
the same spot that differ by one base can come from two related B cells
after somatic hypermutation. For BCRs, use `identity = 1.0`, so that
SIRUP does not merge them. TCRs do not mutate, so a one-base difference
in a spot is more likely a sequencing error. For TCRs, use `identity
= 0.99`. In the SIRUP paper, `consensus_freq = 0.51` gave the most
high-quality sequences for both chains.

To compare settings, run `processSpatialSample()` with each value.
Clustering and alignment depend on `identity` only, so with
`skip_existing = TRUE` you can try several `consensus_freq` values and
reuse the clustering.

## Filter the consensus table

The pipeline writes the consensus `sequence_id` as
`<cluster>,<barcode>`. The example data include the consensus table of
MG thymus section 2. It is the `db_pass_umi` table that
`processSpatialSample()` wrote for the SIRUP paper:

``` r
consensus <- airr::read_rearrangement(
  file.path(extdata, "mg_thymus_slice2_bcr.tsv.gz")
)
head(consensus$sequence_id)
#> [1] "9,TGATTCCCGGTTACCT"  "11,TGATTCCCGGTTACCT" "2,TGATTCCCGGTTACCT" 
#> [4] "6,TGATTCCCGGTTACCT"  "5,CCTACATTCACAGACG"  "3,CCTACATTCACAGACG"
```

`extractSpatialIds()` copies the part of an ID that a regular expression
matches into a new column. The `suffix` argument adds `-1`, so the
barcode matches the spot names that Space Ranger writes:

``` r
consensus <- extractSpatialIds(
  consensus,
  id = "sequence_id",
  spatial_id = "barcode_10x",
  regex = ",([ACGT]{16})$",
  suffix = "-1"
)
head(consensus$barcode_10x)
#> [1] "TGATTCCCGGTTACCT-1" "TGATTCCCGGTTACCT-1" "TGATTCCCGGTTACCT-1"
#> [4] "TGATTCCCGGTTACCT-1" "CCTACATTCACAGACG-1" "CCTACATTCACAGACG-1"
```

`addSequenceQualityMetrics()` adds three columns that help you set
thresholds: `non_N_nucleotide_prop` (the fraction of called bases),
`max_gap_length` (the longest gap in the alignment), and
`opening_gap_length`:

``` r
consensus <- addSequenceQualityMetrics(
  consensus,
  sequence = "sequence",
  alignment = "sequence_alignment"
)
summary(consensus[, c("non_N_nucleotide_prop", "max_gap_length", "consensus_count")])
#>  non_N_nucleotide_prop max_gap_length   consensus_count 
#>  Min.   :0.9510        Min.   :  6.00   Min.   : 1.000  
#>  1st Qu.:1.0000        1st Qu.: 12.00   1st Qu.: 1.000  
#>  Median :1.0000        Median : 12.00   Median : 1.000  
#>  Mean   :0.9982        Mean   : 63.68   Mean   : 1.441  
#>  3rd Qu.:1.0000        3rd Qu.:131.00   3rd Qu.: 1.000  
#>  Max.   :1.0000        Max.   :220.00   Max.   :35.000
```

Then filter. These settings match the definition of a quality sequence
that we used to tune the pipeline in the SIRUP paper: a productive
rearrangement, a V alignment of at least 250 nucleotides, a J of at
least 30, a junction of at least 10, and no gap longer than 15
nucleotides.

``` r
bcrs <- filterAirr(
  consensus,
  locus_val = "IGH",
  min_v_length = 250,
  min_j_length = 30,
  min_junction_length = 10,
  productive_only = TRUE,
  min_consensus_count = NULL,
  max_gap_length = 15,
  spatial_id = "barcode_10x"
)
#> filterAirr kept 2078 of 3801 rows (dropped 1723).

nrow(consensus)
#> [1] 3801
nrow(bcrs)
#> [1] 2078
```

This table is `db_pass_umi`, so every row already has at least two
reads. In UMI mode, `consensus_count` counts molecules, not reads, so
set `min_consensus_count = NULL` and use `min_umi_support` for the read
floor. Without UMI support, `consensus_count` counts reads. Then use
`min_consensus_count = 2`.

`max_gap_length` depends on the locus. IMGT numbering pads TRBV
alignments with longer gaps than IGHV alignments, so a limit of 15
removes valid TRB sequences. The SIRUP paper used 15 for IGH and 24 for
TRB.

The table is now ready for `createSIRUP()`.

## Visium HD

Visium HD reads carry a 2 µm barcode, such as `s_002um_00001_00002-1`.
Most analyses use 8 µm or 16 µm bins. Space Ranger writes a file that
maps each 2 µm barcode to its larger bins (`barcode_mappings.parquet`).
`readSpatialIdMapping()` reads it. Parquet files need the arrow package.
You can also save the file as CSV once and read that.

``` r
# bcrs$barcode_2um holds the 2 µm barcode of each consensus sequence.
mapping <- readSpatialIdMapping("outs/barcode_mappings.parquet")
mean(bcrs$barcode_2um %in% mapping$square_002um)

bcrs <- bcrs |>
  left_join(
    select(mapping, square_002um, square_008um, square_016um),
    by = c("barcode_2um" = "square_002um")
  ) |>
  rename(barcode_8um = square_008um, barcode_16um = square_016um)
```

Check the match rate before the join. If few 2 µm barcodes match,
compare a few values from each table. The two tables must write the row
and column of the bin in the same order (`s_002um_<row>_<col>-1`).

Load the 8 µm bins with `Seurat::Load10X_Spatial(bin.size = 8)`, and
give `spatial_id = "barcode_8um"` to `createSIRUP()`.

`collapseTopConsensusPerSpatialId()` keeps the consensus sequence with
the most reads in each spatial ID and locus. Use it when you expect one
cell for each spatial ID, for example at single-cell resolution.

### Identical sequences in neighboring spots

Transcripts can diffuse a short distance before capture. The same
molecule can then appear in neighboring spots. On Visium HD, where 2 µm
bins sit next to each other, this spatial bleed is common.
`collapseIdenticalNeighbors()` finds groups of neighboring spots that
hold the same exact sequence. It keeps one spot per group: the spot with
the highest consensus count, or the spot nearest the group center if
there is a tie.

The function also merges rows with the same sequence in the same spot
and adds their counts. It keeps rows whose spot is not in the Seurat
object, such as bulk rows, without changes.

**Use it on Visium HD data only.** Standard Visium spots are 55 µm wide,
100 µm apart, and hold several cells each. Two neighboring Visium spots
with the same sequence are more likely two related cells than spatial
bleed.

The function estimates the neighbor distance from the spot coordinates.
Set `num_neighbors_non_boundary` to the number of neighbors that a spot
in the middle of the tissue has: 6 on the standard Visium hexagonal
grid, 8 on the Visium HD square grid. The example below runs it on the
Visium thymus data only to show the call and the output:

``` r
thymus <- readRDS(file.path(extdata, "mg_thymus_slice2_spatial.rds"))
sirup <- createSIRUP(
  spatial = thymus,
  receptors = bcrs,
  spatial_id = "barcode_10x",
  image_scale = "lowres",
  swap_coords = FALSE
)

collapsed <- collapseIdenticalNeighbors(
  sirup,
  sequence = "sequence_alignment",
  num_neighbors_non_boundary = 6
)
#> ℹ Ensuring one row per spatial_id/sequence...
#> ✔ Ensuring one row per spatial_id/sequence... [1.6s]
#> 
#> ℹ Estimating distance threshold from full tissue positions...
#> ✔ Estimating distance threshold from full tissue positions... [2.8s]
#> 
#> ℹ Distance threshold: 10.1179
#> ℹ Filtering to spots with observed sequences...
#> ✔ Filtering to spots with observed sequences... [6ms]
#> 
#> ℹ Sequence-positive spots: 802 (of 3,203 total)
#> ℹ Building adjacency matrix over spots with observed sequences...
#> ✔ Building adjacency matrix over spots with observed sequences... [9ms]
#> 
#> ℹ Collapsing components to representative spots...
#> ✔ Collapsing components to representative spots... [189ms]
#> 
#> ℹ Filtering data to representative spots...
#> ✔ Filtering data to representative spots... [10ms]
#> 
#> ✔ Done! Retained 1,951 rows from 2,076 (125 bleedover sequences removed)

c(before = nrow(getReceptors(sirup)), after = nrow(getReceptors(collapsed)))
#> before  after 
#>   2078   1951
```

For Visium HD, give the 8 µm object and `num_neighbors_non_boundary
= 8`.

## Session information

``` r
sessionInfo()
#> R version 4.5.1 (2025-06-13)
#> Platform: x86_64-pc-linux-gnu
#> Running under: Red Hat Enterprise Linux 8.10 (Ootpa)
#> 
#> Matrix products: default
#> BLAS/LAPACK: /usr/lib64/libopenblasp-r0.3.15.so;  LAPACK version 3.9.0
#> 
#> locale:
#>  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
#>  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
#>  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
#>  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
#>  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
#> [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
#> 
#> time zone: America/New_York
#> tzcode source: system (glibc)
#> 
#> attached base packages:
#> [1] stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> other attached packages:
#> [1] scoper_1.5.0     ggplot2_4.0.3    dplyr_1.2.1      SIRUP_0.0.0.9001
#> [5] testthat_3.3.2  
#> 
#> loaded via a namespace (and not attached):
#>   [1] fs_2.1.0                    matrixStats_1.5.0          
#>   [3] spatstat.sparse_3.2-0       bitops_1.1-0               
#>   [5] devtools_2.5.2              httr_1.4.8                 
#>   [7] RColorBrewer_1.1-3          doParallel_1.0.17          
#>   [9] tools_4.5.1                 sctransform_0.4.3          
#>  [11] utf8_1.2.6                  R6_2.6.1                   
#>  [13] lazyeval_0.2.3              uwot_0.2.5                 
#>  [15] withr_3.0.3                 sp_2.2-3                   
#>  [17] prettyunits_1.2.0           gridExtra_2.3.1            
#>  [19] progressr_1.0.0             cli_3.6.6                  
#>  [21] Biobase_2.70.0              phylotate_1.3              
#>  [23] spatstat.explore_3.8-1      fastDummies_1.7.6          
#>  [25] sass_0.4.10                 labeling_0.4.3             
#>  [27] diptest_0.77-2              Seurat_5.5.1               
#>  [29] S7_0.2.2                    spatstat.data_3.1-9        
#>  [31] readr_2.2.0                 ggridges_0.5.7             
#>  [33] pbapply_1.7-4               Rsamtools_2.26.0           
#>  [35] systemfonts_1.3.2           yulab.utils_0.2.5          
#>  [37] dichromat_2.0-1             kedd_1.0.4                 
#>  [39] parallelly_1.48.0           sessioninfo_1.2.4          
#>  [41] maps_3.4.3                  pals_1.10                  
#>  [43] rstudioapi_0.19.0           generics_0.1.4             
#>  [45] gridGraphics_0.5-1          ica_1.0-3                  
#>  [47] spatstat.random_3.5-0       vroom_1.7.1                
#>  [49] Matrix_1.7-3                S4Vectors_0.48.1           
#>  [51] abind_1.4-8                 lifecycle_1.0.5            
#>  [53] yaml_2.3.12                 SummarizedExperiment_1.40.0
#>  [55] SparseArray_1.8.0           Rtsne_0.17                 
#>  [57] grid_4.5.1                  alakazam_1.4.3             
#>  [59] promises_1.5.0              crayon_1.5.3               
#>  [61] pwalign_1.6.0               miniUI_0.1.2               
#>  [63] lattice_0.22-7              cowplot_1.2.0              
#>  [65] cigarillo_1.0.0             mapproj_1.2.12             
#>  [67] pillar_1.11.1               knitr_1.51                 
#>  [69] GenomicRanges_1.62.1        future.apply_1.20.2        
#>  [71] codetools_0.2-20            fastmatch_1.1-8            
#>  [73] glue_1.8.1                  ggiraph_0.9.6              
#>  [75] ggfun_0.2.1                 spatstat.univar_3.2-0      
#>  [77] fontLiberation_0.1.0        data.table_1.18.4          
#>  [79] vctrs_0.7.3                 png_0.1-9                  
#>  [81] treeio_1.34.0               spam_2.11-4                
#>  [83] gtable_0.3.6                cachem_1.1.0               
#>  [85] xfun_0.60                   S4Arrays_1.10.1            
#>  [87] mime_0.13                   Seqinfo_1.0.0              
#>  [89] coda_0.19-4.1               survival_3.8-9             
#>  [91] iterators_1.0.14            ellipsis_0.3.3             
#>  [93] fitdistrplus_1.2-6          ROCR_1.0-12                
#>  [95] nlme_3.1-168                ggtree_4.0.5               
#>  [97] usethis_3.2.1               bit64_4.8.2                
#>  [99] fontquiver_0.2.1            progress_1.2.3             
#> [101] RcppAnnoy_0.0.23            rprojroot_2.1.1            
#> [103] bslib_0.12.0                irlba_2.3.7                
#> [105] airr_2.0.0                  KernSmooth_2.23-26         
#> [107] otel_0.2.0                  colorspace_2.1-3           
#> [109] BiocGenerics_0.56.0         ade4_1.7-24                
#> [111] phangorn_2.12.1             tidyselect_1.2.1           
#> [113] bit_4.6.0                   compiler_4.5.1             
#> [115] desc_1.4.3                  fontBitstreamVera_0.1.1    
#> [117] DelayedArray_0.36.1         plotly_4.12.1              
#> [119] scales_1.4.0                lmtest_0.9-40              
#> [121] quadprog_1.5-8              rappdirs_0.3.4             
#> [123] stringr_1.6.0               digest_0.6.39              
#> [125] goftest_1.2-3               spatstat.utils_3.2-4       
#> [127] dowser_2.5.1.999            rmarkdown_2.31             
#> [129] XVector_0.50.0              htmltools_0.5.9            
#> [131] pkgconfig_2.0.3             MatrixGenerics_1.22.0      
#> [133] fastmap_1.2.0               rlang_1.3.0                
#> [135] htmlwidgets_1.6.4           shiny_1.14.0               
#> [137] jquerylib_0.1.4             farver_2.1.2               
#> [139] zoo_1.9-0                   jsonlite_2.0.0             
#> [141] BiocParallel_1.38.0         magrittr_2.0.5             
#> [143] ggplotify_0.1.3             dotCall64_1.2              
#> [145] patchwork_1.3.2             Rcpp_1.1.2                 
#> [147] ape_5.8-1                   gdtools_0.5.1              
#> [149] reticulate_1.46.0           stringi_1.8.9              
#> [151] brio_1.1.5                  MASS_7.3-65                
#> [153] plyr_1.8.9                  pkgbuild_1.4.8             
#> [155] parallel_4.5.1              listenv_1.0.0              
#> [157] ggrepel_0.9.8               deldir_2.0-4               
#> [159] Biostrings_2.78.0           shazam_1.3.2               
#> [161] splines_4.5.1               tensor_1.5.1               
#> [163] hms_1.1.4                   igraph_2.3.3               
#> [165] fastcluster_1.3.0           spatstat.geom_3.8-1        
#> [167] seqinr_4.2-44               markdown_2.0               
#> [169] RcppHNSW_0.7.0              reshape2_1.4.5             
#> [171] stats4_4.5.1                pkgload_1.5.3              
#> [173] evaluate_1.0.5              SeuratObject_5.4.0         
#> [175] tzdb_0.5.0                  foreach_1.5.2              
#> [177] httpuv_1.6.17               RANN_2.6.3                 
#> [179] tidyr_1.3.2                 purrr_1.2.2                
#> [181] polyclip_1.10-7             future_1.75.0              
#> [183] scattermore_1.2             xtable_1.8-8               
#> [185] RSpectra_0.16-2             tidytree_0.4.8             
#> [187] later_1.4.8                 viridisLite_0.4.3          
#> [189] tibble_3.3.1                aplot_0.3.1                
#> [191] memoise_2.0.1               GenomicAlignments_1.46.0   
#> [193] IRanges_2.44.0              cluster_2.1.8.3            
#> [195] globals_0.19.1
```
