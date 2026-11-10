Process demultiplexed reads
================

<div markdown="1"
style="background: #F7E5BC; padding: 12px; margin-bottom: 24px; line-height: 24px;">

<b>Time</b>: ~30 minutes. <b>Data</b>: example MG thymus read-level data
in SIRUP.

</div>

Here, we will take demultiplexed long-read sequences and process them
into an [AIRR
rearrangement](https://docs.airr-community.org/en/latest/datarep/rearrangements.html).

## Demultiplexing long-read sequences

SIRUP does not demultiplex long-read sequencing data on its own. At
present, users should demultiplex their data with

1.  RAD, the read-structure agnostic demultiplexer, which can be used
    [via the command line](https://github.com/indianewok/rad) or
    [natively in R](https://github.com/indianewok/rrad). SIRUP can then
    process the demultiplexed long-reads.
2.  [LongAIRR](https://github.com/AGImkeller/LongAIRR), which
    demultiplexes and processes long-read sequencing data.

If using LongAIRR, the output is an AIRR rearrangement that can be used
directly, without further processing in SIRUP. In this case, the rest of
this vignette is not needed.

## SIRUP processing pipeline

Spatial AIRR-seq libraries are sets of long-read sequences of B and T
cell receptor transcripts. Each read carries a spatial ID (or barcode),
which is a short DNA tag that identifies the spot where the molecule
came from. Many protocols also include a UMI (unique molecular
identifier), a tag that identifies the original molecule.

Processing in SIRUP is done with the `processSpatialSample()` function.
It consists of the following steps.

|  | Step | What it does |
|----|----|----|
| 1 | IgBLAST (`airrblast::assignGenes()`) | Annotates each read with V(D)J genes |
| 2 | `filterAirr()` | Removes reads that are not receptors or are too short to use |
| 3 | `ClusterSets.py` (pRESTO) | Groups the reads of each spot by sequence identity (using vsearch) |
| 4 | `AlignSets.py` (pRESTO) | Aligns the reads in each group (using MUSCLE) |
| 5 | `BuildConsensus.py` (pRESTO) | Builds one consensus sequence per group |
| 6 | `airrblast::assignGenes()` and `airrblast::addGaps()` | Annotates the consensus sequences, adds the IMGT gaps, and writes the AIRR table |

Long-read sequencing has a high error rate. A consensus of several reads
removes most random errors. The spot barcode stays in the sequence
header through every step, so each consensus sequence keeps its spot.

## Install and set up dependencies

### IgBLAST and IMGT germline references

See more detailed instructions in the [airrblast
documentation](https://github.com/hoehnlab/airrblast/tree/main).

Briefly, download and decompress
[IgBLAST](https://ftp.ncbi.nih.gov/blast/executables/igblast/release/LATEST/).

<div markdown="1"
style="background: #9e584a; padding: 12px; margin-bottom: 24px; line-height: 24px;">

<p style="margin: -12px -12px 12px; padding: 6px 12px; background: #8b3620; color: #ffffff; font-weight: 700;">

<b>WARNING</b>
</p>

There is a patch to the latest stable release of IgBLAST. Download the
patched .aux files
[here](https://ftp.ncbi.nih.gov/blast/executables/igblast/release/patch/)
(in `optional_file`). Then swap them into your `optional_file` directory
in your IgBLAST installation.

</div>

Then download the IMGT germline references and set up the IgBLAST
database with `airrblast`.

``` r
library(airrblast)
library(airr)

# name/location of the uncompressed IgBLAST directory
igblast_dir <- "ncbi-igblast-1.22.0"

# Download IMGT GENE-DB into directory "imgt"
# Do this only once
downloadIMGT("imgt", organism="human", verbose=1)

# Set up search databases for IgBlast using imgt database
# Do this only once
imgtToIgblast("imgt", "igblast", igblast_dir)
```

################ 

# PICK UP HERE

################ 

Install these programs on your computer or cluster:

- **IgBLAST** (`igblastn`) from NCBI, with its `internal_data` and
  `optional_file` folders.
- **IMGT germline references**, both as IgBLAST databases and as gapped
  FASTA files. The [airrblast](https://github.com/hoehnlab/airrblast)
  package downloads and converts them (see below).
- **pRESTO** (`pip install presto`).
- **vsearch**, which `ClusterSets.py` calls.
- **MUSCLE**, which `AlignSets.py` calls. The SIRUP paper used MUSCLE
  3.8.31.

> **Tip.** The [Immcantation Docker
> image](https://immcantation.readthedocs.io) contains most of these
> tools. Run R inside the container and install airrblast.

> **Note.** Paths such as `~/share/igblast` and
> `~/bin/muscle3.8.31_i86linux64` in this page are examples. Replace
> them with the paths on your system.

Check that R can find the programs. This chunk shows code only:

``` r
Sys.which(c(
  "ClusterSets.py", "ParseHeaders.py", "AlignSets.py", "BuildConsensus.py",
  "vsearch"
))
file.exists(file.path("~/share/igblast", "bin", "igblastn"))
```

An empty string means R cannot find that program. Add its folder to
`PATH`, or give its full path to `processSpatialSample()`
(`muscle_exec`). SIRUP runs IgBLAST through airrblast. It uses
`<igblast_base>/bin/igblastn`.

### Build the IgBLAST databases

SIRUP imports the airrblast package. Run this code once to download the
IMGT germlines and build the IgBLAST databases:

``` r
library(airrblast)
downloadIMGT("~/share/imgt", organism = "human")
imgtToIgblast("~/share/imgt", "~/share/igblast_refs", igblast = "~/share/igblast")
```

`imgtToIgblast()` writes the `database/` and `fasta/` folders into
`~/share/igblast_refs`. Give that folder to `processSpatialSample()` as
`igblast_refs`. If you leave `igblast_refs` out, SIRUP looks for those
folders in `igblast_base`. Give the IMGT folder as `imgt_dir`. airrblast
supports human, mouse, and rhesus monkey.

## Demultiplex with rrad

rrad is a separate package. Install it from GitHub:

``` r
if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
if (!requireNamespace("rrad", quietly = TRUE)) remotes::install_github("indianewok/rrad")
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
`processSpatialSample()` calls it for you. Use it yourself to check a
header. This chunk runs:

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

SIRUP includes 2,000 read-level IgBLAST rows from one section of human
thymus. The file stores the barcode and UMI in the columns
`cell_barcode` and `umi`, not in the header. The sample is enriched for
usable reads. In this sample, 750 rows pass the IGH filter below. In the
full section, 0.4% of the reads pass. The chunks in this section and the
next run on the example.

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

> **Note.** Keep `productive_only = FALSE` at the read level. A
> sequencing error can make a single read look unproductive, and the
> consensus step corrects that error. Filter on productivity after the
> consensus step. With `umi_support = TRUE`, SIRUP builds each molecule
> consensus from the productive reads of that molecule.

`plotLengthQc()` shows the V, J, and junction lengths. Use it to choose
the thresholds:

``` r
plotLengthQc(igh_reads, bin_width = 10, title = "IGH reads, thymus section")
```

![](processing_files/figure-gfm/length-qc-1.jpeg)

## Run the consensus pipeline

`processSpatialSample()` runs every step from the table above. It prints
each command before it runs it. If a command fails, it stops and shows
the command. The chunks in this section show code only. They need the
external tools and full-size inputs.

You can start from three types of input:

| `input` | `input_type` | First step that runs |
|----|----|----|
| FASTA or FASTQ from rrad | `"auto"` or `"raw_reads"` | IgBLAST on each read |
| AIRR table (path or data frame) from a first IgBLAST pass | `"auto"` or `"airr"` | `filterAirr()` |
| FASTA that is already filtered, with the spot in the header | `"spatial_fasta"` | pRESTO clustering |

### From rrad FASTA

``` r
result <- processSpatialSample(
  input = "prepared/sample-igblast.fasta",
  output_dir = "processed/bcr",
  receptor = "ig",
  organism = "human",
  igblast_base = "~/share/igblast",
  igblast_refs = "~/share/igblast_refs",
  imgt_dir = "~/share/imgt",
  muscle_exec = "~/bin/muscle3.8.31_i86linux64",
  identity = 1.0,
  consensus_freq = 0.51,
  consensus_min_count = 2,
  nproc = 8
)

result$db_pass
```

Use `receptor = "tr"` for TCRs. Run BCRs and TCRs as separate calls.

`runIgblast()` calls `airrblast::assignGenes()`, which does not search a
constant-region database. A table made from raw reads has no `c_call`
column.

### With UMI support

Reads of one molecule can differ by sequencing errors. At
`identity = 1.0`, those reads fall into different clusters, and each
cluster can then fail the two-read minimum. With `umi_support = TRUE`,
SIRUP first groups the reads of each spot and UMI at a looser identity
(`umi_identity`, default 0.95). It builds one consensus per molecule
from the productive reads only. The spot-level step then combines the
molecules and sums their reads.

`umi_support` needs an AIRR table with `umi` and `productive` columns.
Make it with a first IgBLAST pass on the rrad FASTA file:

``` r
runIgblast(
  input = "prepared/sample-igblast.fasta",
  output_tsv = "processed/reads_igblast.tsv",
  receptor = "ig",
  igblast_base = "~/share/igblast",
  igblast_refs = "~/share/igblast_refs"
)
igblast_reads <- airr::read_rearrangement("processed/reads_igblast.tsv")
igblast_reads <- parseRadSequenceIds(igblast_reads)
```

The rrad FASTA header holds the barcode and UMI. `parseRadSequenceIds()`
copies them into the `barcode` and `umi` columns. The SIRUP paper used
this call for the thymus BCRs:

``` r
result <- processSpatialSample(
  input = igblast_reads,           # read-level AIRR data frame
  input_type = "airr",
  output_dir = "processed/bcr_umi",
  receptor = "ig",
  identity = 1.0,
  consensus_freq = 0.51,
  cluster_field = "barcode",
  igblast_base = "~/share/igblast",
  igblast_refs = "~/share/igblast_refs",
  imgt_dir = "~/share/imgt",
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
`c_call` into the FASTA header. Then give
`consensus_copy_fields = "C_CALL"` and
`consensus_copy_action = "majority"`. Each consensus sequence gets the
constant region of most of its reads, in `c_call`, and the fraction of
reads that agree, in `c_call_freq`.

### Many samples

`processSpatialSamples()` runs `processSpatialSample()` on each input,
with the same settings:

``` r
processSpatialSamples(
  inputs = c("prepared/slice1.fasta", "prepared/slice2.fasta"),
  output_dirs = c("processed/slice1", "processed/slice2"),
  receptor = "ig",
  igblast_base = "~/share/igblast",
  igblast_refs = "~/share/igblast_refs",
  imgt_dir = "~/share/imgt",
  nproc = 8
)
```

> **Tip.** On a cluster, submit one job per sample. Each job calls
> `processSpatialSample()` once.

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
in a spot is more likely a sequencing error. For TCRs, use
`identity = 0.99`. In the SIRUP paper, `consensus_freq = 0.51` gave the
most high-quality sequences for both chains.

To compare settings, run `processSpatialSample()` with each value.
Clustering and alignment depend on `identity` only, so with
`skip_existing = TRUE` you can try several `consensus_freq` values and
reuse the clustering.

## Filter the consensus table

The pipeline writes the consensus `sequence_id` as
`<cluster>,<barcode>`. The example data include the consensus table of
one thymus section. It is the `db_pass_umi` table that
`processSpatialSample()` wrote for the SIRUP paper. The chunks in this
section run on it:

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
#>  Mean   :0.9982        Mean   : 63.71   Mean   : 1.442  
#>  3rd Qu.:1.0000        3rd Qu.:131.00   3rd Qu.: 1.000  
#>  Max.   :1.0000        Max.   :220.00   Max.   :35.000
```

Then filter. These settings follow the definition of a quality sequence
in the SIRUP paper. A quality sequence is a productive rearrangement.
Its V alignment has at least 250 nucleotides, its J at least 30, and its
junction at least 10. It has no gap longer than 15 nucleotides.

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
#> filterAirr kept 2073 of 3793 rows (dropped 1720).

nrow(consensus)
#> [1] 3793
nrow(bcrs)
#> [1] 2073
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

This chunk shows code only. It needs Space Ranger output.

``` r
# bcrs$barcode_2um holds the 2 µm barcode of each consensus sequence.
mapping <- readSpatialIdMapping("outs/barcode_mappings.parquet")
mean(bcrs$barcode_2um %in% mapping$square_002um)

bcrs <- bcrs %>%
  left_join(
    select(mapping, square_002um, square_008um, square_016um),
    by = c("barcode_2um" = "square_002um")
  ) %>%
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

> **Note.** Use it on Visium HD data only. Standard Visium spots are 55
> µm wide, 100 µm apart, and hold several cells each. Two neighboring
> Visium spots with the same sequence are more likely two related cells
> than spatial bleed.

The function estimates the neighbor distance from the spot coordinates.
Set `num_neighbors_non_boundary` to the number of neighbors of a spot in
the middle of the tissue. Use 6 on the standard Visium hexagonal grid
and 8 on the Visium HD square grid. The next chunk runs the function on
the standard Visium example to show the call and its output. Do not use
it on standard Visium data in practice:

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
#> ✔ Ensuring one row per spatial_id/sequence... [2s]
#> 
#> ℹ Estimating distance threshold from full tissue positions...
#> ✔ Estimating distance threshold from full tissue positions... [3.4s]
#> 
#> ℹ Distance threshold: 10.0214
#> ℹ Filtering to spots with observed sequences...
#> ✔ Filtering to spots with observed sequences... [7ms]
#> 
#> ℹ Sequence-positive spots: 800 (of 3,203 total)
#> ℹ Building adjacency matrix over spots with observed sequences...
#> ✔ Building adjacency matrix over spots with observed sequences... [10ms]
#> 
#> ℹ Collapsing components to representative spots...
#> ✔ Collapsing components to representative spots... [220ms]
#> 
#> ℹ Filtering data to representative spots...
#> ✔ Filtering data to representative spots... [11ms]
#> 
#> ✔ Done! Retained 1,946 rows from 2,071 (125 bleedover sequences removed)

c(before = nrow(getReceptors(sirup)), after = nrow(getReceptors(collapsed)))
#> before  after 
#>   2073   1946
```

For Visium HD, give the 8 µm object and
`num_neighbors_non_boundary = 8`.

## Next steps

- Build the object from the table in
  [Get started with SIRUP](get-started.md).
- Draw tissue regions and measure clonal dominance in
  [Draw tissue regions and measure clonal dominance and diversity](regions.md).
- Call clones and build lineage trees in
  [Call clones and build lineage trees on the tissue](clones-and-trees.md).
- See what the object stores in
  [How the SIRUP object works](sirup-object.md).
