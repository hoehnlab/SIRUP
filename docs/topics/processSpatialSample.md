# processSpatialSample

## Process one spatial sample from reads to an AIRR table

## Description

`processSpatialSample`processes demultiplexed reads for one sample and one receptor group (`ig` or `tr`). It outputs an AIRR table of filtered consensus sequences, with the V(D)J genes assigned and IMGT-gaps added.

## Usage

``` R
processSpatialSample(
  input,
  output_dir = NULL,
  receptor = c("ig", "tr"),
  organism = "human",
  input_type = c("auto", "raw_reads", "airr", "spatial_fasta"),
  filter_args = list(),
  invalid_spatial_ids = c("abort", "drop"),
  fasta_sequence_col = "sequence_alignment",
  sample_name = NULL,
  identity = 0.975,
  consensus_freq = 0.51,
  germline_glob = NULL,
  igblast_refs = NULL,
  imgt_dir = NULL,
  igblast_base = "~/software/igblast",
  muscle_exec = "~/software/muscle3.8.31_i86linux64",
  nproc = 8,
  cluster_field = "barcode",
  delimiters = c("|", "=", ","),
  parse_fields = c(cluster_field, "cluster"),
  parse_key = paste0(cluster_field, "_cluster"),
  parse_act = "set",
  consensus_min_count = 1,
  umi_support = FALSE,
  umi_weight = TRUE,
  umi_field = "umi",
  umi_identity = 0.95,
  min_umi_support = 2,
  consensus_copy_fields = NULL,
  consensus_copy_action = "majority",
  consensus_max_error = 0.1,
  consensus_max_gap = 0.3,
  skip_existing = FALSE,
  align_log = NULL,
  consensus_log = NULL
)
```

## Arguments

`input`

:   Raw reads (FASTA or FASTQ, optionally gzipped), the path of an AIRR TSV, an AIRR data frame, or a FASTA that is ready for pRESTO.

`output_dir`

:   Directory for all outputs. The default is the directory of the input file. You must give it if `input` is a data frame.

`receptor`

:   `"ig"` for BCR loci or `"tr"` for TCR loci.

`organism`

:   Organism of the sample (default `"human"`).

`input_type`

:   How to treat `input`: `"auto"` (default), `"raw_reads"`, `"airr"`, or `"spatial_fasta"`. With `"auto"`, the function finds the type from the class of the object and the file extension. A FASTA or FASTQ path becomes `"raw_reads"`. Use `"spatial_fasta"` to send a FASTA that is ready for pRESTO straight to the consensus workflow. For a new rrad FASTA, use `"auto"` or `"raw_reads"`, because `"spatial_fasta"` skips the first annotation of V(D)J genes and the filtering of the reads.

`filter_args`

:   Named list that changes arguments of [`filterAirr()`](filterAirr.md). The pipeline defaults are `complete_vdj = TRUE`, `min_consensus_count = NULL`, `spatial_id = "barcode"`, `add_spatial_sequence_id = TRUE`, and `spatial_id_key = cluster_field`. A `NULL` value turns off an optional filter (`locus_val`, `min_consensus_count`, `max_opening_gap`, or `spatial_id`). The other arguments do not accept `NULL`.

`invalid_spatial_ids`

:   What to do with filtered rows whose spatial barcode is missing or holds a delimiter: `"abort"` (default) stops with a count, and `"drop"` removes the rows with a warning.

`fasta_sequence_col`

:   Name of the AIRR column to write to the spatial FASTA (default `"sequence_alignment"`, the aligned V(D)J region). This keeps extra read sequence out of the clustering and consensus steps. Use `"sequence"` to write the full read.

`sample_name`

:   Stem of the intermediate file names. The default is the input file name without the extension, or `"sample"` for a data frame.

`identity`

:   Identity threshold for the clusters, passed to `ClusterSets.py` (default 0.975).

`consensus_freq`

:   Frequency threshold for the consensus, passed to `BuildConsensus.py` (default 0.51).

`germline_glob`

:   Pattern of the germline files, such as `⁠<imgt_dir>/human/vdj/imgt_human_*.fasta⁠`. The function uses it only to find `imgt_dir` if `imgt_dir` is `NULL`.

`igblast_refs`

:   Folder with the IgBLAST `⁠database/⁠` and `⁠fasta/⁠` folders (the `refs` argument of `airrblast::assignGenes()`). If `NULL`, the function uses `igblast_base`.

`imgt_dir`

:   IMGT folder with the gapped germlines, at `⁠<imgt_dir>/<organism>/vdj⁠` (the `gapdb` argument of `airrblast::addGaps()`). If `NULL`, the function finds it from `germline_glob`, which must point into a `⁠<organism>/vdj⁠` folder.

`igblast_base`

:   IgBLAST installation folder (default `"~/software/igblast"`). The function passes it to `airrblast::assignGenes()` as `igblast`. The folder must hold `bin/igblastn`.

`muscle_exec`

:   Path of the MUSCLE program, passed to `AlignSets.py --exec` (default `"~/software/muscle3.8.31_i86linux64"`).

`nproc`

:   Number of processors for the IgBLAST step on raw reads, for pRESTO, and for the gene assignment (default 8).

`cluster_field`

:   FASTA header field for the clustering within a spot (default `"barcode"`, the parsed rrad `BARCODE` field). The pipeline writes the headers with this field name. Change it for other spatial ID fields, such as `barcode_8um`.

`delimiters`

:   Character vector of the pRESTO delimiters, passed after `--delim` (default `c("|", "=", ",")`).

`parse_fields`

:   Header fields for `⁠ParseHeaders.py merge -f⁠` (default `c(cluster_field, "cluster")`).

`parse_key`

:   Name of the header field that `⁠ParseHeaders.py merge -k⁠` makes (default `paste0(cluster_field, "_cluster")`).

`parse_act`

:   Action for `⁠ParseHeaders.py merge --act⁠` (default `"set"`). With `"set"`, the function joins the spatial ID and the cluster ID with a comma, so that you can recover the spatial ID from the consensus header. With `"cat"`, it joins them with no separator, so that you cannot recover the spatial ID of a barcode that has a variable length. Both give the same groups, so use `"set"`.

`consensus_min_count`

:   Smallest number of sequences for a consensus, passed to `BuildConsensus.py -n` (default 1). Use at least 2 for data that do not have single-cell resolution. This value decides which groups the function writes. It does not change the consensus of a group that passes. So you can run with `1` and filter to a higher value later with `filterAirr(min_consensus_count = )`. `BuildConsensus.py -q` is fixed at 20 and has no effect here, because the pipeline input is FASTA and has no base qualities.

`umi_support`

:   Group the reads by spatial barcode and UMI before the barcode step (default `FALSE`)? Reads of one molecule that differ by a sequencing error end up in different clusters at `identity = 1.0`, and the function then drops them as singletons. This stage first joins them at `umi_identity`, so that they support each other. Each barcode-UMI cluster gives one consensus. The consensus carries `cc_bc_umi` (the number of reads in the cluster, counted before the productive subset) and `n_productive_reads` (the number of reads that the consensus uses). The function discards a cluster that has no such reads. The barcode step then runs with `⁠BuildConsensus -n 1⁠` and sums both fields, and `min_umi_support` replaces the read minimum. With `umi_weight = TRUE`, the function writes each consensus once for each productive read, so that one read is one vote in the barcode consensus. A short read from one molecule then cannot trim the consensus of many full-length reads. For this reason, `consensus_count` counts reads and not molecules. The option needs an AIRR table with the columns `umi` and `productive`.

`umi_weight`

:   Write each molecule consensus once for each productive read (default `TRUE`)? The function uses it only with `umi_support`. If `FALSE`, it writes each molecule once, so that one UMI cluster is one vote.

`umi_field`

:   Name of the column that holds the UMI (default `"umi"`). The function uses it only with `umi_support`.

`umi_identity`

:   Identity threshold for the barcode-UMI clusters (default 0.95). It is lower than `identity` on purpose, so that reads of one molecule with errors form one group.

`min_umi_support`

:   Smallest summed `cc_bc_umi` that a barcode consensus needs (default 2). The function applies it after the gene assignment. It writes a filtered copy next to `db_pass`, so that you can change the threshold without a new run of the pipeline.

`consensus_copy_fields`

:   Header fields for `BuildConsensus.py --cf`. With `umi_support`, the function adds `CC_BC_UMI` and `N_PRODUCTIVE_READS` at the start of the list and sums them.

`consensus_copy_action`

:   Action for `BuildConsensus.py --act`, if you give `consensus_copy_fields` (default `"majority"`).

`consensus_max_error`

:   Largest error rate, passed to `BuildConsensus.py --maxerror` (default 0.1).

`consensus_max_gap`

:   Largest gap fraction, passed to `BuildConsensus.py --maxgap` (default 0.3).

`skip_existing`

:   Skip a pipeline command if its output already exists and is not empty (default `FALSE`)? Use it to restart an interrupted run. You can also use it to reuse one clustering and alignment pass for several `consensus_freq` values, because those steps depend only on `identity`.

`align_log`

:   Path of a log file for `AlignSets.py --log`.

`consensus_log`

:   Path of a log file for `BuildConsensus.py --log`.

## Details

The function calls external programs: IgBLAST, pRESTO, vsearch, and MUSCLE. It prints each command before it runs it. The function runs all commands with `system2()` and no shell. A non-zero exit status stops the function, and the error message shows the command that failed. The final output is the `db_pass` AIRR TSV, which is in the returned list of paths.

The `input` can be one of these:

- Raw reads: a FASTA or FASTQ file, optionally gzipped, from RRAD (or RAD). Export the rrad FASTA with `to_fasta = TRUE` and `header_format = "presto"`, so that the barcode and UMI tags stay in the IgBLAST identifier. The function first annotates the reads with [`runIgblast()`](runIgblast.md) (which calls `airrblast::assignGenes()`). This step needs `igblast_base` and `igblast_refs`.

- An AIRR table: a `.tsv` path or a data frame. The function skips the IgBLAST step. A data frame needs `output_dir`.

- A FASTA that is ready for pRESTO (`input_type = "spatial_fasta"`). The function skips the IgBLAST and filter steps, and starts with the consensus workflow.

For raw reads and AIRR tables, the function first does these steps:

1.  It splits the rrad tags `⁠|BARCODE=...|UMI=...⁠` (or the legacy colon tags) out of `sequence_id` with [`parseRadSequenceIds()`](parseRadSequenceIds.md).

2.  It filters sequences with [`filterAirr()`](filterAirr.md). The pipeline defaults are `complete_vdj = TRUE`, `spatial_id = "barcode"`, and a combined header ID with the key `cluster_field`. Change them with `filter_args`.

3.  It writes the raw and the filtered AIRR tables to `output_dir`, so that you can inspect them.

4.  It writes the filtered sequences to a spatial FASTA file. The headers of the file hold the spatial ID.

The spatial FASTA then goes through a number of Immcantation functions:

1.  `⁠ClusterSets.py set⁠` groups the reads within each spatial ID. It uses `cluster_field`, the identity threshold `identity`, and vsearch.

2.  `⁠ParseHeaders.py merge⁠` writes a new header field, `parse_key`. The field combines the spatial ID and the cluster ID (`parse_fields`) with the action `parse_act`. `delimiters` sets the delimiters of the FASTA header.

3.  `⁠AlignSets.py muscle⁠` aligns the reads within each `parse_key` group. `muscle_exec` is the path of the MUSCLE program (<https://github.com/rcedgar/muscle>). The function can write a log to `align_log`.

4.  `BuildConsensus.py` makes a consensus sequence for each `parse_key` group. `consensus_min_count` is the smallest number of sequences that give a valid consensus. `consensus_freq` is the fraction of occurrences below which the function assigns an ambiguous character. The function can write a log to `consensus_log`. To copy more header fields to the consensus file, give their names in `consensus_copy_fields` and the way to copy them in `consensus_copy_action` (`majority`, `first`, `last`, or `all`).

5.  `airrblast::assignGenes()` assigns the V(D)J genes with IgBLAST. It replaces `AssignGenes.py`. `receptor` is `ig` (BCR) or `tr` (TCR). `organism` is the species (human, mouse, or rhesus_monkey). `igblast_base` is the IgBLAST installation, which must hold `bin/igblastn`. `igblast_refs` is the folder with the IgBLAST `⁠database/⁠` and `⁠fasta/⁠` folders.

6.  `airrblast::addGaps()` adds the IMGT-gapped sequences from the germlines in `imgt_dir`. SIRUP then reads the pRESTO header fields into columns, as `MakeDb.py --extended` did. It keeps the sequences that have a V call, a J call, a productive flag, an aligned sequence, and a junction. The result is the `db_pass` table.

## Value

A named list, invisibly. It holds `raw_airr`, `filtered_airr`, `spatial_fasta`, `n_raw`, and `n_filtered`, which are `NULL` for the steps that did not run. It also holds the output paths of pRESTO and IgBLAST. `db_pass` is the final AIRR TSV. With `umi_support`, `db_pass_umi` is the same table, filtered to `min_umi_support`.

## See Also

[`processSpatialSamples()`](processSpatialSamples.md), [`runIgblast()`](runIgblast.md), [`filterAirr()`](filterAirr.md), [`parseRadSequenceIds()`](parseRadSequenceIds.md)

## Examples

``` R
## Not run: 
# This example needs IgBLAST, pRESTO, vsearch, and MUSCLE. See
# vignette("processing", package = "SIRUP").
result <- processSpatialSample(
  "sample1_reads.fasta.gz",
  output_dir = "sample1_output",
  receptor = "ig",
  organism = "human",
  igblast_base = "~/software/igblast",
  imgt_dir = "~/share/germlines/imgt",
  muscle_exec = "~/software/muscle",
  nproc = 8
)
bcrs <- airr::read_rearrangement(result$db_pass)

## End(Not run)
```
