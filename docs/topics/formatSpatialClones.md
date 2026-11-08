# formatSpatialClones

## Format the receptors of a SIRUP object as in dowser\'s `formatClones()`

## Description

`formatSpatialClones` runs [`dowser::formatClones()`](https://dowser.readthedocs.io/en/latest/topics/formatClones/) on the receptor table of a SIRUP object. It stores the result as an analysis.

## Usage

``` R
formatSpatialClones(
  object,
  analysis_name = "f_clones",
  seq = "sequence_alignment",
  clone = "clone_id",
  subgroup = "clone_subgroup",
  id = "sequence_id",
  germ = "germline_alignment_d_mask",
  v_call = "v_call",
  j_call = "j_call",
  junc_len = "junction_length",
  mask_char = "N",
  max_mask = 0,
  pad_end = TRUE,
  text_fields = NULL,
  num_fields = NULL,
  seq_fields = NULL,
  add_count = TRUE,
  verbose = FALSE,
  collapse = TRUE,
  cell = "cell_id",
  locus = "locus",
  traits = NULL,
  mod3 = TRUE,
  randomize = TRUE,
  use_regions = TRUE,
  dup_singles = FALSE,
  nproc = 1,
  chain = "H",
  heavy = "IGH",
  filterstop = TRUE,
  minseq = 2,
  split_light = FALSE,
  light_traits = FALSE,
  majoronly = FALSE,
  columns = NULL,
  validate = TRUE
)
```

## Arguments

`object`

:   A SIRUP object.

`analysis_name`

:   Name for the stored formatted clones (default `"f_clones"`).

`seq`

:   Name of the column that holds the observed DNA sequences (default `"sequence_alignment"`). The sequences must be aligned to each other.

`clone`

:   Name of the column that holds the clone ID (default `"clone_id"`).

`subgroup`

:   Name of the column that holds the subgroup ID (default `"clone_subgroup"`).

`id`

:   Name of the column that holds the sequence IDs (default `"sequence_id"`).

`germ`

:   Name of the column that holds the germline DNA sequences (default `"germline_alignment_d_mask"`). Each clone must have one germline. It must be aligned to the sequences in `seq`.

`v_call`

:   Name of the column that holds the V-segment alleles (default `"v_call"`). Within a clone, the entries must be the same at the gene level.

`j_call`

:   Name of the column that holds the J-segment alleles (default `"j_call"`). Within a clone, the entries must be the same at the gene level.

`junc_len`

:   Name of the numeric column that holds the junction length (default `"junction_length"`). Within a clone, the entries must be the same.

`mask_char`

:   Character for masking and padding (default `"N"`).

`max_mask`

:   Largest number of characters to mask at the start and at the end of the sequences (default 0). With `0`, the function does not mask. With `NULL`, the function uses the largest number of leading or trailing Ns in the sequences.

`pad_end`

:   Pad the end of each sequence with `mask_char` so that all sequences have the same length (default `TRUE`)?

`text_fields`

:   Names of the text columns to keep and merge when the function removes duplicates.

`num_fields`

:   Names of the numeric columns to keep and sum when the function removes duplicates.

`seq_fields`

:   Names of the sequence columns to keep and collapse when the function removes duplicates. Do not repeat `seq` or `germ` here.

`add_count`

:   Add the column `COLLAPSE_COUNT`, which holds the number of sequences that the function collapsed (default `TRUE`)?

`verbose`

:   Report the numbers of input, dropped, and output sequences (default `FALSE`)? The function passes it to `collapseDuplicates`.

`collapse`

:   Collapse identical sequences (default `TRUE`)?

`cell`

:   Name of the column that holds the cell IDs (default `"cell_id"`).

`locus`

:   Name of the column that holds the locus (default `"locus"`).

`traits`

:   Names of the columns to keep separate when the function collapses sequences. Include the spatial ID column if you want to plot the tips on the tissue.

`mod3`

:   Pad the sequences to a length that is a multiple of three (default `TRUE`)?

`randomize`

:   Put the sequences in random order (default `TRUE`)? This is important for PHYLIP.

`use_regions`

:   Assign the CDR and FWR regions (default `TRUE`)?

`dup_singles`

:   Duplicate the sequence of a clone that has one sequence, so that the clone gives a tree (default `FALSE`)?

`nproc`

:   Number of cores for the formatting (default 1).

`chain`

:   Use `"HL"` to include the light chain, if it exists. The default is `"H"`.

`heavy`

:   Name of the heavy chain locus (default `"IGH"`).

`filterstop`

:   Use only sequences with no in-frame stop codon (default `TRUE`)?

`minseq`

:   Smallest number of sequences in a clone (default 2).

`split_light`

:   Split the subgroups, or lump them together (default `FALSE`)? See `resolveLightChains` in dowser.

`light_traits`

:   Include the traits of the light chain when the function joins and collapses the sequences (default `FALSE`)?

`majoronly`

:   Return only the largest subgroup and the sequences with no light chain (default `FALSE`)?

`columns`

:   Names of other columns to keep in the result.

`validate`

:   Validate the updated SIRUP object (default `TRUE`)?

## Details

The function puts clones into a standardized format for downstream analysis. See [`dowser::formatClones()`](https://dowser.readthedocs.io/en/latest/topics/formatClones/) for details. The function stores the result as an analysis in the SIRUP object. The default analysis name is `"f_clones"`.

[`dowser::formatClones()`](https://dowser.readthedocs.io/en/latest/topics/formatClones/) takes a data frame and returns a tibble of `airrClone` objects. `formatSpatialClones` calls it with the receptor table of a SIRUP object and adds four things that spatial analysis needs: This function automatically does not collapse identical sequences across distinct samples. To plot the tips on the tissue, add the spatial ID column to `traits`.

## Value

The updated SIRUP object, with the formatted clones as an analysis.

## See Also

[`getSpatialTrees()`](getSpatialTrees.md), [`dowser::formatClones()`](https://dowser.readthedocs.io/en/latest/topics/formatClones/)

## Examples

``` R
x <- readSIRUP(system.file("extdata", "sirup_synthetic.rds", package = "SIRUP"))
x <- formatSpatialClones(
  x,
  id = "sequence_id",
  germ = "germline_alignment_d_mask",
  traits = "spatial_id"
)
getAnalysis(x, "f_clones")
```
