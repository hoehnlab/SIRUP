# downloadSIRUPExample

## Download the SIRUP example data

## Description

Download the myasthenia gravis (MG) thymus example files that the guides on the SIRUP website use. The files are not part of the package, to keep the package small. The function skips files that are already in `dest` and checks every download against a stored MD5 sum.

## Usage

``` R
downloadSIRUPExample(
  dest = file.path(tempdir(), "SIRUP_example"),
  files = NULL,
  overwrite = FALSE,
  base_url = getOption("SIRUP.example_url", .sirup_example_url),
  quiet = FALSE
)
```

## Arguments

`dest`

:   Folder to save the files in. The function creates it if needed. The default is a folder in the R session\'s temporary directory, so the files are lost when you close R. Give a permanent folder to keep them.

`files`

:   Names of the files to download. If `NULL` (default), the function downloads all files. Use this to download only what a guide needs.

`overwrite`

:   Download a file again even if it is already in `dest` (default `FALSE`)?

`base_url`

:   Web address of the folder that holds the files. Change it only to use a mirror. A `⁠file://⁠` address of a local folder also works. The default is the option `SIRUP.example_url` if you set it.

`quiet`

:   Hide the download progress messages (default `FALSE`)?

## Details

The files are:

- `mg_thymus_slice1_spatial.rds`, `mg_thymus_slice2_spatial.rds`: spatial Seurat objects for the two tissue sections.

- `mg_thymus_slice1_bcr.tsv.gz`, `mg_thymus_slice2_bcr.tsv.gz`: AIRR rearrangement tables of IGH sequences with spatial IDs.

- `mg_thymus_slice2_tcr.tsv.gz`, `mg_thymus_slice2_reads.tsv.gz`: TCR sequences and the demultiplexed reads of section 2.

- `mg_thymus_slice1_gc_polygons.csv`, `mg_thymus_slice2_gc_polygons.csv`: annotated polygons of the ectopic germinal centers.

- `mg_thymus_slice2_qupath_gcs.geojson`: six germinal center outlines that a user drew on section 2 in QuPath. Use it with [`addQuPathPolygons()`](addQuPathPolygons.md).

- `mg_thymus_cell_types.csv.gz`: cell-type proportions for each spot.

- `mg_thymus_slice2.sirup.rds`: a saved SIRUP object for section 2, with BCR and TCR slots, cell-type proportions, and BCR clones. Read it with [`readSIRUP()`](sirup-io.md).

- `mg_thymus_slice1_spatial_brighter.rds`, `mg_thymus_slice2_spatial_brighter.rds`, `mg_thymus_slice2_brighter.sirup.rds`: the same objects with a brighter H&E image. The guides use these files.

- `mg_thymus_bulk_*.tsv.gz`: bulk repertoire tables from thymus and blood.

## Value

The path of `dest`, invisibly. Use `file.path()` with the path and a file name to read a file.

## See Also

[`readSIRUP()`](sirup-io.md)

## Examples

``` R

extdata <- downloadSIRUPExample(
  files = c("mg_thymus_slice2_spatial.rds", "mg_thymus_slice2_bcr.tsv.gz")
)
thymus <- readRDS(file.path(extdata, "mg_thymus_slice2_spatial.rds"))
```
