# SIRUP: Spatial Immune Receptor Utilities Package

> **Under construction.** This site is not finished. Some pages can be
> incomplete or wrong. SIRUP is in active development, and function names and
> arguments can change before the first release.

SIRUP is an R package for spatial adaptive immune receptor repertoire
sequencing (spAIRR-seq). It works with B cell receptor (BCR) and T cell receptor
(TCR) sequences that carry a spatial barcode, together with spatial gene
expression from the same tissue section.

SIRUP covers three stages:

1. **Process.** Turn demultiplexed long reads into consensus receptor
   sequences. IgBLAST and pRESTO build an AIRR rearrangement table with one
   spatial ID for each sequence.
2. **Integrate.** Put the receptor table and the spatial Seurat object into one
   validated SIRUP object. One object can hold several tissue sections, bulk
   AIRR-seq, and blood.
3. **Analyze and visualize.** Draw tissue regions, measure clonal dominance and
   diversity, test spatial patterns, and show lineage trees on the tissue.

## Getting started

Start with [Get started with SIRUP](vignettes/get-started.md). It takes about 15
minutes and uses example data that comes with the package. You build a SIRUP
object, plot receptors on the tissue, and measure clonal dominance in a region.

Install SIRUP from GitHub with BiocManager. BiocManager also installs the
Bioconductor packages that SIRUP needs, such as ggtree:

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("hoehnlab/SIRUP")
```

Use R 4.1 or later.

## Contents

**Tutorials**

* [Process reads into spatial AIRR tables](vignettes/processing.md). Turn
  demultiplexed reads into a table of consensus BCR or TCR sequences.
* [Draw tissue regions and measure clonal dominance and
  diversity](vignettes/regions.md). Compare the repertoires of regions such as
  germinal centers.
* [Call clones and build lineage trees on the
  tissue](vignettes/clones-and-trees.md). Build trees with Dowser and show them
  on the H&E image.
* [Test spatial patterns of genes, cells, and
  clones](vignettes/spatial-statistics.md). Run Moran's I, Lee's L, Ripley's
  cross K, and join counts.
* [Analyze several sections with bulk AIRR-seq and
  blood](vignettes/multisample.md). Combine sections, bulk data, and blood in
  one object.

**Guides**

* [How the SIRUP object works](vignettes/sirup-object.md). The five parts of the
  object, and the rules that SIRUP checks.
* [Choose a spatial statistic](vignettes/choose-spatial-statistic.md). Match
  your question to a function.

**Function reference**

The Function reference section lists every function, with its arguments and
examples.

## Contact

If you need help or have a question, contact the Hoehn Lab. To report a bug or
ask for a feature, open an issue on the
[issue tracker](https://github.com/hoehnlab/SIRUP/issues).

## Authors

Hunter J. Melton and Chinmay Vaidya.

## Citing SIRUP

Melton HJ, Vaidya CM, Panda A, Fielding JJ, Sergison E, Kolling FW, O'Connor KC,
Hoehn KB. SIRUP, a toolkit for spatial B and T cell receptor repertoires, reveals
differences between normal and ectopic germinal centers in myasthenia gravis.
bioRxiv.

## License

MIT. See the
[LICENSE file](https://github.com/hoehnlab/SIRUP/blob/main/LICENSE).
