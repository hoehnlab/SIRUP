# SIRUP: Spatial Immune Receptor Utilities Package

SIRUP is an R package for spatial adaptive immune receptor repertoire
sequencing (spAIRR-seq). It works with B cell receptor (BCR) and T cell
receptor (TCR) sequences that carry a spatial barcode, together with spatial
gene expression from the same tissue section.

SIRUP covers three stages:

1. **Process.** Turn demultiplexed long reads into consensus receptor
   sequences. IgBLAST and pRESTO build an AIRR rearrangement table with one
   spatial ID per sequence.
2. **Integrate.** Put the receptor table and the spatial Seurat object into one
   validated SIRUP object. One object can hold several tissue sections, bulk
   AIRR-seq, and blood.
3. **Analyze and visualize.** Draw tissue regions, measure clonal dominance and
   diversity, test spatial patterns, and show lineage trees on the tissue.

SIRUP is in active development. Function names and arguments can change before
the first release.

## Citation

Melton, H.J., Vaidya, C.M., Panda, A., Fielding, J.J., Sergison, E., Kolling, F.W., O'Connor, K.C. Hoehn, K.B. *SIRUP, a toolkit for spatial B and T cell receptor repertoires, reveals differences between normal and ectopic germinal centers in myasthenia gravis.*  bioRxiv, date, doi.

## Getting started

Start with the [get-started guide](TODO). It shows how to build
the SIRUP object and a few plots and analyses.

A full tutorial is available at [TODO].


## Install SIRUP and dependencies

Use R 4.1 or later. For now, install it from GitHub with BiocManager. This also
installs the Bioconductor packages that SIRUP needs, such as ggtree:

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("hjm2217/SIRUP")
```

TODO: CRAN

### Other dependencies

SIRUP needs these programs only to process reads with `processSpatialSample()`
and `runIgblast()`. You do not need them to analyze an AIRR table that you
already have.

| Program | Use | Download | Install |
|---|---|---|---|
| [pRESTO](https://presto.readthedocs.io) | Clusters the reads of each spot, aligns them, and builds one consensus sequence for each group (`ClusterSets.py`, `ParseHeaders.py`, `AlignSets.py`, `BuildConsensus.py`) | [PyPI](https://pypi.org/project/presto/) or [GitHub](https://github.com/immcantation/presto) | `pip3 install presto --user`. It needs Python 3.7 or later. |
| [vsearch](https://github.com/torognes/vsearch) | Clusters reads by sequence identity. Used with `ClusterSets.py` during `processSpatialSample()`. | [GitHub releases](https://github.com/torognes/vsearch/releases) | Download the binary for Linux, macOS, or Windows. Add its `bin` folder to your `PATH`. You can also use `conda install -c bioconda vsearch` or `brew install vsearch`. |
| [MUSCLE](https://drive5.com/muscle/) v3.8.31 | Aligns the reads in each group. `AlignSets.py` calls it. | [drive5.com, MUSCLE v3 downloads](https://drive5.com/muscle/downloads_v3.htm) | Download the binary for your system. Unpack it with `tar -xvzf`. Give the path of the file to `muscle_exec`. |
| [IgBLAST](https://ncbi.github.io/igblast/) | Assigns V(D)J genes to reads and to consensus sequences. Used during `processSpatialSample()` through the airrblast package. | [NCBI FTP](https://ftp.ncbi.nih.gov/blast/executables/igblast/release/) | Follow the instructions in the [airrblast repo](https://github.com/hoehnlab/airrblast/tree/main) to set IgBLAST up. |

Notes:

- **MUSCLE version.** SIRUP uses muscle through [pRESTO](https://presto.readthedocs.io), which requires MUSCLE v3. MUSCLE v5 is not yet supported to our knowledge.


## License

MIT. See [LICENSE](LICENSE).
