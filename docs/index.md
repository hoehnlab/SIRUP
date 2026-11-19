# SIRUP: Spatial Immune Receptor Utilities Package

<div markdown="1" style="background: #9e584a; padding: 12px; margin-bottom: 24px; line-height: 24px;">
  <p style="margin: -12px -12px 12px; padding: 6px 12px; background: #8b3620; color: #ffffff; font-weight: 700;">Under construction</p>
  This site is still under construction. Be aware that some information may be missing or incorrect.
</div>

<img src="img/sirup-logo.png" alt="SIRUP logo" align="right" width="130" style="margin: 0 0 12px 24px;">
SIRUP, the Spatial Immune Receptor Utilities Package, is an R package for working with
spatial adaptive immune receptor repertoire sequencing (spAIRR-seq) data. 
It works with B cell receptor (BCR) and T cell receptor
(TCR) sequences that have a spatial barcode, together with spatial gene
expression from the same tissue section. Isn't that sweet?

SIRUP helps you:

1. **Process.** Convert demultiplexed long-read spatial immune receptor sequences 
   into high-quality consensus receptor sequences for each spatial location.
   Outputs an [AIRR rearrangement style table](https://docs.airr-community.org/en/latest/datarep/rearrangements.html).
2. **Integrate.** Combine spatially-resolved immune receptors and a Seurat spatial
   transcriptomics object into one unified SIRUP object. The SIRUP object can hold
   multiple spatial tissue sections, as well as associated non-spatial AIRR-seq data
   from bulk or single-cell experiments.
3. **Analyze and visualize.** Select regions of interest through *in silico* micro-dissection,
   measure clonal dominance and diversity within ROIs, test spatial patterns,
   build lineage trees and visualize them over tissue, and more!

<div markdown="1" style="background: #ffedcc; padding: 12px; margin-bottom: 24px; line-height: 24px;">
  <p style="margin: -12px -12px 12px; padding: 6px 12px; background: #d6b676; color: #ffffff; font-weight: 700;">Did you know?</p>
  Maple syrup has been spelled [maple sirup](https://maplesyruphistory.com/2022/09/04/the-origins-of-the-spelling-of-maple-sirup-with-an-i/) at different points of U.S. history.
</div>

## Getting started

Start with [Get started with SIRUP](vignettes/get-started.md). It takes about 15
minutes and uses example data that comes with the package. You build a SIRUP
object, plot receptors on the tissue, and measure clonal dominance in a region.

Install SIRUP from GitHub with BiocManager, which will also install dependencies:

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("hoehnlab/SIRUP")
```

Use R 4.1 or later.

## Contents

**Tutorials**

* [Process demultiplexed reads](vignettes/processing.md). Turn
  demultiplexed reads into a table of consensus BCR or TCR sequences.
* [Annotate tissue ROIs and analyze local clonal
  structure](vignettes/regions.md). Compare the repertoires of regions such as
  germinal centers.
* [Call clones and build lineage trees on the
  tissue](vignettes/clones-and-trees.md). Build trees with Dowser and show them
  on the H&E image.
* [Analyze several sections with bulk AIRR-seq and
  blood](vignettes/multisample.md). Combine sections and combine spatial and bulk receptor data in
  one object.

**Guides**

* [How the SIRUP object works](vignettes/sirup-object.md). 
* [Choose the right spatial statistic for your
  question](vignettes/choose-spatial-statistic.md). Match your question to a
  function.
* [Use QuPath annotations to define ROIs](vignettes/qupath-annotations.md). Export outlines from QuPath and add them
  to a SIRUP object.
* [Use SIRUP with your LLM agent](vignettes/llm-agents.md).

## For AI agents and LLMs

SIRUP provides the following resources for AI agents. See the [guide](vignettes/llm-agents.md) for more details.

* [llms.txt](https://github.com/hoehnlab/SIRUP/blob/main/llms.txt) is a short Markdown index of the SIRUP documentation, with links.
* [llms-full.txt](https://github.com/hoehnlab/SIRUP/blob/main/llms-full.txt) holds the whole documentation in one text file.
* [llms.json](https://github.com/hoehnlab/SIRUP/blob/main/llms.json) is the same index as JSON. It lists each tutorial, guide, and function, with a summary and links.
* [llms-full.json](https://github.com/hoehnlab/SIRUP/blob/main/llms-full.json) holds the whole documentation as JSON, with the full text of each page.
* The **SIRUP skill** teaches an agent how to use SIRUP: how to build an object, which function answers which question, and which checks to run. Find the skill in `inst/skills/sirup`, and install it as shown below.

**Set up the skill**

Install SIRUP. Then use the `installSIRUPSkill()` function to install the skill for your agent. Currently,  `installSIRUPSkill` is set up to work for Claude Code, Codex, or Antigravity. Other types of agents can be set up by manually installing the skill.

```r
library(SIRUP)
installSIRUPSkill("claude")   # or "codex" or "antigravity"
```

By default, this installs the skill for the agent for all projects. To install the skill for one project only, give 
the skills folder of the project, for example `installSIRUPSkill("claude", install_dir = ".claude/skills")`.
Then start a new session and ensure that the SIRUP skill is available.

## Contact

If you need help or have a question, contact the [Hoehn Lab](https://sites.dartmouth.edu/hoehn/).
To report a bug or ask for a feature, open an issue on the
[issue tracker](https://github.com/hoehnlab/SIRUP/issues).

## Dependencies

**Depends:** R (>= 4.1.0)  
**Imports:** airr, airrblast, alakazam, cli, cowplot, data.table, dowser, dplyr, ggplot2, ggtree, grDevices, grid, igraph, kedd, KernSmooth, magrittr, Matrix, methods, parallel, patchwork, pbapply, RANN, rlang, scales, scoper, Seurat, SeuratObject, sf, shiny, spdep, stats, stringr, tibble, tidyr, tools, utils  
**Suggests:** ape, arrow, BiocNeighbors, knitr, nabor, pals, qpdf, rgl, rmarkdown, shazam, spatstat.explore, spatstat.geom, spatstat.random, testthat (>= 3.0.0), withr

## Authors

[Hunter J. Melton](mailto:Hunter.J.Melton@dartmouth.edu)<br>
Chinmay M. Vaidya<br>
Apekshya Panda<br>
Jessie J. Fielding<br>
Elizabeth Sergison<br>
Fred W. Kolling<br>
Kevin C. O'Connor<br>
[Kenneth B. Hoehn](mailto:Kenneth.B.Hoehn@dartmouth.edu)<br>


## Citing SIRUP

Melton HJ, Vaidya CM, Panda A, Fielding JJ, Sergison E, Kolling FW, O'Connor KC,
Hoehn KB. SIRUP, a toolkit for spatial B and T cell receptor repertoires, reveals
differences between normal and ectopic germinal centers in myasthenia gravis.
bioRxiv.

## License

MIT. See the
[LICENSE file](https://github.com/hoehnlab/SIRUP/blob/main/LICENSE).
