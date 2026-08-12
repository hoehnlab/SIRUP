# SIRAP Genome Research Spatial Omics Manuscript Plan

Working target for the Genome Research spatial omics special issue.

- Current planning date: 2026-07-08
- First rough draft target: 2026-09-01
- Submission target: 2026-11-01
- CFP link to verify manually: https://www.cshlpress.com/email_news/GRcfp626/spatial_omics.html
- PBMC/MG study DOI to verify manually: https://www.pnas.org/doi/10.1073/pnas.2007206117
- SPaCeD preprint to verify manually: https://www.biorxiv.org/content/10.64898/2026.04.19.719512v1.full
- RAD preprint to verify manually: https://www.biorxiv.org/content/10.64898/2026.05.29.728810v1

## One-Sentence Pitch

SIRAP is an R toolkit for spatial AIRR-seq (spAIRR) data that connects receptor sequences to tissue coordinates, then supports spatial visualization and analysis of BCR/TCR clones, germinal center structure, receptor lineages, and cross-compartment clonal tracking.

## Method Framing

SIRAP should be framed as a two-part spAIRR toolkit:

1. Process spAIRR data.
   - SIRAP interfaces with RAD through R wrappers.
   - RAD is the demultiplexing/read-structure package that works with SIRAP, not a method we should present as internal to SIRAP.
   - SIRAP provides pre-built downstream processing pipelines after RAD, including pRESTO/IgBLAST-based AIRR generation, QC, barcode extraction, Visium HD barcode mapping, and consensus/top-sequence filtering.
2. Visualize and analyze spAIRR data.
   - SIRAP provides spatial overlays, receptor/gene colocalization, clone dominance summaries, region/polygon workflows, lineage tree/spatial plot pairing, and spatial statistics.
   - SPaCeD is its own method for spatial heavy/light pairing. In SIRAP, we re-implemented/integrated the SPaCeD-style analysis surface so it can be used alongside the rest of the spAIRR workflow.

Manuscript language should keep these boundaries clear:

- RAD: companion demultiplexing package used by SIRAP wrappers.
- SIRAP: R package for processing, integrating, visualizing, and analyzing spAIRR data.
- SPaCeD: separate spatial pairing method re-implemented/integrated in SIRAP.

## Core Stories

### Story 1: Ectopic Versus Canonical Germinal Centers

Compare ectopic thymic germinal centers from the Visium MG thymus data against canonical tonsil germinal centers from the Engblom tonsil data using matched spatial transcriptomic and receptor-level summaries.

Primary questions:

- Are ectopic thymic GCs transcriptionally organized like canonical tonsil GCs?
- Do thymic GCs show BCR clonal expansion, SHM/lineage structure, and spatial clonal dominance consistent with local antigen-driven selection?
- Are TCRs enriched near B-cell/GC regions, AChR-expressing or antigen-presenting epithelial populations, or other disease-relevant niches?

Minimum viable results:

- Harmonized GC annotations or region definitions for tonsil and thymus.
- BCR clone abundance, diversity, SHM, V/J usage, and spatial dominance summaries in GC versus non-GC regions.
- Spatial colocalization statistics linking receptors, GC markers, and relevant tissue states.
- A clear statement of similarities and differences between ectopic thymic and canonical tonsil GCs.

Nice-to-have results:

- TCR clone localization or diversity summaries if quality and counts support it.
- Heavy/light pairing for major thymic BCR clones using the SIRAP re-implementation/integration of SPaCeD.
- BCR lineage trees with spatially clustered subclades for representative clones.

### Story 2: Thymus-PBMC Clonal Connections in MG

Connect spatial thymic BCR clones to PBMC repertoires sampled at 0 and 12 months using the Visium MG data plus the BCR data from the MG PNAS study.

Primary questions:

- Do BCR clones observed in thymic GCs also appear in peripheral blood?
- Are thymus-linked clones persistent, expanded, class-switched, or mutated in PBMCs at 0 and/or 12 months?
- Do thymus-localized clonal lineages provide evidence for exchange between ectopic lymphoid tissue and circulation?

Minimum viable results:

- Reprocessed or imported PBMC AIRR tables in a format comparable to SIRAP thymus output.
- Shared clone calling across thymus and PBMC using consistent heavy-chain criteria.
- Clone-sharing tables for thymus versus PBMC month 0 and thymus versus PBMC month 12.
- Representative lineage trees showing tissue/timepoint labels for shared clones.

Nice-to-have results:

- Quantification of directional lineage structure, such as thymus-interior nodes versus PBMC-derived tips.
- Separate summaries by isotype, mutation burden, and GC region of origin.
- Comparison against permutation or donor-mismatched clone-sharing controls.

## Figure Skeleton

### Figure 1: SIRAP Workflow and Data Overview

Purpose: establish the method and the datasets.

Panels:

- Experimental/computational schematic: Visium HD spatial transcriptome plus long-read BCR/TCR processing.
- SIRAP processing arm: R wrappers around RAD demultiplexing followed by pRESTO/IgBLAST-based AIRR processing, QC, barcode mapping, and consensus filtering.
- SIRAP analysis arm: receptor overlays, spatial statistics, GC region summaries, lineage/spatial plots, cross-compartment clone tracking, and SPaCeD-style pairing.
- Dataset table: MG thymus samples, Engblom tonsil samples, PBMC timepoints, receptor modalities, and available spatial transcriptomics.
- Representative spatial tissue view with BCR/TCR-positive spots.

Existing anchors:

- `SIRAP/inst/scripts/visium_hd_bcr_tcr_walkthrough.R`
- `SIRAP/docs/readthedocs/visium-hd-processing.Rmd`
- `spatialExpressionPlot()`
- `plotReceptorsExpression()`

### Figure 2: Ectopic Thymic GCs Show Spatially Localized BCR Clones

Purpose: show that thymic GCs contain expanded, spatially structured BCR clones.

Panels:

- GC/mantle/other region polygons over thymus tissue.
- Clone dominance per GC polygon.
- Join-count or Moran-style clone localization statistics.
- Representative BCR clone spatial map and lineage tree.

Existing anchors:

- `selectSpatialPolygons()`
- `spatialCloneDominance()`
- `spatialCloneJoinCount()`
- `plotSpatialCloneTreePairs()`
- `SIRAP/docs/readthedocs/spatial-clonal-dominance.Rmd`
- `SIRAP/docs/readthedocs/spatial-clones-and-trees.Rmd`

### Figure 3: Canonical Tonsil Versus Ectopic Thymic GC Comparison

Purpose: make the comparative biology claim.

Panels:

- Matched GC marker/module scores in tonsil and thymus.
- BCR clone diversity/dominance distributions across GC regions.
- SHM and isotype/V-gene summaries, if available in both datasets.
- Spatial colocalization of BCRs with GC markers and local cell states.

Existing anchors:

- `spatialMoran()`
- `spatialLocalMoran()`
- `spatialRipleyCrossK()`
- `plotSpatialLocalMoran()`
- `plotRipleyCrossK()`

Needed:

- Import/harmonization script for Engblom tonsil data.
- Shared region/annotation vocabulary for tonsil and thymus.
- Shared receptor QC table definitions.

### Figure 4: TCR and Tissue-Niche Context of Ectopic GCs

Purpose: add immune-context depth without letting TCR become the schedule risk.

Panels:

- TCR-positive spots relative to GC B-cell markers, mTEC/cTEC/CHRNA1-positive regions, or antigen-presentation markers.
- TCR clone localization summary if counts are sufficient.
- BCR/TCR overlap or neighborhood enrichment around GC regions.

Existing anchors:

- `plotReceptorsExpression()`
- `spatialRipleyCrossK()`
- `spatialExpressionBlendPlot()`
- `SIRAP/docs/readthedocs/spatial-expression-and-colocalization.Rmd`

Decision gate:

- By 2026-08-05, decide whether TCR supports a main figure, supplement-only result, or is dropped from the first manuscript story.

### Figure 5: Thymus-PBMC Clonal Connections Over Time

Purpose: show disease-relevant exchange between ectopic thymic GCs and circulation.

Panels:

- Clone-sharing network or UpSet-style summary across thymus, PBMC month 0, and PBMC month 12.
- Abundance and mutation/isotype summaries for shared versus thymus-only/PBMC-only clones.
- Representative shared-clone lineage tree labeled by compartment/timepoint and spatial GC origin.
- Spatial map of the thymic GC location for the shared representative clone.

Existing anchors:

- `clusterTreesTips()`
- `plotSpatialCloneTreePairs()`
- AIRR processing helpers in `processing.R` and clone matrix helpers in `clones.R`

Needed:

- PBMC data acquisition and import.
- Cross-sample clone-calling script.
- Shared-clone QC and false-match controls.

## First-Draft Outline

1. Abstract
2. Introduction
   - Why immune receptor identity matters in spatial omics.
   - Why ectopic GCs in MG are a compelling biological test case.
   - What spAIRR adds beyond transcriptomics-only spatial analysis.
   - How SIRAP sits between RAD-based processing and spatial receptor analysis, with SPaCeD treated as a distinct pairing method.
3. Results
   - SIRAP processes spAIRR data by wrapping RAD and downstream pRESTO/IgBLAST workflows.
   - SIRAP maps BCR/TCR sequences into Visium HD tissue coordinates for visualization and analysis.
   - Ectopic thymic GCs contain spatially localized BCR clonal expansions.
   - Ectopic thymic GCs resemble and differ from canonical tonsil GCs.
   - TCRs and tissue niches provide local immune context.
   - Thymic GC clones connect to circulating PBMC repertoires over time.
4. Discussion
   - Spatial receptor data as lineage evidence for local immune responses.
   - Ectopic lymphoid tissue in MG as a source or reservoir of autoreactive clones.
   - Methodological limitations: capture efficiency, clone-calling thresholds, spatial resolution, cross-sample clone matching.
5. Methods
   - Visium HD processing and RAD demultiplexing through SIRAP wrappers.
   - AIRR processing, pRESTO/IgBLAST, clone calling.
   - Spatial barcode mapping and Seurat integration.
   - Spatial statistics and GC polygon analysis.
   - SIRAP re-implementation/integration of SPaCeD-style heavy/light pairing.
   - BCR lineage tree reconstruction.
   - PBMC repertoire import and cross-compartment clone matching.

## Analysis Checklist

### Data Inventory

- Confirm MG thymus samples and output directories for M1 and M5.
- Confirm which sample(s) contain high-confidence ectopic GCs.
- Locate/import Engblom tonsil spatial transcriptomic and BCR/TCR data.
- Locate/import PBMC AIRR/BCR data from the MG PNAS study.
- Create a table with sample IDs, tissue, disease/control status, timepoint, spatial platform, receptor modality, and file paths.

### Thymus SIRAP Outputs

- Re-run or validate final BCR/TCR AIRR tables.
- Confirm BCR/TCR barcode mapping from 2um to 8um/16um.
- Generate receptor-positive spatial overlays.
- Call BCR clones with a fixed, documented threshold.
- Call TCR clones if supported.
- Create GC polygon/annotation tables.
- Calculate clone dominance and join-count localization in GC regions.
- Build representative BCR lineage trees.

### Tonsil Comparison

- Standardize tonsil cell/spot annotations with thymus annotations.
- Define canonical GC regions or GC-positive spots.
- Compute matched receptor and transcriptomic summaries.
- Compare GC versus non-GC distributions within each tissue.
- Compare canonical tonsil GC versus ectopic thymic GC effect sizes.

### PBMC Connection

- Convert PBMC repertoires to AIRR-compatible schema.
- Harmonize V/J gene names and junction definitions with thymus data.
- Decide clone definition for cross-compartment matching.
- Quantify shared clones by timepoint.
- Test robustness to stricter/laxer clone thresholds.
- Build trees for top shared clones.
- Link shared clones back to thymic spatial regions.

## Timeline to Rough Draft

### 2026-07-08 to 2026-07-15: Inventory and Reproducibility

- Freeze file paths, sample IDs, and available metadata.
- Create manuscript analysis scripts rather than relying only on walkthroughs.
- Produce a dataset inventory table.
- Identify which existing SIRAP functions need small hardening for paper outputs.

Deliverable:

- `data_inventory.tsv`
- Empty or first-pass figure output directories.
- Reproducible script skeletons for Figures 1, 2, and 5.

### 2026-07-16 to 2026-07-31: Thymus GC Story

- Regenerate thymus BCR/TCR spatial overlays.
- Define/validate GC polygons or annotations.
- Produce clone dominance, join-count, and representative lineage-tree outputs.
- Pick the first two candidate Figure 2 clones.

Deliverable:

- Draft Figure 2 panels.
- Results notes for thymic GC clonal localization.

### 2026-08-01 to 2026-08-14: Tonsil Comparison

- Import and harmonize Engblom tonsil data.
- Run matched GC/non-GC receptor and transcriptomic summaries.
- Decide whether TCR is main text or supplement.

Deliverable:

- Draft Figure 3 panels.
- Results notes comparing ectopic and canonical GCs.

### 2026-08-15 to 2026-08-25: PBMC Linkage Story

- Import MG PBMC repertoires.
- Run cross-compartment clone matching.
- Build clone-sharing summaries and representative trees.
- Pick the strongest Figure 5 examples.

Deliverable:

- Draft Figure 5 panels.
- Results notes for thymus-PBMC clonal connections.

### 2026-08-26 to 2026-09-01: Assemble Bad First Draft

- Write all Results subsections around existing panels.
- Add placeholder Introduction, Discussion, and Methods.
- Mark every missing number and figure panel as `TODO`.
- Circulate rough draft.

Deliverable:

- First complete manuscript draft, even if ugly.

## Immediate Next Actions

1. Create `SIRAP/inst/scripts/manuscript/` for figure-specific reproducible scripts.
2. Create a data inventory TSV template.
3. Start Figure 1 and Figure 2 scripts from the existing Visium HD walkthrough.
4. Determine where Engblom tonsil and MG PBMC data live locally or on the cluster.
5. Decide the canonical cross-sample BCR clone definition before running PBMC matching.

## Risks and Decisions

- CFP deadline and article type details need manual verification because the CFP page blocked automated retrieval here.
- The PBMC DOI/source needs manual verification because the DOI page blocked automated retrieval here.
- The SPaCeD and RAD preprints need manual metadata verification because the bioRxiv pages blocked automated retrieval here.
- Engblom tonsil data availability may determine whether Story 1 is a true matched comparison or a more limited external benchmark.
- TCR should not block the first draft; it graduates to main text only if counts, clone calls, and spatial signal are strong by 2026-08-05.
- Cross-sample clone matching can overcall convergent public clones; we need sensitivity analyses with stricter matching and donor/sample controls where possible.
