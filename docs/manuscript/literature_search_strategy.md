# Literature Search Strategy for the SIRAP Introduction

Goal: build the introduction quickly by searching for claims, not papers.

The output should be a short claim matrix that supports 4-6 introduction
paragraphs. Each paragraph gets one job, 2-4 high-quality citations, and a
clear handoff into the SIRAP manuscript.

## Fastest Successful Workflow

### 1. Define the Introduction Claim Ladder

Use this ladder unless the paper's story changes:

1. Spatial omics preserves tissue context, but transcriptomics alone does not
   reveal the receptor lineages that define adaptive immune responses.
2. AIRR-seq reveals BCR/TCR clonal structure, selection, somatic mutation,
   class switching, and lineage history, but conventional AIRR-seq usually
   loses tissue location.
3. Spatial AIRR-seq (spAIRR) should connect receptor sequence identity to
   tissue architecture, local gene expression, and cellular neighborhoods.
4. Ectopic germinal centers in myasthenia gravis are a strong biological test
   case because thymic immune niches may shape autoreactive B-cell responses.
5. Existing workflows leave a practical gap between raw spAIRR reads,
   coordinate-aware AIRR tables, and spatial receptor analysis.
6. SIRAP fills this gap in two parts: processing spAIRR data with RAD
   interfaces and pRESTO/IgBLAST pipelines, then visualizing/analyzing spAIRR
   data with spatial clone, lineage, region, and pairing tools.

### 2. Search by Evidence Buckets

Search each bucket separately. Stop once each bucket has a few strong papers
plus one or two recent reviews.

#### Bucket A: Spatial Omics Context

Purpose: establish why spatial context matters.

Search strings:

- `"spatial transcriptomics" immune tissue review`
- `"spatial omics" immunology review`
- `"spatial transcriptomics" "germinal center"`
- `"Visium HD" spatial transcriptomics immune`

Preferred citations:

- Method/review papers from Nature Methods, Nature Biotechnology, Cell,
  Science, Genome Research, Nature Reviews, Immunity.
- Papers that explicitly discuss tissue architecture, immune neighborhoods, or
  spatial cell-cell interactions.

#### Bucket B: AIRR-Seq and Receptor Lineage Context

Purpose: explain what receptor sequencing adds.

Search strings:

- `"AIRR-seq" review BCR TCR repertoire`
- `"adaptive immune receptor repertoire sequencing" clonal lineage review`
- `"B cell receptor repertoire" somatic hypermutation lineage germinal center`
- `"T cell receptor repertoire" spatial transcriptomics`

Preferred citations:

- AIRR standards/repertoire review papers.
- BCR lineage and germinal center repertoire papers.
- Methods papers for clone calling or lineage reconstruction only if directly
  needed.

#### Bucket C: Spatial Receptor Sequencing / spAIRR

Purpose: establish the gap SIRAP addresses.

Search strings:

- `"spatial" "B cell receptor" sequencing`
- `"spatial" "T cell receptor" sequencing`
- `"spatial immune repertoire" sequencing`
- `"spatial AIRR" OR "spAIRR"`
- `"Visium" "BCR" "TCR"`
- `"spatial transcriptomics" "immune receptor"`

Preferred citations:

- Papers that recover paired gene expression and immune receptor sequence with
  spatial coordinates.
- Methods that explicitly handle BCR/TCR in tissue.
- SIRAP, RAD, and SPaCeD preprints once metadata is verified.

#### Bucket D: Ectopic Germinal Centers in MG

Purpose: justify the disease biology.

Search strings:

- `"myasthenia gravis" thymus "germinal center"`
- `"myasthenia gravis" ectopic germinal center thymus`
- `"myasthenia gravis" thymus B cell receptor repertoire`
- `"acetylcholine receptor" thymus germinal center myasthenia gravis`
- `"myasthenia gravis" thymectomy B cells germinal centers`

Preferred citations:

- Foundational MG thymus ectopic GC papers.
- Recent single-cell/spatial MG thymus papers.
- Repertoire/PBMC longitudinal work, including the PNAS study already noted.

#### Bucket E: Canonical Tonsil GC Comparator

Purpose: justify tonsil as a regular/canonical GC benchmark.

Search strings:

- `"human tonsil" germinal center spatial transcriptomics`
- `"tonsil" "B cell receptor" spatial`
- `"human tonsil" "germinal center" "single-cell"`
- `"Engblom" tonsil germinal center BCR`

Preferred citations:

- The Engblom tonsil paper/data source.
- Human tonsil spatial or single-cell GC references with strong annotation of
  dark zone, light zone, follicular dendritic cells, Tfh, and plasma cells.

### 3. Use a Claim Matrix, Not a Bibliography

For each paper, record only what we need:

- The claim it supports.
- Why it matters for SIRAP.
- Whether it is background, method precedent, biological precedent, or direct
  comparator.
- One sentence we might write in the introduction.

Use `intro_claim_matrix.tsv` for this.

### 4. Triage Rules

Read in this order:

1. Title and abstract.
2. Figures and captions.
3. Introduction last paragraph and Discussion first paragraph.
4. Methods only if the paper is a direct technical comparator.
5. Full paper only for the 10-15 citations likely to survive into the draft.

Reject or deprioritize papers that:

- Are only transcriptomics and do not help frame the receptor gap.
- Are only repertoire sequencing and do not help with tissue/spatial context.
- Are disease-adjacent but not MG, ectopic GCs, or germinal center biology.
- Are methodologically impressive but not necessary for the introduction.

### 5. Introduction Assembly

After the matrix has enough evidence, write a rough introduction in this shape:

Paragraph 1: spatial omics changed tissue biology, but spatial transcriptomics
alone does not recover adaptive receptor lineages.

Paragraph 2: AIRR-seq captures clonal identity and lineage history, but standard
AIRR workflows lose spatial context.

Paragraph 3: spAIRR data should connect receptor identity, local gene
expression, tissue regions, and lineage relationships, but practical processing
and analysis tooling is fragmented.

Paragraph 4: MG thymic ectopic GCs provide a biological case where spatial
receptor analysis can test local clonal expansion and tissue-periphery
connections.

Paragraph 5: SIRAP contribution: process spAIRR data with RAD interfaces and
downstream AIRR pipelines; visualize/analyze spAIRR data with spatial receptor,
clone, lineage, region, and SPaCeD-style pairing analyses.

## What Codex Can Do Well

Codex can successfully:

- Run targeted web literature searches.
- Build and maintain the claim matrix.
- Summarize abstracts and papers into manuscript-relevant claims.
- Draft introduction paragraphs from the claim matrix.
- Keep RAD, SIRAP, and SPaCeD framing distinct.
- Track uncertainty, missing metadata, and papers that require manual download.

## What Would Make It Better

Best setup:

- A Zotero export (`.bib` or RIS) for papers the group already trusts.
- PDFs or links for paywalled/key papers in a local folder.
- Permission to use command-line literature APIs if needed, such as PubMed,
  Europe PMC, Crossref, Semantic Scholar, or bioRxiv APIs.
- A short list of "must cite" papers from Ken/Chinmay/Hunter.

Minimum setup:

- Keep using web search plus manually provided PDFs/links for blocked pages.
- Add search results and summaries to `intro_claim_matrix.tsv`.

