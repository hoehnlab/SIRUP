# Render the package vignettes into the read-the-docs folders.
#
# The vignettes in `vignettes/` are the single source. This script writes a
# GitHub markdown copy (with figure folders) to `markdown/` and a standalone
# HTML copy to `html/`. Run it from the package root after you change a
# vignette:
#
#   Rscript docs/readthedocs/render_docs.R
#
# The vignettes need the suggested packages (airr, Seurat, scoper, shazam,
# dowser, ggtree, spdep, spatstat, alakazam, sf, rgl). The tree sections also
# need IMGT germline references: set SIRUP_IMGT_DIR to the human VDJ folder,
# or those chunks show code only. The pages for rrad and SPaCeD are not
# vignettes. Their sources stay in this folder, and this script renders them
# too.

devtools::load_all(".", quiet = TRUE)

docs_dir <- file.path("docs", "readthedocs")
sources <- c(
  list.files("vignettes", pattern = "[.]Rmd$", full.names = TRUE),
  list.files(docs_dir, pattern = "[.]Rmd$", full.names = TRUE)
)
# Each page knits in its own folder, as in R CMD build, so relative paths such
# as figures/polygon-app.jpg resolve.

# Screenshots that the vignettes include. The markdown copies need them beside
# them; the HTML copies embed them.
dir.create(file.path(docs_dir, "markdown", "figures"), recursive = TRUE, showWarnings = FALSE)
file.copy(
  list.files(file.path("vignettes", "figures"), full.names = TRUE),
  file.path(docs_dir, "markdown", "figures"),
  overwrite = TRUE
)

for (input in sources) {
  name <- tools::file_path_sans_ext(basename(input))
  message("Rendering ", name)

  rmarkdown::render(
    input,
    output_format = rmarkdown::github_document(html_preview = FALSE),
    output_dir = file.path(docs_dir, "markdown"),
    envir = new.env(),
    quiet = TRUE
  )

  rmarkdown::render(
    input,
    output_dir = file.path(docs_dir, "html"),
    envir = new.env(),
    quiet = TRUE
  )
}
