# docs/render_vignettes.R - vignettes/*.Rmd -> docs/vignettes/*.md (+ figures).
#
# STATUS: run on a small fixture package (R 4.6.1, Nov 2026): figures, local images and
# vignette() links work. Not yet run on a package with heavy vignettes (Seurat and so on).
# It renders rmarkdown::github_document, which has been used in practice.
#
# vignettes/*.Rmd stays the single source. Chunks that call tools the machine lacks should
# already carry eval = FALSE in the Rmd; this script never edits the Rmd.

# only: NULL renders every vignette. A character vector of file stems renders just those.
render_vignettes <- function(config, only = NULL) {
  pkg <- config$package
  out_dir <- file.path("docs", "vignettes"); dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  if (requireNamespace("devtools", quietly = TRUE)) devtools::load_all(".", quiet = TRUE) else library(pkg, character.only = TRUE)

  sources <- unlist(lapply(config$vignette_sources, function(d) list.files(d, pattern = "[.](Rmd|qmd)$", full.names = TRUE)))
  # Leave out vignettes that must not become site pages (vignette_exclude in docs-config.yml).
  sources <- sources[!(tools::file_path_sans_ext(basename(sources)) %in% unlist(config$vignette_exclude))]
  titles <- vapply(sources, function(s) {
    t <- rmarkdown::yaml_front_matter(s)$title
    if (is.null(t)) tools::file_path_sans_ext(basename(s)) else t
  }, "")
  names(titles) <- tools::file_path_sans_ext(basename(sources))
  # Titles of all vignettes stay available, so links between vignettes resolve in a partial run.
  if (!is.null(only)) sources <- sources[names(titles) %in% only]

  for (src in sources) {
    stem <- tools::file_path_sans_ext(basename(src))
    message("Rendering ", stem)
    # Each page knits in its own folder (as R CMD build does) so relative paths resolve.
    rmarkdown::render(src, output_format = rmarkdown::github_document(html_preview = FALSE),
                      output_dir = out_dir, envir = new.env(), quiet = TRUE)
    md_path <- file.path(out_dir, paste0(stem, ".md"))
    md <- paste(readLines(md_path, warn = FALSE), collapse = "\n")
    # rmarkdown writes figure links as absolute paths when output_dir is set. Make them
    # relative to the page, or the links break on the built site.
    md <- gsub(paste0(normalizePath(out_dir), "/"), "", md, fixed = TRUE)
    md <- rewrite_html_images(md)
    if (!isFALSE(config$extras$link_rewrite)) md <- rewrite_vignette_refs(md, pkg, titles)
    writeLines(md, md_path)
    copy_referenced_assets(md, dirname(src), out_dir)
  }
  # remove pages whose source vignette is gone (full runs only)
  if (!is.null(only)) return(invisible(NULL))
  stale <- setdiff(list.files(out_dir, pattern = "[.]md$"), paste0(names(titles), ".md"))
  if (length(stale)) {
    message("removing stale vignette pages: ", paste(stale, collapse = ", "))
    file.remove(file.path(out_dir, stale))
    unlink(file.path(out_dir, paste0(sub("[.]md$", "", stale), "_files")), recursive = TRUE)
  }
  invisible(NULL)
}

# knitr writes a raw HTML <img> tag when a chunk sets out.width. MkDocs does not rewrite the
# path in raw HTML, so the image breaks on the built site (a page at /vignettes/x/ looks for
# the file in the wrong folder). Change each relative <img> to Markdown image syntax, which
# MkDocs does rewrite. The width attribute is dropped. The theme limits images to the page width.
rewrite_html_images <- function(md) {
  replace_all <- get("replace_all", mode = "function")   # defined in build_topics.R
  replace_all(md, "(<img\\s[^>]*>)", function(c) {
    tag <- c[1]
    src <- sub('^.*\\ssrc="([^"]*)".*$', "\\1", tag)
    if (identical(src, tag) || grepl("^(https?:|data:|/)", src)) return(tag)
    alt <- if (grepl('\\salt="', tag)) sub('^.*\\salt="([^"]*)".*$', "\\1", tag) else ""
    sprintf("![%s](%s)", alt, src)
  })
}

# `vignette("x", package = "pkg")` -> [Title](x.md). Unknown names are left alone and reported.
rewrite_vignette_refs <- function(md, pkg, titles) {
  pat <- paste0("`?vignette\\(\\s*[\"']([^\"']+)[\"']\\s*(?:,\\s*package\\s*=\\s*[\"']", pkg, "[\"']\\s*)?\\)`?")
  replace_all <- get("replace_all", mode = "function")   # defined in build_topics.R
  replace_all(md, pat, function(c) {
    if (c[1] %in% names(titles)) sprintf("[%s](%s.md)", titles[[c[1]]], c[1])
    else { message("  unresolved vignette reference: ", c[1]); sprintf("`vignette(\"%s\")`", c[1]) }
  })
}

# Static images referenced as ![](figures/x.png) or <img src="figures/x.png"> sit next to the Rmd;
# copy them beside the generated page so the relative links keep working.
copy_referenced_assets <- function(md, src_dir, out_dir) {
  refs <- c(regmatches(md, gregexpr("!\\[[^]]*\\]\\(([^) ]+)", md, perl = TRUE))[[1]],
            regmatches(md, gregexpr("<img[^>]*src=\"([^\"]+)\"", md, perl = TRUE))[[1]])
  refs <- unique(sub("^.*?(?:\\(|src=\")", "", refs, perl = TRUE))
  refs <- sub("[\")]$", "", refs)
  for (r in refs[!grepl("^(https?:|data:)", refs)]) {
    from <- file.path(src_dir, r); to <- file.path(out_dir, r)
    if (file.exists(from) && !file.exists(to)) { dir.create(dirname(to), recursive = TRUE, showWarnings = FALSE); file.copy(from, to) }
  }
}
