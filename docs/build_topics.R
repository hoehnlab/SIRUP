# docs/build_topics.R - man/*.Rd -> docs/topics/*.md using only base R tools + pandoc.
#
# STATUS: run against the 71 Rd files of SIRUP with R 4.6.1 and pandoc 3.12 (Nov 2026).
# Still read two or three generated pages after the first run on a new package.
# See references/validation.md.
#
# Pipeline per Rd file:
#   1. rewrite \link{} / \link[=x]{} / \link[pkg:x]{} / \code{\link{}} into \href{url}{text}
#      (this is what markr could not do; \href renders to a plain <a> in every R version)
#   2. tools::Rd2HTML()  -> full HTML page. Do NOT pass fragment = TRUE: in R that flag
#      means "the input is an Rd fragment", so parse_Rd rejects \name and \alias.
#   3. tidy the HTML (keep <main> only, drop R's header table and version footer,
#      arguments table -> <dl>)
#   4. pandoc HTML -> Markdown; MkDocs `def_list` renders <dl> terms like Dowser's pages

BASE_PKGS <- c("base", "stats", "utils", "graphics", "grDevices", "methods", "tools", "grid", "parallel")

# The base-R manual names each page after the Rd FILE, not the alias. For example, the topic
# saveRDS lives in readRDS.html. Look up the Rd file that holds an alias in a base package.
.base_alias_cache <- new.env()
base_rd_stem <- function(pkg, topic) {
  if (is.null(.base_alias_cache[[pkg]])) {
    db <- tools::Rd_db(pkg)
    stems <- character(0)
    for (f in names(db)) {
      rd <- db[[f]]
      tags <- vapply(rd, function(x) { t <- attr(x, "Rd_tag"); if (is.null(t)) "" else t }, "")
      aliases <- vapply(which(tags == "\\alias"), function(i) paste(unlist(rd[[i]]), collapse = ""), "")
      stems[aliases] <- sub("[.]Rd$", "", f)
    }
    .base_alias_cache[[pkg]] <- stems
  }
  stem <- .base_alias_cache[[pkg]][topic]
  if (is.na(stem)) topic else unname(stem)
}
base_manual_url <- function(pkg, topic) {
  sprintf("https://stat.ethz.ch/R-manual/R-devel/library/%s/html/%s.html", pkg, base_rd_stem(pkg, topic))
}

replace_all <- function(txt, pattern, fn) {
  m <- gregexpr(pattern, txt, perl = TRUE)[[1]]
  if (m[1] == -1L) return(txt)
  starts <- as.integer(m); lens <- attr(m, "match.length")
  cs <- attr(m, "capture.start"); cl <- attr(m, "capture.length")
  out <- character(0); pos <- 1L
  for (i in seq_along(starts)) {
    caps <- vapply(seq_len(ncol(cs)), function(j) substr(txt, cs[i, j], cs[i, j] + cl[i, j] - 1L), "")
    out <- c(out, substr(txt, pos, starts[i] - 1L), fn(caps))
    pos <- starts[i] + lens[i]
  }
  paste0(c(out, substr(txt, pos, nchar(txt))), collapse = "")
}

# Where does a link to <topic> in <pkg> go?  Returns a URL, or NA if it cannot be resolved.
resolve_link <- function(topic, pkg, own_pkg, alias_map, external) {
  if (is.na(pkg) || pkg == "" || pkg == own_pkg) {
    stem <- alias_map[[topic]]
    if (!is.null(stem)) return(paste0(stem, ".md"))   # same folder: docs/topics/
    # Not an alias of this package. R resolves an unqualified \link{} on the search path,
    # so look for a base-R topic (for example saveRDS) before giving up.
    hit <- sub("^package:", "", utils::find(topic))
    base_hit <- intersect(hit, BASE_PKGS)
    if (length(base_hit)) return(base_manual_url(base_hit[1], topic))
    return(NA_character_)
  }
  tmpl <- external[[pkg]]
  if (!is.null(tmpl)) return(gsub("{topic}", topic, tmpl, fixed = TRUE))
  if (pkg %in% BASE_PKGS) return(base_manual_url(pkg, topic))
  sprintf("https://rdrr.io/cran/%s/man/%s.html", pkg, topic)   # best effort; override in docs-config.yml
}

rewrite_rd_links <- function(txt, own_pkg, alias_map, external, warn_fn) {
  link <- function(topic, pkg, label) {
    url <- resolve_link(topic, pkg, own_pkg, alias_map, external)
    if (is.na(url)) {
      warn_fn(sprintf("unresolved link: %s%s", if (is.na(pkg)) "" else paste0(pkg, ":"), topic))
      return(label)
    }
    sprintf("\\href{%s}{%s}", url, label)
  }
  # order matters: most specific first. \code{\link{x}} needs no special case: the \href
  # simply sits inside \code{}.
  txt <- replace_all(txt, "\\\\link\\[=([^\\]]+)\\]\\{([^}]*)\\}",          function(c) link(c[1], NA, c[2]))
  txt <- replace_all(txt, "\\\\link\\[([^\\]=:]+):([^\\]]+)\\]\\{([^}]*)\\}", function(c) link(c[2], c[1], c[3]))
  txt <- replace_all(txt, "\\\\link\\[([^\\]=:]+)\\]\\{([^}]*)\\}",          function(c) link(c[2], c[1], c[2]))
  txt <- replace_all(txt, "\\\\link\\{([^}]*)\\}",                           function(c) link(c[1], NA, c[1]))
  txt
}

tidy_html <- function(html) {
  html <- paste(html, collapse = "\n")
  # Rd2HTML returns a whole page (DOCTYPE, head, scripts). Keep only the <main> content.
  html <- sub("(?s)^.*?<main>", "", html, perl = TRUE)
  html <- sub("(?s)</main>.*$", "", html, perl = TRUE)
  html <- sub("(?s)^\\s*<table[^>]*>.*?</table>", "", html, perl = TRUE)                 # "topic {pkg} / R Documentation"
  html <- sub("(?s)<hr\\s*/?>\\s*<div[^>]*>\\s*\\[Package.*?</div>\\s*$", "", html, perl = TRUE)  # version footer
  # Arguments: <table><tr><td><code>x</code></td><td>desc</td></tr></table> -> <dl>
  html <- replace_all(html, "(?s)<table[^>]*>(.*?)</table>", function(c) {
    rows <- regmatches(c[1], gregexpr("(?s)<tr[^>]*>.*?</tr>", c[1], perl = TRUE))[[1]]
    items <- vapply(rows, function(r) {
      cells <- regmatches(r, gregexpr("(?s)<td[^>]*>.*?</td>", r, perl = TRUE))[[1]]
      cells <- gsub("(?s)^<td[^>]*>|</td>$", "", cells, perl = TRUE)
      if (length(cells) < 2) return("")
      sprintf("<dt>%s</dt><dd>%s</dd>", trimws(cells[1]), trimws(paste(cells[-1], collapse = " ")))
    }, "")
    paste0("<dl>", paste(items, collapse = "\n"), "</dl>")
  })
  html <- gsub("<h3", "<h2", html, fixed = TRUE); html <- gsub("</h3>", "</h2>", html, fixed = TRUE)
  html
}

build_topics <- function(config) {
  pkg <- config$package
  rd_files <- sort(list.files("man", pattern = "[.]Rd$", full.names = TRUE))
  out_dir <- file.path("docs", "topics"); dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  pandoc <- rmarkdown::pandoc_exec()
  # --markdown-headings needs pandoc 2.11.2 or later. Stop early with a clear message.
  if (rmarkdown::pandoc_version() < "2.11.2")
    stop("pandoc ", as.character(rmarkdown::pandoc_version()), " is too old. Use 2.11.2 or later ",
         "(set RSTUDIO_PANDOC to the folder of a newer pandoc).")

  # alias -> page stem, over ALL Rd files (a link may point at any alias)
  alias_map <- list()
  for (f in rd_files) {
    a <- regmatches(readLines(f, warn = FALSE), gregexpr("\\\\alias\\{[^}]*\\}", readLines(f, warn = FALSE)))
    for (al in gsub("^\\\\alias\\{|\\}$", "", unlist(a))) alias_map[[al]] <- tools::file_path_sans_ext(basename(f))
  }
  warnings_seen <- character(0)
  for (f in rd_files) {
    stem <- tools::file_path_sans_ext(basename(f))
    raw <- paste(readLines(f, warn = FALSE), collapse = "\n")
    if (!isTRUE(config$include_internal) && grepl("\\\\keyword\\{internal\\}", raw)) next
    if (!isFALSE(config$extras$link_rewrite))
      raw <- rewrite_rd_links(raw, pkg, alias_map, config$external_link_targets,
                              function(m) warnings_seen <<- c(warnings_seen, paste0(stem, ": ", m)))
    tmp_rd <- tempfile(fileext = ".Rd"); tmp_html <- tempfile(fileext = ".html"); tmp_md <- tempfile(fileext = ".md")
    writeLines(raw, tmp_rd)
    tools::Rd2HTML(tmp_rd, out = tmp_html, package = pkg)
    writeLines(tidy_html(readLines(tmp_html, warn = FALSE)), tmp_html)
    status <- system2(pandoc, c("-f", "html", "-t", "markdown-inline_code_attributes-fenced_divs", "--wrap=none", "--markdown-headings=atx",
                                shQuote(tmp_html), "-o", shQuote(tmp_md)))
    if (status != 0) stop("pandoc failed for ", f)
    body <- readLines(tmp_md, warn = FALSE)
    name <- sub("^.*\\\\name\\{([^}]*)\\}.*$", "\\1", grep("\\\\name\\{", readLines(f, warn = FALSE), value = TRUE)[1])
    writeLines(c(paste0("# ", name), "", body), file.path(out_dir, paste0(stem, ".md")))
  }
  # remove pages whose Rd is gone, so stale topics do not linger
  keep <- paste0(tools::file_path_sans_ext(basename(rd_files)), ".md")
  stale <- setdiff(list.files(out_dir, pattern = "[.]md$"), keep)
  if (length(stale)) { message("removing stale pages: ", paste(stale, collapse = ", ")); file.remove(file.path(out_dir, stale)) }
  if (length(warnings_seen)) message("Link warnings (", length(warnings_seen), "):\n  ", paste(unique(warnings_seen), collapse = "\n  "))
  invisible(NULL)
}
