# NCPM (library-size) scaling on sparse
#' @export
scaleNcpm <- function(M, per = 1e6) {
  M <- Matrix::forceSymmetric(M, uplo = "U")
  s <- sum(M@x); if (s > 0) M@x <- M@x / (s / per)
  M
}

# Observed/Expected by band distance on a DENSE matrix (upper & lower together)
#' @export
obsExpDense <- function(A) {
  n <- nrow(A)
  for (k in 1:(n - 1L)) {
    idx_i <- 1:(n - k); idx_j <- idx_i + k
    v <- A[cbind(idx_i, idx_j)]
    v <- v[is.finite(v) & v > 0]
    if (!length(v)) next
    expk <- stats::median(v)     # robust typical value for that distance
    if (!is.finite(expk) || expk <= 0) next
    A[cbind(idx_i, idx_j)] <- A[cbind(idx_i, idx_j)] / expk
    A[cbind(idx_j, idx_i)] <- A[cbind(idx_j, idx_i)] / expk
  }
  A
}

# 1D Gaussian kernel in *bins*
#' @export
gaussianKernel <- function(sigma_bins, radius = ceiling(3 * sigma_bins)) {
  if (!is.finite(sigma_bins) || sigma_bins <= 0) return(1)
  x <- -radius:radius
  k <- exp(-0.5 * (x / sigma_bins)^2)
  k / sum(k)
}

# separable Gaussian blur with simple reflect padding
#' @export
gaussBlurDense <- function(mat, krow) {
  if (length(krow) == 1L) return(mat)
  r <- (length(krow) - 1L) %/% 2L
  nr <- nrow(mat); nc <- ncol(mat)
  
  pad_reflect <- function(v, r) {
    if (r == 0) return(v)
    c(v[pmin.int(seq(r, 1), length(v))], v, v[pmax.int(length(v) - seq(0, r - 1), 1L)])
  }
  
  # row pass
  B <- matrix(0, nr, nc)
  for (i in seq_len(nr)) {
    pv <- pad_reflect(mat[i, ], r)
    conv <- stats::filter(pv, krow, sides = 2, circular = FALSE)
    B[i, ] <- conv[(r + 1):(r + nc)]
  }
  # col pass
  C <- matrix(0, nr, nc)
  for (j in seq_len(nc)) {
    pv <- pad_reflect(B[, j], r)
    conv <- stats::filter(pv, krow, sides = 2, circular = FALSE)
    C[, j] <- conv[(r + 1):(r + nr)]
  }
  C
}

#' @export
binLabelsFromCoords <- function(
    n_bins, 
    region_start, 
    bin_bp, 
    coords, 
    annotate_on = "name"
) {
  stopifnot(all(c("start","end", annotate_on) %in% names(coords)))
  region_start <- as.integer(region_start); bin_bp <- as.integer(bin_bp)
  
  bins <- data.table::data.table(
    bin   = 0:(n_bins-1L),
    start = region_start + (0:(n_bins-1L))*bin_bp,
    end   = region_start + (1:n_bins)*bin_bp - 1L
  )
  data.table::setkey(bins, start, end)
  
  cd <- data.table::as.data.table(coords)[
    , .(start = as.integer(start),
        end   = as.integer(end),
        label = as.character(get(annotate_on)))]
  cd[, `:=`(start = pmin(start,end), end = pmax(start,end))]
  data.table::setkey(cd, start, end)
  
  lab <- rep("intergenic", n_bins)
  ov  <- data.table::foverlaps(bins, cd, nomatch = 0L)
  if (nrow(ov)) {
    ov[, ovlen := pmin(end, i.end) - pmax(start, i.start) + 1L]
    best <- ov[ovlen > 0, .SD[which.max(ovlen)], by = bin]
    if (nrow(best)) lab[best$bin + 1L] <- best$label
  }
  lab
}

#' @export
geneContactMatrixSparse <- function(
    M,                      # dsCMatrix/dgCMatrix or base matrix (n_bins x n_bins)
    lab,                    # character vector length n_bins (one label per bin)
    drop_labels = "intergenic",  # bins with these labels are excluded
    weights = NULL          # optional per-bin weights (same length as lab)
) {
  if (!inherits(M, "Matrix")) M <- Matrix::Matrix(M, sparse = TRUE)
  n <- ncol(M); stopifnot(length(lab) == n)
  keep <- rep(TRUE, n)
  if (!is.null(drop_labels) && length(drop_labels))
    keep <- !(lab %in% drop_labels)
  
  # keep genomic order (first appearance)
  f <- factor(lab[keep], levels = unique(lab[keep]))
  ng <- nlevels(f); if (ng == 0) stop("No gene labels after filtering.")
  
  # sparse incidence B (n_bins x n_genes), one 1 per kept bin (or weight)
  B <- Matrix::sparseMatrix(
    i = which(keep),
    j = as.integer(f),
    x = if (is.null(weights)) rep(1, sum(keep)) else weights[keep],
    dims = c(n, ng),
    dimnames = list(NULL, levels(f))
  )
  
  # aggregate: G = t(B) %*% M %*% B
  # (works well even when M is triangular/symmetric sparse)
  G <- Matrix::t(B) %*% (M %*% B)
  G <- Matrix::forceSymmetric(G, uplo = "U")
  dimnames(G) <- list(levels(f), levels(f))
  G
}

#' @export
normalizeGeneContacts <- function(G, method = c("none","ncpm","oe"), scale = 1e6) {
  method <- match.arg(method)
  if (method == "none") return(G)
  if (method == "ncpm") {
    tot <- sum(G)
    if (tot > 0) return((G / tot) * scale) else return(G)
  }
  # observed/expected by marginals
  r <- Matrix::rowSums(G); T <- sum(r)
  E <- (r %o% r) / max(T, 1)            # symmetric expected
  out <- as.matrix(G) / E
  out[!is.finite(out)] <- 0
  Matrix::Matrix(out, sparse = FALSE)
}
#' @export
blurGeneMatrix <- function(G, sigma_genes = 0) {
  if (!is.finite(sigma_genes) || sigma_genes <= 0) return(G)
  k <- .gauss_k(sigma_genes)
  X <- if (inherits(G,"Matrix")) as.matrix(G) else G
  X <- t(apply(X, 1, .conv1d_reflect, k = k))
  X <- apply(X, 2, .conv1d_reflect, k = k)
  dim(X) <- dim(G); dimnames(X) <- dimnames(G)
  Matrix::Matrix(X, sparse = FALSE)
}

#' @export
plotGeneContactHeatmap <- function(
    G, palette = NULL, cap_q = 0.995, gamma = 1,
    legend_title = "signal", cluster = FALSE
) {
  X <- if (inherits(G,"Matrix")) as.matrix(G) else G
  # transform
  X <- log1p(pmax(X, 0))
  if (!is.null(cap_q)) {
    vals <- X[is.finite(X) & X > 0]
    if (length(vals)) {
      qq <- stats::quantile(vals, cap_q, names = FALSE)
      if (is.finite(qq) && qq > 0) X[X > qq] <- qq
    }
  }
  if (is.finite(gamma) && gamma != 1) {
    mx <- max(X, na.rm = TRUE); if (mx > 0) X <- (X/mx)^gamma * mx
  }
  if (cluster) {
    ord <- stats::hclust(stats::dist(X + t(X)))$order
    X <- X[ord, ord, drop = FALSE]
  }
  
  rn <- if (!is.null(rownames(X))) rownames(X) else as.character(seq_len(nrow(X)))
  cn <- if (!is.null(colnames(X))) colnames(X) else as.character(seq_len(ncol(X)))
  df <- data.frame(
    gene_i = factor(rep(rn, times = ncol(X)), levels = rn),
    gene_j = factor(rep(cn, each  = nrow(X)), levels = cn),
    value  = as.vector(X),
    stringsAsFactors = FALSE
  )

  p <- ggplot2::ggplot(df, ggplot2::aes(gene_j, gene_i, fill = value)) +
    ggplot2::geom_tile() +
    ggplot2::coord_fixed() +
    ggplot2::scale_fill_gradient(low = "white", high = "#B2182B", name = legend_title) +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(
      axis.title = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5),
      panel.grid = ggplot2::element_blank(),
      panel.background = ggplot2::element_rect(fill = "white", colour = NA),
      plot.background  = ggplot2::element_rect(fill = "white", colour = NA)
    )
  
  # color axis tick labels with your palette (if provided)
  if (!is.null(palette) && requireNamespace("ggtext", quietly = TRUE)) {
    labs_x <- function(x) sprintf("<span style='color:%s'>%s</span>", palette[x], x)
    labs_y <- function(x) sprintf("<span style='color:%s'>%s</span>", palette[x], x)
    p <- p +
      ggplot2::scale_x_discrete(labels = labs_x) +
      ggplot2::scale_y_discrete(labels = labs_y) +
      ggplot2::theme(axis.text.x = ggtext::element_markdown(size = 9),
                     axis.text.y = ggtext::element_markdown(size = 9))
  }
  p
}

#' @export
locusBlocksDf <- function(
    coords_db, 
    region_start, 
    bin_bp, 
    n_bins,
    annotate_on = "name",
    col_start = "start", 
    col_end = "end",
    collapse = TRUE
) {
  stopifnot(all(c(col_start, col_end, annotate_on) %in% names(coords_db)))
  rs <- as.integer(region_start); bb <- as.integer(bin_bp); nb <- as.integer(n_bins)
  reg_end <- rs + nb*bb - 1L
  
  cd <- data.table::as.data.table(coords_db)[
    , .(start = as.integer(get(col_start)),
        end   = as.integer(get(col_end)),
        label = as.character(get(annotate_on)))]
  # sanitize, clip to region
  cd[, `:=`(start = pmin(start, end), end = pmax(start, end))]
  cd <- cd[end >= rs & start <= reg_end]
  if (!nrow(cd)) return(data.table::data.table(label=character(), xmin=integer(), xmax=integer()))
  
  # convert bp -> bins (1-based, inclusive)
  cd[, xmin := ((pmax(start, rs)    - rs) %/% bb) + 1L]
  cd[, xmax := ((pmin(end,   reg_end) - rs) %/% bb) + 1L]
  cd[, `:=`(
    xmin = pmax(1L, pmin(nb, as.integer(xmin))),
    xmax = pmax(1L, pmin(nb, as.integer(xmax)))
  )]
  cd <- cd[xmin <= xmax]
  data.table::setorder(cd, label, xmin, xmax)
  
  if (!collapse) return(cd[, .(label, xmin, xmax)])
  
  # INTEGER-SAFE collapse of adjacent/overlapping ranges per label
  # treat [..xmax] and next [xmin..] as same group if gap <= 1 bin
  cd[, grp := cumsum(
    data.table::fifelse(xmin > data.table::shift(xmax, fill = 0L) + 1L, 1L, 0L)
  ), by = label]
  cd[, .(xmin = min(xmin), xmax = max(xmax)), by = .(label, grp)][, grp := NULL][]
}

#' @export
prepHeatmapMatrix <- function(
    M,
    bin_bp,
    norm   = c("ncpm", "oe", "none"),
    sigma_bp = 0,
    cap_q = 0.995,
    gamma = 1,
    upper_only = TRUE
) {
  norm <- match.arg(norm)
  
  # Coerce to symmetric sparse or dense
  Msym <-
    if (inherits(M, "Matrix")) Matrix::forceSymmetric(M, uplo = "U")
  else if (is.matrix(M))       if (isSymmetric(M)) M else (M + t(M))/2
  else stop("M must be a Matrix or base matrix")
  
  # Keep dimnames
  rn <- if (inherits(Msym, "Matrix")) rownames(Msym) else rownames(M)
  cn <- if (inherits(Msym, "Matrix")) colnames(Msym) else colnames(M)
  
  # Normalization before densify
  if (inherits(Msym, "Matrix") && norm == "ncpm") Msym <- scaleNcpm(Msym)
  
  # Densify + O/E if requested
  A <- if (inherits(Msym, "Matrix")) Matrix::as.matrix(Msym) else as.matrix(Msym)
  nr <- nrow(A); nc <- ncol(A); dim(A) <- c(nr, nc)
  
  if (norm == "oe") {
    A <- obsExpDense(A)           # ensures matrix + dims
  }
  
  # KDE smoothing
  if (is.finite(sigma_bp) && sigma_bp > 0) {
    sigma_bins <- sigma_bp / bin_bp
    krow <- gaussianKernel(sigma_bins)
    A <- gaussBlurDense(A, krow)  # ensures matrix + dims
  } else {
    dim(A) <- c(nr, nc)             # reassert
  }
  
  # log1p, cap, gamma (preserve dims)
  A <- log1p(pmax(A, 0)); dim(A) <- c(nr, nc)
  if (!is.null(cap_q)) {
    vals <- A[is.finite(A) & A > 0]
    if (length(vals)) {
      qq <- stats::quantile(vals, cap_q, names = FALSE)
      if (is.finite(qq) && qq > 0) A[A > qq] <- qq
    }
  }
  if (!is.null(gamma) && is.finite(gamma) && gamma != 1) {
    mx <- max(A[is.finite(A)], 0)
    if (mx > 0) { A <- (A / mx) ^ gamma * mx; dim(A) <- c(nr, nc) }
  }
  
  # mask lower triangle if desired
  if (upper_only) {
    lt <- lower.tri(A)
    A[lt] <- NA_real_
    dim(A) <- c(nr, nc)
  }
  
  # reattach dimnames if present
  if (!is.null(rn)) rownames(A) <- rn
  if (!is.null(cn)) colnames(A) <- cn
  
  A  # guaranteed matrix
}

#' @export
plotSplicingMatrix<-function(
    large_matrix,
    coords_db,
    annotate_on = "name",
    bin_bp      = 1000L,
    region_start = NULL,
    tick_bp     = 1e6,
    show_ticks  = c("none","mb"),
    upper_only  = TRUE,
    title       = "Splice map",
    norm = "oe",
    sigma_bp = 0,
    cap_q = NULL,
    gamma = 1
) {
  stopifnot(all(c("start","end", annotate_on) %in% names(coords_db)))
  show_ticks <- match.arg(show_ticks)
  
  # --- symmetrize, log1p, (optionally) mask lower triangle ---
  mat<-prepHeatmapMatrix(
    M = large_matrix,
    bin_bp = bin_bp, 
    upper_only = upper_only, 
    norm = norm, 
    sigma_bp = sigma_bp, 
    cap_q = cap_q, 
    gamma = gamma
    )
  n <-nrow(mat)
  # derive region so bins line up with matrix
  if (is.null(region_start)) region_start <- as.integer(min(coords_db$start, na.rm = TRUE))
  region_start   <- as.integer(region_start)
  region_end_eff <- region_start + n * as.integer(bin_bp) - 1L
  
  # ticks
  if (show_ticks == "mb") {
    ticks_pos <- seq(ceiling(region_start / tick_bp) * tick_bp, region_end_eff, by = tick_bp)
    tick_bins <- ((ticks_pos - region_start) %/% bin_bp) + 1L
    tick_labs <- paste0(round(ticks_pos / 1e6, 1), " Mb")
    axis_labs <- rep("", n); axis_labs[pmax(1L, pmin(n, tick_bins))] <- tick_labs
  } else {
    axis_labs <- rep("", n)
  }
  
  # bins (valid intervals)
  bins <- data.table::data.table(
    bin   = 0:(n - 1L),
    start = region_start + (0:(n - 1L)) * bin_bp,
    end   = region_start + (1:n) * bin_bp - 1L
  )
  bins[end > region_end_eff, end := region_end_eff]
  bins <- bins[start <= end]
  data.table::setkey(bins, start, end)
  
  # coords -> label per bin
  cd <- data.table::as.data.table(coords_db)[
    , .(start = as.integer(start),
        end   = as.integer(end),
        label = as.character(get(annotate_on)))]
  cd[, `:=`(start = pmin(start, end), end = pmax(start, end))]
  cd <- cd[end >= region_start & start <= region_end_eff]
  if (nrow(cd)) data.table::setkey(cd, start, end)
  
  ann_by_bin <- rep("intergenic", n)
  if (nrow(cd)) {
    ov <- data.table::foverlaps(bins, cd, nomatch = 0L)
    if (nrow(ov)) {
      ov[, ovlen := pmin(end, i.end) - pmax(start, i.start) + 1L]
      best <- ov[ovlen > 0, .SD[which.max(ovlen)], by = bin]
      if (nrow(best)) ann_by_bin[best$bin + 1L] <- best$label
    }
  }
  
  ids <- paste0("b", seq_len(n))
  rownames(mat) <- colnames(mat) <- ids
  ann_df <- data.frame(label = ann_by_bin, row.names = ids)
  lev    <- unique(ann_df$label)
  pal    <- setNames(hcl.colors(length(lev), "Spectral"), lev)
  
  pheatmap::pheatmap(
    mat,
    cluster_rows = FALSE, cluster_cols = FALSE,
    labels_row = axis_labs, labels_col = axis_labs,
    annotation_row = ann_df, annotation_col = ann_df,
    annotation_colors = list(label = pal),
    na_col = "white",           # blanks for lower triangle
    border_color = NA,
    main = sprintf("%s (bin = %s bp)", title, bin_bp)
  )
}

#' @export
heatmapWr <- function(
    H, 
    lab, 
    title = "Splice map"
) {
  stopifnot(is.matrix(H) || inherits(H, "Matrix"))
  n <- ncol(H); stopifnot(length(lab) == n)
  
  vmax <- max(H[is.finite(H)], na.rm = TRUE)
  
  pal_lab <- setNames(hcl.colors(length(unique(lab)), "Spectral"),
                      unique(lab))
  
  ggheatmap(H, filling = "raster") +                      # alias of heatmap_layout
    ggplot2::scale_fill_gradient(low = "white", high = "#B2182B",  # 1-D scale
                        limits = c(0, vmax),
                        oob = scales::squish,
                        na.value = "white",
                        name = "signal") +
    anno_top(size = 6) +                                  # add a top annotation band
    ggalign(                                              # draw a colored strip from `lab`
      data = data.frame(label = lab)
    ) +
    ggplot2::geom_tile(ggplot2::aes(x = .x, y = 1, fill = label),           # .x comes from the layout
              height = 1) +
    ggplot2::scale_fill_manual(values = pal_lab, guide = "none") + # use same colors for row/col bars
    ggplot2::ggtitle(title) &
    ggplot2::theme_minimal(base_size = 9) &
    ggplot2::theme(axis.text = ggplot2::element_blank(),
          panel.grid = ggplot2::element_blank(),
          plot.margin = ggplot2::margin(2, 2, 2, 2))
}

#' @export
plotSpliceArcs <- function(
    large_matrix,
    k_top          = 200,
    min_dist_bp    = 0,
    bin_bp         = 1000L,
    title          = "Top arcs"
) {
  # prefer raw/symmetric counts for ranking
  Msym <- Matrix::forceSymmetric(large_matrix, uplo = "U")
  S <- Matrix::summary(Matrix::triu(Msym))  # i, j, x (1-based)
  S <- S[S$i < S$j & is.finite(S$x), , drop = FALSE]
  if (min_dist_bp > 0) {
    min_bins <- ceiling(min_dist_bp / bin_bp)
    S <- S[(S$j - S$i) >= min_bins, , drop = FALSE]
  }
  if (!nrow(S)) {
    return(ggplot2::ggplot() + ggplot2::theme_void() + ggplot2::ggtitle("Top arcs (none)"))
  }
  S <- S[order(S$x, decreasing = TRUE), , drop = FALSE]
  if (nrow(S) > k_top) S <- S[seq_len(k_top), , drop = FALSE]
  n <- max(S$j, na.rm = TRUE)
  
  ggplot2::ggplot(as.data.frame(S)) +
    ggplot2::geom_curve(
      ggplot2::aes(x = i, xend = j, y = 0, yend = 0,
                   alpha = x, linewidth = x),
      curvature = 0.25, lineend = "round"
    ) +
    ggplot2::scale_x_continuous(limits = c(1, n), expand = c(0, 0)) +
    ggplot2::scale_alpha_continuous(range = c(0.2, 0.8), guide = "none") +
    ggplot2::scale_linewidth_continuous(range = c(0.2, 1.6), guide = "none") +
    ggplot2::theme_void() +
    ggplot2::ggtitle(title)
}

#' @export
arcsFromMatrix <- function(
    H,
    k = 200,
    min_dist_bins = 1L,
    y = 1
) {
  Hsym <- Matrix::forceSymmetric(H, uplo = "U")
  S <- Matrix::summary(Matrix::triu(Hsym))
  S <- S[S$i < S$j & (S$j - S$i) >= min_dist_bins & is.finite(S$x), , drop = FALSE]
  if (!nrow(S)) return(data.frame())
  S <- S[order(S$x, decreasing = TRUE), , drop = FALSE]
  if (nrow(S) > k) S <- S[seq_len(k), , drop = FALSE]

  data.frame(
    x    = S$i,
    xend = S$j,
    y    = y,
    yend = y,
    w    = scales::rescale(S$x, to = c(0.2, 0.95))
  )
}

#' @export
arcsPanel <- function(
    H, 
    k = 200, 
    min_dist_bins = 5L
) {
  n <- ncol(H)
  S <- arcsFromMatrix(H, k = k, min_dist_bins = min_dist_bins)
  if (!nrow(S)) {
    return(ggplot2::ggplot() + 
             ggplot2::xlim(1, n) + ggplot2::ylim(0, 1) +
             ggplot2::theme_void() + ggplot2::ggtitle("Top arcs (none)"))
  }
  ggplot2::ggplot(S) +
    ggplot2::geom_curve(ggplot2::aes(x = x, xend = xend, y = 0, yend = 0,
                   alpha = ..curv.., size = ..curv..),
               curvature = 0.25, lineend = "round") +
    ggplot2::scale_x_continuous(limits = c(1, n), expand = c(0, 0)) +
    ggplot2::scale_alpha_continuous(range = c(0.15, 0.8), guide = "none") +
    ggplot2::scale_size_continuous(range = c(0.2, 1.5),  guide = "none") +
    ggplot2::theme_void() +
    ggplot2::ggtitle("Top arcs")
}

#' @export
plotLocusTrackBins <- function(
    blocks, 
    n_bins, 
    pal = NULL, 
    show_labels = TRUE, 
    label_min_bins = 3L,
    label_size = 3, 
    label_color = "black", 
    title = "Locus"
) {
  if (!nrow(blocks)) return(ggplot2::ggplot() + ggplot2::theme_void() + ggplot2::ggtitle(paste0(title, " (none)")))
  labs <- unique(blocks$label)
  if (is.null(pal)) pal <- setNames(hcl.colors(length(labs), "Spectral"), labs)
  
  p <- ggplot2::ggplot(blocks) +
    ggplot2::geom_rect(ggplot2::aes(xmin = xmin, xmax = xmax + 1L, ymin = 0, ymax = 1, fill = label),
              color = NA, alpha = 0.6) +
    ggplot2::scale_fill_manual(values = pal, guide = "none") +
    ggplot2::scale_x_continuous(limits = c(1, n_bins), expand = c(0, 0)) +
    ggplot2::coord_cartesian(ylim = c(0, 1)) +
    ggplot2::theme_void() + ggplot2::ggtitle(title)
  
  if (show_labels) {
    lbl <- blocks[(xmax - xmin + 1L) >= as.integer(label_min_bins),
                  .(x = (xmin + xmax + 1)/2, label = label)]
    if (nrow(lbl)) {
      p <- p + ggplot2::geom_text(data = lbl, ggplot2::aes(x = x, y = 0.5, label = label),
                         size = label_size, color = label_color, check_overlap = TRUE)
    }
  }
  p
}

#' @export
plotSpliceArcsColored <- function(
    large_matrix,
    k_top        = 200,
    min_dist_bp  = 0,
    bin_bp       = 1000L,
    title        = "Top arcs",
    lab          = NULL,                # character vector length = ncol(M); bin -> gene label
    coords_db    = NULL,                # if lab is NULL, build from coords_db + region_start
    annotate_on  = "name",
    region_start = NULL,
    color_by     = c("pair", "left", "right", "same-diff"),
    pal_lab      = NULL,                # named vector gene -> color (optional, for left/right)
    curvature    = 0.25,
    show_legend  = TRUE
) {
  color_by <- match.arg(color_by)
  
  # arcs from raw/symmetric counts
  Msym <- Matrix::forceSymmetric(large_matrix, uplo = "U")
  S <- Matrix::summary(Matrix::triu(Msym))  # i, j, x  (1-based)
  S <- S[S$i < S$j & is.finite(S$x), , drop = FALSE]
  if (min_dist_bp > 0) {
    min_bins <- ceiling(min_dist_bp / bin_bp)
    S <- S[(S$j - S$i) >= min_bins, , drop = FALSE]
  }
  if (!nrow(S)) {
    return(ggplot2::ggplot() + ggplot2::theme_void() + ggplot2::ggtitle("Top arcs (none)"))
  }
  S <- S[order(S$x, decreasing = TRUE), , drop = FALSE]
  if (nrow(S) > k_top) S <- S[seq_len(k_top), , drop = FALSE]
  S <- as.data.frame(S)
  n <- max(S$j, na.rm = TRUE)
  
  # --- get bin -> gene labels ---
  if (is.null(lab)) {
    if (is.null(coords_db) || is.null(region_start))
      stop("Provide `lab` or (`coords_db` + `region_start`) to color arcs by genes.")
    lab <- binLabelsFromCoords(n_bins = max(n, ncol(Msym)),
                                  region_start = region_start,
                                  bin_bp = bin_bp,
                                  coords = coords_db,
                                  annotate_on = annotate_on)
  }
  if (length(lab) < n) stop("`lab` must have length >= number of bins in matrix.")
  
  S$g_i <- lab[S$i]
  S$g_j <- lab[S$j]
  S$pair <- ifelse(S$g_i <= S$g_j,
                   paste0(S$g_i, "<->", S$g_j),
                   paste0(S$g_j, "<->", S$g_i))
  S$class <- ifelse(S$g_i == S$g_j, "intra", "inter")
  
  # --- aesthetics: weight by count; precompute alpha/width ---
  S$w <- scales::rescale(log1p(S$x), to = c(0.2, 1.0))  # visual weight
  
  # --- choose color mapping ---
  if (color_by == "pair") {
    colvar <- "pair"
    pairs <- unique(S$pair)
    pal <- setNames(hcl.colors(length(pairs), "Spectral"), pairs)
    sc_color <- ggplot2::scale_color_manual(values = pal, guide = if (show_legend) "legend" else "none")
    legend_title <- "Gene pair"
  } else if (color_by == "left") {
    colvar <- "g_i"
    genes <- unique(S$g_i)
    if (is.null(pal_lab)) pal_lab <- setNames(hcl.colors(length(genes), "Spectral"), genes)
    pal <- pal_lab
    sc_color <- ggplot2::scale_color_manual(values = pal, guide = if (show_legend) "legend" else "none")
    legend_title <- "Left gene"
  } else if (color_by == "right") {
    colvar <- "g_j"
    genes <- unique(S$g_j)
    if (is.null(pal_lab)) pal_lab <- setNames(hcl.colors(length(genes), "Spectral"), genes)
    pal <- pal_lab
    sc_color <- ggplot2::scale_color_manual(values = pal, guide = if (show_legend) "legend" else "none")
    legend_title <- "Right gene"
  } else { # same-diff
    colvar <- "class"
    pal <- c(intra = "#1b9e77", inter = "#d95f02")
    sc_color <- ggplot2::scale_color_manual(values = pal, guide = if (show_legend) "legend" else "none")
    legend_title <- "Class"
  }
  
  ggplot2::ggplot(S) +
    ggplot2::geom_curve(
      ggplot2::aes(x = i, xend = j, y = 0, yend = 0,
          color = .data[[colvar]],
          alpha = w, linewidth = w),
      curvature = curvature, lineend = "round"
    ) +
    sc_color + ggplot2::guides(color = ggplot2::guide_legend(title = legend_title, override.aes = list(alpha = 1, linewidth = 1))) +
    ggplot2::scale_alpha_identity() +
    ggplot2::scale_linewidth_identity() +
    ggplot2::scale_x_continuous(limits = c(1, n), expand = c(0, 0)) +
    ggplot2::theme_void() +
    ggplot2::ggtitle(title)
}

#' @export
composeHeatmapWithArcs <- function(
    heatmap_obj,
    arcs_plot,
    heights = c(1, 4)
) {
  ggalign::align_plots(arcs_plot, heatmap_obj, ncol = 1, heights = heights)
}

#' @export
composeLocusAndArcs <- function(
    arcs_plot, 
    coords_db, 
    region_start, 
    bin_bp,
    n_bins, 
    annotate_on = "name", 
    pal = NULL, 
    locus_height = 1, 
    arcs_height = 3,
    locus_title = "Locus", 
    arcs_title = NULL) 
{
  blocks <- locusBlocksDf(coords_db, region_start, bin_bp, n_bins,
                            annotate_on = annotate_on, collapse = TRUE)
  p_locus <- plotLocusTrackBins(blocks, n_bins, pal = pal, title = locus_title)
  
  # apply optional title override to arcs
  if (!is.null(arcs_title)) arcs_plot <- arcs_plot + ggplot2::ggtitle(arcs_title)
  
  # Use ggalign's align_plots if you're already in that ecosystem; otherwise patchwork
  # With base ggplot2, patchwork is simplest:
  suppressPackageStartupMessages(requireNamespace("patchwork"))
  patchwork::wrap_plots(p_locus, arcs_plot, ncol = 1, heights = c(locus_height, arcs_height))
}

#' @export
plotSplicingMatrixGgalign2 <- function(
    large_matrix,
    coords_db,
    annotate_on      = "name",
    bin_bp           = 1000L,
    region_start     = NULL,
    # display / normalization
    upper_only       = TRUE,
    norm             = c("ncpm", "oe", "none"),
    sigma_bp         = 0,
    cap_q            = 0.995,
    gamma            = 1,
    # readability helpers
    mask_diagonal_bp = 0,          # blank a diagonal band (half-width in bp; 0 = off)
    show_boundaries  = TRUE,       # thin gridlines at gene borders
    shade_by_gene    = FALSE,      # vertical background strips per gene
    shade_alpha      = 0.06,
    # diagonal ribbon
    diag_band_bp     = 0,          # diagonal ribbon half-width (bp; 0 = off)
    diag_band_alpha  = 0.25,
    diag_color_by    = c("row","col"),
    # palette & labels
    palette          = NULL,       # named vector label->color; NULL = auto Spectral
    label_min_bins   = 3L,         # only label top blocks with width >= this
    title            = "Splice map"
) {
  stopifnot(all(c("start","end", annotate_on) %in% names(coords_db)))
  norm <- match.arg(norm)
  diag_color_by <- match.arg(diag_color_by)
  
  # --- 1) prep display matrix ---
  H <- prepHeatmapMatrix(
    M = large_matrix, bin_bp = bin_bp, norm = norm,
    sigma_bp = sigma_bp, cap_q = cap_q, gamma = gamma, upper_only = upper_only
  )
  H <- as.matrix(H)
  n <- nrow(H); if (n == 0) stop("Empty matrix after prep.")
  
  na_df <- as.data.frame(which(is.na(H), arr.ind = TRUE))
  colnames(na_df) <- c("i","j")
  # triangle mask (white via NA)
  triangle<-"lower"
  if(upper_only) triangle<-"upper"
  if (triangle == "upper") {
    H[lower.tri(H, diag = FALSE)] <- NA_real_
  } else if (triangle == "lower") {
    H[upper.tri(H, diag = FALSE)] <- NA_real_   # <- upper white
  } # "full" shows both
  
  # optional near-diagonal mask
  if (is.numeric(mask_diagonal_bp) && mask_diagonal_bp > 0) {
    k <- ceiling(mask_diagonal_bp / bin_bp)
    ii <- row(H); jj <- col(H)
    H[abs(ii - jj) <= k] <- NA_real_
  }
  
  # --- 2) per-bin labels from coords ---
  if (is.null(region_start)) region_start <- as.integer(min(coords_db$start, na.rm = TRUE))
  region_start   <- as.integer(region_start)
  region_end_eff <- region_start + n * as.integer(bin_bp) - 1L
  
  bins <- data.table::data.table(
    bin   = 0:(n - 1L),
    start = region_start + (0:(n - 1L)) * bin_bp,
    end   = pmin(region_start + (1:n) * bin_bp - 1L, region_end_eff)
  )
  data.table::setkey(bins, start, end)
  
  cd <- data.table::as.data.table(coords_db)[
    , .(start = as.integer(start),
        end   = as.integer(end),
        label = as.character(get(annotate_on)))]
  cd[, `:=`(start = pmin(start, end), end = pmax(start, end))]
  cd <- cd[end >= region_start & start <= region_end_eff]
  if (nrow(cd)) data.table::setkey(cd, start, end)
  
  lab <- rep("intergenic", n)
  if (nrow(cd)) {
    ov <- data.table::foverlaps(bins, cd, nomatch = 0L)
    if (nrow(ov)) {
      ov[, ovlen := pmin(end, i.end) - pmax(start, i.start) + 1L]
      best <- ov[ovlen > 0, .SD[which.max(ovlen)], by = bin]
      if (nrow(best)) lab[best$bin + 1L] <- best$label
    }
  }
  
  # contiguous blocks (for labels, shading, boundaries)
  r      <- rle(lab)
  ends   <- cumsum(r$lengths)
  starts <- c(1L, head(ends, -1L) + 1L)
  blocks <- data.table::data.table(label = r$values, xmin = starts, xmax = ends)
  
  # per-bin strip data (n rows)  <<-- key for ggalign top annotation
  strip_df <- data.frame(label = lab)
  
  # palette
  lev_all <- unique(blocks$label)
  pal <- if (is.null(palette)) {
    setNames(hcl.colors(length(lev_all), "Spectral"), lev_all)
  } else {
    missing <- setdiff(lev_all, names(palette))
    c(palette, setNames(hcl.colors(length(missing), "Spectral"), missing))
  }
  
  vmax <- max(H[is.finite(H)], na.rm = TRUE); if (!is.finite(vmax) || vmax <= 0) vmax <- 1
  
  # --- 3) base heatmap ---
  p <- ggalign::ggheatmap(H, filling = "raster") +
    ggplot2::scale_fill_gradient(
      low = "white", high = "#B2182B",
      limits = c(0, vmax), oob = scales::squish,
      na.value = "white", name = "signal"
    )
  
  # boundaries
  if (isTRUE(show_boundaries)) {
    b_edges <- which(lab[-length(lab)] != lab[-1]) + 0.5
    if (length(b_edges)) {
      p <- p +
        ggplot2::geom_vline(xintercept = b_edges, color = "grey85", linewidth = 0.2) +
        ggplot2::geom_hline(yintercept = b_edges, color = "grey85", linewidth = 0.2)
    }
  }
  
  # optional vertical shading by gene (main panel) - separate fill scale
  shade_only_non_intergenic<-TRUE
  if (isTRUE(shade_by_gene) && nrow(blocks)) {
    bg_blocks <- if (shade_only_non_intergenic) blocks[label != "intergenic"] else blocks
    p <- tryCatch(
      p + ggnewscale::new_scale_fill() +
        ggplot2::geom_rect(
          data = bg_blocks,
          ggplot2::aes(xmin = xmin, xmax = xmax + 1, ymin = 0.5, ymax = n + 0.5, fill = label),
          inherit.aes = FALSE, alpha = shade_alpha, color = NA
        ) +
        ggplot2::scale_fill_manual(values = pal, guide = "none"),
      error = function(e) {
        warning("plotSplicingMatrixGgalign2: ggnewscale incompatible with the loaded ggalign; skipping gene-shading overlay.", call. = FALSE)
        p
      }
    )
  }

  # optional gene-colored diagonal ribbon (main panel) - separate fill scale
  if (is.numeric(diag_band_bp) && diag_band_bp > 0) {
    k <- max(1L, ceiling(diag_band_bp / bin_bp))
    idx <- seq_len(n)
    overlay <- data.table::CJ(i = idx, d = -k:k)[, j := i + d][j >= 1 & j <= n & abs(i - j) <= k]
    overlay[, label := if (diag_color_by == "row") lab[i] else lab[j]]
    p <- tryCatch(
      p + ggnewscale::new_scale_fill() +
        ggplot2::geom_tile(
          data = overlay,
          ggplot2::aes(x = j, y = i, fill = label),
          inherit.aes = FALSE, alpha = diag_band_alpha, width = 1, height = 1
        ) +
        ggplot2::scale_fill_manual(values = pal, guide = "none"),
      error = function(e) {
        warning("plotSplicingMatrixGgalign2: ggnewscale incompatible with the loaded ggalign; skipping diagonal-ribbon overlay.", call. = FALSE)
        p
      }
    )
  }
  
  if (triangle != "full") {
    if (triangle == "upper") {
      # cover the upper triangle (y > x)
      mask_poly <- data.frame(
        x = c(0.5, n + 0.5, 0.5),
        y = c(n + 0.5, n + 0.5, 0.5)
      )
    } else { # triangle == "upper"
      # cover the lower triangle (y < x)
      mask_poly <- data.frame(
        x = c(0.5, n + 0.5, n + 0.5),
        y = c(0.5, 0.5, n + 0.5)
      )
    }
    p <- p + ggplot2::geom_polygon(
      data = mask_poly,
      ggplot2::aes(x = x, y = y),
      inherit.aes = FALSE, fill = "white", color = NA
    )
  }
  # --- 4) top gene strip (per-bin; matches nobs) + labels from blocks ---
  # ggalign 1.x gives each annotation panel its own scale; the legacy
  # ggnewscale::new_scale_fill() reset is unnecessary (and incompatible with
  # ggalign's S7 layout objects). We try it for older ggalign versions, then
  # fall through if it errors.
  p <- tryCatch(p + ggnewscale::new_scale_fill(), error = function(e) p)
  p <- p +
    ggalign::anno_top(size = 0.16) +
    ggalign::ggalign(data = strip_df) +
    ggplot2::geom_tile(
      ggplot2::aes(x = .data$.x, y = 1, fill = label),
      height = 1, width = 1, show.legend = TRUE
    ) +
    ggplot2::scale_fill_manual(
      values = pal, name = annotate_on,
      guide = ggplot2::guide_legend(nrow = 1, byrow = TRUE, override.aes = list(alpha = 1))
    )
  if (requireNamespace("ggrepel", quietly = TRUE)) {
    label_exclude<-c("intergenic")
    lab_blocks <- blocks[(xmax - xmin + 1L) >= label_min_bins]
    if (!is.null(label_exclude) && length(label_exclude)) {
      lab_blocks <- lab_blocks[!(label %in% label_exclude)]
    }
    if (nrow(lab_blocks)) {
      lab_blocks[, x := (xmin + xmax + 1) / 2]
      p <- p +
        ggrepel::geom_text_repel(
          data = lab_blocks, ggplot2::aes(
            x = x,
            y = 1.05,
            label = label),
          inherit.aes = FALSE, direction = "x",
          box.padding = 0.2, point.padding = 0.1,
          min.segment.length = 0, max.overlaps = Inf,
          size = 3, color = "black", segment.color = "grey60"
        )
    }
  }

  # --- 5) cosmetics ---
  p <- p +
    ggplot2::ggtitle(sprintf("%s (bin = %s bp)", title, bin_bp)) +
    ggplot2::theme_minimal(base_size = 9) +
    ggplot2::theme(
      axis.text = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(2, 2, 2, 2),
      legend.position = "top"
    )
  p <- p +
    ggplot2::theme(
      panel.background = ggplot2::element_rect(fill = "white", colour = NA_real_),
      plot.background  = ggplot2::element_rect(fill = "white", colour = NA_real_)
    )
  
  invisible(list(
    plot      = p,
    labels    = lab,
    palette   = pal,
    blocks    = blocks,
    region    = c(start = region_start, end = region_end_eff),
    bin_bp    = bin_bp
  ))
  p
}

#' @export
plotGeneContactsGg <- function(
    large_matrix,                 # dsCMatrix/dgCMatrix or base matrix (bins x bins)
    coords_db,                    # data.frame/data.table with start, stop, and annotate_on
    annotate_on   = "name",
    bin_bp        = 1000L,
    region_start  = NULL,         # if NULL uses min(coords_db$start)
    # gene filtering
    drop_labels   = "intergenic", # set NULL to keep all
    # normalization at gene level
    norm          = c("ncpm","oe","none"),
    # small blur in *gene index* space (NOT bin space)
    sigma_genes   = 0,            # e.g. 0.6-1.0 genes
    # display shaping
    cap_q         = 0.995,
    gamma         = 1,
    cluster       = FALSE,        # cluster genes instead of genomic order
    palette       = NULL,         # named vector gene->color; NULL auto "Spectral"
    title         = "Gene x gene splice contacts"
) {
  stopifnot(all(c("start","end", annotate_on) %in% names(coords_db)))
  norm <- match.arg(norm)
  
  # --- 0) prepare matrix (symmetric sparse) ---
  M <- large_matrix
  if (!inherits(M, "Matrix")) M <- Matrix::Matrix(M, sparse = TRUE)
  M <- Matrix::forceSymmetric(M, uplo = "U")
  n_bins <- ncol(M)
  if (n_bins != nrow(M)) stop("Matrix must be square (bins x bins).")
  
  # --- 1) per-bin labels from coords (like your other fn) ---
  if (is.null(region_start)) region_start <- as.integer(min(coords_db$start, na.rm = TRUE))
  region_start <- as.integer(region_start)
  region_end_eff <- region_start + n_bins * as.integer(bin_bp) - 1L
  
  bins <- data.table::data.table(
    bin   = 0:(n_bins - 1L),
    start = region_start + (0:(n_bins - 1L)) * bin_bp,
    end   = pmin(region_start + (1:n_bins) * bin_bp - 1L, region_end_eff)
  )
  data.table::setkey(bins, start, end)
  
  cd <- data.table::as.data.table(coords_db)[
    , .(start = as.integer(start),
        end   = as.integer(end),
        label = as.character(get(annotate_on)))]
  cd[, `:=`(start = pmin(start, end), end = pmax(start, end))]
  cd <- cd[end >= region_start & start <= region_end_eff]
  if (nrow(cd)) data.table::setkey(cd, start, end)
  
  lab <- rep("intergenic", n_bins)
  if (nrow(cd)) {
    ov <- data.table::foverlaps(bins, cd, nomatch = 0L)
    if (nrow(ov)) {
      ov[, ovlen := pmin(end, i.end) - pmax(start, i.start) + 1L]
      best <- ov[ovlen > 0, .SD[which.max(ovlen)], by = bin]
      if (nrow(best)) lab[best$bin + 1L] <- best$label
    }
  }
  
  # --- 2) build sparse incidence B and aggregate: G = t(B) %*% M %*% B ---
  keep <- rep(TRUE, n_bins)
  if (!is.null(drop_labels) && length(drop_labels)) keep <- !(lab %in% drop_labels)
  if (!any(keep)) stop("No bins retained after drop_labels filter.")
  
  # maintain genomic order of first appearance
  f <- factor(lab[keep], levels = unique(lab[keep]))
  ng <- nlevels(f)
  B <- Matrix::sparseMatrix(
    i = which(keep),
    j = as.integer(f),
    x = 1,
    dims = c(n_bins, ng),
    dimnames = list(NULL, levels(f))
  )
  G <- Matrix::t(B) %*% (M %*% B)
  G <- Matrix::forceSymmetric(G, uplo = "U")
  dimnames(G) <- list(levels(f), levels(f))
  
  # --- 3) normalization at gene level ---
  if (norm == "ncpm") {
    tot <- sum(G)
    if (tot > 0) G <- (G / tot) * 1e6
  } else if (norm == "oe") {
    r <- Matrix::rowSums(G)
    Tt <- sum(r)
    E <- (r %o% r) / max(Tt, 1)          # expected under independence
    X <- as.matrix(G) / E
    X[!is.finite(X)] <- 0
    G <- Matrix::Matrix(X, sparse = FALSE)
  }
  
  # --- 4) optional small blur in gene-index space (Nebulosa-ish) ---
  if (is.finite(sigma_genes) && sigma_genes > 0) {
    k_gauss <- function(s) { r <- ceiling(3*s); x <- -r:r; w <- exp(-0.5*(x/s)^2); w/sum(w) }
    conv1d_reflect <- function(v, k) {
      m <- length(v); if (length(k) == 1L) return(v)
      r <- (length(k)-1L)%/%2L; out <- numeric(m)
      for (p in seq_len(m)) {
        acc <- 0.0
        for (t in -r:r) {
          idx <- p + t
          if (idx < 1L) idx <- 2L - idx else if (idx > m) idx <- 2L*m - idx
          acc <- acc + v[idx] * k[t + r + 1L]
        }
        out[p] <- acc
      }
      out
    }
    K <- k_gauss(sigma_genes)
    X <- if (inherits(G, "Matrix")) as.matrix(G) else G
    X <- t(apply(X, 1, conv1d_reflect, k = K))
    X <- apply(X, 2, conv1d_reflect, k = K)
    dim(X) <- dim(G); dimnames(X) <- dimnames(G)
    G <- Matrix::Matrix(X, sparse = FALSE)
  }
  
  # --- 5) transform for display (log1p + cap + gamma) ---
  X <- if (inherits(G, "Matrix")) as.matrix(G) else G
  X <- log1p(pmax(X, 0))
  if (!is.null(cap_q)) {
    vals <- X[is.finite(X) & X > 0]
    if (length(vals)) {
      qq <- stats::quantile(vals, cap_q, names = FALSE)
      if (is.finite(qq) && qq > 0) X[X > qq] <- qq
    }
  }
  if (is.finite(gamma) && gamma != 1) {
    mx <- max(X, na.rm = TRUE); if (mx > 0) X <- (X/mx)^gamma * mx
  }
  
  # --- 6) ordering & palette (ensure dimnames exist) ---
  if (is.null(rownames(G)) || is.null(colnames(G))) {
    # fall back to factor levels (genes_in_order) if available
    genes_in_order <- if (exists("f")) levels(f) else seq_len(nrow(G))
    dimnames(G) <- list(genes_in_order, genes_in_order)
  }
  X <- if (inherits(G, "Matrix")) as.matrix(G) else G
  
  genes_i <- rownames(X)
  genes_j <- colnames(X)
  if (length(genes_i) == 0L || length(genes_j) == 0L)
    stop("Gene x gene matrix has no dimnames; likely no genes remained after filtering (drop_labels?).")
  
  if (cluster) {
    ord <- stats::hclust(stats::dist(X + t(X)))$order
    X <- X[ord, ord, drop = FALSE]
    genes_i <- rownames(X); genes_j <- colnames(X)
  }
  
  pal <- if (is.null(palette)) {
    setNames(hcl.colors(length(genes_i), "Spectral"), genes_i)
  } else {
    missing <- setdiff(genes_i, names(palette))
    c(palette, setNames(hcl.colors(length(missing), "Spectral"), missing))[genes_i]
  }
  
  # --- 7) robust long-format consistent with column-major as.vector() ---
  df <- data.frame(
    gene_j = rep(genes_j, each = nrow(X)),   # column-major: j changes slowest
    gene_i = rep(genes_i, times = ncol(X)),  # i runs within each column
    value  = as.vector(X),
    stringsAsFactors = FALSE
  )
  # lock plotting order
  df$gene_i <- factor(df$gene_i, levels = genes_i)
  df$gene_j <- factor(df$gene_j, levels = genes_j)
  
  # --- 8) plot (white->red; colored tick labels if ggtext present) ---
  p <- ggplot2::ggplot(df, ggplot2::aes(x = gene_j, y = gene_i, fill = value)) +
    ggplot2::geom_tile() +
    ggplot2::coord_fixed() +
    ggplot2::scale_fill_gradient(low = "white", high = "#B2182B", name = "signal") +
    ggplot2::ggtitle(title) +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(
      axis.title = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5),
      panel.grid = ggplot2::element_blank(),
      panel.background = ggplot2::element_rect(fill = "white", colour = NA),
      plot.background  = ggplot2::element_rect(fill = "white", colour = NA),
      legend.position  = "right"
    )
  
  #if (requireNamespace("ggtext", quietly = TRUE)) {
  #  lab_x <- function(x) sprintf("<span style='color:%s'>%s</span>", pal[as.character(x)], x)
  #  lab_y <- function(x) sprintf("<span style='color:%s'>%s</span>", pal[as.character(x)], x)
  #  p <- p +
  #    ggplot2::scale_x_discrete(labels = lab_x) +
  #    ggplot2::scale_y_discrete(labels = lab_y) +
  #    ggplot2::theme(axis.text.x = ggtext::element_markdown(size = 9),
  #                   axis.text.y = ggtext::element_markdown(size = 9))
  #}
  p <- p + ggplot2::scale_x_discrete() + ggplot2::scale_y_discrete()
  
  invisible(list(
    plot       = p,
    G          = G,
    display    = X,
    genes      = genes_i,
    palette    = pal,
    bin_labels = lab,
    kept_bins  = which(keep),
    region     = c(start = region_start, end = region_end_eff),
    bin_bp     = bin_bp
  ))
  p
  
}

# tiny blur across the *gene index* (not bins); sigma in "genes"
.gauss_k <- function(s) { if (s<=0) return(1); r <- ceiling(3*s); x <- -r:r; k <- exp(-0.5*(x/s)^2); k/sum(k) }

.conv1d_reflect <- function(v, k) {
  m <- length(v); if (length(k)==1L) return(v)
  r <- (length(k)-1L)%/%2L; out <- numeric(m)
  for (p in seq_len(m)) {
    acc <- 0.0
    for (t in -r:r) {
      idx <- p + t
      if (idx < 1L) idx <- 2L - idx else if (idx > m) idx <- 2L*m - idx
      acc <- acc + v[idx] * k[t + r + 1L]
    }
    out[p] <- acc
  }
  out
}
