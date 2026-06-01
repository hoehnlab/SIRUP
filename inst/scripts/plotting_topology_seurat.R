plot_joint_density <- function(seurat_object, gene_1, gene_2, reduction,
                               show_contours = TRUE,
                               contour_dist_threshold = 100, pt.size = 1) {
  library(FNN)
  library(spdep)
  library(igraph)
  library(akima)
  library(patchwork)
  
  # --- Single Nebulosa call for all three plots ---
  int_all <- Nebulosa::plot_density(
    object = seurat_object, features = c(gene_1, gene_2),
    joint = TRUE, combine = FALSE,
    reduction = reduction, slot = 'data', size = pt.size)
  
  d1     <- int_all[[1]]$data
  d2     <- int_all[[2]]$data
  djoint <- int_all[[3]]$data
  
  setDT(d1,     keep.rownames = 'barcode')
  setDT(d2,     keep.rownames = 'barcode')
  setDT(djoint, keep.rownames = 'barcode')
  
  setnames(d1,     'feature', 'density_gene_1')
  setnames(d2,     'feature', 'density_gene_2')
  setnames(djoint, 'feature', 'joint_density')
  
  tdd <- merge(djoint, d1[, .(barcode, density_gene_1)], by = 'barcode')
  tdd <- merge(tdd,   d2[, .(barcode, density_gene_2)], by = 'barcode')
  
  # --- Fetch normalized counts ---
  norm_counts <- FetchData(seurat_object, vars = c(gene_1, gene_2), layer = 'data')
  setDT(norm_counts, keep.rownames = 'barcode')
  setnames(norm_counts, old = c(gene_1, gene_2), new = c('norm_gene_1', 'norm_gene_2'))
  tdd <- merge(tdd, norm_counts, by = 'barcode')
  
  # --- Main connected component ---
  coords   <- as.matrix(tdd[, .(STX_1, STX_2)])
  nb_full  <- spdep::knn2nb(spdep::knearneigh(coords, k = 6))
  nb_graph <- igraph::graph_from_adj_list(nb_full)
  components <- igraph::components(nb_graph)
  main_idx <- which(components$membership == which.max(components$csize))
  tdd_main <- tdd[main_idx, ]
  
  # --- Spatial neighbor weights ---
  coords_main <- as.matrix(tdd_main[, .(STX_1, STX_2)])
  nb_main <- spdep::knn2nb(spdep::knearneigh(coords_main, k = 6))
  lw_main <- spdep::nb2listw(nb_main, style = 'W')
  
  # --- Lee's L: KDE densities (supracellular, smoothed) ---
  lee_kde  <- spdep::lee.test(tdd_main$density_gene_1, tdd_main$density_gene_2, lw_main)
  moran_g1_kde <- spdep::moran.test(tdd_main$density_gene_1, lw_main)
  moran_g2_kde <- spdep::moran.test(tdd_main$density_gene_2, lw_main)
  
  # --- Lee's L: normalized counts (spot-level, unsmoothed) ---
  lee_raw  <- spdep::lee.test(tdd_main$norm_gene_1, tdd_main$norm_gene_2, lw_main)
  moran_g1_raw <- spdep::moran.test(tdd_main$norm_gene_1, lw_main)
  moran_g2_raw <- spdep::moran.test(tdd_main$norm_gene_2, lw_main)
  
  message(sprintf("--- Spatial statistics: %s vs %s ---", gene_1, gene_2))
  message("-- KDE density (supracellular) --")
  message(sprintf("Moran's I (%s): %.3f, p = %.4f",            gene_1, moran_g1_kde$estimate[1], moran_g1_kde$p.value))
  message(sprintf("Moran's I (%s): %.3f, p = %.4f",            gene_2, moran_g2_kde$estimate[1], moran_g2_kde$p.value))
  message(sprintf("Lee's L: %.3f, p = %.4f",                           lee_kde$estimate[1],      lee_kde$p.value))
  message("-- Normalized counts (spot-level) --")
  message(sprintf("Moran's I (%s): %.3f, p = %.4f",            gene_1, moran_g1_raw$estimate[1], moran_g1_raw$p.value))
  message(sprintf("Moran's I (%s): %.3f, p = %.4f",            gene_2, moran_g2_raw$estimate[1], moran_g2_raw$p.value))
  message(sprintf("Lee's L: %.3f, p = %.4f",                           lee_raw$estimate[1],      lee_raw$p.value))
  
  # --- Contour helper ---
  make_contour_df <- function(data, z_col, nx = 200, ny = 200) {
    d <- data[is.finite(get(z_col))]
    interp_result <- akima::interp(
      x = d$STX_1, y = d$STX_2, z = d[[z_col]],
      nx = nx, ny = ny, duplicate = 'mean')
    cdf <- expand.grid(STX_1 = interp_result$x, STX_2 = interp_result$y)
    cdf$z <- as.vector(interp_result$z)
    nn_dist <- FNN::knnx.dist(
      data  = as.matrix(data[, .(STX_1, STX_2)]),
      query = as.matrix(cdf[, c('STX_1', 'STX_2')]),
      k = 1)
    cdf$z[nn_dist > contour_dist_threshold] <- NA
    cdf
  }
  
  # --- Build three plots ---
  lee_label_kde <- sprintf("Lee's L (KDE) = %.3f, p = %.4f",  lee_kde$estimate[1], lee_kde$p.value)
  lee_label_raw <- sprintf("Lee's L (counts) = %.3f, p = %.4f", lee_raw$estimate[1], lee_raw$p.value)
  
  p1 <- ggscatter(tdd, x = 'STX_1', y = 'STX_2', color = 'density_gene_1',
                  size = pt.size, legend = 'right') +
    scale_color_continuous(palette = 'plasma', na.value = 'grey99', name = gene_1) +
    coord_fixed() + ggtitle(gene_1)
  
  p2 <- ggscatter(tdd, x = 'STX_1', y = 'STX_2', color = 'density_gene_2',
                 size = pt.size, legend = 'right') +
    scale_color_continuous(palette = 'plasma', na.value = 'grey99', name = gene_2) +
    coord_fixed() + ggtitle(gene_2)
  
  pj <- ggscatter(tdd, x = 'STX_1', y = 'STX_2', color = 'joint_density',
                  size = pt.size, legend = 'right') +
    scale_color_continuous(palette = 'plasma', na.value = 'grey99',
                           name = paste0(gene_1, '\n+ ', gene_2)) +
    coord_fixed() +
    ggtitle(paste0('Joint density\n', lee_label_kde, '\n', lee_label_raw))
  
  if (show_contours) {
    contour_g1 <- make_contour_df(tdd_main, 'density_gene_1')
    contour_g2 <- make_contour_df(tdd_main, 'density_gene_2')
    
    p1 <- p1 + geom_contour(data = contour_g1, aes(x = STX_1, y = STX_2, z = z),
                            color = 'white', linewidth = 0.6, na.rm = TRUE)
    p2 <- p2 + geom_contour(data = contour_g2, aes(x = STX_1, y = STX_2, z = z),
                            color = 'white', linewidth = 0.6, na.rm = TRUE)
   # pj <- pj +
    #  geom_contour(data = contour_g1, aes(x = STX_1, y = STX_2, z = z),
    #               color = 'white', linewidth = 0.6, na.rm = TRUE) +
    #  geom_contour(data = contour_g2, aes(x = STX_1, y = STX_2, z = z),
    #               color = 'red', linewidth = 0.6, na.rm = TRUE)
  }
  
  p1 + p2 + pj + patchwork::plot_layout(ncol = 3)
}