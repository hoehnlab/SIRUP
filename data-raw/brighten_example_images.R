# Make example spatial objects with the brightened H&E image. The existing example
# objects stay as they are.
#
# make_example_data.R builds example-data/mg_thymus_slice{1,2}_spatial.rds with the
# standard H&E image. The updated thymus SIRUP object (analysis/fig_4_..._11_19) was
# built from Seurat objects that use tissue_hires_image_brighter.png. That image has
# the same size and scale factors as the standard hires image. This script reduces
# the brightened hires image to Space Ranger lowres size, with the same block
# averaging as make_example_data.R. It puts the image into a copy of each example
# object and saves the copy as mg_thymus_slice{1,2}_spatial_brighter.rds. The
# expression, the spot coordinates, and the polygons do not change.
#
# After this script, run make_example_sirup_objects.R, and add the MD5 sums of the new
# files to R/utils.R (.sirup_example_files).
#
# Run on Discovery (the updated object is 2 GB; request about 48 GB of memory):
#   module load R/4.5.1; Rscript data-raw/brighten_example_images.R

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
})

repo <- "/dartfs/rc/lab/H/HoehnK/Hunter/spatial_bcr_tcr_package"
updated_path <- file.path(
  repo, "analysis/fig_4_bcr_publici_tonsil_11_19/02_build_sirup_objects/analysis/thymus_multisample_sirup.rds"
)
example_dir <- file.path(repo, "SIRUP", "example-data")
samples <- c("mg_thymus_slice1", "mg_thymus_slice2")

updated <- readRDS(updated_path)
message("sections in the updated object: ", paste(names(updated$spatial), collapse = ", "))

for (s in samples) {
  # The updated object names its sections thymus_slice1 and thymus_slice2.
  key <- grep(paste0(sub("^mg_", "", s), "$"), names(updated$spatial), value = TRUE)
  stopifnot(length(key) == 1L)
  new_image <- updated$spatial[[key]]@images[[1]]

  example_file <- file.path(example_dir, paste0(s, "_spatial.rds"))
  example <- readRDS(example_file)
  image_name <- Images(example)[[1]]
  img <- example@images[[image_name]]

  # Reduce the brightened hires raster to Space Ranger lowres size by block averaging.
  hires <- new_image@image
  shrink <- new_image@scale.factors$hires / new_image@scale.factors$lowres
  row_group <- ceiling(seq_len(dim(hires)[[1]]) / shrink)
  col_group <- ceiling(seq_len(dim(hires)[[2]]) / shrink)
  lowres <- array(0, dim = c(max(row_group), max(col_group), 3))
  for (channel in 1:3) {
    m <- rowsum(hires[, , channel], row_group) / as.vector(table(row_group))
    m <- t(rowsum(t(m), col_group) / as.vector(table(col_group)))
    lowres[, , channel] <- m
  }

  # The new image must have the same size as the old one, so the coordinates still fit.
  stopifnot(identical(dim(lowres), dim(img@image)))
  message(
    s, ": mean intensity ", round(mean(img@image), 4), " -> ", round(mean(lowres), 4),
    " (", paste(dim(lowres), collapse = " x "), ")"
  )

  img@image <- round(lowres * 255) / 255
  example[[image_name]] <- img
  saveRDS(example, file.path(example_dir, paste0(s, "_spatial_brighter.rds")), compress = "xz")
}
