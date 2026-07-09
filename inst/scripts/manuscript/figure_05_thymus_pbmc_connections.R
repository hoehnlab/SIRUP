# Figure 5: thymus-PBMC clonal connections

library(SIRAP)
library(dplyr)
library(airr)
library(scoper)

sample_id <- Sys.getenv("SIRAP_SAMPLE_ID", "M5")
sample_lower <- tolower(sample_id)
nproc <- as.integer(Sys.getenv("SIRAP_NPROC", "8"))

thymus_dir <- Sys.getenv(
  "SIRAP_THYMUS_DIR",
  "/dartfs/rc/lab/H/HoehnK/Hunter/Thymus_VisiumHD"
)
walkthrough_dir <- file.path(thymus_dir, "SIRAP_walkthrough", sample_id)
final_airr_dir <- file.path(walkthrough_dir, "03_final_spatial_airr")
manuscript_out_dir <- file.path(thymus_dir, "SIRAP_manuscript", "figure_05", sample_id)
dir.create(manuscript_out_dir, recursive = TRUE, showWarnings = FALSE)

thymus_bcr_path <- file.path(final_airr_dir, sprintf("%s_bcrs_heavy_final_spatial_airr.tsv", sample_lower))
pbmc_month0_path <- Sys.getenv("SIRAP_PBMC_MONTH0_AIRR", "")
pbmc_month12_path <- Sys.getenv("SIRAP_PBMC_MONTH12_AIRR", "")

if (!file.exists(thymus_bcr_path)) {
  stop("Missing thymus BCR AIRR file: ", thymus_bcr_path)
}
if (!nzchar(pbmc_month0_path) || !file.exists(pbmc_month0_path)) {
  stop("Set SIRAP_PBMC_MONTH0_AIRR to an AIRR TSV for PBMC month 0.")
}
if (!nzchar(pbmc_month12_path) || !file.exists(pbmc_month12_path)) {
  stop("Set SIRAP_PBMC_MONTH12_AIRR to an AIRR TSV for PBMC month 12.")
}

thymus_bcr <- airr::read_rearrangement(thymus_bcr_path) |>
  mutate(compartment = "thymus", timepoint = "thymus", source_sample = sample_id)
pbmc_month0 <- airr::read_rearrangement(pbmc_month0_path) |>
  mutate(compartment = "PBMC", timepoint = "0mo", source_sample = "PBMC_0mo")
pbmc_month12 <- airr::read_rearrangement(pbmc_month12_path) |>
  mutate(compartment = "PBMC", timepoint = "12mo", source_sample = "PBMC_12mo")

combined_bcr <- bind_rows(thymus_bcr, pbmc_month0, pbmc_month12) |>
  filter(productive %in% TRUE)

required_cols <- c("sequence_id", "v_call", "j_call", "junction", "junction_length")
missing_cols <- setdiff(required_cols, colnames(combined_bcr))
if (length(missing_cols) > 0L) {
  stop("Combined AIRR table is missing required columns: ", paste(missing_cols, collapse = ", "))
}

combined_clones <- scoper::hierarchicalClones(
  combined_bcr,
  threshold = 0.1,
  cell_id = "sequence_id",
  summarize_clones = FALSE,
  nproc = nproc
)

saveRDS(combined_clones, file.path(manuscript_out_dir, "combined_thymus_pbmc_bcr_clones.rds"))

clone_presence <- combined_clones |>
  filter(!is.na(clone_id)) |>
  distinct(clone_id, compartment, timepoint, source_sample) |>
  count(clone_id, compartment, timepoint, name = "present_records") |>
  tidyr::pivot_wider(
    names_from = c(compartment, timepoint),
    values_from = present_records,
    values_fill = 0
  )

utils::write.table(
  clone_presence,
  file.path(manuscript_out_dir, "clone_presence_by_compartment.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

shared_clone_summary <- combined_clones |>
  filter(!is.na(clone_id)) |>
  mutate(is_thymus = compartment == "thymus", is_pbmc = compartment == "PBMC") |>
  group_by(clone_id) |>
  summarize(
    n_sequences = n(),
    n_thymus = sum(is_thymus),
    n_pbmc_0mo = sum(compartment == "PBMC" & timepoint == "0mo"),
    n_pbmc_12mo = sum(compartment == "PBMC" & timepoint == "12mo"),
    shared_thymus_pbmc = n_thymus > 0 & (n_pbmc_0mo > 0 | n_pbmc_12mo > 0),
    .groups = "drop"
  ) |>
  arrange(desc(shared_thymus_pbmc), desc(n_sequences))

utils::write.table(
  shared_clone_summary,
  file.path(manuscript_out_dir, "shared_clone_summary.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
