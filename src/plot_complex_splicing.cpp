// [[Rcpp::plugins(cpp17)]]
// [[Rcpp::plugins(openmp)]]
// [[Rcpp::depends(Rcpp)]]

#include <Rcpp.h>

#include "paf_contact_core.hpp"

using namespace Rcpp;

static inline bool has_col(DataFrame df, const std::string& nm) {
  CharacterVector nms = df.names();
  for (R_xlen_t i = 0; i < nms.size(); ++i) {
    if (nms[i] == NA_STRING) {
      continue;
    }
    if (std::string_view(CHAR(nms[i])) == std::string_view(nm)) {
      return true;
    }
  }
  return false;
}

static inline NumericVector get_num_col(DataFrame df, const std::string& nm) {
  SEXP col = df[nm];
  if (Rf_isNull(col)) {
    stop("Column '%s' not found", nm.c_str());
  }
  if (TYPEOF(col) == INTSXP) {
    return as<NumericVector>(IntegerVector(col));
  }
  if (TYPEOF(col) == REALSXP) {
    return NumericVector(col);
  }
  stop("Column '%s' must be integer or numeric", nm.c_str());
  return NumericVector();
}

static inline CharacterVector get_chr_col(DataFrame df, const std::string& nm) {
  SEXP col = df[nm];
  if (Rf_isNull(col)) {
    stop("Column '%s' not found", nm.c_str());
  }
  if (TYPEOF(col) != STRSXP) {
    stop("Column '%s' must be character", nm.c_str());
  }
  return CharacterVector(col);
}

static inline std::string_view normalize_cg(std::string_view cg) {
  if (cg.size() >= 5 && cg.substr(0, 5) == std::string_view("cg:Z:")) {
    return cg.substr(5);
  }
  return cg;
}

//' @export
// [[Rcpp::export]]
SEXP spliceContactsFromPaf(std::string paf_path,
                              std::string target = "",
                              double region_start = 1.0,
                              double region_end = 0.0,
                              int bin_bp = 250,
                              int min_mapq = 0,
                              int chunk_lines = 200000,
                              bool symmetric = true,
                              int threads = 0,
                              bool verbose = false) {
  if (bin_bp <= 0) {
    stop("bin_bp must be > 0");
  }
  if (region_start < 1.0) {
    region_start = 1.0;
  }

  rad_paf::configure_openmp_threads(threads, verbose);

  rad_paf::StreamConfig cfg;
  cfg.target = target;
  cfg.region_start1b = static_cast<int64_t>(region_start);
  cfg.region_end1b = (region_end <= 0.0) ? 0 : static_cast<int64_t>(region_end);
  cfg.bin_bp = bin_bp;
  cfg.min_mapq = min_mapq;
  cfg.min_anchor_bp = 1;
  cfg.chunk_lines = std::max(1, chunk_lines);
  cfg.threads = threads;
  cfg.verbose = verbose;

  return rad_paf::stream_paf_contacts_to_sparse(paf_path, cfg, symmetric);
}

//' @export
// [[Rcpp::export]]
SEXP spliceContactsFromPafDt(
    DataFrame df,
    std::string target = "",
    double region_start = 1.0,
    double region_end = 0.0,
    int bin_bp = 250,
    int min_mapq = 0,
    int min_anchor_bp = 1,
    int chunk_rows = 1000000,
    bool symmetric = true,
    int threads = 0,
    bool verbose = false,
    std::string col_tname = "tname",
    std::string col_tstart = "tstart",
    std::string col_tend = "tend",
    std::string col_cg = "cg",
    std::string col_mapq = "mapq",
    std::string col_strand = "strand") {
  if (bin_bp <= 0) {
    stop("bin_bp must be > 0");
  }
  if (region_start < 1.0) {
    region_start = 1.0;
  }

  rad_paf::configure_openmp_threads(threads, verbose);

  const int64_t n = df.nrows();
  if (n == 0) {
    stop("Empty data.frame");
  }

  CharacterVector v_tname = get_chr_col(df, col_tname);
  NumericVector v_tstart = get_num_col(df, col_tstart);
  CharacterVector v_cg = get_chr_col(df, col_cg);

  const bool have_mapq = has_col(df, col_mapq);
  NumericVector v_mapq = have_mapq ? get_num_col(df, col_mapq) : NumericVector(n, 60.0);

  (void)col_tend;
  (void)col_strand;

  rad_paf::StreamConfig cfg;
  cfg.target = target;
  cfg.region_start1b = static_cast<int64_t>(region_start);
  cfg.region_end1b = (region_end <= 0.0) ? 0 : static_cast<int64_t>(region_end);
  cfg.bin_bp = bin_bp;
  cfg.min_mapq = min_mapq;
  cfg.min_anchor_bp = std::max(1, min_anchor_bp);
  cfg.chunk_lines = 0;
  cfg.threads = threads;
  cfg.verbose = verbose;

  std::vector<rad_paf::Pair> global_agg;
  std::vector<rad_paf::Pair> merge_buf;
  uint32_t max_bin_seen = 0;

  const bool filter_target = !target.empty();
  const std::string_view target_view(target);

  auto process_range = [&](int64_t a, int64_t b) {
    if (a >= b) {
      return;
    }

    int nt = 1;
#ifdef _OPENMP
#pragma omp parallel
    {
#pragma omp single
      { nt = omp_get_num_threads(); }
    }
#endif

    std::vector<std::vector<uint64_t>> tls_keys(static_cast<size_t>(nt));

#ifdef _OPENMP
#pragma omp parallel
#endif
    {
#ifdef _OPENMP
      const int tid = omp_get_thread_num();
#else
      const int tid = 0;
#endif
      std::vector<uint64_t>& keys = tls_keys[static_cast<size_t>(tid)];
      keys.reserve(static_cast<size_t>((b - a) / (nt * 2)) + 1024);

      std::vector<rad_paf::Block> blocks;
      blocks.reserve(8);

#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
      for (int64_t i = a; i < b; ++i) {
        if (filter_target) {
          SEXP tname = v_tname[i];
          if (tname == NA_STRING) {
            continue;
          }
          if (std::string_view(CHAR(tname)) != target_view) {
            continue;
          }
        }

        if (have_mapq && min_mapq > 0) {
          const double mq = v_mapq[i];
          if (NumericVector::is_na(mq) || mq < min_mapq) {
            continue;
          }
        }

        SEXP cg = v_cg[i];
        if (cg == NA_STRING) {
          continue;
        }
        std::string_view cgsv = normalize_cg(std::string_view(CHAR(cg)));
        if (cgsv.empty()) {
          continue;
        }

        const double d_tstart = v_tstart[i];
        if (NumericVector::is_na(d_tstart)) {
          continue;
        }
        const int64_t tstart1b = static_cast<int64_t>(d_tstart) + 1;

        rad_paf::parse_cigar_blocks(tstart1b, cgsv, blocks);
        rad_paf::append_junction_pairs(blocks, cfg, keys);
      }
    }

    size_t total = 0;
    for (const auto& v : tls_keys) {
      total += v.size();
    }
    if (total == 0) {
      return;
    }

    std::vector<uint64_t> all;
    all.reserve(total);
    for (const auto& v : tls_keys) {
      all.insert(all.end(), v.begin(), v.end());
    }

    for (uint64_t packed : all) {
      const uint32_t j = static_cast<uint32_t>(packed & 0xffffffffu);
      if (j > max_bin_seen) {
        max_bin_seen = j;
      }
    }

    std::sort(all.begin(), all.end());
    std::vector<rad_paf::Pair> local_agg;
    local_agg.reserve(all.size() / 2 + 1);
    rad_paf::reduce_sorted_keys(all, local_agg);

    if (global_agg.empty()) {
      global_agg.swap(local_agg);
    } else {
      rad_paf::merge_counts(global_agg, local_agg, merge_buf);
      global_agg.swap(merge_buf);
      merge_buf.clear();
    }
  };

  if (chunk_rows <= 0) {
    chunk_rows = 1000000;
  }

  for (int64_t a = 0; a < n; a += chunk_rows) {
    const int64_t b = std::min<int64_t>(n, a + chunk_rows);
    if (verbose) {
      Rcpp::Rcout << "Processing rows [" << a << "," << b << ")\n";
    }
    process_range(a, b);
  }

  if (global_agg.empty()) {
    stop("No edges produced (check cg/mapq/target filters).");
  }

  const uint32_t n_bins = rad_paf::infer_n_bins(cfg, max_bin_seen);
  return rad_paf::pairs_to_sparse_matrix(global_agg, n_bins, symmetric);
}

