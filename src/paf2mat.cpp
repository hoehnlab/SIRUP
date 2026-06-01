// [[Rcpp::plugins(cpp17)]]
// [[Rcpp::plugins(openmp)]]
// [[Rcpp::depends(Rcpp)]]

#include <Rcpp.h>

#include "paf_contact_core.hpp"

using namespace Rcpp;

//' @export
// [[Rcpp::export]]
SEXP streamPaf2mat(std::string path,
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

  return rad_paf::stream_paf_contacts_to_sparse(path, cfg, symmetric);
}

//' @export
// [[Rcpp::export]]
SEXP streamPaf2matAnch(std::string path,
                         std::string target = "",
                         double region_start = 1.0,
                         double region_end = 0.0,
                         int bin_bp = 250,
                         int min_mapq = 0,
                         int min_anchor_bp = 1,
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
  cfg.min_anchor_bp = std::max(1, min_anchor_bp);
  cfg.chunk_lines = std::max(1, chunk_lines);
  cfg.threads = threads;
  cfg.verbose = verbose;

  return rad_paf::stream_paf_contacts_to_sparse(path, cfg, symmetric);
}

