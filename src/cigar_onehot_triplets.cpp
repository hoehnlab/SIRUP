// [[Rcpp::plugins(cpp17)]]
// [[Rcpp::depends(Rcpp)]]

#include <Rcpp.h>
#include <algorithm>
#include <string>
#include <vector>

using namespace Rcpp;

inline void emit_cigar_triplets(const std::string& cigar,
                                int aln_start,
                                int region_start,
                                int region_end,
                                int row1,
                                std::vector<int>& i_out,
                                std::vector<int>& j_out) {
  long long ref_pos = aln_start;
  const char* p = cigar.c_str();
  const char* const pend = p + cigar.size();

  while (p < pend) {
    long long len = 0;
    if (*p < '0' || *p > '9') {
      ++p;
      continue;
    }
    while (p < pend && *p >= '0' && *p <= '9') {
      len = len * 10 + (*p - '0');
      ++p;
    }
    if (p >= pend) {
      break;
    }

    const char op = *p++;
    switch (op) {
    case 'M':
    case '=': {
      const long long seg_start = ref_pos;
      const long long seg_end = ref_pos + len;
      const long long cs = std::max(seg_start, static_cast<long long>(region_start));
      const long long ce = std::min(seg_end, static_cast<long long>(region_end));
      if (ce > cs) {
        const int base_col0 = static_cast<int>(cs - region_start);
        const int n = static_cast<int>(ce - cs);
        i_out.reserve(i_out.size() + n);
        j_out.reserve(j_out.size() + n);
        for (int k = 0; k < n; ++k) {
          i_out.push_back(row1);
          j_out.push_back(base_col0 + k + 1);
        }
      }
      ref_pos += len;
      break;
    }
    case 'D':
    case 'N':
    case 'X':
      ref_pos += len;
      break;
    default:
      break;
    }
  }
}

//' @export
// [[Rcpp::export]]
Rcpp::List cigarOnehotTriplets(Rcpp::List cigars,
                                 Rcpp::IntegerVector aln_starts,
                                 Rcpp::IntegerVector aln_ends,
                                 int region_start,
                                 int region_end) {
  const int n = cigars.size();
  if (n != aln_starts.size() || n != aln_ends.size()) {
    Rcpp::stop("Lengths of cigars, aln_starts, aln_ends must match.");
  }
  if (region_end < region_start) {
    Rcpp::stop("region_end must be >= region_start.");
  }

  std::vector<int> i_out;
  std::vector<int> j_out;
  i_out.reserve(1024);
  j_out.reserve(1024);

  for (int r = 0; r < n; ++r) {
    const std::string cigar = Rcpp::as<std::string>(cigars[r]);
    const int aln_start = aln_starts[r];
    emit_cigar_triplets(cigar, aln_start, region_start, region_end, r + 1, i_out, j_out);
  }

  return Rcpp::List::create(
    Rcpp::Named("i") = Rcpp::IntegerVector(i_out.begin(), i_out.end()),
    Rcpp::Named("j") = Rcpp::IntegerVector(j_out.begin(), j_out.end()),
    Rcpp::Named("nrow") = n,
    Rcpp::Named("ncol") = region_end - region_start
  );
}

