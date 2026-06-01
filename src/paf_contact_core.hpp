#pragma once

#include <Rcpp.h>
#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <memory>
#include <string>
#include <string_view>
#include <utility>
#include <vector>
#include <zlib.h>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace rad_paf {

struct Block {
  int64_t s;
  int64_t e;
};

struct Pair {
  uint64_t key;
  uint32_t cnt;
};

struct Entry {
  uint32_t j;
  uint32_t i;
  double x;
};

struct PafLite {
  std::string_view tname;
  int64_t tstart = -1;  // 0-based
  int64_t tend = -1;    // 0-based end
  char strand = '+';
  int mapq = 0;
  std::string_view cg;
};

struct StreamConfig {
  std::string target;
  int64_t region_start1b = 1;
  int64_t region_end1b = 0;
  int bin_bp = 250;
  int min_mapq = 0;
  int min_anchor_bp = 1;
  int chunk_lines = 200000;
  int threads = 0;
  bool verbose = false;
};

inline void configure_openmp_threads(int threads, bool verbose) {
#ifdef _OPENMP
  if (threads > 0) {
    omp_set_num_threads(threads);
  }
#else
  if (threads > 0 && verbose) {
    Rcpp::Rcout << "OpenMP not available; running single-threaded.\n";
  }
#endif
}

inline void parse_cigar_blocks(int64_t ref_start1b,
                               std::string_view cig,
                               std::vector<Block>& out) {
  out.clear();
  int64_t rpos = ref_start1b;
  int64_t bs = -1;
  int64_t be = -1;
  auto flush = [&]() {
    if (bs >= 0) {
      out.push_back({bs, be});
      bs = -1;
      be = -1;
    }
  };

  uint32_t n = 0;
  for (char c : cig) {
    if (c >= '0' && c <= '9') {
      n = n * 10 + static_cast<uint32_t>(c - '0');
      continue;
    }
    if (n == 0) {
      continue;
    }

    switch (c) {
    case 'M':
    case '=':
    case 'X':
      if (bs < 0) {
        bs = rpos;
      }
      be = rpos + static_cast<int64_t>(n) - 1;
      rpos += static_cast<int64_t>(n);
      break;
    case 'D':
      rpos += static_cast<int64_t>(n);
      break;
    case 'N':
      flush();
      rpos += static_cast<int64_t>(n);
      break;
    case 'I':
    case 'S':
    case 'H':
    case 'P':
      break;
    default:
      break;
    }
    n = 0;
  }
  flush();
}

inline int32_t bin_of(int64_t pos1b, int64_t region_start1b, int32_t bin_bp) {
  if (pos1b < region_start1b) {
    pos1b = region_start1b;
  }
  int64_t zero = pos1b - region_start1b;
  if (zero < 0) {
    zero = 0;
  }
  return static_cast<int32_t>(zero / bin_bp);
}

inline uint64_t pack_key(uint32_t i, uint32_t j) {
  return (static_cast<uint64_t>(i) << 32) | static_cast<uint64_t>(j);
}

inline void reduce_sorted_keys(const std::vector<uint64_t>& keys,
                               std::vector<Pair>& out) {
  out.clear();
  if (keys.empty()) {
    return;
  }

  uint64_t cur = keys[0];
  uint32_t c = 1;
  for (size_t k = 1; k < keys.size(); ++k) {
    if (keys[k] == cur) {
      ++c;
    } else {
      out.push_back({cur, c});
      cur = keys[k];
      c = 1;
    }
  }
  out.push_back({cur, c});
}

inline void merge_counts(const std::vector<Pair>& a,
                         const std::vector<Pair>& b,
                         std::vector<Pair>& dst) {
  dst.clear();
  dst.reserve(a.size() + b.size());

  size_t i = 0;
  size_t j = 0;
  while (i < a.size() && j < b.size()) {
    if (a[i].key < b[j].key) {
      dst.push_back(a[i++]);
    } else if (b[j].key < a[i].key) {
      dst.push_back(b[j++]);
    } else {
      dst.push_back({a[i].key, a[i].cnt + b[j].cnt});
      ++i;
      ++j;
    }
  }
  while (i < a.size()) {
    dst.push_back(a[i++]);
  }
  while (j < b.size()) {
    dst.push_back(b[j++]);
  }
}

inline int anchor_overlap_bp_donor(const Block& a, int64_t region_start1b, int32_t bin_bp) {
  uint32_t bi = static_cast<uint32_t>(bin_of(a.e, region_start1b, bin_bp));
  int64_t bs = region_start1b + static_cast<int64_t>(bi) * bin_bp;
  int64_t be = bs + bin_bp - 1;
  int64_t ov = std::min<int64_t>(a.e, be) - std::max<int64_t>(a.s, bs) + 1;
  return static_cast<int>(std::max<int64_t>(0, ov));
}

inline int anchor_overlap_bp_acceptor(const Block& b, int64_t region_start1b, int32_t bin_bp) {
  uint32_t bj = static_cast<uint32_t>(bin_of(b.s, region_start1b, bin_bp));
  int64_t bs = region_start1b + static_cast<int64_t>(bj) * bin_bp;
  int64_t be = bs + bin_bp - 1;
  int64_t ov = std::min<int64_t>(b.e, be) - std::max<int64_t>(b.s, bs) + 1;
  return static_cast<int>(std::max<int64_t>(0, ov));
}

inline bool parse_paf_line(std::string& line, PafLite& out) {
  size_t n = line.size();
  size_t p = 0;

  auto next_field = [&](size_t& pos) -> std::pair<size_t, size_t> {
    if (pos >= n) {
      return {n, n};
    }
    size_t s = pos;
    while (pos < n && line[pos] != '\t' && line[pos] != '\n' && line[pos] != '\r') {
      ++pos;
    }
    size_t e = pos;
    if (pos < n && line[pos] == '\t') {
      ++pos;
    }
    return {s, e};
  };

  std::vector<std::pair<size_t, size_t>> fields;
  fields.reserve(16);
  while (p < n) {
    auto se = next_field(p);
    if (se.first == se.second && se.first == n) {
      break;
    }
    fields.push_back(se);
  }

  if (fields.size() < 12) {
    return false;
  }

  out = PafLite();

  if (fields[4].second <= fields[4].first) {
    return false;
  }
  out.strand = line[fields[4].first];

  out.tname = std::string_view(line.data() + fields[5].first,
                               fields[5].second - fields[5].first);

  out.tstart = std::strtoll(line.c_str() + fields[7].first, nullptr, 10);
  out.tend = std::strtoll(line.c_str() + fields[8].first, nullptr, 10);
  out.mapq = std::atoi(line.c_str() + fields[11].first);

  for (size_t k = 12; k < fields.size(); ++k) {
    const size_t s = fields[k].first;
    const size_t e = fields[k].second;
    if (e > s + 5 &&
        line[s + 0] == 'c' && line[s + 1] == 'g' && line[s + 2] == ':' &&
        line[s + 3] == 'Z' && line[s + 4] == ':') {
      out.cg = std::string_view(line.data() + s + 5, e - (s + 5));
      break;
    }
  }

  return !out.cg.empty();
}

struct ILineReader {
  virtual bool getline(std::string& out) = 0;
  virtual ~ILineReader() {}
};

struct IfsReader : ILineReader {
  std::ifstream f;
  explicit IfsReader(const std::string& path) : f(path.c_str()) {}

  bool getline(std::string& out) override {
    if (!std::getline(f, out)) {
      return false;
    }
    if (!out.empty() && out.back() == '\r') {
      out.pop_back();
    }
    return true;
  }

  bool ok() const {
    return f.good();
  }
};

struct GzReader : ILineReader {
  gzFile f = nullptr;
  std::vector<char> buf;

  explicit GzReader(const std::string& path, size_t buf_sz = (1 << 16))
      : f(gzopen(path.c_str(), "rb")), buf(buf_sz, 0) {}

  bool getline(std::string& out) override {
    out.clear();
    if (!f) {
      return false;
    }

    for (;;) {
      char* r = gzgets(f, buf.data(), static_cast<int>(buf.size()));
      if (!r) {
        return !out.empty();
      }

      const size_t len = std::strlen(r);
      if (len > 0) {
        out.append(r, len);
      }

      if (!out.empty() && (out.back() == '\n' || out.back() == '\r')) {
        while (!out.empty() && (out.back() == '\n' || out.back() == '\r')) {
          out.pop_back();
        }
        return true;
      }

      if (len < buf.size() - 1) {
        return true;
      }
    }
  }

  bool ok() const {
    return f != nullptr;
  }

  ~GzReader() {
    if (f) {
      gzclose(f);
    }
  }
};

inline bool ends_with(const std::string& s, const std::string& suffix) {
  return s.size() >= suffix.size() &&
    std::equal(suffix.rbegin(), suffix.rend(), s.rbegin());
}

inline void append_junction_pairs(const std::vector<Block>& blocks,
                                  const StreamConfig& cfg,
                                  std::vector<uint64_t>& keys) {
  if (blocks.size() < 2) {
    return;
  }

  const int64_t R1 = cfg.region_start1b;
  const int64_t R2 = cfg.region_end1b;

  for (size_t k = 0; k + 1 < blocks.size(); ++k) {
    const Block& a = blocks[k];
    const Block& b = blocks[k + 1];

    if (R2 > 0) {
      if (a.e < R1 || a.e > R2) {
        continue;
      }
      if (b.s < R1 || b.s > R2) {
        continue;
      }
    }

    if (cfg.min_anchor_bp > 1) {
      const int left_ov = anchor_overlap_bp_donor(a, R1, cfg.bin_bp);
      const int right_ov = anchor_overlap_bp_acceptor(b, R1, cfg.bin_bp);
      if (left_ov < cfg.min_anchor_bp || right_ov < cfg.min_anchor_bp) {
        continue;
      }
    }

    uint32_t bi = static_cast<uint32_t>(bin_of(a.e, R1, cfg.bin_bp));
    uint32_t bj = static_cast<uint32_t>(bin_of(b.s, R1, cfg.bin_bp));
    if (bi > bj) {
      std::swap(bi, bj);
    }
    keys.push_back(pack_key(bi, bj));
  }
}

inline void process_paf_chunk(std::vector<std::string>& chunk,
                              const StreamConfig& cfg,
                              std::vector<Pair>& global_agg,
                              std::vector<Pair>& merge_buf,
                              uint32_t& max_bin_seen) {
  if (chunk.empty()) {
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
  const bool filter_target = !cfg.target.empty();
  const std::string_view target_view(cfg.target);

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
    keys.reserve(chunk.size() / (static_cast<size_t>(nt) * 2) + 1024);

    std::vector<Block> blocks;
    blocks.reserve(8);

#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
    for (int i = 0; i < static_cast<int>(chunk.size()); ++i) {
      PafLite pl;
      if (!parse_paf_line(chunk[static_cast<size_t>(i)], pl)) {
        continue;
      }
      if (cfg.min_mapq > 0 && pl.mapq < cfg.min_mapq) {
        continue;
      }
      if (filter_target && pl.tname != target_view) {
        continue;
      }

      const int64_t tstart1b = pl.tstart + 1;
      parse_cigar_blocks(tstart1b, pl.cg, blocks);
      append_junction_pairs(blocks, cfg, keys);
    }
  }

  size_t total = 0;
  for (const auto& v : tls_keys) {
    total += v.size();
  }

  if (total == 0) {
    chunk.clear();
    return;
  }

  std::vector<uint64_t> all;
  all.reserve(total);
  for (auto& v : tls_keys) {
    all.insert(all.end(), v.begin(), v.end());
  }

  for (uint64_t packed : all) {
    uint32_t j = static_cast<uint32_t>(packed & 0xffffffffu);
    if (j > max_bin_seen) {
      max_bin_seen = j;
    }
  }

  std::sort(all.begin(), all.end());
  std::vector<Pair> local_agg;
  local_agg.reserve(all.size() / 2 + 1);
  reduce_sorted_keys(all, local_agg);

  if (global_agg.empty()) {
    global_agg.swap(local_agg);
  } else {
    merge_counts(global_agg, local_agg, merge_buf);
    global_agg.swap(merge_buf);
    merge_buf.clear();
  }

  chunk.clear();
}

inline void stream_paf_to_pairs(const std::string& path,
                                const StreamConfig& cfg,
                                std::vector<Pair>& global_agg,
                                std::vector<Pair>& merge_buf,
                                uint32_t& max_bin_seen) {
  std::unique_ptr<ILineReader> reader;
  if (ends_with(path, ".gz")) {
    std::unique_ptr<GzReader> gz(new GzReader(path));
    if (!gz->ok()) {
      Rcpp::stop("Failed to open gzip PAF: %s", path.c_str());
    }
    reader = std::move(gz);
  } else {
    std::unique_ptr<IfsReader> ifs(new IfsReader(path));
    if (!ifs->ok()) {
      Rcpp::stop("Failed to open PAF: %s", path.c_str());
    }
    reader = std::move(ifs);
  }

  std::vector<std::string> lines;
  lines.reserve(static_cast<size_t>(std::max(1, cfg.chunk_lines)));

  std::string line;
  while (reader->getline(line)) {
    lines.emplace_back(std::move(line));
    if (static_cast<int>(lines.size()) >= cfg.chunk_lines) {
      if (cfg.verbose) {
        Rcpp::Rcout << "Processing " << lines.size() << " lines...\n";
      }
      process_paf_chunk(lines, cfg, global_agg, merge_buf, max_bin_seen);
    }
  }

  if (!lines.empty()) {
    if (cfg.verbose) {
      Rcpp::Rcout << "Processing " << lines.size() << " lines...\n";
    }
    process_paf_chunk(lines, cfg, global_agg, merge_buf, max_bin_seen);
  }
}

inline uint32_t infer_n_bins(const StreamConfig& cfg, uint32_t max_bin_seen) {
  if (cfg.region_end1b > 0) {
    const int64_t region_len = cfg.region_end1b - cfg.region_start1b + 1;
    if (region_len <= 0) {
      Rcpp::stop("region_end must be >= region_start");
    }
    return static_cast<uint32_t>((region_len + cfg.bin_bp - 1) / cfg.bin_bp);
  }
  return max_bin_seen + 1;
}

inline SEXP pairs_to_sparse_matrix(const std::vector<Pair>& pairs,
                                   uint32_t n_bins,
                                   bool symmetric) {
  std::vector<Entry> entries;
  entries.reserve(pairs.size());
  for (const Pair& p : pairs) {
    const uint32_t i = static_cast<uint32_t>(p.key >> 32);
    const uint32_t j = static_cast<uint32_t>(p.key & 0xffffffffu);
    entries.push_back({j, i, static_cast<double>(p.cnt)});
  }

  std::sort(entries.begin(), entries.end(), [](const Entry& a, const Entry& b) {
    if (a.j != b.j) {
      return a.j < b.j;
    }
    return a.i < b.i;
  });

  if (symmetric) {
    std::vector<int> i;
    std::vector<int> p(n_bins + 1, 0);
    std::vector<double> x;
    i.reserve(entries.size());
    x.reserve(entries.size());

    size_t k = 0;
    for (uint32_t j = 0; j < n_bins; ++j) {
      p[j] = static_cast<int>(i.size());
      while (k < entries.size() && entries[k].j == j) {
        i.push_back(static_cast<int>(entries[k].i));
        x.push_back(entries[k].x);
        ++k;
      }
    }
    p[n_bins] = static_cast<int>(i.size());

    Rcpp::S4 mat("dsCMatrix");
    mat.slot("i") = Rcpp::IntegerVector(i.begin(), i.end());
    mat.slot("p") = Rcpp::IntegerVector(p.begin(), p.end());
    mat.slot("x") = Rcpp::NumericVector(x.begin(), x.end());
    mat.slot("Dim") = Rcpp::IntegerVector::create(n_bins, n_bins);
    mat.slot("Dimnames") = Rcpp::List::create(R_NilValue, R_NilValue);
    mat.slot("uplo") = Rcpp::CharacterVector::create("U");
    mat.slot("factors") = Rcpp::List::create();
    return mat;
  }

  std::vector<int> col_nnz(n_bins, 0);
  for (const auto& e : entries) {
    ++col_nnz[e.j];
  }
  for (const auto& e : entries) {
    if (e.i != e.j) {
      ++col_nnz[e.i];
    }
  }

  std::vector<int> p(n_bins + 1, 0);
  for (uint32_t j = 0; j < n_bins; ++j) {
    p[j + 1] = p[j] + col_nnz[j];
  }

  const int nnz = p[n_bins];
  std::vector<int> i(nnz);
  std::vector<double> x(nnz);
  std::vector<int> cursor = p;

  for (const auto& e : entries) {
    const int idx = cursor[e.j]++;
    i[idx] = static_cast<int>(e.i);
    x[idx] = e.x;
  }
  for (const auto& e : entries) {
    if (e.i != e.j) {
      const int idx = cursor[e.i]++;
      i[idx] = static_cast<int>(e.j);
      x[idx] = e.x;
    }
  }

  Rcpp::S4 mat("dgCMatrix");
  mat.slot("i") = Rcpp::IntegerVector(i.begin(), i.end());
  mat.slot("p") = Rcpp::IntegerVector(p.begin(), p.end());
  mat.slot("x") = Rcpp::NumericVector(x.begin(), x.end());
  mat.slot("Dim") = Rcpp::IntegerVector::create(n_bins, n_bins);
  mat.slot("Dimnames") = Rcpp::List::create(R_NilValue, R_NilValue);
  mat.slot("factors") = Rcpp::List::create();
  return mat;
}

inline SEXP stream_paf_contacts_to_sparse(const std::string& path,
                                          const StreamConfig& cfg,
                                          bool symmetric) {
  std::vector<Pair> global_agg;
  std::vector<Pair> merge_buf;
  uint32_t max_bin_seen = 0;

  stream_paf_to_pairs(path, cfg, global_agg, merge_buf, max_bin_seen);

  if (global_agg.empty()) {
    Rcpp::stop("No edges produced (check cg/mapq/target filters).");
  }

  const uint32_t n_bins = infer_n_bins(cfg, max_bin_seen);
  return pairs_to_sparse_matrix(global_agg, n_bins, symmetric);
}

}  // namespace rad_paf

