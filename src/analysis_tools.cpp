#include "rad.h"
using namespace std;
using namespace Rcpp;

KSEQ_INIT(gzFile, gzread)

//' @export
// [[Rcpp::export]]
Rcpp::DataFrame read_seqs_cpp(std::string filename, std::string type, bool full_id = false, int max_reads = -1) {
  std::vector<std::string> ids;
  std::vector<std::string> full_ids;
  std::vector<std::string> seqs;
  std::vector<std::string> quals;
  bool is_fastq = (type == "fq" || type == "fastq");
  
  // Open the file
  gzFile fp = gzopen(filename.c_str(), "r");
  if (!fp) {
    Rcpp::stop("Could not open file: " + filename);
  }
  
  // Initialize kseq
  kseq_t *seq = kseq_init(fp);
  int read_counter = 0;
  int result;
  
  // Read sequences
  while ((result = kseq_read(seq)) >= 0) {
    // Handle full ID vs shortened ID
    std::string full_id_str = std::string(seq->name.s);
    if (seq->comment.l) {
      full_id_str += " " + std::string(seq->comment.s);
    }
    
    std::string short_id = full_id_str;
    size_t space_pos = short_id.find(' ');
    if (space_pos != std::string::npos) {
      short_id = short_id.substr(0, space_pos);
    }
    
    // Store sequence data
    ids.push_back(short_id);
    if (full_id) {
      full_ids.push_back(full_id_str);
    }
    seqs.push_back(std::string(seq->seq.s));
    
    // Handle quality scores for FASTQ
    if (is_fastq) {
      if (seq->qual.l == 0) {
        Rcpp::stop("Missing quality scores in FASTQ file");
      }
      quals.push_back(std::string(seq->qual.s));
    }
    
    read_counter++;
    if (max_reads > 0 && read_counter >= max_reads) {
      break;
    }
  }
  
  // Clean up
  kseq_destroy(seq);
  gzclose(fp);
  
  // Return appropriate DataFrame based on format and full_id option
  if (is_fastq) {
    if (full_id) {
      return Rcpp::DataFrame::create(
        Rcpp::Named("id") = ids,
        Rcpp::Named("full_id") = full_ids,
        Rcpp::Named("seq") = seqs,
        Rcpp::Named("qual") = quals
      );
    } else {
      return Rcpp::DataFrame::create(
        Rcpp::Named("id") = ids,
        Rcpp::Named("seq") = seqs,
        Rcpp::Named("qual") = quals
      );
    }
  } else {
    if (full_id) {
      return Rcpp::DataFrame::create(
        Rcpp::Named("id") = ids,
        Rcpp::Named("full_id") = full_ids,
        Rcpp::Named("seq") = seqs
      );
    } else {
      return Rcpp::DataFrame::create(
        Rcpp::Named("id") = ids,
        Rcpp::Named("seq") = seqs
      );
    }
  }
}

//' @export
// [[Rcpp::export]]
SEXP stream_fastqas(SEXP fn, std::string type = "fq", bool full_id = false, int max_reads = -1) {
  Rcpp::Function file_exists("file.exists");
  Rcpp::Function setDT("setDT", "data.table");
  Rcpp::Function pblapply("pblapply", "pbapply");
  Rcpp::Function rbindlist("rbindlist", "data.table");
  
  // Input validation
  if (type != "fa" && type != "fasta" && type != "fq" && type != "fastq") {
    Rcpp::stop("type must be one of: fa, fasta, fq, fastq");
  }
  
  // Normalize type to internal representation
  type = (type == "fasta") ? "fa" : (type == "fastq" ? "fq" : type);
  
  // Handle single file case
  auto read_single = [&](const std::string& filename) {
    if (!Rcpp::as<bool>(file_exists(filename))) {
      Rcpp::stop("File not found: " + filename);
    }
    
    Rcpp::DataFrame res = read_seqs_cpp(filename, type, full_id, max_reads);
    return setDT(res);
  };
  
  // Check if input is a vector/list
  if (TYPEOF(fn) == STRSXP && Rf_length(fn) > 1) {
    // Handle multiple files
    Rcpp::Function identity("identity");
    Rcpp::List out = pblapply(fn, identity);
    
    std::vector<SEXP> results;
    for (int i = 0; i < out.length(); i++) {
      results.push_back(read_single(Rcpp::as<std::string>(out[i])));
    }
    
    return rbindlist(Rcpp::wrap(results));
  } else {
    // Handle single file
    return read_single(Rcpp::as<std::string>(fn));
  }
}

//' @export
// [[Rcpp::export]]
NumericMatrix cmv_consensusMatrix(SEXP sequences) {
  CharacterVector seqs(sequences);
  int n_seq = seqs.length();
  if (n_seq == 0) return NumericMatrix(0);
  
  // Get length of first sequence and convert to std::string
  int seq_length = 0;
  for (int i = 0; i < n_seq; i++) {
    std::string curr_seq = std::string(seqs[i]);
    seq_length = std::max(seq_length, (int)curr_seq.length());
  }  
  // Initialize 4x seq_length matrix
  NumericMatrix consensus(4, seq_length);
  rownames(consensus) = CharacterVector::create("A", "C", "G", "T");
  // Count occurrences of each base at each position
  for (int i = 0; i < n_seq; i++) {
    std::string curr_seq = std::string(seqs[i]);
    
    for (int j = 0; j < seq_length; j++) {
      char base = curr_seq[j];
      switch(base) {
      case 'A':
      case 'a':
        consensus(0, j)++;
        break;
      case 'C':
      case 'c':
        consensus(1, j)++;
        break;
      case 'G':
      case 'g':
        consensus(2, j)++;
        break;
      case 'T':
      case 't':
        consensus(3, j)++;
        break;
      }
    }
  }
  return consensus;
}

//' @export
// [[Rcpp::export]]
Rcpp::IntegerVector find_longest_matches(Rcpp::StringVector cigars) {
  std::vector<int> longest_matches;
  longest_matches.reserve(cigars.length());
  
  for(const auto& cigar : cigars) {
    std::string cigar_str = Rcpp::as<std::string>(cigar);
    int current_length = 0;
    int max_length = 0;
    std::string number_buffer;
    
    for(size_t i = 0; i < cigar_str.length(); i++) {
      char c = cigar_str[i];
      
      if(std::isdigit(c)) {
        number_buffer += c;
      } else {
        // Process operation
        if(c == '=') {
          // Found a match operation
          current_length = std::stoi(number_buffer);
          max_length = std::max(max_length, current_length);
        }
        number_buffer.clear();
      }
    }
    
    longest_matches.push_back(max_length);
  }
  
  return Rcpp::wrap(longest_matches);
}

//' @export
// [[Rcpp::export]]
Rcpp::List analyze_cigar_patterns(Rcpp::StringVector cigars) {
  std::vector<int> longest_matches;
  std::vector<int> total_matches;
  std::vector<int> match_blocks;
  longest_matches.reserve(cigars.length());
  total_matches.reserve(cigars.length());
  match_blocks.reserve(cigars.length());
  
  for(const auto& cigar : cigars) {
    std::string cigar_str = Rcpp::as<std::string>(cigar);
    int current_length = 0;
    int max_length = 0;
    int total_match = 0;
    int blocks = 0;
    std::string number_buffer;
    
    for(size_t i = 0; i < cigar_str.length(); i++) {
      char c = cigar_str[i];
      
      if(std::isdigit(c)) {
        number_buffer += c;
      } else {
        int num = std::stoi(number_buffer);
        if(c == '=') {
          // Found a match operation
          current_length = num;
          max_length = std::max(max_length, current_length);
          total_match += num;
          blocks++;
        }
        number_buffer.clear();
      }
    }
    
    longest_matches.push_back(max_length);
    total_matches.push_back(total_match);
    match_blocks.push_back(blocks);
  }
  
  return Rcpp::List::create(
    Rcpp::Named("longest_match") = Rcpp::wrap(longest_matches),
    Rcpp::Named("total_matches") = Rcpp::wrap(total_matches),
    Rcpp::Named("match_blocks") = Rcpp::wrap(match_blocks)
  );
}

//' @export
// [[Rcpp::export]]
bool has_minimum_match_length(const std::string& cigar, int min_length) {
  std::string number_buffer;
  
  for(size_t i = 0; i < cigar.length(); i++) {
    char c = cigar[i];
    
    if(std::isdigit(c)) {
      number_buffer += c;
    } else {
      if(c == '=' && std::stoi(number_buffer) >= min_length) {
        return true;
      }
      number_buffer.clear();
    }
  }
  
  return false;
}

//' @export
// [[Rcpp::export]]
List edlib_search(std::vector<std::string> queries, 
                                 std::vector<std::string> targets,
                                 int max_distance = 3,
                                 int num_threads = 0) {
  
  // Set number of threads (0 = auto-detect)
  if (num_threads <= 0) {
    num_threads = omp_get_max_threads();
  }
  omp_set_num_threads(num_threads);
  
  Rcout << "Starting multithreaded search with " << num_threads << " threads..." << std::endl;
  Rcout << "Processing " << queries.size() << " queries against " 
        << targets.size() << " targets (max distance: " << max_distance << ")" << std::endl;
  
  // Pre-allocate result list
  List results(queries.size());
  CharacterVector query_names(queries.size());
  
  // Progress tracking variables (need to be thread-safe)
  std::vector<int> progress_counter(num_threads, 0);
  int total_matches = 0;
  
  // Parallel loop over queries
#pragma omp parallel for schedule(dynamic)
  for (int q = 0; q < (int)queries.size(); q++) {
    
    int thread_id = omp_get_thread_num();
    progress_counter[thread_id]++;
    
    // Progress reporting (only from thread 0 to avoid race conditions)
    const std::string& query = queries[q];
    const char* query_cstr = query.c_str();
    int query_len = query.length();
    
    std::vector<std::string> matched_targets;
    std::vector<int> matched_distances;
    
    // Pre-create edlib config (thread-local)
    EdlibAlignConfig config = edlibNewAlignConfig(max_distance, 
                                                  EDLIB_MODE_NW, 
                                                  EDLIB_TASK_DISTANCE, 
                                                  NULL, 0);
    
    // Inner loop over targets (not parallelized to avoid nested parallelism)
    for (size_t t = 0; t < targets.size(); t++) {
      
      const std::string& target = targets[t];
      const char* target_cstr = target.c_str();
      int target_len = target.length();
      
      // Quick length-based filteringd
      if (abs(query_len - target_len) > max_distance) {
        continue;
      }
      
      EdlibAlignResult result = edlibAlign(query_cstr, query_len,
                                           target_cstr, target_len,
                                           config);
      if (result.editDistance != -1 && result.editDistance <= max_distance) {
        matched_targets.push_back(target);
        matched_distances.push_back(result.editDistance);
        total_matches++;
      }
      edlibFreeAlignResult(result);
    }
    // Create DataFrame for this query's matches (thread-safe since each thread writes to different index)
    if (matched_targets.size() > 0) {
      DataFrame query_matches = DataFrame::create(Named("target") = matched_targets,
                                                  Named("distance") = matched_distances);
      results[q] = query_matches;
    } else {
      // Empty DataFrame for queries with no matches
      results[q] = DataFrame::create(Named("target") = CharacterVector(0),
                                     Named("distance") = IntegerVector(0));
    }
    query_names[q] = query;
  }
  
  results.names() = query_names;
  Rcout << "Multithreaded search complete!" << std::endl;
  Rcout << "Found " << total_matches << " total matches across all queries." << std::endl;
  return results;
}
