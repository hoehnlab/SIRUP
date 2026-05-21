#include "rad.h"
#include "static_alignment_parity.h"
using namespace std;
using namespace Rcpp;

KSEQ_INIT(gzFile, gzread)

std::string extract_corrected_barcode(const std::string& full_header) {
  size_t start_pos = full_header.find("CB:Z:");
  if (start_pos == std::string::npos) return "";
  
  size_t barcode_start = start_pos + 5; // Skip "CB:Z:"
  
  // Find the next space or tab after the barcode, stopping before any UB:Z: field
  size_t barcode_end = full_header.find_first_of("\t ", barcode_start);
  if (barcode_end == std::string::npos) {
    barcode_end = full_header.length(); // In case there is no space/tab after the barcode
  }
  
  return full_header.substr(barcode_start, barcode_end - barcode_start);
}

std::string extract_blaze_id(const std::string& header) {
  size_t start_pos = header.find("#");
  size_t end_pos = header.find("_", start_pos);
  
  if (start_pos != std::string::npos && end_pos != std::string::npos) {
    return header.substr(start_pos + 1, end_pos - start_pos - 1);
  }
  
  // Return the full header if the pattern is not found
  return header;
}

//' @export
// [[Rcpp::export]]
Rcpp::DataFrame extract_blaze_barcode(std::string fastq_file, int print_every = 10000) {
  gzFile fp = gzopen(fastq_file.c_str(), "r");
  if (!fp) Rcpp::stop("Failed to open file.");
  
  kseq_t *seq = kseq_init(fp);
  int l;
  
  // Hash map to store the counts of corrected barcodes
  std::unordered_map<std::string, int> barcode_counts;
  
  int sequence_count = 0;
  
  // Read through the fastq file with kseq
  while ((l = kseq_read(seq)) >= 0) {
    sequence_count++;
    
    // Capture both the sequence name (header) and the comment (which might contain CB:Z)
    std::string header = "@" + std::string(seq->name.s);  // Add '@' symbol to the header
    std::string full_header = header;
    if (seq->comment.l > 0) {
      full_header += "\t" + std::string(seq->comment.s); // Append the comment section (CB:Z might be here)
    }
    
    // Extract the corrected barcode (CB:Z: field)
    std::string corrected_barcode = extract_corrected_barcode(full_header);
    
    // If the corrected barcode is missing, skip
    if (corrected_barcode.empty()) continue;
    
    // Increment the count for the corrected barcode
    barcode_counts[corrected_barcode]++;
    
    // Print progress every 'print_every' sequences
    if (sequence_count % print_every == 0) {
      std::cout << "Processed " << sequence_count << " sequences so far." << std::endl;
    }
  }
  
  kseq_destroy(seq);
  gzclose(fp);
  
  std::cout << "Total sequences processed: " << sequence_count << std::endl;
  
  // Prepare vectors for the DataFrame
  Rcpp::CharacterVector barcodes;
  Rcpp::IntegerVector total_counts;
  for (const auto& entry : barcode_counts) {
    barcodes.push_back(entry.first);              // Corrected barcode
    total_counts.push_back(entry.second);  // Total count
  }
  
  // Return the DataFrame with barcode and total counts
  return Rcpp::DataFrame::create(Rcpp::Named("barcode") = barcodes,
    Rcpp::Named("total") = total_counts);
}

//' @export
// [[Rcpp::export]]
void extract_blaze_id_bc(std::string fastq_file, std::string output_file, int print_every = 10000) {
  gzFile fp = gzopen(fastq_file.c_str(), "r");
  if (!fp) Rcpp::stop("Failed to open file.");
  
  kseq_t *seq = kseq_init(fp);
  int l;
  
  // Vectors to store the extracted IDs and barcodes
  std::vector<std::string> ids;
  std::vector<std::string> barcodes;
  
  int sequence_count = 0;
  
  // Read through the fastq file with kseq
  while ((l = kseq_read(seq)) >= 0) {
    sequence_count++;
    
    // Capture both the sequence name (header) and the comment (which might contain CB:Z)
    std::string header = "@" + std::string(seq->name.s);  // Add '@' symbol to the header
    std::string full_header = header;
    if (seq->comment.l > 0) {
      full_header += "\t" + std::string(seq->comment.s); // Append the comment section (CB:Z might be here)
    }
    
    // Extract the corrected barcode (CB:Z: field)
    std::string corrected_barcode = extract_corrected_barcode(full_header);
    
    // If the corrected barcode is missing, skip
    if (corrected_barcode.empty()) continue;
    
    // Extract the read ID part between '#' and '_'
    std::string extracted_id = extract_blaze_id(header);
    
    // Add the extracted ID and corrected barcode to the vectors
    ids.push_back(extracted_id);
    barcodes.push_back(corrected_barcode);
    
    // Print progress every 'print_every' sequences
    if (sequence_count % print_every == 0) {
      std::cout << "Processed " << sequence_count << " sequences so far." << std::endl;
    }
  }
  
  kseq_destroy(seq);
  gzclose(fp);
  
  std::cout << "Total sequences processed: " << sequence_count << std::endl;
  
  // Write the results to a CSV file
  std::ofstream out_file(output_file);
  if (!out_file.is_open()) {
    Rcpp::stop("Failed to open output file.");
  }
  
  // Write the header
  out_file << "id,barcode\n";
  
  // Write the data
  for (size_t i = 0; i < ids.size(); ++i) {
    out_file << ids[i] << "," << barcodes[i] << "\n";
  }
  
  out_file.close();
  
  std::cout << "Data successfully written to " << output_file << std::endl;
}

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
Rcpp::CharacterVector masked_sigalign(
    Rcpp::CharacterVector adapters,
    std::vector<std::string> sequences,
    std::vector<std::string> ids,
    const Rcpp::DataFrame& misalignment_threshold, 
    int nthreads) {
  
  std::vector<std::string> queries = Rcpp::as<std::vector<std::string>>(adapters);
  std::vector<std::string> query_names = Rcpp::as<std::vector<std::string>>(adapters.names());
  
  // Create masked versions of adapters
  std::vector<std::pair<std::string, int>> masked_adapters;
  for(const auto& query : queries) {
    // Store original length for position correction
    int orig_length = query.length();
    
    // Create masked version - keeping last 10 bases unmasked
    // Adjust these numbers based on your needs
    int unmask_length = std::min(10, (int)query.length());
    int mask_length = query.length() - unmask_length;
    
    std::string masked_query = std::string(mask_length, 'N') + 
      query.substr(query.length() - unmask_length);
    
    masked_adapters.push_back({masked_query, orig_length});
  }
  
  std::map<int, std::string> signature_map;
  std::map<std::string, std::pair<double, double>> null_dist_map;
  Rcpp::StringVector query_id = misalignment_threshold["query_id"];
  Rcpp::NumericVector misal_threshold = misalignment_threshold["misal_threshold"];
  Rcpp::NumericVector misal_sd = misalignment_threshold["misal_sd"];
  
  for(int i = 0; i < query_id.size(); ++i) {
    null_dist_map[Rcpp::as<std::string>(query_id[i])] = 
      std::make_pair((double) misal_threshold[i], (double) misal_sd[i]);
  }
  
  omp_set_num_threads(nthreads);
#pragma omp parallel for schedule(dynamic)
  for (int i = 0; i < sequences.size(); ++i) {
    const auto& sequence = sequences[i];
    std::vector<std::string> temp_signature_parts;
    int length = sequence.size();
    
    for (int j = 0; j < queries.size(); ++j) {
      const auto& query_name = query_names[j];
      const auto& masked_adapter = masked_adapters[j];
      
      if (query_name == "poly_a" || query_name == "poly_t") {
        char poly_base = (query_name == "poly_a") ? 'A' : 'T';
        std::string match_str = findPolyTails(sequence, poly_base, 14, 12);
        if (!match_str.empty()) {
#pragma omp critical
          temp_signature_parts.push_back(match_str);
        }
      } else {
        // Configure edlib to treat N as equivalent to any base
        EdlibEqualityPair additionalEqualities[4] = {
          {'N', 'A'}, {'N', 'C'}, {'N', 'G'}, {'N', 'T'}
        };
        
        EdlibAlignConfig config = edlibNewAlignConfig(
          -1,                    // no max edit distance
          EDLIB_MODE_HW,        // semi-global alignment
          EDLIB_TASK_LOC,       // location tasks
          additionalEqualities,  // N matches any base
          4                      // number of equality pairs
        );
        
        char* cquery = const_cast<char*>(masked_adapter.first.c_str());
        char* csequence = const_cast<char*>(sequence.c_str());
        
        EdlibAlignResult cresult = edlibAlign(
          cquery, 
          masked_adapter.first.size(), 
          csequence, 
          sequence.size(), 
          config
        );
        
        // Process alignment results with position correction
        std::vector<int> start_positions(
            cresult.startLocations, 
            cresult.startLocations + cresult.numLocations
        );
        std::vector<int> end_positions(
            cresult.endLocations, 
            cresult.endLocations + cresult.numLocations
        );
        
        // Adjust positions based on original adapter length
        for (int& pos : start_positions) {
          pos = pos + 1 - (masked_adapter.first.length() - masked_adapter.second);
        }
        for (int& pos : end_positions) { ++pos; }
        
        // Rest of your existing alignment processing code...
        std::sort(start_positions.begin(), start_positions.end());
        std::sort(end_positions.begin(), end_positions.end());
        
        auto last = std::unique(start_positions.begin(), start_positions.end());
        start_positions.erase(last, start_positions.end());
        
        std::vector<int> unique_end_positions;
        for (const auto& start : start_positions) {
          auto pos = std::find(start_positions.begin(), start_positions.end(), start) - 
            start_positions.begin();
          unique_end_positions.push_back(end_positions[pos]);
        }
        
        int edit_distance = cresult.editDistance;
        std::set<UniqueAlignment> uniqueAlignments;
        
        for (size_t k = 0; k < start_positions.size(); ++k) {
          UniqueAlignment ua = {
            edit_distance, 
            start_positions[k], 
            unique_end_positions[k]
          };
          uniqueAlignments.insert(ua);
        }
        
        // Extract best alignments and process results
        auto it = uniqueAlignments.begin();
        UniqueAlignment best = (it != uniqueAlignments.end()) ? 
        *it : UniqueAlignment{-1, -1, -1}; 
        ++it;
        UniqueAlignment secondBest = (it != uniqueAlignments.end()) ? 
        *it : UniqueAlignment{-1, -1, -1};
        
        bool skip_second_best = false;
        auto null_data = null_dist_map.find(query_name);
        if (null_data != null_dist_map.end()) {
          double threshold = null_data->second.first - null_data->second.second;
          if (secondBest.edit_distance >= threshold || secondBest.edit_distance == -1) {
            skip_second_best = true;
          }
        }
        
#pragma omp critical
{
  std::string best_signature_part = query_name + ":" + 
    std::to_string(best.edit_distance) + ":" + 
    std::to_string(best.start_position) + ":" + 
    std::to_string(best.end_position);
  temp_signature_parts.push_back(best_signature_part);
  
  if (!skip_second_best) {
    std::string second_best_signature_part = query_name + ":" + 
      std::to_string(secondBest.edit_distance) + ":" + 
      std::to_string(secondBest.start_position) + ":" + 
      std::to_string(secondBest.end_position);
    temp_signature_parts.push_back(second_best_signature_part);
  }
}
edlibFreeAlignResult(cresult);
      }
    }
    
    // Add your hardcoded signature parts
    temp_signature_parts.push_back("seq_start:0:1:1");
    temp_signature_parts.push_back("rc_seq_start:0:1:1");
    temp_signature_parts.push_back("seq_stop:0:" + std::to_string(length) + ":" + 
      std::to_string(length));
    temp_signature_parts.push_back("rc_seq_stop:0:" + std::to_string(length) + ":" + 
      std::to_string(length));
    
    // Sort and join signature parts
    std::sort(temp_signature_parts.begin(), temp_signature_parts.end(),
      [](const std::string &a, const std::string &b) -> bool {
        int start_pos_a = std::stoi(a.substr(a.find_last_of(":") + 1));
        int start_pos_b = std::stoi(b.substr(b.find_last_of(":") + 1));
        return start_pos_a < start_pos_b;
      }
    );
    
    std::string final_signature = std::accumulate(
      temp_signature_parts.begin(), 
      temp_signature_parts.end(),
      std::string(),
      [](const std::string& a, const std::string& b) -> string {
        return a + (a.length() > 0 ? "|" : "") + b;
      }
    );
    
#pragma omp critical
{
  final_signature += "<" + std::to_string(length) + ":" + ids[i] + ":undecided>";
  signature_map[i + 1] = final_signature;
  temp_signature_parts.clear();
}
  }
  
  Rcpp::CharacterVector signature_strings(signature_map.size());
  for (const auto& pair : signature_map) {
    signature_strings[pair.first - 1] = pair.second;
  }
  return signature_strings;
}
//' @export
// [[Rcpp::export]]
Rcpp::DataFrame n_masked_misalignment_cigar_string(
    const std::string& input_path,
    Rcpp::CharacterVector adapters,
    int nthreads = 1,
    int max_sequences = 500000,
    size_t chunk_size = 50000) {
  std::vector<std::string> input_files = get_fastq_files(input_path);
  if (input_files.empty()) {
    Rcpp::warning("No FASTQ files found in the specified path.");
    return Rcpp::DataFrame::create();
  }
  std::vector<std::string> queries = Rcpp::as<std::vector<std::string>>(adapters);
  Rcpp::CharacterVector query_names = adapters.names();
  
  std::vector<std::string> filtered_queries;
  Rcpp::CharacterVector filtered_query_names;
  for (int i = 0; i < query_names.size(); ++i) {
    std::string name = Rcpp::as<std::string>(query_names[i]);
    if (name.find("poly") == std::string::npos) {
      filtered_queries.push_back(queries[i]);
      filtered_query_names.push_back(query_names[i]);
    }
  }
  
  std::vector<int> ids;
  std::vector<std::string> query_ids;
  std::vector<int> best_edit_distances;
  std::vector<int> best_start_positions;
  std::vector<int> best_stop_positions;
  std::vector<int> second_edit_distances;
  std::vector<std::string> aligned_sequences;
  std::vector<int> sequence_lengths;
  std::vector<std::string> cigar_strings;
  
  Rcpp::Rcout << "Generating misalignment threshold data!\n";
  
  omp_set_num_threads(nthreads);
  
  auto process_chunk = [&](const std::vector<std::string>& chunk_sequences, const std::vector<std::string>& chunk_ids, int sequence_counter, int remaining) {
#pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < std::min(static_cast<int>(chunk_sequences.size()), remaining); ++i) {
      const auto& sequence = chunk_sequences[i];
      std::vector<AlignmentInfo> alignments(filtered_queries.size());
      std::vector<std::string> local_aligned_seqs(filtered_queries.size());
      std::vector<std::string> local_cigars(filtered_queries.size());
      
      for (size_t j = 0; j < filtered_queries.size(); ++j) {
        const auto& query = filtered_queries[j];
        
        EdlibEqualityPair additionalEqualities[4] = {
          {'N', 'A'}, {'N', 'C'}, {'N', 'G'}, {'N', 'T'}
        };
        
        EdlibAlignConfig config = edlibNewAlignConfig(
          -1,
          EDLIB_MODE_HW,
          EDLIB_TASK_PATH,
          additionalEqualities,
          4
        );
        
        EdlibAlignResult result = edlibAlign(query.c_str(), query.size(), sequence.c_str(), sequence.size(), config);
        
        alignments[j] = {
          result.editDistance,
          result.startLocations[0] + 1,
          result.endLocations[0] + 1,
          (result.numLocations > 1) ? result.editDistance : -1
        };
        
        int start_pos = result.startLocations[0];
        int end_pos = result.endLocations[0] + 1;
        if (start_pos >= 0 && end_pos <= sequence.length()) {
          local_aligned_seqs[j] = sequence.substr(start_pos, end_pos - start_pos);
          
          // Generate CIGAR string using edlib's function
          if (result.alignment != nullptr) {
            char* cigar = edlibAlignmentToCigar(result.alignment, result.alignmentLength, EDLIB_CIGAR_EXTENDED);
            if (cigar != nullptr) {
              local_cigars[j] = std::string(cigar);
              free(cigar);
            } else {
              local_cigars[j] = "";
            }
          } else {
            local_cigars[j] = "";
          }
        } else {
          local_aligned_seqs[j] = "";
          local_cigars[j] = "";
        }
        
        edlibFreeAlignResult(result);
      }
      
#pragma omp critical
{
  for (size_t j = 0; j < alignments.size(); ++j) {
    ids.push_back(sequence_counter + i);
    query_ids.push_back(Rcpp::as<std::string>(filtered_query_names[j]));
    best_edit_distances.push_back(alignments[j].best_edit_distance);
    best_start_positions.push_back(alignments[j].best_start_pos);
    best_stop_positions.push_back(alignments[j].best_stop_pos);
    second_edit_distances.push_back(alignments[j].second_edit_distance);
    aligned_sequences.push_back(local_aligned_seqs[j]);
    sequence_lengths.push_back(sequence.length());
    cigar_strings.push_back(local_cigars[j]);
  }
}
    }
  };
  
  int total_sequences = 0;
  for (const auto& file_path : input_files) {
    gzFile file = gzopen(file_path.c_str(), "rb");
    if (!file) {
      Rcpp::warning("Failed to open file: " + file_path);
      continue;
    }
    kseq_t* seq = kseq_init(file);
    std::vector<std::string> chunk_sequences;
    std::vector<std::string> chunk_ids;
    while (kseq_read(seq) >= 0 && (max_sequences == -1 || total_sequences < max_sequences)) {
      chunk_sequences.push_back(seq->seq.s);
      chunk_ids.push_back(seq->name.s);
      if (chunk_sequences.size() >= chunk_size || (max_sequences != -1 && total_sequences + chunk_sequences.size() >= max_sequences)) {
        int remaining = (max_sequences == -1) ? chunk_sequences.size() : std::min(static_cast<int>(chunk_sequences.size()), max_sequences - total_sequences);
        process_chunk(chunk_sequences, chunk_ids, total_sequences, remaining);
        total_sequences += chunk_sequences.size();
        chunk_sequences.clear();
        chunk_ids.clear();
        if (max_sequences != -1 && total_sequences >= max_sequences) break;
      }
    }
    
    if (!chunk_sequences.empty() && (max_sequences == -1 || total_sequences < max_sequences)) {
      int remaining = (max_sequences == -1) ? chunk_sequences.size() : std::min(static_cast<int>(chunk_sequences.size()), max_sequences - total_sequences);
      process_chunk(chunk_sequences, chunk_ids, total_sequences, remaining);
      total_sequences += chunk_sequences.size();
    }
    
    kseq_destroy(seq);
    gzclose(file);
    
    if (max_sequences != -1 && total_sequences >= max_sequences) break;
  }
  
  if (total_sequences == 0) {
    Rcpp::warning("No sequences were processed.");
    return Rcpp::DataFrame::create();
  }
  
  return Rcpp::DataFrame::create(
    Rcpp::Named("id") = ids,
    Rcpp::Named("query_id") = query_ids,
    Rcpp::Named("best_edit_distance") = best_edit_distances,
    Rcpp::Named("best_start_pos") = best_start_positions,
    Rcpp::Named("best_stop_pos") = best_stop_positions,
    Rcpp::Named("second_edit_distance") = second_edit_distances,
    Rcpp::Named("aligned_sequence") = aligned_sequences,
    Rcpp::Named("sequence_length") = sequence_lengths,
    Rcpp::Named("cigar") = cigar_strings
  );
}
//' @export
// [[Rcpp::export]]
Rcpp::CharacterVector n_masked_sigalign(
    Rcpp::CharacterVector adapters,
    std::vector<std::string> sequences,
    std::vector<std::string> ids,
    const Rcpp::DataFrame& misalignment_threshold,
    int nthreads,
    bool verbose = false) {

  std::vector<std::string> queries = Rcpp::as<std::vector<std::string>>(adapters);
  std::vector<std::string> query_names = Rcpp::as<std::vector<std::string>>(adapters.names());

  Rcpp::StringVector query_id = misalignment_threshold["query_id"];
  Rcpp::NumericVector misal_threshold = misalignment_threshold["misal_threshold"];
  Rcpp::NumericVector misal_sd = misalignment_threshold["misal_sd"];
  Rcpp::NumericVector align_mean_start = misalignment_threshold["align_mean_start"];
  Rcpp::NumericVector align_mean_stop = misalignment_threshold["align_mean_stop"];
  Rcpp::NumericVector align_sd_start = misalignment_threshold["align_sd_start"];
  Rcpp::NumericVector align_sd_stop = misalignment_threshold["align_sd_stop"];

  struct AlignmentBounds {
    double misal_threshold;
    double misal_sd;
    double start_mean;
    double stop_mean;
    double start_sd;
    double stop_sd;
  };

  std::map<std::string, AlignmentBounds> alignment_bounds;
  for (int i = 0; i < query_id.size(); ++i) {
    alignment_bounds[Rcpp::as<std::string>(query_id[i])] = {
      misal_threshold[i],
      misal_sd[i],
      align_mean_start[i],
      align_mean_stop[i],
      align_sd_start[i],
      align_sd_stop[i]
    };
  }

  size_t total_adapter_length = 0;
  for (size_t j = 0; j < query_names.size(); ++j) {
    if (query_names[j] != "poly_a" && query_names[j] != "poly_t") {
      total_adapter_length += queries[j].length();
    }
  }

  std::map<int, std::string> signature_map;

  omp_set_num_threads(nthreads);
#pragma omp parallel for schedule(dynamic)
  for (int i = 0; i < static_cast<int>(sequences.size()); ++i) {
    const std::string& original_sequence = sequences[i];
    int length = static_cast<int>(original_sequence.size());

    if (length < static_cast<int>(total_adapter_length)) {
#pragma omp critical
      {
        signature_map[i + 1] = "<" + std::to_string(length) + ":" + ids[i] + ":undecided>";
      }
      continue;
    }

    std::string mutable_sequence = original_sequence;
    std::vector<std::string> temp_signature_parts;

    for (int j = 0; j < static_cast<int>(queries.size()); ++j) {
      const std::string& query = queries[j];
      const std::string& query_name = query_names[j];

      if (query_name == "poly_a" || query_name == "poly_t") {
        const char poly_base = (query_name == "poly_a") ? 'A' : 'T';
        const std::string match_str = findPolyTails(mutable_sequence, poly_base, 14, 12);
        if (!match_str.empty()) {
          const std::vector<std::string> poly_parts = rad_static::split_sig_parts(match_str);
          for (const auto& part : poly_parts) {
            temp_signature_parts.push_back(part);
            int start_pos = -1;
            int stop_pos = -1;
            if (rad_static::parse_sig_span(part, start_pos, stop_pos)) {
              rad_static::mask_span(mutable_sequence, start_pos, stop_pos);
            }
          }
        }
        continue;
      }

      auto bounds_it = alignment_bounds.find(query_name);
      if (bounds_it == alignment_bounds.end()) {
        continue;
      }

      const AlignmentBounds& bounds = bounds_it->second;
      const int max_edit_distance = std::max(0, static_cast<int>(std::floor(bounds.misal_threshold)));

      int expected_start = 1;
      int expected_end = length;
      if (bounds.start_mean < 50.0) {
        expected_end = std::min(length, static_cast<int>(((bounds.start_mean + (bounds.start_sd * 1.25)) / 100.0) * length));
      } else {
        expected_start = std::max(1, static_cast<int>(((bounds.start_mean - (bounds.start_sd * 1.25)) / 100.0) * length));
      }

      rad_static::StaticAlignment alignment = rad_static::align_static_element(
        query,
        mutable_sequence,
        verbose,
        max_edit_distance,
        expected_start,
        expected_end
      );

      if (!alignment.success) {
        continue;
      }

      for (const auto& pos : alignment.positions) {
        const int start_pos = std::max(1, pos.first);
        const int stop_pos = std::min(length, pos.second);
        temp_signature_parts.push_back(
          query_name + ":" + std::to_string(alignment.edit_distance) + ":" +
            std::to_string(start_pos) + ":" + std::to_string(stop_pos)
        );
        rad_static::mask_span(mutable_sequence, start_pos, stop_pos);
      }
    }

    temp_signature_parts.push_back("seq_start:0:1:1");
    temp_signature_parts.push_back("rc_seq_start:0:1:1");
    temp_signature_parts.push_back("seq_stop:0:" + std::to_string(length) + ":" + std::to_string(length));
    temp_signature_parts.push_back("rc_seq_stop:0:" + std::to_string(length) + ":" + std::to_string(length));

    std::sort(temp_signature_parts.begin(), temp_signature_parts.end(),
      [](const std::string& a, const std::string& b) {
        int start_pos_a = std::stoi(a.substr(a.find_last_of(":") + 1));
        int start_pos_b = std::stoi(b.substr(b.find_last_of(":") + 1));
        return start_pos_a < start_pos_b;
      }
    );

    const std::string final_signature = std::accumulate(
      temp_signature_parts.begin(),
      temp_signature_parts.end(),
      std::string(),
      [](const std::string& left, const std::string& right) {
        return left + (left.empty() ? "" : "|") + right;
      }
    );

#pragma omp critical
    {
      signature_map[i + 1] = final_signature + "<" + std::to_string(length) + ":" + ids[i] + ":undecided>";
    }
  }

  Rcpp::CharacterVector signature_strings(signature_map.size());
  for (const auto& pair : signature_map) {
    signature_strings[pair.first - 1] = pair.second;
  }
  return signature_strings;
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
