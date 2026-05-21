#include "rad.h"
#include "static_alignment_parity.h"
using namespace std;
using namespace Rcpp;
namespace fs = std::filesystem;

KSEQ_INIT(gzFile, gzread)
  
  std::mutex fasta_mutex;
  std::mutex sig_map_mutex;

struct BarcodeCount {
    std::string barcode;
    int count;
  };
  
struct SearchResult {
  std::string read_id;
  std::string sequence;
  std::string qual_scores;
};

struct AlignmentBounds {
  double misal_threshold;
  double misal_sd;
  double start_mean;
  double stop_mean;
  double start_sd;
  double stop_sd;
};

bool ends_with(const std::string& str, const std::string& suffix) {
  return str.size() >= suffix.size() && str.compare(str.size() - suffix.size(), suffix.size(), suffix) == 0;
}

bool is_regular_file(const std::string& path) {
  struct stat path_stat;
  if (stat(path.c_str(), &path_stat) != 0) {
    return false;
  }
  return S_ISREG(path_stat.st_mode);
}

bool is_directory(const std::string& path) {
  struct stat path_stat;
  if (stat(path.c_str(), &path_stat) != 0) {
    return false;
  }
  return S_ISDIR(path_stat.st_mode);
}

std::string format_duration(std::chrono::seconds seconds) {
  auto hours = std::chrono::duration_cast<std::chrono::hours>(seconds);
  seconds -= std::chrono::duration_cast<std::chrono::seconds>(hours);
  auto minutes = std::chrono::duration_cast<std::chrono::minutes>(seconds);
  seconds -= std::chrono::duration_cast<std::chrono::seconds>(minutes);
  
  std::stringstream ss;
  ss << std::setfill('0') << std::setw(2) << hours.count() << ":"
     << std::setfill('0') << std::setw(2) << minutes.count() << ":"
     << std::setfill('0') << std::setw(2) << seconds.count();
  return ss.str();
}

std::vector<SearchResult> readFastqChunk(
    kseq_t* seq, 
    gzFile fp, 
    size_t chunkSize, 
    bool& is_fastq, 
    int& read_count, 
    int max_sequences) {
  std::vector<SearchResult> records;
  for (size_t i = 0; i < chunkSize && kseq_read(seq) >= 0; ++i) {
    if (max_sequences != -1 && read_count >= max_sequences) break;
    SearchResult sr;
    sr.read_id = seq->name.s;
    sr.sequence = seq->seq.s;
    if (seq->qual.l) {
      sr.qual_scores = seq->qual.s;
      is_fastq = true;
    } else {
      is_fastq = false;
    }
    records.push_back(sr);
    ++read_count;
  }
  return records;
}

std::vector<std::string> get_fastq_files(const std::string& path) {
  std::vector<std::string> fastq_files;
  
  if (is_regular_file(path)) {
    fastq_files.push_back(path);
  } else if (is_directory(path)) {
    DIR* dir = opendir(path.c_str());
    if (dir == nullptr) {
      Rcpp::stop("Could not open directory: " + path);
    }
    
    struct dirent* entry;
    while ((entry = readdir(dir)) != nullptr) {
      std::string filename = entry->d_name;
      if (ends_with(filename, ".fastq.gz") || ends_with(filename, ".fq.gz")) {
        fastq_files.push_back(path + "/" + filename);
      }
    }
    
    closedir(dir);
  } else {
    Rcpp::stop("Invalid path: " + path);
  }
  
  return fastq_files;
}

std::pair<std::vector<std::string>,
          std::vector<std::string>> parseAdapters(
              const Rcpp::CharacterVector& adapters) {
            std::vector<std::string> queries(adapters.size());
            std::vector<std::string> query_names(adapters.size());
            Rcpp::CharacterVector names = adapters.names();
            for (int i = 0; i < adapters.size(); ++i) {
              queries[i] = Rcpp::as<std::string>(adapters[i]);
              query_names[i] = Rcpp::as<std::string>(names[i]);
            }
            return std::make_pair(queries, query_names);
          }

std::map<std::string, AlignmentBounds> parseThresholds(
    const Rcpp::DataFrame& misalignment_threshold) {
  std::map<std::string, AlignmentBounds> alignment_bounds;
  
  Rcpp::StringVector query_id = misalignment_threshold["query_id"];
  Rcpp::NumericVector misal_threshold = misalignment_threshold["misal_threshold"];
  Rcpp::NumericVector misal_sd = misalignment_threshold["misal_sd"];
  Rcpp::NumericVector align_mean_start = misalignment_threshold["align_mean_start"];
  Rcpp::NumericVector align_mean_stop = misalignment_threshold["align_mean_stop"];
  Rcpp::NumericVector align_sd_start = misalignment_threshold["align_sd_start"];
  Rcpp::NumericVector align_sd_stop = misalignment_threshold["align_sd_stop"];
  
  for(int i = 0; i < query_id.size(); ++i) {
    std::string id = Rcpp::as<std::string>(query_id[i]);
    alignment_bounds[id] = {
      misal_threshold[i],
                     misal_sd[i],
                             align_mean_start[i],
                                             align_mean_stop[i],
                                                            align_sd_start[i],
                                                                          align_sd_stop[i]
    };
  }
  return alignment_bounds;
}

std::vector<std::string> n_masked_sigalign(
    const std::vector<std::string>& queries,
    const std::vector<std::string>& query_names,
    const std::vector<std::string>& sequences,
    const std::vector<std::string>& ids,
    const std::map<std::string, AlignmentBounds>& alignment_bounds,
    int nthreads,
    bool verbose = false) {

  std::map<int, std::string> signature_map;
  std::mutex local_sig_map_mutex;

  size_t total_adapter_length = 0;
  for (size_t j = 0; j < query_names.size(); ++j) {
    if (query_names[j] != "poly_a" && query_names[j] != "poly_t") {
      total_adapter_length += queries[j].length();
    }
  }

  omp_set_num_threads(nthreads);
#pragma omp parallel for schedule(dynamic)
  for (int i = 0; i < static_cast<int>(sequences.size()); ++i) {
    const std::string& original_sequence = sequences[i];
    const int length = static_cast<int>(original_sequence.size());

    if (length < static_cast<int>(total_adapter_length)) {
      std::lock_guard<std::mutex> lock(local_sig_map_mutex);
      signature_map[i + 1] = "<" + std::to_string(length) + ":" + ids[i] + ":undecided>";
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
              });

    const std::string final_signature = std::accumulate(
      temp_signature_parts.begin(),
      temp_signature_parts.end(),
      std::string(),
      [](const std::string& left, const std::string& right) {
        return left + (left.empty() ? "" : "|") + right;
      }
    );

    std::lock_guard<std::mutex> lock(local_sig_map_mutex);
    signature_map[i + 1] = final_signature + "<" + std::to_string(length) + ":" + ids[i] + ":undecided>";
  }

  std::vector<std::string> signature_strings(signature_map.size());
  for (const auto& pair : signature_map) {
    signature_strings[pair.first - 1] = pair.second;
  }
  return signature_strings;
}

std::vector<std::string> sigrun_cpp(
    const std::vector<std::string>& sigstrings,
    const Rcpp::DataFrame& read_layout,
    const Rcpp::DataFrame& misalignment_threshold,
    int nthreads = 1,
    bool verbose = false) {
  container = prep_read_layout_cpp(read_layout, misalignment_threshold, verbose);
  VarScan(container, verbose);
  PositionFuncMap positionFuncMap = createPositionFunctionMap(container, verbose);
  return process_sigstrings_cpp(container, sigstrings, positionFuncMap, verbose, nthreads);
}

void writeSequenceFile(
    const std::string& fileName, 
    const std::vector<std::pair<std::string, 
                                std::pair<std::string, std::string>>>& sequences, 
                                bool write_fastq, 
                                bool compress) {
  std::string actualFileName = fileName + (compress ? ".gz" : "");
  gzFile compressedFile = NULL;
  std::ofstream normalFile;
  
  if (compress) {
    compressedFile = gzopen(actualFileName.c_str(), "ab");  // Open in append mode
    if (!compressedFile) {
      Rcpp::stop("Failed to open compressed output file.");
    }
  } else {
    normalFile.open(actualFileName, std::ios::out | std::ios::app | std::ios::binary);  // Open in append mode
    if (!normalFile.is_open()) {
      Rcpp::stop("Failed to open output file.");
    }
  }
  
  for (const auto& seq : sequences) {
    if (!seq.second.first.empty()) {
      std::string output;
      if (write_fastq) {
        output = "@" + seq.first + "\n" + seq.second.first + "\n+\n" + seq.second.second + "\n";
      } else {
        output = ">" + seq.first + "\n" + seq.second.first + "\n";
      }
      
      if (compress) {
        gzwrite(compressedFile, output.c_str(), output.length());
      } else {
        normalFile << output;
      }
    }
  }
  
  if (compress) {
    gzclose(compressedFile);
  } else {
    normalFile.close();
  }
}

void writeSigstrings(
    const std::string& sigstringsFilePath,
    const std::vector<std::string>& sigstrings,
    bool compress = true) {
  std::string actualFileName = sigstringsFilePath + (compress ? ".gz" : "");
  gzFile compressedFile = NULL;
  std::ofstream normalFile;
  
  if (compress) {
    compressedFile = gzopen(actualFileName.c_str(), "ab");  // Always open in append mode
    if (!compressedFile) {
      throw std::runtime_error("Failed to open compressed sigstrings file.");
    }
  } else {
    normalFile.open(actualFileName, std::ios::out | std::ios::app | std::ios::binary);  // Always open in append mode
    if (!normalFile.is_open()) {
      throw std::runtime_error("Failed to open sigstrings file.");
    }
  }
  
  for (const auto& sig : sigstrings) {
    std::string output = sig + "\n";
    if (compress) {
      gzwrite(compressedFile, output.c_str(), output.length());
    } else {
      normalFile << output;
    }
  }
  
  if (compress) {
    gzclose(compressedFile);
  } else {
    normalFile.close();
  }
}

void sig_extraction(
    const Rcpp::DataFrame& read_layout,
    const Rcpp::DataFrame& misalignment_threshold,
    const std::vector<std::string>& original_ids,
    const std::vector<std::string>& original_sequences,
    const std::vector<std::string>& original_quality_scores,
    std::vector<std::string>& processed_sigstrings,
    const std::string& outputFilePath,
    bool write_fastq,
    bool compress,
    int min_length,
    int nthreads,
    bool verbose) {
  
  // Initial setup
  ReadLayout readLayout = prep_read_layout_cpp(read_layout, misalignment_threshold, verbose);
  VarScan(readLayout, verbose);
  PositionFuncMap positionFuncMap = createPositionFunctionMap(readLayout, verbose);
  
  // Create sequence map (similar to Rcpp version)
  std::unordered_map<std::string, std::vector<std::string>> seq_map;
  for (size_t i = 0; i < original_ids.size(); ++i) {
    seq_map[original_ids[i]].push_back(original_sequences[i]);
  }
  
  // Extract all unique class IDs
  std::set<std::string> uniqueClassIds;
  auto& class_id_index = readLayout.get<id_tag>();
  for (const auto& elem : class_id_index) {
    if (elem.type == "variable") {
      if (elem.class_id.substr(0, 3) == "rc_") {
        uniqueClassIds.insert(elem.class_id.substr(3));
      } else {
        uniqueClassIds.insert(elem.class_id);
      }
    }
  }
  
  // Store results
  std::vector<std::pair<std::string, std::pair<std::string, std::string>>> result_sequences;
  std::mutex result_mutex;
  
#pragma omp parallel for num_threads(nthreads)
  for (size_t i = 0; i < processed_sigstrings.size(); ++i) {
    std::string signature = processed_sigstrings[i];
    std::map<std::string, std::string> extracted_elements;
    std::string read_sequence, read_quality;
    
    // Parse signature info
    auto read_info_pos = signature.find_last_of('<');
    std::string lengthTypeStr = signature.substr(read_info_pos + 1, signature.length() - read_info_pos - 2);
    std::stringstream lengthTypeStream(lengthTypeStr);
    std::string lengthStr, read_id, type;
    std::getline(lengthTypeStream, lengthStr, ':');
    std::getline(lengthTypeStream, read_id, ':');
    std::getline(lengthTypeStream, type);
    
    if (type == "undecided") continue;
    
    // Remove concatenation suffix if present
    auto plus_pos = read_id.find('+');
    if (plus_pos != std::string::npos) {
      read_id = read_id.substr(0, plus_pos);
    }
    
    // Process each element in the signature
    std::stringstream ss(signature.substr(0, read_info_pos));
    std::string token;
    while (std::getline(ss, token, '|')) {
      std::stringstream tokenStream(token);
      std::string element_id;
      int edit_distance, start_pos, end_pos;
      
      std::getline(tokenStream, element_id, ':');
      tokenStream >> edit_distance;
      tokenStream.ignore(1);
      tokenStream >> start_pos;
      tokenStream.ignore(1);
      tokenStream >> end_pos;
      
      if (!tokenStream) continue;
      auto it = class_id_index.find(element_id);
      if (it != class_id_index.end() && it->type == "variable") {
        auto seq_it = seq_map.find(read_id);
        if (seq_it == seq_map.end()) continue;
        
        for (const auto& seq : seq_it->second) {
          if (start_pos > 0 && end_pos >= start_pos && 
              static_cast<size_t>(end_pos) <= seq.length()) {
            try {
              std::string extracted_seq = seq.substr(
                static_cast<size_t>(start_pos - 1),
                static_cast<size_t>(end_pos - start_pos + 1)
              );
              
              if (it->direction != "forward") {
                extracted_seq = revcomp_cpp(extracted_seq);
              }
              
              std::string classId = element_id;
              if (classId.substr(0, 3) == "rc_") {
                classId = classId.substr(3);
              }
              
              if (it->global_class == "read") {
                read_sequence = extracted_seq;
                // Get corresponding quality scores
                auto qual_it = std::find(original_ids.begin(), original_ids.end(), read_id);
                if (qual_it != original_ids.end()) {
                  size_t idx = qual_it - original_ids.begin();
                  read_quality = original_quality_scores[idx].substr(
                    static_cast<size_t>(start_pos - 1),
                    static_cast<size_t>(end_pos - start_pos + 1)
                  );
                  if (it->direction != "forward") {
                    std::reverse(read_quality.begin(), read_quality.end());
                  }
                }
              } else {
                extracted_elements[classId] = extracted_seq;
              }
            } catch (const std::exception& e) {
              if (verbose) {
#pragma omp critical
{
  std::cout << "Error extracting sequence for " << element_id 
            << ": " << e.what() << std::endl;
}
              }
              continue;
            }
          }
        }
      }
    }
    
    // Check if we have all required elements
    bool has_barcode = extracted_elements.find("barcode") != extracted_elements.end();
    bool has_umi = extracted_elements.find("umi") != extracted_elements.end();
    bool has_read = !read_sequence.empty();
    
    if (has_barcode && has_umi && has_read && 
        read_sequence.length() >= min_length) {
      
      // Construct new ID with extracted elements
      std::string new_id = read_id;
      for (const auto& elem : extracted_elements) {
        new_id += "|" + elem.first + ":" + elem.second;
      }
      
      std::lock_guard<std::mutex> lock(result_mutex);
      result_sequences.emplace_back(
        new_id,
        std::make_pair(read_sequence, read_quality)
      );
      
      if (verbose) {
#pragma omp critical
{
  std::cout << "Accepted sequence " << new_id 
            << " with length " << read_sequence.length() << std::endl;
}
      }
    } else {
      if (verbose) {
#pragma omp critical
{
  std::cout << "Rejected sequence " << read_id 
            << ": barcode=" << has_barcode 
            << ", umi=" << has_umi 
            << ", read=" << has_read 
            << ", length=" << read_sequence.length() << std::endl;
}
      }
    }
  }
  // Sort results by ID to maintain consistency
  std::sort(result_sequences.begin(), result_sequences.end(),
            [](const auto& a, const auto& b) {
              return a.first < b.first;
            });
  
  if (verbose) {
    std::cout << "Total sequences processed: " << processed_sigstrings.size() << std::endl;
    std::cout << "Sequences passing all filters: " << result_sequences.size() << std::endl;
  }
  
  // Write output files
  if (!result_sequences.empty()) {
    if (write_fastq) {
      std::string output_mode = compress ? "wb" : "w";
      gzFile outfile = gzopen(outputFilePath.c_str(), output_mode.c_str());
      
      if (outfile != nullptr) {
        for (const auto& result : result_sequences) {
          std::string fastq_entry = "@" + result.first + "\n" +
            result.second.first + "\n" +
            "+\n" +
            result.second.second + "\n";
          gzwrite(outfile, fastq_entry.c_str(), fastq_entry.length());
        }
        gzclose(outfile);
      } else {
        std::cerr << "Error: Could not open output file " << outputFilePath << std::endl;
      }
    } else {
      std::ofstream outfile(outputFilePath);
      if (outfile.is_open()) {
        for (const auto& result : result_sequences) {
          outfile << ">" << result.first << "\n" << result.second.first << "\n";
        }
        outfile.close();
      } else {
        std::cerr << "Error: Could not open output file " << outputFilePath << std::endl;
      }
    }
  } else {
    if (verbose) {
      std::cout << "Warning: No sequences passed all filters" << std::endl;
    }
    // Create empty output file to indicate processing completed
    std::ofstream outfile(outputFilePath);
    outfile.close();
  }
}

//' @export
// [[Rcpp::export]]
void sigstream(
    const std::string& input_path,
    const std::string& sigstringsFilePath,
    Rcpp::CharacterVector adapters,
    const Rcpp::DataFrame& read_layout,
    const Rcpp::DataFrame& misalignment_threshold,
    int chunkSize, int nthreads,
    const std::string& outputPath,
    bool write_fastq = true,
    bool compress = true,
    int max_sequences = -1,
    int min_length = 100,
    bool verbose = false) {
  auto start_time = std::chrono::high_resolution_clock::now();
  
  std::cout << "Starting sigstream processing..." << std::endl;
  std::vector<std::string> input_files = get_fastq_files(input_path);
  if (input_files.empty()) {
    Rcpp::warning("No FASTQ files found in the specified path.");
    return;
  }
  std::cout << "Found " << input_files.size() << " input files." << std::endl;
  auto parsed_adapters = parseAdapters(adapters);
  std::vector<std::string> queries = parsed_adapters.first;
  std::vector<std::string> query_names = parsed_adapters.second;
  
  std::cout << "Parsing thresholds..." << std::endl;
  auto alignment_bounds = parseThresholds(misalignment_threshold);

  std::cout << "Preparing read layout..." << std::endl;
  ReadLayout container = prep_read_layout_cpp(read_layout, misalignment_threshold, true);
  VarScan(container, true);
  PositionFuncMap positionFuncMap;
  positionFuncMap = createPositionFunctionMap(container, true);
  
  int total_sequences = 0;
  int file_count = 0;
  
  std::cout << "Setting up " << nthreads << " thread(s) for processing." << std::endl;
  omp_set_num_threads(nthreads);
  
  auto process_chunk = [&](const std::vector<SearchResult>& chunk, int sequence_counter) {
    auto chunk_start_time = std::chrono::high_resolution_clock::now();
    
#pragma omp critical
{
  std::cout << "Processing chunk starting at sequence " << sequence_counter << " with " << chunk.size() << " sequences." << std::endl;
}
std::vector<std::string> chunk_sequences;
std::vector<std::string> chunk_ids;
std::vector<std::string> chunk_quality_scores;
for (const auto& sr : chunk) {
  if (sr.sequence.length() >= min_length) {
    chunk_sequences.push_back(sr.sequence);
    chunk_ids.push_back(sr.read_id);
    chunk_quality_scores.push_back(sr.qual_scores);
  }
}
std::vector<std::string> sigalign_vector = n_masked_sigalign(queries, query_names, chunk_sequences, chunk_ids, alignment_bounds, nthreads);
std::vector<std::string> sigrun_results = sigrun_cpp(sigalign_vector, read_layout, misalignment_threshold, nthreads, false);
if (verbose) {
#pragma omp critical
{
  std::cout << "Sigalign Vector Size: " << sigalign_vector.size() << std::endl;
  std::cout << "Sigrun Results Size: " << sigrun_results.size() << std::endl;
}
}

std::vector<std::string> chunk_processed_signatures;
std::vector<std::string> chunk_chunk_ids;
std::vector<std::string> chunk_chunk_sequences;
std::vector<std::string> chunk_chunk_quality_scores;

#pragma omp critical
{
  std::cout << "Storing processed results..." << std::endl;
}

for (size_t i = 0; i < sigrun_results.size(); ++i) {
  chunk_processed_signatures.push_back(sigrun_results[i]);
  size_t original_index = i % chunk_ids.size();
  // Only add sequences to output if they meet the criteria
  if (chunk_sequences[original_index].length() >= min_length &&
    sigrun_results[i].find(":undecided>") == std::string::npos) {
    chunk_chunk_ids.push_back(chunk_ids[original_index]);
    chunk_chunk_sequences.push_back(chunk_sequences[original_index]);
    chunk_chunk_quality_scores.push_back(chunk_quality_scores[original_index]);
  }
}

#pragma omp critical
{
  writeSigstrings(sigstringsFilePath, chunk_processed_signatures, compress);
}

sig_extraction(
  read_layout,
  misalignment_threshold,
  chunk_chunk_ids,
  chunk_chunk_sequences,
  chunk_chunk_quality_scores,
  chunk_processed_signatures,
  outputPath,
  write_fastq,
  compress,
  min_length,
  nthreads,
  false);

auto chunk_end_time = std::chrono::high_resolution_clock::now();
auto chunk_duration = std::chrono::duration_cast<std::chrono::seconds>(chunk_end_time - chunk_start_time);
auto total_duration = std::chrono::duration_cast<std::chrono::seconds>(chunk_end_time - start_time);

#pragma omp critical
{
  std::cout << "Chunk processing completed. Chunk time: " << format_duration(chunk_duration) << " | Total runtime: " << format_duration(total_duration) << std::endl;
}
  };
  
  for (const auto& file_path : input_files) {
#pragma omp critical
{
  std::cout << "Processing file " << ++file_count << "/" << input_files.size() << ": " << file_path << std::endl;
}
    gzFile file = gzopen(file_path.c_str(), "rb");
    if (!file) {
      Rcpp::warning("Failed to open file: " + file_path);
      continue;
    }
    kseq_t* seq = kseq_init(file);
    std::vector<SearchResult> chunk;
    int chunk_count = 0;
    bool is_fastq;
    int file_sequences = 0;
    while (true) {
      int prev_total = total_sequences;
      chunk = readFastqChunk(seq, file, chunkSize, is_fastq, total_sequences, max_sequences);
      if (chunk.empty()) break;
      
#pragma omp critical
{
  std::cout << "Processing chunk " << ++chunk_count << " of file " << file_count << std::endl;
  std::cout << "Chunk size: " << chunk.size() << ", Total sequences before: " << prev_total << ", after: " << total_sequences << std::endl;
}
process_chunk(chunk, file_sequences);
file_sequences += chunk.size();
if (max_sequences != -1 && total_sequences >= max_sequences) break;
    }
    
    kseq_destroy(seq);
    gzclose(file);
#pragma omp critical
{
  std::cout << "Finished processing file " << file_count << ". Sequences in this file: " << file_sequences << ", Total sequences processed: " << total_sequences << std::endl;
}
if (max_sequences != -1 && total_sequences >= max_sequences) {
#pragma omp critical
{
  std::cout << "Reached maximum number of sequences. Stopping processing." << std::endl;
}
  break;
}
  }
  if (total_sequences == 0) {
    Rcpp::warning("No sequences were processed.");
    return;
  }
  
  auto end_time = std::chrono::high_resolution_clock::now();
  auto total_duration = std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time);
  
#pragma omp critical
{
  std::cout << "Read extraction & demultiplexing completed. Total sequences processed: " << total_sequences << std::endl;
  std::cout << "Total runtime: " << format_duration(total_duration) << std::endl;
}
}

//' @export
// [[Rcpp::export]]
Rcpp::DataFrame misalignment_stream(
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
  
  std::cout << "Generating misalignment threshold data!\n";
  
  omp_set_num_threads(nthreads);
  
  auto process_chunk = [&](const std::vector<std::string>& chunk_sequences, const std::vector<std::string>& chunk_ids, int sequence_counter, int remaining) {
#pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < std::min(static_cast<int>(chunk_sequences.size()), remaining); ++i) {
      const auto& sequence = chunk_sequences[i];
      std::vector<AlignmentInfo> alignments(filtered_queries.size());
      
      for (size_t j = 0; j < filtered_queries.size(); ++j) {
        const auto& query = filtered_queries[j];
        EdlibAlignConfig config = edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_LOC, NULL, 0);
        EdlibAlignResult result = edlibAlign(query.c_str(), query.size(), sequence.c_str(), sequence.size(), config);
        alignments[j] = {
          result.editDistance,
          result.startLocations[0] + 1,
          result.endLocations[0] + 1,
          (result.numLocations > 1) ? result.editDistance : -1
        };
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
    Rcpp::Named("second_edit_distance") = second_edit_distances
  );
}

//' @export
// [[Rcpp::export]]
void tabulate_variable_sequences(const std::vector<std::string>& input_files, 
  const std::string& output_prefix, 
  const std::vector<std::string>& identity_elements,
  bool compress = true, 
  size_t batch_size = 50000) {
  auto write_output = [compress](gzFile& gz_file, const std::string& content) {
    if (compress) {
      gzwrite(gz_file, content.c_str(), content.length());
    } else {
      std::cout << "Compression should be true if using gzFile." << std::endl;
    }
  };
  
  std::vector<std::string> column_order = identity_elements;  // Use the provided identity elements as column names
  
  // Initialize gzFile for the master sequence output
  gzFile master_gz = gzopen((output_prefix + "mvs.csv.gz").c_str(), "wb");
  if (!master_gz) {
    Rcpp::stop("Failed to open the master sequence file.");
  }
  
  // Write the header for the master sequence file
  std::stringstream master_ss;
  master_ss << "id";  // Start with the ID column
  for (const auto& col : column_order) {
    master_ss << "," << col;
  }
  master_ss << "\n";
  write_output(master_gz, master_ss.str());  // Write the header to the file
  
  // Hash maps to store barcode counts for each barcode column
  // Each barcode map will have: {barcode => {non_concatenate_count, concatenate_count}}
  std::vector<std::unordered_map<std::string, std::pair<int, int>>> barcode_counts;
  for (const auto& col : column_order) {
    if (col.find("barcode") != std::string::npos) {
      barcode_counts.emplace_back();  // Initialize count for each barcode-related column
    }
  }
  
  // Process each input file
  for (const auto& file_path : input_files) {
    gzFile fp = gzopen(file_path.c_str(), "r");
    if (!fp) {
      Rcpp::warning("Error opening file: " + file_path);
      continue;
    }
    
    kseq_t* seq = kseq_init(fp);  // Initialize kseq for sequence reading
    int l;
    size_t processed = 0;
    
    // Process reads in batches
    while ((l = kseq_read(seq)) >= 0) {
      std::string header(seq->name.s);  // Read header
      std::istringstream iss(header);
      std::string token, id, barcode_seq;
      bool concatenate_flag = header.find("+FR_RF") != std::string::npos;  // Check for "+FR_RF" in the ID
      
      // Preallocate a vector to hold values for each identity element in the correct order
      std::vector<std::string> identity_values(identity_elements.size(), "");
      
      // Extract ID
      std::getline(iss, id, '|');  // Extract the read ID
      
      // Build the master sequence file line
      master_ss.str("");  // Clear the stringstream
      master_ss << id;  // Start the line with the read ID
      
      // Track which barcode column we are currently processing
      size_t barcode_column_index = 0;
      
      // Extract barcode, UMI, or other identity elements using find
      while (std::getline(iss, token, '|')) {
        for (size_t i = 0; i < identity_elements.size(); ++i) {
          const auto& elem = identity_elements[i];
          // Look for "elem:" (e.g., "barcode:", "umi:")
          size_t pos = token.find(elem + ":");
          if (pos != std::string::npos) {
            barcode_seq = token.substr(pos + elem.size() + 1);  // Extract the sequence part
            identity_values[i] = barcode_seq;  // Store the extracted sequence in the correct position
            
            // If this element is a "barcode", increment the count
            if (elem.find("barcode") != std::string::npos) {
              std::string revcomp_seq = revcomp_cpp(barcode_seq);  // Compute the reverse complement
              if (barcode_counts[barcode_column_index].count(revcomp_seq)) {
                // If the reverse complement exists, collapse the counts into the forward complement
                barcode_counts[barcode_column_index][barcode_seq].first += barcode_counts[barcode_column_index][revcomp_seq].first;
                barcode_counts[barcode_column_index][barcode_seq].second += barcode_counts[barcode_column_index][revcomp_seq].second;
                barcode_counts[barcode_column_index].erase(revcomp_seq);  // Remove the reverse complement entry
              }
              
              if (concatenate_flag) {
                barcode_counts[barcode_column_index][barcode_seq].second++;  // Increment concatenate count
              } else {
                barcode_counts[barcode_column_index][barcode_seq].first++;   // Increment non-concatenate count
              }
              barcode_column_index++;
            }
            break;  // No need to check the remaining identity elements for this token
          }
        }
      }
      
      // Append all values in the correct order (matching the header)
      for (const auto& value : identity_values) {
        master_ss << "," << value;
      }
      
      // Complete the line and write it to the file
      master_ss << "\n";
      write_output(master_gz, master_ss.str());
      
      processed++;
      
      if (processed % batch_size == 0) {
        std::cout << "Processed " << processed << " reads from " << file_path << std::endl;
      }
    }
    
    kseq_destroy(seq);  // Clean up kseq resources
    gzclose(fp);  // Close the input file
  }
  
  gzclose(master_gz);  // Close the master sequence file at the end
  
  // Write the barcode counts to separate CSV files
  for (size_t i = 0; i < barcode_counts.size(); ++i) {
    std::string output_file = output_prefix + "barcode";
    if (i > 0) {
      output_file += "_" + std::to_string(i);  // Add a suffix for subsequent barcode columns
    }
    output_file += ".csv.gz";
    
    gzFile barcode_gz = gzopen(output_file.c_str(), "wb");
    if (!barcode_gz) {
      Rcpp::warning("Error opening file: " + output_file);
      continue;
    }
    
    std::stringstream ss;
    ss << "barcode,non_concatenate_count,concatenate_count\n";  // Write header
    for (const auto& entry : barcode_counts[i]) {
      ss << entry.first << "," << entry.second.first << "," << entry.second.second << "\n";  // Write barcode and counts
    }
    write_output(barcode_gz, ss.str());
    gzclose(barcode_gz);  // Close the barcode output file
  }
}

//' @export
// [[Rcpp::export]]
void tabulate_sigs(const std::string& file_path,
                   const std::string& output_prefix,
                   int chunk_size = 50000,
                   int max_sequences = -1,
                   int nthreads = 1,
                   bool compress = true,
                   const std::string& output_type = "summary") {
  omp_set_num_threads(nthreads);
  auto write_output = [compress](const std::string& filename, const std::string& content) {
    if (compress) {
      gzFile file = gzopen(filename.c_str(), "wb");
      if (!file) {
        Rcpp::stop("Error opening output file: " + filename);
      }
      gzwrite(file, content.c_str(), content.length());
      gzclose(file);
    } else {
      std::ofstream file(filename);
      if (!file.is_open()) {
        Rcpp::stop("Error opening output file: " + filename);
      }
      file << content;
      file.close();
    }
  };
  auto extract_id = [](const std::string& sigstring) {
    std::regex id_regex("<\\d+:(.*?):");
    std::smatch match;
    if (std::regex_search(sigstring, match, id_regex)) {
      return match[1].str();
    }
    return std::string("");
  };
  auto extract_direction = [](const std::string& sigstring) {
    std::regex direction_regex(":([^:]*?)>");
    std::smatch match;
    if (std::regex_search(sigstring, match, direction_regex)) {
      return match[1].str();
    }
    return std::string("");
  };
  std::string data_output_file = output_prefix + "/sigparsed_info.csv" + (compress ? ".gz" : "");
  std::string summary_output_file = output_prefix + "/summary_table.csv" + (compress ? ".gz" : "");
  std::stringstream data_ss;
  data_ss << "id,element,edit_distance,start_pos,stop_pos,concatenate,direction\n";
  std::unordered_map<std::string, int> summary_map;
  int sequence_count = 0;
  gzFile file = gzopen(file_path.c_str(), "r");
  if (!file) {
    Rcpp::stop("Failed to open file: " + file_path);
  }
  kstream_t* ks = ks_init(file);
  kstring_t str = {0, 0, 0};
  std::regex element_regex(R"(([^:]+):(\d+):(\d+):(\d+))");
  while ((max_sequences == -1 || sequence_count < max_sequences) && ks_getuntil(ks, '\n', &str, 0) >= 0) {
    std::vector<std::string> sigstrings;
    sigstrings.push_back(str.s);
    for (int i = 1; i < chunk_size && (max_sequences == -1 || sequence_count + i < max_sequences) && ks_getuntil(ks, '\n', &str, 0) >= 0; ++i) {
      sigstrings.push_back(str.s);
    }
#pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < sigstrings.size(); ++i) {
      std::string sigstring = sigstrings[i];
      std::string id = extract_id(sigstring);
      std::string direction = extract_direction(sigstring);
      std::string concatenate = (sigstring.find("+FR_RF") != std::string::npos) ? "TRUE" : "FALSE";
      std::smatch match;
      auto start = std::sregex_iterator(sigstring.begin(), sigstring.end(), element_regex);
      auto end = std::sregex_iterator();
      std::stringstream local_data_ss;
      for (std::sregex_iterator it = start; it != end; ++it) {
        std::smatch match = *it;
        std::string element = "|" + match[1].str();
        int edit_distance = std::stoi(match[2].str());
        int start_pos = std::stoi(match[3].str());
        int stop_pos = std::stoi(match[4].str());
        local_data_ss << id << "," << element << "," << edit_distance << "," << start_pos << 
          "," << stop_pos << "," << concatenate << "," << direction << "\n";
      }
#pragma omp critical
{
  if (output_type == "diagnostic") {
    data_ss << local_data_ss.str();
  }
  summary_map[concatenate + ":" + direction]++;
}
    }
    
    sequence_count += sigstrings.size();
    if (sequence_count % chunk_size == 0) {
      Rcpp::checkUserInterrupt();
#pragma omp critical
{
  std::cout << "Processed " << sequence_count << " sigstrings." << std::endl;
}
    }
  }
  ks_destroy(ks);
  gzclose(file);
  free(str.s);
  if (output_type == "diagnostic") {
    write_output(data_output_file, data_ss.str());
  }
  std::stringstream summary_ss;
  summary_ss << "concatenate,direction,unique_id_count\n";
  for (const auto& entry : summary_map) {
    std::stringstream ss(entry.first);
    std::string concatenate_str, direction_str;
    std::getline(ss, concatenate_str, ':');
    std::getline(ss, direction_str, ':');
    summary_ss << concatenate_str << "," << direction_str << "," << entry.second << "\n";
  }
  write_output(summary_output_file, summary_ss.str());
  std::cout << "Completed processing " << sequence_count << " sigstrings." << std::endl;
}

std::string find_and_replace_barcode(
  const std::string& header, 
  const std::string& barcode_id,
  const std::unordered_map<std::string, std::string>& barcode_map,
  const std::vector<std::string>& correct_barcodes,
  bool &barcode_found,
  bool &barcode_corrected) {
  // Look for the barcode identifier (e.g., "barcode:") in the header
  size_t barcode_start = header.find(barcode_id + ":");
  if (barcode_start == std::string::npos) {
    barcode_found = false;
    barcode_corrected = false;
    return header;
  }
  
  // Find where the barcode sequence starts (just after 'barcode_id:')
  barcode_start += barcode_id.size() + 1;
  
  // Find where the barcode sequence ends (next pipe `|` or end of line)
  size_t barcode_end = header.find_first_of("|", barcode_start);
  if (barcode_end == std::string::npos) {
    barcode_end = header.size();
  }
  
  // Extract the barcode sequence
  std::string barcode = header.substr(barcode_start, barcode_end - barcode_start);
  
  // Check if the barcode is in the barcode_map (incorrect barcodes)
  auto it = barcode_map.find(barcode);
  if (it != barcode_map.end()) {
    // Replace the incorrect barcode with the corrected one
    std::string corrected_header = header;
    corrected_header.replace(barcode_start, barcode.size(), it->second);
    barcode_found = true;
    barcode_corrected = true;
    return corrected_header;
  }
  
  // Check if the barcode is in the correct_barcodes list
  if (std::find(correct_barcodes.begin(), correct_barcodes.end(), barcode) != correct_barcodes.end()) {
    barcode_found = true;
    barcode_corrected = false;
    return header;
  }
  
  // If the barcode is neither incorrect nor correct, return the original header and mark it as not found
  barcode_found = false;
  barcode_corrected = false;
  return header;
}

//' @export
// [[Rcpp::export]]
Rcpp::DataFrame fastq_correction(
    std::string input_file,
    std::string output_file,
    Rcpp::StringVector incorrect_bcs,
    Rcpp::StringVector correct_bcs,
    std::string barcode_id = "barcode",
    int chunk_size = 50000,
    int print_freq = 10000,
    int nthreads = 4,
    int downsample = -1) {
  // Set the number of OpenMP threads
  omp_set_num_threads(nthreads);
  // Create a map for barcode correction
  std::unordered_map<std::string, std::string> barcodeMap;
  for (int i = 0; i < incorrect_bcs.size(); ++i) {
    barcodeMap[Rcpp::as<std::string>(incorrect_bcs[i])] = Rcpp::as<std::string>(correct_bcs[i]);
  }
  std::vector<std::string> correct_bcsVec = Rcpp::as<std::vector<std::string>>(correct_bcs);
  // Open input file
  gzFile fp = gzopen(input_file.c_str(), "r");
  if (!fp) {
    Rcpp::stop("Failed to open input file");
  }
  // Open output file
  gzFile outFp = gzopen(output_file.c_str(), "wb");
  if (!outFp) {
    gzclose(fp);
    Rcpp::stop("Failed to open output file");
  }
  // Initialize kseq for reading fastq
  kseq_t* seq = kseq_init(fp);
  // Variables to track read counts
  std::atomic<size_t> totalReads{0};
  std::atomic<size_t> keptReads{0};
  std::atomic<size_t> correctedReads{0};
  std::atomic<size_t> discardedReads{0};
  // Chunk storage for parallel processing
  std::vector<std::string> chunkHeaders;
  std::vector<std::string> chunkSequences;
  std::vector<std::string> chunkQuals;
  std::vector<std::string> chunkComments;
  // Process each record in the fastq file
  while (kseq_read(seq) >= 0) {
    if (downsample != -1 && totalReads >= static_cast<size_t>(downsample)) {
      break;
    }
    ++totalReads;
    // Capture fastq elements (header, sequence, quality, etc.)
    std::string header = "@" + std::string(seq->name.s);
    std::string seqStr(seq->seq.s, seq->seq.l);
    std::string qualStr(seq->qual.s, seq->qual.l);
    std::string commentStr = (seq->comment.l > 0) ? std::string(seq->comment.s, seq->comment.l) : "+";
    
    // Store elements in chunk
    chunkHeaders.push_back(header);
    chunkSequences.push_back(seqStr);
    chunkQuals.push_back(qualStr);
    chunkComments.push_back(commentStr);
    
    // If chunk reaches the chunk size, process it
    if (chunkHeaders.size() >= static_cast<size_t>(chunk_size)) {
#pragma omp parallel
{
  // Local thread storage to avoid critical section contention
  std::vector<std::string> localHeaders;
  std::vector<std::string> localSequences;
  std::vector<std::string> localQuals;
  std::vector<std::string> localComments;
  size_t localDiscarded = 0;
  
#pragma omp for
  for (size_t i = 0; i < chunkHeaders.size(); ++i) {
    bool barcodeFound = false;
    bool barcodeCorrected = false;
    
    // Find and replace the barcode in the header
    std::string processedHeader = find_and_replace_barcode(
      chunkHeaders[i],
      barcode_id,
      barcodeMap,
      correct_bcsVec,
      barcodeFound,
      barcodeCorrected);
    
    if (barcodeFound) {
      localHeaders.push_back(processedHeader);
      localSequences.push_back(chunkSequences[i]);
      localQuals.push_back(chunkQuals[i]);
      localComments.push_back(chunkComments[i]);
      ++keptReads;
      if (barcodeCorrected) {
        ++correctedReads;
      }
    } else {
      ++localDiscarded;
    }
  }
  
#pragma omp critical
{
  // Write out local data to file
  for (size_t i = 0; i < localHeaders.size(); ++i) {
    gzwrite(outFp, localHeaders[i].c_str(), localHeaders[i].size());
    gzwrite(outFp, "\n", 1);
    gzwrite(outFp, localSequences[i].c_str(), localSequences[i].size());
    gzwrite(outFp, "\n", 1);
    gzwrite(outFp, localComments[i].c_str(), localComments[i].size());
    gzwrite(outFp, "\n", 1);
    gzwrite(outFp, localQuals[i].c_str(), localQuals[i].size());
    gzwrite(outFp, "\n", 1);
  }
  discardedReads += localDiscarded;
}
}
      
      // Clear chunk after processing
      chunkHeaders.clear();
      chunkSequences.clear();
      chunkQuals.clear();
      chunkComments.clear();
    }
    
    // Log progress
    if (totalReads % print_freq == 0) {
      std::cout << "\rProcessed " << totalReads << " reads, kept " << keptReads 
                << ", corrected " << correctedReads << ", discarded " << discardedReads << std::flush;
    }
  }
  
  // Process any remaining reads
  if (!chunkHeaders.empty()) {
#pragma omp parallel
{
  // Local thread storage for processing remaining chunk
  std::vector<std::string> localHeaders;
  std::vector<std::string> localSequences;
  std::vector<std::string> localQuals;
  std::vector<std::string> localComments;
  size_t localDiscarded = 0;
  
#pragma omp for
  for (size_t i = 0; i < chunkHeaders.size(); ++i) {
    bool barcodeFound = false;
    bool barcodeCorrected = false;
    
    // Process barcode replacement
    std::string processedHeader = find_and_replace_barcode(
      chunkHeaders[i],
      barcode_id,
      barcodeMap,
      correct_bcsVec,
      barcodeFound,
      barcodeCorrected);
    
    if (barcodeFound) {
      localHeaders.push_back(processedHeader);
      localSequences.push_back(chunkSequences[i]);
      localQuals.push_back(chunkQuals[i]);
      localComments.push_back(chunkComments[i]);
      ++keptReads;
      if (barcodeCorrected) {
        ++correctedReads;
      }
    } else {
      ++localDiscarded;
    }
  }
  
#pragma omp critical
{
  // Write remaining data to file
  for (size_t i = 0; i < localHeaders.size(); ++i) {
    gzwrite(outFp, localHeaders[i].c_str(), localHeaders[i].size());
    gzwrite(outFp, "\n", 1);
    gzwrite(outFp, localSequences[i].c_str(), localSequences[i].size());
    gzwrite(outFp, "\n", 1);
    gzwrite(outFp, localComments[i].c_str(), localComments[i].size());
    gzwrite(outFp, "\n", 1);
    gzwrite(outFp, localQuals[i].c_str(), localQuals[i].size());
    gzwrite(outFp, "\n", 1);
  }
  discardedReads += localDiscarded;
}
}
  }
  
  // Final log
  std::cout << "\nFinished processing. Total reads: " << totalReads 
            << ", kept: " << keptReads << ", corrected: " << correctedReads 
            << ", discarded: " << discardedReads << "." << std::endl;
  
  // Cleanup
  kseq_destroy(seq);
  gzclose(fp);
  gzclose(outFp);
  
  return Rcpp::DataFrame::create(
    Rcpp::Named("total_reads") = totalReads.load(),
    Rcpp::Named("kept_reads") = keptReads.load(),
    Rcpp::Named("corrected_reads") = correctedReads.load(),
    Rcpp::Named("discarded_reads") = discardedReads.load()
  );
}
