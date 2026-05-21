#include "rad.h"

using namespace std;
using namespace Rcpp;

std::string revcomp_cpp(const std::string& sequence) {
  if (sequence.empty()) {
    return sequence;
  }

  std::string rc(sequence.rbegin(), sequence.rend());
  std::transform(rc.begin(), rc.end(), rc.begin(), [](char c) {
    switch (c) {
    case 'A': return 'T';
    case 'T': return 'A';
    case 'G': return 'C';
    case 'C': return 'G';
    case 'N': return 'N';
    case 'n': return 'n';
    default: return c;
    }
  });
  return rc;
}

//' @export
// [[Rcpp::export]]
Rcpp::StringVector revcomp(Rcpp::StringVector sequences) {
  int n = sequences.size();
  Rcpp::StringVector results(n);
  int j = 0; // Counter for results
  for (int i = 0; i < n; ++i) {
    if (sequences[i] == NA_STRING) {
      results[j] = NA_STRING;
      continue;
    }
    results[j] = revcomp_cpp(Rcpp::as<std::string>(sequences[i]));
    ++j;
  }
  return results;
}

std::vector<int64_t> sequence_to_bits_cpp(const std::string& sequence) {
  if (sequence.empty()) {
    return {0};
  }

  const int chunk_size = 32; // Number of bases per chunk
  std::vector<int64_t> results;
  for (size_t i = 0; i < sequence.size(); i += chunk_size) {
    int64_t result = 0;
    size_t end = std::min(i + chunk_size, sequence.size());
    for (size_t j = i; j < end; ++j) {
      char c = sequence[j];
      result <<= 2;
      switch (c) {
      case 'A': break;
      case 'C': result |= 1; break;
      case 'T': result |= 2; break;
      case 'G': result |= 3; break;
      default:
        return {0};
      }
    }
    results.push_back(result);
  }
  return results;
}

//' @export
// [[Rcpp::export]]
Rcpp::NumericVector sequence_to_bits(Rcpp::StringVector sequences) {
  int n = sequences.size();
  std::vector<int64_t> results;
  for(int i = 0; i < n; ++i) {
    std::string sequence = Rcpp::as<std::string>(sequences[i]);
    std::vector<int64_t> chunk_results = sequence_to_bits_cpp(sequence);
    results.insert(results.end(), chunk_results.begin(), chunk_results.end());
  }
  return Rcpp::toInteger64(results);
}

//' @export
// [[Rcpp::export]]
List sequence_to_bitlist(Rcpp::StringVector sequences) {
  int n = sequences.size();
  List results(n);
  for(int i = 0; i < n; ++i) {
    std::string sequence = as<std::string>(sequences[i]);
    std::vector<int64_t> chunk_results = sequence_to_bits_cpp(sequence);
    results[i] = wrap(chunk_results);
  }
  return results;
}

std::string bits_to_sequence_cpp(const std::vector<int64_t>& input, int sequence_length) {
  int n = input.size();
  std::string final_sequence;
  int chunk_size = 32;
  for (int i = 0; i < n; ++i) {
    int64_t int64_code = input[i];
    std::string result;
    int bases_to_decode = (i < n - 1) ? chunk_size : sequence_length % chunk_size;
    bases_to_decode = (bases_to_decode == 0) ? chunk_size : bases_to_decode;
    for (int j = 0; j < bases_to_decode; ++j) {
      switch (int64_code & 3) {
      case 0: result.insert(result.begin(), 'A'); break;
      case 1: result.insert(result.begin(), 'C'); break;
      case 2: result.insert(result.begin(), 'T'); break;
      case 3: result.insert(result.begin(), 'G'); break;
      }
      int64_code >>= 2;
    }
    final_sequence += result;
  }
  return final_sequence;
}

std::vector<unsigned int> bits_to_uint_cpp(const std::vector<int64_t>& input, int sequence_length) {
  int n = input.size();
  std::vector<unsigned int> final_sequence;
  int chunk_size = 32; // Assuming 32 bases per int64_t value
  for (int i = 0; i < n; ++i) {
    int64_t int64_code = input[i];
    int bases_to_decode = (i < n - 1) ? chunk_size : sequence_length % chunk_size;
    bases_to_decode = (bases_to_decode == 0) ? chunk_size : bases_to_decode;
    for (int j = 0; j < bases_to_decode; ++j) {
      unsigned int base = int64_code & 3; // Extract the last two bits
      final_sequence.push_back(base); // Add the base as an unsigned int to the sequence
      int64_code >>= 2; // Move to the next base
    }
  }
  return final_sequence;
}

//' @export
// [[Rcpp::export]]
Rcpp::StringVector bits_to_sequence(Rcpp::NumericVector input, Rcpp::IntegerVector sequence_length = 16) {
  // Convert the entire NumericVector to std::vector<int64_t>
  std::vector<int64_t> cpp_input_list = Rcpp::as<std::vector<int64_t>>(input);
  int n = cpp_input_list.size(); // Number of sequences to process
  bool use_sl = sequence_length.size() == 1;
  int seq_length = use_sl ? sequence_length[0] : 0; // Use this if only one length is provided
  Rcpp::StringVector final_sequences(n);
  for (int i = 0; i < n; ++i) {
    std::vector<int64_t> bits_input = {cpp_input_list[i]}; // Single element vector
    int seqlength = use_sl ? seq_length : sequence_length[i]; // Use the single length or individual lengths
    // Process the single sequence
    std::string sequence = bits_to_sequence_cpp(bits_input, seqlength);
    final_sequences[i] = sequence;
  }
  return final_sequences;
}

std::string bits_to_hex_cpp(const std::vector<int64_t>& input, int sequence_length) {
  std::string final_sequence;
  const int chunk_size = 32;
  for (size_t i = 0; i < input.size(); ++i) {
    int64_t int64_code = input[i];
    std::string result;
    
    int k_mers_to_decode = (i < input.size() - 1) ? chunk_size : sequence_length % chunk_size;
    k_mers_to_decode = (k_mers_to_decode == 0) ? chunk_size : k_mers_to_decode;
    
    for (int j = 0; j < k_mers_to_decode; ++j) {
      int k_mer = int64_code & 0xF;
      if (k_mer < 10) {
        result.insert(result.begin(), '0' + k_mer);
      } else {
        result.insert(result.begin(), 'A' + (k_mer - 10));
      }
      int64_code >>= 4;
    }
    final_sequence += result;
  }
  return final_sequence;
}

//' @export
// [[Rcpp::export]]
StringVector bits_to_hex(Rcpp::NumericVector input, int sequence_length = 16) {
  Rcpp::StringVector hex_sequences(input.size());
  std::vector<int64_t> cpp_input = Rcpp::as<std::vector<int64_t>>(input);
  int slen = sequence_length/2; 
  //for consistency and compression's sake, we convert from a rep of 2 bits->one char to 4 bits->one char
  for(size_t i = 0; i < cpp_input.size(); ++i) {
    std::vector<int64_t> single_input = {cpp_input[i]};
    std::string hex_sequence = bits_to_hex_cpp(single_input, slen);
    hex_sequences[i] = hex_sequence;
  }
  return hex_sequences;
}

std::vector<int64_t> revcomp_bits_cpp(const std::vector<int64_t>& input, int sequence_length = 16) {
  std::vector<int64_t> reversed(input.rbegin(), input.rend()); // Reverse the vector
  for (int64_t& chunk : reversed) {
    int64_t result = 0;
    for (int i = 0; i < sequence_length; i++) {
      result <<= 2; // Shift result to make room for the new bits
      int64_t base = chunk & 3; // Extract the two rightmost bits
      switch (base) {
      case 0: result |= 2; break; // A (00) becomes T (10)
      case 1: result |= 3; break; // C (01) becomes G (11)
      case 2: result |= 0; break; // T (10) becomes A (00)
      case 3: result |= 1; break; // G (11) becomes C (01)
      }
      chunk >>= 2; // Move to the next base in the chunk
    }
    // Assign the reversed and complemented chunk back
    chunk = result;
  }
  return reversed;
}

//' @export
// [[Rcpp::export]]
Rcpp::NumericVector revcomp_bits(Rcpp::NumericVector input, Rcpp::IntegerVector sequence_length = 16) {
  std::vector<int64_t> cpp_input_list = Rcpp::as<std::vector<int64_t>>(input);
  int n = cpp_input_list.size(); // Number of sequences to process
  bool use_sl = sequence_length.size() == 1;
  int seq_length = use_sl ? sequence_length[0] : 0; // Use this if only one length is provided
  std::vector<int64_t> results;
  for (int i = 0; i < n; ++i) {
    std::vector<int64_t> bits_input = {cpp_input_list[i]}; // Single element vector
    int seqlength = use_sl ? seq_length : sequence_length[i]; // Use the single length or individual lengths
    std::vector<int64_t> subresult = revcomp_bits_cpp(bits_input, seqlength);
    results.insert(results.end(), subresult.begin(), subresult.end());
  }
  return Rcpp::toInteger64(results);
}

// Function to calculate the number of differing bits for int64_t sequences
int hamming_distance(int64_t a, int64_t b) {
  int bit_diff_count = std::bitset<64>(a ^ b).count();
  return std::ceil(bit_diff_count / 2.0);  // Base Hamming distance
}

//' @export
// [[Rcpp::export]]
int hamming_distance_rcpp(Rcpp::NumericVector a, Rcpp::NumericVector b) {
  // Convert Rcpp::NumericVector (integer64) to std::vector<int64_t>
  std::vector<int64_t> a_int = Rcpp::fromInteger64(a);
  std::vector<int64_t> b_int = Rcpp::fromInteger64(b);
  
  if (a_int.size() != 1 || b_int.size() != 1) {
    Rcpp::stop("Expecting single integer64 values.");
  }
  
  // Return the Hamming distance between the first (and only) elements
  return hamming_distance(a_int[0], b_int[0]);
}

// Function to generate mutations for a single sequence
std::vector<int64_t> generate_quaternary_mutations_cpp_heap(int64_t sequence, int sequence_length) {
  // Allocate mutations vector on the heap
  auto mutations = std::make_unique<std::vector<int64_t>>();
  for (int position = 0; position < sequence_length; ++position) {
    int64_t original_nucleotide = (sequence >> (2 * position)) & 0b11;
    for (int64_t new_nucleotide = 0; new_nucleotide < 4; ++new_nucleotide) {
      if (new_nucleotide != original_nucleotide) {
        int64_t mutation_mask = new_nucleotide << (2 * position);
        int64_t clear_mask = ~(0b11LL << (2 * position));
        int64_t mutated_sequence = (sequence & clear_mask) | mutation_mask;
        mutations->push_back(mutated_sequence);
      }
    }
  }
  // Return the vector stored on the heap
  return *mutations;
}

// Function to generate mutations for a single sequence
std::vector<int64_t> generate_quaternary_mutations_cpp(int64_t sequence, int sequence_length) {
  std::vector<int64_t> mutations;
  for (int position = 0; position < sequence_length; ++position) {
    int64_t original_nucleotide = (sequence >> (2 * position)) & 0b11;
    for (int64_t new_nucleotide = 0; new_nucleotide < 4; ++new_nucleotide) {
      if (new_nucleotide != original_nucleotide) {
        int64_t mutation_mask = new_nucleotide << (2 * position);
        int64_t clear_mask = ~(0b11LL << (2 * position));
        int64_t mutated_sequence = (sequence & clear_mask) | mutation_mask;
        mutations.push_back(mutated_sequence);
      }
    }
  }
  return mutations;
}

//' @export
// [[Rcpp::export]]
Rcpp::NumericVector generate_quaternary_mutations(Rcpp::NumericVector sequences, Rcpp::IntegerVector sequence_length = 32) {
  std::vector<int64_t> cpp_input_list = Rcpp::as<std::vector<int64_t>>(sequences);
  bool use_sl = sequence_length.size() == 1;
  int seq_length = use_sl ? sequence_length[0] : 0;
  std::vector<int64_t> results;
  for (size_t i = 0; i < cpp_input_list.size(); ++i) {
    int64_t quaternary_input = cpp_input_list[i];
    int seqlength = use_sl ? seq_length : sequence_length[i];
    std::vector<int64_t> subresult = generate_quaternary_mutations_cpp(quaternary_input, seqlength);
    results.insert(results.end(), subresult.begin(), subresult.end());
  }
  return Rcpp::wrap(results);
}

// Function to generate recursive mutations
std::vector<int64_t> generate_recursive_quaternary_mutations_cpp(
    int64_t sequence, 
  int mutation_rounds, 
  int sequence_length = 32) {
  std::unordered_set<int64_t> all_mutations{sequence};
  std::vector<int64_t> current_round{sequence};
  for (int round = 0; round < mutation_rounds; ++round) {
    std::vector<int64_t> next_round;
    for (const auto& seq : current_round) {
      auto mutations = generate_quaternary_mutations_cpp(seq, sequence_length);
      for (const auto& mutation : mutations) {
        if (all_mutations.insert(mutation).second) {
          next_round.push_back(mutation);
        }
      }
    }
    if (next_round.empty()) break;
    current_round = std::move(next_round);
  }
  all_mutations.erase(sequence);
  return std::vector<int64_t>(all_mutations.begin(), all_mutations.end());
}

std::vector<int64_t> generate_recursive_quaternary_mutations_cpp_v2(
    int64_t sequence,
    int mutation_rounds,
    int sequence_length = 32) {
  std::unordered_set<int64_t> all_mutations{sequence};
  std::queue<std::pair<int64_t, int>> sequences_to_mutate;
  sequences_to_mutate.push({sequence, 0});
  
  while (!sequences_to_mutate.empty()) {
    auto current = sequences_to_mutate.front();
    int64_t current_seq = current.first;
    int current_round = current.second;
    sequences_to_mutate.pop();
    
    if (current_round < mutation_rounds) {
      std::vector<int64_t> new_mutations = generate_quaternary_mutations_cpp(current_seq, sequence_length);
      for (const auto& mutation : new_mutations) {
        if (all_mutations.insert(mutation).second) {
          sequences_to_mutate.push(std::make_pair(mutation, current_round + 1));
        }
      }
    }
  }
  return std::vector<int64_t>(all_mutations.begin(), all_mutations.end());
}

std::vector<int64_t> generate_recursive_quaternary_mutations_cpp_v3(
    int64_t sequence,
    int mutation_rounds,
    int sequence_length = 32) {
  // Use unique_ptr to allocate containers on the heap
  auto all_mutations = std::make_unique<std::unordered_set<int64_t>>();
  all_mutations->insert(sequence);
  
  auto sequences_to_mutate = std::make_unique<std::queue<std::pair<int64_t, int>>>();
  sequences_to_mutate->push({sequence, 0});
  
  while (!sequences_to_mutate->empty()) {
    auto current = sequences_to_mutate->front();
    int64_t current_seq = current.first;
    int current_round = current.second;
    sequences_to_mutate->pop();
    
    if (current_round < mutation_rounds) {
      // Store new mutations on the heap
      auto new_mutations = std::make_unique<std::vector<int64_t>>(generate_quaternary_mutations_cpp(current_seq, sequence_length));
      for (const auto& mutation : *new_mutations) {
        if (all_mutations->insert(mutation).second) {
          sequences_to_mutate->push(std::make_pair(mutation, current_round + 1));
        }
      }
    }
  }
  
  // Return vector created from the heap-allocated unordered_set
  return std::vector<int64_t>(all_mutations->begin(), all_mutations->end());
}

//' @export
// [[Rcpp::export]]
Rcpp::NumericVector generate_recursive_quaternary_mutations(
    Rcpp::NumericVector sequences,
    int mutation_rounds = 2, 
  Rcpp::IntegerVector sequence_length = 32) {
  std::vector<int64_t> cpp_input_list = Rcpp::as<std::vector<int64_t>>(sequences);
  std::vector<int64_t> results;
  bool use_sl = sequence_length.size() == 1;
  int seq_length = use_sl ? sequence_length[0] : 0;
  
  for (size_t i = 0; i < cpp_input_list.size(); ++i) {
    int64_t quaternary_input = cpp_input_list[i];
    int seqlength = use_sl ? seq_length : sequence_length[i];
    std::vector<int64_t> subresult = generate_recursive_quaternary_mutations_cpp(quaternary_input, mutation_rounds, seqlength); 
    results.insert(results.end(), subresult.begin(), subresult.end());
  }
  return Rcpp::wrap(results);
}
