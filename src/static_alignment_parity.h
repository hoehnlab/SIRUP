#pragma once

#include "rad.h"
#include "include/external/ssw/ssw_cpp.h"

namespace rad_static {

struct StaticAlignment {
  std::vector<std::pair<int, int>> positions;
  int edit_distance = -1;
  bool success = false;
  std::string seq;
  std::string cigar;
};

inline int cigar_max_consecutive_matches(const char* cigar) {
  int max_matches = 0;
  int current_matches = 0;
  int number = 0;
  for (const char* p = cigar; *p; ++p) {
    if (std::isdigit(*p)) {
      number = number * 10 + (*p - '0');
      continue;
    }
    if (*p == 'M' || *p == '=') {
      current_matches = number;
    } else {
      current_matches = 0;
    }
    max_matches = std::max(max_matches, current_matches);
    number = 0;
  }
  return max_matches;
}

inline int cigar_max_consecutive_matches_with_indels(const char* cigar) {
  std::vector<std::pair<int, char>> segments;
  while (*cigar) {
    int count = 0;
    while (std::isdigit(*cigar)) {
      count = count * 10 + (*cigar - '0');
      ++cigar;
    }
    if (!*cigar) {
      break;
    }
    segments.push_back({count, *cigar});
    ++cigar;
  }

  int best = 0;
  for (size_t i = 0; i < segments.size(); ++i) {
    if (segments[i].second != '=' && segments[i].second != 'M') {
      continue;
    }
    int current = 0;
    bool used_indel = false;
    for (size_t j = i; j < segments.size(); ++j) {
      const char op = segments[j].second;
      const int count = segments[j].first;
      if (op == '=' || op == 'M') {
        current += count;
      } else if ((op == 'I' || op == 'D' || op == 'X') && count <= 2) {
        if (used_indel) {
          break;
        }
        used_indel = true;
      } else {
        break;
      }
      best = std::max(best, current);
    }
  }
  return best;
}

inline int cigar_edit_distance(const char* cigar) {
  int edit_distance = 0;
  int number = 0;
  while (*cigar) {
    if (std::isdigit(*cigar)) {
      number = number * 10 + (*cigar - '0');
    } else {
      if (*cigar == 'X' || *cigar == 'I' || *cigar == 'D') {
        edit_distance += number;
      }
      number = 0;
    }
    ++cigar;
  }
  return edit_distance;
}

inline std::vector<std::pair<int, int>> collapse_intervals(
    std::vector<std::pair<int, int>> intervals) {
  std::vector<std::pair<int, int>> out;
  if (intervals.empty()) {
    return out;
  }
  std::sort(intervals.begin(), intervals.end(),
            [](const auto& a, const auto& b) { return a.first < b.first; });
  std::pair<int, int> current = intervals.front();
  for (size_t i = 1; i < intervals.size(); ++i) {
    if (intervals[i].first <= current.second) {
      current.second = std::max(current.second, intervals[i].second);
    } else {
      out.push_back(current);
      current = intervals[i];
    }
  }
  out.push_back(current);
  return out;
}

inline StaticAlignment align_static_element(
    const std::string& query,
    const std::string& target,
    bool verbose,
    int max_edit_distance = -1,
    int expected_start = 1,
    int expected_end = -1) {
  StaticAlignment alignment;
  const int min_match_bases = 5;
  if (expected_end < 0) {
    expected_end = static_cast<int>(target.size());
  }

  const int num_eq = 13;
  EdlibEqualityPair custom_equalities[num_eq] = {
    {'A', 'a'}, {'C', 'c'}, {'T', 't'}, {'G', 'g'},
    {'A', 'x'}, {'C', 'x'}, {'T', 'x'}, {'G', 'x'},
    {'N', 'A'}, {'N', 'C'}, {'N', 'T'}, {'N', 'G'},
    {'N', 'x'}
  };

  EdlibAlignConfig config = edlibNewAlignConfig(
    max_edit_distance,
    EDLIB_MODE_HW,
    EDLIB_TASK_LOC,
    custom_equalities,
    num_eq
  );

  EdlibAlignResult edlib_result = edlibAlign(
    query.c_str(),
    static_cast<int>(query.size()),
    target.c_str(),
    static_cast<int>(target.size()),
    config
  );

  std::vector<std::pair<int, int>> candidates;
  if (edlib_result.status == EDLIB_STATUS_OK &&
      edlib_result.numLocations > 0 &&
      edlib_result.editDistance > -1) {
    std::vector<std::pair<int, int>> intervals;
    intervals.reserve(edlib_result.numLocations);
    for (int i = 0; i < edlib_result.numLocations; ++i) {
      intervals.push_back({edlib_result.startLocations[i] + 1, edlib_result.endLocations[i] + 1});
    }
    candidates = collapse_intervals(std::move(intervals));
  }

  if (candidates.size() == 1) {
    const auto cand = candidates.front();
    if (cand.first >= expected_start && cand.second <= expected_end) {
      alignment.positions.push_back(cand);
      alignment.success = true;
      alignment.edit_distance = edlib_result.editDistance;
      alignment.seq = target.substr(cand.first - 1, cand.second - cand.first + 1);
      edlibFreeAlignResult(edlib_result);
      return alignment;
    }
  }
  edlibFreeAlignResult(edlib_result);

  std::string query_ssw = query;
  std::string target_ssw = target;
  std::transform(query_ssw.begin(), query_ssw.end(), query_ssw.begin(), ::toupper);
  std::transform(target_ssw.begin(), target_ssw.end(), target_ssw.begin(), ::toupper);

  int32_t mask_len = static_cast<int32_t>(query_ssw.size() / 2);
  if (mask_len < 15) {
    mask_len = 15;
  }

  StripedSmithWaterman::Aligner ssw_aligner(2, 2, 3, 2);
  StripedSmithWaterman::Filter ssw_filter;
  ssw_filter.report_cigar = true;
  StripedSmithWaterman::Alignment ssw_alignment;
  ssw_aligner.Align(
    target_ssw.c_str(),
    query_ssw.c_str(),
    static_cast<int>(query_ssw.size()),
    ssw_filter,
    &ssw_alignment,
    mask_len
  );

  if (ssw_alignment.query_begin < 0 || ssw_alignment.query_end < ssw_alignment.query_begin) {
    return alignment;
  }

  const int ssw_start = ssw_alignment.query_begin + 1;
  const int ssw_end = ssw_alignment.query_end + 1;
  const int ssw_length = ssw_end - ssw_start + 1;
  const int ssw_max_matches = cigar_max_consecutive_matches(ssw_alignment.cigar_string.c_str());
  const int ssw_max_matches_ind = cigar_max_consecutive_matches_with_indels(ssw_alignment.cigar_string.c_str());

  int deviation = 0;
  if (ssw_start < expected_start) {
    deviation += expected_start - ssw_start;
  }
  if (ssw_end > expected_end) {
    deviation += ssw_end - expected_end;
  }

  if (verbose) {
#pragma omp critical
    {
      std::ostringstream oss;
      oss << "SSW alignment: region [" << ssw_start << ", " << ssw_end << "], length = " << ssw_length
          << ", max_matches_with_indels = " << ssw_max_matches_ind
          << ", max_matches = " << ssw_max_matches
          << ", deviation = " << deviation
          << ", cigar = " << ssw_alignment.cigar_string << "\n";
      std::cout << oss.str();
    }
  }

  const bool pass_near =
    ssw_max_matches > min_match_bases &&
    deviation <= 100 &&
    ssw_length >= 10 &&
    ssw_max_matches_ind >= 8 &&
    ssw_max_matches >= 8;
  const bool pass_far =
    deviation > 100 &&
    ssw_length >= 10 &&
    ssw_max_matches_ind >= 10;

  if (!(pass_near || pass_far)) {
    return alignment;
  }

  alignment.positions.push_back({ssw_start, ssw_end});
  alignment.edit_distance = cigar_edit_distance(ssw_alignment.cigar_string.c_str());
  alignment.success = alignment.edit_distance <= max_edit_distance;
  alignment.cigar = ssw_alignment.cigar_string;
  alignment.seq = target.substr(ssw_start - 1, ssw_end - ssw_start + 1);
  return alignment;
}

inline std::vector<std::string> split_sig_parts(const std::string& sig) {
  std::vector<std::string> out;
  if (sig.empty()) {
    return out;
  }
  size_t start = 0;
  while (start < sig.size()) {
    size_t end = sig.find('|', start);
    if (end == std::string::npos) {
      out.push_back(sig.substr(start));
      break;
    }
    out.push_back(sig.substr(start, end - start));
    start = end + 1;
  }
  return out;
}

inline bool parse_sig_span(const std::string& part, int& start_pos, int& stop_pos) {
  size_t last_colon = part.find_last_of(':');
  if (last_colon == std::string::npos) {
    return false;
  }
  size_t second_last_colon = part.find_last_of(':', last_colon - 1);
  if (second_last_colon == std::string::npos) {
    return false;
  }
  try {
    start_pos = std::stoi(part.substr(second_last_colon + 1, last_colon - second_last_colon - 1));
    stop_pos = std::stoi(part.substr(last_colon + 1));
  } catch (...) {
    return false;
  }
  return true;
}

inline void mask_span(std::string& sequence, int start_pos, int stop_pos) {
  if (sequence.empty()) {
    return;
  }
  int start = std::max(1, start_pos);
  int stop = std::min(static_cast<int>(sequence.size()), stop_pos);
  if (start > stop) {
    return;
  }
  std::fill(sequence.begin() + (start - 1), sequence.begin() + stop, 'X');
}

}  // namespace rad_static

