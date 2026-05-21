// ============================================================================
// rad_conversion.h - Conversion utilities between rad C++ types and R data.table
// ============================================================================
#pragma once

#include <Rcpp.h>
#include "include/rad/rad_headers.h"

namespace rad_conversion {

// ============================================================================
// data.table creation helper
// Creates a proper data.table from column vectors
// ============================================================================
inline SEXP make_datatable(Rcpp::List columns) {
    // Set data.frame attributes first
    columns.attr("class") = Rcpp::CharacterVector::create("data.table", "data.frame");

    // Get number of rows from first column
    if (columns.size() > 0) {
        SEXP first_col = columns[0];
        R_xlen_t nrows = Rf_xlength(first_col);
        columns.attr("row.names") = Rcpp::seq(1, nrows);
    } else {
        columns.attr("row.names") = Rcpp::IntegerVector::create();
    }

    // Call data.table::setDT to properly initialize internal structure
    Rcpp::Environment dt_env = Rcpp::Environment::namespace_env("data.table");
    Rcpp::Function setDT = dt_env["setDT"];
    setDT(columns);

    return columns;
}

// ============================================================================
// SigElement (multi-index container) to data.table
// ============================================================================
inline SEXP sigelement_to_dt(const SigElement& sig) {
    size_t n = sig.size();

    Rcpp::CharacterVector class_ids(n);
    Rcpp::CharacterVector global_classes(n);
    Rcpp::CharacterVector types(n);
    Rcpp::CharacterVector directions(n);
    Rcpp::IntegerVector orders(n);
    Rcpp::IntegerVector pos_starts(n);
    Rcpp::IntegerVector pos_stops(n);
    Rcpp::IntegerVector edit_distances(n);
    Rcpp::LogicalVector element_pass(n);
    Rcpp::CharacterVector seqs(n);
    Rcpp::CharacterVector quals(n);

    size_t i = 0;
    for (const auto& elem : sig) {
        class_ids[i] = elem.class_id;
        global_classes[i] = elem.global_class;
        types[i] = elem.type;
        directions[i] = elem.direction;
        orders[i] = elem.order;
        pos_starts[i] = elem.position.first;
        pos_stops[i] = elem.position.second;

        if (elem.edit_distance) {
            edit_distances[i] = *elem.edit_distance;
        } else {
            edit_distances[i] = NA_INTEGER;
        }

        if (elem.element_pass) {
            element_pass[i] = *elem.element_pass ? TRUE : FALSE;
        } else {
            element_pass[i] = NA_LOGICAL;
        }

        seqs[i] = elem.seq.value_or("");
        quals[i] = elem.qual.value_or("");

        ++i;
    }

    Rcpp::List dt = Rcpp::List::create(
        Rcpp::Named("class_id") = class_ids,
        Rcpp::Named("global_class") = global_classes,
        Rcpp::Named("type") = types,
        Rcpp::Named("direction") = directions,
        Rcpp::Named("order") = orders,
        Rcpp::Named("position_start") = pos_starts,
        Rcpp::Named("position_stop") = pos_stops,
        Rcpp::Named("edit_distance") = edit_distances,
        Rcpp::Named("element_pass") = element_pass,
        Rcpp::Named("seq") = seqs,
        Rcpp::Named("qual") = quals
    );

    return make_datatable(dt);
}

// ============================================================================
// Layout_Struct (multi-index container) to data.table
// ============================================================================
inline SEXP layout_to_dt(const Layout_Struct& layout) {
    size_t n = layout.size();

    Rcpp::CharacterVector class_ids(n);
    Rcpp::CharacterVector seqs(n);
    Rcpp::CharacterVector masked_seqs(n);
    Rcpp::IntegerVector expected_lengths(n);
    Rcpp::CharacterVector types(n);
    Rcpp::IntegerVector orders(n);
    Rcpp::CharacterVector directions(n);
    Rcpp::LogicalVector unidirectional(n);
    Rcpp::CharacterVector global_classes(n);
    Rcpp::CharacterVector whitelist_paths(n);

    size_t i = 0;
    for (const auto& elem : layout) {
        class_ids[i] = elem.class_id;
        seqs[i] = elem.seq;
        masked_seqs[i] = elem.masked_seq;

        if (elem.expected_length) {
            expected_lengths[i] = *elem.expected_length;
        } else {
            expected_lengths[i] = NA_INTEGER;
        }

        types[i] = elem.type;
        orders[i] = elem.order;
        directions[i] = elem.direction;
        unidirectional[i] = elem.unidirectional ? TRUE : FALSE;
        global_classes[i] = elem.global_class;
        whitelist_paths[i] = elem.whitelist_path;

        ++i;
    }

    Rcpp::List dt = Rcpp::List::create(
        Rcpp::Named("class_id") = class_ids,
        Rcpp::Named("seq") = seqs,
        Rcpp::Named("masked_seq") = masked_seqs,
        Rcpp::Named("expected_length") = expected_lengths,
        Rcpp::Named("type") = types,
        Rcpp::Named("order") = orders,
        Rcpp::Named("direction") = directions,
        Rcpp::Named("unidirectional") = unidirectional,
        Rcpp::Named("global_class") = global_classes,
        Rcpp::Named("whitelist_path") = whitelist_paths
    );

    return make_datatable(dt);
}

// ============================================================================
// Sequences data.table (id, seq, qual) creation
// ============================================================================
inline SEXP sequences_to_dt(
    const std::vector<std::string>& ids,
    const std::vector<std::string>& seqs,
    const std::vector<std::string>& quals = {}
) {
    size_t n = ids.size();

    Rcpp::CharacterVector r_ids(n);
    Rcpp::CharacterVector r_seqs(n);

    for (size_t i = 0; i < n; ++i) {
        r_ids[i] = ids[i];
        r_seqs[i] = seqs[i];
    }

    Rcpp::List dt;
    if (quals.empty()) {
        dt = Rcpp::List::create(
            Rcpp::Named("id") = r_ids,
            Rcpp::Named("seq") = r_seqs
        );
    } else {
        Rcpp::CharacterVector r_quals(n);
        for (size_t i = 0; i < n; ++i) {
            r_quals[i] = quals[i];
        }
        dt = Rcpp::List::create(
            Rcpp::Named("id") = r_ids,
            Rcpp::Named("seq") = r_seqs,
            Rcpp::Named("qual") = r_quals
        );
    }

    return make_datatable(dt);
}

// ============================================================================
// read_streaming::sequence vector to data.table
// ============================================================================
inline SEXP read_sequences_to_dt(const std::vector<read_streaming::sequence>& sequences) {
    size_t n = sequences.size();

    Rcpp::CharacterVector ids(n);
    Rcpp::CharacterVector seqs(n);
    Rcpp::CharacterVector quals(n);
    Rcpp::CharacterVector comments(n);

    bool has_qual = false;

    for (size_t i = 0; i < n; ++i) {
        ids[i] = sequences[i].id;
        seqs[i] = sequences[i].seq;
        if (!sequences[i].qual.empty()) {
            quals[i] = sequences[i].qual;
            has_qual = true;
        } else {
            quals[i] = NA_STRING;
        }
        comments[i] = sequences[i].comment;
    }

    Rcpp::List dt;
    if (has_qual) {
        dt = Rcpp::List::create(
            Rcpp::Named("id") = ids,
            Rcpp::Named("seq") = seqs,
            Rcpp::Named("qual") = quals,
            Rcpp::Named("comment") = comments
        );
    } else {
        dt = Rcpp::List::create(
            Rcpp::Named("id") = ids,
            Rcpp::Named("seq") = seqs,
            Rcpp::Named("comment") = comments
        );
    }

    return make_datatable(dt);
}

// ============================================================================
// Barcode correction results to data.table
// ============================================================================
inline SEXP barcode_correction_to_dt(
    const std::vector<std::string>& original,
    const std::vector<std::string>& corrected,
    const std::vector<int>& distances,
    const std::vector<int>& counts
) {
    size_t n = original.size();

    Rcpp::CharacterVector r_original(n);
    Rcpp::CharacterVector r_corrected(n);
    Rcpp::IntegerVector r_distances(n);
    Rcpp::IntegerVector r_counts(n);

    for (size_t i = 0; i < n; ++i) {
        r_original[i] = original[i];
        r_corrected[i] = corrected[i];
        r_distances[i] = distances[i];
        r_counts[i] = counts[i];
    }

    Rcpp::List dt = Rcpp::List::create(
        Rcpp::Named("original") = r_original,
        Rcpp::Named("corrected") = r_corrected,
        Rcpp::Named("distance") = r_distances,
        Rcpp::Named("count") = r_counts
    );

    return make_datatable(dt);
}

// ============================================================================
// Mutation results to data.table
// ============================================================================
inline SEXP mutations_to_dt(
    const std::vector<std::string>& original,
    const std::vector<std::string>& mutated,
    const std::vector<int>& rounds
) {
    size_t n = original.size();

    Rcpp::CharacterVector r_original(n);
    Rcpp::CharacterVector r_mutated(n);
    Rcpp::IntegerVector r_rounds(n);

    for (size_t i = 0; i < n; ++i) {
        r_original[i] = original[i];
        r_mutated[i] = mutated[i];
        r_rounds[i] = rounds[i];
    }

    Rcpp::List dt = Rcpp::List::create(
        Rcpp::Named("original") = r_original,
        Rcpp::Named("mutated") = r_mutated,
        Rcpp::Named("round") = r_rounds
    );

    return make_datatable(dt);
}

// ============================================================================
// Levenshtein distance results to data.table
// ============================================================================
inline SEXP levenshtein_to_dt(
    const std::vector<std::string>& queries,
    const std::vector<std::string>& targets,
    const std::vector<int>& distances
) {
    size_t n = queries.size();

    Rcpp::CharacterVector r_queries(n);
    Rcpp::CharacterVector r_targets(n);
    Rcpp::IntegerVector r_distances(n);

    for (size_t i = 0; i < n; ++i) {
        r_queries[i] = queries[i];
        r_targets[i] = targets[i];
        r_distances[i] = distances[i];
    }

    Rcpp::List dt = Rcpp::List::create(
        Rcpp::Named("query") = r_queries,
        Rcpp::Named("target") = r_targets,
        Rcpp::Named("distance") = r_distances
    );

    return make_datatable(dt);
}

// ============================================================================
// Extract columns from data.table
// ============================================================================
inline std::vector<std::string> dt_get_string_column(SEXP dt, const std::string& col_name) {
    Rcpp::List dt_list(dt);
    Rcpp::CharacterVector col = dt_list[col_name];

    std::vector<std::string> result(col.size());
    for (R_xlen_t i = 0; i < col.size(); ++i) {
        result[i] = Rcpp::as<std::string>(col[i]);
    }
    return result;
}

inline std::vector<int> dt_get_int_column(SEXP dt, const std::string& col_name) {
    Rcpp::List dt_list(dt);
    Rcpp::IntegerVector col = dt_list[col_name];

    std::vector<int> result(col.size());
    for (R_xlen_t i = 0; i < col.size(); ++i) {
        result[i] = col[i];
    }
    return result;
}

} // namespace rad_conversion
