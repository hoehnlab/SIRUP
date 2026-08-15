#' Cluster tree tips by cutting tree branches
#'
#' Partitions a phylogenetic tree into connected components by cutting `k - 1`
#' branches. A connected component is one disconnected tree piece left after
#' branch cuts; tips in the same component receive the same cluster label.
#'
#' With `method = "cut_largest_branches"`, the function cuts the longest
#' eligible branches globally. This is the graph-theoretic tree clustering
#' heuristic used for fixed-k minimum spanning tree clustering, applied to a
#' phylogenetic tree.
#'
#' With `method = "split_largest_subtrees"`, the function recursively splits the
#' current component with the largest genetic span, measured as the largest
#' patristic distance between any two tips in that component. Within that
#' component, candidate cuts must split it into two non-empty tip-containing
#' components. The selected cut maximizes the size of the smaller daughter
#' component, and branch length is used as a tie-breaker. Setting
#' `cut_terminal_edges = FALSE` makes the function prefer internal branches,
#' which usually avoids peeling off long single-tip branches before larger clades
#' have been split.
#'
#' More explicitly, `method = "split_largest_subtrees"` works as follows:
#' 1. Start with the full tree as one connected component.
#' 2. Repeat until `k - 1` branches have been cut, or no valid cut remains.
#' 3. Calculate the genetic span of each current component as the largest
#'    patristic distance between any two tips in that component.
#' 4. Choose the component with the largest genetic span, using tip count as a
#'    tie-breaker.
#' 5. Consider eligible branches inside that component. If
#'    `cut_terminal_edges = FALSE`, internal branches are tried first and
#'    terminal branches are used only as a fallback.
#' 6. Keep only candidate cuts that split the component into two non-empty
#'    tip-containing daughter components.
#' 7. Choose the candidate cut that maximizes the size of the smaller daughter
#'    component, using branch length as the tie-breaker.
#' 8. Cut that branch and repeat.
#' 9. Label the final connected components as `tree_cluster_1`,
#'    `tree_cluster_2`, and so on.
#'
#' @param tree A phylogenetic tree object with `edge`, `edge.length`, and
#'   `tip.label` fields.
#' @param k Number of connected tree components to form before ignored tips are
#'   removed. If larger than the number of tips in the tree, a warning is
#'   raised and `k` is set to the number of tips for that tree.
#' @param id Name of the returned tip identifier column.
#' @param cluster Name of the returned cluster column.
#' @param ignore_tips Optional tip labels to drop from the returned data frame,
#'   for example `"Germline"`.
#' @param prefix Prefix used for cluster labels.
#' @param method Tree clustering method. `"cut_largest_branches"` removes the
#'   `k - 1` longest eligible branches globally. `"split_largest_subtrees"`
#'   repeatedly splits the largest remaining tip-containing subtree.
#' @param cut_terminal_edges Cut terminal branches when choosing the longest
#'   branches to remove. Defaults to `TRUE`; set to `FALSE` to prioritize
#'   internal branches that split larger clades.
#' @return Data frame with one row per retained tip and a tree-cluster label.
#' @export
clusterTreeTips <- function(tree,
                            k,
                            id = "sequence_id",
                            cluster = "tree_cluster",
                            ignore_tips = "Germline",
                            prefix = "tree_cluster_",
                            method = c("split_largest_subtrees", "cut_largest_branches"),
                            cut_terminal_edges = TRUE) {
  method <- match.arg(method)
  if (!is.list(tree) || !all(c("edge", "edge.length", "tip.label") %in% names(tree))) {
    rlang::abort("`tree` must have `edge`, `edge.length`, and `tip.label` fields.")
  }
  if (!is.numeric(k) || length(k) != 1L || is.na(k) || k < 1) {
    rlang::abort("`k` must be a positive integer.")
  }
  k <- as.integer(k)
  n_tips <- length(tree$tip.label)
  if (n_tips == 0L) {
    rlang::abort("`tree` has no tips.")
  }
  if (k > n_tips) {
    warning(
      "`k` is larger than the number of tips in this tree; using k = ",
      n_tips,
      " instead.",
      call. = FALSE
    )
    k <- n_tips
  }

  edge <- as.matrix(tree$edge)
  edge_length <- as.numeric(tree$edge.length)
  if (nrow(edge) != length(edge_length)) {
    rlang::abort("`tree$edge` and `tree$edge.length` must have matching lengths.")
  }
  if (any(is.na(edge_length))) {
    rlang::abort("`tree$edge.length` contains missing values.")
  }

  keep_edge <- rep(TRUE, nrow(edge))
  if (k > 1L) {
    cut_edges <- .sirap_tree_cut_edges(
      edge = edge,
      edge_length = edge_length,
      n_tips = n_tips,
      n_cut = k - 1L,
      cut_terminal_edges = cut_terminal_edges,
      method = method
    )
    keep_edge[cut_edges] <- FALSE
  }

  nodes <- sort(unique(as.integer(edge)))
  adjacency <- stats::setNames(vector("list", length(nodes)), as.character(nodes))
  for (i in which(keep_edge)) {
    parent <- as.character(edge[i, 1])
    child <- as.character(edge[i, 2])
    adjacency[[parent]] <- c(adjacency[[parent]], child)
    adjacency[[child]] <- c(adjacency[[child]], parent)
  }

  visited <- stats::setNames(rep(FALSE, length(nodes)), as.character(nodes))
  component_id <- stats::setNames(rep(NA_integer_, length(nodes)), as.character(nodes))
  component <- 0L
  for (node in names(adjacency)) {
    if (visited[[node]]) {
      next
    }
    component <- component + 1L
    queue <- node
    visited[[node]] <- TRUE
    while (length(queue) > 0L) {
      current <- queue[[1]]
      queue <- queue[-1]
      component_id[[current]] <- component
      for (neighbor in adjacency[[current]]) {
        if (!visited[[neighbor]]) {
          visited[[neighbor]] <- TRUE
          queue <- c(queue, neighbor)
        }
      }
    }
  }

  tip_nodes <- as.character(seq_len(n_tips))
  out <- data.frame(
    .tip_id = tree$tip.label,
    .component = component_id[tip_nodes],
    stringsAsFactors = FALSE
  )
  if (!is.null(ignore_tips)) {
    out <- out[!out$.tip_id %in% ignore_tips, , drop = FALSE]
  }
  cluster_levels <- unique(out$.component)
  cluster_labels <- stats::setNames(
    paste0(prefix, seq_along(cluster_levels)),
    as.character(cluster_levels)
  )
  out[[id]] <- out$.tip_id
  out[[cluster]] <- cluster_labels[as.character(out$.component)]
  out[, c(id, cluster), drop = FALSE]
}

#' Choose k for tree tip clustering
#'
#' Evaluates fixed-k tree partitions from `clusterTreeTips()` and chooses a
#' useful number of connected tip clusters. The default rule picks the previous
#' k when the relative improvement in within-cluster patristic distance falls
#' below `improvement_threshold`.
#'
#' @param tree A phylogenetic tree object with `edge`, `edge.length`, and
#'   `tip.label` fields.
#' @param k_range Candidate k values.
#' @param improvement_threshold Relative improvement threshold. The previous
#'   `k` is chosen once the improvement in mean within-cluster patristic
#'   distance falls below this value.
#' @param ignore_tips Optional tip labels to drop before scoring, for example
#'   `"Germline"`.
#' @param method Tree clustering method passed to `clusterTreeTips()`.
#' @param cut_terminal_edges Cut terminal branches when choosing branches to
#'   remove, passed to `clusterTreeTips()`.
#' @return List with selected `k` and a candidate score data frame.
#' @export
chooseTreeTipK <- function(tree,
                           k_range = 1:10,
                           improvement_threshold = 0.10,
                           ignore_tips = "Germline",
                           method = c("cut_largest_branches", "split_largest_subtrees"),
                           cut_terminal_edges = TRUE) {
  method <- match.arg(method)
  retained_tips <- setdiff(tree$tip.label, ignore_tips)
  max_k <- length(retained_tips)
  k_range <- sort(unique(as.integer(k_range)))
  k_range <- k_range[!is.na(k_range) & k_range >= 1L & k_range <= max_k]
  if (length(k_range) == 0L) {
    rlang::abort("`k_range` does not contain any valid k values for this tree.")
  }

  scores <- do.call(rbind, lapply(k_range, function(k) {
    clusters <- clusterTreeTips(
      tree,
      k = k,
      id = "sequence_id",
      cluster = "tree_cluster",
      ignore_tips = ignore_tips,
      method = method,
      cut_terminal_edges = cut_terminal_edges
    )
    cluster_score <- .sirap_tree_cluster_score(
      tree,
      clusters,
      id_col = "sequence_id",
      cluster_col = "tree_cluster"
    )
    data.frame(
      k = k,
      score = cluster_score,
      stringsAsFactors = FALSE
    )
  }))

  scores$relative_improvement <- NA_real_
  if (nrow(scores) > 1L) {
    previous <- scores$score[-nrow(scores)]
    current <- scores$score[-1L]
    scores$relative_improvement[-1L] <- ifelse(
      previous == 0,
      0,
      (previous - current) / previous
    )
  }

  low_gain <- which(!is.na(scores$relative_improvement) & scores$relative_improvement < improvement_threshold)
  chosen_k <- if (length(low_gain) > 0L) {
    scores$k[max(1L, low_gain[[1]] - 1L)]
  } else {
    scores$k[nrow(scores)]
  }

  list(k = as.integer(chosen_k), scores = scores)
}

#' Cluster tips across a dowser trees object
#'
#' Applies `clusterTreeTips()` to each tree in a dowser `getTrees()` result and
#' adds the resulting cluster column back to each tree's clone data.
#'
#' @param trees Dowser `getTrees()` output, or a list containing a `trees`
#'   element with that output.
#' @param k Number of connected tree components to form for each tree, or
#'   `"auto"` to choose k separately for each tree with `chooseTreeTipK()`.
#' @param id Column in each clone data frame matching tree tip labels.
#' @param cluster Name of the cluster column to add.
#' @param ignore_tips Optional tip labels to drop before joining clusters.
#' @param prefix Prefix used for cluster labels.
#' @param keep Optional column in each clone data frame used to choose
#'   whether a tree has the expected retained tips. Clustering is run before
#'   any keep-filtering; tips where `keep != keep_val` then have their
#'   computed cluster labels overwritten with their `keep` value.
#' @param keep_val Value in `keep` that indicates a tip should be kept.
#' @param method Tree clustering method passed to `clusterTreeTips()`.
#' @param cut_terminal_edges Cut terminal branches when choosing branches to
#'   remove. Defaults to `FALSE`, which prioritizes internal clade splits; set
#'   to `TRUE` to allow terminal tips to be peeled off by branch cuts.
#' @param k_range Candidate k values passed to `chooseTreeTipK()` when
#'   `k = "auto"`.
#' @param improvement_threshold Relative improvement threshold passed to
#'   `chooseTreeTipK()` when `k = "auto"`.
#' @return `trees`, with `cluster` added to each clone data object.
#' @export
clusterTreesTips <- function(trees,
                             k = "auto",
                             id = "sequence_id",
                             cluster = "tree_cluster",
                             prefix = "tree_cluster_",
                             ignore_tips = "Germline",
                             keep = NULL,
                             keep_val = NULL,
                             method = c("split_largest_subtrees", "cut_largest_branches"),
                             cut_terminal_edges = FALSE,
                             k_range = 1:8,
                             improvement_threshold = 0.10) {
  method <- match.arg(method)
  auto_k <- identical(k, "auto")
  if (is.list(trees) && !is.data.frame(trees) && "trees" %in% names(trees)) {
    trees$trees <- clusterTreesTips(
      trees$trees,
      k = k,
      id = id,
      cluster = cluster,
      ignore_tips = ignore_tips,
      prefix = prefix,
      keep = keep,
      keep_val = keep_val,
      method = method,
      cut_terminal_edges = cut_terminal_edges,
      k_range = k_range,
      improvement_threshold = improvement_threshold
    )
    return(trees)
  }

  if (!is.data.frame(trees) || !all(c("trees", "data") %in% colnames(trees))) {
    rlang::abort("`trees` must be a dowser trees data frame or a list containing one in `trees`.")
  }

  out <- trees
  chosen_k <- integer(nrow(out))
  k_scores <- vector("list", nrow(out))
  use_keep_filter <- .sirap_check_tree_keep_args(keep, keep_val)
  for (i in seq_len(nrow(out))) {
    tree_for_clustering <- out$trees[[i]]
    extra_tip_clusters <- NULL
    if (use_keep_filter) {
      extra_tip_clusters <- .sirap_tree_non_keep_clusters(
        tree = out$trees[[i]],
        clone_data = out$data[[i]],
        id_col = id,
        keep_col = keep,
        keep_value = keep_val,
        cluster_col = cluster,
        ignore_tips = ignore_tips
      )
    }

    tree_k <- k
    if (auto_k) {
      k_choice <- chooseTreeTipK(
        tree_for_clustering,
        k_range = k_range,
        improvement_threshold = improvement_threshold,
        ignore_tips = ignore_tips,
        method = method,
        cut_terminal_edges = cut_terminal_edges
      )
      tree_k <- k_choice$k
      k_scores[[i]] <- k_choice$scores
    }
    chosen_k[[i]] <- as.integer(tree_k)
    if (!auto_k) {
      n_tree_tips <- length(tree_for_clustering$tip.label)
      if (as.integer(tree_k) > n_tree_tips) {
        warning(
          .sirap_tree_row_label(out, i),
          " requested k = ",
          as.integer(tree_k),
          " but has only ",
          n_tree_tips,
          " tips; using k = ",
          n_tree_tips,
          " for this tree.",
          call. = FALSE
        )
        tree_k <- n_tree_tips
        chosen_k[[i]] <- as.integer(tree_k)
      }
    }
    tree_tip_clusters <- clusterTreeTips(
      tree_for_clustering,
      k = tree_k,
      id = id,
      cluster = cluster,
      ignore_tips = ignore_tips,
      prefix = prefix,
      method = method,
      cut_terminal_edges = cut_terminal_edges
    )
    if (!is.null(extra_tip_clusters) && nrow(extra_tip_clusters) > 0L) {
      tree_tip_clusters <- tree_tip_clusters[
        !tree_tip_clusters[[id]] %in% extra_tip_clusters[[id]],
        ,
        drop = FALSE
      ]
      tree_tip_clusters <- rbind(tree_tip_clusters, extra_tip_clusters)
    }
    tree_tip_clusters[[cluster]] <- .sirap_compact_tree_cluster_labels(
      tree_tip_clusters[[cluster]],
      prefix = prefix
    )
    out$data[[i]] <- .sirap_add_tree_tip_clusters(
      out$data[[i]],
      tree_tip_clusters,
      id_col = id,
      cluster_col = cluster
    )
  }
  out[[paste0(cluster, "_k")]] <- chosen_k
  attr(out, paste0(cluster, "_scores")) <- k_scores
  out
}

#' Format SIRAP receptors as dowser clones
#'
#' Runs [dowser::formatClones()] on the SIRAP receptor table and stores the
#' result as a named analysis. The default germline column follows dowser's
#' default; if it is missing, add germlines to the receptor table first with
#' `updateReceptors()` and [dowser::createGermlines()].
#'
#' @param object A SIRAP object.
#' @param analysis_name Name used to store the formatted clone analysis.
#' @param seq dowser: name of the column containing observed DNA sequences. All
#'   sequences in this column must be multiple aligned.
#' @param clone dowser: name of the column containing the identifier for the
#'   clone. All entries in this column should be identical.
#' @param subgroup dowser: name of the column containing the identifier for the
#'   subgroup.
#' @param id dowser: name of the column containing sequence identifiers.
#' @param germ dowser: name of the column containing germline DNA sequences.
#'   All entries in this column should be identical for any given clone, and
#'   they must be multiple aligned with the data in the seq column.
#' @param v_call dowser: name of the column containing V-segment allele
#'   assignments. All entries in this column should be identical to the gene
#'   level.
#' @param j_call dowser: name of the column containing J-segment allele
#'   assignments. All entries in this column should be identical to the gene
#'   level.
#' @param junc_len dowser: name of the column containing the length of the
#'   junction as a numeric value. All entries in this column should be identical
#'   for any given clone.
#' @param mask_char dowser: character to use for masking and padding.
#' @param max_mask dowser: maximum number of characters to mask at the leading
#'   and trailing sequence ends. If NULL then the upper masking bound will be
#'   automatically determined from the maximum number of observed leading or
#'   trailing Ns amongst all sequences. If set to 0 (default) then masking will
#'   not be performed.
#' @param pad_end dowser: if TRUE pad the end of each sequence with mask_char to
#'   make every sequence the same length.
#' @param text_fields dowser: text annotation columns to retain and merge during
#'   duplicate removal.
#' @param num_fields dowser: numeric annotation columns to retain and sum during
#'   duplicate removal.
#' @param seq_fields dowser: sequence annotation columns to retain and collapse
#'   during duplicate removal. Note, this is distinct from the seq and germ
#'   arguments, which contain the primary sequence data for the clone and should
#'   not be repeated in this argument.
#' @param add_count dowser: if TRUE add an additional annotation column called
#'   COLLAPSE_COUNT during duplicate removal that indicates the number of
#'   sequences that were collapsed.
#' @param verbose dowser: passed on to collapseDuplicates. If TRUE, report the
#'   numbers of input, discarded and output sequences; otherwise, process
#'   sequences silently.
#' @param collapse dowser: collapse identical sequences?
#' @param cell dowser: name of the column containing cell assignment information
#' @param locus dowser: name of the column containing locus information
#' @param traits dowser: column ids to keep distinct during sequence collapse
#' @param mod3 dowser: pad sequences to length multiple three?
#' @param randomize dowser: randomize sequence order? Important if using PHYLIP
#' @param use_regions dowser: assign CDR/FWR regions?
#' @param dup_singles dowser: Duplicate sequences in singleton clones to include
#'   them as trees?
#' @param nproc dowser: number of cores to parallelize formatting over.
#' @param chain dowser: if HL, include light chain information if available.
#' @param heavy dowser: name of heavy chain locus (default = "IGH")
#' @param filterstop dowser: only use sequences that do not contain an in-frame
#'   stop codon
#' @param minseq dowser: minimum number of sequences per clone
#' @param split_light dowser: split or lump subgroups? See resolveLightChains.
#' @param light_traits dowser: Include the traits from the light chain when
#'   concatenating and collapsing trees?
#' @param majoronly dowser: only return largest subgroup and sequences without
#'   light chains
#' @param columns dowser: additional data columns to include in output
#' @param validate Whether to validate the updated SIRAP object.
#' @return An updated SIRAP object with the formatted clones stored as an analysis.
#' @export
formatSpatialClones <- function(object,
                                analysis_name = "f_clones",
                                seq = "sequence_alignment",
                                clone = "clone_id",
                                subgroup = "clone_subgroup",
                                id = "sequence_id",
                                germ = "germline_alignment_d_mask",
                                v_call = "v_call",
                                j_call = "j_call",
                                junc_len = "junction_length",
                                mask_char = "N",
                                max_mask = 0,
                                pad_end = TRUE,
                                text_fields = NULL,
                                num_fields = NULL,
                                seq_fields = NULL,
                                add_count = TRUE,
                                verbose = FALSE,
                                collapse = TRUE,
                                cell = "cell_id",
                                locus = "locus",
                                traits = NULL,
                                mod3 = TRUE,
                                randomize = TRUE,
                                use_regions = TRUE,
                                dup_singles = FALSE,
                                nproc = 1,
                                chain = "H",
                                heavy = "IGH",
                                filterstop = TRUE,
                                minseq = 2,
                                split_light = FALSE,
                                light_traits = FALSE,
                                majoronly = FALSE,
                                columns = NULL,
                                validate = TRUE) {
  .require_namespace("dowser")
  validateSIRAP(object)
  receptors <- getReceptors(object)
  if (!germ %in% colnames(receptors)) {
    rlang::abort(paste0(
      "The receptor table does not contain `", germ,
      "`. Run `updateReceptors()` on this SIRAP object with ",
      "`dowser::createGermlines()` before `formatSpatialClones()`."
    ))
  }
  # Carry a verbatim copy of the receptor id (and the sample column, when the
  # SIRAP `sample` field is set) through tree building. dowser sanitizes tip
  # labels for Newick output (":" and "," become "_"), so the tree `id` column
  # can no longer be matched against the raw receptor id. `.sirap_uid` is passed
  # as a text field, which dowser preserves verbatim (and comma-joins when it
  # collapses identical sequences), so `plotSpatialCloneTreePairs()` can match
  # tree tips back to this object's receptors. See that function for the match.
  uid_col <- ".sirap_uid"
  receptors[[uid_col]] <- as.character(receptors[[id]])
  text_fields <- union(text_fields, uid_col)
  sample_field <- .sirap_receptor_fields(object)$sample
  if (!is.null(sample_field) && sample_field %in% colnames(receptors)) {
    traits <- union(traits, sample_field)
  }
  f_clones <- dowser::formatClones(
    data = receptors,
    seq = seq,
    clone = clone,
    subgroup = subgroup,
    id = id,
    germ = germ,
    v_call = v_call,
    j_call = j_call,
    junc_len = junc_len,
    mask_char = mask_char,
    max_mask = max_mask,
    pad_end = pad_end,
    text_fields = text_fields,
    num_fields = num_fields,
    seq_fields = seq_fields,
    add_count = add_count,
    verbose = verbose,
    collapse = collapse,
    cell = cell,
    locus = locus,
    traits = traits,
    mod3 = mod3,
    randomize = randomize,
    use_regions = use_regions,
    dup_singles = dup_singles,
    nproc = nproc,
    chain = chain,
    heavy = heavy,
    filterstop = filterstop,
    minseq = minseq,
    split_light = split_light,
    light_traits = light_traits,
    majoronly = majoronly,
    columns = columns
  )
  setAnalysis(object, analysis_name, f_clones, validate = validate)
}

#' Infer dowser trees from formatted SIRAP clone analyses
#'
#' Runs [dowser::getTrees()] on a stored formatted-clone analysis and stores the
#' resulting trees as another SIRAP analysis.
#'
#' @param object A SIRAP object.
#' @param f_clones_analysis Name of the stored formatted-clone analysis.
#' @param analysis_name Name used to store the tree analysis.
#' @param trait dowser: trait to use for parsimony models (required if igphyml
#'   specified)
#' @param id dowser: unique identifier for this analysis (required if igphyml or
#'   dnapars specified)
#' @param dir dowser: directory where temporary files will be placed.
#' @param modelfile dowser: file specifying parsimony model to use
#' @param build dowser: program to use for tree building (pratchet, pml,
#'   dnapars, dnaml, igphyml, raxml)
#' @param exec dowser: location of desired phylogenetic executable
#' @param igphyml dowser: optional location of igphyml executable for parsimony
#' @param fixtrees dowser: if TRUE, use supplied tree topologies
#' @param nproc dowser: number of cores to parallelize computations
#' @param quiet dowser: amount of rubbish to print to console
#' @param rm_temp dowser: remove temporary files (default=TRUE)
#' @param palette dowser: deprecated
#' @param seq dowser: column name containing sequence information
#' @param collapse dowser: Collapse internal nodes with identical sequences?
#' @param check_divergence dowser: run checkDivergence on trees?
#' @param ... dowser: Additional arguments passed to tree building programs
#' @param validate Whether to validate the updated SIRAP object.
#' @return An updated SIRAP object with dowser trees stored as an analysis.
#' @export
getSpatialTrees <- function(object,
                            f_clones_analysis = "f_clones",
                            analysis_name = "trees",
                            trait = NULL,
                            id = NULL,
                            dir = NULL,
                            modelfile = NULL,
                            build = "pratchet",
                            exec = NULL,
                            igphyml = NULL,
                            fixtrees = FALSE,
                            nproc = 1,
                            quiet = 0,
                            rm_temp = TRUE,
                            palette = NULL,
                            seq = NULL,
                            collapse = FALSE,
                            check_divergence = TRUE,
                            ...,
                            validate = TRUE) {
  .require_namespace("dowser")
  validateSIRAP(object)
  clones <- getAnalysis(object, f_clones_analysis)
  trees <- dowser::getTrees(
    clones = clones,
    trait = trait,
    id = id,
    dir = dir,
    modelfile = modelfile,
    build = build,
    exec = exec,
    igphyml = igphyml,
    fixtrees = fixtrees,
    nproc = nproc,
    quiet = quiet,
    rm_temp = rm_temp,
    palette = palette,
    seq = seq,
    collapse = collapse,
    check_divergence = check_divergence,
    ...
  )
  setAnalysis(object, analysis_name, trees, validate = validate)
}

#' Assign stable colors within selected clones
#'
#' Adds a per-tip color identifier and color value that can be used consistently
#' in tree plots and matched spatial plots.
#'
#' @param clones Clone-assigned AIRR data frame.
#' @param clone_val Optional clone IDs to retain before assigning colors.
#' @param clone Column name for clone identifiers.
#' @param id Column name for sequence/tip identifiers.
#' @param shared_color Optional column whose shared values should receive
#'   the same tip color within each clone, for example a spatial ID column.
#' @param palette Optional color vector. Defaults to `pals::glasbey()` when
#'   available, otherwise base R qualitative colors.
#' @return Data frame with `tip_color_id` and `tip_color` columns.
#' @export
assignTipColors <- function(clones,
                              clone_val = NULL,
                              clone = "clone_id",
                              id = "sequence_id",
                              shared_color = NULL,
                              palette = NULL) {
  .check_required_cols(clones, c(clone, id), "clones")
  if (!is.null(shared_color)) {
    .check_required_cols(clones, shared_color, "clones")
  }
  out <- .as_data_frame(clones)
  if (!is.null(clone_val)) {
    out <- out[out[[clone]] %in% clone_val, , drop = FALSE]
  }
  out <- out[!is.na(out[[clone]]) & !is.na(out[[id]]), , drop = FALSE]
  if (is.null(shared_color)) {
    out <- out %>%
      dplyr::group_by(dplyr::across(dplyr::all_of(clone))) %>%
      dplyr::arrange(.data[[id]], .by_group = TRUE) %>%
      dplyr::mutate(tip_color_id = as.character(dplyr::row_number())) %>%
      dplyr::ungroup()
  } else {
    out <- out %>%
      dplyr::group_by(dplyr::across(dplyr::all_of(clone))) %>%
      dplyr::arrange(.data[[shared_color]], .data[[id]], .by_group = TRUE) %>%
      dplyr::mutate(tip_color_id = as.character(dplyr::dense_rank(.data[[shared_color]]))) %>%
      dplyr::ungroup()
  }

  n_colors <- max(as.integer(out$tip_color_id), na.rm = TRUE)
  if (is.null(palette)) {
    if (requireNamespace("pals", quietly = TRUE)) {
      palette <- .sirap_extend_palette(pals::glasbey(32), n_colors)
    } else {
      palette <- grDevices::hcl.colors(max(12, n_colors), palette = "Dark 3")
    }
  } else if (length(palette) < n_colors) {
    palette <- .sirap_extend_palette(palette, n_colors)
  }
  out$tip_color <- palette[as.integer(out$tip_color_id)]
  out
}

# Dispatch to the requested strategy for choosing which tree edges to remove.
.sirap_tree_cut_edges <- function(edge,
                                  edge_length,
                                  n_tips,
                                  n_cut,
                                  cut_terminal_edges,
                                  method) {
  if (identical(method, "cut_largest_branches")) {
    return(.sirap_tree_cut_largest_branches(
      edge = edge,
      edge_length = edge_length,
      n_tips = n_tips,
      n_cut = n_cut,
      cut_terminal_edges = cut_terminal_edges
    ))
  }

  .sirap_tree_split_largest_subtrees(
    edge = edge,
    edge_length = edge_length,
    n_tips = n_tips,
    n_cut = n_cut,
    cut_terminal_edges = cut_terminal_edges
  )
}

# Choose the globally longest eligible branches for fixed-k tree clustering.
.sirap_tree_cut_largest_branches <- function(edge,
                                             edge_length,
                                             n_tips,
                                             n_cut,
                                             cut_terminal_edges) {
  edge_order <- order(edge_length, decreasing = TRUE)
  if (isTRUE(cut_terminal_edges)) {
    return(edge_order[seq_len(n_cut)])
  }

  internal_edges <- edge_order[edge[edge_order, 2] > n_tips]
  terminal_edges <- edge_order[edge[edge_order, 2] <= n_tips]
  c(internal_edges, terminal_edges)[seq_len(n_cut)]
}

# Recursively split the current component with the largest tip-to-tip span.
.sirap_tree_split_largest_subtrees <- function(edge,
                                               edge_length,
                                               n_tips,
                                               n_cut,
                                               cut_terminal_edges) {
  keep_edge <- rep(TRUE, nrow(edge))
  cut_edges <- integer()
  nodes <- sort(unique(as.integer(edge)))

  for (step in seq_len(n_cut)) {
    components <- .sirap_tree_components(edge, keep_edge, nodes)
    tip_counts <- vapply(
      components,
      function(component) sum(component <= n_tips),
      integer(1)
    )
    weighted_adjacency <- .sirap_tree_weighted_adjacency(edge, edge_length, keep_edge, nodes)
    component_scores <- vapply(
      components,
      .sirap_tree_component_diameter,
      numeric(1),
      n_tips = n_tips,
      weighted_adjacency = weighted_adjacency
    )
    component_order <- order(component_scores, tip_counts, decreasing = TRUE)

    best_edge <- NA_integer_
    for (component_idx in component_order) {
      if (tip_counts[[component_idx]] < 2L) {
        next
      }

      component_nodes <- components[[component_idx]]
      component_edges <- which(
        keep_edge &
          edge[, 1] %in% component_nodes &
          edge[, 2] %in% component_nodes
      )

      candidate_edges <- if (isTRUE(cut_terminal_edges)) {
        component_edges
      } else {
        component_edges[edge[component_edges, 2] > n_tips]
      }
      best_edge <- .sirap_best_tree_cut_edge(
        candidate_edges,
        edge = edge,
        edge_length = edge_length,
        keep_edge = keep_edge,
        nodes = nodes,
        n_tips = n_tips,
        component_nodes = component_nodes
      )

      if (is.na(best_edge) && !isTRUE(cut_terminal_edges)) {
        terminal_edges <- component_edges[edge[component_edges, 2] <= n_tips]
        best_edge <- .sirap_best_tree_cut_edge(
          terminal_edges,
          edge = edge,
          edge_length = edge_length,
          keep_edge = keep_edge,
          nodes = nodes,
          n_tips = n_tips,
          component_nodes = component_nodes
        )
      }

      if (!is.na(best_edge)) {
        break
      }
    }

    if (is.na(best_edge)) {
      break
    }
    keep_edge[[best_edge]] <- FALSE
    cut_edges <- c(cut_edges, best_edge)
  }

  cut_edges
}

# Score candidate branch cuts and choose the most balanced split.
.sirap_best_tree_cut_edge <- function(candidate_edges,
                                      edge,
                                      edge_length,
                                      keep_edge,
                                      nodes,
                                      n_tips,
                                      component_nodes) {
  if (length(candidate_edges) == 0L) {
    return(NA_integer_)
  }

  candidate_scores <- lapply(candidate_edges, function(edge_idx) {
    keep_after_cut <- keep_edge
    keep_after_cut[[edge_idx]] <- FALSE
    components <- .sirap_tree_components(edge, keep_after_cut, nodes)
    split_components <- components[
      vapply(components, function(component) all(component %in% component_nodes), logical(1))
    ]
    tip_counts <- vapply(
      split_components,
      function(component) sum(component <= n_tips),
      integer(1)
    )
    if (length(tip_counts) < 2L || any(tip_counts == 0L)) {
      return(NULL)
    }
    data.frame(
      edge_idx = edge_idx,
      min_tip_count = min(tip_counts),
      edge_length = edge_length[[edge_idx]],
      stringsAsFactors = FALSE
    )
  })

  candidate_scores <- do.call(rbind, candidate_scores[!vapply(candidate_scores, is.null, logical(1))])
  if (is.null(candidate_scores) || nrow(candidate_scores) == 0L) {
    return(NA_integer_)
  }
  candidate_scores <- candidate_scores[
    order(candidate_scores$min_tip_count, candidate_scores$edge_length, decreasing = TRUE),
    ,
    drop = FALSE
  ]
  candidate_scores$edge_idx[[1]]
}

# Return connected components after dropping selected tree edges.
.sirap_tree_components <- function(edge, keep_edge, nodes) {
  kept_edge <- edge[keep_edge, , drop = FALSE]
  adjacency <- stats::setNames(vector("list", length(nodes)), as.character(nodes))
  if (nrow(kept_edge) > 0L) {
    for (i in seq_len(nrow(kept_edge))) {
      parent <- as.character(kept_edge[i, 1])
      child <- as.character(kept_edge[i, 2])
      adjacency[[parent]] <- c(adjacency[[parent]], child)
      adjacency[[child]] <- c(adjacency[[child]], parent)
    }
  }

  visited <- stats::setNames(rep(FALSE, length(nodes)), as.character(nodes))
  components <- list()
  for (node in names(adjacency)) {
    if (visited[[node]]) {
      next
    }
    queue <- node
    visited[[node]] <- TRUE
    component <- integer()
    while (length(queue) > 0L) {
      current <- queue[[1]]
      queue <- queue[-1]
      component <- c(component, as.integer(current))
      for (neighbor in adjacency[[current]]) {
        if (!visited[[neighbor]]) {
          visited[[neighbor]] <- TRUE
          queue <- c(queue, neighbor)
        }
      }
    }
    components[[length(components) + 1L]] <- component
  }
  components
}

# Build a branch-length-weighted adjacency list for retained tree edges.
.sirap_tree_weighted_adjacency <- function(edge, edge_length, keep_edge, nodes) {
  kept_idx <- which(keep_edge)
  adjacency <- stats::setNames(vector("list", length(nodes)), as.character(nodes))
  for (i in kept_idx) {
    parent <- as.character(edge[i, 1])
    child <- as.character(edge[i, 2])
    weight <- edge_length[[i]]
    adjacency[[parent]] <- c(adjacency[[parent]], stats::setNames(weight, child))
    adjacency[[child]] <- c(adjacency[[child]], stats::setNames(weight, parent))
  }
  adjacency
}

# Measure a component's genetic span as its largest tip-to-tip distance.
.sirap_tree_component_diameter <- function(component, n_tips, weighted_adjacency) {
  tips <- component[component <= n_tips]
  if (length(tips) < 2L) {
    return(0)
  }

  max_distance <- 0
  for (tip in tips) {
    distances <- .sirap_tree_distances_from_node(tip, weighted_adjacency)
    tip_distances <- distances[as.character(tips)]
    max_distance <- max(max_distance, tip_distances, na.rm = TRUE)
  }
  max_distance
}

# Compute shortest branch-length distances from one tree node.
.sirap_tree_distances_from_node <- function(start, weighted_adjacency) {
  start <- as.character(start)
  distances <- stats::setNames(rep(NA_real_, length(weighted_adjacency)), names(weighted_adjacency))
  distances[[start]] <- 0
  queue <- start
  while (length(queue) > 0L) {
    current <- queue[[1]]
    queue <- queue[-1]
    neighbors <- weighted_adjacency[[current]]
    for (neighbor in names(neighbors)) {
      if (is.na(distances[[neighbor]])) {
        distances[[neighbor]] <- distances[[current]] + unname(neighbors[[neighbor]])
        queue <- c(queue, neighbor)
      }
    }
  }
  distances
}

# Score clustering quality for one candidate k in automatic k selection.
.sirap_tree_cluster_score <- function(tree,
                                      clusters,
                                      id_col,
                                      cluster_col) {
  dist_mat <- .sirap_tip_patristic_matrix(tree, clusters[[id_col]])
  total <- 0
  n_pairs <- 0L
  for (cluster in unique(clusters[[cluster_col]])) {
    ids <- clusters[[id_col]][clusters[[cluster_col]] == cluster]
    if (length(ids) < 2L) {
      next
    }
    sub <- dist_mat[ids, ids, drop = FALSE]
    values <- sub[upper.tri(sub, diag = FALSE)]
    total <- total + sum(values, na.rm = TRUE)
    n_pairs <- n_pairs + length(values)
  }
  if (n_pairs == 0L) {
    0
  } else {
    total / n_pairs
  }
}

#' Calculate patristic distances among selected tree tips
#'
#' Returns a pairwise patristic-distance matrix for selected tip labels,
#' preserving the requested tip order. This is the public wrapper around
#' SIRAP's tree-tip distance helper used by tree-tip clustering workflows.
#'
#' @param tree A phylogenetic tree with `edge`, `edge.length`, and `tip.label`.
#' @param tip_labels Character vector of tip labels to include.
#' @return A numeric matrix of pairwise patristic distances.
#' @export
tipPatristicMatrix <- function(tree, tip_labels) {
  .sirap_tip_patristic_matrix(tree = tree, tip_labels = tip_labels)
}

# Return patristic distances among selected tips, preserving tip order.
.sirap_tip_patristic_matrix <- function(tree, tip_labels) {
  edge <- as.matrix(tree$edge)
  edge_length <- as.numeric(tree$edge.length)
  n_tips <- length(tree$tip.label)
  nodes <- sort(unique(as.integer(edge)))
  parent <- stats::setNames(rep(NA_integer_, length(nodes)), as.character(nodes))
  branch_length <- stats::setNames(rep(0, length(nodes)), as.character(nodes))
  children <- edge[, 2]
  parent[as.character(children)] <- edge[, 1]
  branch_length[as.character(children)] <- edge_length
  root <- setdiff(edge[, 1], edge[, 2])[[1]]

  depth <- stats::setNames(rep(NA_real_, length(nodes)), as.character(nodes))
  depth[as.character(root)] <- 0
  for (node in names(depth)) {
    .sirap_node_depth(as.integer(node), parent, branch_length, depth)
  }
  for (node in names(depth)) {
    depth[[node]] <- .sirap_node_depth(as.integer(node), parent, branch_length, depth)
  }

  tip_nodes <- stats::setNames(seq_len(n_tips), tree$tip.label)
  tip_nodes <- tip_nodes[tip_labels]
  out <- matrix(0, nrow = length(tip_labels), ncol = length(tip_labels), dimnames = list(tip_labels, tip_labels))
  ancestors <- lapply(tip_nodes, .sirap_node_ancestors, parent = parent)
  names(ancestors) <- tip_labels
  for (i in seq_along(tip_labels)) {
    for (j in seq_along(tip_labels)) {
      if (j <= i) {
        next
      }
      shared <- intersect(ancestors[[i]], ancestors[[j]])
      lca <- shared[which.max(depth[as.character(shared)])]
      distance <- depth[as.character(tip_nodes[[i]])] + depth[as.character(tip_nodes[[j]])] - 2 * depth[as.character(lca)]
      out[i, j] <- distance
      out[j, i] <- distance
    }
  }
  out
}

# Calculate root-to-node depth from parent and branch-length maps.
.sirap_node_depth <- function(node, parent, branch_length, depth) {
  node_name <- as.character(node)
  if (!is.na(depth[[node_name]])) {
    return(depth[[node_name]])
  }
  parent_node <- parent[[node_name]]
  if (is.na(parent_node)) {
    return(0)
  }
  .sirap_node_depth(parent_node, parent, branch_length, depth) + branch_length[[node_name]]
}

# Walk from a node to the root and return its ancestor chain.
.sirap_node_ancestors <- function(node, parent) {
  ancestors <- integer()
  current <- node
  while (!is.na(current)) {
    ancestors <- c(ancestors, current)
    current <- parent[[as.character(current)]]
  }
  ancestors
}

# Relabel arbitrary component IDs as compact, stable cluster labels.
.sirap_compact_tree_cluster_labels <- function(labels, prefix) {
  labels <- as.character(labels)
  prefix_re <- gsub("([][{}()+*^$|\\\\?.])", "\\\\\\1", prefix)
  pattern <- paste0("^", prefix_re, "[0-9]+$")
  tree_labels <- unique(labels[!is.na(labels) & grepl(pattern, labels)])
  if (length(tree_labels) == 0L) {
    return(labels)
  }
  ids <- as.integer(sub(pattern = paste0("^", prefix_re), replacement = "", tree_labels))
  tree_labels <- tree_labels[order(ids)]
  label_map <- stats::setNames(paste0(prefix, seq_along(tree_labels)), tree_labels)
  replace <- !is.na(labels) & labels %in% names(label_map)
  labels[replace] <- label_map[labels[replace]]
  labels
}

# Build a human-readable label for a row in a dowser-style tree table.
.sirap_tree_row_label <- function(trees, i) {
  label_cols <- intersect(
    c("clone_id", "clone", "clone_id_rank", "sample_id", "sample", "locus"),
    colnames(trees)
  )
  for (col in label_cols) {
    value <- trees[[col]][[i]]
    if (length(value) == 1L && !is.na(value) && nzchar(as.character(value))) {
      return(paste0("Tree ", col, "=", value, " at row ", i))
    }
  }
  paste0("Tree at row ", i)
}

# Validate the paired arguments used to ignore reference/germline tips.
.sirap_check_tree_keep_args <- function(keep_col, keep_value) {
  has_keep_col <- !is.null(keep_col)
  has_keep_value <- !is.null(keep_value)
  if (xor(has_keep_col, has_keep_value)) {
    rlang::abort("`keep_col` and `keep_value` must be supplied together.")
  }
  if (!has_keep_col) {
    return(FALSE)
  }
  if (!is.character(keep_col) || length(keep_col) != 1L || is.na(keep_col) || keep_col == "") {
    rlang::abort("`keep_col` must be a single non-empty column name.")
  }
  if (length(keep_value) != 1L || is.na(keep_value)) {
    rlang::abort("`keep_value` must be a single non-missing value.")
  }
  TRUE
}

# Identify clusters that contain no retained non-reference tips.
.sirap_tree_non_keep_clusters <- function(tree,
                                          clone_data,
                                          id_col,
                                          keep_col,
                                          keep_value,
                                          cluster_col,
                                          ignore_tips) {
  data <- .sirap_tree_clone_data_frame(clone_data)
  .check_required_cols(data, c(id_col, keep_col), "tree clone data")

  data[[id_col]] <- as.character(data[[id_col]])
  data[[keep_col]] <- as.character(data[[keep_col]])
  tip_labels <- setdiff(as.character(tree$tip.label), ignore_tips)
  tip_data <- data[data[[id_col]] %in% tip_labels, , drop = FALSE]
  tip_data <- tip_data[!duplicated(tip_data[[id_col]]), , drop = FALSE]

  keep_value <- as.character(keep_value)
  keep_tip_ids <- tip_data[[id_col]][!is.na(tip_data[[keep_col]]) & tip_data[[keep_col]] == keep_value]
  if (length(keep_tip_ids) == 0L) {
    rlang::abort(paste0(
      "No tree tips had `", keep_col, " == ", keep_value, "` for clustering."
    ))
  }

  non_keep_tip_data <- tip_data[
    is.na(tip_data[[keep_col]]) | tip_data[[keep_col]] != keep_value,
    c(id_col, keep_col),
    drop = FALSE
  ]
  names(non_keep_tip_data) <- c(id_col, cluster_col)
  non_keep_tip_data[[cluster_col]] <- as.character(non_keep_tip_data[[cluster_col]])
  non_keep_tip_data
}

# Extract the clone metadata data frame from common dowser tree containers.
.sirap_tree_clone_data_frame <- function(clone_data) {
  if (inherits(clone_data, "airrClone")) {
    return(.as_data_frame(clone_data@data))
  }
  if (is.list(clone_data) && "data" %in% names(clone_data) && is.data.frame(clone_data$data)) {
    return(.as_data_frame(clone_data$data))
  }
  if (is.data.frame(clone_data)) {
    return(.as_data_frame(clone_data))
  }
  rlang::abort("Could not find a clone data frame in `trees$data`.")
}

# Add computed tip cluster labels back into the original clone-data container.
.sirap_add_tree_tip_clusters <- function(clone_data,
                                         tree_tip_clusters,
                                         id_col,
                                         cluster_col) {
  add_clusters <- function(data) {
    .check_required_cols(data, id_col, "tree clone data")
    data <- .as_data_frame(data)
    if (cluster_col %in% colnames(data)) {
      data <- data[, colnames(data) != cluster_col, drop = FALSE]
    }
    data %>% dplyr::left_join(tree_tip_clusters, by = id_col)
  }

  if (inherits(clone_data, "airrClone")) {
    clone_data@data <- add_clusters(clone_data@data)
    return(clone_data)
  }

  if (is.list(clone_data) && "data" %in% names(clone_data) && is.data.frame(clone_data$data)) {
    clone_data$data <- add_clusters(clone_data$data)
    return(clone_data)
  }

  if (is.data.frame(clone_data)) {
    return(add_clusters(clone_data))
  }

  rlang::abort("Could not find a clone data frame to add tree tip clusters to.")
}

# Extend a user palette with hue colors when there are more groups than colors.
.sirap_extend_palette <- function(palette, n) {
  palette <- unique(as.character(stats::na.omit(palette)))
  if (length(palette) == 0L) {
    palette <- grDevices::hcl.colors(max(2L, n), palette = "Dark 3")
  } else if (length(palette) == 1L && n > 1L) {
    palette <- unique(c(palette, grDevices::hcl.colors(max(2L, n), palette = "Dark 3")))
  }
  if (length(palette) >= n) {
    return(palette[seq_len(n)])
  }

  rgb <- grDevices::col2rgb(palette)
  averaged <- character()
  for (i in seq_len(ncol(rgb))) {
    for (j in seq_len(ncol(rgb))) {
      if (i >= j) {
        next
      }
      avg <- round(rowMeans(cbind(rgb[, i], rgb[, j])))
      averaged <- c(averaged, grDevices::rgb(avg[1], avg[2], avg[3], maxColorValue = 255))
      if (length(unique(c(palette, averaged))) >= n) {
        return(unique(c(palette, averaged))[seq_len(n)])
      }
    }
  }

  grDevices::colorRampPalette(palette)(n)
}
