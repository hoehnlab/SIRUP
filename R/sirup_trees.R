#' Cluster the tips of a tree
#'
#' `clusterTreeTips` groups the tips of a phylogenetic tree into `k` connected
#' clusters.
#'
#' @details
#' The function groups the tips of the tree into `k` connected clusters. This
#' helps with visualization: tips in the same cluster get the same label, and you
#' can give them the same color. Each cluster is a subtree. The function makes
#' the `k` subtrees by cutting `k - 1` branches. The `method` argument sets how
#' the function chooses the branches to cut:
#' \itemize{
#'   \item `split_largest_subtrees` (default): The function splits one subtree
#'     at each step. It picks the subtree with the largest genetic span, which
#'     is the largest patristic distance between two of its tips. If two subtrees
#'     tie, it picks the subtree with more tips. In that subtree, a candidate
#'     branch must leave two subtrees that each have at least one tip. The
#'     function cuts the branch that makes the smaller of the two subtrees as
#'     large as possible (by tip count). If branches tie, it cuts the longest
#'     branch. If `cut_terminal_edges` is `FALSE`, the function considers
#'     internal branches first. It considers terminal branches (the branches
#'     that end in a tip) only if no internal branch is valid. This stops the
#'     function from cutting off long terminal branches before it splits the
#'     larger clades. The function stops after `k - 1` cuts, or when no valid
#'     cut remains.
#'   \item `cut_largest_branches`: The function cuts the `k - 1` longest
#'     branches of the whole tree. If `cut_terminal_edges` is `FALSE`, it cuts
#'     the longest internal branches first. It cuts terminal branches only if
#'     the tree has fewer than `k - 1` internal branches. This is the heuristic
#'     for fixed-k minimum spanning tree clustering, applied to a phylogenetic
#'     tree.
#'   \item `monophyletic_clades`: Each cluster is one complete clade, which is
#'     a node and all of its descendants. No cluster is a paraphyletic
#'     remainder. Use this method if each cluster must match a real clade, for
#'     example if each tip color in a figure must mark one monophyletic group.
#'     The function starts with the whole tree as one clade. At each step it
#'     picks the clade with the largest genetic span (the clade needs at least
#'     two tips). It replaces that clade with all of its child clades, which are
#'     the subtrees that start at the children of its root node. The function
#'     stops when `k` clades exist, or when no clade can split.
#'     `cut_terminal_edges` has no effect. For a bifurcating tree, the result
#'     has exactly `k` clades. At a multifurcation, one step adds all child
#'     clades, so the result can have more than `k` clades. If an ignored tip
#'     (such as the germline) forms its own clade, the function removes that
#'     clade. The result can then have fewer than `k` clusters.
#' }
#' The function labels the clusters with `prefix` and a number, such as
#' `tree_cluster_1`.
#'
#' @param tree Phylogenetic tree with the fields `edge`, `edge.length`, and
#'   `tip.label`.
#' @param k Number of subtrees to make, before the function removes the ignored
#'   tips. If `k` is larger than the number of tips, the function gives a
#'   warning and uses the number of tips.
#' @param id Name of the tip ID column in the result (default `"sequence_id"`).
#' @param cluster Name of the cluster column in the result (default
#'   `"tree_cluster"`).
#' @param ignore_tips Tip labels to remove **before** the clustering (default
#'   `"Germline"`). These tips do not form a cluster, do not use up a cluster,
#'   and do not appear in the result. The function always removes a tip named
#'   `"Germline"` (in any case). It prunes the tree with `ape::drop.tip()`. If
#'   `ape` is not installed, the function removes the tips from the result table
#'   after the clustering.
#' @param prefix Prefix of the cluster labels (default `"tree_cluster_"`).
#' @param method Clustering method: `"split_largest_subtrees"` (default),
#'   `"cut_largest_branches"`, or `"monophyletic_clades"`. See Details.
#' @param cut_terminal_edges Allow cuts of terminal branches (default `TRUE`)? If
#'   `FALSE`, the function cuts internal branches first, so that it splits larger
#'   clades before it cuts off single tips.
#'
#' @return A data frame with one row for each kept tip. It has the tip ID column
#'   (named by `id`) and the cluster label column (named by `cluster`).
#' @seealso [chooseTreeTipK()], [clusterTreesTips()]
#'
#' @examples
#' x <- readSIRUP(system.file("extdata", "sirup_synthetic.rds", package = "SIRUP"))
#' x <- formatSpatialClones(
#'   x,
#'   id = "sequence_id",
#'   germ = "germline_alignment_d_mask",
#'   traits = "spatial_id"
#' )
#' x <- getSpatialTrees(x)
#' trees <- getAnalysis(x, "trees")
#' tree <- trees$trees[[1]]
#'
#' # Cut the tree into 3 clusters.
#' clusterTreeTips(tree, k = 3)
#'
#' # Make each cluster one complete clade.
#' clusterTreeTips(tree, k = 3, method = "monophyletic_clades")
#' @export
clusterTreeTips <- function(tree,
                            k,
                            id = "sequence_id",
                            cluster = "tree_cluster",
                            ignore_tips = "Germline",
                            prefix = "tree_cluster_",
                            method = c("split_largest_subtrees", "cut_largest_branches", "monophyletic_clades"),
                            cut_terminal_edges = TRUE) {
  method <- match.arg(method)
  if (!is.list(tree) || !all(c("edge", "edge.length", "tip.label") %in% names(tree))) {
    rlang::abort("`tree` must have `edge`, `edge.length`, and `tip.label` fields.")
  }
  if (!is.numeric(k) || length(k) != 1L || is.na(k) || k < 1) {
    rlang::abort("`k` must be a positive integer.")
  }
  k <- as.integer(k)

  # Exclude the germline (and any explicitly ignored tips) BEFORE clustering, so
  # it never forms or consumes a cluster. It is absent from the returned table
  # and therefore stays uncoloured (grey) downstream. Matches a tip named
  # "Germline"/"germline" (any case) as well as anything in `ignore_tips`.
  drop_now <- union(
    if (!is.null(ignore_tips)) intersect(ignore_tips, tree$tip.label) else character(0),
    tree$tip.label[grepl("^germline$", tree$tip.label, ignore.case = TRUE)]
  )
  if (length(drop_now) > 0L &&
      (length(tree$tip.label) - length(drop_now)) >= 1L &&
      requireNamespace("ape", quietly = TRUE)) {
    # `ape::drop.tip()` is an S3 generic for `phylo`, but this function accepts a
    # bare edge/tip list (the documented input). Coerce to a minimal valid
    # `phylo` (supplying `Nnode` if missing) before pruning; the rest of the
    # function stays class-agnostic.
    if (!inherits(tree, "phylo")) {
      if (is.null(tree$Nnode)) {
        tree$Nnode <- length(unique(as.integer(tree$edge))) - length(tree$tip.label)
      }
      class(tree) <- "phylo"
    }
    tree <- ape::drop.tip(tree, drop_now)
  }

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
    cut_edges <- .sirup_tree_cut_edges(
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

#' Choose the number of tip clusters for a tree
#'
#' `chooseTreeTipK` picks a number of tip clusters for one tree.
#'
#' @details
#' For each value in `k_range`, the function runs [clusterTreeTips()]. It scores
#' the result as the mean within cluster patristic distance. It finds the
#' first `k` where the relative improvement of the score is below
#' `improvement_threshold`, and picks the `k` before it. If the improvement
#' never falls below the threshold, it picks the largest `k`.
#'
#' @param tree Phylogenetic tree with the fields `edge`, `edge.length`, and
#'   `tip.label`.
#' @param k_range Values of `k` to test (default `1:10`).
#' @param improvement_threshold Smallest relative improvement of the score that
#'   justifies a larger `k` (default 0.10).
#' @param ignore_tips Tip labels to remove before the scoring (default
#'   `"Germline"`).
#' @param method Clustering method for [clusterTreeTips()]. The default is
#'   `"cut_largest_branches"`.
#' @param cut_terminal_edges Allow cuts of terminal branches (default `TRUE`)?
#'   The function passes the value to [clusterTreeTips()].
#'
#' @return A list with the chosen `k` and `scores`, a data frame with the columns
#'   `k`, `score`, and `relative_improvement`.
#' @seealso [clusterTreeTips()], [clusterTreesTips()]
#'
#' @examples
#' x <- readSIRUP(system.file("extdata", "sirup_synthetic.rds", package = "SIRUP"))
#' x <- formatSpatialClones(
#'   x,
#'   id = "sequence_id",
#'   germ = "germline_alignment_d_mask",
#'   traits = "spatial_id"
#' )
#' x <- getSpatialTrees(x)
#' trees <- getAnalysis(x, "trees")
#' tree <- trees$trees[[1]]
#'
#' chooseTreeTipK(tree)
#' @export
chooseTreeTipK <- function(tree,
                           k_range = 1:10,
                           improvement_threshold = 0.10,
                           ignore_tips = "Germline",
                           method = c("cut_largest_branches", "split_largest_subtrees", "monophyletic_clades"),
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
    cluster_score <- .sirup_tree_cluster_score(
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

#' Cluster the tips of all trees in a dowser result
#'
#' `clusterTreesTips` runs [clusterTreeTips()] on each tree of a dowser result.
#' It adds the cluster column to the clone data of each tree.
#'
#' @param trees Result of [dowser::getTrees()], or a list with an element
#'   `trees` that holds this result.
#' @param k Number of subtrees to make in each tree. Use `"auto"`
#'   (default) to choose `k` for each tree with [chooseTreeTipK()].
#' @param id Name of the column in each clone data frame that matches the tree
#'   tip labels (default `"sequence_id"`).
#' @param cluster Name of the cluster column to add (default `"tree_cluster"`).
#' @param prefix Prefix of the cluster labels (default `"tree_cluster_"`).
#' @param ignore_tips Tip labels to remove before the clustering (default
#'   `"Germline"`).
#' @param keep Name of a column in each clone data frame (default `NULL`). The
#'   function runs the clustering first. Then, for each tip where `keep` differs
#'   from `keep_val`, it replaces the cluster label with the value of `keep`.
#' @param keep_val Value of `keep` that marks a tip to keep.
#' @param method Clustering method for [clusterTreeTips()]. The default is
#'   `"split_largest_subtrees"`.
#' @param cut_terminal_edges Allow cuts of terminal branches (default `FALSE`)?
#'   With `FALSE`, the function prefers internal splits of clades. With `TRUE`,
#'   branch cuts can remove single terminal tips.
#' @param k_range Values of `k` to test, if `k = "auto"` (default `1:8`). The
#'   function passes it to [chooseTreeTipK()].
#' @param improvement_threshold Smallest relative improvement that justifies a
#'   larger `k`, if `k = "auto"` (default 0.10). The function passes it to
#'   [chooseTreeTipK()].
#'
#' @return `trees`, with the `cluster` column added to the clone data of each
#'   tree.
#' @seealso [clusterTreeTips()], [chooseTreeTipK()], [plotSpatialCloneTreePairs()]
#'
#' @examples
#' x <- readSIRUP(system.file("extdata", "sirup_synthetic.rds", package = "SIRUP"))
#' x <- formatSpatialClones(
#'   x,
#'   id = "sequence_id",
#'   germ = "germline_alignment_d_mask",
#'   traits = "spatial_id"
#' )
#' x <- getSpatialTrees(x)
#' trees <- getAnalysis(x, "trees")
#' tree <- trees$trees[[1]]
#'
#' # Use 3 clusters for every tree.
#' clustered <- clusterTreesTips(trees, k = 3)
#' table(clustered$data[[1]]@data$tree_cluster)
#'
#' # Choose k for each tree.
#' clustered <- clusterTreesTips(trees, k = "auto")
#' @export
clusterTreesTips <- function(trees,
                             k = "auto",
                             id = "sequence_id",
                             cluster = "tree_cluster",
                             prefix = "tree_cluster_",
                             ignore_tips = "Germline",
                             keep = NULL,
                             keep_val = NULL,
                             method = c("split_largest_subtrees", "cut_largest_branches", "monophyletic_clades"),
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
  use_keep_filter <- .sirup_check_tree_keep_args(keep, keep_val)
  for (i in seq_len(nrow(out))) {
    tree_for_clustering <- out$trees[[i]]
    extra_tip_clusters <- NULL
    if (use_keep_filter) {
      extra_tip_clusters <- .sirup_tree_non_keep_clusters(
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
          .sirup_tree_row_label(out, i),
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
    tree_tip_clusters[[cluster]] <- .sirup_compact_tree_cluster_labels(
      tree_tip_clusters[[cluster]],
      prefix = prefix
    )
    out$data[[i]] <- .sirup_add_tree_tip_clusters(
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

#' Format the receptors of a SIRUP object as in dowser's `formatClones()`
#'
#' `formatSpatialClones` runs [dowser::formatClones()] on the receptor table of a
#' SIRUP object. It stores the result as an analysis.
#'
#' @details
#' The function puts clones into a standardized format for downstream analysis.
#' See [dowser::formatClones()] for details. The function stores the result as an
#' analysis in the SIRUP object. The default analysis name is `"f_clones"`.
#'
#' [dowser::formatClones()] takes a data frame and returns a tibble of
#' `airrClone` objects. `formatSpatialClones` calls it with the receptor table of
#' a SIRUP object and adds four things that spatial analysis needs:
#
#' This function automatically does not collapse identical sequences across distinct samples.
#' To plot the tips on the tissue, add the spatial ID column to `traits`.
#'
#'
#' @param object A SIRUP object.
#' @param analysis_name Name for the stored formatted clones (default
#'   `"f_clones"`). If the object has several receptor slots, the default also
#'   ends with `_` and the slot name.
#' @param slot Name of the receptor slot to use. Leave it `NULL` if the SIRUP
#'   object has one receptor table. You must give it if the object has several
#'   slots. 

#' @param seq Name of the column that holds the observed DNA sequences (default
#'   `"sequence_alignment"`). The sequences must be aligned to each other.
#' @param clone Name of the column that holds the clone ID (default
#'   `"clone_id"`).
#' @param subgroup Name of the column that holds the subgroup ID (default
#'   `"clone_subgroup"`).
#' @param id Name of the column that holds the sequence IDs (default
#'   `"sequence_id"`).
#' @param germ Name of the column that holds the germline DNA sequences (default
#'   `"germline_alignment_d_mask"`). Each clone must have one germline. It must
#'   be aligned to the sequences in `seq`.
#' @param v_call Name of the column that holds the V-segment alleles (default
#'   `"v_call"`). Within a clone, the entries must be the same at the gene level.
#' @param j_call Name of the column that holds the J-segment alleles (default
#'   `"j_call"`). Within a clone, the entries must be the same at the gene level.
#' @param junc_len Name of the numeric column that holds the junction length
#'   (default `"junction_length"`). Within a clone, the entries must be the same.
#' @param mask_char Character for masking and padding (default `"N"`).
#' @param max_mask Largest number of characters to mask at the start and at the
#'   end of the sequences (default 0). With `0`, the function does not mask. With
#'   `NULL`, the function uses the largest number of leading or trailing Ns in
#'   the sequences.
#' @param pad_end Pad the end of each sequence with `mask_char` so that all
#'   sequences have the same length (default `TRUE`)?
#' @param text_fields Names of the text columns to keep and merge when the
#'   function removes duplicates.
#' @param num_fields Names of the numeric columns to keep and sum when the
#'   function removes duplicates.
#' @param seq_fields Names of the sequence columns to keep and collapse when the
#'   function removes duplicates. Do not repeat `seq` or `germ` here.
#' @param add_count Add the column `COLLAPSE_COUNT`, which holds the number of
#'   sequences that the function collapsed (default `TRUE`)?
#' @param verbose Report the numbers of input, dropped, and output sequences
#'   (default `FALSE`)? The function passes it to `collapseDuplicates`.
#' @param collapse Collapse identical sequences (default `TRUE`)?
#' @param cell Name of the column that holds the cell IDs (default `"cell_id"`).
#' @param locus Name of the column that holds the locus (default `"locus"`).
#' @param traits Names of the columns to keep separate when the function
#'   collapses sequences. Include the spatial ID column if you want to plot the
#'   tips on the tissue.
#' @param mod3 Pad the sequences to a length that is a multiple of three
#'   (default `TRUE`)?
#' @param randomize Put the sequences in random order (default `TRUE`)? This is
#'   important for PHYLIP.
#' @param use_regions Assign the CDR and FWR regions (default `TRUE`)?
#' @param dup_singles Duplicate the sequence of a clone that has one sequence, so
#'   that the clone gives a tree (default `FALSE`)?
#' @param nproc Number of cores for the formatting (default 1).
#' @param chain Use `"HL"` to include the light chain, if it exists. The default
#'   is `"H"`.
#' @param heavy Name of the heavy chain locus (default `"IGH"`).
#' @param filterstop Use only sequences with no in-frame stop codon (default
#'   `TRUE`)?
#' @param minseq Smallest number of sequences in a clone (default 2).
#' @param split_light Split the subgroups, or lump them together (default
#'   `FALSE`)? See `resolveLightChains` in dowser.
#' @param light_traits Include the traits of the light chain when the function
#'   joins and collapses the sequences (default `FALSE`)?
#' @param majoronly Return only the largest subgroup and the sequences with no
#'   light chain (default `FALSE`)?
#' @param columns Names of other columns to keep in the result.
#' @param validate Validate the updated SIRUP object (default `TRUE`)?
#' @param slot Name of the receptor slot to use. Leave it `NULL` if the SIRUP
#'   object has one receptor table. You must give it if the object has several
#'   slots.
#'
#' @return The updated SIRUP object, with the formatted clones as an analysis.
#' @seealso [getSpatialTrees()], [dowser::formatClones()]
#'
#' @examples
#' x <- readSIRUP(system.file("extdata", "sirup_synthetic.rds", package = "SIRUP"))
#' x <- formatSpatialClones(
#'   x,
#'   id = "sequence_id",
#'   germ = "germline_alignment_d_mask",
#'   traits = "spatial_id"
#' )
#' getAnalysis(x, "f_clones")
#' @export
formatSpatialClones <- function(object,
                                analysis_name = NULL,
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
                                slot = NULL,
                                validate = TRUE) {
  validateSIRUP(object)
  receptors <- getReceptors(object, slot = slot)
  if (!germ %in% colnames(receptors)) {
    rlang::abort(paste0(
      "The receptor table does not contain `", germ,
      "`. Run `updateReceptors()` on this SIRUP object with ",
      "`dowser::createGermlines()` before `formatSpatialClones()`."
    ))
  }
  # Carry a verbatim copy of the receptor id (and the sample column, when the
  # SIRUP `sample` field is set) through tree building. dowser sanitizes tip
  # labels for Newick output (":" and "," become "_"), so the tree `id` column
  # can no longer be matched against the raw receptor id. `.sirup_uid` is passed
  # as a text field, which dowser preserves verbatim (and comma-joins when it
  # collapses identical sequences), so `plotSpatialCloneTreePairs()` can match
  # tree tips back to this object's receptors. See that function for the match.
  uid_col <- ".sirup_uid"
  receptors[[uid_col]] <- as.character(receptors[[id]])
  text_fields <- union(text_fields, uid_col)
  sample_field <- .sirup_receptor_fields(object)$sample
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
  analysis_name <- analysis_name %||% .sirup_slot_analysis_name(object, "f_clones", slot)
  setAnalysis(object, analysis_name, f_clones, validate = validate)
}

#' Build trees from the formatted clones of a SIRUP object
#'
#' `getSpatialTrees` runs [dowser::getTrees()] on the stored formatted clones. It
#' stores the trees as an analysis.
#'
#' @details
#' [dowser::getTrees()] takes a tibble of formatted clones and returns a tibble
#' with a `trees` column. `getSpatialTrees` does not change how dowser builds the
#' trees. It reads the formatted clones from the SIRUP object by name
#' (`f_clones_analysis`) and stores the trees in the object as an analysis. The
#' trees and the clone data of each tree then stay with the spatial data.
#' [plotSpatialCloneTreePairs()] and [clusterTreesTips()] can use them without
#' more input.
#'
#' By default, the function builds maximum parsimony trees with
#' `phangorn::pratchet()`. It estimates the branch lengths with `ape::acctran()`.
#' If you use `build = "igphyml"`, the function predicts the values of the
#' internal nodes with maximum parsimony. In this case, give `dir` (a directory
#' for the temporary files) and `id` (a unique name for this analysis). Use a
#' different `id` for each call that shares a `dir`, so that the calls do not
#' overwrite each other's files. The function writes `modelfile` if you do not
#' give it. The file has no constraints. The function deletes the temporary files
#' unless you set `rm_temp = FALSE`.
#'
#' The arguments from `trait` to `...` go to [dowser::getTrees()].
#'
#' @param object A SIRUP object.
#' @param f_clones_analysis Name of the stored formatted clones. If `NULL`, the
#'   function uses `"f_clones"`, or `"f_clones_<slot>"` if the object has several
#'   receptor slots.
#' @param analysis_name Name for the stored trees (default `"trees"`). If the
#'   object has several receptor slots, the default also ends with `_` and the
#'   slot name.
#' @param slot Name of the receptor slot that the trees come from. The function
#'   uses it only to find the default analysis names. Leave it `NULL` if the
#'   object has one receptor table. You must give it if the object has several
#'   slots and you use the default names.
#' @param trait Trait for the parsimony models. You need it if you use
#'   `igphyml`.
#' @param id Unique name for this analysis. You need it if you use `igphyml` or
#'   `dnapars`.
#' @param dir Directory for the temporary files.
#' @param modelfile File that sets the parsimony model.
#' @param build Program that builds the trees: `pratchet` (default), `pml`,
#'   `dnapars`, `dnaml`, `igphyml`, or `raxml`.
#' @param exec Location of the phylogenetic program.
#' @param igphyml Location of the `igphyml` program, for parsimony.
#' @param fixtrees Use the tree topologies that you supply (default `FALSE`)?
#' @param nproc Number of cores for the calculations (default 1).
#' @param quiet Amount of text that the function prints to the console (default
#'   0).
#' @param rm_temp Delete the temporary files (default `TRUE`)?
#' @param palette Deprecated.
#' @param seq Name of the column that holds the sequences.
#' @param collapse Collapse the internal nodes that have identical sequences
#'   (default `FALSE`)?
#' @param check_divergence Run `checkDivergence` on the trees (default `TRUE`)?
#' @param ... Other arguments for the tree-building programs.
#' @param validate Validate the updated SIRUP object (default `TRUE`)?
#'
#' @return The updated SIRUP object, with the trees as an analysis.
#' @seealso [formatSpatialClones()], [clusterTreesTips()],
#'   [plotSpatialCloneTreePairs()]
#'
#' @examples
#' x <- readSIRUP(system.file("extdata", "sirup_synthetic.rds", package = "SIRUP"))
#' x <- formatSpatialClones(
#'   x,
#'   id = "sequence_id",
#'   germ = "germline_alignment_d_mask",
#'   traits = "spatial_id"
#' )
#' x <- getSpatialTrees(x)
#' getAnalysis(x, "trees")
#' @export
getSpatialTrees <- function(object,
                            f_clones_analysis = NULL,
                            analysis_name = NULL,
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
                            slot = NULL,
                            validate = TRUE) {
  validateSIRUP(object)
  f_clones_analysis <- f_clones_analysis %||% .sirup_slot_analysis_name(object, "f_clones", slot)
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
  analysis_name <- analysis_name %||% .sirup_slot_analysis_name(object, "trees", slot)
  setAnalysis(object, analysis_name, trees, validate = validate)
}

#' Assign colors to the tips of clones
#'
#' `assignTipColors` gives each tip of each clone a color ID and a color. Tree
#' plots and the matching spatial plots use the same colors.
#'
#' @details
#' The function numbers the tips of each clone in the order of `id`. If you give
#' `shared_color`, tips with the same value of that column get the same number
#' within a clone. For example, tips in the same spot get the same color.
#'
#' @param clones Data frame of AIRR data with clone IDs.
#' @param clone_val Clone IDs to keep before the function assigns the colors. If
#'   `NULL`, the function keeps all clones.
#' @param clone Name of the column that holds the clone IDs (default
#'   `"clone_id"`).
#' @param id Name of the column that holds the sequence or tip IDs (default
#'   `"sequence_id"`).
#' @param shared_color Name of a column. Tips that share a value in this column
#'   get the same color within a clone, for example a spatial ID column.
#' @param palette Vector of colors. If `NULL`, the function uses
#'   `pals::glasbey()` if the package is installed. Otherwise, it uses the base R
#'   colors from `grDevices::hcl.colors()`.
#'
#' @return A data frame with the columns `tip_color_id` and `tip_color` added.
#' @seealso [plotSpatialCloneTreePairs()]
#'
#' @examples
#' x <- readSIRUP(system.file("extdata", "sirup_synthetic.rds", package = "SIRUP"))
#' colors <- assignTipColors(getReceptors(x), shared_color = "spatial_id")
#' head(colors[, c("clone_id", "spatial_id", "tip_color_id", "tip_color")])
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
      palette <- .sirup_extend_palette(pals::glasbey(32), n_colors)
    } else {
      palette <- grDevices::hcl.colors(max(12, n_colors), palette = "Dark 3")
    }
  } else if (length(palette) < n_colors) {
    palette <- .sirup_extend_palette(palette, n_colors)
  }
  out$tip_color <- palette[as.integer(out$tip_color_id)]
  out
}

# Dispatch to the requested strategy for choosing which tree edges to remove.
.sirup_tree_cut_edges <- function(edge,
                                  edge_length,
                                  n_tips,
                                  n_cut,
                                  cut_terminal_edges,
                                  method) {
  if (identical(method, "cut_largest_branches")) {
    return(.sirup_tree_cut_largest_branches(
      edge = edge,
      edge_length = edge_length,
      n_tips = n_tips,
      n_cut = n_cut,
      cut_terminal_edges = cut_terminal_edges
    ))
  }

  if (identical(method, "monophyletic_clades")) {
    return(.sirup_tree_split_clades(
      edge = edge,
      edge_length = edge_length,
      n_tips = n_tips,
      n_cut = n_cut
    ))
  }

  .sirup_tree_split_largest_subtrees(
    edge = edge,
    edge_length = edge_length,
    n_tips = n_tips,
    n_cut = n_cut,
    cut_terminal_edges = cut_terminal_edges
  )
}

# Choose the globally longest eligible branches for fixed-k tree clustering.
.sirup_tree_cut_largest_branches <- function(edge,
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
.sirup_tree_split_largest_subtrees <- function(edge,
                                               edge_length,
                                               n_tips,
                                               n_cut,
                                               cut_terminal_edges) {
  keep_edge <- rep(TRUE, nrow(edge))
  cut_edges <- integer()
  nodes <- sort(unique(as.integer(edge)))

  for (step in seq_len(n_cut)) {
    components <- .sirup_tree_components(edge, keep_edge, nodes)
    tip_counts <- vapply(
      components,
      function(component) sum(component <= n_tips),
      integer(1)
    )
    weighted_adjacency <- .sirup_tree_weighted_adjacency(edge, edge_length, keep_edge, nodes)
    component_scores <- vapply(
      components,
      .sirup_tree_component_diameter,
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
      best_edge <- .sirup_best_tree_cut_edge(
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
        best_edge <- .sirup_best_tree_cut_edge(
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
# Partition tips into strictly monophyletic clades. Unlike the branch-cutting
# methods, this never leaves a paraphyletic remainder: it starts from the root
# clade and repeatedly replaces the clade with the largest genetic span by ALL of
# its child clades, so every returned group is the complete descendant set of a
# single node. The cut set is the parent edges of the final clade roots; the
# tip-free backbone above the clades drops out during component labelling. For a
# bifurcating tree this cuts exactly `n_cut` edges (k = n_cut + 1 clades); at a
# multifurcation a split adds all children at once, so the result can have more
# than `n_cut + 1` clades.
.sirup_tree_split_clades <- function(edge, edge_length, n_tips, n_cut) {
  k_target <- n_cut + 1L
  nodes <- sort(unique(as.integer(edge)))

  children <- stats::setNames(vector("list", length(nodes)), as.character(nodes))
  parent_edge <- stats::setNames(rep(NA_integer_, length(nodes)), as.character(nodes))
  for (i in seq_len(nrow(edge))) {
    parent <- as.character(edge[i, 1])
    child <- as.character(edge[i, 2])
    children[[parent]] <- c(children[[parent]], as.integer(edge[i, 2]))
    parent_edge[[child]] <- i
  }
  root <- setdiff(as.integer(edge[, 1]), as.integer(edge[, 2]))[[1]]

  # Descendant node set (self + all nodes below) for every node, in post-order so
  # each node is filled after its children.
  descendants <- stats::setNames(vector("list", length(nodes)), as.character(nodes))
  for (nd in .sirup_tree_postorder(root, children)) {
    kids <- children[[as.character(nd)]]
    if (is.null(kids)) {
      descendants[[as.character(nd)]] <- nd
    } else {
      desc <- nd
      for (kid in kids) desc <- c(desc, descendants[[as.character(kid)]])
      descendants[[as.character(nd)]] <- desc
    }
  }
  n_desc_tips <- function(nd) sum(descendants[[as.character(nd)]] <= n_tips)

  weighted_adjacency <- .sirup_tree_weighted_adjacency(
    edge, edge_length, rep(TRUE, nrow(edge)), nodes
  )
  clade_diameter <- function(nd) {
    .sirup_tree_component_diameter(descendants[[as.character(nd)]], n_tips, weighted_adjacency)
  }

  frontier <- root
  while (length(frontier) < k_target) {
    splittable <- frontier[vapply(frontier, function(nd) {
      length(children[[as.character(nd)]]) > 0L && n_desc_tips(nd) >= 2L
    }, logical(1))]
    if (length(splittable) == 0L) {
      break
    }
    diams <- vapply(splittable, clade_diameter, numeric(1))
    tip_counts <- vapply(splittable, n_desc_tips, integer(1))
    best <- splittable[[order(diams, tip_counts, decreasing = TRUE)[[1]]]]
    frontier <- c(setdiff(frontier, best), children[[as.character(best)]])
  }

  cut_nodes <- frontier[frontier != root]
  cut_edges <- unname(parent_edge[as.character(cut_nodes)])
  cut_edges[!is.na(cut_edges)]
}

# Iterative post-order node ordering (descendants before ancestors) for a rooted
# tree given a parent -> children list.
.sirup_tree_postorder <- function(root, children) {
  out <- integer()
  stack <- root
  while (length(stack) > 0L) {
    nd <- stack[[length(stack)]]
    stack <- stack[-length(stack)]
    out <- c(nd, out)
    kids <- children[[as.character(nd)]]
    if (!is.null(kids)) {
      stack <- c(stack, kids)
    }
  }
  out
}

.sirup_best_tree_cut_edge <- function(candidate_edges,
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
    components <- .sirup_tree_components(edge, keep_after_cut, nodes)
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
.sirup_tree_components <- function(edge, keep_edge, nodes) {
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
.sirup_tree_weighted_adjacency <- function(edge, edge_length, keep_edge, nodes) {
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
.sirup_tree_component_diameter <- function(component, n_tips, weighted_adjacency) {
  tips <- component[component <= n_tips]
  if (length(tips) < 2L) {
    return(0)
  }

  max_distance <- 0
  for (tip in tips) {
    distances <- .sirup_tree_distances_from_node(tip, weighted_adjacency)
    tip_distances <- distances[as.character(tips)]
    max_distance <- max(max_distance, tip_distances, na.rm = TRUE)
  }
  max_distance
}

# Compute shortest branch-length distances from one tree node.
.sirup_tree_distances_from_node <- function(start, weighted_adjacency) {
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
.sirup_tree_cluster_score <- function(tree,
                                      clusters,
                                      id_col,
                                      cluster_col) {
  dist_mat <- .sirup_tip_patristic_matrix(tree, clusters[[id_col]])
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

#' Calculate patristic distances between tips in a tree
#'
#' `tipPatristicMatrix` returns the matrix of patristic distances between the
#' selected tips of a tree. The rows and columns follow the order of
#' `tip_labels`.
#'
#' @param tree Phylogenetic tree with the fields `edge`, `edge.length`, and
#'   `tip.label`.
#' @param tip_labels Character vector of the tip labels to use.
#'
#' @return A numeric matrix of patristic distances.
#' @seealso [clusterTreeTips()]
#'
#' @examples
#' x <- readSIRUP(system.file("extdata", "sirup_synthetic.rds", package = "SIRUP"))
#' x <- formatSpatialClones(
#'   x,
#'   id = "sequence_id",
#'   germ = "germline_alignment_d_mask",
#'   traits = "spatial_id"
#' )
#' x <- getSpatialTrees(x)
#' trees <- getAnalysis(x, "trees")
#' tree <- trees$trees[[1]]
#'
#' tipPatristicMatrix(tree, tree$tip.label[1:4])
#' @export
tipPatristicMatrix <- function(tree, tip_labels) {
  .sirup_tip_patristic_matrix(tree = tree, tip_labels = tip_labels)
}

# Return patristic distances among selected tips, preserving tip order.
.sirup_tip_patristic_matrix <- function(tree, tip_labels) {
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
    .sirup_node_depth(as.integer(node), parent, branch_length, depth)
  }
  for (node in names(depth)) {
    depth[[node]] <- .sirup_node_depth(as.integer(node), parent, branch_length, depth)
  }

  tip_nodes <- stats::setNames(seq_len(n_tips), tree$tip.label)
  tip_nodes <- tip_nodes[tip_labels]
  out <- matrix(0, nrow = length(tip_labels), ncol = length(tip_labels), dimnames = list(tip_labels, tip_labels))
  ancestors <- lapply(tip_nodes, .sirup_node_ancestors, parent = parent)
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
.sirup_node_depth <- function(node, parent, branch_length, depth) {
  node_name <- as.character(node)
  if (!is.na(depth[[node_name]])) {
    return(depth[[node_name]])
  }
  parent_node <- parent[[node_name]]
  if (is.na(parent_node)) {
    return(0)
  }
  .sirup_node_depth(parent_node, parent, branch_length, depth) + branch_length[[node_name]]
}

# Walk from a node to the root and return its ancestor chain.
.sirup_node_ancestors <- function(node, parent) {
  ancestors <- integer()
  current <- node
  while (!is.na(current)) {
    ancestors <- c(ancestors, current)
    current <- parent[[as.character(current)]]
  }
  ancestors
}

# Relabel arbitrary component IDs as compact, stable cluster labels.
.sirup_compact_tree_cluster_labels <- function(labels, prefix) {
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
.sirup_tree_row_label <- function(trees, i) {
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
.sirup_check_tree_keep_args <- function(keep_col, keep_value) {
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
.sirup_tree_non_keep_clusters <- function(tree,
                                          clone_data,
                                          id_col,
                                          keep_col,
                                          keep_value,
                                          cluster_col,
                                          ignore_tips) {
  data <- .sirup_tree_clone_data_frame(clone_data)
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
.sirup_tree_clone_data_frame <- function(clone_data) {
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
.sirup_add_tree_tip_clusters <- function(clone_data,
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
.sirup_extend_palette <- function(palette, n) {
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
