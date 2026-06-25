#' Build a clone-by-spatial-barcode count matrix
#'
#' @param clones Clone-assigned AIRR data frame.
#' @param clone_col,barcode_col Column names.
#' @param sparse Return a sparse Matrix when the optional `Matrix` package exists.
#' @return Clone-by-barcode count matrix.
#' @export
buildCloneCountMatrix <- function(clones,
                                     clone_col = "clone_id",
                                     barcode_col = "barcode",
                                     sparse = TRUE) {
  .check_required_cols(clones, c(clone_col, barcode_col), "clones")
  x <- clones[!is.na(clones[[clone_col]]) & !is.na(clones[[barcode_col]]), , drop = FALSE]
  mat <- table(x[[clone_col]], x[[barcode_col]]) %>% as.matrix()
  if (sparse && requireNamespace("Matrix", quietly = TRUE)) {
    mat <- Matrix::Matrix(mat, sparse = TRUE)
  }
  mat
}

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
#' @param id_col Name of the returned tip identifier column.
#' @param cluster_col Name of the returned cluster column.
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
                            id_col = "sequence_id",
                            cluster_col = "tree_cluster",
                            ignore_tips = "Germline",
                            prefix = "tree_cluster_",
                            method = c("cut_largest_branches", "split_largest_subtrees"),
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
  out[[id_col]] <- out$.tip_id
  out[[cluster_col]] <- cluster_labels[as.character(out$.component)]
  out[, c(id_col, cluster_col), drop = FALSE]
}

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
#' @param criterion Selection rule. `"relative_improvement"` chooses the
#'   previous k once the score improvement falls below `improvement_threshold`.
#'   `"elbow"` chooses the point farthest from the line connecting the first and
#'   last candidate scores.
#' @param score Score to minimize. `"mean_pairwise"` uses the mean within-cluster
#'   patristic distance across all same-cluster tip pairs. `"sum_pairwise"` uses
#'   the total within-cluster patristic distance.
#' @param improvement_threshold Relative improvement threshold used when
#'   `criterion = "relative_improvement"`.
#' @param ignore_tips Optional tip labels to drop before scoring, for example
#'   `"Germline"`.
#' @param method,cut_terminal_edges Arguments passed to `clusterTreeTips()`.
#' @return List with selected `k` and a candidate score data frame.
#' @export
chooseTreeTipK <- function(tree,
                           k_range = 1:10,
                           criterion = c("relative_improvement", "elbow"),
                           score = c("mean_pairwise", "sum_pairwise"),
                           improvement_threshold = 0.10,
                           ignore_tips = "Germline",
                           method = c("cut_largest_branches", "split_largest_subtrees"),
                           cut_terminal_edges = TRUE) {
  criterion <- match.arg(criterion)
  score <- match.arg(score)
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
      id_col = "sequence_id",
      cluster_col = "tree_cluster",
      ignore_tips = ignore_tips,
      method = method,
      cut_terminal_edges = cut_terminal_edges
    )
    cluster_score <- .sirap_tree_cluster_score(
      tree,
      clusters,
      id_col = "sequence_id",
      cluster_col = "tree_cluster",
      score = score
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

  if (criterion == "relative_improvement") {
    low_gain <- which(!is.na(scores$relative_improvement) & scores$relative_improvement < improvement_threshold)
    chosen_k <- if (length(low_gain) > 0L) {
      scores$k[max(1L, low_gain[[1]] - 1L)]
    } else {
      scores$k[nrow(scores)]
    }
  } else {
    chosen_k <- .sirap_choose_elbow_k(scores)
  }

  list(k = as.integer(chosen_k), scores = scores)
}

.sirap_tree_cluster_score <- function(tree,
                                      clusters,
                                      id_col,
                                      cluster_col,
                                      score) {
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
  if (score == "sum_pairwise") {
    total
  } else if (n_pairs == 0L) {
    0
  } else {
    total / n_pairs
  }
}

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

.sirap_node_ancestors <- function(node, parent) {
  ancestors <- integer()
  current <- node
  while (!is.na(current)) {
    ancestors <- c(ancestors, current)
    current <- parent[[as.character(current)]]
  }
  ancestors
}

.sirap_choose_elbow_k <- function(scores) {
  if (nrow(scores) < 3L) {
    return(scores$k[nrow(scores)])
  }
  x <- scores$k
  y <- scores$score
  x <- (x - min(x)) / (max(x) - min(x))
  y <- if (max(y) == min(y)) rep(0, length(y)) else (y - min(y)) / (max(y) - min(y))
  start <- c(x[[1]], y[[1]])
  end <- c(x[[length(x)]], y[[length(y)]])
  line <- end - start
  distances <- vapply(seq_along(x), function(i) {
    point <- c(x[[i]], y[[i]])
    abs(line[[1]] * (start[[2]] - point[[2]]) - (start[[1]] - point[[1]]) * line[[2]]) / sqrt(sum(line^2))
  }, numeric(1))
  scores$k[which.max(distances)]
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
#' @param id_col Column in each clone data frame matching tree tip labels.
#' @param cluster_col Name of the cluster column to add.
#' @param ignore_tips Optional tip labels to drop before joining clusters.
#' @param prefix Prefix used for cluster labels.
#' @param keep_col Optional column in each clone data frame used to choose
#'   whether a tree has the expected retained tips. Clustering is run before
#'   any keep-filtering; tips where `keep_col != keep_value` then have their
#'   computed cluster labels overwritten with their `keep_col` value.
#' @param keep_value Expected value in `keep_col`.
#' @param method Tree clustering method passed to `clusterTreeTips()`.
#' @param cut_terminal_edges Cut terminal branches when choosing branches to
#'   remove. Defaults to `TRUE`; set to `FALSE` to prioritize internal clade
#'   splits.
#' @param k_range,criterion,score,improvement_threshold Arguments passed to
#'   `chooseTreeTipK()` when `k = "auto"`.
#' @return `trees`, with `cluster_col` added to each clone data object.
#' @export
clusterTreesTips <- function(trees,
                             k,
                             id_col = "sequence_id",
                             cluster_col = "tree_cluster",
                             ignore_tips = "Germline",
                             prefix = "tree_cluster_",
                             keep_col = NULL,
                             keep_value = NULL,
                             method = c("cut_largest_branches", "split_largest_subtrees"),
                             cut_terminal_edges = TRUE,
                             k_range = 1:8,
                             criterion = c("relative_improvement", "elbow"),
                             score = c("mean_pairwise", "sum_pairwise"),
                             improvement_threshold = 0.10) {
  method <- match.arg(method)
  auto_k <- identical(k, "auto")
  if (auto_k) {
    criterion <- match.arg(criterion)
    score <- match.arg(score)
  }
  if (is.list(trees) && !is.data.frame(trees) && "trees" %in% names(trees)) {
    trees$trees <- clusterTreesTips(
      trees$trees,
      k = k,
      id_col = id_col,
      cluster_col = cluster_col,
      ignore_tips = ignore_tips,
      prefix = prefix,
      keep_col = keep_col,
      keep_value = keep_value,
      method = method,
      cut_terminal_edges = cut_terminal_edges,
      k_range = k_range,
      criterion = criterion,
      score = score,
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
  use_keep_filter <- .sirap_check_tree_keep_args(keep_col, keep_value)
  for (i in seq_len(nrow(out))) {
    tree_for_clustering <- out$trees[[i]]
    extra_tip_clusters <- NULL
    if (use_keep_filter) {
      extra_tip_clusters <- .sirap_tree_non_keep_clusters(
        tree = out$trees[[i]],
        clone_data = out$data[[i]],
        id_col = id_col,
        keep_col = keep_col,
        keep_value = keep_value,
        cluster_col = cluster_col,
        ignore_tips = ignore_tips
      )
    }

    tree_k <- k
    if (auto_k) {
      k_choice <- chooseTreeTipK(
        tree_for_clustering,
        k_range = k_range,
        criterion = criterion,
        score = score,
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
      id_col = id_col,
      cluster_col = cluster_col,
      ignore_tips = ignore_tips,
      prefix = prefix,
      method = method,
      cut_terminal_edges = cut_terminal_edges
    )
    if (!is.null(extra_tip_clusters) && nrow(extra_tip_clusters) > 0L) {
      tree_tip_clusters <- tree_tip_clusters[
        !tree_tip_clusters[[id_col]] %in% extra_tip_clusters[[id_col]],
        ,
        drop = FALSE
      ]
      tree_tip_clusters <- rbind(tree_tip_clusters, extra_tip_clusters)
    }
    tree_tip_clusters[[cluster_col]] <- .sirap_compact_tree_cluster_labels(
      tree_tip_clusters[[cluster_col]],
      prefix = prefix
    )
    out$data[[i]] <- .sirap_add_tree_tip_clusters(
      out$data[[i]],
      tree_tip_clusters,
      id_col = id_col,
      cluster_col = cluster_col
    )
  }
  out[[paste0(cluster_col, "_k")]] <- chosen_k
  attr(out, paste0(cluster_col, "_scores")) <- k_scores
  out
}

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

#' Assign stable colors within selected clones
#'
#' Adds a per-tip color identifier and color value that can be used consistently
#' in tree plots and matched spatial plots.
#'
#' @param clones Clone-assigned AIRR data frame.
#' @param clone_ids Optional clone IDs to retain before assigning colors.
#' @param clone_col,id_col Column names for clone and sequence/tip identifiers.
#' @param shared_color_col Optional column whose shared values should receive
#'   the same tip color within each clone, for example a spatial barcode column.
#' @param palette Optional color vector. Defaults to `pals::glasbey()` when
#'   available, otherwise base R qualitative colors.
#' @return Data frame with `tip_color_id` and `tip_color` columns.
#' @export
assignTipColors <- function(clones,
                              clone_ids = NULL,
                              clone_col = "clone_id",
                              id_col = "sequence_id",
                              shared_color_col = NULL,
                              palette = NULL) {
  .check_required_cols(clones, c(clone_col, id_col), "clones")
  if (!is.null(shared_color_col)) {
    .check_required_cols(clones, shared_color_col, "clones")
  }
  out <- .as_data_frame(clones)
  if (!is.null(clone_ids)) {
    out <- out[out[[clone_col]] %in% clone_ids, , drop = FALSE]
  }
  out <- out[!is.na(out[[clone_col]]) & !is.na(out[[id_col]]), , drop = FALSE]
  if (is.null(shared_color_col)) {
    out <- out %>%
      dplyr::group_by(dplyr::across(dplyr::all_of(clone_col))) %>%
      dplyr::arrange(.data[[id_col]], .by_group = TRUE) %>%
      dplyr::mutate(tip_color_id = as.character(dplyr::row_number())) %>%
      dplyr::ungroup()
  } else {
    out <- out %>%
      dplyr::group_by(dplyr::across(dplyr::all_of(clone_col))) %>%
      dplyr::arrange(.data[[shared_color_col]], .data[[id_col]], .by_group = TRUE) %>%
      dplyr::mutate(tip_color_id = as.character(dplyr::dense_rank(.data[[shared_color_col]]))) %>%
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
