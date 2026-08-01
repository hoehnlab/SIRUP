test_that("clusterTreeTips cuts the longest branches into connected tip clusters", {
  tree <- list(
    edge = matrix(
      c(
        7, 5,
        7, 6,
        5, 1,
        5, 2,
        6, 3,
        6, 4
      ),
      ncol = 2,
      byrow = TRUE
    ),
    edge.length = c(10, 1, 1, 1, 1, 1),
    tip.label = c("A", "B", "C", "D")
  )

  clusters <- clusterTreeTips(tree, k = 2, id = "tip", cluster = "cluster")
  cluster_lookup <- stats::setNames(clusters$cluster, clusters$tip)

  expect_equal(cluster_lookup[["A"]], cluster_lookup[["B"]])
  expect_equal(cluster_lookup[["C"]], cluster_lookup[["D"]])
  expect_false(identical(cluster_lookup[["A"]], cluster_lookup[["C"]]))
  expect_equal(length(unique(clusters$cluster)), 2)
})

test_that("clusterTreeTips can ignore germline tips", {
  tree <- list(
    edge = matrix(
      c(
        6, 5,
        5, 1,
        5, 2,
        6, 3,
        6, 4
      ),
      ncol = 2,
      byrow = TRUE
    ),
    edge.length = c(5, 1, 1, 1, 1),
    tip.label = c("A", "B", "C", "Germline")
  )

  clusters <- clusterTreeTips(tree, k = 2, id = "tip", cluster = "cluster")

  expect_false("Germline" %in% clusters$tip)
  expect_equal(nrow(clusters), 3)
})

test_that("chooseTreeTipK returns a selected k and score table", {
  tree <- list(
    edge = matrix(
      c(
        7, 5,
        7, 6,
        5, 1,
        5, 2,
        6, 3,
        6, 4
      ),
      ncol = 2,
      byrow = TRUE
    ),
    edge.length = c(10, 1, 1, 1, 1, 1),
    tip.label = c("A", "B", "C", "D")
  )

  choice <- chooseTreeTipK(tree, k_range = 1:4)

  expect_true(choice$k %in% 1:4)
  expect_true(all(c("k", "score", "relative_improvement") %in% colnames(choice$scores)))
  expect_false(any(c("criterion", "score") %in% names(formals(chooseTreeTipK))))
})

test_that("clusterTreesTips adds clusters to dowser-style tree data frames", {
  tree <- list(
    edge = matrix(
      c(
        7, 5,
        7, 6,
        5, 1,
        5, 2,
        6, 3,
        6, 4
      ),
      ncol = 2,
      byrow = TRUE
    ),
    edge.length = c(10, 1, 1, 1, 1, 1),
    tip.label = c("A", "B", "C", "D")
  )
  trees <- data.frame(clone_id = "clone_1", stringsAsFactors = FALSE)
  trees$trees <- list(tree)
  trees$data <- list(data.frame(sequence_id = c("A", "B", "C", "D")))

  out <- clusterTreesTips(trees, k = "auto", k_range = 1:3)

  expect_true("tree_cluster" %in% colnames(out$data[[1]]))
  expect_true("tree_cluster_k" %in% colnames(out))
  expect_false(any(c("criterion", "score") %in% names(formals(clusterTreesTips))))
})
