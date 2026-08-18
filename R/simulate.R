#' Simulate network data with spillovers
#'
#' @description
#' Simulates covariates, a network, treatment assignments, and outcomes under a
#' first-order network-interference model.
#'
#' @param parameters_graph A list with entries `n`, `type_graph`, and `neighb`.
#' @param parameters_model A list with entries `mu1`, `mu2`, `mu3`, `mu4`,
#'   `heterogeneity`, `p`, `variance`, and `type`.
#' @param D Optional binary treatment vector. If `NULL`, treatment is assigned at
#'   random with probability 1/2.
#' @param W Optional adjacency matrix. If `NULL`, a graph is simulated.
#' @param seed Optional random seed.
#'
#' @return A list containing `y`, `X`, `D`, `W`, `W_0`, `overall_effect`,
#'   `estimator`, and `expected_estimator`.
#' @export
simulate_network_data <- function(
    parameters_graph = list(n = 100, type_graph = "geometric", neighb = 2),
    parameters_model = list(
      mu1 = NA_real_, mu2 = NA_real_, mu3 = NA_real_, mu4 = NA_real_,
      heterogeneity = FALSE, p = 4, variance = 1, type = "linear"
    ),
    D = NULL,
    W = NULL,
    seed = NULL) {

  if (!is.null(seed)) {
    set.seed(seed)
  }

  n <- parameters_graph$n
  if (!is.null(D)) {
    n <- length(D)
  }

  p <- parameters_model$p
  if (!is.numeric(p) || length(p) != 1L || p < 1) {
    stop("`parameters_model$p` must be a positive integer.")
  }
  p <- as.integer(p)

  coeff <- draw_model_coefficients(p, isTRUE(parameters_model$heterogeneity))

  X <- matrix(stats::runif(n * p, min = -1, max = 1), nrow = n, ncol = p)

  if (is.null(W)) {
    graph_type <- parameters_graph$type_graph
    neighb <- parameters_graph$neighb

    W <- switch(
      graph_type,
      geometric = build_geometric_graph(X),
      `Erdos-Renyi` = build_erdos_renyi_graph(n, neighb),
      Barabasi = build_barabasi_graph(n, neighb),
      stop("Unsupported graph type. Use 'geometric', 'Erdos-Renyi', or 'Barabasi'.")
    )
  }

  W <- validate_adjacency_matrix(W, binary = TRUE)
  W_0 <- W

  if (is.null(D)) {
    D <- stats::rbinom(nrow(W), size = 1, prob = 0.5)
  }
  D <- as.integer(D)

  mu1 <- parameters_model$mu1
  mu2 <- parameters_model$mu2
  mu3 <- parameters_model$mu3
  mu4 <- parameters_model$mu4

  if (is.na(mu1)) {
    mu1 <- coeff$mu
    mu2 <- coeff$mu
    mu3 <- coeff$mu
    mu4 <- coeff$mu
  }

  spill1 <- as.vector(mu1 + X %*% coeff$b3)
  spill2 <- as.vector(mu2 + X %*% coeff$b5)
  cate <- as.vector(mu3 + X %*% coeff$b4)

  errors <- stats::rnorm(n, sd = sqrt(parameters_model$variance))
  max_deg <- max(1, max(rowSums(W_0)))
  errors <- errors / sqrt(2) + as.vector(W_0 %*% errors) / sqrt(2 * max_deg)

  neighbors <- pmax(rowSums(W_0), 1)
  D_all <- rep(1L, nrow(W_0))

  model_type <- parameters_model$type
  if (identical(model_type, "linear")) {
    y <- as.vector(W_0 %*% D) * spill2 / neighbors +
      D * cate +
      D * as.vector(W_0 %*% D) * spill1 / neighbors +
      errors

    overall_effect <- mean(
      as.vector(W_0 %*% D_all) * spill2 / neighbors +
        D_all * cate +
        D_all * as.vector(W_0 %*% D_all) * spill1 / neighbors
    )
  } else if (identical(model_type, "endogenous")) {
    y_helper <- as.vector(W_0 %*% D) * spill2 / neighbors +
      D * cate +
      D * as.vector(W_0 %*% D) * spill1 / neighbors +
      errors

    y_helper_all <- as.vector(W_0 %*% D_all) * spill2 / neighbors +
      D_all * cate +
      D_all * as.vector(W_0 %*% D_all) * spill1 / neighbors +
      errors

    system_matrix <- diag(nrow(W_0)) - mu4 * W_0
    y <- as.vector(solve(system_matrix, y_helper))
    y_all <- as.vector(solve(system_matrix, y_helper_all))
    overall_effect <- mean(y_all)
  } else {
    stop("Unsupported outcome model type. Use 'linear' or 'endogenous'.")
  }

  estimator <- safe_difference_in_means(y, D)
  expected_estimator <- safe_difference_in_means(y - errors, D)

  list(
    y = y,
    X = X,
    D = D,
    W = W_0,
    W_0 = W_0,
    overall_effect = overall_effect,
    estimator = estimator,
    expected_estimator = expected_estimator
  )
}

#' Assign treatment at the cluster level
#'
#' @param W Adjacency matrix.
#' @param cluster_design Logical; if `TRUE`, assign one binary treatment per
#'   cluster. If `FALSE`, draw a cluster-specific saturation probability.
#' @param clusters Integer vector of cluster labels.
#' @param seed Optional random seed.
#'
#' @return A list with entries `assignments` and `clusters`.
#' @export
assign_cluster_treatment <- function(W, cluster_design = TRUE, clusters, seed = NULL) {
  W <- validate_adjacency_matrix(W, binary = TRUE)
  if (length(clusters) != nrow(W)) {
    stop("`clusters` must have length equal to nrow(W).")
  }

  if (!is.null(seed)) {
    set.seed(seed)
  }

  assignments <- integer(nrow(W))
  cluster_ids <- unique(clusters)

  if (isTRUE(cluster_design)) {
    for (cl in cluster_ids) {
      assignments[clusters == cl] <- stats::rbinom(1, size = 1, prob = 0.5)
    }
  } else {
    probs <- stats::runif(length(cluster_ids))
    for (i in seq_along(cluster_ids)) {
      cl <- cluster_ids[i]
      assignments[clusters == cl] <- stats::rbinom(sum(clusters == cl), size = 1, prob = probs[i])
    }
  }

  list(assignments = assignments, clusters = clusters)
}
