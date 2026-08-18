# Graph and clustering helpers ---------------------------------------------

#' Compute a graph Laplacian
#'
#' @param W Symmetric matrix used to construct the graph Laplacian.
#' @param normalized If `TRUE`, return the normalized Laplacian.
#' @param binary If `TRUE`, use the binary support of `W`.
#'
#' @return A square matrix.
#' @export
graph_laplacian <- function(W, normalized = TRUE, binary = TRUE) {
  W <- validate_adjacency_matrix(W, binary = binary)
  deg <- pmax(rowSums(W != 0), 1)
  n <- nrow(W)

  if (isTRUE(normalized)) {
    d_half <- diag(1 / sqrt(deg), nrow = n, ncol = n)
    diag(n) - d_half %*% W %*% d_half
  } else {
    diag(deg, nrow = n, ncol = n) - W
  }
}

#' Compute the left-normalized adjacency matrix
#'
#' @description Returns \eqn{L = V^{-1} A}, where \eqn{A} is the binary
#' adjacency matrix and \eqn{V_{ii} = max(1, |N_i|)}.
#'
#' @param W Symmetric adjacency matrix; non-zero entries are interpreted as edges.
#'
#' @return A square matrix.
#' @export
left_normalized_adjacency <- function(W) {
  W <- validate_adjacency_matrix(W, binary = TRUE)
  deg <- pmax(rowSums(W != 0), 1)
  diag(1 / deg, nrow = nrow(W), ncol = nrow(W)) %*% W
}

paper_left_normalized_adjacency <- left_normalized_adjacency

validate_xi <- function(xi) {
  if (!is.numeric(xi) || length(xi) != 1L || is.na(xi) || !is.finite(xi) || xi < 0) {
    stop("`xi` must be a single finite non-negative number.", call. = FALSE)
  }
  as.numeric(xi)
}

validate_positive_xi <- function(xi, name = "xi") {
  if (!is.numeric(xi) || length(xi) != 1L || is.na(xi) ||
      !is.finite(xi) || xi <= 0) {
    stop("`", name, "` must be a single finite positive number.", call. = FALSE)
  }
  as.numeric(xi)
}

validate_objective_type <- function(objective_type = "squared") {
  if (is.null(objective_type)) {
    return("squared")
  }
  if (!is.character(objective_type) || length(objective_type) != 1L ||
      is.na(objective_type) || !identical(objective_type, "squared")) {
    stop("`objective_type` must be \"squared\".", call. = FALSE)
  }
  "squared"
}

validate_gamma_bar <- function(gamma_bar) {

  if (!is.numeric(gamma_bar) || length(gamma_bar) != 1L ||
      is.na(gamma_bar) || !is.finite(gamma_bar) || gamma_bar < 1) {
    stop("`gamma_bar` must be a single finite number at least 1.", call. = FALSE)
  }
  as.numeric(gamma_bar)
}

validate_box_constraints <- function(box_constraints) {
  if (!is.logical(box_constraints) || length(box_constraints) != 1L || is.na(box_constraints)) {
    stop("`box_constraints` must be TRUE or FALSE.", call. = FALSE)
  }
  isTRUE(box_constraints)
}

validate_k_range <- function(min_k, max_k, n) {
  min_k <- max(1L, as.integer(min_k))
  max_k <- min(as.integer(max_k), n)
  if (min_k > max_k) {
    stop("`min_k` must be less than or equal to `max_k`.", call. = FALSE)
  }
  seq.int(min_k, max_k)
}

validate_nonnegative_grid <- function(x, name) {
  if (!is.numeric(x) || length(x) == 0L || anyNA(x) || any(!is.finite(x)) || any(x < 0)) {
    stop("`", name, "` must be a non-empty numeric vector of finite non-negative values.", call. = FALSE)
  }
  as.numeric(x)
}

validate_positive_grid <- function(x, name) {
  if (!is.numeric(x) || length(x) == 0L || anyNA(x) || any(!is.finite(x)) || any(x <= 0)) {
    stop("`", name, "` must be a non-empty numeric vector of finite positive values.", call. = FALSE)
  }
  as.numeric(x)
}

#' Available discretization methods
#'
#' @return Character vector of method names.
#' @export
available_discretization_methods <- function() {
  c(
    "kmeans",
    "hierarchical",
    "spectral_norm_kmeans",
    "spectral_unnorm_kmeans",
    "spectral_norm_hierarchical",
    "spectral_unnorm_hierarchical"
  )
}

standardize_discretization_methods <- function(methods) {
  if (is.null(methods) || length(methods) == 0L) {
    stop("`methods` must contain at least one discretization method.", call. = FALSE)
  }

  aliases <- c(
    search_kmeans = "kmeans",
    search_hierarchical = "hierarchical",
    spectral_normalized_kmeans = "spectral_norm_kmeans",
    spectral_unnormalized_kmeans = "spectral_unnorm_kmeans",
    spectral_normalized_hierarchical = "spectral_norm_hierarchical",
    spectral_unnormalized_hierarchical = "spectral_unnorm_hierarchical"
  )

  methods <- as.character(methods)
  mapped <- unname(ifelse(methods %in% names(aliases), aliases[methods], methods))
  choices <- available_discretization_methods()
  bad <- setdiff(mapped, choices)

  if (length(bad) > 0L) {
    stop(
      "Unknown discretization method(s): ", paste(bad, collapse = ", "),
      ". Valid methods are: ", paste(choices, collapse = ", "), ".",
      call. = FALSE
    )
  }

  unique(mapped)
}

spectral_embedding <- function(W, n_eig = 2, normalized = TRUE, binary = TRUE) {
  W <- validate_adjacency_matrix(W, binary = binary)
  n <- nrow(W)
  n_eig <- max(1L, min(as.integer(n_eig), n))

  L <- graph_laplacian(W, normalized = normalized, binary = FALSE)
  eig <- eigen(L, symmetric = TRUE)
  idx <- seq.int(from = n - n_eig + 1L, to = n)
  eig$vectors[, idx, drop = FALSE]
}

#' Epsilon-net clustering
#'
#' @param W Symmetric binary adjacency matrix.
#' @param epsilon Radius parameter. Only `epsilon = 3` is supported.
#'
#' @return A list with entry `clusters`.
#' @export
cluster_epsilon_net <- function(W, epsilon = 3) {
  W <- validate_adjacency_matrix(W, binary = TRUE)
  if (!identical(epsilon, 3)) {
    stop("Only `epsilon = 3` is supported in this implementation.", call. = FALSE)
  }

  n <- nrow(W)
  excluded <- rep(FALSE, n)
  net <- rep(FALSE, n)

  while (!all(excluded)) {
    candidate <- sample(which(!excluded), size = 1)
    first_hop <- which(W[candidate, ] != 0)
    second_hop <- unique(unlist(lapply(first_hop, function(j) which(W[j, ] != 0))))
    remove_set <- unique(c(candidate, first_hop, second_hop))
    net[candidate] <- TRUE
    excluded[remove_set] <- TRUE
  }

  clusters <- integer(n)
  net_nodes <- which(net)
  clusters[net_nodes] <- net_nodes

  non_net_nodes <- which(!net)
  if (length(non_net_nodes) > 0L) {
    for (i in non_net_nodes) {
      first_hop <- which(W[i, ] != 0)
      net_neighbors <- intersect(first_hop, net_nodes)
      if (length(net_neighbors) > 0L) {
        clusters[i] <- sample(net_neighbors, size = 1)
        next
      }

      second_hop_net <- unique(unlist(
        lapply(first_hop, function(j) intersect(which(W[j, ] != 0), net_nodes))
      ))
      clusters[i] <- if (length(second_hop_net) > 0L) {
        sample(second_hop_net, size = 1)
      } else {
        sample(net_nodes, size = 1)
      }
    }
  }

  list(clusters = clusters)
}

#' Louvain community detection
#'
#' @param W Symmetric binary adjacency matrix.
#'
#' @return Integer cluster-membership vector.
#' @export
cluster_louvain_membership <- function(W) {
  W <- validate_adjacency_matrix(W, binary = TRUE)
  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop("Package 'igraph' is required for Louvain clustering.", call. = FALSE)
  }
  g <- igraph::graph_from_adjacency_matrix(W, mode = "undirected", diag = FALSE, weighted = NULL)
  igraph::cluster_louvain(g)$membership
}

#' Spectral clustering baseline
#'
#' @param W Symmetric matrix used to construct the graph Laplacian.
#' @param num_clusters Number of clusters.
#' @param normalized Whether to use the normalized graph Laplacian.
#' @param n_eig Number of eigenvectors.
#' @param hierarchical If `TRUE`, use hierarchical clustering instead of k-means.
#' @param seed Optional k-means seed.
#' @param binary If `TRUE`, use the binary support of `W`.
#'
#' @return Integer cluster-membership vector.
#' @export
cluster_spectral <- function(W,
                             num_clusters,
                             normalized = TRUE,
                             n_eig = num_clusters,
                             hierarchical = FALSE,
                             seed = NULL,
                             binary = TRUE) {
  W <- validate_adjacency_matrix(W, binary = binary)
  n <- nrow(W)
  num_clusters <- as.integer(num_clusters)
  if (num_clusters < 1L || num_clusters > n) {
    stop("`num_clusters` must be between 1 and nrow(W).", call. = FALSE)
  }
  if (num_clusters == 1L) {
    return(rep(1L, n))
  }

  emb <- spectral_embedding(W, n_eig = n_eig, normalized = normalized, binary = FALSE)
  if (isTRUE(hierarchical)) {
    hc <- stats::hclust(stats::dist(emb), method = "ward.D2")
    return(stats::cutree(hc, k = num_clusters))
  }

  if (!is.null(seed)) set.seed(seed)
  stats::kmeans(emb, centers = num_clusters, iter.max = 50)$cluster
}

# Relaxations ---------------------------------------------------------------

spectral_relaxation_matrix <- function(W, xi) {
  W <- validate_adjacency_matrix(W, binary = TRUE)
  xi <- validate_xi(xi)
  n <- nrow(W)
  B <- matrix(xi, nrow = n, ncol = n) - left_normalized_adjacency(W)
  0.5 * (B + t(B))
}

symmetric_entry_matrix <- function(n, i, j) {
  M <- Matrix::Matrix(0, nrow = n, ncol = n, sparse = TRUE)
  if (i == j) {
    M[i, j] <- 1
  } else {
    M[i, j] <- 0.5
    M[j, i] <- 0.5
  }
  M
}

empty_linear_block <- function(m) {
  Matrix::Matrix(0, nrow = m, ncol = 1, sparse = TRUE)
}

linear_unit_block <- function(m, idx, value = 1) {
  out <- empty_linear_block(m)
  out[idx, 1] <- value
  out
}

sdp_maxcut <- function(B) {
  check_sdp_dependencies()
  if (!is.matrix(B) || !is.numeric(B) || nrow(B) != ncol(B)) {
    stop("`B` must be a square numeric matrix.", call. = FALSE)
  }
  B <- 0.5 * (B + t(B))
  n <- nrow(B)

  blk <- matrix(list(), 1, 2)
  blk[[1, 1]] <- "s"
  blk[[1, 2]] <- n

  C <- matrix(list(), 1, 1)
  C[[1]] <- B

  A <- matrix(list(), 1, n)
  b <- rep(1, n)
  for (k in seq_len(n)) {
    A[[1, k]] <- Matrix::Matrix(0, nrow = n, ncol = n, sparse = TRUE)
    A[[1, k]][k, k] <- 1
  }

  svec_fun <- utils::getFromNamespace("svec", "sdpt3r")
  sqlp_base_fun <- utils::getFromNamespace("sqlp_base", "sdpt3r")
  Avec <- svec_fun(blk, M = A, isspx = matrix(0, nrow(blk), 1))
  out <- sqlp_base_fun(blk, Avec, C, b)
  dim(out$X) <- NULL
  dim(out$Z) <- NULL
  out
}

sdp_objective_value <- function(solution, xi, n) {
  J <- matrix(1, nrow = n, ncol = n)
  X <- as.matrix(solution$X[[1]])
  Y <- as.matrix(solution$X[[2]])
  as.numeric((xi / n^2) * sum(J * X) + Y[1, 1])
}

#' SDP relaxation for Equation (11)
#'
#' @description Solves the semidefinite relaxation used by Algorithm 1.
#' The returned objective is in Equation (9) normalization:
#' `xi * tr(11'X) / n^2 + t`.
#'
#' @param W Symmetric binary adjacency matrix.
#' @param xi Non-negative calibration value.
#' @param num_clusters Number of clusters K; required when `k_constraint = TRUE`.
#' @param gamma_bar Maximum proportional cluster-size constant.
#' @param k_constraint Whether to impose the optional K-dependent constraints.
#' @param box_constraints Whether to impose `0 <= X_ij <= 1` when
#'   `k_constraint = TRUE`.
#'
#' @return An SDPT3 solution object with added entries `X_matrix`, `objective`,
#'   and `lower_bound`.
sdp_quadratic_relaxation <- function(W,
                                     xi,
                                     num_clusters = NULL,
                                     gamma_bar = 10,
                                     k_constraint = FALSE,
                                     box_constraints = TRUE) {
  check_sdp_dependencies()
  W <- validate_adjacency_matrix(W, binary = TRUE)
  xi <- validate_xi(xi)
  gamma_bar <- validate_gamma_bar(gamma_bar)
  box_constraints <- validate_box_constraints(box_constraints)

  n <- nrow(W)
  J <- matrix(1, nrow = n, ncol = n)
  L <- left_normalized_adjacency(W)
  L_sym <- 0.5 * (L + t(L))

  if (isTRUE(k_constraint)) {
    if (is.null(num_clusters)) {
      stop("`num_clusters` must be supplied when `k_constraint = TRUE`.", call. = FALSE)
    }
    num_clusters <- as.integer(num_clusters)
    if (num_clusters < 1L || num_clusters > n) {
      stop("`num_clusters` must be between 1 and nrow(W).", call. = FALSE)
    }
  }

  n_slacks <- 0L
  if (isTRUE(k_constraint)) {
    n_slacks <- n_slacks + 2L
    if (isTRUE(box_constraints)) {
      n_slacks <- n_slacks + 2L * as.integer(n * (n + 1L) / 2L)
    }
  }

  n_blocks <- if (n_slacks > 0L) 3L else 2L
  blk <- matrix(list(), n_blocks, 2)
  blk[[1, 1]] <- "s"; blk[[1, 2]] <- n
  blk[[2, 1]] <- "s"; blk[[2, 2]] <- 2
  if (n_slacks > 0L) {
    blk[[3, 1]] <- "l"; blk[[3, 2]] <- n_slacks
  }

  n_constraints <- n + 2L
  if (isTRUE(k_constraint)) {
    n_constraints <- n_constraints + 2L
    if (isTRUE(box_constraints)) {
      n_constraints <- n_constraints + 2L * as.integer(n * (n + 1L) / 2L)
    }
  }

  C <- matrix(list(), n_blocks, 1)
  A <- matrix(list(), n_blocks, n_constraints)
  b <- numeric(n_constraints)

  C[[1]] <- Matrix::Matrix((xi / n^2) * J, sparse = TRUE)
  C[[2]] <- Matrix::Matrix(c(1, 0, 0, 0), nrow = 2, ncol = 2, byrow = TRUE, sparse = TRUE)
  if (n_slacks > 0L) C[[3]] <- empty_linear_block(n_slacks)

  zero_X <- function() Matrix::Matrix(0, nrow = n, ncol = n, sparse = TRUE)
  zero_Y <- function() Matrix::Matrix(0, nrow = 2, ncol = 2, sparse = TRUE)
  zero_S <- function() empty_linear_block(n_slacks)

  for (k in seq_len(n)) {
    A[[1, k]] <- zero_X()
    A[[1, k]][k, k] <- 1
    A[[2, k]] <- zero_Y()
    if (n_slacks > 0L) A[[3, k]] <- zero_S()
    b[k] <- 1
  }

  row_y22 <- n + 1L
  A[[1, row_y22]] <- zero_X()
  A[[2, row_y22]] <- zero_Y()
  A[[2, row_y22]][2, 2] <- 1
  if (n_slacks > 0L) A[[3, row_y22]] <- zero_S()
  b[row_y22] <- 1

  row_z <- n + 2L
  A[[1, row_z]] <- Matrix::Matrix(L_sym / n, sparse = TRUE)
  A[[2, row_z]] <- zero_Y()
  A[[2, row_z]][1, 2] <- 0.5
  A[[2, row_z]][2, 1] <- 0.5
  if (n_slacks > 0L) A[[3, row_z]] <- zero_S()
  b[row_z] <- sum(L) / n

  row <- row_z
  slack <- 0L

  if (isTRUE(k_constraint)) {
    lower <- n^2 / num_clusters
    upper <- gamma_bar * n^2 / num_clusters

    row <- row + 1L; slack <- slack + 1L
    A[[1, row]] <- Matrix::Matrix(J, sparse = TRUE)
    A[[2, row]] <- zero_Y()
    A[[3, row]] <- linear_unit_block(n_slacks, slack, -1)
    b[row] <- lower

    row <- row + 1L; slack <- slack + 1L
    A[[1, row]] <- Matrix::Matrix(J, sparse = TRUE)
    A[[2, row]] <- zero_Y()
    A[[3, row]] <- linear_unit_block(n_slacks, slack, 1)
    b[row] <- upper

    if (isTRUE(box_constraints)) {
      for (i in seq_len(n)) {
        for (j in i:n) {
          Eij <- symmetric_entry_matrix(n, i, j)

          row <- row + 1L; slack <- slack + 1L
          A[[1, row]] <- Eij
          A[[2, row]] <- zero_Y()
          A[[3, row]] <- linear_unit_block(n_slacks, slack, -1)
          b[row] <- 0

          row <- row + 1L; slack <- slack + 1L
          A[[1, row]] <- Eij
          A[[2, row]] <- zero_Y()
          A[[3, row]] <- linear_unit_block(n_slacks, slack, 1)
          b[row] <- 1
        }
      }
    }
  }

  svec_fun <- utils::getFromNamespace("svec", "sdpt3r")
  sqlp_base_fun <- utils::getFromNamespace("sqlp_base", "sdpt3r")
  Avec <- svec_fun(blk, M = A, isspx = matrix(0, nrow(blk), 1))
  out <- sqlp_base_fun(blk, Avec, C, b)

  dim(out$X) <- NULL
  dim(out$Z) <- NULL
  out$X_matrix <- as.matrix(out$X[[1]])
  out$objective <- sdp_objective_value(out, xi, n)
  out$lower_bound <- out$objective
  out$xi <- xi
  out$num_clusters <- num_clusters
  out$k_constraint <- isTRUE(k_constraint)
  out$gamma_bar <- gamma_bar
  out$box_constraints <- isTRUE(box_constraints)
  out
}

#' SDP endpoint-regret relaxation for Algorithm 2
#'
#' @param W Symmetric binary adjacency matrix.
#' @param xi_lower Positive lower endpoint of the xi range.
#' @param xi_upper Positive upper endpoint of the xi range.
#' @param lower_bound_lower SDP lower bound at `xi_lower`.
#' @param lower_bound_upper SDP lower bound at `xi_upper`.
#'
#' @return An SDPT3 solution object with entries `X_matrix` and `rho`.
sdp_regret_relaxation <- function(W,
                                  xi_lower,
                                  xi_upper,
                                  lower_bound_lower,
                                  lower_bound_upper) {
  check_sdp_dependencies()
  W <- validate_adjacency_matrix(W, binary = TRUE)
  xi_lower <- validate_positive_xi(xi_lower, "xi_lower")
  xi_upper <- validate_positive_xi(xi_upper, "xi_upper")

  if (xi_lower > xi_upper) {
    tmp <- xi_lower; xi_lower <- xi_upper; xi_upper <- tmp
    tmp <- lower_bound_lower; lower_bound_lower <- lower_bound_upper; lower_bound_upper <- tmp
  }
  if (!is.numeric(lower_bound_lower) || length(lower_bound_lower) != 1L ||
      !is.finite(lower_bound_lower) || lower_bound_lower <= 0) {
    stop("`lower_bound_lower` must be a positive finite number.", call. = FALSE)
  }
  if (!is.numeric(lower_bound_upper) || length(lower_bound_upper) != 1L ||
      !is.finite(lower_bound_upper) || lower_bound_upper <= 0) {
    stop("`lower_bound_upper` must be a positive finite number.", call. = FALSE)
  }

  n <- nrow(W)
  J <- matrix(1, nrow = n, ncol = n)
  L <- left_normalized_adjacency(W)
  L_sym <- 0.5 * (L + t(L))

  blk <- matrix(list(), 3, 2)
  blk[[1, 1]] <- "s"; blk[[1, 2]] <- n
  blk[[2, 1]] <- "s"; blk[[2, 2]] <- 2
  blk[[3, 1]] <- "l"; blk[[3, 2]] <- 3

  n_constraints <- n + 4L
  C <- matrix(list(), 3, 1)
  A <- matrix(list(), 3, n_constraints)
  b <- numeric(n_constraints)

  C[[1]] <- Matrix::Matrix(0, nrow = n, ncol = n, sparse = TRUE)
  C[[2]] <- Matrix::Matrix(0, nrow = 2, ncol = 2, sparse = TRUE)
  C[[3]] <- Matrix::Matrix(c(1, 0, 0), nrow = 3, ncol = 1, sparse = TRUE)

  zero_X <- function() Matrix::Matrix(0, nrow = n, ncol = n, sparse = TRUE)
  zero_Y <- function() Matrix::Matrix(0, nrow = 2, ncol = 2, sparse = TRUE)
  zero_S <- function() Matrix::Matrix(0, nrow = 3, ncol = 1, sparse = TRUE)

  for (k in seq_len(n)) {
    A[[1, k]] <- zero_X()
    A[[1, k]][k, k] <- 1
    A[[2, k]] <- zero_Y()
    A[[3, k]] <- zero_S()
    b[k] <- 1
  }

  row_y22 <- n + 1L
  A[[1, row_y22]] <- zero_X()
  A[[2, row_y22]] <- zero_Y()
  A[[2, row_y22]][2, 2] <- 1
  A[[3, row_y22]] <- zero_S()
  b[row_y22] <- 1

  row_z <- n + 2L
  A[[1, row_z]] <- Matrix::Matrix(L_sym / n, sparse = TRUE)
  A[[2, row_z]] <- zero_Y()
  A[[2, row_z]][1, 2] <- 0.5
  A[[2, row_z]][2, 1] <- 0.5
  A[[3, row_z]] <- zero_S()
  b[row_z] <- sum(L) / n

  row_low <- n + 3L
  A[[1, row_low]] <- Matrix::Matrix(-(xi_lower / n^2) * J, sparse = TRUE)
  A[[2, row_low]] <- Matrix::Matrix(c(-1, 0, 0, 0), nrow = 2, ncol = 2, byrow = TRUE, sparse = TRUE)
  A[[3, row_low]] <- Matrix::Matrix(c(lower_bound_lower, -1, 0), nrow = 3, ncol = 1, sparse = TRUE)
  b[row_low] <- 0

  row_up <- n + 4L
  A[[1, row_up]] <- Matrix::Matrix(-(xi_upper / n^2) * J, sparse = TRUE)
  A[[2, row_up]] <- Matrix::Matrix(c(-1, 0, 0, 0), nrow = 2, ncol = 2, byrow = TRUE, sparse = TRUE)
  A[[3, row_up]] <- Matrix::Matrix(c(lower_bound_upper, 0, -1), nrow = 3, ncol = 1, sparse = TRUE)
  b[row_up] <- 0

  svec_fun <- utils::getFromNamespace("svec", "sdpt3r")
  sqlp_base_fun <- utils::getFromNamespace("sqlp_base", "sdpt3r")
  Avec <- svec_fun(blk, M = A, isspx = matrix(0, nrow(blk), 1))
  out <- sqlp_base_fun(blk, Avec, C, b)

  dim(out$X) <- NULL
  dim(out$Z) <- NULL
  out$X_matrix <- as.matrix(out$X[[1]])
  out$rho <- as.numeric(out$X[[3]][1, 1])
  out$objective <- out$rho
  out$xi_lower <- xi_lower
  out$xi_upper <- xi_upper
  out$lower_bound_lower <- lower_bound_lower
  out$lower_bound_upper <- lower_bound_upper
  out
}

sdp_relaxation_solution <- function(W,
                                    xi,
                                    objective_type = "squared",
                                    num_clusters = NULL,
                                    gamma_bar = 10,
                                    k_constraint = FALSE,
                                    box_constraints = TRUE) {
  W <- validate_adjacency_matrix(W, binary = TRUE)
  xi <- validate_xi(xi)
  validate_objective_type(objective_type)

  sdp_quadratic_relaxation(
    W = W,
    xi = xi,
    num_clusters = num_clusters,
    gamma_bar = gamma_bar,
    k_constraint = k_constraint,
    box_constraints = box_constraints
  )
}

sdp_relaxation_matrix <- function(W,
                                  xi,
                                  objective_type = "squared",
                                  num_clusters = NULL,
                                  gamma_bar = 10,
                                  k_constraint = FALSE,
                                  box_constraints = TRUE) {
  sdp_relaxation_solution(
    W = W,
    xi = xi,
    objective_type = objective_type,
    num_clusters = num_clusters,
    gamma_bar = gamma_bar,
    k_constraint = k_constraint,
    box_constraints = box_constraints
  )$X_matrix
}

search_matrix <- function(W,
                          xi,
                          engine = c("sdp", "spectral"),
                          objective_type = "squared",
                          num_clusters = NULL,
                          gamma_bar = 10,
                          k_constraint = FALSE,
                          box_constraints = TRUE) {
  W <- validate_adjacency_matrix(W, binary = TRUE)
  xi <- validate_xi(xi)
  engine <- match.arg(engine)
  objective_type <- validate_objective_type(objective_type)

  if (identical(engine, "spectral")) {
    if (isTRUE(k_constraint)) {
      stop("`k_constraint = TRUE` requires `engine = 'sdp'`.", call. = FALSE)
    }
    return(spectral_relaxation_matrix(W, xi))
  }

  sdp_relaxation_matrix(
    W = W,
    xi = xi,
    objective_type = objective_type,
    num_clusters = num_clusters,
    gamma_bar = gamma_bar,
    k_constraint = k_constraint,
    box_constraints = box_constraints
  )
}

# Rounding and objective helpers -------------------------------------------

embedding_from_search_matrix <- function(my_matrix,
                                         n_eig,
                                         engine = c("sdp", "spectral")) {
  engine <- match.arg(engine)
  n <- nrow(my_matrix)
  n_eig <- max(1L, min(as.integer(n_eig), n))

  eig <- if (identical(engine, "spectral")) {
    eigen(-my_matrix, symmetric = TRUE)
  } else {
    eigen(my_matrix, symmetric = TRUE)
  }

  eig$vectors[, seq_len(n_eig), drop = FALSE]
}

cluster_from_search_matrix <- function(my_matrix,
                                       num_clusters,
                                       n_eig = num_clusters,
                                       engine = c("sdp", "spectral"),
                                       method = c("kmeans", "hierarchical"),
                                       seed = NULL) {
  engine <- match.arg(engine)
  method <- match.arg(method)
  num_clusters <- as.integer(num_clusters)
  n <- nrow(my_matrix)

  if (num_clusters < 1L || num_clusters > n) {
    stop("`num_clusters` must be between 1 and nrow(my_matrix).", call. = FALSE)
  }
  if (num_clusters == 1L) {
    return(rep(1L, n))
  }

  emb <- embedding_from_search_matrix(my_matrix, n_eig = n_eig, engine = engine)

  if (identical(method, "hierarchical")) {
    hc <- stats::hclust(stats::dist(emb), method = "ward.D2")
    return(stats::cutree(hc, k = num_clusters))
  }

  if (!is.null(seed)) set.seed(seed)
  stats::kmeans(emb, centers = num_clusters, iter.max = 50)$cluster
}

is_valid_cluster_vector <- function(clusters, n) {
  length(clusters) == n && !anyNA(clusters) && length(unique(clusters)) >= 1L
}

cluster_size_constraint_ok <- function(clusters,
                                       num_clusters,
                                       n,
                                       gamma_bar,
                                       tol = 1e-8) {
  if (!is_valid_cluster_vector(clusters, n)) return(FALSE)
  cluster_ids <- unique(clusters)
  if (length(cluster_ids) != as.integer(num_clusters)) return(FALSE)
  sizes <- tabulate(match(clusters, cluster_ids), nbins = length(cluster_ids))
  max(sizes) <= gamma_bar * n / as.integer(num_clusters) + tol
}

#' Compute objective components for a clustering
#'
#' @description Computes `variance = n^{-2} sum_k n_k^2` and
#' `bias = n^{-1} sum_i |N_i|^{-1} |{j in N_i: c(i) != c(j)}|`.
#'
#' @param W Symmetric binary adjacency matrix.
#' @param clusters Cluster labels.
#'
#' @return List with `variance`, `bias`, `num_clusters`, and `cluster_sizes`.
#' @export
clustering_objective_components <- function(W, clusters) {
  W <- validate_adjacency_matrix(W, binary = TRUE)
  n <- nrow(W)

  if (length(clusters) != n) {
    stop("`clusters` must have length equal to nrow(W).", call. = FALSE)
  }
  if (anyNA(clusters)) {
    stop("`clusters` cannot contain missing values.", call. = FALSE)
  }

  cluster_ids <- unique(clusters)
  cluster_index <- match(clusters, cluster_ids)
  cluster_sizes <- as.integer(tabulate(cluster_index, nbins = length(cluster_ids)))
  names(cluster_sizes) <- as.character(cluster_ids)

  adjacency <- W != 0
  diag(adjacency) <- FALSE
  deg <- pmax(rowSums(adjacency), 1)

  bias_sum <- 0
  for (i in seq_len(n)) {
    neigh <- which(adjacency[i, ])
    if (length(neigh) > 0L) {
      bias_sum <- bias_sum + sum(clusters[neigh] != clusters[i]) / deg[i]
    }
  }

  list(
    variance = sum(cluster_sizes^2) / n^2,
    bias = bias_sum / n,
    num_clusters = length(cluster_ids),
    cluster_sizes = cluster_sizes
  )
}

#' Compute the causal-clustering objective
#'
#' @description Computes the design objective
#' `xi * n^{-2} sum_k n_k^2 + b_n(C)^2`.
#'
#' @param W Symmetric binary adjacency matrix.
#' @param xi Non-negative tuning parameter.
#' @param clusters Cluster labels.
#' @param objective_type Must be `"squared"`.
#'
#' @return Numeric scalar.
#' @export
compute_causal_clustering_objective <- function(W,
                                                xi,
                                                clusters,
                                                objective_type = "squared") {
  xi <- validate_xi(xi)
  validate_objective_type(objective_type)
  components <- clustering_objective_components(W, clusters)
  xi * components$variance + components$bias^2
}

#' Compute the clustering objective
#'
#' @description Wrapper around [compute_causal_clustering_objective()].
#'
#' @param W Symmetric binary adjacency matrix.
#' @param xi Non-negative tuning parameter.
#' @param clusters Cluster labels with length `nrow(W)`.
#' @param squared_bias Logical. Must be `TRUE`.
#' @param objective_type Must be `"squared"`.
#'
#' @return A numeric scalar objective value.
#' @export
compute_objective <- function(W,
                              xi,
                              clusters,
                              squared_bias = TRUE,
                              objective_type = NULL) {
  if (!isTRUE(squared_bias)) {
    stop("`squared_bias` must be TRUE for the design criterion.", call. = FALSE)
  }
  objective_type <- if (is.null(objective_type)) "squared" else objective_type
  objective_type <- validate_objective_type(objective_type)
  compute_causal_clustering_objective(W, xi, clusters, objective_type = objective_type)
}

safe_candidate_objective <- function(W,
                                     xi,
                                     clusters,
                                     objective_type,
                                     num_clusters = NULL,
                                     k_constraint = FALSE,
                                     gamma_bar = 10) {
  n <- nrow(W)
  if (!is_valid_cluster_vector(clusters, n)) return(Inf)
  if (!is.null(num_clusters) && length(unique(clusters)) != as.integer(num_clusters)) return(Inf)
  if (isTRUE(k_constraint) && !cluster_size_constraint_ok(clusters, num_clusters, n, gamma_bar)) return(Inf)

  tryCatch(
    compute_causal_clustering_objective(W, xi, clusters, objective_type = objective_type),
    error = function(e) Inf
  )
}

candidate_clusterings_from_matrix <- function(W,
                                              my_matrix,
                                              num_clusters,
                                              n_eig = num_clusters,
                                              engine = c("sdp", "spectral"),
                                              methods = available_discretization_methods(),
                                              seed = NULL,
                                              try_sign_flip = TRUE,
                                              k_constraint = FALSE,
                                              gamma_bar = 10) {
  W <- validate_adjacency_matrix(W, binary = TRUE)
  engine <- match.arg(engine)
  methods <- standardize_discretization_methods(methods)
  gamma_bar <- validate_gamma_bar(gamma_bar)

  n <- nrow(W)
  num_clusters <- as.integer(num_clusters)
  candidates <- list()

  add_candidate <- function(name, fun) {
    clusters <- tryCatch(fun(), error = function(e) NULL)
    if (is.null(clusters) || !is_valid_cluster_vector(clusters, n)) return(invisible(NULL))
    if (length(unique(clusters)) != num_clusters) return(invisible(NULL))
    if (isTRUE(k_constraint) && !cluster_size_constraint_ok(clusters, num_clusters, n, gamma_bar)) {
      return(invisible(NULL))
    }
    candidates[[name]] <<- clusters
    invisible(NULL)
  }

  signs <- if (isTRUE(try_sign_flip)) c(plus = 1, minus = -1) else c(plus = 1)

  for (sign_name in names(signs)) {
    signed_matrix <- signs[[sign_name]] * my_matrix

    if ("kmeans" %in% methods) {
      add_candidate(paste0("search_kmeans_", sign_name), function() {
        cluster_from_search_matrix(signed_matrix, num_clusters, n_eig, engine, "kmeans", seed)
      })
    }
    if ("hierarchical" %in% methods) {
      add_candidate(paste0("search_hierarchical_", sign_name), function() {
        cluster_from_search_matrix(signed_matrix, num_clusters, n_eig, engine, "hierarchical", seed)
      })
    }
    if ("spectral_norm_kmeans" %in% methods) {
      add_candidate(paste0("spectral_norm_kmeans_", sign_name), function() {
        cluster_spectral(signed_matrix, num_clusters, normalized = TRUE, n_eig = n_eig, hierarchical = FALSE, seed = seed, binary = FALSE)
      })
    }
    if ("spectral_unnorm_kmeans" %in% methods) {
      add_candidate(paste0("spectral_unnorm_kmeans_", sign_name), function() {
        cluster_spectral(signed_matrix, num_clusters, normalized = FALSE, n_eig = n_eig, hierarchical = FALSE, seed = seed, binary = FALSE)
      })
    }
    if ("spectral_norm_hierarchical" %in% methods) {
      add_candidate(paste0("spectral_norm_hierarchical_", sign_name), function() {
        cluster_spectral(signed_matrix, num_clusters, normalized = TRUE, n_eig = n_eig, hierarchical = TRUE, seed = seed, binary = FALSE)
      })
    }
    if ("spectral_unnorm_hierarchical" %in% methods) {
      add_candidate(paste0("spectral_unnorm_hierarchical_", sign_name), function() {
        cluster_spectral(signed_matrix, num_clusters, normalized = FALSE, n_eig = n_eig, hierarchical = TRUE, seed = seed, binary = FALSE)
      })
    }
  }

  candidates
}

best_search_clustering <- function(W,
                                   my_matrix,
                                   xi,
                                   num_clusters,
                                   n_eig = num_clusters,
                                   engine = c("sdp", "spectral"),
                                   methods = available_discretization_methods(),
                                   seed = NULL,
                                   objective_type = "squared",
                                   try_sign_flip = TRUE,
                                   k_constraint = FALSE,
                                   gamma_bar = 10) {
  W <- validate_adjacency_matrix(W, binary = TRUE)
  xi <- validate_xi(xi)
  engine <- match.arg(engine)
  methods <- standardize_discretization_methods(methods)
  objective_type <- validate_objective_type(objective_type)
  gamma_bar <- validate_gamma_bar(gamma_bar)

  candidates <- candidate_clusterings_from_matrix(
    W = W,
    my_matrix = my_matrix,
    num_clusters = num_clusters,
    n_eig = n_eig,
    engine = engine,
    methods = methods,
    seed = seed,
    try_sign_flip = try_sign_flip,
    k_constraint = k_constraint,
    gamma_bar = gamma_bar
  )

  if (length(candidates) == 0L) {
    stop("No valid clustering candidates were produced for K = ", num_clusters, ".", call. = FALSE)
  }

  objectives <- vapply(
    candidates,
    function(cl) safe_candidate_objective(
      W = W,
      xi = xi,
      clusters = cl,
      objective_type = objective_type,
      num_clusters = num_clusters,
      k_constraint = k_constraint,
      gamma_bar = gamma_bar
    ),
    numeric(1)
  )
  best_name <- names(which.min(objectives))

  list(
    clusters = candidates[[best_name]],
    objective = unname(objectives[[best_name]]),
    method = best_name,
    candidates = candidates,
    candidate_objectives = objectives,
    objective_type = objective_type
  )
}

# Algorithm 1 ---------------------------------------------------------------

#' Search causal clustering over K for fixed xi
#'
#' @description Implements Algorithm 1 for a fixed value of `xi`. With
#' `k_constraint = FALSE`, the relaxation is solved once and rounded over the
#' K grid. With `k_constraint = TRUE`, a separate K-specific SDP with optional
#' constraints is solved for every K.
#'
#' @param W Symmetric binary adjacency matrix.
#' @param xi Non-negative tuning parameter.
#' @param min_k Minimum K.
#' @param max_k Maximum K.
#' @param n_eig Optional number of eigenvectors used for rounding.
#' @param seed Optional k-means seed.
#' @param engine `"sdp"` is the SDP relaxation; `"spectral"` is a heuristic.
#' @param methods Character vector of rounding/discretization methods. All requested methods are run and the rounded candidate with the smallest realized objective is selected.
#' @param include_bernoulli Whether to include the all-singleton design.
#' @param objective_type Must be `"squared"` for Equation (9).
#' @param k_constraint Boolean from Algorithm 1.
#' @param gamma_bar Cluster-size constant for optional constraints.
#' @param box_constraints Whether to impose `0 <= X_ij <= 1` in the K-specific SDP.
#' @param try_sign_flip Whether to try rounding both X and -X. This is useful when comparing multiple discretizations after a relaxation step.
#' @param keep_sdp_solutions Whether to return full SDP objects.
#'
#' @return List containing the selected clustering, candidate objectives, SDP
#'   lower bounds, and `Gamma_n` when available.
#' @export
search_causal_clustering <- function(W,
                                     xi = 0,
                                     min_k = 2,
                                     max_k = floor(nrow(W) / 2),
                                     n_eig = NULL,
                                     seed = NULL,
                                     engine = c("sdp", "spectral"),
                                     methods = available_discretization_methods(),
                                     include_bernoulli = FALSE,
                                     objective_type = "squared",
                                     k_constraint = FALSE,
                                     gamma_bar = 10,
                                     box_constraints = TRUE,
                                     try_sign_flip = TRUE,
                                     keep_sdp_solutions = FALSE) {
  W <- validate_adjacency_matrix(W, binary = TRUE)
  xi <- validate_xi(xi)
  engine <- match.arg(engine)
  methods <- standardize_discretization_methods(methods)
  objective_type <- validate_objective_type(objective_type)
  gamma_bar <- validate_gamma_bar(gamma_bar)
  box_constraints <- validate_box_constraints(box_constraints)

  if (isTRUE(k_constraint) && identical(engine, "spectral")) {
    stop("`k_constraint = TRUE` requires `engine = 'sdp'`.", call. = FALSE)
  }
  n <- nrow(W)
  k_grid <- validate_k_range(min_k, max_k, n)
  solutions <- vector("list", length(k_grid))
  objectives <- rep(Inf, length(k_grid))
  sdp_lower_bounds <- rep(NA_real_, length(k_grid))
  sdp_solutions <- vector("list", length(k_grid))

  shared_search_mat <- NULL
  shared_sdp_solution <- NULL
  if (!isTRUE(k_constraint)) {
    if (identical(engine, "sdp")) {
      shared_sdp_solution <- sdp_relaxation_solution(
        W = W,
        xi = xi,
        objective_type = objective_type,
        k_constraint = FALSE
      )
      shared_search_mat <- shared_sdp_solution$X_matrix
      sdp_lower_bounds[] <- shared_sdp_solution$objective
    } else {
      shared_search_mat <- spectral_relaxation_matrix(W, xi)
    }
  }

  for (i in seq_along(k_grid)) {
    k <- k_grid[i]
    eig_k <- if (is.null(n_eig)) k else n_eig

    if (isTRUE(k_constraint) && identical(engine, "sdp")) {
      sdp_sol <- sdp_relaxation_solution(
        W = W,
        xi = xi,
        objective_type = objective_type,
        num_clusters = k,
        gamma_bar = gamma_bar,
        k_constraint = TRUE,
        box_constraints = box_constraints
      )
      search_mat <- sdp_sol$X_matrix
      sdp_lower_bounds[i] <- sdp_sol$objective
      if (isTRUE(keep_sdp_solutions)) sdp_solutions[[i]] <- sdp_sol
    } else {
      search_mat <- shared_search_mat
      if (isTRUE(keep_sdp_solutions) && !is.null(shared_sdp_solution)) {
        sdp_solutions[[i]] <- shared_sdp_solution
      }
    }

    sol <- tryCatch(
      best_search_clustering(
        W = W,
        my_matrix = search_mat,
        xi = xi,
        num_clusters = k,
        n_eig = eig_k,
        engine = engine,
        methods = methods,
        seed = seed,
        objective_type = objective_type,
        try_sign_flip = try_sign_flip,
        k_constraint = isTRUE(k_constraint),
        gamma_bar = gamma_bar
      ),
      error = function(e) list(error = conditionMessage(e), clusters = NULL, objective = Inf)
    )

    solutions[[i]] <- c(sol, list(num_clusters = k))
    objectives[i] <- sol$objective
  }

  objective_names <- paste0("k_", k_grid)
  all_objectives <- objectives
  bernoulli_clusters <- NULL

  if (isTRUE(include_bernoulli)) {
    bernoulli_clusters <- seq_len(n)
    bernoulli_objective <- compute_causal_clustering_objective(
      W = W,
      xi = xi,
      clusters = bernoulli_clusters,
      objective_type = objective_type
    )
    all_objectives <- c(objectives, bernoulli_objective)
    objective_names <- c(objective_names, "bernoulli")
  }

  if (all(!is.finite(all_objectives))) {
    stop("No finite objective was obtained for any candidate clustering.", call. = FALSE)
  }

  best_index <- which.min(all_objectives)
  best_clusters <- if (best_index <= length(k_grid)) solutions[[best_index]]$clusters else bernoulli_clusters
  selected_solution <- if (best_index <= length(k_grid)) solutions[[best_index]] else NULL
  selected_method <- if (!is.null(selected_solution) && !is.null(selected_solution$method)) {
    selected_solution$method
  } else if (best_index > length(k_grid)) {
    "bernoulli"
  } else {
    NA_character_
  }
  selected_candidate_objectives <- if (!is.null(selected_solution) &&
      !is.null(selected_solution$candidate_objectives)) {
    selected_solution$candidate_objectives
  } else {
    NULL
  }

  lower_bound <- suppressWarnings(min(sdp_lower_bounds, na.rm = TRUE))
  if (!is.finite(lower_bound)) lower_bound <- NA_real_
  gamma_n <- if (is.na(lower_bound) || lower_bound <= 0) NA_real_ else unname(all_objectives[best_index]) / lower_bound

  certificate_valid <- identical(engine, "sdp") &&
    identical(objective_type, "squared") &&
    !isTRUE(include_bernoulli) &&
    is.finite(gamma_n)

  list(
    clusters = best_clusters,
    solutions = solutions,
    objectives = stats::setNames(all_objectives, objective_names),
    k_grid = k_grid,
    best_index = best_index,
    selected_k = if (best_index <= length(k_grid)) k_grid[best_index] else n,
    selected_method = selected_method,
    selected_candidate_objectives = selected_candidate_objectives,
    objective = unname(all_objectives[best_index]),
    objective_type = objective_type,
    engine = engine,
    methods = methods,
    include_bernoulli = include_bernoulli,
    k_constraint = isTRUE(k_constraint),
    gamma_bar = gamma_bar,
    box_constraints = isTRUE(box_constraints),
    try_sign_flip = try_sign_flip,
    search_matrix = if (isTRUE(k_constraint)) NULL else shared_search_mat,
    sdp_lower_bounds = stats::setNames(sdp_lower_bounds, paste0("k_", k_grid)),
    sdp_lower_bound = lower_bound,
    approximation_error = gamma_n,
    Gamma_n = gamma_n,
    certificate_valid = certificate_valid,
    sdp_solutions = if (isTRUE(keep_sdp_solutions)) sdp_solutions else NULL
  )
}

#' Algorithm 1: causal clustering for one calibration value
#'
#' @description
#' Implements Algorithm 1 for a single value of `xi`. The `k_constraint`
#' argument implements the Boolean branch in Algorithm 1: `FALSE` solves one relaxation and rounds over
#' all K values, while `TRUE` solves a K-specific SDP relaxation for every K.
#'
#' @param W Symmetric binary adjacency matrix.
#' @param xi Non-negative Algorithm 1 calibration, equal to
#'   `(lambda * phibar_n^2 / psibar)^(-1)` in the paper.
#' @param min_k Minimum number of clusters to consider.
#' @param max_k Maximum number of clusters to consider.
#' @param n_eig Optional number of eigenvectors used for rounding. If `NULL`,
#'   the number of eigenvectors equals the candidate K.
#' @param seed Optional seed used by k-means rounding and randomized baselines.
#' @param engine Optimization engine. Use `"sdp"` for the semidefinite relaxation and
#'   `"spectral"` for a spectral heuristic.
#' @param methods Character vector of rounding/discretization methods. See
#'   [available_discretization_methods()]. All requested methods are run and the rounded candidate with the smallest realized objective is selected.
#' @param include_bernoulli Logical; if `TRUE`, also compare against the
#'   all-singleton Bernoulli design. This option is outside the approximation
#'   certificate in Algorithm 1.
#' @param objective_type Must be `"squared"` for Equation (9).
#' @param k_constraint Logical Boolean from Algorithm 1.
#' @param gamma_bar Upper cluster-size multiplier used by the optional
#'   K-specific constraints.
#' @param box_constraints Logical; whether to impose the optional constraints
#'   `0 <= X_ij <= 1` when `k_constraint = TRUE`.
#' @param try_sign_flip Logical; if `TRUE`, round both the relaxation matrix
#'   and its negative.
#' @param keep_search_matrices Logical; if `TRUE`, keep the relaxation matrix in
#'   the returned object.
#' @param keep_sdp_solutions Logical; if `TRUE`, keep raw SDP solution objects.
#' @param xi_grid Optional scalar alias for `xi`. Use
#'   [causal_clustering_algorithm2()] for calibration ranges.
#' @param calibration Optional positive scalar equal to
#'   `lambda * phibar_n^2 / psibar`; converted to `xi = 1 / calibration`.
#'
#' @return A list containing the selected `clusters`, selected K, objective
#'   value, objective components, SDP lower bounds, and the approximation
#'   certificate `Gamma_n` when available.
#' @export
causal_clustering_algorithm1 <- function(W,
                                         xi = NULL,
                                         min_k = 2,
                                         max_k = floor(nrow(W) / 2),
                                         n_eig = NULL,
                                         seed = NULL,
                                         engine = c("sdp", "spectral"),
                                         methods = available_discretization_methods(),
                                         include_bernoulli = FALSE,
                                         objective_type = "squared",
                                         k_constraint = FALSE,
                                         gamma_bar = 10,
                                         box_constraints = TRUE,
                                         try_sign_flip = TRUE,
                                         keep_search_matrices = FALSE,
                                         keep_sdp_solutions = FALSE,
                                         xi_grid = NULL,
                                         calibration = NULL) {
  W <- validate_adjacency_matrix(W, binary = TRUE)
  engine <- match.arg(engine)
  methods <- standardize_discretization_methods(methods)
  objective_type <- validate_objective_type(objective_type)

  if (is.null(xi)) {
    if (!is.null(calibration)) {
      calibration <- validate_positive_grid(calibration, "calibration")
      if (length(calibration) != 1L) {
        stop("`calibration` must be a single positive number when used to define `xi = 1 / calibration`.", call. = FALSE)
      }
      xi <- 1 / calibration
    } else if (!is.null(xi_grid) && length(xi_grid) == 1L) {
      xi <- xi_grid
    } else {
      stop("Algorithm 1 requires a single `xi`. For a range, use `causal_clustering_algorithm2()`.", call. = FALSE)
    }
  }
  xi <- validate_xi(xi)

  if (!is.null(xi_grid) && length(xi_grid) > 1L) {
    stop("Algorithm 1 does not optimize over an xi grid; use `causal_clustering_algorithm2()`.", call. = FALSE)
  }

  sol <- search_causal_clustering(
    W = W,
    xi = xi,
    min_k = min_k,
    max_k = max_k,
    n_eig = n_eig,
    seed = seed,
    engine = engine,
    methods = methods,
    include_bernoulli = include_bernoulli,
    objective_type = objective_type,
    k_constraint = k_constraint,
    gamma_bar = gamma_bar,
    box_constraints = box_constraints,
    try_sign_flip = try_sign_flip,
    keep_sdp_solutions = keep_sdp_solutions
  )

  components <- clustering_objective_components(W, sol$clusters)
  objective <- compute_causal_clustering_objective(W, xi, sol$clusters, objective_type = objective_type)
  if (!isTRUE(keep_search_matrices)) sol$search_matrix <- NULL

  list(
    clusters = sol$clusters,
    xi = xi,
    selected_k = sol$selected_k,
    selected_method = sol$selected_method,
    selected_candidate_objectives = sol$selected_candidate_objectives,
    best_index = sol$best_index,
    objective = objective,
    objective_type = objective_type,
    components = components,
    search_result = sol,
    engine = engine,
    methods = methods,
    include_bernoulli = include_bernoulli,
    k_constraint = isTRUE(k_constraint),
    gamma_bar = gamma_bar,
    box_constraints = isTRUE(box_constraints),
    sdp_lower_bounds = sol$sdp_lower_bounds,
    sdp_lower_bound = sol$sdp_lower_bound,
    approximation_error = sol$approximation_error,
    Gamma_n = sol$Gamma_n,
    certificate_valid = sol$certificate_valid,
    min_k = min(sol$k_grid),
    max_k = max(sol$k_grid)
  )
}

# Algorithm 2 ---------------------------------------------------------------

safe_ratio <- function(numerator, denominator) {
  if (!is.numeric(denominator) || length(denominator) != 1L || is.na(denominator)) {
    stop("`denominator` must be a single numeric value.", call. = FALSE)
  }
  if (denominator > 0) return(numerator / denominator)
  ifelse(numerator == 0, 1, Inf)
}

canonicalize_clusters <- function(clusters) {
  as.integer(match(clusters, unique(clusters)))
}

#' Algorithm 2: causal clustering over a range of spillover calibrations
#'
#' @description Implements Algorithm 2. It solves endpoint SDP lower
#' bounds at the two endpoints of the xi range, solves the endpoint-regret SDP,
#' rounds the resulting matrix for each K, and selects the clustering with the
#' smallest endpoint regret.
#'
#' @param W Symmetric binary adjacency matrix.
#' @param xi_lower Positive lower endpoint for xi.
#' @param xi_upper Positive upper endpoint for xi.
#' @param xi_range Optional length-2 positive range for xi.
#' @param calibration_range Optional length-2 range for
#'   `lambda * phibar_n^2 / psibar`; converted internally to `xi = 1/calibration`.
#' @param min_k Minimum K.
#' @param max_k Maximum K.
#' @param n_eig Optional number of eigenvectors used for rounding.
#' @param seed Optional k-means seed.
#' @param engine Must be `"sdp"` for Algorithm 2.
#' @param methods Character vector of rounding/discretization methods. All requested methods are run and the rounded candidate with the smallest realized objective is selected.
#' @param include_bernoulli Whether to include the all-singleton design.
#' @param k_constraint Whether to filter rounded candidates by the size constraint.
#' @param gamma_bar Cluster-size constant for optional filtering.
#' @param try_sign_flip Whether to try rounding both X and -X. This is useful when comparing multiple discretizations after a relaxation step.
#' @param keep_sdp_solutions Whether to return full SDP solution objects.
#'
#' @return List containing the selected clustering, endpoint lower bounds, and
#'   endpoint regret values.
#' @export
causal_clustering_algorithm2 <- function(W,
                                         xi_lower = NULL,
                                         xi_upper = NULL,
                                         xi_range = NULL,
                                         calibration_range = NULL,
                                         min_k = 2,
                                         max_k = floor(nrow(W) / 2),
                                         n_eig = NULL,
                                         seed = NULL,
                                         engine = "sdp",
                                         methods = available_discretization_methods(),
                                         include_bernoulli = FALSE,
                                         k_constraint = FALSE,
                                         gamma_bar = 10,
                                         try_sign_flip = TRUE,
                                         keep_sdp_solutions = FALSE) {
  W <- validate_adjacency_matrix(W, binary = TRUE)
  engine <- match.arg(engine)
  methods <- standardize_discretization_methods(methods)
  gamma_bar <- validate_gamma_bar(gamma_bar)

  if (!identical(engine, "sdp")) {
    stop("Algorithm 2 requires `engine = 'sdp'`.", call. = FALSE)
  }

  if (!is.null(xi_range)) {
    if (!is.numeric(xi_range) || length(xi_range) != 2L || anyNA(xi_range) || any(!is.finite(xi_range)) || any(xi_range <= 0)) {
      stop("`xi_range` must be a length-2 positive numeric vector.", call. = FALSE)
    }
    xi_lower <- min(xi_range)
    xi_upper <- max(xi_range)
  } else if (!is.null(calibration_range)) {
    calibration_range <- validate_positive_grid(calibration_range, "calibration_range")
    if (length(calibration_range) != 2L) {
      stop("`calibration_range` must be a length-2 positive numeric vector.", call. = FALSE)
    }
    xi_values <- 1 / calibration_range
    xi_lower <- min(xi_values)
    xi_upper <- max(xi_values)
  }

  xi_lower <- validate_positive_xi(xi_lower, "xi_lower")
  xi_upper <- validate_positive_xi(xi_upper, "xi_upper")
  if (xi_lower > xi_upper) {
    tmp <- xi_lower; xi_lower <- xi_upper; xi_upper <- tmp
  }

  n <- nrow(W)
  k_grid <- validate_k_range(min_k, max_k, n)

  lower_sdp <- sdp_quadratic_relaxation(W = W, xi = xi_lower, k_constraint = FALSE)
  upper_sdp <- sdp_quadratic_relaxation(W = W, xi = xi_upper, k_constraint = FALSE)

  regret_sdp <- sdp_regret_relaxation(
    W = W,
    xi_lower = xi_lower,
    xi_upper = xi_upper,
    lower_bound_lower = lower_sdp$objective,
    lower_bound_upper = upper_sdp$objective
  )

  X_hat <- regret_sdp$X_matrix
  solutions <- vector("list", length(k_grid))
  rho_values <- rep(Inf, length(k_grid))

  for (i in seq_along(k_grid)) {
    k <- k_grid[i]
    eig_k <- if (is.null(n_eig)) k else n_eig

    candidates <- candidate_clusterings_from_matrix(
      W = W,
      my_matrix = X_hat,
      num_clusters = k,
      n_eig = eig_k,
      engine = "sdp",
      methods = methods,
      seed = seed,
      try_sign_flip = try_sign_flip,
      k_constraint = k_constraint,
      gamma_bar = gamma_bar
    )

    if (length(candidates) == 0L) {
      solutions[[i]] <- list(num_clusters = k, error = "No valid rounded candidates", rho = Inf)
      next
    }

    candidate_summary <- data.frame(
      method = names(candidates),
      objective_lower = NA_real_,
      objective_upper = NA_real_,
      rho = NA_real_,
      stringsAsFactors = FALSE
    )

    for (j in seq_along(candidates)) {
      obj_lower <- compute_causal_clustering_objective(W, xi_lower, candidates[[j]], objective_type = "squared")
      obj_upper <- compute_causal_clustering_objective(W, xi_upper, candidates[[j]], objective_type = "squared")
      rho_j <- max(
        safe_ratio(obj_lower, lower_sdp$objective),
        safe_ratio(obj_upper, upper_sdp$objective)
      )
      candidate_summary$objective_lower[j] <- obj_lower
      candidate_summary$objective_upper[j] <- obj_upper
      candidate_summary$rho[j] <- rho_j
    }

    best_j <- which.min(candidate_summary$rho)
    rho_values[i] <- candidate_summary$rho[best_j]
    solutions[[i]] <- list(
      clusters = candidates[[best_j]],
      num_clusters = k,
      method = candidate_summary$method[best_j],
      objective_lower = candidate_summary$objective_lower[best_j],
      objective_upper = candidate_summary$objective_upper[best_j],
      rho = candidate_summary$rho[best_j],
      candidates = candidates,
      candidate_summary = candidate_summary
    )
  }

  names(rho_values) <- paste0("k_", k_grid)
  all_rho <- rho_values
  bernoulli_clusters <- NULL

  if (isTRUE(include_bernoulli)) {
    bernoulli_clusters <- seq_len(n)
    bernoulli_lower <- compute_causal_clustering_objective(W, xi_lower, bernoulli_clusters, objective_type = "squared")
    bernoulli_upper <- compute_causal_clustering_objective(W, xi_upper, bernoulli_clusters, objective_type = "squared")
    all_rho <- c(all_rho, bernoulli = max(
      safe_ratio(bernoulli_lower, lower_sdp$objective),
      safe_ratio(bernoulli_upper, upper_sdp$objective)
    ))
  }

  if (all(!is.finite(all_rho))) {
    stop("No finite endpoint-regret objective was obtained for any rounded candidate.", call. = FALSE)
  }

  best_index <- which.min(all_rho)
  best_clusters <- if (best_index <= length(k_grid)) solutions[[best_index]]$clusters else bernoulli_clusters
  selected_solution <- if (best_index <= length(k_grid)) solutions[[best_index]] else NULL
  selected_method <- if (!is.null(selected_solution) && !is.null(selected_solution$method)) {
    selected_solution$method
  } else if (best_index > length(k_grid)) {
    "bernoulli"
  } else {
    NA_character_
  }
  selected_candidate_summary <- if (!is.null(selected_solution) &&
      !is.null(selected_solution$candidate_summary)) {
    selected_solution$candidate_summary
  } else {
    NULL
  }

  list(
    clusters = best_clusters,
    selected_k = if (best_index <= length(k_grid)) k_grid[best_index] else n,
    selected_index = best_index,
    selected_method = selected_method,
    selected_candidate_summary = selected_candidate_summary,
    objective = unname(all_rho[best_index]),
    rho_values = all_rho,
    xi_lower = xi_lower,
    xi_upper = xi_upper,
    xi_range = c(xi_lower, xi_upper),
    lower_bound_lower = lower_sdp$objective,
    lower_bound_upper = upper_sdp$objective,
    endpoint_lower_bounds = c(lower = lower_sdp$objective, upper = upper_sdp$objective),
    components = clustering_objective_components(W, best_clusters),
    solutions = solutions,
    sdp_regret_objective = regret_sdp$objective,
    regret_sdp = if (isTRUE(keep_sdp_solutions)) regret_sdp else NULL,
    endpoint_sdp = if (isTRUE(keep_sdp_solutions)) list(lower = lower_sdp, upper = upper_sdp) else NULL,
    methods = methods,
    include_bernoulli = include_bernoulli,
    k_constraint = isTRUE(k_constraint),
    min_k = min(k_grid),
    max_k = max(k_grid),
    gamma_bar = gamma_bar,
    engine = engine,
    objective_type = "squared"
  )
}


# Unified public interface --------------------------------------------------

#' Causal clustering algorithm
#'
#' @description
#' Unified interface for the causal clustering algorithms. Supply exactly one
#' calibration input. A scalar `xi` or scalar `calibration` calls Algorithm 1
#' for a fixed calibration. A grid or range of `xi` or `calibration` values
#' calls Algorithm 2 and uses only the endpoint range for the endpoint-regret
#' criterion.
#'
#' @param W Symmetric binary adjacency matrix.
#' @param xi Optional scalar or grid of `xi = (lambda * phibar_n^2 / psibar)^(-1)` values.
#' @param calibration Optional scalar or grid of `lambda * phibar_n^2 / psibar` values.
#' @param xi_grid Optional positive grid of `xi` values.
#' @param xi_range Optional length-two positive range of `xi` values.
#' @param calibration_grid Optional positive grid of `lambda * phibar_n^2 / psibar` values.
#' @param calibration_range Optional length-two positive range of `lambda * phibar_n^2 / psibar` values.
#' @param min_k Minimum number of clusters to consider.
#' @param max_k Maximum number of clusters to consider.
#' @param n_eig Optional number of eigenvectors used for rounding. If `NULL`,
#'   the number of eigenvectors equals the candidate K.
#' @param seed Optional seed used by k-means rounding.
#' @param engine Optimization engine. Use `"sdp"` for the semidefinite relaxation and
#'   `"spectral"` for a spectral heuristic. Calibration grids/ranges require `"sdp"`.
#' @param methods Character vector of rounding/discretization methods. See
#'   [available_discretization_methods()]. All requested methods are run and
#'   the rounded candidate with the smallest realized objective is selected.
#' @param include_bernoulli Logical; if `TRUE`, also compare against the
#'   all-singleton Bernoulli design.
#' @param objective_type Must be `"squared"`.
#' @param k_constraint Logical Boolean from Algorithm 1. For Algorithm 2 it
#'   filters rounded candidates by the cluster-size constraint.
#' @param gamma_bar Upper cluster-size multiplier used by optional constraints or filtering.
#' @param box_constraints Logical; whether to impose optional `0 <= X_ij <= 1`
#'   constraints when Algorithm 1 is run with `k_constraint = TRUE`.
#' @param try_sign_flip Whether to try rounding both X and -X.
#' @param keep_search_matrices Logical; if `TRUE`, keep the relaxation matrix in
#'   Algorithm 1 results.
#' @param keep_sdp_solutions Logical; if `TRUE`, keep raw SDP solution objects.
#'
#' @return A list. For scalar calibrations it contains the Algorithm 1 output,
#'   including `clusters`, `selected_k`, `selected_method`, `objective`,
#'   `components`, SDP lower bounds, and `Gamma_n` when available. For grids or
#'   ranges it contains the Algorithm 2 output, including `clusters`,
#'   `selected_k`, `selected_method`, endpoint lower bounds, endpoint regret
#'   values, and objective components.
#' @export
causal_clustering_algorithm <- function(W,
                                        xi = NULL,
                                        calibration = NULL,
                                        xi_grid = NULL,
                                        xi_range = NULL,
                                        calibration_grid = NULL,
                                        calibration_range = NULL,
                                        min_k = 2,
                                        max_k = floor(nrow(W) / 2),
                                        n_eig = NULL,
                                        seed = NULL,
                                        engine = c("sdp", "spectral"),
                                        methods = available_discretization_methods(),
                                        include_bernoulli = FALSE,
                                        objective_type = "squared",
                                        k_constraint = FALSE,
                                        gamma_bar = 10,
                                        box_constraints = TRUE,
                                        try_sign_flip = TRUE,
                                        keep_search_matrices = FALSE,
                                        keep_sdp_solutions = FALSE) {
  W <- validate_adjacency_matrix(W, binary = TRUE)
  engine <- match.arg(engine)
  methods <- standardize_discretization_methods(methods)
  objective_type <- validate_objective_type(objective_type)

  supplied <- c(
    xi = !is.null(xi),
    calibration = !is.null(calibration),
    xi_grid = !is.null(xi_grid),
    xi_range = !is.null(xi_range),
    calibration_grid = !is.null(calibration_grid),
    calibration_range = !is.null(calibration_range)
  )

  if (sum(supplied) != 1L) {
    stop(
      "Supply exactly one of `xi`, `calibration`, `xi_grid`, `xi_range`, ",
      "`calibration_grid`, or `calibration_range`.",
      call. = FALSE
    )
  }

  fixed_call <- function(xi_value = NULL, calibration_value = NULL) {
    out <- causal_clustering_algorithm1(
      W = W,
      xi = xi_value,
      calibration = calibration_value,
      min_k = min_k,
      max_k = max_k,
      n_eig = n_eig,
      seed = seed,
      engine = engine,
      methods = methods,
      include_bernoulli = include_bernoulli,
      objective_type = objective_type,
      k_constraint = k_constraint,
      gamma_bar = gamma_bar,
      box_constraints = box_constraints,
      try_sign_flip = try_sign_flip,
      keep_search_matrices = keep_search_matrices,
      keep_sdp_solutions = keep_sdp_solutions
    )
    out$algorithm <- "algorithm1_fixed_calibration"
    out
  }

  range_call <- function(xi_range_value = NULL, calibration_range_value = NULL) {
    if (!identical(engine, "sdp")) {
      stop("Calibration grids/ranges use Algorithm 2 and require `engine = 'sdp'`.", call. = FALSE)
    }
    out <- causal_clustering_algorithm2(
      W = W,
      xi_range = xi_range_value,
      calibration_range = calibration_range_value,
      min_k = min_k,
      max_k = max_k,
      n_eig = n_eig,
      seed = seed,
      engine = engine,
      methods = methods,
      include_bernoulli = include_bernoulli,
      k_constraint = k_constraint,
      gamma_bar = gamma_bar,
      try_sign_flip = try_sign_flip,
      keep_sdp_solutions = keep_sdp_solutions
    )
    out$algorithm <- "algorithm2_endpoint_regret"
    out
  }

  if (!is.null(xi)) {
    xi <- as.numeric(xi)
    if (length(xi) == 1L) return(fixed_call(xi_value = xi))
    xi <- validate_positive_grid(xi, "xi")
    return(range_call(xi_range_value = range(xi)))
  }

  if (!is.null(calibration)) {
    calibration <- as.numeric(calibration)
    if (length(calibration) == 1L) return(fixed_call(calibration_value = calibration))
    calibration <- validate_positive_grid(calibration, "calibration")
    return(range_call(calibration_range_value = range(calibration)))
  }

  if (!is.null(xi_grid)) {
    xi_grid <- validate_positive_grid(xi_grid, "xi_grid")
    if (length(xi_grid) == 1L) return(fixed_call(xi_value = xi_grid))
    return(range_call(xi_range_value = range(xi_grid)))
  }

  if (!is.null(xi_range)) {
    xi_range <- validate_positive_grid(xi_range, "xi_range")
    if (length(xi_range) != 2L) stop("`xi_range` must be length two.", call. = FALSE)
    return(range_call(xi_range_value = range(xi_range)))
  }

  if (!is.null(calibration_grid)) {
    calibration_grid <- validate_positive_grid(calibration_grid, "calibration_grid")
    if (length(calibration_grid) == 1L) return(fixed_call(calibration_value = calibration_grid))
    return(range_call(calibration_range_value = range(calibration_grid)))
  }

  calibration_range <- validate_positive_grid(calibration_range, "calibration_range")
  if (length(calibration_range) != 2L) {
    stop("`calibration_range` must be length two.", call. = FALSE)
  }
  range_call(calibration_range_value = range(calibration_range))
}

#' Adaptive causal clustering
#'
#' @description
#' Public wrapper for Algorithm 2. Supply either `xi_range`,
#' `xi_grid`, `calibration_range`, or `calibration_grid`; the function converts
#' the input to an endpoint range and calls [causal_clustering_algorithm2()].
#'
#' @param W Symmetric binary adjacency matrix.
#' @param calibration_grid Optional positive grid for
#'   `lambda * phibar_n^2 / psibar`. The endpoint range of `1 / calibration_grid`
#'   is used.
#' @param xi_grid Optional positive grid of `xi` values. The endpoint range
#'   is used.
#' @param xi_range Optional length-two positive range for `xi`.
#' @param calibration_range Optional length-two positive range for
#'   `lambda * phibar_n^2 / psibar`; converted to a `xi` range.
#' @param min_k Minimum number of clusters to consider.
#' @param max_k Maximum number of clusters to consider.
#' @param n_eig Optional number of eigenvectors used for rounding.
#' @param seed Optional seed used by k-means rounding.
#' @param engine Must be `"sdp"` for Algorithm 2.
#' @param methods Character vector of rounding/discretization methods.
#' @param include_bernoulli Logical; if `TRUE`, also compare the all-singleton
#'   Bernoulli design.
#' @param objective_type Must be `"squared"` for Algorithm 2.
#' @param k_constraint Logical; if `TRUE`, rounded candidates are filtered by
#'   the cluster-size constraint.
#' @param gamma_bar Upper cluster-size multiplier used when `k_constraint = TRUE`.
#' @param try_sign_flip Logical; if `TRUE`, round both the relaxation matrix and
#'   its negative.
#' @param keep_sdp_solutions Logical; if `TRUE`, keep raw SDP solution objects.
#'
#' @return The list returned by [causal_clustering_algorithm2()].
#' @export
adaptive_causal_clustering <- function(W,
                                       calibration_grid = NULL,
                                       xi_grid = NULL,
                                       xi_range = NULL,
                                       calibration_range = NULL,
                                       min_k = 2,
                                       max_k = floor(nrow(W) / 2),
                                       n_eig = NULL,
                                       seed = NULL,
                                       engine = "sdp",
                                       methods = available_discretization_methods(),
                                       include_bernoulli = FALSE,
                                       objective_type = "squared",
                                       k_constraint = FALSE,
                                       gamma_bar = 10,
                                       try_sign_flip = TRUE,
                                       keep_sdp_solutions = FALSE) {
  engine <- match.arg(engine)
  objective_type <- validate_objective_type(objective_type)

  if (is.null(xi_range)) {
    if (!is.null(xi_grid)) {
      xi_grid <- validate_positive_grid(xi_grid, "xi_grid")
      xi_range <- range(xi_grid)
    } else if (!is.null(calibration_range)) {
      calibration_range <- validate_positive_grid(calibration_range, "calibration_range")
      xi_range <- range(1 / calibration_range)
    } else if (!is.null(calibration_grid)) {
      calibration_grid <- validate_positive_grid(calibration_grid, "calibration_grid")
      xi_range <- range(1 / calibration_grid)
    } else {
      stop("Supply `xi_range`, `xi_grid`, `calibration_range`, or `calibration_grid`.", call. = FALSE)
    }
  }

  causal_clustering_algorithm2(
    W = W,
    xi_range = xi_range,
    min_k = min_k,
    max_k = max_k,
    n_eig = n_eig,
    seed = seed,
    engine = engine,
    methods = methods,
    include_bernoulli = include_bernoulli,
    k_constraint = k_constraint,
    gamma_bar = gamma_bar,
    try_sign_flip = try_sign_flip,
    keep_sdp_solutions = keep_sdp_solutions
  )
}

# Additional objective wrappers --------------------------------------------

#' Compute the frontier objective from Equation (9)
#'
#' @param W Symmetric binary adjacency matrix.
#' @param xi Non-negative tuning parameter.
#' @param clusters Cluster labels with length `nrow(W)`.
#' @param objective_type Must be `"squared"` for Equation (9).
#'
#' @return A numeric scalar objective value.
#' @export
compute_frontier_objective <- function(W,
                                       xi,
                                       clusters,
                                       objective_type = "squared") {
  compute_causal_clustering_objective(W, xi, clusters, objective_type = objective_type)
}

#' Compute the mean-squared-error objective from Equation (8)
#'
#' @description
#' Computes `variance + calibration * bias^2`. When `calibration = 1 / xi`,
#' this scalarization has the same minimizers as Equation (9).
#'
#' @param W Symmetric binary adjacency matrix.
#' @param calibration Non-negative value of `lambda * phibar_n^2 / psibar`.
#' @param clusters Cluster labels with length `nrow(W)`.
#'
#' @return A numeric scalar objective value.
#' @export
compute_mse_clustering_objective <- function(W, calibration, clusters) {
  if (!is.numeric(calibration) || length(calibration) != 1L ||
      is.na(calibration) || !is.finite(calibration) || calibration < 0) {
    stop("`calibration` must be a single finite non-negative number.", call. = FALSE)
  }
  components <- clustering_objective_components(W, clusters)
  components$variance + calibration * components$bias^2
}
