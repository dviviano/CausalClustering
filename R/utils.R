# Internal helpers ---------------------------------------------------------

validate_adjacency_matrix <- function(W, binary = TRUE) {
  if (!is.matrix(W)) {
    stop("`W` must be a matrix.", call. = FALSE)
  }
  if (nrow(W) != ncol(W)) {
    stop("`W` must be square.", call. = FALSE)
  }
  if (!is.numeric(W)) {
    stop("`W` must be numeric.", call. = FALSE)
  }
  if (!is.logical(binary) || length(binary) != 1L || is.na(binary)) {
    stop("`binary` must be TRUE or FALSE.", call. = FALSE)
  }
  if (anyNA(W) || any(!is.finite(W))) {
    stop("`W` must contain finite numeric entries.", call. = FALSE)
  }
  W <- 0.5 * (W + t(W))
  if (isTRUE(binary)) {
    W <- 1 * (W != 0)
  }
  diag(W) <- 0
  storage.mode(W) <- "double"
  W
}

safe_difference_in_means <- function(y, D) {
  treated_n <- sum(D == 1)
  control_n <- sum(D == 0)

  if (treated_n == 0 || control_n == 0) {
    return(NA_real_)
  }

  mean(y[D == 1]) - mean(y[D == 0])
}

sample_sparse_sign_vector <- function(p, strong = FALSE) {
  out <- rep(0, p)
  if (p < 2) {
    out[] <- sample(c(-1, 1), size = p, replace = TRUE)
    return(out)
  }

  idx <- sample.int(p, size = 2, replace = FALSE)
  out[idx[1]] <- if (strong) 1.5 else 1
  out[idx[2]] <- if (strong) -1.5 else -1
  out
}

draw_model_coefficients <- function(p, heterogeneity) {
  mu <- sample(c(-1, 1), size = 1)

  if (!isTRUE(heterogeneity)) {
    return(list(
      mu = mu,
      b1 = rep(0, p),
      b2 = rep(0, p),
      b3 = rep(0, p),
      b4 = rep(0, p),
      b5 = rep(0, p)
    ))
  }

  if (p <= 10) {
    return(list(
      mu = mu,
      b1 = sample(c(-1, 1), size = p, replace = TRUE),
      b2 = sample(c(-1, 1), size = p, replace = TRUE),
      b3 = sample(c(-1, 1), size = p, replace = TRUE),
      b4 = sample(c(-1.5, 1.5), size = p, replace = TRUE),
      b5 = sample(c(-1, 1), size = p, replace = TRUE)
    ))
  }

  list(
    mu = mu,
    b1 = sample_sparse_sign_vector(p),
    b2 = sample_sparse_sign_vector(p),
    b3 = sample_sparse_sign_vector(p),
    b4 = sample_sparse_sign_vector(p, strong = TRUE),
    b5 = sample_sparse_sign_vector(p)
  )
}

build_geometric_graph <- function(X) {
  n <- nrow(X)
  if (ncol(X) < 2) {
    stop("Geometric graph generation requires at least 2 covariates.")
  }

  radius <- sqrt(4 / (2.75 * n))
  W <- matrix(0, nrow = n, ncol = n)

  for (i in seq_len(n - 1L)) {
    d <- abs(X[i, 1] - X[(i + 1L):n, 1]) / 2 +
      abs(X[i, 2] - X[(i + 1L):n, 2]) / 2
    W[(i + 1L):n, i] <- as.integer(d <= radius)
  }

  W <- W + t(W)
  diag(W) <- 0
  W
}

build_erdos_renyi_graph <- function(n, neighb) {
  prob <- min(max(neighb / max(n - 1, 1), 0), 1)
  W <- matrix(0, nrow = n, ncol = n)
  idx <- upper.tri(W)
  W[idx] <- stats::rbinom(sum(idx), size = 1, prob = prob)
  W <- W + t(W)
  diag(W) <- 0
  W
}

build_barabasi_graph <- function(n, neighb) {
  if (n < 3) {
    W <- matrix(0, nrow = n, ncol = n)
    diag(W) <- 0
    return(W)
  }

  m0 <- max(2L, min(n, floor(n / 5)))
  m_attach <- max(1L, min(m0, round(neighb)))

  W <- build_erdos_renyi_graph(m0, neighb)
  if (sum(W) == 0) {
    W[1, 2] <- 1
    W[2, 1] <- 1
  }

  for (new_node in seq.int(m0 + 1L, n)) {
    old_n <- nrow(W)
    deg <- rowSums(W)
    prob <- (deg + 1) / sum(deg + 1)
    targets <- sample.int(old_n, size = min(m_attach, old_n), replace = FALSE, prob = prob)

    W_new <- matrix(0, nrow = old_n + 1L, ncol = old_n + 1L)
    W_new[seq_len(old_n), seq_len(old_n)] <- W
    W_new[new_node, targets] <- 1
    W_new[targets, new_node] <- 1
    W <- W_new
  }

  diag(W) <- 0
  W
}

mean_result_frames <- function(result_list) {
  keep <- Filter(Negate(is.null), result_list)
  keep <- Filter(function(x) !inherits(x, "try-error") && !all(is.na(x)), keep)

  if (length(keep) == 0) {
    return(NULL)
  }

  numeric_cols <- vapply(keep[[1]], is.numeric, logical(1))
  total <- keep[[1]][, numeric_cols, drop = FALSE]

  if (length(keep) > 1) {
    for (i in 2:length(keep)) {
      total <- total + keep[[i]][, numeric_cols, drop = FALSE]
    }
  }

  total / length(keep)
}

check_sdp_dependencies <- function() {
  pkgs <- c("Matrix", "sdpt3r")
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  
  if (length(missing) > 0) {
    stop(
      "The SDP engine requires these packages: ",
      paste(missing, collapse = ", "),
      ". Please install them before using engine = 'sdp'.",
      call. = FALSE
    )
  }
  
  invisible(TRUE)
}