#' Run one simulated network instance
#'
#' @param seed Random seed.
#' @param parameters_graph Graph-generation parameters.
#' @param W Optional adjacency matrix. If `NULL`, a network is simulated.
#' @param xi_seq Sequence of xi values.
#' @param min_k Minimum number of clusters in the causal search.
#' @param max_k Maximum number of clusters in the causal search.
#' @param include_louvain Logical; whether to compute the Louvain baseline.
#' @param engine Search engine for the causal clustering step. Use `"sdp"` for
#'   the semidefinite relaxation and `"spectral"` for a spectral heuristic.
#' @param methods Discretization methods used after the relaxation step. All requested methods are tried and the lowest realized objective is selected.
#' @param objective_type Must be `"squared"`.
#' @param k_constraint Logical; passed to `search_causal_clustering()`.
#' @param gamma_bar Upper cluster-size multiplier used when `k_constraint = TRUE`.
#' @param box_constraints Logical; whether to include optional box constraints
#'   when `k_constraint = TRUE`.
#'
#' @return A list with `results`, `W`, and `clusters`.
#' @export
run_single_network <- function(seed,
                               parameters_graph = list(type_graph = "geometric", n = 100, neighb = 2),
                               W = NULL,
                               xi_seq = seq(from = 0.1, to = 10, by = 1),
                               min_k = 2,
                               max_k = NULL,
                               include_louvain = FALSE,
                               engine = c("sdp", "spectral"),
                               methods = available_discretization_methods(),
                               objective_type = "squared",
                               k_constraint = FALSE,
                               gamma_bar = 10,
                               box_constraints = TRUE) {
  engine <- match.arg(engine)
  methods <- standardize_discretization_methods(methods)
  objective_type <- validate_objective_type(objective_type)

  if (!is.null(seed)) {
    set.seed(seed)
  }

  if (is.null(W)) {
    simulated_data <- simulate_network_data(parameters_graph = parameters_graph, seed = seed)
    W <- simulated_data$W
  }
  W <- validate_adjacency_matrix(W, binary = TRUE)

  n <- nrow(W)
  if (is.null(max_k)) {
    max_k <- max(min_k, floor(n / 2))
  }

  epsilon_clusters <- cluster_epsilon_net(W)$clusters
  spectral_clusters <- cluster_spectral(
    W,
    num_clusters = max(2L, floor(n / 3)),
    normalized = TRUE,
    n_eig = max(2L, floor(n / 3)),
    hierarchical = FALSE,
    seed = seed
  )

  louvain_clusters <- NULL
  if (isTRUE(include_louvain)) {
    louvain_clusters <- cluster_louvain_membership(W)
  }

  causal_solutions <- vector("list", length(xi_seq))
  objective_minimax <- numeric(length(xi_seq))
  objective_epsilon <- numeric(length(xi_seq))
  objective_spectral <- numeric(length(xi_seq))
  objective_louvain <- rep(NA_real_, length(xi_seq))
  num_clusters_causal <- integer(length(xi_seq))

  for (i in seq_along(xi_seq)) {
    xi <- xi_seq[i]
    causal_solutions[[i]] <- search_causal_clustering(
      W = W,
      xi = xi,
      min_k = min_k,
      max_k = max_k,
      seed = seed,
      engine = engine,
      methods = methods,
      objective_type = objective_type,
      k_constraint = k_constraint,
      gamma_bar = gamma_bar,
      box_constraints = box_constraints
    )

    objective_minimax[i] <- compute_causal_clustering_objective(
      W, xi, causal_solutions[[i]]$clusters, objective_type = objective_type
    )
    objective_epsilon[i] <- compute_causal_clustering_objective(
      W, xi, epsilon_clusters, objective_type = objective_type
    )
    objective_spectral[i] <- compute_causal_clustering_objective(
      W, xi, spectral_clusters, objective_type = objective_type
    )
    num_clusters_causal[i] <- length(unique(causal_solutions[[i]]$clusters))

    if (!is.null(louvain_clusters)) {
      objective_louvain[i] <- compute_causal_clustering_objective(
        W, xi, louvain_clusters, objective_type = objective_type
      )
    }
  }

  data_results <- data.frame(
    xi_seq = xi_seq,
    num_clusters_causal = num_clusters_causal,
    objective_minimax = objective_minimax,
    objective_epsilon_net = objective_epsilon,
    objective_spectral = objective_spectral,
    objective_louvain = objective_louvain,
    objective_type = objective_type,
    stringsAsFactors = FALSE
  )

  list(
    results = data_results,
    W = W,
    clusters = list(
      epsilon_net = epsilon_clusters,
      louvain = louvain_clusters,
      spectral = spectral_clusters,
      causal = causal_solutions
    ),
    engine = engine,
    methods = methods,
    objective_type = objective_type,
    k_constraint = k_constraint,
    gamma_bar = gamma_bar,
    box_constraints = box_constraints
  )
}

#' Simulate a grid of networks
#'
#' @param params_sim_networks A list with entries `n`, `type_graph`, `n_sim`,
#'   `xi_seq`, and optionally `neighb`, `min_k`, `max_k`.
#' @param engine Search engine for the causal clustering step.
#' @param methods Discretization methods used after the relaxation step. All requested methods are tried and the lowest realized objective is selected.
#' @param objective_type Must be `"squared"`.
#' @param k_constraint Logical; passed to `run_single_network()`.
#' @param gamma_bar Upper cluster-size multiplier used when `k_constraint = TRUE`.
#' @param box_constraints Logical; whether to include optional box constraints.
#'
#' @return A nested list of simulation results.
#' @export
simulate_network_grid <- function(params_sim_networks,
                                  engine = c("sdp", "spectral"),
                                  methods = available_discretization_methods(),
                                  objective_type = "squared",
                                  k_constraint = FALSE,
                                  gamma_bar = 10,
                                  box_constraints = TRUE) {
  engine <- match.arg(engine)
  methods <- standardize_discretization_methods(methods)
  objective_type <- validate_objective_type(objective_type)

  n_values <- params_sim_networks$n
  graph_types <- params_sim_networks$type_graph
  xi_seq <- params_sim_networks$xi_seq
  n_sim <- params_sim_networks$n_sim
  neighb <- if (!is.null(params_sim_networks$neighb)) params_sim_networks$neighb else 2
  min_k <- if (!is.null(params_sim_networks$min_k)) params_sim_networks$min_k else 2
  max_k <- params_sim_networks$max_k

  grid <- expand.grid(sample_size = n_values, network = graph_types, stringsAsFactors = FALSE)
  simulations <- vector("list", nrow(grid))

  for (j in seq_len(nrow(grid))) {
    parameters_graph <- list(
      n = grid$sample_size[j],
      type_graph = grid$network[j],
      neighb = neighb
    )

    simulations[[j]] <- lapply(seq_len(n_sim), function(i) {
      run_single_network(
        seed = i,
        parameters_graph = parameters_graph,
        xi_seq = xi_seq,
        min_k = min_k,
        max_k = if (is.null(max_k)) floor(parameters_graph$n / 2) else max_k,
        engine = engine,
        methods = methods,
        objective_type = objective_type,
        k_constraint = k_constraint,
        gamma_bar = gamma_bar,
        box_constraints = box_constraints
      )
    })
  }

  simulations
}
