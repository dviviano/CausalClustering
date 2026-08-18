# Plotting helpers ----------------------------------------------------------

check_ggplot2 <- function() {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plotting functions.", call. = FALSE)
  }
  invisible(TRUE)
}

#' Plot objective and cluster-count paths for one network run
#'
#' @param data_results Data frame produced by `run_single_network()$results`.
#'
#' @return A named list with two ggplot objects: `plot_objective` and
#'   `plot_num_clusters`.
#' @export
plot_objective_path <- function(data_results) {
  check_ggplot2()

  n_rows <- nrow(data_results)
  long_df <- data.frame(
    xi_seq = rep(data_results$xi_seq, 4),
    objective = c(
      data_results$objective_minimax,
      data_results$objective_epsilon_net,
      data_results$objective_spectral,
      data_results$objective_louvain
    ),
    clustering_type = c(
      rep("Causal clustering", n_rows),
      rep("Epsilon net", n_rows),
      rep("Spectral clustering", n_rows),
      rep("Louvain", n_rows)
    ),
    stringsAsFactors = FALSE
  )

  p_objective <- ggplot2::ggplot(
    long_df,
    ggplot2::aes(
      x = .data$xi_seq,
      y = .data$objective,
      colour = .data$clustering_type,
      group = .data$clustering_type
    )
  ) +
    ggplot2::geom_line(na.rm = TRUE) +
    ggplot2::geom_point(na.rm = TRUE) +
    ggplot2::theme_bw()

  p_clusters <- ggplot2::ggplot(
    data_results,
    ggplot2::aes(x = .data$xi_seq, y = .data$num_clusters_causal)
  ) +
    ggplot2::geom_line(na.rm = TRUE) +
    ggplot2::geom_point(na.rm = TRUE) +
    ggplot2::theme_bw()

  list(plot_objective = p_objective, plot_num_clusters = p_clusters)
}

#' Plot MSE, objective, bias, variance, and cluster-count paths
#'
#' @param data_results Data frame containing the columns produced by the MSE
#'   simulation workflow.
#'
#' @return A named list of ggplot objects.
#' @export
plot_mse_path <- function(data_results) {
  check_ggplot2()

  n_rows <- nrow(data_results)
  long_df <- data.frame(
    xi_seq = rep(data_results$xi_seq, 4),
    objective = c(
      data_results$objective_minimax,
      data_results$objective_epsilon_net,
      data_results$objective_spectral,
      data_results$objective_louvain
    ),
    weighted_mse = c(
      data_results$weighted_mse_causal,
      data_results$weighted_mse_epsilon_net,
      data_results$weighted_mse_spectral,
      data_results$weighted_mse_louvain
    ),
    bias = c(
      data_results$bias_causal,
      data_results$bias_epsilon_net,
      data_results$bias_spectral,
      data_results$bias_louvain
    ),
    variance = c(
      data_results$variance_causal,
      data_results$variance_epsilon_net,
      data_results$variance_spectral,
      data_results$variance_louvain
    ),
    clustering_type = c(
      rep("Causal clustering", n_rows),
      rep("Epsilon net", n_rows),
      rep("Spectral clustering", n_rows),
      rep("Louvain", n_rows)
    ),
    stringsAsFactors = FALSE
  )

  base_aes <- ggplot2::aes(
    x = .data$xi_seq,
    colour = .data$clustering_type,
    group = .data$clustering_type
  )
  base_theme <- ggplot2::theme_bw()

  plot_mse <- ggplot2::ggplot(long_df, base_aes) +
    ggplot2::aes(y = .data$weighted_mse) +
    ggplot2::geom_line(na.rm = TRUE) + ggplot2::geom_point(na.rm = TRUE) + base_theme

  plot_objective <- ggplot2::ggplot(long_df, base_aes) +
    ggplot2::aes(y = .data$objective) +
    ggplot2::geom_line(na.rm = TRUE) + ggplot2::geom_point(na.rm = TRUE) + base_theme

  plot_bias <- ggplot2::ggplot(long_df, base_aes) +
    ggplot2::aes(y = .data$bias) +
    ggplot2::geom_line(na.rm = TRUE) + ggplot2::geom_point(na.rm = TRUE) + base_theme

  plot_variance <- ggplot2::ggplot(long_df, base_aes) +
    ggplot2::aes(y = .data$variance) +
    ggplot2::geom_line(na.rm = TRUE) + ggplot2::geom_point(na.rm = TRUE) + base_theme

  plot_num_clusters <- ggplot2::ggplot(
    data_results,
    ggplot2::aes(x = .data$xi_seq, y = .data$num_clusters_causal)
  ) + ggplot2::geom_line(na.rm = TRUE) + ggplot2::geom_point(na.rm = TRUE) + base_theme

  list(
    plot_mse = plot_mse,
    plot_objective = plot_objective,
    plot_variance = plot_variance,
    plot_num_clusters = plot_num_clusters,
    plot_bias = plot_bias
  )
}

#' Plot number of clusters by sample size
#'
#' @param extracted_results Data frame produced by an extraction workflow over
#'   repeated network simulations.
#'
#' @return A list of ggplot objects, one per network type.
#' @export
plot_cluster_count_by_sample_size <- function(extracted_results) {
  check_ggplot2()

  nets <- unique(extracted_results$network)
  extracted_results$sample_size <- as.character(extracted_results$sample_size)

  lapply(nets, function(net_name) {
    ggplot2::ggplot(
      extracted_results[extracted_results$network == net_name, , drop = FALSE],
      ggplot2::aes(
        x = .data$xi_seq,
        y = .data$num_clusters_causal,
        colour = .data$sample_size,
        group = .data$sample_size
      )
    ) +
      ggplot2::geom_line(na.rm = TRUE) +
      ggplot2::geom_point(na.rm = TRUE) +
      ggplot2::ggtitle(net_name) +
      ggplot2::xlab("xi") +
      ggplot2::ylab("Number of clusters") +
      ggplot2::theme_bw()
  })
}

#' Plot average objective curves across sample sizes and network types
#'
#' @param extracted_results Data frame produced by an objective extraction
#'   workflow over repeated simulations.
#'
#' @return A list of ggplot objects.
#' @export
plot_objective_grid <- function(extracted_results) {
  check_ggplot2()

  n_rows <- nrow(extracted_results)
  final_dataset <- data.frame(
    log_objective = log(c(
      extracted_results$objective_minimax,
      extracted_results$objective_epsilon_net,
      extracted_results$objective_spectral,
      extracted_results$objective_louvain
    )),
    xi_seq = rep(extracted_results$xi_seq, 4),
    sample_size = rep(extracted_results$sample_size, 4),
    network = rep(as.character(extracted_results$network), 4),
    clustering_type = c(
      rep("Causal clustering", n_rows),
      rep("Epsilon net", n_rows),
      rep("Spectral", n_rows),
      rep("Louvain", n_rows)
    ),
    stringsAsFactors = FALSE
  )

  unique_n <- unique(final_dataset$sample_size)
  unique_network <- unique(final_dataset$network)
  plots <- list()
  acc <- 1L

  for (n_val in unique_n) {
    for (net_name in unique_network) {
      plot_data <- final_dataset[
        final_dataset$sample_size == n_val & final_dataset$network == net_name,
        ,
        drop = FALSE
      ]
      plots[[acc]] <- ggplot2::ggplot(
        plot_data,
        ggplot2::aes(
          x = .data$xi_seq,
          y = .data$log_objective,
          colour = .data$clustering_type,
          group = .data$clustering_type
        )
      ) +
        ggplot2::geom_line(na.rm = TRUE) +
        ggplot2::geom_point(na.rm = TRUE) +
        ggplot2::ggtitle(paste0("N = ", n_val, ", ", net_name)) +
        ggplot2::xlab("xi") +
        ggplot2::ylab("log objective") +
        ggplot2::theme_bw()
      acc <- acc + 1L
    }
  }

  plots
}

#' Plot average weighted-MSE, bias, and variance curves
#'
#' @param extracted_results Data frame produced by an MSE extraction workflow
#'   over repeated simulations.
#'
#' @return A list containing MSE plots, variance plots, bias plots, and the long
#'   plotting data set.
#' @export
plot_mse_grid <- function(extracted_results) {
  check_ggplot2()

  n_rows <- nrow(extracted_results)
  final_dataset <- data.frame(
    objective = c(
      extracted_results$objective_minimax,
      extracted_results$objective_epsilon_net,
      extracted_results$objective_spectral,
      extracted_results$objective_louvain
    ),
    population_objective = c(
      extracted_results$weighted_mse_causal,
      extracted_results$weighted_mse_epsilon_net,
      extracted_results$weighted_mse_spectral,
      extracted_results$weighted_mse_louvain
    ),
    bias = c(
      extracted_results$bias_causal,
      extracted_results$bias_epsilon_net,
      extracted_results$bias_spectral,
      extracted_results$bias_louvain
    ),
    variance = c(
      extracted_results$variance_causal,
      extracted_results$variance_epsilon_net,
      extracted_results$variance_spectral,
      extracted_results$variance_louvain
    ),
    xi_seq = rep(extracted_results$xi_seq, 4),
    sample_size = rep(extracted_results$sample_size, 4),
    network = rep(as.character(extracted_results$network), 4),
    clustering_type = c(
      rep("Causal", n_rows),
      rep("Epsilon net", n_rows),
      rep("Spectral", n_rows),
      rep("Louvain", n_rows)
    ),
    stringsAsFactors = FALSE
  )

  unique_n <- unique(final_dataset$sample_size)
  unique_network <- unique(final_dataset$network)

  make_plot_set <- function(y_var, y_lab) {
    plots <- list()
    acc <- 1L
    for (n_val in unique_n) {
      for (net_name in unique_network) {
        plot_data <- final_dataset[
          final_dataset$sample_size == n_val & final_dataset$network == net_name,
          ,
          drop = FALSE
        ]
        plots[[acc]] <- ggplot2::ggplot(
          plot_data,
          ggplot2::aes(
            x = .data$xi_seq,
            y = .data[[y_var]],
            colour = .data$clustering_type,
            group = .data$clustering_type
          )
        ) +
          ggplot2::geom_line(na.rm = TRUE) +
          ggplot2::geom_point(na.rm = TRUE) +
          ggplot2::ggtitle(paste0("N = ", n_val, ", Network = ", net_name)) +
          ggplot2::xlab("xi") +
          ggplot2::ylab(y_lab) +
          ggplot2::theme_bw()
        acc <- acc + 1L
      }
    }
    plots
  }

  list(
    plot_mse = make_plot_set("population_objective", "Weighted MSE"),
    plot_variance = make_plot_set("variance", "Variance"),
    plot_bias = make_plot_set("bias", "Bias"),
    dataset = final_dataset
  )
}
