#' causalClust: causal clustering for network experiments
#'
#' Tools for constructing and evaluating cluster randomized designs under
#' network interference. The main design function is
#' [causal_clustering_algorithm()], which dispatches to the fixed-calibration
#' or endpoint-regret algorithm depending on whether the calibration input is
#' scalar or a grid/range.
#'
#' @keywords internal
"_PACKAGE"

if (getRversion() >= "2.15.1") {
  utils::globalVariables(c(".data"))
}
