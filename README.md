# causalClust

`causalClust` provides tools for designing and evaluating cluster-randomized experiments under network interference. The package implements causal clustering procedures for choosing experimental clusters when units interact through a network and the target estimand is a global treatment effect.

The package is organized around a single main workflow:

1. provide an adjacency matrix for the experimental network;
2. choose a spillover-calibration value `xi`, or provide a grid/range of calibration values;
3. run the causal clustering algorithm;
4. inspect the selected number of clusters, rounding method, bias, variance, objective value, and optimization diagnostics;
5. optionally compare the selected design with baseline clusterings such as epsilon-net, spectral clustering, or Louvain clustering.

## Reference

This package implements methods from:

Viviano, D., Lei, L., Imbens, G., Karrer, B., Schrijvers, O., & Shi, L. (2026). *Causal clustering: design of cluster experiments under network interference*. Manuscript, May 28, 2026.

## Installation

Install the development version from GitHub with:

```r
install.packages("remotes")
remotes::install_github("ostasovskyi/causalClustering")
```

Load the package with:

```r
library(causalClust)
```

Most basic utilities use base R and `stats`. Some optional functionality requires additional packages:

```r
install.packages(c("igraph", "ggplot2", "Matrix"))
```

The recommended engine for the main causal clustering algorithm is the SDP engine. This engine requires `Matrix` and `sdpt3r`.

If `sdpt3r` is not available from CRAN in your R setup, install it from GitHub:

```r
install.packages("remotes")
remotes::install_github("AdamRahman/sdpt3r")
```

On Windows, installing `sdpt3r` from source may require Rtools.

## Quick start

The following example simulates a small network and runs the main clustering algorithm for a single calibration value.

```r
library(causalClust)

sim <- simulate_network_data(
  parameters_graph = list(
    n = 60,
    type_graph = "geometric",
    neighb = 2
  )
)

W <- sim$W

fit <- causal_clustering_algorithm(
  W = W,
  xi = 1,
  min_k = 2,
  max_k = 12,
  engine = "sdp",
  methods = available_discretization_methods(),
  seed = 123
)

fit$algorithm
fit$selected_k
fit$selected_method
fit$objective
fit$components
fit$Gamma_n
fit$certificate_valid
head(fit$clusters)
```

The returned object contains the selected clustering, the selected number of clusters, the selected rounding method, the realized value of the design objective, and diagnostic information from the optimization step.

## Calibration grids

When the calibration value is uncertain, pass a grid or range instead of a single scalar. The unified function then uses the endpoint-regret version of the algorithm.

```r
fit_grid <- causal_clustering_algorithm(
  W = W,
  xi_grid = c(0.2, 0.5, 1, 2),
  min_k = 2,
  max_k = 12,
  engine = "sdp",
  methods = available_discretization_methods(),
  seed = 123
)

fit_grid$algorithm
fit_grid$xi_range
fit_grid$selected_k
fit_grid$selected_method
fit_grid$objective
fit_grid$rho_values
fit_grid$components
```

You can also specify the calibration in the paper's scale using `calibration`, `calibration_grid`, or `calibration_range`. Internally, the package uses

```text
xi = 1 / calibration.
```

## Objective components

For a clustering `clusters`, the package evaluates the design objective

```text
xi * variance + bias^2,
```

where `variance` is the normalized sum of squared cluster sizes and `bias` is the average fraction of neighbors assigned to different clusters.

```r
W_small <- matrix(
  c(
    0, 1, 1, 0,
    1, 0, 0, 1,
    1, 0, 0, 1,
    0, 1, 1, 0
  ),
  nrow = 4,
  byrow = TRUE
)

clusters <- c(1, 1, 2, 2)

components <- clustering_objective_components(W_small, clusters)
objective <- compute_causal_clustering_objective(
  W = W_small,
  xi = 1,
  clusters = clusters
)

components
objective
all.equal(objective, components$variance + components$bias^2)
```

## Core functions

### `causal_clustering_algorithm()`

Main user-facing function. It accepts either a single calibration value or a calibration grid/range. With a scalar calibration, it runs the fixed-calibration algorithm. With a grid or range, it runs the endpoint-regret algorithm.

| Argument group | Inputs |
|---|---|
| Network | `W`, a numeric square adjacency matrix. Nonzero entries define the binary adjacency support. The diagonal is set to zero internally, and nonsymmetric input is symmetrized. |
| Calibration | Exactly one of `xi`, `calibration`, `xi_grid`, `xi_range`, `calibration_grid`, or `calibration_range`. |
| Search range | `min_k`, `max_k`, and optionally `n_eig`. These define the candidate number of clusters and the embedding dimension used for rounding. |
| Engine | `engine = "sdp"` for the semidefinite relaxation, or `engine = "spectral"` for a lightweight scalar-calibration heuristic. |
| Rounding | `methods`, usually from `available_discretization_methods()`. Requested methods are tried and the rounded candidate with the smallest realized objective is selected. |
| Optional constraints | `k_constraint`, `gamma_bar`, and `box_constraints`. These control optional size and box constraints for the SDP engine. |
| Reproducibility and storage | `seed`, `try_sign_flip`, `keep_search_matrices`, and `keep_sdp_solutions`. |

Main outputs:

| Output | Meaning |
|---|---|
| `clusters` | Integer cluster labels for the selected design. |
| `selected_k` | Selected number of clusters. |
| `selected_method` | Rounding or discretization method selected after candidate evaluation. |
| `objective` | Selected value of the design criterion. For grid/range calls, this is the selected endpoint-regret criterion. |
| `components` | Bias, variance, number of clusters, and cluster sizes for the selected clustering. |
| `Gamma_n` | SDP approximation certificate, when available. |
| `certificate_valid` | Whether the approximation certificate is available and valid for the call. |
| `sdp_lower_bound` or `endpoint_lower_bounds` | SDP lower bound diagnostics, when available. |
| `search_result`, `solutions`, `rho_values` | Candidate-level optimization output. |

### `clustering_objective_components()`

Computes the two parts of the design objective for an existing clustering.

Inputs:

- `W`: adjacency matrix;
- `clusters`: cluster labels with length `nrow(W)`.

Outputs:

- `variance`: `n^{-2} sum_k n_k^2`;
- `bias`: average fraction of neighbors assigned to different clusters;
- `num_clusters`: number of distinct clusters;
- `cluster_sizes`: cluster-size vector.

### `compute_causal_clustering_objective()`

Computes the realized scalar design objective for an existing clustering.

Inputs:

- `W`: adjacency matrix;
- `xi`: nonnegative calibration value;
- `clusters`: cluster labels;
- `objective_type`: currently must be `"squared"`.

Output:

- numeric scalar equal to `xi * variance + bias^2`.

### `simulate_network_data()`

Simulates covariates, a network, treatment assignments, and outcomes under a first-order network-interference model.

Inputs:

- `parameters_graph`: graph-generation settings;
- `parameters_model`: outcome-model settings;
- `D`: optional treatment vector;
- `W`: optional adjacency matrix;
- `seed`: optional random seed.

Output:

- a list containing simulated outcomes, covariates, treatment assignments, network matrices, and estimand summaries.

## Secondary functions

| Function | Purpose |
|---|---|
| `available_discretization_methods()` | Lists supported rounding and discretization methods. |
| `causal_clustering_algorithm1()` | Fixed-calibration wrapper. |
| `causal_clustering_algorithm2()` | Endpoint-regret wrapper for calibration ranges. |
| `adaptive_causal_clustering()` | Convenience wrapper for calibration grids or ranges. |
| `cluster_epsilon_net()` | Epsilon-net baseline clustering. |
| `cluster_louvain_membership()` | Louvain community-detection baseline. |
| `cluster_spectral()` | Spectral-clustering baseline. |
| `assign_cluster_treatment()` | Assigns treatment at the cluster level. |
| `run_single_network()` | Runs one network simulation and compares clustering procedures. |
| `simulate_network_grid()` | Runs a grid of network simulations. |
| `plot_objective_path()` | Plots objective and selected-cluster paths. |
| `plot_mse_path()` | Plots weighted MSE, bias, variance, and cluster-count paths. |

## Input requirements

The main input is an adjacency matrix `W`:

- `W` must be a numeric square matrix;
- nonzero entries define the binary adjacency support;
- diagonal entries are set to zero internally;
- nonsymmetric input is symmetrized internally;
- missing or non-finite entries are not allowed.

## Calibration parameters

The tuning parameter `xi` controls the tradeoff between the two pieces of the design criterion:

```text
xi * variance + bias^2.
```

Larger values of `xi` place more weight on the variance component and typically favor more clusters. Smaller values place more weight on reducing cross-cluster exposure and typically favor fewer, larger clusters.

The package also accepts `calibration`, where

```text
xi = 1 / calibration.
```

## Discretization methods

Available methods can be inspected with:

```r
available_discretization_methods()
```

Currently supported methods are:

```r
c(
  "kmeans",
  "hierarchical",
  "spectral_norm_kmeans",
  "spectral_unnorm_kmeans",
  "spectral_norm_hierarchical",
  "spectral_unnorm_hierarchical"
)
```

When more than one method is supplied, each method is applied after the relaxation step. The package evaluates the realized design objective for each rounded candidate and selects the best one.

## SDP engine diagnostics

With `engine = "sdp"`, returned objects may include SDP diagnostics:

```r
fit$Gamma_n
fit$certificate_valid
fit$sdp_lower_bound
```

The argument `k_constraint` controls whether the SDP relaxation is solved once and rounded over the candidate K grid, or solved separately with K-specific constraints.

- `k_constraint = FALSE`: solve one relaxation and round it over the K grid.
- `k_constraint = TRUE`: solve a K-specific SDP relaxation for each candidate K, with optional size and box constraints.

## Lightweight spectral approximation

The package can also be used without an SDP solver for exploratory work. The spectral engine evaluates the same squared-bias objective after rounding, but it does not solve the semidefinite relaxation and does not provide the SDP lower bound or approximation certificate.

```r
fit_spectral <- causal_clustering_algorithm(
  W = W,
  xi = 1,
  min_k = 2,
  max_k = 12,
  engine = "spectral",
  methods = available_discretization_methods(),
  seed = 123
)

fit_spectral$selected_k
fit_spectral$selected_method
fit_spectral$objective
fit_spectral$components
```

## Comparing clustering methods

The following example uses the spectral engine so that it can run without an SDP solver.

```r
out <- run_single_network(
  seed = 123,
  parameters_graph = list(
    type_graph = "geometric",
    n = 100,
    neighb = 2
  ),
  xi_seq = seq(0.1, 5, by = 0.5),
  min_k = 2,
  max_k = 20,
  include_louvain = TRUE,
  engine = "spectral",
  objective_type = "squared",
  methods = available_discretization_methods()
)

head(out$results)

plots <- plot_objective_path(out$results)
plots$plot_objective
plots$plot_num_clusters
```