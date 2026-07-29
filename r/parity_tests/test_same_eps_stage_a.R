source(file.path("r", "chronodiag", "R", "core_metrics.R"))
source(file.path("r", "chronodiag", "R", "composition_diagnostics.R"))
source(file.path("r", "chronodiag", "R", "influence_diagnostics.R"))

output_dir <- file.path("analysis", "same_eps_stage_a_v01")
exact <- read.csv(
  gzfile(file.path(output_dir, "exact_same_eps_long_v01.csv.gz")),
  stringsAsFactors = FALSE
)
key <- read.csv(
  file.path(output_dir, "key_results_v01.csv"),
  stringsAsFactors = FALSE
)
bootstrap_membership <- read.csv(
  gzfile(file.path(output_dir, "bootstrap_multiplicities_v01.csv.gz")),
  stringsAsFactors = FALSE
)
bootstrap_expected <- read.csv(
  gzfile(file.path(output_dir, "bootstrap_results_v01.csv.gz")),
  stringsAsFactors = FALSE
)
loto_expected <- read.csv(
  file.path(output_dir, "loto_results_v01.csv"),
  stringsAsFactors = FALSE
)

scenario_data <- function(name) {
  out <- exact[exact$scenario == name, c("year", "tree_id", "value")]
  rownames(out) <- NULL
  out
}

# Exchangeable subset-to-full correlation must equal sqrt(SSS_9|18).
exchangeable <- scenario_data("exchangeable")
subset_membership <- data.frame(
  replicate_id = rep("S09", 9),
  unit_id = sprintf("T%02d", 1:9),
  stringsAsFactors = FALSE
)
subset_result <- subset_stability(
  exchangeable,
  subset_membership,
  min_overlap = 100L
)
expected_exchangeable <- key[
  key$scenario == "exchangeable",
  "subset_n9_theoretical_correlation"
]
stopifnot(
  abs(subset_result$correlation - expected_exchangeable) < 1e-10
)

# A cluster-aligned 9-vs-9 split is the exact lower-tail two-cluster case.
two_clusters <- scenario_data("two_clusters")
split_membership <- rbind(
  data.frame(
    replicate_id = "H_CLUSTER",
    group = "A",
    unit_id = sprintf("T%02d", 1:9),
    stringsAsFactors = FALSE
  ),
  data.frame(
    replicate_id = "H_CLUSTER",
    group = "B",
    unit_id = sprintf("T%02d", 10:18),
    stringsAsFactors = FALSE
  )
)
split_result <- split_half_stability(
  two_clusters,
  split_membership,
  min_overlap = 100L
)
expected_cluster_min <- key[
  key$scenario == "two_clusters",
  "split_half_min"
]
stopifnot(abs(split_result$correlation - expected_cluster_min) < 1e-10)

# The first ten frozen bootstrap draws agree with the covariance-level analysis.
first_replicates <- sprintf("B%04d", 1:10)
membership_10 <- bootstrap_membership[
  bootstrap_membership$replicate_id %in% first_replicates,
]
bootstrap_actual <- tree_bootstrap(
  exchangeable,
  membership_10,
  min_overlap = 100L
)
bootstrap_reference <- bootstrap_expected[
  bootstrap_expected$scenario == "exchangeable" &
    bootstrap_expected$replicate_id %in% first_replicates,
]
bootstrap_reference <- bootstrap_reference[
  order(bootstrap_reference$replicate_id),
]
bootstrap_actual <- bootstrap_actual[
  order(bootstrap_actual$replicate_id),
]
stopifnot(identical(
  as.character(bootstrap_actual$replicate_id),
  as.character(bootstrap_reference$replicate_id)
))
for (column in c(
  "correlation",
  "standardized_rms_difference",
  "variance_ratio"
)) {
  stopifnot(max(abs(
    bootstrap_actual[[column]] - bootstrap_reference[[column]]
  )) < 1e-10)
}

# LOTO influence is metric-specific in the high-/lower-loading scenario.
influential <- scenario_data("influential_trees")
loto_actual <- leave_one_tree_out(influential, min_overlap = 100L)
loto_reference <- loto_expected[
  loto_expected$scenario == "influential_trees",
]
loto_reference <- loto_reference[
  order(loto_reference$omitted_unit_id),
]
loto_actual <- loto_actual[order(loto_actual$omitted_unit_id),]
stopifnot(identical(
  as.character(loto_actual$omitted_unit_id),
  as.character(loto_reference$omitted_unit_id)
))
for (column in c(
  "correlation",
  "standardized_rms_difference",
  "delta_rbar",
  "delta_eps",
  "delta_effective_dimension"
)) {
  stopifnot(max(abs(
    loto_actual[[column]] - loto_reference[[column]]
  )) < 1e-9)
}
stopifnot(
  loto_actual$delta_eps[loto_actual$omitted_unit_id == "T01"] <
    loto_actual$delta_eps[loto_actual$omitted_unit_id == "T18"]
)
stopifnot(
  loto_actual$standardized_rms_difference[
    loto_actual$omitted_unit_id == "T18"
  ] >
    loto_actual$standardized_rms_difference[
      loto_actual$omitted_unit_id == "T01"
    ]
)

cat("R same-EPS Stage A parity checks passed.\n")
