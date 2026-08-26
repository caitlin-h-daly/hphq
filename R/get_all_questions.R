#' Get all credible hierarchies
#'
#' @description
#' `get_all_questions` runs the inputs through `get_arrangements()`,
#' `get_partial_hierarchies()`, and `get_all_questions()` to produce a list of
#' all hierarchies with empirical probabilities greater than or equal to the
#' threshold. These hierarchies are then run through `find_redundancies()` to
#' determine which of them are redundant.
#'
#' @param inputs the output from `prep_data()`, which consists of a list of
#'   `hierarchy_matrix`, `effects_matrix`, and `ranking_df`.
#' @param n_trt a numeric value indicating the number of treatments in the
#'   network.
#' @param larger_better a logical value indicating whether larger relative
#'   effects are better (TRUE) or not (FALSE).
#' @param thresholds A numeric vector containing three proportions between 0 and
#'   1 for which a combinatorial hierarchy, partial hierarchy, and individual
#'   ranking probabilities would be credible.
#' @param mid a numeric value indicating the absolute minimally important
#'   difference. Default is 0.
#' @param print_plot a logical value indicating whether the rankograms should be
#'   printed (TRUE) or not (FALSE, the default).
#' @param trim_redundant a logical value indicating whether the redundant
#'   hierarchies should be trimmed from the output (TRUE) or not (FALSE, the
#'   default).
#'
#' @return A list of data frames containing the credible hierarchies for ranked
#' permutations, permutations, ranked combinations, combinations, partial
#' hierarchies, and HDR sets.
#' @export
#'
#' @examples
#' inputs <- prep_data(effects_matrix = dat_Thijs2008[, -1],
#'                     reference = "Placebo",
#'                     larger_better = FALSE)
#' get_all_questions(inputs = inputs,
#'                   larger_better = FALSE,
#'                   thresholds = c(0.95, 0.95, 0.95),
#'                   mid = 0,
#'                   print_plot = FALSE,
#'                   trim_redundant = FALSE)
get_all_questions <- function(inputs, n_trt, larger_better, thresholds,
                              mid = 0, print_plot = FALSE,
                              trim_redundant = FALSE) {

  if(max(thresholds) > 1 || min(thresholds) < 0) {
    stop("Please ensure threshold values are between 0 and 1")
  }

  if(length(unique(thresholds)) > 1) {
    stop("To find redundancies, please ensure thresholds are the same")
  } else {
    threshold <- unique(thresholds)
  }

  treatments <- colnames(inputs[[2]])
  n_trt <- length(treatments)

  arrangements <- get_arrangements(inputs$hierarchy_matrix,
                                   thresholds[[1]])
  phier <- get_partial_hierarchies(inputs$effects_matrix, mid,
                                   thresholds[[2]], larger_better)
  single <- get_ranks_by_treatment(inputs$ranking_df, thresholds[[3]],
                                   print_plot)

  all_outputs <- find_redundancies(algo_1 = arrangements,
                                   algo_2 = phier,
                                   algo_3 = single,
                                   n_trt = n_trt,
                                   threshold = threshold,
                                   trim_redundant = trim_redundant)

  all_outputs[[7]] <- single$`Alternative CDR of same size as HDR`

  names(all_outputs) <- c("Ranked Permutations", "Permutations",
                          "Ranked Combinations", "Combinations",
                          "Partial Hierarchies", "HDR",
                          "Alternative CDR of same size as HDR")

  return(all_outputs)
}
