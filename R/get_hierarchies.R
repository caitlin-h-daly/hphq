#' Get all credible hierarchies
#'
#' @description
#' `get_hierarchies` runs the inputs through `consecutive()`, `get_hier()`, and
#' `singular_treatment()` to produce a list of all hierarchies with relative
#' frequencies greater than or equal to the threshold. These hierarchies are
#' then run through `find_supersets()` to determine which of them are supersets.
#'
#' @param inputs the output from `prep_data()`, which consists of a list of
#'   `hierarchy_matrix`, `effects_matrix`, and `ranking_df`.
#' @param largerBetter a logical value indicating whether larger relative
#'   effects are better (TRUE) or not (FALSE).
#' @param thresholds A numeric vector containing three proportions between 0 and
#'   1 for which a combinatorial hierarchy, partial hierarchy, and individual
#'   ranking probabilities would be credible.
#' @param MID a numeric value indicating the absolute minimally important
#'   difference.
#' @param printPlot a logical value indicating whether the rankograms should be
#' printed (TRUE) or not (FALSE, the default).
#'
#' @return A list of data frames containing the credible hierarchies for Ranked Permutations, Permutations, Ranked Combinations, Combinations, Posets, Individual Treatments, and HPD Intervals
#' @export
#'
#' @examples
#'get_hierarchies(inputs, "Waitlist", TRUE, c(0.5,0.6,0.7), 10)

get_hierarchies <- function(inputs, largerBetter, thresholds, MID, printPlot = FALSE) {

  if(max(thresholds) > 1 || min(thresholds) < 0) {
    stop("Please ensure threshold values are between 0 and 1")
  }

  treatments <- colnames(inputs[[2]])
  n_trt <- length(treatments)

  consec_outputs <- consecutive(inputs[[1]], thresholds[[1]])
  phier <- get_cred_hier(inputs[[2]], MID, thresholds[[2]], largerBetter)
  single <- singular_treatment(inputs[[3]], thresholds[[3]], printPlot)

  first_outputs <- find_supersets(consec_outputs, phier, n_trt)

  all_outputs <- append(first_outputs, single)
  names(all_outputs) <- c("Ranked Permutations", "Permutations",
                          "Ranked Combinations", "Combinations",
                          "Partial Hierarchies", "Individual Ranks", "HPD")
  return(all_outputs)
}
