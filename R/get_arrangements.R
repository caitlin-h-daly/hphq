#' Get credible arrangement hierarchies
#'
#' @description
#' `get_arrangements()` finds all ranked permutations, permutations, ranked
#' combinations, and combinations with relative frequencies greater than or
#' equal to a given threshold.
#'
#' @details
#' Note that the progress of the permutations compilation will be printed for
#' each permutation size. The user may suppress these messages using
#' `suppressMessages(get_arrangements())`.
#'
#' @param hierarchy_matrix a matrix where column headers are ranks and each row
#'   displays the treatments assigned to each rank for that iteration.
#' @param threshold a proportion between 0 and 1 for which a hierarchy must be
#'   observed in order to be credible.
#'
#' @return A list of data frames containing the credible hierarchies for ranked
#' permutations, permutations, ranked combinations, and combinations.
#'
#' @importFrom stats aggregate
#' @export
#'
#' @examples
#' inputs <- prep_data(effects_matrix = dat_Thijs2008[, -1], reference = "Placebo", largerbetter = FALSE)
#' get_arrangements(hierarchy_matrix = inputs$hierarchy_matrix, threshold = 0.9)
get_arrangements <- function(hierarchy_matrix, threshold) {

  if(threshold > 1 || threshold < 0) {
    stop("Please ensure threshold value is between 0 and 1")
  }

  treatments <- hierarchy_matrix[1, ]
  n_trt <- ncol(hierarchy_matrix)
  tolerance <- .Machine$double.eps ^ 0.5
  threshold <- threshold - tolerance
  consec_output <- vector("list", length = 4)

  # total number of permutations of sizes 2 to n_trt covering all ranks
  # (e.g., 1:2, 2:3, ..., (n_trt-1):n_trt, 1:3, 2:4, ..., (n_trt-2):n_trt, ...,  1:n_trt)
  n_perm_grps <- ((n_trt - 1) * n_trt) / 2

  # tabulate all observed permutations and their ranks from size 2 to size n_trt
  all_ranked_perm_list <- vector("list", length = n_perm_grps)
  index <- 1
  for(size in 2:n_trt) {
    # start is index of first rank to consider
    for(start in 1:(n_trt - size + 1)) {
      # end is index of last rank to consider
      end <- start + size - 1
      # tabulate all permutations observed between rank `start` and rank `end`
      all_ranked_perm_list[[index]] <- get_ranked_perm(hierarchy_matrix, start:end)
      index <- index + 1
    }
    message(paste0("Permutations of size ", size ," completed"))
  }

  # all observed ranked permutations
  all_ranked_perm <- do.call(rbind, all_ranked_perm_list)
  # credible ranked permutations
  cred_ranked_perm <- subset(all_ranked_perm, all_ranked_perm$Freq > threshold)
  names(cred_ranked_perm)[names(cred_ranked_perm) == 'Var1'] <- 'Ranked Permutations'
  # order credible ranked permutations by frequency
  cred_ranked_perm <- cred_ranked_perm[order(cred_ranked_perm$Freq, decreasing = TRUE), ]
  # remove row names
  rownames(cred_ranked_perm) <- NULL
  consec_output[[1]] <- cred_ranked_perm
  message("Ranked permutations completed")

  # all observed permutations
  all_perm <- get_perm(all_ranked_perm)
  # credible permutations
  cred_perm <- subset(all_perm, all_perm$Freq > threshold)
  names(cred_perm)[names(cred_perm) == 'Var1'] <- 'Permutations'
  # order credible permutations by frequency
  cred_perm <- cred_perm[order(cred_perm$Freq, decreasing = TRUE), ]
  # remove row names
  rownames(cred_perm) <- NULL
  consec_output[[2]] <- cred_perm
  message("Permutations completed")

  # all ranked permutations up to size n_trt - 3
  #all_ranked_perm_combo <- do.call(rbind, all_ranked_perm_list[1:(n_perm_grps - 3)])
  # all ranked combinations
  #all_ranked_combo <- indexed_combo(all_ranked_perm_combo, treatments)
  all_ranked_combo <- get_ranked_comb(all_ranked_perm, treatments)
  # all credible ranked combinations
  cred_ranked_combo <- subset(all_ranked_combo, all_ranked_combo$Freq > threshold)
  # order credible ranked combinations by frequency
  cred_ranked_combo <- cred_ranked_combo[order(cred_ranked_combo$Freq, decreasing = TRUE), ]
  # remove row names
  rownames(cred_ranked_combo) <- NULL
  consec_output[[3]] <- cred_ranked_combo
  names(consec_output[[3]])[names(consec_output[[3]]) == 'Combinations'] <- 'Ranked Combinations'
  message("Ranked combinations completed")

  # all combinations
  all_combo <- get_combo(all_ranked_combo)
  # all credible combinations
  cred_combo <- subset(all_combo, all_combo$Freq > threshold)
  # order credible ranked combinations by frequency
  cred_combo <- cred_combo[order(cred_combo$Freq, decreasing = TRUE), ]
  # remove row names
  rownames(cred_combo) <- NULL
  consec_output[[4]] <- cred_combo
  message("Combinations completed")

  # add appropriate brackets for combinatorial type
  if(nrow(consec_output[[1]]) > 0) {
    consec_output[[1]]$`Ranked Permutations` <- paste0("(", consec_output[[1]]$`Ranked Permutations`, ")")
  }
  if(nrow(consec_output[[2]]) > 0) {
    consec_output[[2]]$Permutations <- paste0("(", consec_output[[2]]$Permutations, ")")
  }
  if(nrow(consec_output[[3]]) > 0) {
    consec_output[[3]]$`Ranked Combinations` <- paste0("{", consec_output[[3]]$`Ranked Combinations`, "}")
  }
  if(nrow(consec_output[[4]]) > 0) {
    consec_output[[4]]$Combinations <- paste0("{", consec_output[[4]]$Combinations, "}")
  }
  names(consec_output) <- c("Ranked Permutations", "Permutations",
                            "Ranked Combinations", "Combinations")

  return(consec_output)
}
