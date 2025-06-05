#' Get ranked permutations
#'
#' @description
#' `get_ranked_perm()` calculates the proportion of samples for which a
#' permutation with ranks `rank_range` was observed in `hierarchy_matrix`.
#'
#' @param hierarchy_matrix a matrix where column headers are ranks and each row
#'   displays the treatments assigned to each rank for that iteration.
#' @param rank_range a numeric vector of sequential ranks to consider.
#'
#' @returns A data frame containing
#'   * `Range`: a string presenting the ranks corresponding to the permutation
#'   in `Var1`, presented as min-max.
#'   * `Var1`: a string of the permutation of treatments.
#'   * `Length`: the number of treatments in the permutation.
#'   * `Freq`: the proportion of samples for which the ranked permutation was
#'   observed.
#'
#' @keywords internal
get_ranked_perm <- function(hierarchy_matrix, rank_range){
  ranked_perm <- as.data.frame(table(apply(hierarchy_matrix[, rank_range], 1,
                                           paste0, collapse = ",")) / nrow(hierarchy_matrix))
  rank_int <- rep(paste0(min(rank_range), "-", max(rank_range)), nrow(ranked_perm))
  ranked_perm <- cbind(rank_int, ranked_perm)
  colnames(ranked_perm) <- c("Range", "Var1", "Freq")
  ranked_perm$Length <- max(rank_range) - min(rank_range) + 1
  ranked_perm <- ranked_perm[, c("Range", "Var1", "Length", "Freq")]
  return(ranked_perm)
}

#' Get permutations
#'
#' @description
#' `get_perm()` groups all ranked permutations by the permutation string
#'   (ignoring ranks), and sums the proportion of samples for which they were
#'   observed.
#'
#' @param all_ranked_perm a data frame consisting of the ranks (`Range`) of all
#'   observed permutations (`Var1`) and the corresponding proportion of samples
#'   for which they were observed (`Freq`).
#'
#' @returns A data frame containing
#'   * `Var1`: a string of the permutation of treatments.
#'   * `Length`: the number of treatments in the permutation.
#'   * `Freq`: the proportion of samples for which the permutation was observed.
#'
#' @keywords internal
get_perm <- function(all_ranked_perm) {
  # takes all ranked permutations, groups it by permutation, and calculates the
  # sum of the frequency
  all_perm <- aggregate(Freq ~ Var1 + Length, data = all_ranked_perm, sum)
  return(all_perm)
}

#' Gets credible ranked combinations
#'
#' @description
#' get_ranked_comb()` first sorts the treatments within permutation strings to
#' create a combination string that ignores order. It then groups all
#' combinations by rank interval and sums the proportion of samples for which
#' they were observed.
#'
#' @param all_ranked_perm a data frame consisting of the ranks (`Range`) of all
#'   observed permutations (`Var1`) and the corresponding proportion of samples
#'   for which they were observed (`Freq`).
#' @param trts a vector of all the treatments strings.
#'
#' @returns A data frame containing
#'   * `Range`: a string presenting the ranks corresponding to the combination
#'   in `Combinations`, presented as min-max.
#'   * `Combinations`: a string of the combination of treatments.
#'   * `Length`: the number of treatments in the combination.
#'   * `Freq`: the proportion of samples for which the ranked combination was
#'   observed.
#'
#' @keywords internal
get_ranked_comb <- function(all_ranked_perm, trts) {
  trts <- sort(trts)
  # create a new column `Combinations`, which sorts within the permutation strings
  all_ranked_perm$Combinations <- vapply(strsplit(as.character(all_ranked_perm$Var1), ','),
                                         function(x) paste(x[match(trts, x, nomatch = 0)], collapse = ','), '')

  # groups combinations by rank interval and calculates sum of the frequency
  all_ranked_combo <- aggregate(Freq ~ Length + Combinations + Range, data = all_ranked_perm, sum)
  all_ranked_combo <- all_ranked_combo[, c( "Range", "Combinations", "Length", "Freq")]
  return(all_ranked_combo)
}

#' Get combinations
#'
#' @description
#' `get_combo()`groups all ranked combinations by the combination string
#'  (ignoring ranks), and sums the proportion of samples for which they were
#'  observed.
#'
#' @param all_ranked_combo a data frame consisting of the ranks (`Range`) of all
#'   observed combinations (`Combinations`) and the corresponding proportion of
#'   samples for which they were observed (`Freq`).
#'
#' @return A data frame containing
#'   * `Combinations`: a string of the combination of treatments.
#'   * `Length`: the number of treatments in the combination.
#'   * `Freq`: the proportion of samples for which the combination was observed.
#'
#' @keywords internal
get_combo <- function(all_ranked_combo) {
  # takes in all ranked combinations, groups it by Combination, and sums the frequencies
  all_combo <- aggregate(Freq ~ Combinations + Length, data = all_ranked_combo, sum)
  return(all_combo)
}

#' Determine redundancy status of a (credible) partial hierarchy within
#' (credible) partial hierarchies
#'
#' @description
#' `is_phier_redundant_within_phier()` determines whether a (credible) partial
#' hierarchy specified by `phier_target` is redundant because of any of the
#' (credible) partial hierarchies specified in `phier_larger_list`.
#'
#' @param phier_target a single character vector of treatment names (without ">"
#'   characters) in order of a (credible) partial hierarchy to assess the
#'   redundancy status of.
#' @param phier_larger_list list of data frames of larger (credible) partial
#'   hierarchies by size, to assess the redundancy status of `phier_target`
#'   against.
#'
#' @return `TRUE` if `phier_target` is a redundant because of any of the
#'   (credible) partial hierarchies listed in `larger_phier_list`.
#'
#' @keywords internal
is_phier_redundant_within_phier <- function(phier_target, larger_phier_list) {
  # to find redundancy status faster, look at smaller larger_phier_list first
  # (smaller are more likely to be credible)
  phier_sizes <- sort(unique(do.call(rbind, larger_phier_list)$Length))
  for (i in sort(phier_sizes)) {
    current_phier_df <- larger_phier_list[[as.character(i)]]
    for(j in 1:nrow(current_phier_df)) {
      # extract jth partial hierarchy in current_phier_df
      phier_larger <- str_split_1(as.character(current_phier_df[j, 1]), " > ")
      # finds position of treatments in phier_target within phier_larger
      positions <- match(phier_target, phier_larger)
      if (!anyNA(positions) && all(diff(positions) >= 0)) {  # checks if position[j+1] >= position[j]
        return(TRUE)
        break
      }
    }
  }
}

#' Creates permutations for `get_cred_phier()` function
#'
#' @description
#' `create_perm()` looks at the credibility of a permutation of size 2 involving
#' a treatment with the worst effect in an existing credible permutation `trts`
#' and a new treatment not in `trts`; if the pair of these two treatments is a
#' credible permutation, the new treatment will be added to `trts` to create a
#' new, larger permutation to be assessed as a potentially credible partial
#' hierarchy.
#'
#' @param trts a data frame consisting of one row of the treatment names
#'   corresponding to those in a permutation that we want to build on.
#' @param trt1 a character string belonging to the treatment with the worst
#'   effect in `trts`.
#' @param new_trts a character vector of the remaining treatments to consider
#'   adding to `trts` to build a new permutation.
#' @param credible a character vector listing credible permutations of size two
#'   (e.g., "trt1_name,trt2_name").
#'
#' @details
#' Note the treatment names should match those in the column names of the
#' `effects_matrix` inputted into `prep_data()`.
#'
#' @return Either a data frame of the new, larger permutations to consider as
#' potentially credible partial hierarchies, or "FALSE" to indicate that no new
#' permutations have been created to be assessed as potentially credible partial
#' hierarchies.
#'
#' @keywords internal
create_perm <- function(trts, trt1, new_trts, credible) {
  actual_new_trts <- c()
  for(trt2 in new_trts) { # iterates through each potential new treatment to add
    pair <- paste0(trt1, ",", trt2) # the potential pair to add
    # if the pair exists in the credible list, we can add it to the permutation
    if(any(pair %in% credible)) {
      actual_new_trts <- c(actual_new_trts, trt2)
    }
  }
  size <- length(actual_new_trts)
  if(size == 0){ # no new treatments to add
    return (FALSE)
  }
  # create vectors where each treatment from trts is repeated size times
  list_trts <- lapply(trts, function(x) {
    rep(x, size)
  })
  list_trts[[length(list_trts) + 1]] <- actual_new_trts
  trt_df <- do.call(cbind, list_trts)  # binds the new treatments to the permutations
  return(trt_df)
}

#' Find high probability density (HPD) set
#'
#' @description
#' `hpd()` determines the subset of ranks with the smallest possible cumulative
#' relative frequency that is at least equal to `threshold`.
#'
#' @param ranks a data frame for a particular treatment, consisting of one
#'   column (`Rank`) of all possible ranks and another column (`Freq`) listing
#'   the proportion of samples for which the treatment was ranked `Rank`.
#' @param threshold a proportion between 0 and 1 for which a hierarchy must be
#'   observed in order to be credible.
#' @param freq_sum a numeric value that should always be 1 (the default).
#'
#' @return A list of 1) a string of the rank(s) in the HPD interval, 2) the
#' corresponding observed relative frequency for the ranks in the HPD interval,
#' 3) a vector of the rank(s) in the HPD interval.
#'
#' @keywords internal
hpd <- function(ranks, threshold, freq_sum = 1) {
  if(threshold == 0) {
    hpd_ranks <- c()
    ranks$Freq <- 0
  } else {
    ranks <- ranks[order(ranks$Freq), ] # sorts in increasing order
    while(freq_sum > threshold && nrow(ranks) > 0) {

      # calculate freq_sum without smallest probability
      freq_sum <- freq_sum - ranks[1, 2]

      # if freq_sum >= threshold, we can drop the first row
      if(freq_sum >= threshold) {
        ranks <- ranks[-1, ]
      }
    }
    ranks <- ranks[order(ranks$Rank), ] # re-orders it in terms of rank
    hpd_ranks <- ranks$Rank
  }

  # formatting
  if(length(hpd_ranks) == 1) { # if there is just one element
    concat_ranks <- hpd_ranks
  } else if(length(hpd_ranks) == 0) {
    concat_ranks <- "N/A"
  } else if(all(diff(as.numeric(as.character(hpd_ranks))) == 1)) { # hpd is an interval
    # formats the ranks into an interval
    concat_ranks <- paste(hpd_ranks[1], hpd_ranks[[length(hpd_ranks)]], sep = "-")
  } else { # hpd is not an interval
    # collapses the ranks into one string
    concat_ranks <- paste(hpd_ranks, collapse = ',')
  }

  hpd_vec <- list(concat_ranks, sum(ranks$Freq), hpd_ranks)

  return(hpd_vec)
}
