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
#'   * `Size`: the number of treatments in the permutation.
#'   * `pi_hat`: the proportion of samples for which the ranked permutation was
#'   observed (i.e., empirical probabilities).
#'
#' @keywords internal
get_ranked_perm <- function(hierarchy_matrix, rank_range){
  ranked_perm <- as.data.frame(table(apply(hierarchy_matrix[, rank_range], 1,
                                           paste0, collapse = ",")) / nrow(hierarchy_matrix))
  rank_int <- rep(paste0(min(rank_range), "-", max(rank_range)), nrow(ranked_perm))
  ranked_perm <- cbind(rank_int, ranked_perm)
  colnames(ranked_perm) <- c("Range", "Var1", "pi_hat")
  ranked_perm$Size <- max(rank_range) - min(rank_range) + 1
  ranked_perm <- ranked_perm[, c("Range", "Var1", "Size", "pi_hat")]
  return(ranked_perm)
}

#' Get permutations
#'
#' @description
#' `get_perm()` groups all ranked permutations by the permutation string
#'   (ignoring ranks), and sums the proportion of samples for which they were
#'   observed (i.e., empirical probabilities).
#'
#' @param all_ranked_perm a data frame consisting of the ranks (`Range`) of all
#'   observed permutations (`Var1`) and the corresponding proportion of samples
#'   for which they were observed (`pi_hat`).
#'
#' @returns A data frame containing
#'   * `Var1`: a string of the permutation of treatments.
#'   * `Size`: the number of treatments in the permutation.
#'   * `pi_hat`: the proportion of samples for which the permutation was
#'   observed  (i.e., empirical probabilities).
#'
#' @keywords internal
get_perm <- function(all_ranked_perm) {
  # takes all ranked permutations, groups it by permutation, and calculates the
  # empirical probability
  all_perm <- aggregate(pi_hat ~ Var1 + Size, data = all_ranked_perm, sum)
  return(all_perm)
}

#' Gets credible ranked combinations
#'
#' @description
#' get_ranked_comb()` first sorts the treatments within permutation strings to
#' create a combination string that ignores order. It then groups all
#' combinations by rank interval and sums the proportion of samples for which
#' they were observed (i.e., empirical probabilities).
#'
#' @param all_ranked_perm a data frame consisting of the ranks (`Range`) of all
#'   observed permutations (`Var1`) and the corresponding proportion of samples
#'   for which they were observed (`pi_hat`).
#' @param trts a vector of all the treatments strings.
#'
#' @returns A data frame containing
#'   * `Range`: a string presenting the ranks corresponding to the combination
#'   in `Combinations`, presented as min-max.
#'   * `Combinations`: a string of the combination of treatments.
#'   * `Size`: the number of treatments in the combination.
#'   * `pi_hat`: the proportion of samples for which the ranked combination was
#'   observed.
#'
#' @keywords internal
get_ranked_comb <- function(all_ranked_perm, trts) {
  trts <- sort(trts)
  # create a new column `Combinations`, which sorts within the permutation strings
  all_ranked_perm$Combinations <- vapply(strsplit(as.character(all_ranked_perm$Var1), ','),
                                         function(x) paste(x[match(trts, x, nomatch = 0)], collapse = ','), '')

  # groups combinations by rank interval and calculates the empirical probability
  all_ranked_combo <- aggregate(pi_hat ~ Size + Combinations + Range, data = all_ranked_perm, sum)
  all_ranked_combo <- all_ranked_combo[, c( "Range", "Combinations", "Size", "pi_hat")]
  return(all_ranked_combo)
}

#' Get combinations
#'
#' @description
#' `get_combo()`groups all ranked combinations by the combination string
#'  (ignoring ranks), and sums the proportion of samples for which they were
#'  observed (i.e., empirical probabilities).
#'
#' @param all_ranked_combo a data frame consisting of the ranks (`Range`) of all
#'   observed combinations (`Combinations`) and the corresponding proportion of
#'   samples for which they were observed (`pi_hat`).
#'
#' @return A data frame containing
#'   * `Combinations`: a string of the combination of treatments.
#'   * `Size`: the number of treatments in the combination.
#'   * `pi_hat`: the proportion of samples for which the combination was
#'   observed  (i.e., empirical probabilities).
#'
#' @keywords internal
get_combo <- function(all_ranked_combo) {
  # takes in all ranked combinations, groups it by Combination, and sums the empirical probabilities
  all_combo <- aggregate(pi_hat ~ Combinations + Size, data = all_ranked_combo, sum)
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
#' @param larger_phier_list list of data frames of larger (credible) partial
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
  phier_sizes <- sort(unique(do.call(rbind, larger_phier_list)$Size))
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

#' Find highest density region (HDR) set
#'
#' @description
#' `hdr()` determines the subset of ranks with the smallest possible cumulative
#' empirical probability that is at least equal to `threshold`.
#'
#' @param ranks a data frame for a particular treatment, consisting of one
#'   column (`Rank`) of all possible ranks and another column (`pi_hat`) listing
#'   the proportion of samples for which the treatment was ranked `Rank`.
#' @param threshold a proportion between 0 and 1 for which a hierarchy must be
#'   observed in order to be credible.
#' @param pi_hat_sum a numeric value that should always be 1 (the default).
#'
#' @return A list of 1) a string of the rank(s) in the HDR set, 2) the
#' corresponding empirical probability for the ranks in the HDR set,
#' 3) a vector of the rank(s) in the HDR set.
#'
#' @keywords internal
hdr <- function(ranks, threshold, pi_hat_sum = 1) {
  if(threshold == 0) {
    hdr_ranks <- c()
    ranks$pi_hat <- 0
  } else {
    alt_cdr <- list()
    excl_ranks <- ranks[0, ]
    ranks <- ranks[order(ranks$pi_hat), ] # sorts in increasing order
    while(pi_hat_sum > threshold && nrow(ranks) > 0) {

      # calculate pi_hat_sum without smallest probability
      pi_hat_sum <- pi_hat_sum - ranks[1, 2]

      # if pi_hat_sum >= threshold, we can drop the first row
      if(pi_hat_sum >= threshold) {
        excl_ranks <- rbind(excl_ranks, ranks[1, ])
        ranks <- ranks[-1, ]
      }
    }
    ranks <- ranks[order(ranks$Rank), ] # re-orders it in terms of rank
    hdr_ranks <- ranks$Rank
  }

  # formatting
  if(length(hdr_ranks) == 1) { # if there is just one element
    concat_ranks <- hdr_ranks
  } else if(length(hdr_ranks) == 0) {
    concat_ranks <- "N/A"
  } else if(all(diff(as.numeric(as.character(hdr_ranks))) == 1)) {
    # hdr is an interval
    # formats the ranks into an interval
    concat_ranks <- paste(hdr_ranks[1], hdr_ranks[[length(hdr_ranks)]], sep = "-")
  } else { # hdr is not an interval
    # collapses the ranks into one string
    concat_ranks <- paste(hdr_ranks, collapse = ',')
  }

  hdr_vec <- list(concat_ranks, sum(ranks$pi_hat), hdr_ranks)

  # alternative credible density regions (CDR) of the same size
  #if(length(hdr_ranks) > 1) {
  hdr_pi_hat <- sum(ranks$pi_hat)
  min_hdr_rank_pi_hat <- min(ranks$pi_hat) # smallest pi_hat contributing to hdr
  req_alt_rank_pi_hat <- threshold - (hdr_pi_hat - min_hdr_rank_pi_hat) # pi_hat required to substitute least likely rank in hdr

  alt_cdr_rank <- excl_ranks[which(excl_ranks$pi_hat >= req_alt_rank_pi_hat), ]

  if(nrow(alt_cdr_rank) > 0) {
    list_ind <- 1

    ranks <- ranks[order(ranks$pi_hat), ] # re-orders it in terms of pi_hat

    ranks_to_sub_ind <- which(ranks$pi_hat == min_hdr_rank_pi_hat) # in case multiple ranks have same pi_hat

    for(i in ranks_to_sub_ind) {
      base_ranks <- ranks[-i, ]

      for(j in 1:nrow(alt_cdr_rank)){
        new_ranks <- rbind(alt_cdr_rank[j, ], base_ranks)
        new_ranks <- new_ranks[order(new_ranks$Rank), ] # re-orders it in terms of rank

        new_cdr_ranks <- new_ranks$Rank

        # format
        if(length(new_cdr_ranks) == 1) { # if there is just one element
          new_concat_ranks <- new_cdr_ranks
        } else if(length(new_cdr_ranks) == 0) {
          new_concat_ranks <- "N/A"
        } else if(all(diff(as.numeric(as.character(new_cdr_ranks))) == 1)) {
          # cdr is an interval
          # formats the ranks into an interval
          new_concat_ranks <- paste(new_cdr_ranks[1],
                                    new_cdr_ranks[[length(new_cdr_ranks)]],
                                    sep = "-")
        } else { # cdr is not an interval
          # collapses the ranks into one string
          new_concat_ranks <- paste(new_cdr_ranks, collapse = ',')
        }

        alt_cdr[[list_ind]] <- list(new_concat_ranks, sum(new_ranks$pi_hat))

        list_ind <- list_ind + 1

      }
    }
  }

  #}

  return(list(hdr = hdr_vec, alt_cdr = alt_cdr))

}
