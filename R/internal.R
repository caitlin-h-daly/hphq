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
#'   * `Freq`: the proportion of samples for which the permutation was observed.
#'
#' @keywords internal
get_perm <- function(all_ranked_perm) {
  # takes all ranked permutations, groups it by permutation, and calculates the
  # sum of the frequency
  all_perm <- aggregate(Freq ~ Var1, data = all_ranked_perm, sum)
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
  all_ranked_combo <- aggregate(Freq ~ Combinations + Range, data = all_ranked_perm, sum)
  all_ranked_combo <- all_ranked_combo[, c("Range", "Combinations", "Freq")]
  return(ranked_combo)
}

#' Get combinations
#'
#' @description
#' `get_combo()`
#'
#' @param all_ranked_combo a data frame consisting of the ranks (`Range`) of all
#'   observed combinations (`Combinations`) and the corresponding proportion of
#'   samples for which they were observed (`Freq`).
#'
#' @return A data frame containing
#'   * `Combinations`: a string of the combination of treatments.
#'   * `Freq`: the proportion of samples for which the combination was observed.
#'
#' @keywords internal
get_combo <- function(all_ranked_combo) {
  # takes in all ranked combinations, groups it by Combination, and sums the frequencies
  all_combo <- aggregate(Freq ~ Combinations, data = all_ranked_combo, sum)
  return(all_combo)
}

#' Determine superset status of a (credible) partial hierarchy within (credible)
#' partial hierarchies
#'
#' @description
#' `is_phier_sup_of_phier()` determines whether a (credible) partial hierachy
#' specified by `phier_target` is a superset of any of the (credible) partial
#' hierarchies specified in `phier_larger_list`.
#'
#' @param phier_target a single character vector of treatment names (without ">"
#'   characters) in order of a (credible) partial hierarchy to assess the
#'   superset status of.
#' @param phier_larger_list list of data frames of larger (credible) partial
#'   hierarchies by size, to assess the superset status of `phier_target`
#'   against.
#'
#' @return `TRUE` if `phier_target` is a superset of any of the (credible)
#'   partial hierarchies listed in `larger_phier_list`.
#'
#' @keywords internal
is_phier_sup_of_phier <- function(phier_target, larger_phier_list) {
  # to find superset status faster, look at smaller larger_phier_list first
  phier_sizes <- sort(unique(do.call(rbind, larger_phier_list)$phier_size))
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

#' Creates permutations for poset function
#' @param trts list of treatment names
#' @param trt1 first treatment in the new pair
#' @param new_trts potential new treatments to add
#' @param perm_size size of the permutation
#' @param credible list of credible size two permutations
#' @keywords internal
create_perm <- function (trts,trt1,new_trts,perm_size,credible) {
  actual_new_trts <- c()
  for (trt2 in new_trts) { # iterates through each potential new treatment to add
    pair <- paste0(trt1,",",trt2) # the potential pair to add
    if (any(pair %in% credible)) { # if the pair exists in the credible list, we can add it to the permutation
      actual_new_trts<-c(actual_new_trts,trt2)
    } # if the pair isn't credible, we do not need to add it to the permutation
  }
  size<-length(actual_new_trts)
  if(size == 0) { # no new treatments to add
    return (FALSE)
  }
  list_trts <- lapply(trts,function(x) { # create vectors where each treatment from trts is repeated size times
    rep(x,size)
  })
  list_trts[[length(list_trts)+1]] <- actual_new_trts
  trt_df <- do.call(cbind, list_trts) # binds the new treatments to the permutations
  return(trt_df)
}

#' Finds HPD intervals
#' @param ranks column of ranks and their frequency
#' @param threshold cutoff for what is considered credible
#' @param freq_sum value is always 1
#' @keywords internal
hpd <- function(ranks,threshold,freq_sum) {
  if (threshold == 0){
    hpd_ranks<-c()
    ranks$Freq<-0
  } else {
    ranks <- ranks[order(ranks$Freq), ] # sorts in increasing order
    while (freq_sum > threshold && nrow(ranks) > 0) {
      freq_sum <- freq_sum - ranks[1,2] # calculates freq_sum without smallest probability

      if (freq_sum >= threshold) { # if >= threshold, we can drop the first row
        ranks <- ranks[-1,]
      }
    }
    ranks <-ranks[order(ranks$Rank),] # re-orders it in terms of rank
    hpd_ranks <- ranks$Rank
  }

  # formatting
  if (length(hpd_ranks) == 1) { # if there is just one element
    concat_ranks <- hpd_ranks
  } else if (length(hpd_ranks) == 0) {
    concat_ranks <- "N/A"
  } else if (all(diff(as.numeric(as.character(hpd_ranks))) == 1)) { # hpd is an interval
    concat_ranks <- paste(hpd_ranks[1], hpd_ranks[[length(hpd_ranks)]], sep = "-")  # formats the ranks into an interval
  } else { # hpd is not an interval
    concat_ranks <- paste(hpd_ranks,collapse=',') # collapses the ranks into one string
  }
  hpd_vec <- list(concat_ranks,sum(ranks$Freq),hpd_ranks)
  return(hpd_vec)
}
