#' Identify and optionally trim redundant hierarchies
#'
#' @description
#' `find_redundancies` compares the outputs from `get_arrangements()`,
#' `get_partial_hiearchies()`, and `get_ranks_by_treatment()` to identify
#' redundant credible hierarchies.
#'
#' @param algo_1 a list of data frames outputted from `get_arrangements()`
#'   containing the credible hierarchies for ranked permutations, permutations,
#'   ranked combinations, and combinations.
#' @param algo_2 a data frame outputted from `get_partial_hierarchies()`
#'   containing credible partial hierarchies.
#' @param algo_3 a data frame outputted from `get_ranks_by_treatments()`
#'   containing credible HDRs.
#' @param n_trt a numeric value indicating the number of treatments in the
#'   network.
#' @param threshold a proportion between 0 and 1 for which a hierarchy must be
#'   observed in order to be credible.
#' @param type a numeric vector indicating what types of hierarchy redundancies
#'   should be trimmed. See details for more information. Default is 1:18.
#' @param trim_redundant a logical value indicating whether the redundant
#'   hierarchies should be trimmed from the output (TRUE) or not (FALSE, the
#'   default).
#'
#' @details
#' Redundancy types can be identified with numbers through 1 through 18 as follows:
#' \enumerate{
#'  \item Within ranked permutations (e.g., \eqn{(A, B)_1^2} is redundant because of \eqn{(A, B, C)_1^3}).
#'  \item Within permutations (e.g., \eqn{(A, B)} is redundant because of \eqn{(A, B, C)}).
#'  \item Within partial hierarchies (e.g., \eqn{A > B} is redundant because of \eqn{A > B > C}).
#'  \item Between permutations and partial hierarchies (e.g., \eqn{A > B} is redundant because of \eqn{(A, B)}).
#'  \item Between ranked permutations and permutations (e.g., \eqn{(A, B)} is redundant because of \eqn{(A, B)_1^2}).
#'  \item Between ranked combinations and combinations (e.g., \eqn{\{A, B\}} is redundant because of \eqn{\{A, B\}_1^2}).
#'  \item Between permutations and combinations (e.g., \eqn{\{A, B\}} is redundant because of \eqn{(A, B)}).
#'  \item Between ranked permutations and ranked combinations (e.g., \eqn{\{A, B\}_1^2} is redundant because of \eqn{(A, B)_1^2}).
#'  \item Between a ranked permutation and multiple single HDR ranks.
#'  \item Between a ranked combination and multiple HDRs.
#'  \item Within top/bottom ranked combinations.
#'  \item Between top/bottom single HDR ranks and a combination ranking from 1:(n_trt-1) or 2:n_trt.
#'  \item Between permutations ranking from 1:(n_trt-1) or 2:n_trt and a top/bottom single HDR rank.
#'  \item Between top/bottom single HDR ranks and a middle combination ranking 2:(n_trt-1).
#'  \item Between (top/bottom single HDR rank + tail-ranked permutation) and a middle ranked combination.
#'  \item Between tail-ranked permutations and a middle ranked combination.
#'  \item Between (ranked permutation + smaller ranked combination) and ranked combination.
#'  \item Between (single HDR rank + smaller ranked combination) and ranked combination.
#' }
#' We recommend these redundancies are identified in increasing order (1 through
#' 18) to ensure a complete assessment and output. Otherwise, a warning will be
#' outputted.
#'
#' @return A list of the credible ranked permutations, permutations, ranked
#' combinations, combinations, partial hierarchies, and HDRs along with their
#' redundancy status.
#'
#' @importFrom stats na.omit
#'
#' @export
#'
#' @examples
#' inputs <- prep_data(effects_matrix = dat_Thijs2008[, -1],
#'                     reference = "Placebo",
#'                     larger_better = FALSE)
#' algo1 <- get_arrangements(hierarchy_matrix = inputs$hierarchy_matrix,
#'                           threshold = 0.95)
#' algo2 <- get_partial_hierarchies(effects_matrix = inputs$effects_matrix,
#'                                  mid = 0,
#'                                  threshold = 0.95,
#'                                  larger_better = FALSE)
#' algo3 <- get_ranks_by_treatment(ranking_df = inputs$ranking_df,
#'                                 threshold = 0.95,
#'                                 print_plot = FALSE)
#' find_redundancies(algo1, algo2, algo3, n_trt = 5, threshold = 0.95)
find_redundancies <- function(algo_1, algo_2, algo_3,
                              n_trt, threshold, type = 1:18,
                              trim_redundant = FALSE) {

  # Verify algo_2 objects correspond to MID = 0
  if(gsub("Treatments at MID = ", "", colnames(algo_2)[1]) != 0){
    warning("Partial hierarchies constructed with a non-zero MID; partial hierarchies that are made redundant by permutations will not be identified.")
    if(2 %in% type) {
      type <- type[-which(type == 2)]  # do not search for partial hierarchies that are redundant because of permutations
    }
  }

  # Verify if inputted types are in increasing order
  if(any(type != 1:18)) {
    warning("`types` is not in increasing order from 1 to 18; some redundancies may be missed.")
  }

  tolerance <- .Machine$double.eps ^ 0.5
  threshold <- threshold - tolerance

  # Extract different hierarchy types
  ranked_perm <- algo_1[[1]]
  ranked_perm[, 2] <- stringr::str_remove_all(ranked_perm[, 2], "[()]")
  perm <- algo_1[[2]]
  perm[, 1] <- stringr::str_remove_all(perm[, 1], "[()]")
  ranked_comb <- algo_1[[3]]
  ranked_comb[, 2] <- stringr::str_remove_all(ranked_comb[, 2], "[{}]")
  comb <- algo_1[[4]]
  comb[, 1] <- stringr::str_remove_all(comb[, 1], "[{}]")
  phier <- algo_2
  hdr <- algo_3[[1]]
  hdr <- hdr[order(hdr$`HDR Rank(s)`, hdr$Treatment),]

  # Now find redundant hierarchies
  for(ind in type) {

    if(ind == 1) {

      # Within ranked permutations (e.g., $(A, B)_1^2$ is redundant because of $(A, B, C)_1^3$).

      if(nrow(ranked_perm) > 0) {

        if(!exists("Redundant", ranked_perm)) {
          ranked_perm$Redundant <- "FALSE"
        }

        # check if credible ranked permutations are redundant within
        ranked_perm_ind <- 1:nrow(ranked_perm)
        for(i in ranked_perm_ind){
          if(ranked_perm[i, "Redundant"] == "FALSE") {
            range_target <- str_split_1(ranked_perm[i, 1], "-")
            ranked_perm_target <- str_split_1(as.character(ranked_perm[i, 2]), ",")
            for(j in ranked_perm_ind[-i]) {
              ranked_perm_range_to_check <- str_split_1(ranked_perm[j, 1], "-")

              # if target permutation ranks are within range of permutation ranks to check
              if(range_target[1] >= ranked_perm_range_to_check[1] &
                 range_target[2] <= ranked_perm_range_to_check[2] ) {

                # find index of ranked_perm_string_to_check that would match positions of treatments to check
                upper <- as.numeric(range_target[2]) - as.numeric(ranked_perm_range_to_check[1]) + 1
                lower <- as.numeric(range_target[1]) - as.numeric(ranked_perm_range_to_check[1]) + 1

                # extract treatments in range of target ranked permutation
                ranked_perm_string_to_check <- str_split_1(as.character(ranked_perm[j, 2]), ",")[lower:upper]

                if(identical(ranked_perm_target, ranked_perm_string_to_check)) {
                  ranked_perm[i, "Redundant"] <- "TRUE"
                  break
                }

              }

            }
          }
        }

      }

    } else if(ind == 2) {

      # Within permutations (e.g., (A, B) is redundant because of (A, B, C)).

      if(nrow(perm) > 0) {

        if(!exists("Redundant", perm)) {
          perm$Redundant <- "FALSE"
        }

        # check if (remaining) credible permutations are redundant within
        perm_ind <- 1:nrow(perm)
        for(i in perm_ind){
          if(perm[i, "Redundant"] == "FALSE") {
            perm_target <- as.character(perm[i, 1])
            for(j in perm_ind[-i]) {
              perm_string <- as.character(perm[j, 1])
              if(str_detect(perm_string, perm_target)) {
                perm[i, "Redundant"] <- "TRUE"
                break
              }
            }
          }
        }

      }

    } else if(ind == 3) {

      # Within partial hierarchies (e.g., A > B is redundant because of A > B > C).

      if(nrow(phier) > 0) {

        if(!exists("Redundant", phier)) {
          phier$Redundant <- "FALSE"
        }

        # order by size of all credible partial hierarchies
        phier <- phier[order(phier$Size), ]

        # create a list of credible hierarchies by size
        phier_list <- split(phier, phier$Size)

        # unique sizes of credible partial hierarchies
        phier_sizes <- sort(unique(phier$Size))

        # check if credible partial hierarchies are redundant within
        if(length(phier_sizes) > 1) {
          # assess whether smaller credible partial hierarchies are within larger credible partial hierarchies
          for(i in phier_sizes[-which.max(phier_sizes)]) {
            # extract credible partial hierarchies of size i
            current_phier_size <- phier_list[[as.character(i)]]
            # list of larger hierarchies to check against
            phier_larger_list <- phier_list[as.character(phier_sizes[which(phier_sizes > i)])]
            for(j in 1:nrow(current_phier_size)) {
              if(current_phier_size[j, "Redundant"] == "FALSE") {
                # extract jth partial hierarchy in current_phier_size
                phier_target <- stringr::str_split_1(as.character(current_phier_size[j, 1]), " > ")
                # check to see if jth partial hierarchy is in larger partial hierarchies
                redundant_check <- is_phier_redundant_within_phier(phier_target,
                                                                   phier_larger_list)
                if(!is.null(redundant_check)){
                  current_phier_size[j, "Redundant"] <- redundant_check
                }
              }
            }
            phier_list[[as.character(i)]] <- current_phier_size
          }
        }
        phier <- do.call(rbind, phier_list)

      }

    } else if(ind == 4) {

      # Between permutations and partial hierarchies (e.g., A > B is redundant
      # because of (A, B)).

      if(nrow(phier) > 0) {

        if(!exists("Redundant", phier)) {
          phier$Redundant <- "FALSE"
        }

        if(nrow(perm) > 0) {

          # check if (remaining) credible partial hierarchies are redundant
          # because of credible permutations
          for (i in 1:nrow(phier)) {
            if(phier[i, "Redundant"] == "FALSE") {
              phier_target <- str_split_1(as.character(phier[i, 1]), " > ")
              for(j in 1:nrow(perm)) {
                perm_string_to_check <- str_split_1(as.character(perm[j, 1]), ",")
                # finds position of treatments in phier_target within perm_string_to_check
                positions <- match(phier_target, perm_string_to_check)
                if (!anyNA(positions) && all(diff(positions) >= 0)) { # checks if position[j+1] >= position[j]
                  phier[i, "Redundant"] <- "TRUE"
                  break
                }
              }
            }
          }

        }
      }

    } else if(ind == 5) {

      # Between ranked permutations and permutations (e.g., $(A, B)$ is redundant because of $(A, B)_1^2$).

      if(nrow(perm) > 0) {

        if(!exists("Redundant", perm)) {
          perm$Redundant <- "FALSE"
        }

        if(nrow(ranked_perm) > 0) {

          # check if (remaining) credible permutations are redundant because
          # of credible ranked permutations
          for(i in 1:nrow(perm)) {
            perm_target <- str_split_1(perm[i, 1], ",")
            if(perm[i, "Redundant"] == "FALSE") {
              for(j in 1:nrow(ranked_perm)) {
                if(setequal(perm_target, str_split_1(ranked_perm[j, 2], ","))) {
                  perm[i, "Redundant"] <- "TRUE"
                  break
                }
              }
            }
          }

        }
      }

    } else if(ind == 6) {

      # Between ranked combinations and combinations (e.g., ${A, B}$
      # is redundant because of \eqn{\{A, B\}_1^2}).

      if(nrow(comb) > 0) {

        if(!exists("Redundant", comb)) {
          comb$Redundant <- "FALSE"
        }

        if(nrow(ranked_comb) > 0) {

          # check if (remaining) credible combinations are redundant because
          # of credible ranked combinations
          for(i in 1:nrow(comb)) {
            comb_target <- str_split_1(comb[i, 1], ",")
            if(comb[i, "Redundant"] == "FALSE") {
              for(j in 1:nrow(ranked_comb)) {
                if(setequal(comb_target, str_split_1(ranked_comb[j, 2], ","))) {
                  comb[i, "Redundant"] <- "TRUE"
                  break
                }
              }
            }
          }

        }
      }

    } else if(ind == 7) {

      # Between permutations and combinations (e.g., ${A, B}$ is redundant
      # because of $(A, B)$).

      if(nrow(comb) > 0) {

        if(!exists("Redundant", comb)) {
          comb$Redundant <- "FALSE"
        }

        if(nrow(perm) > 0) {

          # check if (remaining) credible combinations are redundant because
          # of credible permutations
          for(i in 1:nrow(comb)) {
            comb_target <- str_split_1(comb[i, 1], ",")
            if(comb[i, "Redundant"] == "FALSE") {
              for(j in 1:nrow(perm)) {
                if(setequal(comb_target, str_split_1(perm[j, 1], ","))) {
                  comb[i, "Redundant"] <- "TRUE"
                  break
                }
              }
            }
          }

        }
      }

    } else if(ind == 8) {

      # Between ranked permutations and ranked combinations
      # (e.g., \eqn{\{A, B\}_1^2} is redundant because of \eqn{(A, B)_1^2}).

      if(nrow(ranked_comb) > 0) {

        if(!exists("Redundant", ranked_comb)) {
          ranked_comb$Redundant <- "FALSE"
        }

        if(nrow(ranked_perm) > 0) {

          # check if (remaining) credible ranked combinations are redundant
          # because of credible ranked permutations
          for(i in 1:nrow(ranked_comb)) {
            range_target <- as.character(ranked_comb[i, 1])
            ranked_comb_target <- str_split_1(ranked_comb[i, 2], ",")
            if(ranked_comb[i, "Redundant"] == "FALSE") {
              for(j in 1:nrow(ranked_perm)) {
                if(range_target == as.character(ranked_perm[j, 1]) &&
                   setequal(ranked_comb_target,
                            str_split_1(ranked_perm[j, 2], ","))) {
                  ranked_comb[i, "Redundant"] <- "TRUE"
                  break
                }
              }
            }
          }

        }
      }

    } else if(ind == 9) {

      # Between a ranked permutation and multiple single HDR ranks.

      # Find HDRs with single rank
      HDR_single_ind <- which(nchar(hdr[, "HDR Rank(s)"]) == 1)

      if(nrow(ranked_perm) > 0 & length(HDR_single_ind) > 0) {

        if(!exists("Redundant", hdr)) {
          hdr$Redundant <- "FALSE"
        }

        if(!exists("Redundant", ranked_perm)) {
          ranked_perm$Redundant <- "FALSE"
        }

        # note indices of non-redundant ranked permutations
        ranked_perm_target <- which( ranked_perm[, "Redundant"] == "FALSE" )

        for(i in ranked_perm_target) {
          perm_ranks <- stringr::str_split_1(ranked_perm[i, "Range"], "-")
          perm_ranks <- as.character(perm_ranks[1]:perm_ranks[2])
          perm_set <- stringr::str_split_1(ranked_perm[i, "Ranked Permutations"],
                                           ",")
          if(all(perm_ranks %in% hdr[HDR_single_ind, "HDR Rank(s)"],
                 stats::na.omit(hdr[match(perm_ranks,
                                          hdr[, "HDR Rank(s)"]),
                                    "Treatment"]) == perm_set)) {
            hdr[stats::na.omit(match(perm_ranks,
                                     hdr[, "HDR Rank(s)"])), "Redundant"] <- "TRUE"
          }
        }
      }

    } else if(ind == 10) {

      # Between a ranked combination and multiple HDRs.

      if(!exists("Redundant", hdr)) {
        hdr$Redundant <- "FALSE"
      }

      # Which treatments have the same HDRs
      dup_ind <- which(duplicated(hdr[, "HDR Rank(s)"]) == TRUE)

      # If any treatments have same set of HDRs...
      if(length(dup_ind) > 0) {

        # Which ranked combinations' ranking ranges == the HDRs?
        rk_comb_ind <- which(ranked_comb[, "Range"] %in% hdr[dup_ind, "HDR Rank(s)"] == TRUE)

        # If any of the ranked combinations' ranking range == the HDRs ...
        if(length(rk_comb_ind) > 0) {

          # Do the ranked combinations include the treatments with the same HDRs?
          for(i in rk_comb_ind) {

            # Treatments involved in ranked combination
            rk_comb_trt <- stringr::str_split_1(ranked_comb[i, "Ranked Combinations"], ",")
            # Ranks involved in ranked combination
            rk_comb_range <- ranked_comb[i, "Range"]

            for(j in dup_ind) {
              # Index of HDRs to assess redundancy status
              hdr_to_check <- which(hdr[, "HDR Rank(s)"] == hdr[j, "HDR Rank(s)"])
              # Range of HDR being checked
              hdr_range <- unique(hdr[hdr_to_check, "HDR Rank(s)"])
              # Treatments with HDR being checked
              hdr_trt <-  hdr[hdr_to_check, "Treatment"]

              # If the treatments with same HDRs are exactly the same treatments
              # in a ranked combination with a ranking range == the HDRs...
              if(all(hdr_trt %in% rk_comb_trt,
                     length(hdr_trt) == length(rk_comb_trt),
                     hdr_range == rk_comb_range)) {
                hdr[hdr_to_check, "Redundant"] <- "TRUE"
              }
            }

          }

        }

      }

    } else if(ind == 11) {

      # Within top/bottom ranked combinations.

      if(nrow(ranked_comb) > 0) {

        if(!exists("Redundant", ranked_comb)) {
          warning("For completeness, ranked combinations should ideally be
                  assessed for within redundancies after comparison with
                  other hierarchy question types.")
          ranked_comb$Redundant <- "FALSE"
        }

        # Note indices of redundant combinations ranking 1:j or (j+1):n
        ranked_comb_target <- which((stringr::str_split(ranked_comb[, "Range"],
                                                        "-", simplify = TRUE)[, 1] == "1" |
                                       stringr::str_split(ranked_comb[, "Range"],
                                                          "-", simplify = TRUE)[, 2] == as.character(n_trt)) &
                                      ranked_comb[, "Redundant"] == "TRUE")

        # Check if any of the redundant ranked combinations complete the
        # (remaining) credible ranked combinations
        for(i in 1:nrow(ranked_comb)) {
          if(ranked_comb[i, "Redundant"] == "FALSE") {
            for(j in ranked_comb_target) {

              # Check if (maximum rank of top-ranked combination  - 1) = (minimum rank of bottom ranked combination)
              check_ranked_comb_start <- stringr::str_split_1(ranked_comb[i, "Range"], "-")[1] == "1"
              check_ranked_comb_target_end <- stringr::str_split_1(ranked_comb[j, "Range"], "-")[2] == as.character(n_trt)
              check_ranked_comb_end_target_start <- as.numeric(stringr::str_split_1(ranked_comb[i, "Range"], "-")[2]) == (as.numeric(stringr::str_split_1(ranked_comb[j, "Range"], "-")[1]) - 1)
              check_unique_trt <- all(!(stringr::str_split_1(ranked_comb[i, "Ranked Combinations"], ",") %in%
                                          stringr::str_split_1(ranked_comb[j, "Ranked Combinations"], ",")))

              if(all(check_ranked_comb_start,
                     check_ranked_comb_target_end,
                     check_ranked_comb_end_target_start,
                     check_unique_trt)) {
                ranked_comb[i, "Redundant"] <- "TRUE"
              } else if(all(stringr::str_split_1(ranked_comb[j, "Range"], "-")[1] == "1",
                            stringr::str_split_1(ranked_comb[i, "Range"], "-")[2] == as.character(n_trt),
                            as.numeric(stringr::str_split_1(ranked_comb[j, "Range"], "-")[2]) == (as.numeric(stringr::str_split_1(ranked_comb[i, "Range"], "-")[1]) - 1),
                            all(!(stringr::str_split_1(ranked_comb[i, "Ranked Combinations"], ",") %in%
                                  stringr::str_split_1(ranked_comb[j, "Ranked Combinations"], ",")))
              )) {
                ranked_comb[i, "Redundant"] <- "TRUE"
              }

            }
          }
        }

      }

    } else if(ind == 12) {

      # Between top/bottom single HDR ranks and a
      # combination ranking from 1:(n_trt-1) or 2:n_trt.

      # Find HDR consisting of single rank 1
      HDR_single_1_ind <- which(hdr[, "HDR Rank(s)"] == "1")
      # Find HDR consisting of single rank n_trt
      HDR_single_n_ind <- which(hdr[, "HDR Rank(s)"] == as.character(n_trt))

      if(nrow(ranked_comb) > 0 & (length(HDR_single_1_ind) > 0 | length(HDR_single_n_ind) > 0)) {

        if(!exists("Redundant", ranked_comb)) {
          ranked_comb$Redundant <- "FALSE"
        }

        # note indices of combinations ranking 1:(n_trt - 1)
        ranked_comb_to_check_1 <- which(ranked_comb$Range == paste0("1-", n_trt - 1))

        # note indices of combinations ranking 2:n_trt
        ranked_comb_to_check_n <- which(ranked_comb$Range == paste0("2-", n_trt))

        # Check if the top single HDR rank completes the ranked combination
        if(length(ranked_comb_to_check_n) > 0 & length(HDR_single_1_ind) > 0) {
          for(i in ranked_comb_to_check_n) {
            if(!(hdr[HDR_single_1_ind, "Treatment"] %in%
                 stringr::str_split_1(ranked_comb[i, "Ranked Combinations"], ","))) {
              ranked_comb[i, "Redundant"] <- "TRUE"
            }
          }
        }

        # Check if the bottom single HDR rank completes the ranked combination
        if(length(ranked_comb_to_check_1) > 0 & length(HDR_single_n_ind) > 0) {
          for(i in ranked_comb_to_check_1) {
            if(!(hdr[HDR_single_n_ind, "Treatment"] %in%
                 stringr::str_split_1(ranked_comb[i, "Ranked Combinations"], ","))) {
              ranked_comb[i, "Redundant"] <- "TRUE"
            }
          }
        }

      }

    } else if(ind == 13) {

      # Between permutations ranking from 1:(n_trt-1) or
      # 2:n_trt and a top/bottom single HDR rank.

      if(!exists("Redundant", hdr)) {
        hdr$Redundant <- "FALSE"
      }

      # Find HDR consisting of single rank 1
      HDR_single_1_ind <- which(hdr[, "HDR Rank(s)"] == "1")
      # Find HDR consisting of single rank n_trt
      HDR_single_n_ind <- which(hdr[, "HDR Rank(s)"] == as.character(n_trt))

      if(nrow(ranked_perm) > 0 & (length(HDR_single_1_ind) > 0 | length(HDR_single_n_ind) > 0)) {

        # note index of permutations ranking 1:(n_trt - 1)
        ranked_perm_to_check_1 <- which(ranked_perm$Range == paste0("1-", n_trt - 1) &
                                          ranked_perm$Redundant == "FALSE")

        # note indices of permutations ranking 2:n_trt
        ranked_perm_to_check_n <- which(ranked_perm$Range == paste0("2-", n_trt) &
                                          ranked_perm$Redundant == "FALSE")

        # Check if the top single HDR rank completes the ranked permutation
        if(length(ranked_perm_to_check_n) > 0 & length(HDR_single_1_ind) > 0) {
          for(i in ranked_perm_to_check_n) {
            if(!(hdr[HDR_single_1_ind, "Treatment"] %in%
                 stringr::str_split_1(ranked_perm[i, "Ranked Permutations"], ","))) {
              hdr[HDR_single_1_ind, "Redundant"] <- "TRUE"
            }
          }
        }

        # Check if the bottom single HDR rank completes the ranked permutation
        if(length(ranked_perm_to_check_1) > 0 & length(HDR_single_n_ind) > 0) {
          for(i in ranked_perm_to_check_1) {
            if(!(hdr[HDR_single_n_ind, "Treatment"] %in%
                 stringr::str_split_1(ranked_perm[i, "Ranked Permutations"], ","))) {
              hdr[HDR_single_n_ind, "Redundant"] <- "TRUE"
            }
          }
        }

      }

    } else if(ind == 14) {

      # Between top/bottom single HDR ranks and a middle combination ranking 2:(n_trt-1).

      # Find HDR consisting of single rank 1
      HDR_single_1_ind <- which(hdr[, "HDR Rank(s)"] == "1" & hdr[, "Redundant"] == "TRUE")
      # Find HDR consisting of single rank n_trt
      HDR_single_n_ind <- which(hdr[, "HDR Rank(s)"] == as.character(n_trt) & hdr[, "Redundant"] == "TRUE")

      if(nrow(ranked_comb) > 0 & length(HDR_single_1_ind) > 0 & length(HDR_single_n_ind) > 0) {

        if(!exists("Redundant", ranked_comb)) {
          ranked_comb$Redundant <- "FALSE"
        }

        # note indices of combinations ranking 2:(n_trt - 1)
        ranked_comb_to_check <- which(ranked_comb$Range == paste0("2-", n_trt - 1))

        if(length(ranked_comb_to_check) > 0) {
          # check if any of the individual treatments with credible ranks complete
          # the ranked combination
          for(i in ranked_comb_to_check) {
            if(all(!(hdr[HDR_single_1_ind, "Treatment"] %in%
                     stringr::str_split_1(ranked_comb[i, "Ranked Combinations"], ",")),
                   !(hdr[HDR_single_n_ind, "Treatment"] %in%
                     stringr::str_split_1(ranked_comb[i, "Ranked Combinations"], ",")),
                   min(as.numeric(stringr::str_split_1(ranked_comb[i, "Range"], "-"))) - 1 == 1,
                   max(as.numeric(stringr::str_split_1(ranked_comb[i, "Range"], "-"))) + 1 == n_trt)) {

              if( (hdr[HDR_single_1_ind, "Sum of pi_hat"] + hdr[HDR_single_n_ind, "Sum of pi_hat"] - 1) >= threshold ) {
                ranked_comb[i, "Redundant"] <- "TRUE"
                break
              } else {

                hdr_string_1 <- paste0(hdr[HDR_single_1_ind, "Treatment"], "_{", 1, "}")
                hdr_string_n <- paste0(hdr[HDR_single_1_ind, "Treatment"], "_{", n_trt, "}")

                ranked_comb[i, "Redundant"] <- paste0("Check P(", hdr_string_1,
                                                      " AND ", hdr_string_n, ")")
                break
              }

            }
          }
        }

      }

    } else if(ind == 15) {

      # Between (top/bottom single HDR rank + tail-ranked permutation) and a middle ranked combination.

      # Find HDR consisting of single rank 1
      HDR_single_1_ind <- which(hdr[, "HDR Rank(s)"] == "1" & hdr[, "Redundant"] == "TRUE")
      # Find HDR consisting of single rank n_trt
      HDR_single_n_ind <- which(hdr[, "HDR Rank(s)"] == as.character(n_trt) & hdr[, "Redundant"] == "TRUE")

      if(nrow(ranked_comb) > 0 & nrow(ranked_perm) > 0 & (length(HDR_single_1_ind) > 0 | length(HDR_single_n_ind) > 0)) {

        if(!exists("Redundant", ranked_comb)) {
          ranked_comb$Redundant <- "FALSE"
        }

        # Note indices of combinations with ranks starting at 2
        ranked_comb_to_check_2 <- which( stringr::str_split(ranked_comb[, "Range"], "-", simplify = TRUE)[, 1] == "2" &
                                           as.numeric(stringr::str_split(ranked_comb[, "Range"], "-", simplify = TRUE)[, 2]) <= (n_trt - 1) )

        # Note indices of combinations with ranks ending at (n-1)
        ranked_comb_to_check_n1 <- which( as.numeric(stringr::str_split(ranked_comb[, "Range"], "-", simplify = TRUE)[, 1]) >= 2 &
                                            stringr::str_split(ranked_comb[, "Range"], "-", simplify = TRUE)[, 2] == as.character(n_trt - 1) )

        # Note indices of ranked permutations with ranks starting at 1 that are redundant
        ranked_perm_target_1 <- which( stringr::str_split(ranked_perm[, "Range"], "-", simplify = TRUE)[, 1] == "1" & ranked_perm[, "Redundant"] == "TRUE")

        # Note indices of ranked permutations with ranks ending at n that are redundant
        ranked_perm_target_n <- which( stringr::str_split(ranked_perm[, "Range"], "-", simplify = TRUE)[, 2] == as.character(n_trt) & ranked_perm[, "Redundant"] == "TRUE")

        # Check if any combinations starting with rank 2 are redundant because of a top single HDR rank and a permutation ending with rank n
        if(length(ranked_comb_to_check_2) > 0 & length(HDR_single_1_ind) > 0 & length(ranked_perm_target_n) > 0) {
          for(i in ranked_comb_to_check_2) {
            for(j in ranked_perm_target_n) {
              if(all(!(hdr[HDR_single_1_ind, "Treatment"] %in%
                       stringr::str_split_1(ranked_comb[i, "Ranked Combinations"], ",")),
                     all(!(stringr::str_split_1(ranked_comb[i, "Ranked Combinations"], ",") %in%
                           stringr::str_split_1(ranked_perm[j, "Ranked Permutations"], ","))),
                     !(hdr[HDR_single_1_ind, "Treatment"] %in%
                       stringr::str_split_1(ranked_perm[j, "Ranked Permutations"], ",")),
                     max(as.numeric(stringr::str_split_1(ranked_comb[i, "Range"], "-"))) ==
                     min(as.numeric(stringr::str_split_1(ranked_perm[j, "Range"], "-"))) - 1)) {

                if( (hdr[HDR_single_1_ind, "Sum of pi_hat"] + ranked_perm[j, "pi_hat"] - 1) >= threshold ) {
                  ranked_comb[i, "Redundant"] <- "TRUE"
                  break
                } else {

                  hdr_single_string <- paste0(hdr[HDR_single_1_ind,
                                                  "Treatment"],
                                              "_{", 1, "}")
                  ranked_perm_string <- paste0("(", ranked_perm[j, "Ranked Permutations"], ")_",
                                               stringr::str_split_1(ranked_perm[j, "Range"], "-")[1],
                                               "^",
                                               stringr::str_split_1(ranked_perm[j, "Range"], "-")[2])


                  ranked_comb[i, "Redundant"] <- paste0("Check P(",
                                                        hdr_single_string,
                                                        " AND ",
                                                        ranked_perm_string, ")")
                  break
                }

              }
            }
          }
        }

        # Check if any combinations ending with rank n-1 are redundant because of a bottom single HDR rank and a permutation starting with rank 1
        if(length(ranked_comb_to_check_n1) > 0 & length(HDR_single_n_ind) > 0 & length(ranked_perm_target_1) > 0) {
          for(i in ranked_comb_to_check_n1) {
            for(j in ranked_perm_target_1) {
              if(all(!(hdr[HDR_single_n_ind, "Treatment"] %in%
                       stringr::str_split_1(ranked_comb[i, "Ranked Combinations"], ",")),
                     all(!(stringr::str_split_1(ranked_comb[i, "Ranked Combinations"], ",") %in%
                           stringr::str_split_1(ranked_perm[j, "Ranked Permutations"], ","))),
                     !(hdr[HDR_single_n_ind, "Treatment"] %in%
                       stringr::str_split_1(ranked_perm[j, "Ranked Permutations"], ",")),
                     min(as.numeric(stringr::str_split_1(ranked_comb[i, "Range"], "-"))) ==
                     max(as.numeric(stringr::str_split_1(ranked_perm[j, "Range"], "-"))) + 1)) {

                if( (hdr[HDR_single_n_ind, "Sum of pi_hat"] + ranked_perm[j, "pi_hat"] - 1) >= threshold ) {
                  ranked_comb[i, "Redundant"] <- "TRUE"
                  break
                } else {

                  hdr_single_string <- paste0(hdr[HDR_single_n_ind, "Treatment"], "_{", as.character(n_trt), "}")
                  ranked_perm_string <- paste0("(", ranked_perm[j, "Ranked Permutations"], ")_",
                                               stringr::str_split_1(ranked_perm[j, "Range"], "-")[1],
                                               "^",
                                               stringr::str_split_1(ranked_perm[j, "Range"], "-")[2])


                  ranked_comb[i, "Redundant"] <- paste0("Check P(", hdr_single_string,
                                                        " AND ", ranked_perm_string, ")")
                  break
                }

              }
            }
          }
        }

      }

    } else if(ind == 16) {

      # Between tail-ranked permutations and a middle ranked combination.

      if(nrow(ranked_comb) > 0) {

        if(!exists("Redundant", ranked_comb)) {
          ranked_comb$Redundant <- "FALSE"
        }

        # Index of ranked combinations to check that have not yet been flagged as redundant
        rk_comb_to_check <- which(ranked_comb[, "Redundant"] == "FALSE")

        # If there are any ranked combinations and it is possible to check against tail ranked permutations (i.e., at least two ranked permutations)...
        if(length(rk_comb_to_check) > 0 & nrow(ranked_perm) > 1) {

          # Index of redundant ranked permutations
          red_ranked_perm <- ranked_perm[which(ranked_perm[, "Redundant"] == TRUE), ]

          # If it is possible to check against REDUNDANT tail ranked permutations (i.e., at least two ranked permutations)...
          if(length(red_ranked_perm) > 1) {

            # Minimum rank of ranked permutations
            ranked_perm_min <- sapply(red_ranked_perm[, "Range"], stringr::str_split_1, pattern = "-", simplify = TRUE)[1, ]
            # Maximum rank of ranked permutations
            ranked_perm_max <- sapply(red_ranked_perm[, "Range"], stringr::str_split_1, pattern = "-", simplify = TRUE)[2, ]

            # Index of ranked permutations which start at rank 1
            ranked_perm_1 <- which(ranked_perm_min == 1)
            # Index of ranked permutations which end at rank n
            ranked_perm_n <- which(ranked_perm_max == n_trt)

            # If there are any ranked combinations to check against tail ranks...
            if(length(ranked_perm_1) > 0 & length(ranked_perm_n) > 0) {

              for(i in rk_comb_to_check) {
                for(j in ranked_perm_1) {
                  for(k in ranked_perm_n) {

                    # Range of ranked combination to check
                    ranked_comb_range <- ranked_comb[i, "Range"]
                    # Treatments of ranked combination to check
                    ranked_comb_trt <-  stringr::str_split_1(ranked_comb[i, "Ranked Combinations"], ",")

                    # Range of top ranked permutation to check against
                    ranked_perm_1_range <- red_ranked_perm[j, "Range"]
                    # Treatments of top ranked permutation to check against
                    ranked_perm_1_trt <-  stringr::str_split_1(red_ranked_perm[j, "Ranked Permutations"], ",")

                    # Range of bottom ranked permutation to check against
                    ranked_perm_n_range <- red_ranked_perm[k, "Range"]
                    # Treatments of bottom ranked permutation to check against
                    ranked_perm_n_trt <- stringr::str_split_1(red_ranked_perm[k, "Ranked Permutations"], ",")

                    # If the ranked combination's ranking range is sandwiched by that of the ranked permutations,
                    # the treatments in the ranked combination are not in either ranked permutation,
                    # and there are no duplicate treatments in the ranked permutations...
                    if(all(ranked_comb_range == paste0((as.numeric(stringr::str_split_1(ranked_perm_1_range, "-")[2]) + 1),
                                                       "-",
                                                       (as.numeric(stringr::str_split_1(ranked_perm_n_range, "-")[1]) - 1)),
                           !(ranked_comb_trt %in% c(ranked_perm_1_trt, ranked_perm_n_trt)),
                           !(ranked_perm_n_trt %in% ranked_perm_n_trt))) {

                      # Verify joint empirical probability of tail-ranked permutations is credible
                      if( (red_ranked_perm[j, "pi_hat"] + red_ranked_perm[k, "pi_hat"] - 1) >= threshold ) {
                        ranked_comb[i, "Redundant"] <- "TRUE"
                        break
                      } else {

                        ranked_perm_1_string <- paste0("(", red_ranked_perm[j, "Ranked Permutations"], ")_",
                                                       stringr::str_split_1(red_ranked_perm[j, "Range"], "-")[1],
                                                       "^",
                                                       stringr::str_split_1(red_ranked_perm[j, "Range"], "-")[2])
                        ranked_perm_n_string <- paste0("(", red_ranked_perm[k, "Ranked Permutations"], ")_",
                                                       stringr::str_split_1(red_ranked_perm[k, "Range"], "-")[1],
                                                       "^",
                                                       stringr::str_split_1(red_ranked_perm[k, "Range"], "-")[2])

                        ranked_comb[i, "Redundant"] <- paste0("Check P(", ranked_perm_1_string,
                                                              " AND ", ranked_perm_n_string, ")")
                        break
                      }

                    }

                  }
                }
              }

            }

          }

        }

      }

    } else if(ind == 17) {

      # Between (ranked permutation + smaller ranked combination) and ranked combination.

      if(nrow(ranked_comb) > 0) {

        if(!exists("Redundant", ranked_comb)) {
          ranked_comb$Redundant <- "FALSE"
        }

        # Index of ranked combinations to check that have not yet been flagged as redundant
        rk_comb_to_check <- which(ranked_comb[, "Redundant"] == "FALSE")

        # If there are any ranked combinations and it is possible to check against ranked permutations...
        if(length(rk_comb_to_check) > 0 & nrow(ranked_perm) > 0) {

          for(i in rk_comb_to_check) {

            # Range of ranked combination to check
            ranked_comb_range <- ranked_comb[i, "Range"]
            # Treatments of ranked combination to check
            ranked_comb_trt <- stringr::str_split_1(ranked_comb[i, "Ranked Combinations"], ",")
            # Size of ranked combination to check
            ranked_comb_size <- ranked_comb[i, "Size"]

            # Find ranked combinations of smaller size
            sm_ranked_comb_ind <- which(ranked_comb[, "Size"] < ranked_comb_size)

            if(length(sm_ranked_comb_ind) > 0) {

              for(j in 1:nrow(ranked_perm)) {
                # Range of ranked permutation to check against
                ranked_perm_range <- ranked_perm[j, "Range"]
                # Treatments of ranked permutation to check against
                ranked_perm_trt <- stringr::str_split_1(ranked_perm[j, "Ranked Permutations"], ",")

                for(k in sm_ranked_comb_ind) {

                  # Range of smaller ranked combination to check against
                  sm_ranked_comb_range <- ranked_comb[k, "Range"]
                  # Treatments of ranked permutation to check against
                  sm_ranked_comb_trt <- stringr::str_split_1(ranked_comb[k, "Ranked Combinations"], ",")

                  # Are ranking ranges of ranked permutation and smaller ranked combination consecutive?
                  perm_sm_comb_range_check <- as.numeric(stringr::str_split_1(ranked_perm_range, "-")[2]) + 1 ==
                    as.numeric(stringr::str_split_1(sm_ranked_comb_range, "-")[1])
                  sm_comb_perm_range_check <- as.numeric(stringr::str_split_1(sm_ranked_comb_range, "-")[2]) + 1 ==
                    as.numeric(stringr::str_split_1(ranked_perm_range, "-")[1])

                  # Denote combined ranking range as appropriate
                  if(perm_sm_comb_range_check) {
                    range_check <- paste0(stringr::str_split_1(ranked_perm_range, "-")[1],
                                          "-",
                                          stringr::str_split_1(sm_ranked_comb_range, "-")[2])
                  } else if(sm_comb_perm_range_check) {
                    range_check <- paste0(stringr::str_split_1(sm_ranked_comb_range, "-")[1],
                                          "-",
                                          stringr::str_split_1(ranked_perm_range, "-")[2])
                  } else {
                    range_check <- "FALSE"
                  }

                  # If the ranked combination's ranking range is the same as the consecutive range of the ranked perm and smaller ranked comb,
                  # and the treatments are the exact same...
                  if(range_check != "FALSE") {
                    if(all(ranked_comb_range == range_check,
                           setequal(c(ranked_perm_trt, sm_ranked_comb_trt), ranked_comb_trt))) {

                      # Verify joint empirical probability of ranked permutation + smaller ranked combination is credible
                      if( (ranked_perm[j, "pi_hat"] + ranked_comb[k, "pi_hat"] - 1) >= threshold ) {
                        ranked_comb[i, "Redundant"] <- "TRUE"
                        break
                      } else {

                        ranked_perm_string <- paste0("(", ranked_perm[j, "Ranked Permutations"], ")_",
                                                     stringr::str_split_1(ranked_perm[j, "Range"], "-")[1],
                                                     "^",
                                                     stringr::str_split_1(ranked_perm[j, "Range"], "-")[2])
                        sm_ranked_comb_string <- paste0("(", ranked_comb[k, "Ranked Combinations"], ")_",
                                                        stringr::str_split_1(ranked_comb[k, "Range"], "-")[1],
                                                        "^",
                                                        stringr::str_split_1(ranked_comb[k, "Range"], "-")[2])

                        ranked_comb[i, "Redundant"] <- paste0("Check P(", ranked_perm_string,
                                                              " AND ", sm_ranked_comb_string, ")")
                        break
                      }

                    }
                  }

                }

              }

            }

          }

        }

      }

    } else if(ind == 18) {

      # Between (single HDR rank + smaller ranked combination) and ranked combination.

      if(nrow(ranked_comb) > 0) {

        if(!exists("Redundant", ranked_comb)) {
          ranked_comb$Redundant <- "FALSE"
        }

        # Index of ranked combinations to check that have not yet been flagged as redundant
        rk_comb_to_check <- which(ranked_comb[, "Redundant"] == "FALSE")

        # If there are any ranked combinations and it is possible to check against ranked permutations...
        if(length(rk_comb_to_check) > 0 & nrow(ranked_perm) > 0) {

          for(i in rk_comb_to_check) {

            # Range of ranked combination to check
            ranked_comb_range <- ranked_comb[i, "Range"]
            # Treatments of ranked combination to check
            ranked_comb_trt <- stringr::str_split_1(ranked_comb[i, "Ranked Combinations"], ",")
            # Size of ranked combination to check
            ranked_comb_size <- ranked_comb[i, "Size"]

            # Find ranked combinations of smaller size
            sm_ranked_comb_ind <- which(ranked_comb[, "Size"] < ranked_comb_size)

            # Find HDRs with single rank
            HDR_single_ind <- which(nchar(hdr[, "HDR Rank(s)"]) == 1)

            if(length(sm_ranked_comb_ind) > 0 & length(HDR_single_ind) > 0) {

              for(j in HDR_single_ind) {

                # Single rank in HDR to check against
                HDR_single_rk <- hdr[j, "HDR Rank(s)"]
                # Treatment with an HDR containing single rank to check against
                HDR_single_trt <- hdr[j, "Treatment"]

                for(k in sm_ranked_comb_ind) {
                  # Range of smaller ranked combination to check against
                  sm_ranked_comb_range <- ranked_comb[k, "Range"]
                  # Treatments of ranked permutation to check against
                  sm_ranked_comb_trt <- stringr::str_split_1(ranked_comb[k, "Ranked Combinations"], ",")

                  # Are ranking ranges of HDR and smaller ranked combination consecutive?
                  HDR_sm_comb_range_check <- as.numeric(HDR_single_rk) + 1 ==
                    as.numeric(stringr::str_split_1(sm_ranked_comb_range, "-")[1])
                  sm_comb_HDR_range_check <- as.numeric(stringr::str_split_1(sm_ranked_comb_range, "-")[2]) + 1 ==
                    as.numeric(HDR_single_rk)

                  # Denote combined ranking range as appropriate
                  if(HDR_sm_comb_range_check) {
                    range_check <- paste0(HDR_single_rk,
                                          "-",
                                          stringr::str_split_1(sm_ranked_comb_range, "-")[2])
                  } else if(sm_comb_HDR_range_check) {
                    range_check <- paste0(stringr::str_split_1(sm_ranked_comb_range, "-")[1],
                                          "-",
                                          HDR_single_rk)
                  } else {
                    range_check <- "FALSE"
                  }

                  # If the ranked combination's ranking range is the same as the consecutive range of the ranked perm and HDR with single rank,
                  # and the treatments are the exact same...
                  if(range_check != "FALSE") {
                    if(all(ranked_comb_range == range_check,
                           setequal(c(HDR_single_trt, sm_ranked_comb_trt), ranked_comb_trt))) {

                      # Verify joint empirical probability of single HDR rank + smaller ranked combination is credible
                      if( (hdr[j, "Sum of pi_hat"] + ranked_comb[k, "pi_hat"] - 1) >= threshold ) {
                        ranked_comb[i, "Redundant"] <- "TRUE"
                        break
                      } else {

                        hdr_single_string <- paste0(HDR_single_trt, "_{", HDR_single_rk, "}")
                        sm_ranked_comb_string <- paste0("(", ranked_comb[k, "Ranked Combinations"], ")_",
                                                        stringr::str_split_1(ranked_comb[k, "Range"], "-")[1],
                                                        "^",
                                                        stringr::str_split_1(ranked_comb[k, "Range"], "-")[2])

                        ranked_comb[i, "Redundant"] <- paste0("Check P(", hdr_single_string,
                                                              " AND ", sm_ranked_comb_string, ")")
                        break
                      }

                    }
                  }

                }

              }

            }

          }

        }

      }

    } else {
      stop(paste0("Invalid type = ", ind, "; should be a numeric value between 1 and 18."))
    }
  }

  if(trim_redundant == "TRUE"){
    ranked_perm <- ranked_perm[ranked_perm$Redundant == "FALSE", ]
    perm <- perm[perm$Redundant == "FALSE", ]
    ranked_comb <- ranked_comb[ranked_comb$Redundant == "FALSE", ]
    comb <- comb[comb$Redundant == "FALSE", ]
    phier <- phier[phier$Redundant == "FALSE", ]
    hdr <- hdr[hdr$Redundant == "FALSE", ]
  }

  # remove non-informative row names
  rownames(ranked_perm) <- NULL
  rownames(perm) <- NULL
  rownames(ranked_comb) <- NULL
  rownames(comb) <- NULL
  rownames(phier) <- NULL
  rownames(hdr) <- NULL

  # add appropriate brackets for combinatorial type
  if(nrow(ranked_perm) > 0) {
    ranked_perm$`Ranked Permutations` <- paste0("(", ranked_perm$`Ranked Permutations`, ")")
  }
  if(nrow(perm) > 0) {
    perm$Permutations <- paste0("(", perm$Permutations, ")")
  }
  if(nrow(ranked_comb) > 0) {
    ranked_comb$`Ranked Combinations` <- paste0("{", ranked_comb$`Ranked Combinations`, "}")
  }
  if(nrow(comb) > 0) {
    comb$Combinations <- paste0("{", comb$Combinations, "}")
  }

  return(list(`Ranked Permutations` = ranked_perm,
              `Permutations` = perm,
              `Ranked Combinations` = ranked_comb,
              `Combinations` = comb,
              `Partial Hierachies` = phier,
              `HDR` = hdr))

}
