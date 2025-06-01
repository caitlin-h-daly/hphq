#' Identify and optionally trim redundant hierarchies
#'
#' @description
#' `find_redundant_hierarchies` compares the outputs from `get_arrangements()`
#' and `get_phier()` to identify redundant credible hierarchies.
#'
#' @param algo_1 a list of data frames outputted from `get_arrangements()`
#'   containing the credible hierarchies for ranked permutations, permutations,
#'   ranked combinations, and combinations.
#' @param algo_2 a data frame outputted from `get_phier()` containing credible
#'   partial hierarchies.
#' @param type a numeric vector indicating what types of hierarchy redundancies
#'   should be trimmed. See details for more information. Default is 1:8.
#' @param trim_redundant a logical value indicating whether the redundant
#'   hierarchies should be trimmed from the output (TRUE) or not (FALSE, the
#'   default).
#'
#' @details
#' Hierarchy redundancy types can be identified with numbers through 1 through 8
#' as follows:
#' 1) Within partial hierarchies (e.g., A > B is redundant because of A > B > C).
#' 2) Between partial hierarchies and permutations (e.g., A > B is redundant
#'   because of (A, B)).
#' 3) Between combinations and ranked combinations (e.g., ${A, B}$ is redundant
#'   because of ${A, B}_1^2$).
#' 4) Between combinations and ranked combinations (e.g., ${A, B}$ is redundant
#'   because of $(A, B)$).
#' 5) Between ranked combinations and ranked permutations (e.g., ${A, B}_1^2$ is
#'   redundant because of $(A, B)_1^2$).
#' 6) Within permutations (e.g., (A, B) is redundant because of (A, B, C)).
#' 7) Between permutations and ranked permutations(e.g., $(A, B)$ is redundant
#'   because of $(A, B)_1^2$).
#' 8) Within ranked permutations (e.g., $(A, B)_1^2$ is redundant because of
#'   $(A, B, C)_1^3$).
#'
#' We recommend these hierarchy redundancies are identified in increasing order
#' (1 through 8). Otherwise, a warning will be outputted.
#'
#' @return A list of the credible ranked permutations, permutations, ranked
#' combinations, combinations, and partial hierarchies along with their
#' redundancy status.
#' @export
#'
#' @examples
#' inputs <- prep_data(effects_matrix = dat_Thijs2008[, -1], reference = "Placebo", largerbetter = FALSE)
#' algo1 <- get_arrangements(hierarchy_matrix = inputs$hierarchy_matrix, threshold = 0.9)
#' algo2 <- get_cred_phier(effects_matrix = inputs$effects_matrix, mid = 0, threshold = 0.9, larger_better = FALSE)
#' find_redundant_hierarchies(algo1, algo2)
find_redundant_hierarchies <- function(algo_1, algo_2, type = 1:8, trim_redundant = FALSE) {

  # Verify algo_1 and algo_2 objects correspond to the same threshold
  # TO DO - add threshold to list in algo_1 and algo_2 outputs

  # Verify algo_2 objects correspond to MID = 0
  if(gsub("Treatments at MID = ", "", colnames(algo_2)[1]) != 0){
    warning("Partial hierarchies constructed with a non-zero MID; partial hierarchies that are made redundant by permutations will not be identified.")
    if(2 %in% type) {
      type <- type[-which(type == 2)]  # do not search for partial hierarchies that are redundant because of permutations
    }
  }

  # Verify if inputted types are in increasing order
  if(any(type != sort(type))) {
    warning("`types` is not in increasing order; some redundant hierarchies may be missed.")
  }

  # Extract different hierarchy types
  ranked_perm <- algo_1[[1]]
  ranked_perm[, 2] <- str_remove_all(ranked_perm[, 2], "[()]")
  perm <- algo_1[[2]]
  perm[, 1] <- str_remove_all(perm[, 1], "[()]")
  ranked_comb <- algo_1[[3]]
  ranked_comb[, 2] <- str_remove_all(ranked_comb[, 2], "[{}]")
  comb <- algo_1[[4]]
  comb[, 1] <- str_remove_all(comb[, 1], "[{}]")
  phier <- algo_2

  # Now find redundant hierarchies
  for(ind in type){

    if(ind == 1 & nrow(phier) > 0) {

      if(!exists("Redundant", phier)) {
        phier$Redundant <- FALSE
      }

      # size of all credible partial hierarchies
      phier$phier_size <- sapply(phier[, 1], function(x) { stringr::str_count(x, ">") + 1 } )
      phier <- phier[order(phier$phier_size), ]

      # create a list of credible hierarchies by size
      phier_list <- split(phier, phier$phier_size)

      # unique sizes of credible partial hierarchies
      phier_sizes <- sort(unique(phier$phier_size))

      # check if credible partial hierarchies are redundant within
      if(length(phier_sizes) > 1) {
        # assess whether smaller credible partial hierarchies are within larger credible partial hierarchies
        for(i in phier_sizes[-which.max(phier_sizes)]) {
          # extract credible partial hierarchies of size i
          current_phier_size <- phier_list[[as.character(i)]]
          # list of larger hierarchies to check against
          phier_larger_list <- phier_list[as.character(phier_sizes[which(phier_sizes > i)])]
          for(j in 1:nrow(current_phier_size)) {
            if(current_phier_size[j, "Redundant"] == FALSE) {
              # extract jth partial hierarchy in current_phier_size
              phier_target <- stringr::str_split_1(as.character(current_phier_size[j, 1]), " > ")
              # check to see if jth partial hierarchy is in larger partial hierarchies
              redundant_check <- nmahierarchies:::is_phier_redundant_within_phier(phier_target,
                                                                                  phier_larger_list)
              if(!is.null(redundant_check)){
                current_phier_size[j, "Redundant"] <- redundant_check
              }
            }
          }
          phier_list[[as.character(i)]] <- current_phier_size
        }
      }
      phier <- do.call(rbind, phier_list)[-4]  # drop the phier_size column

    } else if(ind == 2 & nrow(phier) > 0 & nrow(perm) > 0) {

      if(!exists("Redundant", phier)) {
        phier$Redundant <- FALSE
      }

      # check if (remaining) credible partial hierarchies are redundant because of credible permutations
      for (i in 1:nrow(phier)) {
        if(phier[i, "Redundant"] == FALSE) {
          phier_target <- str_split_1(as.character(phier[i, 1]), " > ")
          for(j in 1:nrow(perm)) {
            perm_string_to_check <- str_split_1(as.character(perm[j, 1]), ",")
            # finds position of treatments in phier_target within perm_string_to_check
            positions <- match(phier_target, perm_string_to_check)
            if (!anyNA(positions) && all(diff(positions) >= 0)) { # checks if position[j+1] >= position[j]
              phier[i, "Redundant"] <- TRUE
              break
            }
          }
        }
      }


    } else if(ind == 3 & nrow(comb) > 0 & nrow(ranked_comb) > 0) {

      if(!exists("Redundant", comb)) {
        comb$Redundant <- FALSE
      }

      # check if (remaining) credible combinations are redundant because of credible ranked combinations
      for(i in 1:nrow(comb)) {
        comb_target <- str_split_1(comb[i, 1], ",")
        if(comb[i, "Redundant"] == FALSE) {
          for(j in 1:nrow(ranked_comb)) {
            if(setequal(comb_target, str_split_1(ranked_comb[j, 2], ","))) {
              comb[i, "Redundant"] <- TRUE
              break
            }
          }
        }
      }

    } else if(ind == 4 & nrow(comb) > 0 & nrow(perm) > 0) {

      if(!exists("Redundant", comb)) {
        comb$Redundant <- FALSE
      }

      # check if (remaining) credible combinations are redundant because of credible permutations
      for(i in 1:nrow(comb)) {
        comb_target <- str_split_1(comb[i, 1], ",")
        if(comb[i, "Redundant"] == FALSE) {
          for(j in 1:nrow(perm)) {
            if(setequal(comb_target, str_split_1(perm[j, 1], ","))) {
              comb[i, "Redundant"] <- TRUE
              break
            }
          }
        }
      }

    } else if(ind == 5 & nrow(ranked_comb) > 0 & nrow(ranked_perm) > 0) {

      if(!exists("Redundant", ranked_comb)) {
        ranked_comb$Redundant <- FALSE
      }

      # check if (remaining) credible ranked combinations are redundant because of credible ranked permutations
      for(i in 1:nrow(ranked_comb)) {
        range_target <- as.character(ranked_comb[i, 1])
        ranked_comb_target <- str_split_1(ranked_comb[i, 2], ",")
        if(ranked_comb[i, "Redundant"] == FALSE) {
          for(j in 1:nrow(ranked_perm)) {
            if(range_target == as.character(ranked_perm[j, 1]) &&
               setequal(ranked_comb_target, str_split_1(ranked_perm[j, 2], ","))) {
              ranked_comb[i, "Redundant"] <- TRUE
              break
            }
          }
        }
      }

    } else if(ind == 6 & nrow(perm) > 0) {

      if(!exists("Redundant", perm)) {
        perm$Redundant <- FALSE
      }

      # check if (remaining) credible permutations are redundant within
      perm_ind <- 1:nrow(perm)
      for(i in perm_ind){
        if(perm[i, "Redundant"] == FALSE) {
          perm_target <- as.character(perm[i, 1])
          for(j in perm_ind[-i]) {
            perm_string <- as.character(perm[j, 1])
            if(str_detect(perm_string, perm_target)) {
              perm[i, "Redundant"] <- TRUE
              break
            }
          }
        }
      }

    } else if(ind == 7 & nrow(perm) > 0 & nrow(ranked_perm) > 0) {

      if(!exists("Redundant", perm)) {
        perm$Redundant <- FALSE
      }

      # check if (remaining) credible permutations are redundant because of credible ranked permutations
      for(i in 1:nrow(perm)) {
        perm_target <- str_split_1(perm[i, 1], ",")
        if(perm[i, "Redundant"] == FALSE) {
          for(j in 1:nrow(ranked_perm)) {
            if(setequal(perm_target, str_split_1(ranked_perm[j, 2], ","))) {
              perm[i, "Redundant"] <- TRUE
              break
            }
          }
        }
      }

    } else if(ind == 8 & nrow(ranked_perm) > 0) {

      if(!exists("Redundant", ranked_perm)) {
        ranked_perm$Redundant <- FALSE
      }

      # check if credible ranked permutations are redundant within
      ranked_perm_ind <- 1:nrow(ranked_perm)
      for(i in ranked_perm_ind){
        if(ranked_perm[i, "Redundant"] == FALSE) {
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
                ranked_perm[i, "Redundant"] <- TRUE
                break
              }

            }

          }
        }
      }

    } else {
      stop(paste0("Invalid type = ", ind, "; should be a numeric value between 1 and 8."))
    }

  }

  if(trim_redundant == TRUE){
    ranked_perm <- ranked_perm[ranked_perm$Redundant == FALSE, ]
    perm <- perm[perm$Redundant == FALSE, ]
    ranked_comb <- ranked_comb[ranked_comb$Redundant == FALSE, ]
    comb <- comb[comb$Redundant == FALSE, ]
    phier <- phier[phier$Redundant == FALSE, ]
  }

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
              `Partial Hierachies` = phier))

}
