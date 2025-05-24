#' Identify and optionally trim superset hierarchies
#'
#' @description
#' `find_supersets` compares the outputs from `get_arrangements()` and
#' `get_phier()` to identify redundant credible hierarchies (i.e., supersets).
#'
#' @param algo1 a list of data frames outputted from `get_arrangements()`
#'   containing the credible hierarchies for ranked permutations, permutations,
#'   ranked combinations, and combinations.
#' @param algo2 a data frame outputted from `get_phier()` containing credible
#'   partial hierarchies.
#' @param type a numeric vector indicating what types of supersets should be
#'   trimmed. See details for more information. Default is 1:8.
#' @param trim a logical value indicating whether the output should trim the
#'   include the identified subsets (TRUE) or not (FALSE, the default).
#'
#' @details
#' Superset types can be identified with numbers through 1 through 8 as follows:
#' 1) Within partial hierarchies (e.g., A > B is a superset of A > B > C).
#' 2) Between partial hierarchies and permutations (e.g., A > B is a superset of
#'   (A, B)).
#' 3) Between combinations and ranked combinations (e.g., ${A, B}$ is a superset
#'   of ${A, B}_1^2$).
#' 4) Between combinations and ranked combinations (e.g., ${A, B}$ is a superset
#'   of $(A, B)$).
#' 5) Between ranked combinations and ranked permutations (e.g., ${A, B}_1^2$ is
#'   a superset of $(A, B)_1^2$).
#' 6) Within permutations (e.g., (A, B) is a superset of (A, B, C)).
#' 7) Between permutations and ranked permutations(e.g., $(A, B)$ is a superset
#'   of $(A, B)_1^2$).
#' 8) Within ranked permutations (e.g., $(A, B)_1^2$ is a superset of
#'   $(A, B, C)_1^3$).
#'
#' We recommend these supersets are identified in increasing order (1 through
#' 8). Otherwise, a warning will be outputted.
#'
#' @return A list of the credible ranked permutations, permutations, ranked
#' combinations, combinations, and partial hierarchies along with their superset
#' status.
#' @export
#'
#' @examples
#' find_supersets(algo1, algo2)
find_supersets <- function(algo1, algo2, type = 1:8, trim = FALSE) {

  # Verify algo1 and algo2 objects correspond to the same threshold
  # TO DO - add threshold to list in algo1 and algo2 outputs

  # Verify algo2 objects correspond to MID = 0
  if(gsub("Treatments at MID = ", "", colnames(algo2)[1]) != 0){
    warning("Partial hierarchies constructed with a non-zero MID; supersets of permutations will not be identified.")
    if(2 %in% type) {
      type <- type[-which(type == 2)]  # do not search for partial hierarchies that are supersets of permutations
    }
  }

  # Verify if inputted types are in increasing order
  if(any(type != sort(type))) {
    warning("`types` is not in increasing order; some supersets may be missed.")
  }

  # Extract different hierarchy types
  ranked_perm <- algo1[[1]]
  ranked_perm[, 2] <- str_remove_all(ranked_perm[, 2], "[()]")
  perm <- algo1[[2]]
  perm[, 1] <- str_remove_all(perm[, 1], "[()]")
  ranked_comb <- algo1[[3]]
  ranked_comb[, 2] <- str_remove_all(ranked_comb[, 2], "[{}]")
  comb <- algo1[[4]]
  comb[, 1] <- str_remove_all(comb[, 1], "[{}]")
  phier <- algo2

  # Now find supersets
  for(ind in type){

    if(ind == 1 & nrow(phier) > 0) {

      if(!exists("Superset", phier)) {
        phier$Superset <- FALSE
      }

      # size of all credible partial hierarchies
      phier$phier_size <- sapply(phier[, 1], function(x) { stringr::str_count(x, ">") + 1 } )
      phier <- phier[order(phier$phier_size), ]

      # create a list of credible hierarchies by size
      phier_list <- split(phier, phier$phier_size)

      # unique sizes of credible partial hierarchies
      phier_sizes <- sort(unique(phier$phier_size))

      # check if credible partial hierarchies are supersets within
      if(length(phier_sizes) > 1) {
        # assess whether smaller credible partial hierarchies are within larger credible partial hierarchies
        for(i in phier_sizes[-which.max(phier_sizes)]) {
          # extract credible partial hierarchies of size i
          current_phier_size <- phier_list[[as.character(i)]]
          # list of larger hierarchies to check against
          phier_larger_list <- phier_list[as.character(phier_sizes[which(phier_sizes > i)])]
          for(j in 1:nrow(current_phier_size)) {
            if(current_phier_size[j, "Superset"] == FALSE) {
              # extract jth partial hierarchy in current_phier_size
              phier_target <- stringr::str_split_1(as.character(current_phier_size[j, 1]), " > ")
              # check to see if jth partial hierarchy is in larger partial hierarchies
              superset_check <- nmahierarchies::is_phier_sup_of_phier(phier_target,
                                                                      phier_larger_list)
              if(!is.null(superset_check)){
                current_phier_size[j, "Superset"] <- superset_check
              }
            }
          }
          phier_list[[as.character(i)]] <- current_phier_size
        }
      }
      phier <- do.call(rbind, phier_list)[-4]  # drop the phier_size column

    } else if(ind == 2 & nrow(phier) > 0 & nrow(perm) > 0) {

      if(!exists("Superset", phier)) {
        phier$Superset <- FALSE
      }

      # check if (remaining) credible partial hierarchies are supersets of credible permutations
      for (i in 1:nrow(phier)) {
        if(phier[i, "Superset"] == FALSE) {
          phier_target <- str_split_1(as.character(phier[i, 1]), " > ")
          for(j in 1:nrow(perm)) {
            perm_string_to_check <- str_split_1(as.character(perm[j, 1]), ",")
            # finds position of treatments in phier_target within perm_string_to_check
            positions <- match(phier_target, perm_string_to_check)
            if (!anyNA(positions) && all(diff(positions) >= 0)) { # checks if position[j+1] >= position[j]
              phier[i, "Superset"] <- TRUE
              break
            }
          }
        }
      }


    } else if(ind == 3 & nrow(comb) > 0 & nrow(ranked_comb) > 0) {

      if(!exists("Superset", comb)) {
        comb$Superset <- FALSE
      }

      # check if (remaining) credible combinations are supersets of credible ranked combinations
      for(i in 1:nrow(comb)) {
        comb_target <- str_split_1(comb[i, 1], ",")
        if(comb[i, "Superset"] == FALSE) {
          for(j in 1:nrow(ranked_comb)) {
            if(setequal(comb_target, str_split_1(ranked_comb[j, 2], ","))) {
              comb[i, "Superset"] <- TRUE
              break
            }
          }
        }
      }

    } else if(ind == 4 & nrow(comb) > 0 & nrow(perm) > 0) {

      if(!exists("Superset", comb)) {
        comb$Superset <- FALSE
      }

      # check if (remaining) credible combinations are supersets of credible permutations
      for(i in 1:nrow(comb)) {
        comb_target <- str_split_1(comb[i, 1], ",")
        if(comb[i, "Superset"] == FALSE) {
          for(j in 1:nrow(perm)) {
            if(setequal(comb_target, str_split_1(perm[j, 1], ","))) {
              comb[i, "Superset"] <- TRUE
              break
            }
          }
        }
      }

    } else if(ind == 5 & nrow(ranked_comb) > 0 & nrow(ranked_perm) > 0) {

      if(!exists("Superset", ranked_comb)) {
        ranked_comb$Superset <- FALSE
      }

      # check if (remaining) credible ranked combinations are supersets of credible ranked permutations
      for(i in 1:nrow(ranked_comb)) {
        range_target <- as.character(ranked_comb[i, 1])
        ranked_comb_target <- str_split_1(ranked_comb[i, 2], ",")
        if(ranked_comb[i, "Superset"] == FALSE) {
          for(j in 1:nrow(ranked_perm)) {
            if(range_target == as.character(ranked_perm[j, 1]) &&
               setequal(ranked_comb_target, str_split_1(ranked_perm[j, 2], ","))) {
              ranked_comb[i, "Superset"] <- TRUE
              break
            }
          }
        }
      }

    } else if(ind == 6 & nrow(perm) > 0) {

      if(!exists("Superset", perm)) {
        perm$Superset <- FALSE
      }

      # check if (remaining) credible permutations are supersets within
      perm_ind <- 1:nrow(perm)
      for(i in perm_ind){
        if(perm[i, "Superset"] == FALSE) {
          perm_target <- as.character(perm[i, 1])
          for(j in perm_ind[-i]) {
            perm_string <- as.character(perm[j, 1])
            if(str_detect(perm_string, perm_target)) {
              perm[i, "Superset"] <- TRUE
              break
            }
          }
        }
      }

    } else if(ind == 7 & nrow(perm) > 0 & nrow(ranked_perm) > 0) {

      if(!exists("Superset", perm)) {
        perm$Superset <- FALSE
      }

      # check if (remaining) credible permutations are supersets of credible ranked permutations
      for(i in 1:nrow(perm)) {
        perm_target <- str_split_1(perm[i, 1], ",")
        if(perm[i, "Superset"] == FALSE) {
          for(j in 1:nrow(ranked_perm)) {
            if(setequal(perm_target, str_split_1(ranked_perm[j, 2], ","))) {
              perm[i, "Superset"] <- TRUE
              break
            }
          }
        }
      }

    } else if(ind == 8 & nrow(ranked_perm) > 0) {

      if(!exists("Superset", ranked_perm)) {
        ranked_perm$Superset <- FALSE
      }

      # check if credible ranked permutations are supersets within
      ranked_perm_ind <- 1:nrow(ranked_perm)
      for(i in ranked_perm_ind){
        if(ranked_perm[i, "Superset"] == FALSE) {
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
                ranked_perm[i, "Superset"] <- TRUE
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

  if(trim == TRUE){
    ranked_perm <- ranked_perm[ranked_perm$Superset == FALSE, ]
    perm <- perm[perm$Superset == FALSE, ]
    ranked_comb <- ranked_comb[ranked_comb$Superset == FALSE, ]
    comb <- comb[comb$Superset == FALSE, ]
    phier <- phier[phier$Superset == FALSE, ]
  }

  return(list(`Ranked Permutations` = ranked_perm,
              `Permutations` = perm,
              `Ranked Combinations` = ranked_comb,
              `Combinations` = comb,
              `Partial Hierachies` = phier))

}
