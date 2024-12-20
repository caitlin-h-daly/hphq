#' Retrieves credible hierarchies for consecutive question types
#'
#' @description
#' Finds all Ranked Permutations, Permutations, Ranked Combinations, and Combinations with frequencies greater than or equal to the threshold
#'
#' @param hierarchy.matrix A matrix where column headers are ranks and each row displays the treatments assigned to each rank for that iteration
#' @param threshold A number that determines which hierarchies are credible
#'
#' @return A list of data frames containing the credible hierarchies for Ranked Permutations, Permutations, Ranked Combinations, and Combinations
#' @importFrom stats aggregate
#' @export
#'
#' @examples
#' consecutive(hierarchy.matrix,0.5)

consecutive <- function (hierarchy.matrix, threshold) {
  if(threshold>1 || threshold<0) {
    stop("Please ensure threshold value is between 0 and 1")
  }
  treatments<-hierarchy.matrix[1,]
  n.treatments <-  ncol(hierarchy.matrix)
  n.rankings <- nrow(hierarchy.matrix)
  tolerance <- .Machine$double.eps^0.5
  threshold <- threshold - tolerance
  consec_output <- vector("list",length=4)

  x <- ((n.treatments-1)*n.treatments)/2
  all_ranked_perm_list <- vector("list", length = x)
  index <- 1

  for(length in 2:(n.treatments)){
    for(start in 1:(n.treatments-length+1)){
      end = start + length - 1
      all_ranked_perm_list[[index]] <- (freq(hierarchy.matrix, start:end))
      index <- index + 1
    }
    print(paste0("Permutations of size ",length," completed"))
  }

  # ranked permutations
  list_length <- length(all_ranked_perm_list)
  all_ranked_perm_combo <-do.call(rbind,all_ranked_perm_list[1:(list_length-3)])
  all_ranked_perm <- rbind(all_ranked_perm_combo, do.call(rbind, all_ranked_perm_list[(list_length-2):list_length]))

  ranked_perm <- subset(all_ranked_perm, all_ranked_perm$Freq >threshold)
  names(ranked_perm)[names(ranked_perm) == 'Var1'] <- 'Ranked Permutations'
  ranked_perm <- ranked_perm[order(ranked_perm$Freq, decreasing = TRUE), ]
  consec_output[[1]]<-ranked_perm
  print("Ranked permutations completed")

  # permutations
  all_perm <- get_perm(all_ranked_perm)
  perm <- subset(all_perm, all_perm$Freq >threshold)
  perm <- perm[order(perm$Freq, decreasing = TRUE), ]
  names(perm)[names(perm) == 'Var1'] <- 'Permutations'
  consec_output[[2]] <-perm
  print("Permutations completed")

  # ranked combinations
  all_ranked_combo <- indexed_combo(all_ranked_perm_combo,treatments)
  ranked_combo <- subset(all_ranked_combo, all_ranked_combo$Freq >threshold)
  ranked_combo <- ranked_combo[order(ranked_combo$Freq, decreasing = TRUE), ]
  consec_output[[3]] <-ranked_combo
  names(consec_output[[3]])[names(consec_output[[3]]) == 'Combinations'] <- 'Ranked Combinations'
  print("Ranked combinations completed")

  # combinations
  all_combo <- get_combo(all_ranked_combo)
  combo <- subset(all_combo, all_combo$Freq >threshold)
  combo <- combo[order(combo$Freq, decreasing = TRUE), ]
  consec_output[[4]] <-combo
  print("Combinations completed")

  if(nrow(consec_output[[1]])> 0) {
    consec_output[[1]]$`Ranked Permutations` <- paste0("(", consec_output[[1]]$`Ranked Permutations`, ")")
  }
  if(nrow(consec_output[[2]])> 0) {
    consec_output[[2]]$Permutations <- paste0("(", consec_output[[2]]$Permutations, ")")
  }
  if(nrow(consec_output[[3]])> 0) {
    consec_output[[3]]$`Ranked Combinations` <- paste0("{", consec_output[[3]]$`Ranked Combinations`, "}")
  }
  if(nrow(consec_output[[4]])> 0) {
    consec_output[[4]]$Combinations <- paste0("{", consec_output[[4]]$Combinations, "}")
  }
  names(consec_output)<-c("Ranked Permutations","Permutations","Ranked Combinations","Combinations")
  return(consec_output)
}
