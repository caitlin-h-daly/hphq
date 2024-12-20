#' Retrieves credible hierarchies
#'
#'@description
#' Runs the inputs through consecutive(), get_posets(), and singular_treatment(), and returns a list of hierarchies with frequencies greater than or equal to the thresholds specified. These hierarchies are then run through get_superset() to determine which of them are supersets
#'
#' @param inputs The output produced by running prep_data(), which consists of a list of the hierarchy.matrix, effects.matrix, and ranking.df
#' @param largerBetter Boolean value where True indicates larger values are better, and False otherwise
#' @param thresholds A numeric vector containing three values that define the thresholds for identifying credible hierarchies across consecutive, poset, and individual treatment question types
#' @param MID A number indicating the absolute minimally important difference
#' @param printPlot Boolean value where True indicates that the rankograms will be printed (default = False)
#'
#' @return A list of data frames containing the credible hierarchies for Ranked Permutations, Permutations, Ranked Combinations, Combinations, Posets, Individual Treatments, and HPD Intervals
#' @export
#'
#' @examples
#'get_hierarchies(inputs,"Waitlist",TRUE,c(0.5,0.6,0.7),10)

get_hierarchies <- function(inputs,largerBetter,thresholds,MID,printPlot=FALSE) {
  if(max(thresholds)>1 || min(thresholds)<0) {
    stop("Please ensure threshold values are between 0 and 1")
  }
  treatments <- colnames(inputs[[2]])
  n.treatments<-length(treatments)
  consec_outputs <-consecutive(inputs[[1]],thresholds[[1]])
  posets <- get_posets(inputs[[2]],MID,thresholds[[2]],largerBetter)
  single <- singular_treatment(inputs[[3]],thresholds[[3]],printPlot)
  first_outputs<-get_superset(consec_outputs,posets,n.treatments)
  all_outputs<-append(first_outputs,single)
  names(all_outputs)<-c("Ranked Permutations","Permutations","Ranked Combinations","Combinations","Strict Posets","Individual Rank","HPD")
  return(all_outputs)
}
