#' Get credible rankings and HPD sets for each treatment
#'
#' @description
#' `get_cred_hier_single()` finds all ranks and high probability density (HPD)
#' sets with relative frequencies greater than or equal to the threshold for
#' each treatment. The HPD intervals provide the subset of ranks with the
#' smallest possible cumulative relative frequency that is at least equal to the
#' threshold.
#'
#' @param ranking_df a data frame of each treatment's ranks and associated
#'   frequencies.
#' @param threshold a proportion between 0 and 1 for which a hierarchy must be
#'   observed in order to be credible.
#' @param printPlot a logical value indicating whether the rankograms should be
#' printed (TRUE) or not (FALSE, the default).
#'
#' @return A list of data frames containing the credible rankings and HPD
#' set for each treatment.
#'
#' @importFrom graphics barplot
#' @importFrom graphics legend
#' @importFrom graphics par
#' @export
#'
#' @examples
#' get_cred_hier_single(df,0.1,c("CBT.exp","Placebo","PCT"))
get_cred_hier_single <- function(ranking_df, threshold, printPlot = FALSE) {

  if(threshold > 1 || threshold < 0) {
    stop("Please ensure threshold value is between 0 and 1")
  }

  df <- ranking_df
  treatments <- unique(df$Treatment)
  n_trt <- length(treatments)
  tolerance <- .Machine$double.eps ^ 0.5
  threshold <- threshold - tolerance
  comparator <- seq_len(n_trt)
  outputs <- vector("list", length = 2)

  # proportion of times each treatment is a specific rank
  filtered_df <- subset(df, df$Freq > threshold)
  filtered_df <- filtered_df[order(filtered_df$Freq, decreasing = TRUE), ]
  outputs[[1]] <- filtered_df

  if(printPlot) {
    sucra <- matrix(, nrow = n_trt, ncol = n_trt - 1,
                    dimnames = list(treatments, seq_len(n_trt - 1)))
    for(i in 1:(n_trt - 1)) {
      temp_df <- subset(df, df$Rank == i)
      sucra[, i] <- (temp_df[match(treatments, temp_df$Treatment), ])$Freq
    }
    sucra_matrix <- t(apply(sucra, 1, cumsum))
    sucra_values <- rowMeans(sucra_matrix)
    sorted_sucra <- sort(sucra_values, decreasing = TRUE)
  } else {
    sorted_sucra <- treatments
    names(sorted_sucra) <- treatments
  }

  # freq for hpd
  hpd_list <- vector("list", length = n_trt)
  rank_list <- vector("list", length = n_trt)
  for(i in 1:n_trt) {
    trt_name <- names(sorted_sucra[i])
    ranks <- (subset(df, df$Treatment == trt_name))[-1]
    hpd_vec <- hpd(ranks, threshold, 1)
    rank_list[[i]] <- hpd_vec[[3]]
    hpd_vec <- data.frame(hpd_vec[[1]], hpd_vec[[2]])
    hpd_list[[i]] <- cbind(trt_name, hpd_vec)
    colnames(hpd_list[[i]]) <- c("Treatment", "HPD Rank(s)", "Sum of Freq")
  }

  hpd_df <- do.call(rbind, hpd_list)
  outputs[[2]]<-hpd_df

  if(printPlot) {
    rows <- ceiling(sqrt(n_trt))
    cols <- ceiling(n_trt / rows)
    par(mfrow = c(rows, cols))

    # rankograms
    for(i in 1:n_trt) {
      rank_vec <- rank_list[[i]]
      colour_vec <- rep("black", n_trt)
      x_ranks <- as.list(comparator)
      for(x in rank_vec) {
        colour_vec[x] <-"lightblue"
        x_ranks[x] <- paste0(x, "*")
      }
      trt_name <- names(sorted_sucra[i])
      trt_df <- (subset(df, df$Treatment == trt_name))
      sucra_val<-round(sorted_sucra[i],3)
      barplot(trt_df$Freq, names.arg = x_ranks,
              col = colour_vec,
              main = paste("Rank of", trt_name),
              sub = paste("SUCRA = ", sucra_val),
              xlab = "Rank", ylab = "Frequency",
              ylim = c(0, 1))
      legend("topright",
             legend = c("Non HPD", "HPD*"),
             fill = c("black", "lightblue"), bty = "n")
    }
    par(mfrow = c(1, 1))
  }

  names(outputs) <- c("Individual Rank", "HPD")
  return(outputs)

}
