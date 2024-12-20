#' Retrieves credible hierarchies for individual treatment question types
#'
#' @description
#' Finds all treatments' ranks with frequencies greater than or equal to the threshold and HPD intervals which includes the subset of ranks with the smallest possible cumulative frequency that is at least equal to the threshold
#'
#' @param ranking.df A data frame of each treatment's ranks and associated frequencies
#' @param threshold A number that determines which hierarchies are credible
#' @param printPlot Boolean value where True indicates that the rankograms will be printed (default = False)
#'
#' @return A list of data frames containing the credible hierarchies for Individual Treatments and HPD intervals
#' @importFrom graphics barplot
#' @importFrom graphics legend
#' @importFrom graphics par
#' @export
#'
#' @examples
#' singular_treatment(df,0.1,c("CBT.exp","Placebo","PCT"))

singular_treatment <- function(ranking.df,threshold,printPlot=FALSE) {
  if(threshold>1 || threshold<0) {
    stop("Please ensure threshold value is between 0 and 1")
  }
  df<-ranking.df
  treatments<-unique(df$Treatment)
  n.treatments <- length(treatments)
  tolerance <- .Machine$double.eps^0.5
  threshold <- threshold - tolerance
  comparator <- seq_len(n.treatments)
  outputs<-vector("list",length=2)

  # proportion of times each treatment is a specific rank
  filtered_df <- subset(df, df$Freq > threshold)
  filtered_df <- filtered_df[order(filtered_df$Freq, decreasing = TRUE), ]
  outputs[[1]]<-filtered_df

  if (printPlot) {
    sucra <- matrix(,nrow=n.treatments,ncol=(n.treatments-1),dimnames = list(treatments, seq_len(n.treatments-1)))
    for (i in 1:(n.treatments-1)) {
      temp_df<-subset(df, df$Rank==i)
      sucra[,i] <- (temp_df[match(treatments, temp_df$Treatment),])$Freq
    }
    sucra_matrix <- t(apply(sucra, 1, cumsum))
    sucra_values<-rowMeans(sucra_matrix)
    sorted_sucra <- sort(sucra_values,decreasing=TRUE)
  } else {
    sorted_sucra<-treatments
    names(sorted_sucra)<-treatments
  }

  # freq for hpd
  hpd_list <- vector("list",length = n.treatments)
  rank_list <- vector("list",length = n.treatments)
  for (i in seq(1:n.treatments)) {
    trt_name <- names(sorted_sucra[i])
    ranks <- (subset(df, df$Treatment == trt_name))[-1]
    hpd_vec <- hpd(ranks,threshold,1)
    rank_list[[i]] <- hpd_vec[[3]]
    hpd_vec <- data.frame(hpd_vec[[1]],hpd_vec[[2]])
    hpd_list[[i]] <- cbind(trt_name,hpd_vec)
    colnames(hpd_list[[i]]) <- c("Treatment", "HPD Rank(s)", "Sum of Freq")
  }

  hpd_df <- do.call(rbind, hpd_list)
  outputs[[2]]<-hpd_df

  if (printPlot) {
    rows <- ceiling(sqrt(n.treatments))
    cols <- ceiling(n.treatments / rows)
    par(mfrow = c(rows, cols))

    # rankograms
    for (i in seq(1:n.treatments)) {
      rank_vec<-rank_list[[i]]
      colour_vec<-rep("black",n.treatments)
      x_ranks <- as.list(comparator)
      for (x in rank_vec) {
        colour_vec[x] <-"lightblue"
        x_ranks[x] <- paste0(x, "*")
      }
      trt_name <- names(sorted_sucra[i])
      trt_df <- (subset(df, df$Treatment == trt_name))
      sucra_val<-round(sorted_sucra[i],3)
      barplot(trt_df$Freq, names.arg = x_ranks, col = colour_vec,
              main = paste("Rank of", trt_name), sub=paste("SUCRA = ",sucra_val),xlab = "Rank", ylab = "Frequency",
              ylim = c(0, 1))
      legend("topright",
             legend = c("Non HPD", "HPD*"),
             fill = c("black","lightblue"),bty = "n")
    }
    par(mfrow = c(1, 1))
  }
  names(outputs)<-c("Individual Rank","HPD")
  return(outputs)
}
