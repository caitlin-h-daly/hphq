#' Get credible HDR sets for each treatment
#'
#' @description
#' `get_ranks_by_treatment()` finds all highest density region (HDR) sets with
#' empirical probabilities greater than or equal to the threshold for each
#' treatment. The HDR sets provide the subset of ranks with the smallest
#' possible cumulative empirical probability that is at least equal to the
#' threshold.
#'
#' @param ranking_df a data frame of each treatment's ranks and associated
#'   empirical probabilities.
#' @param threshold a proportion between 0 and 1 for which a hierarchy must be
#'   observed in order to be credible.
#' @param print_plot a logical value indicating whether the rankograms should be
#' printed (TRUE) or not (FALSE, the default).
#'
#' @return A data frame containing the credible HDR set for each treatment.
#'
#' @importFrom graphics barplot
#' @importFrom graphics legend
#' @importFrom graphics par
#' @export
#'
#' @examples
#' inputs <- prep_data(effects_matrix = dat_Thijs2008[, -1], reference = "Placebo", larger_better = FALSE)
#' get_ranks_by_treatment(ranking_df = inputs$ranking_df, threshold = 0.9, print_plot = FALSE)
get_ranks_by_treatment <- function(ranking_df, threshold, print_plot = FALSE) {

  if(threshold > 1 || threshold < 0) {
    stop("Please ensure threshold value is between 0 and 1")
  }

  df <- ranking_df
  treatments <- unique(df$Treatment)
  n_trt <- length(treatments)
  tolerance <- .Machine$double.eps ^ 0.5
  threshold <- threshold - tolerance
  comparator <- seq_len(n_trt)

  if(print_plot) {
    sucra <- matrix(, nrow = n_trt, ncol = n_trt - 1,
                    dimnames = list(treatments, seq_len(n_trt - 1)))
    for(i in 1:(n_trt - 1)) {
      temp_df <- subset(df, df$Rank == i)
      sucra[, i] <- (temp_df[match(treatments, temp_df$Treatment), ])$pi_hat
    }
    sucra_matrix <- t(apply(sucra, 1, cumsum))
    sucra_values <- rowMeans(sucra_matrix)
    sorted_sucra <- sort(sucra_values, decreasing = TRUE)
  } else {
    sorted_sucra <- treatments
    names(sorted_sucra) <- treatments
  }

  # pi_hat for hdr
  hdr_list <- vector("list", length = n_trt)
  rank_list <- vector("list", length = n_trt)
  cdr_list <- vector("list", length = n_trt)
  for(i in 1:n_trt) {
    trt_name <- names(sorted_sucra[i])
    ranks <- (subset(df, df$Treatment == trt_name))[-1]
    hdr_cdr_list <- hdr(ranks, threshold, 1)
    hdr_vec <- hdr_cdr_list$hdr
    rank_list[[i]] <- hdr_vec[[3]]
    hdr_vec <- data.frame(hdr_vec[[1]], hdr_vec[[2]])
    hdr_list[[i]] <- cbind(trt_name, hdr_vec)
    colnames(hdr_list[[i]]) <- c("Treatment", "HDR Rank(s)", "Sum of pi_hat")

    # other credible sets of the same size
    if(length(hdr_cdr_list$alt_cdr) > 0) {
      all_cdr <- hdr_cdr_list$alt_cdr
      cdr_list[[i]] <- vector("list", length = length(hdr_cdr_list$alt_cdr))

      for(j in 1:length(all_cdr)){
        cdr_temp <- all_cdr[[j]]
        cdr_temp <- data.frame(cdr_temp[[1]], cdr_temp[[2]])
        cdr_list[[i]][[j]] <- cbind(trt_name, cdr_temp)
        colnames(cdr_list[[i]][[j]]) <- c("Treatment", "CDR Rank(s)", "Sum of pi_hat")
      }
      cdr_list[[i]] <- do.call(rbind, cdr_list[[i]])
    }
  }

  hdr_df <- do.call(rbind, hdr_list)
  row.names(hdr_df) <- NULL

  cdr_df <- do.call(rbind, cdr_list)
  row.names(cdr_df) <- NULL

  if(is.null(cdr_df)) {
    cdr_df <- hdr_df[FALSE, ]
    colnames(cdr_df) <- c("Treatment", "CDR Rank(s)", "Sum of pi_hat")
  }

  if(print_plot) {
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
      barplot(trt_df$pi_hat, names.arg = x_ranks,
              col = colour_vec,
              main = paste("Rank of", trt_name),
              sub = paste("SUCRA = ", sucra_val),
              xlab = "Rank", ylab = "Empirical Probability",
              ylim = c(0, 1))
      legend("topright",
             legend = c("Non HDR", "HDR*"),
             fill = c("black", "lightblue"), bty = "n")
    }
    par(mfrow = c(1, 1))
  }

  return(list(HDR = hdr_df, "Alternative CDR of same size as HDR" = cdr_df))
}
