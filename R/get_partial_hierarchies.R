#' Get credible partial hierarchies
#'
#' @description
#' `get_partial_hierarchies()` finds all partial hierarchies reflecting treatment
#' differences greater than or equal to `mid` and with relative frequencies
#' greater than or equal to the threshold.
#'
#' @param effects_matrix a matrix where the column headers are treatment names
#'   and each row displays each treatment’s sampled relative effect for that
#'   iteration.
#' @param mid a numeric value indicating the absolute minimally important
#'   difference. Default is 0.
#' @param threshold a proportion between 0 and 1 for which a hierarchy must be
#'   observed in order to be credible.
#' @param larger_better a logical value indicating whether larger relative
#'   effects are better (TRUE) or not (FALSE).
#' @param freq_digits a numeric value indicating the desired number of decimal
#'   places for which the relative frequencies will be rounded. Default is 4.
#'
#' @return A data frame containing the credible partial hierarchies.
#' @importFrom utils combn
#' @export
#'
#' @examples
#' inputs <- prep_data(effects_matrix = dat_Thijs2008[, -1], reference = "Placebo", largerbetter = FALSE)
#' get_partial_hierarchies(effects_matrix = inputs$effects_matrix, mid = 0, threshold = 0.9, larger_better = FALSE)
get_partial_hierarchies <- function(effects_matrix, mid = 0, threshold, larger_better, freq_digits = 4) {

  if(threshold > 1 || threshold < 0) {
    stop("Please ensure threshold value is between 0 and 1")
  }

  treatments <- colnames(effects_matrix)
  output_list <- vector('list')
  output_list_index <- 1
  n_trt <- length(treatments)
  n_iter <- nrow(effects_matrix)
  tolerance <- .Machine$double.eps ^ 0.5
  mid <- abs(mid)
  if (larger_better) {
    mid_t <- mid - tolerance
  } else {
    mid_t <- -mid + tolerance
    mid <- -mid
  }
  threshold <- threshold - tolerance

  # permutations of size = 2
  combinations_matrix <- combn(treatments, 2)
  combos <- as.data.frame(t(combinations_matrix))
  colnames(combos) <- c("V1", "V2")
  other_combos <- data.frame(combos$V2, combos$V1)
  colnames(other_combos) <- c("V1", "V2")
  perms <- rbind(combos, other_combos)

  # finds frequencies for permutations of size 2
  Freq <- apply(perms, 1, function(row) {
    trt1 <- row['V1']
    trt2 <- row['V2']
    if(larger_better) {
      sum((effects_matrix[,trt1] - effects_matrix[,trt2]) > mid_t) / n_iter
    } else {
      sum((effects_matrix[,trt1] - effects_matrix[,trt2]) < mid_t)/ n_iter
    }
  })
  perms$Freq <- Freq

  # formatting
  filtered_perms <- subset(perms, perms$Freq > threshold)
  total_rows <- nrow(filtered_perms)
  if (total_rows == 0) {
    finished_perms <- data.frame(filtered_perms$V1, filtered_perms$Freq)
  } else {
    credible <- paste0(filtered_perms$V1, ",", filtered_perms$V2)
    formatted_perms <- paste0(filtered_perms$V1, " > ", filtered_perms$V2)
    finished_perms <- data.frame(formatted_perms, filtered_perms$Freq)
  }
  heading <- paste0("Treatments at MID = ", mid)
  colnames(finished_perms) <- c(heading, "Freq")
  output_list[[output_list_index]] <- finished_perms
  output_list_index <- output_list_index + 1
  perm_size <- 3

  # size > 2, loop repeats while there are still permutations left after
  # filtering and permutation size <= number of treatments
  while(total_rows > 0 && perm_size <= n_trt) {
    filtered_perms <- filtered_perms[, 1:ncol(filtered_perms) - 1]

    # creates permutations to look for
    df_list <- vector("list", length = total_rows)
    for (i in seq(1:total_rows)) {
      trts <- filtered_perms[i, ]
      trt1 <- trts[perm_size - 1]
      new_trts <- treatments[!(treatments %in% trts)]
      perm_df <- create_perm(trts, trt1, new_trts, perm_size, credible)
      if (!is.logical(perm_df)) {
        df_list[[i]] <- perm_df
      }
    }
    perms <- do.call(rbind, df_list)
    if(length(perms) == 0) {
      break
    }

    # calculates frequency of each permutation
    Freq <- apply(perms, 1, function(row) {
      trts <- row[]
      compare_matrix <- matrix(, nrow = n_iter, ncol = (length(trts) - 1))
      for (i in seq (1:(length(trts) - 1))) {
        trt1 <- row[i]
        trt2 <- row[i + 1]
        if (larger_better) {
          compare_matrix[,i] <- (effects_matrix[, trt1] - effects_matrix[, trt2]) > mid_t
        } else {
          compare_matrix[,i] <- (effects_matrix[, trt1] - effects_matrix[, trt2]) < mid_t
        }
      }
      count <- sum(apply(compare_matrix, 1, function(row) {
        if(row[1] && length(unique(row)) == 1) {
          return(1)
        } else {
          return(0)
        }
      }))
      count / n_iter
    })

    # formatting data
    filtered_perms <- data.frame(perms, Freq)
    filtered_perms <- subset(filtered_perms, filtered_perms$Freq > threshold)
    just_perms <- filtered_perms[, 1:ncol(filtered_perms) - 1]
    hierarchies <- apply(just_perms, 1, function(x) {
      paste(x, collapse = " > ")
    })
    all_perms <- data.frame(hierarchies, filtered_perms$Freq)
    colnames(all_perms) <- c(heading, "Freq")
    output_list[[output_list_index]] <- all_perms
    output_list_index <- output_list_index + 1
    total_rows <- nrow(filtered_perms)
    perm_size <- perm_size + 1
  }

  all_output <- do.call(rbind, output_list)
  all_output <- all_output[order(all_output$Freq, decreasing = TRUE), ]
  rownames(all_output) <- NULL
  all_output$Freq <- round(all_output$Freq, digits = freq_digits)
  return(all_output)
}
