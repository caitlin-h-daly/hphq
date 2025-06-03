#' Prepare data to determine credible hierarchies
#'
#' @description
#' `prep_data()` produces a list of inputs to be used altogether in
#' `get_hierarchies()`, or individually in `get_arranagements()`,
#' `get_cred_phier()`, and `singular_treatment()`.
#'
#' @param effects_matrix a data frame where the column headers are treatment
#'   names and each row displays each treatment’s sampled relative effect for
#'   that iteration.
#' @param reference a character string of the reference treatment's name.
#' @param larger_better a logical value indicating whether larger relative
#'   effects are better (TRUE) or not (FALSE).
#'
#' @return \code{hierarchy_matrix} - A matrix where column headers are ranks and
#' each row displays the treatments assigned to each rank for that iteration.
#' @return \code{effects_matrix} - A matrix where the column headers are
#' treatment names and each row displays each treatment’s sampled relative
#' effect for that iteration.
#' @return \code{ranking_df} - A data frame of each treatment's ranks and
#' associated frequencies.
#' @importFrom stats reshape
#' @importFrom stringr str_remove_all
#' @importFrom stringr str_split_1
#' @importFrom stringr str_detect
#' @importFrom stringr str_count
#' @export
#'
#' @examples
#' inputs <- prep_data(effects_matrix = dat_Thijs2008[, -1], reference = "Placebo", largerbetter = FALSE)
#' head(inputs$hierarchy_matrix)
#' head(inputs$effects_matrix)
#' head(inputs$ranking_df)
prep_data <- function(effects_matrix, reference, larger_better) {
  treatments <- colnames(effects_matrix)
  if (reference %in% treatments) {
    x <- effects_matrix
  } else {
    warning("Relative effects for reference treatment not detected. A vector of 0's for the reference treatment has been added, assuming the input is an MCMC sample from a Bayesian framework. Ensure sampled relative effects are on the additive scale such that the null effect is 0.")
    x <- data.frame(0, effects_matrix)
    colnames(x)[1] <- reference
    treatments <- colnames(x)
  }
  for (trt in treatments) {
    if (str_detect(trt, ",")) {
      stop("Comma detected in treatment name. Please remove all commas from treatment names.")
    }
  }

  n_trt <- length(treatments)
  value_matrix <- as.matrix(x)
  x2 <- x
  n_iter <- nrow(x2)
  x2$iteration <- seq_len(n_iter)
  x3 <- reshape(x2,
                varying = names(x2)[-ncol(x2)],
                v.names = "value",
                timevar = "trt",
                times = names(x2)[-ncol(x2)],
                direction = "long")
  x3 <- x3[, -which(names(x3) == "id")]

  if (larger_better == TRUE) {
    x3_sorted <- x3[order(x3$iteration, -x3$value), ]
  } else {
    x3_sorted <- x3[order(x3$iteration, x3$value), ]
  }

  # prep for arrangements
  strs <- sapply(split(x3_sorted$trt, x3_sorted$iteration), function(x) paste(x, collapse = ","))
  tbl<- t(sapply(strs, str_split_1, ","))

  inputs <- vector("list", length=3)
  inputs[[1]] <- tbl
  inputs[[2]] <- value_matrix

  # prep for single treatment
  table_list <- list()

  for (i in seq(1:n_trt)) {
    temp_table <- table(tbl[,i]) / n_iter
    temp_df <- data.frame(names(temp_table), i, as.numeric(temp_table))
    colnames(temp_df) <- c("Treatment","Rank","Freq")
    if (nrow(temp_df) != n_trt) {
      x <- temp_df$Treatment
      missing <- setdiff(treatments, x)
      vec <- data.frame(missing, i, 0)
      colnames(vec) <- c("Treatment", "Rank", "Freq")
      temp_df <- rbind(temp_df, vec)
    }
    table_list[[i]] <- temp_df
  }
  inputs[[3]] <- do.call(rbind, table_list)
  names(inputs) <- c("hierarchy_matrix", "effects_matrix", "ranking_df")
  return (inputs)
}
