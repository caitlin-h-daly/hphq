#' Data preparation
#'
#' @description
#' Produces a list of inputs to be used altogether in get_hierarchies(), or individually in consecutive(), get_posets(), and singular_treatment()
#'
#' @param df A data frame where the column headers are treatment names and row headers are the iteration number. Each row displays each treatment’s sampled relative effect for that iteration
#' @param reference A string of the reference treatment's name
#' @param largerBetter Boolean value where True indicates larger values are better, and False otherwise
#'
#' @return \code{hierarchy.matrix} - A matrix where column headers are ranks and each row displays the treatments assigned to each rank for that iteration
#' @return \code{effects.matrix} - A matrix where the column headers are treatment names and each row displays each treatment’s sampled relative effect for that iteration
#' @return \code{ranking.df} - A data frame of each treatment's ranks and associated frequencies
#' @importFrom stats reshape
#' @export
#'
#' @examples
#' prep_data(df,"waitlist",TRUE)
prep_data <- function(df,reference,largerBetter) {
  treatments<-colnames(df)
  if (reference %in% treatments) {
    x <- df
  } else {
    warning("Relative effects for reference treatment not detected. A vector of 0's for the reference treatment has been added, assuming the input is an MCMC sample from a Bayesian framework. Ensure sampled relative effects are on the additive scale such that the null effect is 0.")
    x<-data.frame(0,df)
    colnames(x)[1] <- reference
    treatments <-colnames(x)
  }
  for (trt in treatments) {
    if (str_detect(trt,",")) {
      stop("Comma detected in treatment name. Please remove all commas from treatment names.")
    }
  }

  n.treatments<-length(treatments)
  value_matrix <- as.matrix(x)
  x2 <- x
  n.iterations <- nrow(x2)
  x2$iteration <- seq_len(n.iterations)
  x3 <- reshape(x2,
                varying = names(x2)[-ncol(x2)],
                v.names = "value",
                timevar = "trt",
                times = names(x2)[-ncol(x2)],
                direction = "long")
  x3 <- x3[, -which(names(x3) == "id")]

  if (largerBetter == TRUE) {
    x3_sorted <- x3[order(x3$iteration, -x3$value), ]
  } else {
    x3_sorted <- x3[order(x3$iteration, x3$value), ]
  }

  # prep for consecutive
  strs <- sapply(split(x3_sorted$trt, x3_sorted$iteration), function(x) paste(x, collapse = ","))
  tbl<- t(sapply(strs, str_split_1, ","))

  inputs <- vector("list",length=3)
  inputs[[1]] <- tbl
  inputs[[2]] <- value_matrix

  # prep for single treatment
  table_list<-list()

  for (i in seq(1:n.treatments)) {
    temp_table <- table(tbl[,i])/n.iterations
    temp_df <- data.frame(names(temp_table),i,as.numeric(temp_table))
    colnames(temp_df) <- c("Treatment","Rank","Freq")
    if (nrow(temp_df) != n.treatments) {
      x <- temp_df$Treatment
      missing <- setdiff(treatments, x)
      vec <- data.frame(missing,i,0)
      colnames(vec) <- c("Treatment","Rank","Freq")
      temp_df <- rbind(temp_df,vec)
    }
    table_list[[i]] <- temp_df
  }
  inputs[[3]]<-do.call(rbind,table_list)
  names(inputs)<-c("hierarchy.matrix","effects.matrix","ranking.df")
  return (inputs)
}
