#' Determines if credible hierarchies are supersets
#'
#' @description
#' Adds an additional column to each data frame, where TRUE indicates it is a superset, and FALSE indicates it is a subset
#'
#' @param consec_outputs Output from consecutive(), which consists of a list of data frames containing the credible hierarchies for Ranked Permutations, Permutations, Ranked Combinations, and Combinations
#' @param posets Output from get_posets(), which consists of a data frame of the credible Posets
#' @param n.treatments Number of treatments
#'
#' @return A list of credible hierarchies with a column indicating which hierarchies are supersets
#' @importFrom stringr str_remove_all
#' @importFrom stringr str_split_1
#' @importFrom stringr str_detect
#' @importFrom stringr str_count
#' @export
#'
#' @examples
#' get_superset(consec_outputs,posets,3)

get_superset <- function(consec_outputs,posets,n.treatments) {
  perm<-consec_outputs[[2]]
  perm[,1]<-str_remove_all(perm[,1], "[()]")
  ranked_perm<-consec_outputs[[1]]
  ranked_perm[,2]<-str_remove_all(ranked_perm[,2], "[()]")
  combo<-consec_outputs[[4]]
  combo[,1]<-str_remove_all(combo[,1], "[{}]")
  ranked_combo<-consec_outputs[[3]]
  ranked_combo[,2]<-str_remove_all(ranked_combo[,2], "[{}]")
  MID <- as.numeric((str_split_1(colnames(posets)[1], " "))[5])

  # Superset for nonconsec
  if (nrow(posets) > 0) {
    size_vec <- sapply(posets[,1], function(x) {
      str_count(x, ">") # gets size of all the permutations
    })
    y<-cbind(posets,size_vec)
    y <- y[order(y$size_vec), ] # order it by size
    nonconsec_list<-split(y, y$size_vec) # split it into a list based on size

    n.size <- length(nonconsec_list)
    if (n.size > 1) { # checks within itself
      for(i in 1:(n.size-1)) {
        current_df <- nonconsec_list[[i]][-3] # drop the Size column
        Superset <- rep(FALSE,nrow(current_df))
        for (j in 1:nrow(current_df)) {
          target <-str_split_1(as.character(current_df[j,1])," > ")
          Superset[j] <-terminal_nonconsec(nonconsec_list[(i+1):n.size],target)
        }
        nonconsec_list[[i]] <-cbind(current_df,Superset)
      }
    }
    if(n.size >= 1) {
      Superset<-FALSE
      nonconsec_list[[n.size]] <-cbind(nonconsec_list[[n.size]][-3],Superset) # biggest size defaults to FALSE
    }
    posets<-do.call(rbind,nonconsec_list)
    if(nrow(perm) > 0 && MID == 0) { # checks against permutations
      for (i in 1:nrow(posets)) {
        if(posets[i,3] == FALSE) {
          target<-str_split_1(as.character(posets[i,1])," > ")
          posets[i,3]<-compare_nonconsec(perm,target)
        }
      }
    }
  } else {
    Superset<-rep(TRUE,0)
    posets<-data.frame(posets,Superset)
  }
  posets <- posets[order(posets$Freq, decreasing = TRUE), ]

  # Superset for permutations
  n.perm <- nrow(perm)
  n.rankedperm <- nrow(ranked_perm)
  Superset<-rep(TRUE,n.perm)
  if (n.perm > 0) {
    for(i in 1:n.perm) {
      target <-as.character(perm[i,1])
      Superset[i] <- calc_terminal_unranked(perm,target,i,n.perm)
      if (!Superset[i] && n.rankedperm > 0) {
        for(j in 1:n.rankedperm) {
          new_string<-as.character(ranked_perm[j,2])
          if (str_detect(new_string,target)) {
            Superset[i]<-TRUE
            break
          }
        }
      }
    }
  }
  consec_outputs[[2]]<- cbind(consec_outputs[[2]],Superset)

  # Superset for combinations
  n.combo <- nrow(combo)
  n.rankedcombo<-nrow(ranked_combo)
  Superset<-rep(FALSE,n.combo)
  if ((n.perm>0 || n.rankedcombo>0) && (n.combo > 0)) {
    for(i in 1:n.combo) {
      target <-str_split_1(combo[i,1],",")
      if(nrow(perm)>0){
        for(j in 1:nrow(perm)) {
          if (setequal(target,str_split_1(perm[j,1],","))) {
            Superset[i] <- TRUE
            break
          }
        }
      }
      if(!Superset[i] && (nrow(ranked_combo)>0)) {
        for(j in 1:nrow(ranked_combo)) {
          if (setequal(target,str_split_1(ranked_combo[j,2],","))) {
            Superset[i] <- TRUE
            break
          }
        }
      }
    }
  }
  consec_outputs[[4]] <- cbind(consec_outputs[[4]],Superset)

  # Superset for ranked permutations
  Superset<-rep(FALSE,n.rankedperm)
  if (n.rankedperm > 0) {
    for(i in 1:n.rankedperm) {
      current_range <- str_split_1(ranked_perm[i,1],"-")
      target <-str_split_1(as.character(ranked_perm[i,2]),",")
      Superset[i] <- calc_terminal_ranked(ranked_perm,current_range,target,i,n.rankedperm)
    }
  }
  consec_outputs[[1]] <- cbind(consec_outputs[[1]],Superset)

  # Superset for ranked combinations
  Superset<-rep(FALSE,n.rankedcombo)
  ranked_combo<-cbind(ranked_combo,Superset)
  if (n.rankedcombo > 0) {
    if (n.rankedperm > 0) {
      for(i in 1:n.rankedcombo) {
        range<-as.character(ranked_combo[i,1])
        target <-str_split_1(ranked_combo[i,2],",")
        for(j in 1:n.rankedperm) {
          if (range == (as.character(ranked_perm[j,1])) && setequal(target,str_split_1(ranked_perm[j,2],","))) {
            ranked_combo[i,4] <- TRUE
            break
          }
        }
      }
    }

    intervals.drop.left <- paste0(1, "-", (floor(n.treatments/2)+1):(n.treatments-2))
    intervals.drop.right <- paste0(3:(floor(n.treatments/2)+1), "-", n.treatments)
    intervals <- c(intervals.drop.left,intervals.drop.right)
    rows_drop <- which(ranked_combo$Range %in% intervals)

    for (index in rows_drop) {
      if(ranked_combo[index,4]) {
        interval <- as.numeric(str_split_1(ranked_combo[index,1],"-"))
        trts<-ranked_combo[index,2]
        if (interval[1] == 1) {
          complement_interval<-paste0((interval[2]+1),"-",n.treatments)
        } else {
          complement_interval<-paste0(1,"-",(interval[1]-1))
        }
        cond1 <- (ranked_combo[,1] == complement_interval)
        cond2 <- sapply(ranked_combo[,2],function(x) {(anyDuplicated(c(x, trts)) == 0)})
        row_complement <- which(cond1 & cond2)
        ranked_combo[row_complement,4] <- TRUE
      }
    }
    ranked_combo <- ranked_combo[-(rows_drop), ]
    ranked_combo$`Ranked Combinations` <- paste0("{", ranked_combo$`Ranked Combinations`, "}")
  }
  consec_outputs[[3]] <- ranked_combo
  output_list <- list(consec_outputs[[1]],consec_outputs[[2]],consec_outputs[[3]],consec_outputs[[4]],posets)
  names(output_list)<-c("Ranked Permutations","Permutations","Ranked Combinations","Combinations","Strict Posets")
  return(output_list)
}
