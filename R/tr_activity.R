# ©2023 Rutgers, The State University of New Jersey, All rights reserved. Do not copy or reproduce without permission.  

library(viper)

#' Estimates single-sample activity level of each transcriptional regulatory program
#'
#' This function computes Normalized Enrichment Score (NES) of each transcriptional regulatory (TR) program for each sample
#'
#' @param exprMat Normalized gene expression matrix where rows are genes and columns are samples
#' @param interactome Object of class regulon
#'
#' @return A matrix of activity level for each TR across samples
#' @export
#'
#' @examples
#' computeTRActivity(exprMat,interactome)
computeTRActivity <- function(exprMat,interactome){
  message("computing TR activity...")
  #transpose expression matrix
  scaled_normalized_exp_t <- t(exprMat)
  # compute NES
  tr_activity <- viper::viper(scaled_normalized_exp_t,interactome,method= c("none"),minsize=5)
  message("finished computing TR activity")
  return(tr_activity)
}
