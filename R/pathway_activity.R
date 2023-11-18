# ©2023 Rutgers, The State University of New Jersey, All rights reserved. Do not copy or reproduce without permission.  

library(doParallel)
library(foreach)

`%dopar%` <- foreach::`%dopar%`
`%do%` <- foreach::`%do%`

ss_pathway_activity<-function(scaledExprMat,pathwayGenelist,nperms=10){
  
  #create and register cluster
  n.cores <- parallel::detectCores() - 1

  pathway_activity_cluster <- parallel::makeCluster(n.cores)
  doParallel::registerDoParallel(cl = pathway_activity_cluster)

  pb <- txtProgressBar(min = 0,      # Minimum value of the progress bar
                       max = length(pathwayGenelist), # Maximum value of the progress bar
                       style = 3,    # Progress bar style (also available style = 1 and style = 2)
                       width = 50,   # Progress bar width. Defaults to getOption("width")
                    char = "=")
 result = list()
 
 for (i in 1:length(pathwayGenelist)) {
    r1<-foreach(j=1:nrow(scaledExprMat),.combine = 'c') %dopar% {
      r<-ssgsea(scaledExprMat[j,], pathwayGenelist[[i]],nperms)$NES
    }
    result[[i]]<-r1
    # update the progress bar
    setTxtProgressBar(pb, i)

  }

  # stop the cluster
  parallel::stopCluster(pathway_activity_cluster)

  pathway_activity_data<- matrix(unlist(result), ncol =length(pathwayGenelist), nrow = nrow(scaledExprMat))
  return(pathway_activity_data)

}


## Estimates single-sample activity level of molecular pathway from geneset provided by user
#'
#' This function computes single-sample activity level of each pathway from C2 collections
#'
#' @param scaledExprMat Scaled gene expression matrix where rows are samples and columns are genes
#'
#' @return Matrix of activity level of each pathway from userdefined genelist across samples
#'
#' @examples
#' scaledExprMat <- scale(t(exprMat))
#' ssPathActivityC2(scaledExprMat)
ssPathActivity <- function(scaledExprMat,pathwayGenelist,nperms) {
  message('computing NES for given pathways....')
  normalized_samples <- scaledExprMat
  pathway_names<-names(pathwayGenelist)
  # compute unique pathway names
  gn <- unique(unlist(pathwayGenelist))
  ## convert gene names to geneids
  colnames(normalized_samples)<-sapply(colnames(normalized_samples),
                                          function(x) m2h$hs[which(m2h$hs.gn==x)])
  normalized_samples <- normalized_samples[,which(colnames(normalized_samples)%in%gn)]
  
  
  genesets <- list()
  for (i in 1:length(pathway_names)) {
    pi<-which(names(pathwayGenelist)==pathway_names[i])
    genesets[[i]] <- pathwayGenelist[[pi]]
  }
  
  
  # estimate activity level of each user pathway in each sample
  pathway_activity <- ss_pathway_activity(normalized_samples,genesets,nperms)
  # match rownames & column names
  
  rownames(pathway_activity) <-   rownames(normalized_samples)
  colnames(pathway_activity) <- pathway_names
  
  message('\n finished computing NES for User defined  pathways')
  return(pathway_activity)
}




#' Estimates single-sample activity level of molecular pathways from C2 collections
#'
#' This function computes single-sample activity level of each pathway from C2 collections
#'
#' @param scaledExprMat Scaled gene expression matrix where rows are samples and columns are genes
#'
#' @return Matrix of activity level of each pathway from C2 across samples
#'
#' @examples
#' scaledExprMat <- scale(t(exprMat))
#' ssPathActivityC2(scaledExprMat)
ssPathActivityC2 <- function(scaledExprMat,nperms=10) {
  #C2 start ------------------
  message('computing NES for C2 pathways....')
  normalized_samples_C2 <- scaledExprMat
  # extracting KEGG, BioCarta and Reactome pathways
  pathways_C2 <- c(grep("KEGG",names(C2)), grep("BIOCARTA", names(C2)),
                grep("REACTOME", names(C2)))
  #extracting names from the pathways
  pathway_names_C2<-names(C2)[pathways_C2]
  # compute unique pathway names
  gn <- unique(unlist(C2))
  ## convert gene names to geneids
  colnames(normalized_samples_C2)<-sapply(colnames(normalized_samples_C2),
                                            function(x) m2h$hs[which(m2h$hs.gn==x)])
  normalized_samples_C2 <- normalized_samples_C2[,which(colnames(normalized_samples_C2)%in%gn)]
  

  genesets <- list()
  for (i in 1:length(pathway_names_C2)) {
    pi<-which(names(C2)==pathway_names_C2[i])
    genesets[[i]] <- C2[[pi]]
  }
  

  # estimate activity level of each C2 pathway in each sample
  pathway_activity_C2 <- ss_pathway_activity(normalized_samples_C2,genesets,nperms)
  # match rownames & column names
  
  rownames(pathway_activity_C2) <-   rownames(normalized_samples_C2)
  colnames(pathway_activity_C2) <- pathway_names_C2

  message('\n finished computing NES for C2 pathways')
  return(pathway_activity_C2)
  #C2 end ------------------
}

#' Estimates single-sample activity level of molecular pathways from Hallmark collections
#'
#' This function estimates single-sample activity level of each pathway from Hallmark collections
#'
#' @param scaledExprMat Scaled gene expression matrix where rows are samples and columns are genes
#'
#' @return Matrix of activity level of each pathway from Hallmark across samples
#'
#' @examples
#' scaledExprMat <- scale(t(exprMat))
#' ssPathActivityHlmrk(scaledExprMat)
ssPathActivityHlmrk <- function(scaledExprMat,nperms=10) {
  #hallmark start ------------------
  normalized_samples_hallmark <- scaledExprMat
  message('computing NES for Hallmark pathways...')
  # compute unique pathway names
  hallmark_gn <- unique(unlist(hallmarksdb))
  # identify hallmark genes
  normalized_samples_hallmark <- normalized_samples_hallmark[,which(colnames(normalized_samples_hallmark)%in%hallmark_gn)]
  # extract hallmark pathways and pathway names
  pathways_hallmark <- c(grep("HALLMARK",names(hallmarksdb)))
  pathway_names_hallmark <-names(hallmarksdb)[pathways_hallmark]
  
  genesets <- list()
  for (i in 1:length(pathway_names_hallmark)) {
    pi<-which(names(hallmarksdb)==pathway_names_hallmark[i])
    genesets[[i]] <- hallmarksdb[[pi]]
  }
  
  
  # estimate activity level of each Hallmark pathway in each sample
  pathway_activity_hallmark <- ss_pathway_activity(normalized_samples_hallmark,genesets,nperms)
  rownames(pathway_activity_hallmark) <- rownames(normalized_samples_hallmark)
  colnames(pathway_activity_hallmark) <- pathway_names_hallmark
  message('\n finished computing NES for Hallmark pathways')
  #hallmark end ------------------
  return(pathway_activity_hallmark)
}


#' Estimates single-sample activity level of molecular pathways from Hallmark collections
#'
#' This function estimates single-sample activity level of each pathway from Hallmark collections
#'
#' @param scaledExprMat Scaled gene expression matrix where rows are samples and columns are genes
#'
#' @return Matrix of activity level of each pathway from Hallmark across samples
#' @export
#' @examples
#' computePathwayActivity(exprMat,pathwaygenelist,nperms=10)
computePathwayActivity<- function(scaledExprMat,pathwayGenelist=NULL,nperms=10){

NES_of_pathway <-NULL
if(is.null(pathwayGenelist)) {
  pathway_activity_C2 <- ssPathActivityC2(scaledExprMat,nperms)
  pathway_activity_hallmark<-ssPathActivityHlmrk(scaledExprMat,nperms)
  
  #### cbind pathway_activity_C2 and pathway_activity_hallmark
  # combine pathway activity
  NES_hallmark_wrt_C2 <- pathway_activity_hallmark[match(rownames(pathway_activity_C2),rownames(pathway_activity_hallmark)),]
  NES_of_pathway = cbind(pathway_activity_C2,NES_hallmark_wrt_C2)
}else{
  
  # discarding mismatched genes
  exprGenes <- rownames(scaledExprMat)
  
  if(is.list(pathwayGenelist)==FALSE){
    stop("The pathway gene list should be a list")
  }
  is_numeric_string <- function(x) {
    !is.na(as.numeric(x))
  }
  
  
  if(all(sapply(unlist(pathwayGenelist), is_numeric_string))==FALSE){
    stop("the pathway related genes in pathwaygenelist should be gene ids")
  }
  # discard unwanted genes 
  
  NES_of_pathway = ssPathActivity(scaledExprMat,pathwayGenelist,nperms)
}

return(NES_of_pathway)
  

}
