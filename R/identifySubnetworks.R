# ©2023 Rutgers, The State University of New Jersey, All rights reserved. Do not copy or reproduce without permission.  

library(viper)

#' @export
identifySubnetworks <- function(weightedNetwork,exp_mat,phenos,interactome,pathwayGenesets=NULL)
{   
  # compute signature phenos
  signaturePhenos<- list()
  # phenos <- unlist(phenos)
  for(i in 1:(length(phenos)-1)){
    p1<- exp_mat[,phenos[[i+1]]]
    p2<- exp_mat[,phenos[[i]]]
    signaturePhenos[[i]]<- ttest(p1,p2)
  }


  # compute marina
  marinaList <-list()
  for(i in 1:length(signaturePhenos)){
    currentSignature <- signaturePhenos[[i]]
    marinaForCurrentSignature <- viper::msviper(currentSignature[,1],interactome,minsize = 5)
    TRs <- summary(marinaForCurrentSignature,length(interactome))
    TRs <- TRs[complete.cases(TRs),]
    marinaList[[i]]<- TRs
  }

  ## determine the leading TRs
  leadingEdgeTRList <- list()
  for(i in 1: (length(marinaList)-1)){
    currentMarina = marinaList[[i]]
    nextMarina = marinaList[[i+1]]
    if((i+1)%%2==0){
      # even case
      sigTRS <- nextMarina[which((nextMarina$NES>0)&(nextMarina$p.value<0.05)),]
    }
    else {
      # odd case
      sigTRS <- nextMarina[which((nextMarina$NES<0)&(nextMarina$p.value<0.05)),]
    }
    query <- sigTRS$Regulon
    reflistTR <- as.matrix(currentMarina[,c("NES")])
    rownames(reflistTR) <- currentMarina$Regulon
    reflistTR <- t(reflistTR)
    r <- ssgsea(reflistTR[1,],query)
    leadingEdgeTRs <- as.matrix(r$LE.genes)
    leadingEdgeTRList[[i]]<- leadingEdgeTRs
  }


  message('computed leading edge')

  # pathway enrichment
  pathwayEnrichmentList =list()

  
  n.cores <- parallel::detectCores() - 1
  my.cluster <- parallel::makeCluster(n.cores)
  doParallel::registerDoParallel(cl = my.cluster)
  
  for (i in 1:(length(signaturePhenos))){
    enrichedPathway <- pathwayEnrichment(signaturePhenos[[i]],pathwayGenesets)
    pathwayEnrichmentList[[i]]<-  enrichedPathway
  }

  
  parallel::stopCluster(my.cluster)
  

  # determine leading edge pathways
  leadEdgePathwayList =list()
  for(i in 1:(length(pathwayEnrichmentList)-1)){
    currentPathway = pathwayEnrichmentList[[i]]
    nextPathway = pathwayEnrichmentList[[i+1]]
    if((i+1)%%2==0){
      # even case
      sigPathway <- nextPathway[which((nextPathway[,1]>0)&(nextPathway[,2]<0.05)),]
    }else{
      # odd case
      sigPathway <- nextPathway[which((nextPathway[,1]<0)&(nextPathway[,2]<0.05)),]
    }
    sigPathway <- row.names(sigPathway)
    reflistPath <- currentPathway
    reflistPath <- t(reflistPath)
    rpath <- gsea(reflistPath[1,],sigPathway)
    leadingEdgePathway <- as.matrix(rpath$LE.genes)
    leadEdgePathwayList[[i]]<- leadingEdgePathway
  }

  

  # intersection  for pathways
  leadingEdgePathwaysFinal <- leadEdgePathwayList[[1]]
  if(length(leadEdgePathwayList)>1 ){
    for(i in 2: (length(leadEdgePathwayList)-1)){
      leadingEdgePathwaysFinal <- intersect(leadingEdgePathwaysFinal,leadEdgePathwayList[[i]] )
    }
  }


  # intersection  for TRs
  leadingEdgeTRSFinal <- leadingEdgeTRList[[1]]
  if(length(leadingEdgeTRList)>1){
    for(i in 2: (length(leadingEdgeTRList)-1)){
      leadingEdgeTRSFinal <- intersect(leadingEdgeTRSFinal,leadingEdgeTRList[[i]])
    }
  }
  
  LEpathwaysinnetwork <- weightedNetwork[match(leadingEdgePathwaysFinal,weightedNetwork$Pathway),]
  LEpathwaysinnetwork = LEpathwaysinnetwork[complete.cases(LEpathwaysinnetwork),]
  subnetworklist = list()
  for(i in 1:nrow(LEpathwaysinnetwork))
  {
    ## identify each leading edge pathway
    LEPathway <- as.character(LEpathwaysinnetwork[i,c("Pathway")])
    ## identify TRs and then leading edge TRs associated to each pathway
    LEPathwayNetwork <- weightedNetwork[which(weightedNetwork$Pathway==LEPathway),]
    subNetwork <- LEPathwayNetwork[match(leadingEdgeTRSFinal,LEPathwayNetwork$TRs),]
    subNetwork <- subNetwork[complete.cases(subNetwork),]
    if(nrow(subNetwork)>1 & !is.null(subNetwork)){
    subnetworklist[[LEPathway]] = subNetwork
    }
  }  
  return(subnetworklist)
}


## Run ttest between two phenotypes
ttest <- function(pheno1Mat,pheno2Mat){
  tval = matrix(0,nrow(pheno1Mat),1)
  pval = matrix(0,nrow(pheno1Mat),1)
  pheno2MatwrtPheno1 = pheno2Mat[match(rownames(pheno1Mat),rownames(pheno2Mat)),]
  for(i in 1:nrow(pheno1Mat))
  {
    r <- t.test(pheno1Mat[i,],pheno2MatwrtPheno1[i,])
    tval[i,1] <- r$statistic
    pval[i,1] <- r$p.value
  }

  signature <- cbind(tval,pval)
  rownames(signature) <- row.names(pheno1Mat)
  colnames(signature) <- c("tval","pval")
  return(signature)
}
