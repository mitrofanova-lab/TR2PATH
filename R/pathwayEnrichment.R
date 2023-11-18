# ©2023 Rutgers, The State University of New Jersey, All rights reserved. Do not copy or reproduce without permission.  

pathwayEnrichment <- function(signaturebtwPhenos,pathwayGenesets){
  
  pathwayEnrichment <- NULL
  
  if(is.null(pathwayGenesets)){
    ## perform pathway enrichment for C2
    pathwayEnrichmentforC2 = pathwayEnrichmentC2(signaturebtwPhenos)
    message('finshied c2 enrich')
    # perform pathway enrichment for hallmark
    pathwayEnrichmentforHallmark = pathwayEnrichmentHallmark(signaturebtwPhenos)
    message('finshed halmark enrich')
    pathwayEnrichment = rbind(pathwayEnrichmentforC2,pathwayEnrichmentforHallmark)
  }else{
   
    if(is.list(pathwayGenesets)==FALSE){
      stop("The pathway gene list should be a list")
    }
    is_numeric_string <- function(x) {
      !is.na(as.numeric(x))
    }
    
    
    if(all(sapply(unlist(pathwayGenesets), is_numeric_string))==FALSE){
      stop("the pathway related genes in pathwaygenelist should be gene ids")
    }
    pathwayEnrichment = pathwayEnrichmentuserIP(signaturebtwPhenos,pathwayGenesets)
    
  }
  
  return(pathwayEnrichment)
}


pathwayEnrichmentuserIP <- function(signaturebtwPhenos,pathwayGenesets){
  reflistbtwPheno = t(signaturebtwPhenos)
  colnames(reflistbtwPheno)<-sapply(colnames(reflistbtwPheno),function(x) m2h$hs[which(m2h$hs.gn==x)])
  gn <- unique(unlist(pathwayGenesets))
  reflistbtwPheno <- reflistbtwPheno[,which(colnames(reflistbtwPheno)%in%gn)]
  pathwayNames<-names(pathwayGenesets)
  
  result <- foreach(
    i = 1:length(pathwayNames), 
    .combine = 'c'
  ) %dopar% {
    pi<-which(names(pathwayGenesets)==pathwayNames[i])
    
    r<-ssgsea(reflistbtwPheno[1,], pathwayGenesets[[pi]])
    NES = r$NES
    pvalue = r$p.value
    gseaVal = cbind(NES,pvalue)
    gseaVal
  }
  
  
  pathwayEnrichUser<- matrix(unlist(result), length(pathwayNames),2)
  
  rownames(pathwayEnrichUser) = pathwayNames
  colnames(pathwayEnrichUser) = c("NES","pvalue")
  return(pathwayEnrichUser)
}



## pathway enrichment for HALLMARK pathways using signature
pathwayEnrichmentHallmark <- function(signaturebtwPhenos){

  reflistbtwPheno = t(signaturebtwPhenos)
  hallmarkgn <- unique(unlist(hallmarksdb))

  reflistbtwPheno <- reflistbtwPheno[,which(colnames(reflistbtwPheno)%in%hallmarkgn)]
  pathways <- c(grep("HALLMARK",names(hallmarksdb)))    ## 50
  pathwayNames<-names(hallmarksdb)[pathways]

  result <- foreach(
    i = 1:length(pathwayNames), 
    .combine = 'c'
  ) %dopar% {
      pi<-which(names(hallmarksdb)==pathwayNames[i])

      r<-ssgsea(reflistbtwPheno[1,], hallmarksdb[[pi]])
      NES = r$NES
      pvalue = r$p.value
      gsea_val = cbind(NES,pvalue)
      gsea_val
  }
  
  pathwayEnrichHallmark<- matrix(unlist(result), length(pathwayNames),2)
  
  rownames(pathwayEnrichHallmark) = pathwayNames
  colnames(pathwayEnrichHallmark) = c("NES","pvalue")
  return(pathwayEnrichHallmark)
}

## pathway enrichment for C2 pathways using signature
pathwayEnrichmentC2 <- function(signaturebtwPhenos){
  reflistbtwPheno = t(signaturebtwPhenos)
  colnames(reflistbtwPheno)<-sapply(colnames(reflistbtwPheno),function(x) m2h$hs[which(m2h$hs.gn==x)])
  gn <- unique(unlist(C2))
  reflistbtwPheno <- reflistbtwPheno[,which(colnames(reflistbtwPheno)%in%gn)]
  pathways <- c(grep("KEGG",names(C2)), grep("BIOCARTA", names(C2)),
                grep("REACTOME", names(C2)))    ## 833
  pathwayNames<-names(C2)[pathways]
  
  result <- foreach(
    i = 1:length(pathwayNames), 
    .combine = 'c'
  ) %dopar% {
    pi<-which(names(C2)==pathwayNames[i])
    
    r<-ssgsea(reflistbtwPheno[1,], C2[[pi]])
    NES = r$NES
    pvalue = r$p.value
    gseaVal = cbind(NES,pvalue)
    gseaVal
  }
  
  
  pathwayEnrichC2<- matrix(unlist(result), length(pathwayNames),2)
  
  rownames(pathwayEnrichC2) = pathwayNames
  colnames(pathwayEnrichC2) = c("NES","pvalue")
  return(pathwayEnrichC2)
}

