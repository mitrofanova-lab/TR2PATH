# ©2023 Rutgers, The State University of New Jersey, All rights reserved. Do not copy or reproduce without permission.  

library(reshape2)
library(doParallel)
library(parallel)



generateNetwork<- function(NESPathway,TRActivity,sampling=FALSE) {
  # checking samples in tr activity vector is similar to samplpes in pathway
  tr_activity_wrt_pathway = TRActivity[,match(rownames(NESPathway),colnames(TRActivity))]
  # initializing slope and pvalues
  slope = matrix(0,nrow(tr_activity_wrt_pathway),ncol(NESPathway))
  pvalue = matrix(0,nrow(tr_activity_wrt_pathway),ncol(NESPathway))
  ind = sample(1:ncol(tr_activity_wrt_pathway),ncol(tr_activity_wrt_pathway),replace= TRUE)
  for(i in 1:ncol(NESPathway)){
    pathway_NES = as.matrix(NESPathway[,i])
    
    if(sampling){
      pathway_NES = as.matrix(NESPathway[ind,i])
    }
    for(j in 1:nrow(tr_activity_wrt_pathway)){
      TR_activity = as.matrix(tr_activity_wrt_pathway[j,])
      TR_activity_wrt_pathwy = as.matrix(TR_activity[match(rownames(pathway_NES),rownames(TR_activity)),])
     
      X<-cbind(1,as.matrix(TR_activity_wrt_pathwy))
      # X<-pathway_NES
      y<-pathway_NES
      
      m<-NULL
      
      m <- .lm.fit(X,y)
      
      rss <- sum(m$residuals^2)
      rdf <- length(y) - ncol(X)
      resvar <- rss/rdf
      R <- chol2inv(m$qr)
      se <- sqrt(diag(R) * resvar)
      pval<-2*pt(abs(m$coef/se),rdf,lower.tail=FALSE)
    
      slope[j,i] = coef(m)[2]
      pvalue[j,i] = pval[2]
     
     rm(list = c("m", "rss", "rdf","resvar","R","se","pval"))
      
    }
  }
  fdr_p = matrix(0,nrow(pvalue),ncol(pvalue))
  for(i in 1:ncol(pvalue))
  {
    pvalue_mat = as.matrix(pvalue[,i])
    fdr_p[,i] = p.adjust(pvalue_mat,method="fdr")
  }
  
  # assign TR names to rows and pathway names to columns of slope/pvalue matrix
  rownames(slope) = rownames(tr_activity_wrt_pathway)
  colnames(slope) = colnames(NESPathway)
  
  rownames(pvalue) = rownames(tr_activity_wrt_pathway)
  colnames(pvalue) = colnames(NESPathway)
  
  rownames(fdr_p) = rownames(tr_activity_wrt_pathway)
  colnames(fdr_p) = colnames(NESPathway)
  
  rVal <-  list(
    slope=slope,
    pvalue=pvalue,
    fdr_p = fdr_p)
 rm(list = c("tr_activity_wrt_pathway", "NESPathway", "slope","pvalue","fdr_p","TRActivity"))
  return ( rVal)
}

  #' Reconstruction of mechanism-centric regulatory network
  #'
  #' The function reconstructs the TR-2-PATH network utilizing activity level of transcriptional regulatory programs (TRs) and activity level of molecular pathways
  #'
  #' @param NESPathway Matrix of activity level of each pathway across samples
  #' @param TRActivity Matrix of activity level of each TR across samples
  #' @param fdr.cutoff Numeric [default 0.05]
  #'
  #' @return A data frame containing following components:
  #' \itemize{
  #'  \item \strong{TRs} Transcriptional regulatory programs
  #'  \item \strong{Pathway} Molecular pathways
  #'  \item\strong{Slope} Beta coefficient from linear regression analysis
  #'  \item \strong{pvalue} Pvalue from linear regression analysis
  #'   \item\strong{fdr} false discovery rate
  #'   }
  #' @export
  #'
  #' @examples
  #'  \dontrun{
  #'  TestPathData <- TR2PATH::TestNESPathwayData
  #'  TestTRData<-TR2PATH::TestTRActivityData
  #' reconstNetwork(TestPathData,TestTRData,fdr.cutoff=0.05)
  #' }
  #'
  reconstNetwork <- function(NESPathway,TRActivity,fdr.cutoff=0.05,bootstrap=FALSE,bootstrap_count=10) {
    message("performing tr2path network reconstruction...")

    orgnl_ntrwk = generateNetwork(NESPathway,TRActivity,sampling=FALSE) # orginal ntrwk

    if(bootstrap){
      random_networks = generateRandomNetowrks(NESPathway,TRActivity,sampleCount=bootstrap_count) # orginal ntrwk

      edg_wts = computeEdgewts(random_networks,orgnl_ntrwk,fdr.cutoff)

      # combining slope|pvalue|fdr
      melted_Slope = melt(orgnl_ntrwk$slope)
      melted_Slope$combined = paste(melted_Slope$Var1,melted_Slope$Var2,sep="_")
      melted_pvalue = melt(orgnl_ntrwk$pvalue)
      melted_pvalue$combined = paste(melted_pvalue$Var1,melted_pvalue$Var2,sep="_")
      melted_fdr = melt(orgnl_ntrwk$fdr_p)
      melted_fdr$combined = paste(melted_fdr$Var1,melted_fdr$Var2,sep="_")

      melted_edges = melt(edg_wts)
      melted_edges$combined = paste(melted_edges$Var1,melted_edges$Var2,sep="_")

      melted_pvalue_wrt_slope = melted_pvalue[match(melted_Slope$combined,melted_pvalue$combined),]
      melted_fdr_wrt_slope = melted_fdr[match(melted_Slope$combined,melted_fdr$combined),]

      melted_edges_wrt_slope = melted_edges[match(melted_Slope$combined,melted_edges$combined),]


      network_reconst = cbind(melted_Slope[,1:3],melted_pvalue_wrt_slope$value,melted_fdr_wrt_slope$value,melted_edges_wrt_slope$value)
      colnames(network_reconst) = c("TRs","Pathway","Slope","pvalue","fdr","edge_weights")
      # computing edges of the network with fdr cutoff
      weighted_tr2_path = network_reconst[which(network_reconst$fdr<fdr.cutoff),]
      weighted_tr2_path = weighted_tr2_path[which(weighted_tr2_path$edge_weights>79),]

      message("finished tr2path network reconstruction")
      return (weighted_tr2_path)
    }else{
      # combining slope|pvalue|fdr
      melted_Slope = melt(orgnl_ntrwk$slope)
      melted_Slope$combined = paste(melted_Slope$Var1,melted_Slope$Var2,sep="_")
      melted_pvalue = melt(orgnl_ntrwk$pvalue)
      melted_pvalue$combined = paste(melted_pvalue$Var1,melted_pvalue$Var2,sep="_")
      melted_fdr = melt(orgnl_ntrwk$fdr_p)
      melted_fdr$combined = paste(melted_fdr$Var1,melted_fdr$Var2,sep="_")

      melted_pvalue_wrt_slope = melted_pvalue[match(melted_Slope$combined,melted_pvalue$combined),]
      melted_fdr_wrt_slope = melted_fdr[match(melted_Slope$combined,melted_fdr$combined),]

      network_reconst = cbind(melted_Slope[,1:3],melted_pvalue_wrt_slope$value,melted_fdr_wrt_slope$value)
      colnames(network_reconst) = c("TRs","Pathway","Slope","pvalue","fdr")
      # computing edges of the network with fdr cutoff
      unwghtd_tr2_path = network_reconst[which(network_reconst$fdr<fdr.cutoff),]
      message("finished tr2path network without bootstrap reconstruction")
      return (unwghtd_tr2_path)
    }
  }

## Generate Random Networks
generateRandomNetowrks<- function(NESPathway,TR_activity,sampleCount=100) {
    n.cores <- parallel::detectCores() - 1
    randomNetworkCluster <- parallel::makeCluster(n.cores, type="PSOCK")
    doParallel::registerDoParallel(cl = randomNetworkCluster)
    randomNetworks <- foreach(i=1:sampleCount) %dopar% generateNetwork(NESPathway,TR_activity,sampling=TRUE)
    # stop the cluster
    parallel::stopCluster(randomNetworkCluster)
    
    
    return(randomNetworks)
}

## compute Edge Weights
computeEdgewts <- function(random_networks,orgnl_ntwrk,fdr.cutoff=0.05){
  bootstrap_count <- length(random_networks)
  test_matrix <- random_networks[[1]]$slope
  coefficients_of_MRS_for_each_pathway = matrix(0,nrow(test_matrix),bootstrap_count)
  fdr_of_MRS_for_each_pathway = matrix(0,nrow(test_matrix),bootstrap_count)
  counts = matrix(0,nrow(test_matrix),ncol(test_matrix))

  for (i in 1:ncol(test_matrix))
  {
    # slope of 100 networks ( 1st column)
    for(j in 1:bootstrap_count)
    {
      coefficients_of_MRS_for_each_pathway[,j] = random_networks[[j]]$slope[,i]
      fdr_of_MRS_for_each_pathway[,j] = random_networks[[j]]$fdr_p[,i]
    }

    for(k in 1:nrow(test_matrix))
    {
      coeff_for_MR = as.matrix(coefficients_of_MRS_for_each_pathway[k,])
      pvalue_for_MR = as.matrix(fdr_of_MRS_for_each_pathway[k,])
      data = cbind(coeff_for_MR,pvalue_for_MR)
      original_coeff_slope =orgnl_ntwrk$slope[k,i]
      original_coeff_pval = orgnl_ntwrk$fdr_p[k,i]
      if((original_coeff_slope >0) & (original_coeff_pval< fdr.cutoff))
      {
        val = data[which((data[,1]>0)&(data[,2]< fdr.cutoff)),]
        val = as.matrix(val)
        counts[k,i] = nrow(val)

      } else if ((original_coeff_slope <0)& (original_coeff_pval< fdr.cutoff)){
        val = data[which((data[,1]<0)&(data[,2]< fdr.cutoff)),]
        val = as.matrix(val)
        counts[k,i] = nrow(val)
      } else{
        counts[k,i] = 0
      }
    }
  }

  rownames(counts) = row.names(test_matrix)
  colnames(counts) = colnames(test_matrix)

  return(counts)

}
