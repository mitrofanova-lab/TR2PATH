# ©2023 Rutgers, The State University of New Jersey, All rights reserved. Do not copy or reproduce without permission.  

#' Estimates single-sample activity level of each transcriptional regulatory program
#'
#' This function computes Normalized Enrichment Score (NES) of each transcriptional regulatory (TR) program for each sample
#'
#' @param exprMat Normalized gene expression matrix where rows are genes and columns are samples
#' @param interactome Object of class regulon
#'
#' @return A matrix of activity level for each TR across samples
#'
#' @examples
#'  \dontrun{
#' computeTRActivity(exprMat,interactome)
#' }
#' 
computePLS <- function(subNetwrks,TRActivity,pathwayActivity,weightedNetwork){
  anglesandCorrelationperPathwayList <- list()
  ## identify each subnetwork
    subNetwork = subNetwrks
    ## determine TR activity for each LE TR
    LEPathway = as.character(unique(subNetwork$Pathway))
    subNetworkTRactivity <- as.matrix(TRActivity[match(subNetwork$TRs,rownames(TRActivity)),])
    ## determine pathway activity for each LE pathway
    subNetworkPathwayActivity <- as.matrix(pathwayActivity[,LEPathway])
    subNetworkPathwayActivity <- t(subNetworkPathwayActivity)
    ## concatenate pathway activity and TR activity
    pathwayActivitywrtTR <- as.matrix(subNetworkPathwayActivity[,match(colnames(subNetworkTRactivity),colnames(subNetworkPathwayActivity))])
    colnames(pathwayActivitywrtTR) <- LEPathway
    subNetworkTRactivity <- t(subNetworkTRactivity)
    activityMat <- cbind(pathwayActivitywrtTR,subNetworkTRactivity)
    ## PLS anlysis
    pls1 <- plsreg1(activityMat[,2:ncol(activityMat)],activityMat[,1,drop=FALSE], comps = 2)
    # calculate Correlation with LV1 and LV2 vectors
    correlatedValues <- pls1$cor.xyt
    ## calculate angle for each TR
    angleAssociatedwithTRs <- computeAngle(correlatedValues)
    ## concatenate correlation with LV1 and LV2 and angle to matrix
    correlationwrtAngle <- correlatedValues[match(rownames(angleAssociatedwithTRs),rownames(correlatedValues)),]
    angleandCorrelation <- as.data.frame(cbind(rownames(correlationwrtAngle),correlationwrtAngle,angleAssociatedwithTRs) )
    angleandCorrelationandPathways <- as.data.frame(cbind(LEPathway,angleandCorrelation))
    rownames(angleandCorrelationandPathways) <- 1:nrow(angleandCorrelationandPathways)
    colnames(angleandCorrelationandPathways) <- c("Pathway","TRs","Correlation with LV1","Correlation with LV2","angles_of_TRS","angle_Wrt_pathway")
    anglesandCorrelationperPathwayList[[LEPathway]] <- angleandCorrelationandPathways
    
    return(anglesandCorrelationperPathwayList)
}


## compute angle
computeAngle <- function(circleofCorrelationvValues){
  angleslist = list()
  corrVal1and2Comp = circleofCorrelationvValues
  a = nrow(corrVal1and2Comp)
  correlationDist = sqrt((corrVal1and2Comp[,1])^2 + corrVal1and2Comp[,2]^2)
  
  correlationDist = as.matrix(correlationDist)
  cosThetanw = abs(corrVal1and2Comp[,1])/correlationDist
  theta = acos(cosThetanw)
  
  thetaupdated = matrix(0,nrow(theta),1)
  for(j in 1:nrow(theta))
  {
    if(corrVal1and2Comp[j,1]>0 & corrVal1and2Comp[j,2]>0)
    {
      thetaupdated[j,] = theta[j,]
    }
    else if (corrVal1and2Comp[j,1]< 0 & corrVal1and2Comp[j,2]>0)
    {
      thetaupdated[j,] = 3.14- theta[j,]
    }
    else if (corrVal1and2Comp[j,1]< 0 & corrVal1and2Comp[j,2]<0)
    {
      thetaupdated[j,] = theta[j,]+3.14
    }
    else if (corrVal1and2Comp[j,1]> 0 & corrVal1and2Comp[j,2]<0)
    {
      thetaupdated[j,] = 6.28 - theta[j,]
    }
  }
  
  rownames(thetaupdated) = rownames(theta)
  angle= rCAT::rad2deg(thetaupdated) ## need to return this 
  angle_for_TRS = as.matrix(angle[1:a-1,])
  anglewithPathway = angle-angle[a,]
  anglewithPathway = anglewithPathway[-a,]
  anglewithPathway= as.matrix(anglewithPathway)
  angles_matrix = cbind(angle_for_TRS,anglewithPathway)
  colnames(angles_matrix) = c("angles_of_TRS","angles_wrt_pathway")
  angleslist = angles_matrix
  return(angleslist)
}
