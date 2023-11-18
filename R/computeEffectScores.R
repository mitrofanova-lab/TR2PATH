# ©2023 Rutgers, The State University of New Jersey, All rights reserved. Do not copy or reproduce without permission.  

library(psych)
library(dplyr)

# function to make angle per quadrants
makeAngleClusters = function(pathwayTRsandAnglemat){
  final_matrix_for_angles_in_quad = list()
  matrix_for_angles_in_quad = list()
  anglesforTRs = as.matrix(as.numeric(pathwayTRsandAnglemat$angles_of_TRS))
  rownames(anglesforTRs) = pathwayTRsandAnglemat$TRs
  anglesforTRs= as.data.frame(anglesforTRs)
  
  
  # Apply categorization to value column
  anglesforTRs$quadrant <- sapply(anglesforTRs$V1, categorize_quadrant)
  
  unique_quads <- unique(anglesforTRs$quadrant)
  
  
  for (quad in unique_quads) {
    # Subset the matrix for the current value
    subset_matrix <- anglesforTRs[anglesforTRs$quadrant == quad, ]
    if(nrow(subset_matrix)<2){
      matrix_of_angles = subset_matrix
    } else{
      matrix_of_angles = matrix(0,nrow(subset_matrix),nrow(subset_matrix))
      rownames(matrix_of_angles) = row.names(subset_matrix)
      colnames(matrix_of_angles) = row.names(subset_matrix)
      
      ## calculate the angles for each quad
      for(i in 1:nrow(matrix_of_angles))
      {
        for(j in 1:nrow(matrix_of_angles))
        {
          matrix_of_angles[i,j] = abs(subset_matrix[i,1] - subset_matrix[j,1])
        }
      }
      
    }
    matrix_for_angles_in_quad[[quad]] <- matrix_of_angles
  }
  
  # change for single list 
  final_matrix_for_angles_in_quad = list(angleMat= matrix_for_angles_in_quad)
  return(final_matrix_for_angles_in_quad)
}

# Function to categorize angles into quadrants
categorize_quadrant <- function(angle) {
  if (angle >= 0 && angle < 90) {
    return("Q1")
  } else if (angle >= 90 && angle < 180) {
    return("Q2")
  } else if (angle >= 180 && angle < 270) {
    return("Q3")
  } else if (angle >= 270 && angle < 360) {
    return("Q4")
  } else {
    return("Unknown")
  }
}

# utility function to find optmial cluster size based on distance metric
find_distance_elbow <- function(wss){
  
  if(length(wss)<3){
    return(1)
  }
  
  # Calculate the line equation
  x1 <- length(wss) # last point
  y1 <- wss[x1]
  x2 <- 1 # 1st point
  y2 <- wss[x2]
  
  # For line Ax + By + C = 0
  A <- y2 - y1
  B <- x1 - x2
  C <- x2*y1 - x1*y2
  
  # Calculate the distance of each point to the line
  distances <- numeric(length(wss))
  
  for (i in 1:length(wss)) {
    x <- i
    y <- wss[i]
    
    # Distance from point to line formula |Ax+By+C|/sqrt(A^2+B^2)
    distances[i] <- abs(A*x + B*y + C) / sqrt(A^2 + B^2)
  }
  
  which.max(distances)
}

# utility for wss score computation
computeWss<- function(matrix_of_angles){
  wss <- function(k) {
    kmeans(matrix_of_angles, k, nstart = 10 )$tot.withinss
  }
  
  kmax = nrow(matrix_of_angles) - 1
  k.values = NULL
  if (kmax ==2){
    k.values <- c(2)
  }else {
    k.values <- 2:nrow(matrix_of_angles) -1
  }
  
  # extract wss
  wss_values <- sapply(k.values, wss)
  return(wss_values)
  
}

#plot single quadrant
plotQuadrant<-function(quad,wss_values){
  path = getwd()
  pathpdf = paste(path,"plot_for_clusters_in",sep="/")
  outputPath = paste0(pathpdf,sep="_",quad,".pdf")
  pdf(outputPath,width=20,height = 15)
  plot(1:length(wss_values), wss_values,
       type="b", pch = 19, frame = FALSE, 
       xlab="Number of clusters K",
       ylab="Total within-clusters sum of squares")
  title("Elbow method-optimal number of clusters")
  dev.off()
}

computeEffectScores <- function(weightedNetwork,PLSoutputpathway,clusters,pathwayname){
  message("\n Calculating effectscores per pathway...")
  effectScoreperPathway <- list()
  trsnamelist <- list()
  colnames(clusters) = c("TRs","Pathway","cluster")
  ##identify pathways for which clusters are provided
  pathwayCluster <- as.matrix(unique(clusters$Pathway))
  for(i in 1:nrow(pathwayCluster)){
    pathwayName <- pathwayCluster[i,]
    ## determine TRs associated to each pathway
    # TRsperPathway <-  PLSoutput[[pathwayName]]
    TRsperPathway<- PLSoutputpathway
    ## determine weights associated to edge between TRs and a pathway
    weightedpathway = weightedNetwork[which(weightedNetwork$Pathway==pathwayName),]
    weightsperTR <- weightedpathway[match(TRsperPathway$TRs,weightedpathway$TRs),]
    weightedSubnetwork <- cbind(TRsperPathway,weightsperTR$edge_weights)
    maxNoofClusters <- max(as.numeric(unique(clusters$cluster)))
    meanValuesforClusters <- matrix(0,maxNoofClusters,4)
    clusterName <- matrix(0,maxNoofClusters,2)
    ## determine mean angle, correlation with LV1, LV2 and edge weight for each cluster
    for(j in 1:maxNoofClusters)
    {
      clusterMat <- as.data.frame(clusters[which(clusters$cluster==j),])
      TRsname <- clusterMat$TRs
      trsnamelist[[j]] = TRsname
      valuesforCluster <- weightedSubnetwork[match(clusterMat$TRs,weightedSubnetwork$TRs),]
      meanAngle <- mean(abs(as.numeric(valuesforCluster$angle_Wrt_pathway)))
      meanValuesforClusters[j,1] <- meanAngle
      meanCorrelationLV1 <- mean(as.numeric(valuesforCluster$`Correlation with LV1`))
      meanValuesforClusters[j,2] <-  meanCorrelationLV1
      meanCorrelationLV2 <- mean(as.numeric(valuesforCluster$`Correlation with LV2`))
      meanValuesforClusters[j,3] <-  meanCorrelationLV2
      meanWeight <- mean(valuesforCluster$`weightsperTR$edge_weights`)
      meanValuesforClusters[j,4] <- meanWeight
      clstName <- paste0("cluster",j)
      clusterName[j,1] <- clstName
    }
    
    rownames(meanValuesforClusters) <- clusterName[,1]
    colnames(meanValuesforClusters) <- c("Angle","correlation_with_LV1","correlation_with_LV2","weight")
    ## calculate effect scores
    effectScores <- EffectScorescalc(meanValuesforClusters)
    effectScoreperPathway[[pathwayName]] <- list(effectScores,trsnamelist)
    
    
  }
  message("finished calculating effect scores")
  effectScoreperPathway = list(pathwayname = pathwayname, cluster = names(effectScoreperPathway[[1]][[1]]),effectscores = effectScoreperPathway[[1]][[1]],TRs = effectScoreperPathway[[1]][[2]])
  return(effectScoreperPathway)
}

## Calculate Effect scores per Pathway
EffectScorescalc <- function(clusterValues){
  
  # first column angle normal rank
  clusterValues[,1] <- dplyr::dense_rank(clusterValues[,1])
  
  # reverse rank other 3 columns
  clusterValues[,2:4] <- apply(clusterValues[,2:4], 2, function(x) dplyr::dense_rank(-x))
  
  # Compute geometric mean for each row
  effectScores <- apply(clusterValues, 1, psych::geometric.mean)
  
  return(effectScores)
  
}


#' Generate clusters per quadrant for a single pathway
#' 
#'
#' @param anglematclusters  
#' @param PLSoutputpathway 
#'
#' @return clustermat
generateClusters = function(anglematclusters,PLSoutputpathway){
  

  merged_cluster = matrix(data= NA,0,2)
  
  quads = names(anglematclusters$angleMat)

  for(quad in quads)
  {
    kval = 1
    
    quadanglemat = anglematclusters$angleMat[[quad]]
    
    if(nrow(quadanglemat)>1){
      # wss
      wss_values<- computeWss(quadanglemat)
      
      # optimal size
      kval <- find_distance_elbow(wss_values)
      
      # plot quadrant
      plotQuadrant(quad,wss_values)
    }
 
    
    if(kval>1){
      
      anglemat = as.matrix(as.numeric(PLSoutputpathway$angle_Wrt_pathway))
      rownames(anglemat) = PLSoutputpathway$TRs
      anglesforquad = anglemat[match(rownames(quadanglemat),rownames(anglemat)),]
      anglesforquad = as.matrix(anglesforquad)
      
      dist_mat <- dist(anglesforquad, method = 'euclidean')
      hclust_avg <- hclust(dist_mat, method = 'average')
      cut_avg <- as.matrix(cutree(hclust_avg, k = kval))
      
     # kmeanclstr = kmeans(quadanglemat,kval)
     # clusters_forquad = as.matrix(kmeanclstr$cluster)
      
      clusters_forquad = as.data.frame(cut_avg)
      clusters_forquad$quadno = quad
      merged_cluster =    rbind(merged_cluster,clusters_forquad)
      
    } else if(kval ==1){
      quadanglemat = anglematclusters$angleMat[[quad]]
      clusters_forquad = matrix(0,nrow(quadanglemat),2)
      rownames(clusters_forquad) = rownames(quadanglemat)
      colnames(clusters_forquad) = c("V1","quadno")
      clusters_forquad[,1]=1
      clusters_forquad[,2] = quad
      merged_cluster = rbind(merged_cluster,clusters_forquad)
    }
  }
  
  merged_cluster = as.data.frame(merged_cluster)
  merged_cluster$rownames = rownames(merged_cluster)
  # add valiation for to check columns and data exist 
  df_unique <- merged_cluster %>%
    distinct(quadno, V1) %>%
    arrange(quadno, V1) %>%
    mutate(unique_number = row_number())
  
  #View(merged_cluster)
  # Merge back to the original data frame
  merged_cluster_nw <- merge(merged_cluster, df_unique, by = c("quadno", "V1"))
  
  rownames(merged_cluster_nw) = merged_cluster_nw$rownames
 # View(merged_cluster_nw)
  
  merged_cluster_nw = merged_cluster_nw[,c("quadno","V1","unique_number")]
  merged_Cluster_with_PLSpathway = merged_cluster_nw[match(PLSoutputpathway$TRs,rownames(merged_cluster_nw)),]
  clustermat = cbind(PLSoutputpathway$TRs,PLSoutputpathway$Pathway,merged_Cluster_with_PLSpathway$unique_number)
  colnames(clustermat) = c("TRS","Pathway","clusters")
  #View(clustermat)
  clustermat = as.data.frame(clustermat)
  return(clustermat)
}


#' @export
computeScores <- function(weightedNetwork,subnetrwklist,pathwayname,TRActivity,pathwayActivity){
  
  pathwayname = as.character(pathwayname)
  subnetrwk = subnetrwklist[[pathwayname]]
  if(length(subnetrwk)==0){
    stop("pathway is not in the subnetwork")
  }
  PLSoutputpathway <- computePLS(subnetrwk,TRActivity,pathwayActivity,weightedNetwork)
  
  PLSoutputpathway <- PLSoutputpathway[[1]]
  #  make quads
  anglematclusters<- makeAngleClusters(PLSoutputpathway)
  
  # genrate  clusters
  # computewss,optmalsize,plot
  clusters <- generateClusters(anglematclusters,PLSoutputpathway)
  
  # compute scores 
  priorityscores <- computeEffectScores(weightedNetwork,PLSoutputpathway,clusters,pathwayname)
  priorityscoresdf <- data.frame(clusters = priorityscores$cluster,
                    TRs = as.character(priorityscores$TRs,quote=FALSE),effectscore= priorityscores$effectscores)
  
  priorityscoresdf$TRs =gsub("^c\\(|\\)$", "", priorityscoresdf$TRs)
  priorityscoresdf$TRs <- gsub("\"", "", priorityscoresdf$TRs)
  
  return(priorityscoresdf)
  
}