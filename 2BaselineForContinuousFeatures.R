# Objective : Extract graph betti signatures from persistence diagrams
# requires results from GraphPersistenceDiagramExtraction.R
# Created by: **
# Created on: 2020-12-18
rm(list = ls())
library(dplyr)
library(plyr)
library(igraph)



subresultsFile <- "PDSubFiltration.txt"

dataPath="IdeaProjects/PeaceCorps/"

# node features whose signatures wil be prepared on graphs
nodeFeatures <- c("betweenness","authority","closeness")
baselineSignatureFile <- "baselineSub.txt"

signatureLength<-100

for(nodeFeature in nodeFeatures){
  whichInputFile<-paste0(dataPath,nodeFeature,subresultsFile)
  
  data <- read.table(file = whichInputFile,header = F, sep = "\t")
  # the first line contains NA due to creating a data frame
  data <- na.omit(data)
  colnames(data) <- c("betti", "birth", "death", "graphId", "dataset")
  
  
  # does every graph have noth betti 0 and 1?
  sanityCheck <- ddply(data, .(dataset, graphId), summarize, l = length(unique(betti)))
  table(sanityCheck$l)
  
  
  
  
  whichOutputFile<-paste0("baseline/",nodeFeature,baselineSignatureFile)
  if (file.exists(whichOutputFile)) {
    #Delete file if it exists
    file.remove(whichOutputFile)
  }
  
  # format persistent shapes (end at infinity)
  data$death[!is.finite(data$death)] <- NA
  
  resultsTable<-data.frame()
  for (datasetName in unique(data$dataset)) {
    message("Processing ", datasetName, " for feature ",nodeFeature)
    dataOfADataset <- data[data$dataset == datasetName,]
    
    
    # Start extracting graph signatures
    for (graphId in unique(dataOfADataset$graphId)) {
      specificGraphData_ <- dataOfADataset[dataOfADataset$graphId == graphId,]
      for (bettiNum in unique(specificGraphData_$betti)) {
        specificGraphData = specificGraphData_[specificGraphData_$betti == bettiNum,]
        
        # Normalize birth and death threshod values to an interval between 0 and 100
        maxThres <- max(specificGraphData[,c("birth","death")],na.rm=T)
        specificGraphData$birth<-as.integer(100*specificGraphData$birth/maxThres)
        specificGraphData$death<-as.integer(100*specificGraphData$death/maxThres)
        if(max(specificGraphData$death,na.rm = TRUE)>100||min(specificGraphData$birth,na.rm = TRUE)<0){
          message("Warning: ", nodeFeature," values for ",datasetName," is not within [0,100]")
        }
        # normalization ends here
        starts <- unique(specificGraphData$birth)
        # count how many births occur at what threshold
        birthData <- ddply(specificGraphData, .(birth), summarise, bs = length(birth))
        # count how many deaths occur at what threshold
        deathData <- ddply(specificGraphData, .(death), summarise, bs = length(birth))
        deaths <- unique(deathData$death)
        nrowResult <- nrow(birthData)
        
        thresholdValues <- sort(unique(c(starts, deaths)))
        thresholdValues <- thresholdValues[!is.na(thresholdValues)]
        bettiSignatureArray = array(0,dim=100)
        
        # create the betti step signature.
        if (nrowResult > 0) {
          latestCountValue <- 0
          minThreshold=min(starts,na.rm = TRUE)
          
          for (threshold in thresholdValues) {
            birthCount <- birthData[birthData$birth == threshold,]$bs[1]
            deathCount <- deathData[deathData$death == threshold,]$bs[1]
            if(is.na(deathCount)) deathCount=0
            if(is.na(birthCount)) birthCount=0
            
            if(threshold==0){
              latestCountValue = birthCount
            }
            else {
              bettiSignatureArray[threshold] <- latestCountValue - deathCount+ birthCount
              latestCountValue = bettiSignatureArray[threshold]
            }
          }
          # in our thresholds scheme, some values will remain NA because
          # no new bettis are born or die around some thresholds.
          # if all bettis die, the signature will end in 0
          # otherwise, if some bettis survive, the signature will end in a value>0
          for (i in seq(1, (1+signatureLength))) {
            if (is.na(bettiSignatureArray[i])) {
              if(i==1){
                bettiSignatureArray[i]=0
              }else {
                bettiSignatureArray[i] <- bettiSignatureArray[i - 1]
              }
            } 
          }
          
          # write signatures to the file
          bettiSignatureArray = cbind(datasetName, graphId, bettiNum, paste(bettiSignatureArray, collapse = " "))
          #resultsTable<-bind_rows(resultsTable,bettiSignatureArray)
          write.table(bettiSignatureArray, file = whichOutputFile, sep = "\t", row.names = FALSE, col.names = FALSE, append = T, quote = FALSE)
        }
      }
    }
  }
}
