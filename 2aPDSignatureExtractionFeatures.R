# Title     : PDSignatureExtractionForContinuousFeatures
# Objective : Extract graph betti signatures from persistence diagrams
# requires results from GraphPersistenceDiagramExtraction.R
# Created by: **
# Created on: 2020-12-18
rm(list = ls())
library(dplyr)
library(plyr)
library(igraph)


dataPath="PeaceCorps/SingleFeature/Filtrations/"
outputPath="PeaceCorps/SingleFeature/SignaturesMinMax/"

dataPath <- "/home/jupiter/GraphML/SingleFeature/Filtrations/"
outputPath="/home/jupiter/GraphML/SingleFeature/SignaturesMinMax/"



inputFile <- "PDSubFiltration.txt"
outputFile <-  "Signature.txt"


sigComputation<-function(graphData,bettiNum,sLen){
  graphBettiData = graphData[graphData$betti == bettiNum,]
  signatureLength<-3*sLen
  if(nrow(graphBettiData)==0){
    return(rep(0,signatureLength))
  }
  
  # Normalize birth and death threshod values to an interval between 0 and sLen
  maxThreshold <- max(graphBettiData[,c("birth","death")],na.rm=T)
  minThreshold <- min(graphBettiData[,c("birth","death")],na.rm=T)
  # signatures from min to max
  graphBettiData$birth<-as.integer(sLen*(graphBettiData$birth-minThreshold)/(maxThreshold-minThreshold))
  graphBettiData$death<-as.integer(sLen*(graphBettiData$death-minThreshold)/(maxThreshold-minThreshold))
  # signatures from 0 to max
  #graphBettiData$birth<-as.integer(sLen*(graphBettiData$birth)/(maxThreshold))
  #graphBettiData$death<-as.integer(sLen*(graphBettiData$death)/(maxThreshold))
  
  if(max(graphBettiData$death,na.rm = TRUE)>sLen||min(graphBettiData$birth,na.rm = TRUE)<0){
    message("Warning: ", nodeFeature," values for ",dataAlias," is not within [0,sLen]")
  }
  starts <- unique(graphBettiData$birth)
  # count how many births occur at what threshold
  birthData <- ddply(graphBettiData, .(birth), summarise, bs = length(birth))
  # count how many deaths occur at what threshold
  deathData <- ddply(graphBettiData, .(death), summarise, bs = length(birth))
  deaths <- unique(deathData$death)
  nrowResult <- nrow(birthData)
  
  thresholdValues <- sort(unique(c(starts, deaths)))
  thresholdValues <- thresholdValues[!is.na(thresholdValues)]
  bettiSignatureArray = array(c(0))
  
  # create the betti step signature.
  # we use an index scheme of 3* threshold
  # index 3*threhold-1 holds bettis that exist before deaths at the threshold
  # index 3*threhold holds bettis that remain after deaths at the threshold
  # index 3*threshold+1 holds bettis that remain after deaths + those that are newly born at the threshold
  if (nrowResult > 0) {
    results<-data.frame()
    
    latestCountValue <- 0
    minThreshold=min(starts, na.rm = TRUE)
    for (threshold in thresholdValues) {
      birthCount <- birthData[birthData$birth == threshold,]$bs[1]
      deathCount <- deathData[deathData$death == threshold,]$bs[1]
      if(is.na(deathCount)) deathCount=0
      if(is.na(birthCount)) birthCount=0
      if(threshold==0){
        bettiSignatureArray[1] <- birthCount-deathCount
        bettiSignatureArray[2] <- birthCount
        latestCountValue = birthCount
      }else{
        bettiSignatureArray[(threshold * 3) - 1] <- latestCountValue
        bettiSignatureArray[threshold * 3] <- latestCountValue - deathCount
        if(threshold!=sLen){
          bettiSignatureArray[threshold * 3 + 1] <- bettiSignatureArray[threshold * 3]+birthCount
          latestCountValue = bettiSignatureArray[threshold * 3 + 1]
        }
      }
    }
    # in our 3* thresholds scheme, some values will remain NA because
    # no new bettis are born or die around some thresholds.
    # if all bettis die, the signature will end in 0
    # otherwise, if some bettis survive, the signature will end in a value>0
    for (i in seq(1, (signatureLength))) {
      if (is.na(bettiSignatureArray[i])) {
        if(i==1){
          bettiSignatureArray[i]=0
        }else {
          bettiSignatureArray[i] <- bettiSignatureArray[i - 1]
        }
      } 
    }
  }
  # write signatures to the file
  return(bettiSignatureArray)
}


computeSaw<-function(dataset,dataAlias,feature,sLen){
  signatureLength<-3*sLen
  whichInputFile<-paste0(dataPath,feature,dataAlias,inputFile)
  if (!file.exists(whichInputFile)) {
     message(dataAlias," persistence diagram file does not exist!")
    return(-1);
  }
  data <- read.table(file = whichInputFile,header = F, sep = "\t")
  colnames(data) <- c("betti", "birth", "death", "graphId", "dataAlias")
  
  
  # does every graph have noth betti 0 and 1?
  #sanityCheck <- ddply(data, .(graphId), summarize, l = length(unique(betti)))
  #table(sanityCheck$l)
  
  whichOutputFile<-paste0(outputPath,dataAlias,feature,outputFile)
  if (file.exists(whichOutputFile)) {
    #Delete file if it exists
    file.remove(whichOutputFile)
  }
  
  # format persistent shapes (end at infinity)
  data$death[!is.finite(data$death)] <- NA
  graphs<-unique(data$graphId)
  message("Processing ",dataAlias," for feature ",feature,". Num of Graphs: ",length(graphs))
  
  bettis<-c(0,1)
  results<-data.frame()
  for(gId in graphs){
    graphData <- data[data$graphId == gId,]
    for (bettiNum in bettis) {
      bettiSignature<-sigComputation(graphData,bettiNum,sLen)
      bettiInfo = cbind(graphId=gId, bettiNum, paste(bettiSignature, collapse = " "))
      results<-rbind(results,bettiInfo)
    }
    
  }
  write.table(results, file = whichOutputFile, sep = "\t", row.names = FALSE, col.names = FALSE, append = T, quote = FALSE)
  
  
}


nodeFeatures <- c("degree","betweenness","closeness")
nodeFeatures2<-c("ricci","forman")
for(f in c(nodeFeatures)){
  computeSaw(dataset="REDDIT-MULTI-5K/REDDIT-MULTI-5K.",dataAlias ="REDDIT5K", feature=f, sLen=100)
  
  
  if(TRUE){
    computeSaw(dataset="REDDIT-BINARY/REDDIT-BINARY.",dataAlias ="RedditBinary", feature=f, sLen=100)
    computeSaw(dataset="IMDB-MULTI/IMDB-MULTI.",dataAlias ="IMDBMulti", feature=f, sLen=100)
    computeSaw(dataset="IMDB-BINARY/IMDB-BINARY.",dataAlias ="IMDBBinary", feature=f, sLen=100)
  computeSaw(dataset="ENZYMES/ENZYMES.",dataAlias ="Enzyme", feature=f, sLen=100)
  computeSaw(dataset="BZR/BZR.",dataAlias ="BZR", feature=f, sLen=100)
  #
  computeSaw(dataset="COX2/COX2.",dataAlias ="COX2", feature=f, sLen=100)
  computeSaw(dataset="DHFR/DHFR.",dataAlias ="DHFR", feature=f, sLen=100)
  computeSaw(dataset="NCI1/NCI1.",dataAlias ="NCI1", feature=f, sLen=100)
  computeSaw(dataset="FRANKENSTEIN/FRANKENSTEIN.",dataAlias ="FRANKENSTEIN", feature=f, sLen=100)
  
  
  #computeSaw(dataset="ENZYMES/ENZYMES.",dataAlias ="Enzyme", feature=f, sLen=100)
  computeSaw(dataset="proteins/proteins.",dataAlias ="Protein", feature=f, sLen=100)
 
  }
}

source("3aGraphClassifTreeMethodsOnSingleSignature.R")


