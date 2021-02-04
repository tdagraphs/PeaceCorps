# Objective : Run Random Forest and XGBoost classifiers on betti signature to
# predict graph labels.
# requires results from 2aPDSignatureExtractionForDiscreteFeatures.R
# and 2bPDSignatureExtractionForContinuousFeatures.R
# Created by: **
# Created on: 2020-12-18
rm(list = ls())
library(plyr)
library(dplyr)
library(randomForest)
library(ggplot2)
library(stringr)
library(xgboost)
library(igraph)

useSubLevel <- TRUE
subSignatureFile <- "baselineSub.txt"
dataPath<-"Documents/GraphML/"
dataPath2<-"PeaceCorps/"
nodeFeatures <- c("betweenness","closeness","authority")#,"eccentricity","degree"

for(nodeFeature in nodeFeatures){
  whichSignatureFile<-paste0(dataPath2,"baseline/",nodeFeature,subSignatureFile)
  
  data <- read.table(whichSignatureFile, header = F, sep = "\t")
  data <- na.omit(data)
  colnames(data) <- c("dataset", "graphId", "betti", "bettisignature")
  
  # does every graph have both betti 0 and 1?
  # at least some must have
  sanityCheck <- ddply(data, .(dataset, graphId), summarize, l = length(unique(betti)))
  table(sanityCheck$l)
  
  trainSize <- 0.8
  
  for (dataset in unique(data$dataset)) {
    datasetData <- data[data$dataset == dataset,]
    idFile <- paste0(dataPath,dataset, "graph_labels")
    labels <- read.table(idFile, quote = "\"", comment.char = "", sep = ",")
    colnames(labels) <- "label"
    # add graph id to the label data
    labels$graphId <- seq(1, nrow(labels))
    # we need label conversion because some ML algs require labels to be numeric and start from 0
    lblIndex<-0
    labels2<-as.data.frame(sort(unique(labels$label)))
    labels2$id<-seq(0,(nrow(labels2)-1))
    colnames(labels2)<-c("oldlabel","newlabel")
    for(ind in seq(0,nrow(labels))){
      oldl<-labels[ind,"label"]
      labels[ind,"label"]<-labels2[labels2$oldlabel==oldl,]$newlabel
    }
    for (bettiNumber in unique(datasetData$betti)) {
      datasetSingleBettiData <- datasetData[datasetData$betti == bettiNumber,]
      if(length(unique(datasetSingleBettiData$bettisignature))==1){
        message(dataset,nodeFeature,bettiNumber,": All data has a single betti signature. Classifier cannot be built on such data.")
      }else{      
        # compute the length of the signature vector that we need to create for graphs of this dataset
        vectorLength <- 1 + max(str_count(datasetSingleBettiData$bettisignature, " "))
      
        graphSignatureVectors <- data.frame()
        for (row in seq(seq_len(nrow(datasetSingleBettiData)))) {
          graphId <- (datasetSingleBettiData[row,]$graphId)
          bfunct <- as.character(datasetSingleBettiData[row,]$bettisignature)
          value <- strsplit(bfunct, split = " ")[[1]]
          # normalizing the signature
          value <- as.integer(value)
          #meanOfValues<-mean(value)
          #stdOfValues<-sd(value)
          #value<-(value-meanOfValues)/stdOfValues
          # convert signature array to our final vector
          signatureArray <- value[1:vectorLength]
          # normalization by max value
          maxVal<-max(signatureArray)
          if(maxVal>0){
          signatureArray <- signatureArray/max(signatureArray)
          }
          formattedGraphSignature <- cbind.data.frame(name = dataset,
                                                      graphId = as.integer(graphId))
          for (ind in seq(1:vectorLength)) {
            formattedGraphSignature <- cbind.data.frame(formattedGraphSignature, signatureArray[ind])
          }
          # add graph signature to the dataset
          graphSignatureVectors <- rbind.data.frame(graphSignatureVectors, formattedGraphSignature)
          
        }
        
        # Classification starts at this point!
        # we do not need dataset name in classification, do not ask me why I put it in the 1st place.
        graphSignatureVectors$name <- NULL
        # assign labels to graphs
        labeledData <- merge(labels, graphSignatureVectors, by = "graphId")
        colnames(labeledData) <- c("graphId", "label", sprintf("betti%s", seq(1:vectorLength)))
        labeledData$graphId<-NULL
        labeledData$label <- as.factor(labeledData$label)
        
        trainingDataSize <- floor(nrow(labeledData) * trainSize)
        # Generate a random sample of "data_set_size" indexes
        indexes <- sample(seq_len(nrow(labeledData)), size = trainingDataSize)
        # Divide the data to the training and test sets
        training <- labeledData[indexes,]
        test <- labeledData[-indexes,]
        # create a random forest classifier
        rf <- randomForest(formula = label ~ ., data = training, ntree = 500, mtry = sqrt(vectorLength), importance = TRUE)
        varImpPlot(rf,type=2)
        # predict labels of test set graphs
        predictedLabels <- predict(rf, test)
        rfAccuracy <- sum(test$label == predictedLabels) / nrow(test)
        message(dataset, "\t",nodeFeature,"\t",bettiNumber,"\trf\t", rfAccuracy)
        
        # xgboost classifier
        # just to prevent data contamination, we recreate trainign and test datasets with the same indices
        # as those used in the random forest
        labeledData <- merge(labels, graphSignatureVectors, by = "graphId")
        colnames(labeledData) <- c("graphId", "label", sprintf("betti%s", seq(1:vectorLength)))
        
        labeledData$label<-as.integer(labeledData$label)
        training <- labeledData[indexes,]
        test <- labeledData[-indexes,]
        numClasses <- length(unique(training$label))
        
        trainMatrix <- as.matrix(training[,3:length(training)])
        
        
        bst <- xgboost(data = trainMatrix,label = training$label,nrounds = 25,
                       objective = "multi:softprob", num_class = numClasses, verbose=F)
        pred <- predict(bst, as.matrix(test[,3:length(test)]))
        
        pred <- matrix(pred, ncol=numClasses, byrow=TRUE)
        # convert the probabilities to softmax labels
        pred_labels <- max.col(pred) - 1
        xgboostAccuracy<- sum(pred_labels == test$label)/length(test$label)
        message(dataset, "\t",nodeFeature,"\t",bettiNumber,"\txgboost\t", xgboostAccuracy)
        
        
      }
    }
  }
}
