# Title     : GraphTreeMethods.R
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
library(mlbench)
library(caret)
library(MLmetrics)


inputPath="/home/jupiter/GraphML/SingleFeature/SignaturesMinMax/"
outputPath = "/home/jupiter/GraphML/SingleFeature/Predictions/"
labelPath<-"/home/jupiter/GraphML/"

 
trainSize = 0.8
kNeighbor=5

classifyWithSaw<-function(dataset ,dataAlias , feature){
  whichSignatureFile<-paste0(inputPath,dataAlias,feature,inputFile)
  if (!file.exists(whichSignatureFile)) {
   message(whichSignatureFile," signature does not exist. ")
    return(-1);
  }
  data <- read.table(whichSignatureFile, header = F, sep = "\t")
  data <- na.omit(data)
  colnames(data) <- c("graphId", "betti", "bettisignature")
  idFile <- paste0(labelPath,dataset, "graph_labels")
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
   
  for(bettiNumber in c(0,1)){
    datasetSingleBettiData <- data[data$betti == bettiNumber,]
    if(length(unique(datasetSingleBettiData$bettisignature))==1){
      message(dataset,feature,bettiNumber,": All data has a single betti signature. Classifier cannot be built on such data.")
    }else{      
      # compute the length of the signature vector that we need to create for graphs of this dataset
      vectorLength <- 1 + max(str_count(datasetSingleBettiData$bettisignature, " "))
      
      graphSignatureVectors <- data.frame()
      for (row in seq(seq_len(nrow(datasetSingleBettiData)))) {
        graphId <- (datasetSingleBettiData[row,]$graphId)
        bfunct <- as.character(datasetSingleBettiData[row,]$bettisignature)
        value <- strsplit(bfunct, split = " ")[[1]]
        value <- as.integer(value)
        if(sum(value)==0) {
          # this should happen sometimes with betti1
          next;
        }
        # convert signature array to our final vector
        signatureArray <- value[1:vectorLength]
        
        formattedGraphSignature <- cbind.data.frame(name = dataAlias,
                                                    graphId = graphId)
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
      numcol <- sqrt(ncol(labeledData))
      # we select mtry values around the sqrt(numcols). This is standard selection
      tunegrid <- expand.grid(.mtry=seq((numcol-3),(numcol+3),by=1))
      control <- trainControl(method="repeatedcv", number=10, repeats=1, search="grid")
      
      rf_default <- train(label~., 
                          data=training, 
                          method='rf', 
                          metric='Accuracy', 
                          tunegrid=tunegrid,
                          trControl=control)
      trainMaxAccuracy= max(rf_default$results$Accuracy)
      testLabels <- predict(rf_default, test)
      # compare predicted outcome and true outcome
      testRfAccuracy <- sum(test$label == testLabels) / nrow(test)
      message(dataAlias,"CrossValidation\t",feature,"\t",bettiNumber,
              "\trf\t", trainMaxAccuracy,"\t",testRfAccuracy,"\t",nrow(training))
      
      # svm classifier
      if(FALSE){
      train_control <- trainControl(method="repeatedcv", number=10, repeats=1)
      svm1 <- train(label ~., data = training, method = "svmRadial", trControl = train_control,  preProcess = c("center","scale"))
      predictedTest <- predict(svm1,test)
      trainingSvmAcc = svm1$results$Accuracy
      testSvmAcc= Accuracy(predictedTest,test$label)
      message(dataAlias,"CrossValidation\t",feature,"\t",bettiNumber,
              "\tsvm\t", trainingSvmAcc,"\t",testSvmAcc,"\t",nrow(training))
      }
      # knn classifier
      if(FALSE){
      matr=dtwDist(data.matrix(test[-1]),data.matrix(training[-1]))
      pred = array(0)
      for( i in 1:nrow(matr)){
       n<- (names(tail(sort.int(matr[i,],decreasing = T), kNeighbor)))
       label = names(which.max(table(training[n,"label"])))
       pred[i]=label
      }
      testKnnAcc= Accuracy(pred,test$label)
      message(dataAlias,"CrossValidation\t",feature,"\t",bettiNumber,
              "\tknn\t", testKnnAcc,"\t",nrow(training))
      }
      
      if(FALSE){
        #xgboost cross validation
        cv.ctrl <- trainControl(method = "repeatedcv", repeats = 3,number = 10, 
                                #summaryFunction = twoClassSummary,
                                classProbs = TRUE,
                                allowParallel=T)
        
        xgb.grid <- expand.grid(nrounds = 1000,
                                eta = c(0.01,0.05,0.1),
                                colsample_bytree = 1,
                                min_child_weight = 100,
                                subsample = 1,
                                gamma=1,
                                max_depth = c(2,4,6,8,10,14)
        )
        set.seed(135)
        # caret library expects data labels to not start with numbers.
        # We will change labels to start with letter P
        training2<-training
        training2$label<-paste0("P",training2$label)
        xgb_tune <-train(x=training2[,which(names(training2) != "label")],
                         y=training2$label,
                         method="xgbTree",
                         trControl=cv.ctrl,
                         tuneGrid=xgb.grid,
                         verbose=T,
                         metric="Kappa",
                         nthread =5
        )
        test2<-test
        test2$label<-paste0("P",test2$label)
        pred <- predict(xgb_tune, as.matrix(test2[,which(names(test2) != "label")]))
        numClasses <- length(unique(training2$label))
        pred <- matrix(pred, ncol=numClasses, byrow=TRUE)
        # convert the probabilities to softmax labels
        xgboostCvAccuracy<- sum(pred == test2$label)/length(test2$label)
        message(dataAlias, "\t",feature,"\t",bettiNumber,"\txgboost\t", xgboostCvAccuracy)
        
        # xgboost classifier
        
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
        message(dataAlias, "\t",feature,"\t",bettiNumber,"\txgboost\t", xgboostAccuracy,"\t",xgboostCvAccuracy)
      }
    }
    
  }
}

features <- c("betweenness","closeness","degree")#"eccentricity","authority"
features2<-c("ricci","forman")
 
for(i in 1:30){
   if(TRUE){
	classifyWithSaw(dataset="REDDIT-MULTI-5K/REDDIT-MULTI-5K.",dataAlias ="REDDIT5K", feature="closeness")
	classifyWithSaw(dataset="REDDIT-MULTI-5K/REDDIT-MULTI-5K.",dataAlias ="REDDIT5K", feature="betweenness") 
 
    classifyWithSaw(dataset="IMDB-MULTI/IMDB-MULTI.",dataAlias ="IMDBMulti", feature="closeness")
    classifyWithSaw(dataset="IMDB-MULTI/IMDB-MULTI.",dataAlias ="IMDBMulti", feature="degree")
    classifyWithSaw(dataset="IMDB-BINARY/IMDB-BINARY.",dataAlias ="IMDBBinary", feature="closeness")
    classifyWithSaw(dataset="IMDB-BINARY/IMDB-BINARY.",dataAlias ="IMDBBinary", feature="betweenness")
    classifyWithSaw(dataset="REDDIT-BINARY/REDDIT-BINARY.",dataAlias ="RedditBinary", feature="degree")
    classifyWithSaw(dataset="REDDIT-BINARY/REDDIT-BINARY.",dataAlias ="RedditBinary", feature="betweenness")
    classifyWithSaw(dataset="proteins/proteins.",dataAlias ="Protein", feature="ricci" )
    classifyWithSaw(dataset="proteins/proteins.",dataAlias ="Protein", feature="forman" )
    classifyWithSaw(dataset="proteins/proteins.",dataAlias ="Protein", feature="betweenness" )
    classifyWithSaw(dataset="BZR/BZR.",dataAlias ="BZR", feature="closeness")
    classifyWithSaw(dataset="BZR/BZR.",dataAlias ="BZR", feature="degree")
    classifyWithSaw(dataset="NCI1/NCI1.",dataAlias ="NCI1", feature="closeness")
    
    classifyWithSaw(dataset="COX2/COX2.",dataAlias ="COX2", feature="ricci")
    classifyWithSaw(dataset="COX2/COX2.",dataAlias ="COX2", feature="betweenness")
    classifyWithSaw(dataset="DHFR/DHFR.",dataAlias ="DHFR", feature="betweenness")
    classifyWithSaw(dataset="DHFR/DHFR.",dataAlias ="DHFR", feature="closeness")
    classifyWithSaw(dataset="FRANKENSTEIN/FRANKENSTEIN.",dataAlias ="FRANKENSTEIN", feature="degree")
    classifyWithSaw(dataset="FRANKENSTEIN/FRANKENSTEIN.",dataAlias ="FRANKENSTEIN", feature="betweenness")
  }   
}
