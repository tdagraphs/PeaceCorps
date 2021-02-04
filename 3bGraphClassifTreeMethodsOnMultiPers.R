rm(list = ls())
library(imager)
library(ggplot2)
library(plyr)
library(dplyr)
library(randomForest)
library(stringr)
library(xgboost)
library(igraph)
library(mlbench)
library(caret)
library(e1071)



#Metric compare model is Accuracy
metric <- "Accuracy"
set.seed(123)
trainSize <- 0.8

inputFile <- "PersistenceGrid.txt"
features <- c("betweenness","eccentricity","closeness","degree","authority")


inputPath <-"/home/jupiter/GraphML/MultiFeature/Grids/"
outputPath<-"/home/jupiter/GraphML/SingleFeature/Predictions/"
labelPath<-"/home/jupiter/GraphML/"

grid1Size=11
grid2Size=12

multiPersClassification<-function(dataset,dataAlias,feature1,feature2, whichBettiToUse=2){
  whichGridFile<-paste0(inputPath,dataAlias,feature1,feature2,inputFile)
  
  
  if(!file.exists(whichGridFile)){
    whichGridFile<-paste0(inputPath,dataAlias,feature2,feature1,inputFile)
  }
  if(!file.exists(whichGridFile)) return()
  data <- read.table(whichGridFile, header = F, sep = "\t")
  colnames(data) <- c("graphId", "dataset", "grid1", "val1",
                      "f1","f2",sprintf("grid2_%d", 1:grid2Size),"betti")
  # we will have betti numbers from a grid starting at he 7th index
  sigStartCol=7
  sigEndCol=sigStartCol+grid2Size-1
  
  graphSignatureVectors <- data.frame()
  f1gridVals=1:(2*grid1Size)
  for(gId in unique(data$graphId)){
    graphImage<-array(0,dim=c(grid1Size,2*grid2Size))
    gData<-data[data$graphId==gId,]
    if(nrow(gData)!=2*grid1Size){
      message("ERROR: Grid 1 values are missing!")
    }
    
    for(thisRowIndex in f1gridVals){
      row = gData[thisRowIndex,]
      # index of the grid, out of grid2Size total
      grid1=row[["grid1"]]
      # betti number of the row, 0 or 1
      betti = row[["betti"]]
      # betti values of the secondary grid
      grid2_vals=as.numeric(row[,sigStartCol:sigEndCol])
      # we need to compute where to insert these values.
      # first, we need to find which column of gridSizeXgridSize matrix
      g2start = (1+betti*grid2Size)
      g2end = g2start+grid2Size-1
      graphImage[grid1,g2start:g2end]=(grid2_vals) 
    }
    flattenedArray = c(gId,c(graphImage))
    graphSignatureVectors <- rbind.data.frame(graphSignatureVectors, flattenedArray)
  }
  
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
  colnames(graphSignatureVectors)<-c("graphId",sprintf("grid_%d", 1:(2*grid1Size*grid2Size)))
  labeledData <- merge(labels, graphSignatureVectors, by = "graphId")
  labeledData$graphId<-NULL
  labeledData$label <- as.factor(labeledData$label)
  
  trainingDataSize <- floor(nrow(labeledData) * trainSize)
  # Generate a random sample of "trainingDataSize" indexes
  indexes <- sample(seq_len(nrow(labeledData)), size = trainingDataSize)
  # Divide the data to the training and test sets
  training <- labeledData[indexes,]
  test <- labeledData[-indexes,]
  # use betti 0 only
  if(whichBettiToUse==0){
    training=training[,1:(1+grid1Size*grid2Size)]
    test=test[,1:(1+grid1Size*grid2Size)]
  }else if(whichBettiToUse==1){
    training[,2:(1+grid1Size*grid2Size)]<-NULL
    test[,2:(1+grid1Size*grid2Size)]<-NULL
  }
  # create a random forest classifier
  
  numcol <- sqrt(ncol(training))
  # we select mtry values around the sqrt(numcols). This is standard selection
  tunegrid <- expand.grid(.mtry=seq(11,numcol+6,by=3))
  control <- trainControl(method="repeatedcv", number=2, repeats=1, search="grid")
  
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
  message(dataAlias,"CrossValidation\t",feature1,"\t",feature2,"\trf\t", trainMaxAccuracy,"\t",testRfAccuracy)
  
}


 
for(i in 1:30){
   
  
  multiPersClassification(dataset="IMDB-MULTI/IMDB-MULTI.",dataAlias ="IMDBMulti", feature1="betweenness",feature2="closeness")
  multiPersClassification(dataset="IMDB-BINARY/IMDB-BINARY.",dataAlias ="IMDBBinary", feature1="degree",feature2="betweenness")
  multiPersClassification(dataset="REDDIT-BINARY/REDDIT-BINARY.",dataAlias ="RedditBinary", feature1="degree",feature2="closeness")
  multiPersClassification(dataset="proteins/proteins.",dataAlias ="Protein", feature1="betweenness",feature2="closeness")
  multiPersClassification("BZR/BZR.","BZR",feature1="degree",feature2="closeness")
  multiPersClassification("NCI1/NCI1.","NCI1",feature1="degree",feature2="betweenness")
  multiPersClassification("REDDIT-MULTI-5K/REDDIT-MULTI-5K.","REDDIT5K",feature1="closeness",feature2="betweenness")
  multiPersClassification("COX2/COX2.","COX2",feature1="degree",feature2="betweenness")
  multiPersClassification("DHFR/DHFR.","DHFR",feature1="closeness",feature2="betweenness")
  multiPersClassification("FRANKENSTEIN/FRANKENSTEIN.","FRANKENSTEIN",feature1="degree",feature2="betweenness")
} 


 
