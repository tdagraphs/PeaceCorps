# Title     : GraphPersistenceDiagramExtraction
# Objective : Extract sub level persistence diagrams from graphs
# requires dataset folders from http://networkrepository.com/labeled.php
# Created by: **
# Created on: 2020-12-18
library(igraph)
library(TDA)
library(dplyr)
rm(list = ls())
options(java.parameters = "-Xmx200g")


# 3 Graph ML datasets
dataPath <- "/home/jupiter/GraphML/"
outputPath="/home/jupiter/GraphML/SingleFeature/Filtrations/"


# max dimension of the homological features to be computed. (e.g. 0 for connected components, 1 for connected components and loops, 2 for connected components, loops, voids, etc.)
maxDimension <- 1

# upper limit on the simplex dimension, size of the cliques to find
maxClique <-2
# node count thresholds to ignore a graph
minNodeCount <- 4

outputFile <- "PDSubFiltration.txt"

# if subLevel=FALSE, consider reducing maxNodeCount because power filtration
# does not work with large graphs
maxNodeCount <- 6000

source("nodeFunctions.R")

extractPD<-function(dataset, dataAlias, feature){
  whichOutputFile<-paste0(outputPath,feature,dataAlias,outputFile)
  # Remove earlier result files, if they exist
  if (file.exists(whichOutputFile)) {
    file.remove(whichOutputFile)
  }
  edgeFile <- paste0(dataPath, dataset, "edges")
  graphFile <- paste0(dataPath, dataset, "graph_idx")
  graphIdData <-
    read.table(graphFile,
               quote = "\"",
               comment.char = "",
               sep = ",")
  allData <-
    read.table(edgeFile,
               quote = "\"",
               comment.char = "",
               sep = ",")
  colnames(allData) <- c("source", "target")
  graphs <- unique(graphIdData$V1)
  message("Processing ",dataAlias," for feature ",feature,". Num of Graphs: ",length(graphs))
  results<-data.frame()
  for (graphId in graphs) {
    # locate nodes of the graph
    thisGraphNodes <- which(graphIdData$V1 == graphId)
    # load edges of the graph
    edgeData <-
      allData[allData$source %in% thisGraphNodes |
                allData$target %in% thisGraphNodes, ]
    graph <- graph.data.frame(edgeData, directed = FALSE)
    nodeCount <- vcount(graph)
    
    # if the graph is not too small nor too big, compute filtration
    if (nodeCount > minNodeCount && nodeCount < maxNodeCount) {
      # below we use sub level filtrations.
      
      #sublevel filtration
      # compute node function values 
      tryCatch({
        F.values = computeNodeVals(feature, graph,dataset=dataset,dataPath=dataPath)
      }, error = function(error_condition) {
        message("Error in extracting PDs for ",dataAlias)
        return(-1);
      })
      # for maxClique=3 below means we are adding 0,1,2 simplices (nodes,edges,triangles) to our complex
      cmplx <- cliques(as.undirected(graph), min = 0, max = maxClique)
      # use sublevel=T for sublevel, sublevel=F for superlevel filtration
      # F.values are node values. At these values, node appear in the complex,
      # and their edges to other active nodes also activate.
      FltRips <- funFiltration(FUNvalues = F.values,
                               cmplx = cmplx,
                               sublevel = T) # Construct filtration using F.values
      
      #extract the persistence diagram
      persistenceDiagram <-
        filtrationDiag(filtration = FltRips, maxdimension = maxDimension)$diagram
      
      #extract betti number - unused
      # b0OrigGraph <- sum(persistenceDiagram[, 1] == 0)
      # b1OrigGraph <- sum(persistenceDiagram[, 1] == 1)
      # b2OrigGraph <- sum(persistenceDiagram[, 1] == 2)
      
      #extract birth and death times
      pd2 <-cbind(persistenceDiagram,GraphId = graphId,dataset = dataAlias)
      results<-rbind(results,pd2)
      
    } else {
      message("Ignoring ",dataAlias," graph ",graphId," Node count:",nodeCount)
    }
  }
  write.table(
    results,
    file = whichOutputFile,
    sep = "\t",
    row.names = FALSE,
    col.names = FALSE,
    append = T,
    quote = FALSE
  )
}

nodeFeatures <- c("degree","betweenness","closeness")#"eccentricity","hub","authority"

nodeFeatures2 <-c("ricci","forman")
for(f in c(nodeFeatures2)){
  # reddit multi pds are missing. 
  if(TRUE){
    extractPD("NCI1/NCI1.","NCI1",feature=f)
    extractPD("ENZYMES/ENZYMES.","Enzyme",feature=f)
    extractPD("BZR/BZR.","BZR",feature=f)
    extractPD("COX2/COX2.","COX2",feature=f)
    extractPD("DHFR/DHFR.","DHFR",feature=f)
    extractPD("FRANKENSTEIN/FRANKENSTEIN.","FRANKENSTEIN",feature=f)
    extractPD("REDDIT-MULTI-5K/REDDIT-MULTI-5K.","REDDIT5K",feature=f)
    extractPD("proteins/proteins.","Protein",feature=f)
    extractPD("REDDIT-BINARY/REDDIT-BINARY.","RedditBinary",feature=f)
    extractPD("IMDB-MULTI/IMDB-MULTI.","IMDBMulti",feature=f)
    extractPD("IMDB-BINARY/IMDB-BINARY.","IMDBBinary",feature=f)
  }
  
}
