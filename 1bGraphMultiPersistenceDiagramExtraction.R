# Objective : Extract sub level "persistence grids" from graphs
# requires dataset folders from http://networkrepository.com/labeled.php
# Created by: **
# Created on: 2021-1-11
library(igraph)
library(TDA)
rm(list = ls())
options(java.parameters = "-Xmx200g")


# 3 Graph ML datasets
dataPath <- "Documents/GraphML/"
outputPath="PeaceCorps/MultiFeature/Grids/"

p1 <-  "ENZYMES/ENZYMES."
p2 <-  "proteins/proteins."
p3 <-  "REDDIT-BINARY/REDDIT-BINARY."
p4 <-  "IMDB-MULTI/IMDB-MULTI."
p5 <-  "NCI1/NCI1."
p6 <-  "IMDB-BINARY/IMDB-BINARY."

# max dimension of the homological features to be computed. (e.g. 0 for connected components, 1 for connected components and loops, 2 for connected components, loops, voids, etc.)
maxDimension <- 1

# upper limit on the simplex dimension, size of the cliques to find
maxClique <-2

# node count thresholds to ignore a graph
minNodeCount <- 4

# how many digits to consider in the grid
featureSensitivity =1

nodeFeatures <- c("closeness","betweenness","hub","degree")#,"eccentricity","authority")

outputFile <- "PersistenceGrid.txt"

# threshold to ignore large graphs
maxNodeCount <- 6000

GRIDCOUNT <- 10

source("nodeFunctions.R")


# function to find the secondary grid id of the value in a row.
get2ndGridId<-function(min,step,birth,death){
  if(birth>death){
    message("birth is later than death?",birth,">",death)
    return(-1)
  }
  gridIdStart<-0
  gridIdEnd<-0
  if(birth<=min){
    gridIdStart<-1
  } else{
    gridIdStart<-1+ceiling((birth-min)/step)
  }
  
  # if the hole has not died, it will have INF as value
  if(!is.finite(death)){
    gridIdEnd<-2+GRIDCOUNT
  }else{
    gridIdEnd<-1+ceiling((death-min)/step)
  }
  if(gridIdStart>gridIdEnd){
    message("An error occurred in the 2nd grid: ",gridIdStart,">",gridIdEnd)
    return(-1)
  }
  return(gridIdStart:gridIdEnd)
}

computePersistentGrid<-function(dataset,dataAlias,feature1,feature2){
  whichOutputFile<-paste0(outputPath,dataAlias,feature1,feature2,outputFile)
  #Check if a previous result file  exists
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
  # I am converting graph vertex ids to strings that starts with V
  # This is done to avoid orig vertex ids where igraph given internal ids
  allData$source<-paste0("v",allData$source)
  allData$target<-paste0("v",allData$target)
  bettiFrame<-data.frame()
  graphs<-unique(graphIdData$V1)
  message("Processing ",dataAlias," for feature ",feature1," and ",feature2,": ",length(graphs)," graphs.")
  
  for(graphId in graphs){
    thisGraphNodes <- paste0("v",which(graphIdData$V1 == graphId))
    edgeData <-
      allData[allData$source %in% thisGraphNodes |
                allData$target %in% thisGraphNodes, ]
    graph <- graph.data.frame(edgeData, directed = FALSE)
    nodeCount <- vcount(graph)
    
    # if the graph is not too small nor too big, compute filtration
    if (nodeCount > minNodeCount && nodeCount < maxNodeCount) {
      # below we use sub level filtrations.
      #sublevel filtration
      # we will compute activation values of nodes based on two features
      
      f1Values = computeNodeVals(feature1, graph,dataset=dataset,dataPath=dataPath)
      f2Values = computeNodeVals(feature2, graph,dataset=dataset,dataPath=dataPath)
      
      # create a 10 step grid for feature1
      grid1Id<-0
      f1Vals <- f1Values[is.finite(f1Values)]
      min<-min(f1Vals)
      max<-max(f1Vals)
      # we cannot do multi persistence if our max and min values are Nan or Inf
      if(!is.finite(min)|!is.finite(max)){
        next
      }
      feature1StepVal = (max-min)/GRIDCOUNT
      val1<-min
      
      # create a grid for feature2
      gridId2<-0
      f2Vals <- f2Values[is.finite(f2Values)]
      min2<-min(f2Vals)
      max2<-max(f2Vals)
      # we cannot do multi persistence if our max and min values are Nan or Inf
      if(!is.finite(min2)|!is.finite(max2)){
        next
      }
      feature2StepVal = (max2-min2)/GRIDCOUNT
      epsilon <- 1e-10
      if(feature1StepVal< epsilon|feature2StepVal< epsilon){
        # we cannot do multi persistence for this graph.
        message("max value is equal to min in ",dataAlias," graph ", graphId,
                ":",feature1," and ",feature2)
        next;
      }
      
      # at this point, not all f1 and f2 values are NaN or Inf
      # we will convert NaN values to min-(step/2) 
      # and Inf values to max+(step/2)
      f1Values[is.na(f1Values)]<-min-(feature1StepVal/2)
      f1Values[!is.finite(f1Values)]<-max+(feature1StepVal/2)
      
      f2Values[is.na(f2Values)]<-min2-feature2StepVal/2
      f2Values[!is.finite(f2Values)]<-max2+feature2StepVal/2
      
      while(grid1Id<=GRIDCOUNT){
        # choose nodes that are active within the grid1
        grid1Id<-grid1Id+1
        vertices<-names(f1Values[f1Values<=val1])
        filteredGraph = induced_subgraph(graph,vertices)
        
        # We will compute multi persistence only if some nodes are activated
        if(vcount(filteredGraph)>0){
          # choose the feature 2 values of active nodes
          filteredF2Values<-f2Values[vertices]
          # for maxClique=3 below means we are adding 0,1,2 simplices (nodes,edges,triangles) to our complex
          cmplx <- cliques(as.undirected(filteredGraph), min = 0, max = maxClique)
          # use sublevel=T for sublevel, sublevel=F for superlevel filtration
          # F.values are node values. At these values, node appear in the complex,
          # and their edges to other active nodes also activate.
          FltRips <- funFiltration(FUNvalues = filteredF2Values,
                                   cmplx = cmplx,
                                   sublevel = T) # Construct filtration using F.values
          
          #extract the persistence diagram
          # if there is a single activated vertex, the code below will give
          # a warning message.
          persistenceDiagram <-
            filtrationDiag(filtration = FltRips, maxdimension = maxDimension)$diagram
           
          bArray<-array(0,dim=c(2,GRIDCOUNT+3))
          
          for(rowIndex in 1:nrow(persistenceDiagram)) {
            row <- persistenceDiagram[rowIndex,]
            # R indexes start from 1. We will add 1 to b0 and b1 
            bettiNum<-1+as.integer(row[["dimension"]])
            birth<-row[["Birth"]]
            death<-row[["Death"]]
            gridId2<-get2ndGridId(min2,feature2StepVal,birth,death)
            if(gridId2!=-1){
              bArray[bettiNum,gridId2]<-(1+bArray[bettiNum,gridId2])
            }
          }
          sumOfHoles <- sum(bArray)
          #  if we had any holes, we will write them to a file
          if(sumOfHoles>0){
            # write the betti number as the last value of the signature 
            bArray[1,GRIDCOUNT+3]<-0
            bArray[2,GRIDCOUNT+3]<-1
            pd2 <-cbind(GraphId = graphId,dataset = dataAlias,
                        Grid1=grid1Id,filVal=val1,f1=feature1,
                        f2=feature2,bArray)
            bettiFrame<-rbind(bettiFrame,pd2)
          }
        }
        
        val1<-val1+feature1StepVal
      }
    }else {
      message("Ignoring ",dataAlias," graph ",graphId," Node count:",nodeCount)
    }
  }
  write.table(
    bettiFrame,
    file = whichOutputFile,
    sep = "\t",
    row.names = FALSE,
    col.names = FALSE,
    append = T,
    quote = FALSE
  )
}




features<-c("closeness","degree","betweenness")#"authority","eccentricity"
features2<-c("ricci","forman")
for(i1 in 1:length(features)){
  for(i2 in (i1+1):length(features)){
    f1 = features[[i1]]
    f2 = features[[i2]]
    
    computePersistentGrid("NCI1/NCI1.","NCI1",feature1=f1,feature2=f2)
    if(FALSE){
    computePersistentGrid("BZR/BZR.","BZR",feature1=f1,feature2=f2)
    computePersistentGrid("REDDIT-MULTI-5K/REDDIT-MULTI-5K.","REDDIT5K",feature1=f1,feature2=f2)
    computePersistentGrid("COX2/COX2.","COX2",feature1=f1,feature2=f2)
    computePersistentGrid("DHFR/DHFR.","DHFR",feature1=f1,feature2=f2)
    computePersistentGrid("FRANKENSTEIN/FRANKENSTEIN.","FRANKENSTEIN",feature1=f1,feature2=f2)
    computePersistentGrid("ENZYMES/ENZYMES.","Enzyme",feature1=f1,feature2=f2)
    computePersistentGrid("proteins/proteins.","Protein",feature1=f1,feature2=f2)
    computePersistentGrid("REDDIT-BINARY/REDDIT-BINARY.","RedditBinary",feature1=f1,feature2=f2)
    computePersistentGrid("IMDB-MULTI/IMDB-MULTI.","IMDBMulti",feature1=f1,feature2=f2)
    
    computePersistentGrid("IMDB-BINARY/IMDB-BINARY.","IMDBBinary",feature1=f1,feature2=f2)
    }
  }
}
 
dataset=  "ENZYMES/ENZYMES."
dataAlias="Enzyme"
feature1="closeness"
feature2="degree"
   






 