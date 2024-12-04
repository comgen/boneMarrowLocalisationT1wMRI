library("ggplot2")
library("pals")

# This script computes key structure stats for each vertex in the input and, for vertices with sufficient total intensity, it generates an intensity plot of all "good" vertices and computes key structure summary stats across these "good" vertices
# Inputs: the intensity matrix and the bounds (computed by the NN)
# Outputs:
# - stats for all vertices
# - stats for calvarial vertices and figure
# - stats across calvarial vertices and figure

# Vertices outside the region of interest have had their intensity array set to 0. We define "calvarial" as vertices that have a sum of intensities above a threshold indicating that they are in the region of the cranium that we wish to analyse.  Note that there can be vertices with high intensities that are outside the region of interest (for example the tops of the ears) and there can be vertices with low intensities that are within the region of interest (relatively unlikely, but could occur).

# At each vertex we sum the intensity of all layers, if this is below the cutoff this indicates that
# the vertex is likely outside the region of interest, and we do not include it in the plots or
# in the summary statistics across "good" vertices.
vertexIntSumCutoff = 1000 
# Expected number of layers in the intensities file
nbLayersPerVertex = 50
# Nb of vertices to put in each panel of the plot (typically 5000 calvarial vertices in a head).
nbVerticesInPanel = 1000
# Maximum intensity of voxel
maxInt = 255
# Mahalanobis cutoff for identifying outliers
MDcutoff = 25

pdf(NULL) # to suppress the creation of the Rplots.pdf

# Variables used for testing
dir="/Users/timothyh/home/proj_bma/010_bmDetectionAlgorithms/test/"
sample="FS_1334172"

# Rscript options.R x y z
# Parse the command line
args = commandArgs(trailingOnly = TRUE)
if(length(args) != 2){
  stop("You must supply the intensity file and the bounds file (in that order)")
}

# Two inputs: intensity file and bounds file
#intensitiesFilePath = paste(dir,sample,"_layers.csv",sep=") # for testing
#boundsFilePath = paste(dir, sample,"_bounds_NN.txt",sep="") # for testing
intensitiesFilePath = args[1];
boundsFilePath = args[2];

cat(paste("Present working directory: ",getwd(), "\n"))
cat(paste("Intensities file: ", intensitiesFilePath, "\n"))
cat(paste("Bounds file: ", boundsFilePath, "\n"))
cat(paste("Assuming separator in both input files is ';'\n"))
cat(paste("Output files will be written to: ", dirname(intensitiesFilePath), "\n"))

# Output files ###########################################################################

# Vertex level - all vertices
outputFilePath_1 = paste(tools::file_path_sans_ext(intensitiesFilePath),"_structData_allVertices.csv",sep="")
# Vertex level - calvarial vertices
outputFilePath_2 = paste(tools::file_path_sans_ext(intensitiesFilePath),"_structData_calvarialVertices.csv",sep="") # also add MD distance
outputFilePath_3 = paste(tools::file_path_sans_ext(intensitiesFilePath),"_structData_calvarialVertices_bmFirstVsBmIntMean.png",sep="") # for checking for outliers (QC)
# Summaries across vertices - scan level data used for QC and for final BMA intensity output
outputFilePath_4a = paste(tools::file_path_sans_ext(intensitiesFilePath),"_structData_calvarialVertices_summaryStats.csv",sep="") # summary stats across vertices
outputFilePath_4b = paste(tools::file_path_sans_ext(intensitiesFilePath),"_structData_calvarialVertices_summaryStats_exclMDoutliers.csv",sep="") # summary stats across vertices
outputFilePath_5 = paste(tools::file_path_sans_ext(intensitiesFilePath),"_structData_calvarialVertices_intensitiesAndBounds.png",sep="") # for checking quality of NN

# Load data ##############################################################################

intensityMatrix = as.matrix(read.table(intensitiesFilePath, header = FALSE, sep = ";"))
bounds = as.matrix(read.table(boundsFilePath, header = FALSE, sep = ";"))

cat("Dimension of intensity matrix: ",dim(intensityMatrix),"\n")
cat("Dimension of bounds matrix: ",dim(bounds),"\n")

# Check size of both matrices: vertices as rows.
# Use layersPerVertex on intensityMatrix to determine if transpose is necessary
if (dim(intensityMatrix)[2] != nbLayersPerVertex){
  if(dim(intensityMatrix)[1] == nbLayersPerVertex){
    cat("Transposing intensity matrix to have vertices as rows\n")
    intensityMatrix = t(intensityMatrix)	
  } else {	
    stop(paste("Neither dimension of the input intensity matrix is equal to ", nbLayersPerVertex, " Make sure first arg is intensities file and second arg is bounds file. Aborting."))
  }
}

if(dim(bounds)[1] != dim(intensityMatrix)[1]){
  if(dim(bounds)[2] == dim(intensityMatrix)[1]){
    cat("Transposing bounds matrix to have rows as vertices\n")	
    bounds = t(bounds)	
  } else {
    stop(paste("Neither dimension of the bounds matrix is equal to the nb of vertices in the input intensity matrix ", dim(intensityMatrix)[1], ". Aborting."))
  }
}
###########################################
# END of loading input data and checking it
###########################################


###########################################################
# Compute the key stats for each vertex and write to file
###########################################################

cat("Computing structure data for each vertex.\n")
# Setup dataframe
vertexMeansAndSd = setNames(data.frame(matrix(ncol = 24, nrow = 0)), nm = c("obIntMean", "bmIntMean", "ibIntMean", "allIntMean", "obIntSum", "bmIntSum", "ibIntSum", "allIntSum", "obIntSd", "bmIntSd", "ibIntSd","allIntSd","obThick", "bmThick", "ibThick", "allThick", "obFirst", "bmFirst", "ibFirst", "allFirst", "obLast", "bmLast", "ibLast", "allLast" ))

for(row in 1:dim(bounds)[1]) {
  ob = intensityMatrix[row,bounds[row,1]:bounds[row,2]]
  bm = intensityMatrix[row,bounds[row,3]:bounds[row,4]]
  ib = intensityMatrix[row,bounds[row,5]:bounds[row,6]]
  all = intensityMatrix[row,1:nbLayersPerVertex]
  
  vertexMeansAndSd[row,]$obFirst = bounds[row,1]
  vertexMeansAndSd[row,]$bmFirst = bounds[row,3]
  vertexMeansAndSd[row,]$ibFirst = bounds[row,5]
  vertexMeansAndSd[row,]$allFirst = 1;
  
  vertexMeansAndSd[row,]$obLast = bounds[row,2]
  vertexMeansAndSd[row,]$bmLast = bounds[row,4]
  vertexMeansAndSd[row,]$ibLast = bounds[row,6]
  vertexMeansAndSd[row,]$allLast = length(all);
  
  vertexMeansAndSd[row,]$obIntMean =  mean(ob)
  vertexMeansAndSd[row,]$bmIntMean = mean(bm)
  vertexMeansAndSd[row,]$ibIntMean = mean(ib)
  vertexMeansAndSd[row,]$allIntMean = mean(all)
  
  vertexMeansAndSd[row,]$obIntSum =  sum(ob)
  vertexMeansAndSd[row,]$bmIntSum = sum(bm)
  vertexMeansAndSd[row,]$ibIntSum = sum(ib)
  vertexMeansAndSd[row,]$allIntSum = sum(all)
  
  vertexMeansAndSd[row,]$obIntSd = sd(ob)
  vertexMeansAndSd[row,]$bmIntSd = sd(bm)
  vertexMeansAndSd[row,]$ibIntSd = sd(ib)
  vertexMeansAndSd[row,]$allIntSd = sd(all)
  
  vertexMeansAndSd[row,]$obThick = length(ob)
  vertexMeansAndSd[row,]$bmThick = length(bm)
  vertexMeansAndSd[row,]$ibThick = length(ib)
  vertexMeansAndSd[row,]$allThick = length(all)

  # SD is NA for scalar (when a layer has a thickness of 1), so need to fix by setting to 0
  vertexMeansAndSd[row,] = ifelse(is.na(vertexMeansAndSd[row,]), 0, vertexMeansAndSd[row,])
  
} # END of loop over vertices

# Write dataframe to file
# Include a header in the output file so that we know what the columns are.
# Use bounds filename as root
cat(">>>>> Outputing to file:", outputFilePath_1, "\n\n")
write.table(vertexMeansAndSd, file=outputFilePath_1, row.names=FALSE, col.names=TRUE, sep=";")

####################################################################################################
# REDUX: reduce number of vertices to those that are in the region of the cranium we wish to analyse.
# For plotting and summary stats across vertices: Only use vertices which have above a certain intensity summed across all layers
####################################################################################################
cat("Reducing to calvarial vertices (vertices with allIntSum greater than cutoff).\n")

useRow = with(vertexMeansAndSd, allIntMean*allThick > vertexIntSumCutoff)

intensityMatrixRedux = intensityMatrix[useRow,]
boundsRedux = bounds[useRow,]
vertexMeansAndSdRedux= vertexMeansAndSd[useRow,]

nbVs = dim(boundsRedux)[1]

####################################################################################################
# 1. On the REDUX dataset, compute the Mahalonobis Distance (MD) based on the 2D distribution (bmFirst vs. bmIntMean) 
# 2. Write this vertex level data to file
# 3. Create a QC plot of the Mahalanobis distribution and write to file
####################################################################################################

cat("Computing Mahalanobis distance on all calvarial vertices .\n")

# Compute Mahalanobis distance
vars = vertexMeansAndSdRedux[,c("bmFirst","bmIntMean")]
m_dist <- mahalanobis(vars, colMeans(vars), cov(vars))
vertexMeansAndSdRedux$MD = round(m_dist,1)

# Write this redux data with MD measure to file
cat(">>>>> Outputing to file:", outputFilePath_2, "\n\n")
write.table(vertexMeansAndSdRedux, file=outputFilePath_2, row.names=FALSE, col.names=TRUE, sep=";")

# Create the QC plot and write to file
cat("Generating point plot for all calvarial vertices of bmFirst vs. bmIntMean (marking MD outliers).\n")
ggplot() + theme_bw() + geom_jitter(data=vertexMeansAndSdRedux, aes(x=bmFirst, y=bmIntMean, colour=MD>25), alpha=0.1, size=0.5) + xlim (0,50) + ylim(0,255) + geom_vline(xintercept=10) + geom_vline(xintercept=30)
cat(">>>>> Outputing to file:", outputFilePath_3, "\n\n")
ggsave(file=outputFilePath_3, height=10, width=10, units="cm")


####################################################################################################
# Compute the summary stats across vertices for all 25 of the vertex variables (4 intensity means, 4 sums, 4 intensity sd, 4 thicknesses, 4 first layers, 4 last layers ,MD). The 4 structures are ob (outerbone), bm (bone marrow), ob (outerbone), all (all layers).
####################################################################################################

computeSummaryStats = function(inputVector){
  result = c(length(inputVector), min(inputVector), median(inputVector), mean(inputVector), max(inputVector), sum(inputVector), sd(inputVector))
  result
}

dfHeader = c("statName", "nbVertices","min", "median", "mean", "max", "sum", "sd")

# First for all calvarial
cat("Computing summary stats across calvarial vertices.\n")
summaryStatsAcrossCalvarialVertices = setNames(data.frame(matrix(ncol = length(dfHeader), nrow = 0)), nm = dfHeader)
for (stat in colnames(vertexMeansAndSdRedux)){
  summaryStatsAcrossCalvarialVertices[nrow(summaryStatsAcrossCalvarialVertices)+1,] = c(stat,computeSummaryStats(vertexMeansAndSdRedux[[stat]]))
}
cat(">>>>> Outputing to file:", outputFilePath_4a, "\n\n")
write.table(summaryStatsAcrossCalvarialVertices, file=outputFilePath_4a, row.names=FALSE, col.names=TRUE, sep=";")

# Second for calvarial but excluding MD outliers
vertexMeansAndSdRedux_ExclOutliers = subset(vertexMeansAndSdRedux, MD < MDcutoff)

cat("Computing summary stats across calvarial vertices **excluding** MD outliers.\n")
summaryStatsAcrossCalvarialVertices_ExclOutliers = setNames(data.frame(matrix(ncol = length(dfHeader), nrow = 0)), nm = dfHeader)
for (stat in colnames(vertexMeansAndSdRedux_ExclOutliers)){
  summaryStatsAcrossCalvarialVertices_ExclOutliers[nrow(summaryStatsAcrossCalvarialVertices_ExclOutliers)+1,] = c(stat,computeSummaryStats(vertexMeansAndSdRedux_ExclOutliers[[stat]]))
}

cat(">>>>> Outputing to file:", outputFilePath_4b, "\n\n")
write.table(summaryStatsAcrossCalvarialVertices_ExclOutliers, file=outputFilePath_4b, row.names=FALSE, col.names=TRUE, sep=";")


###############################################################################################################
## Define function to plot intensities of all vertices and all layers with bound lines
###############################################################################################################

## both inputs are matrices (and not dataframes) so first step of function is to make them into dataframes
# Intensities are in a matrix bcse there are 50 layers
# Bounds are held in both a matrix and a dataframe (vertexMeansAndSdRedux), but use matrix here to match intensities
# Dataframe created for intensities has to have a special form (see below)
plotMatrix = function(verticesByLayers, bounds, firstVertexToPlot, lastVertexToPlot, title){
  
  # Create a matrix with 3 columns (vertexID, layerID, intensity) that we will fill with data and transform to a dataframe
  dfMatrix = matrix(nrow=nrow(verticesByLayers)*ncol(verticesByLayers),ncol=3)
  # Fill the matrix with data
  for(row in 1:nrow(verticesByLayers)) {
    for(col in 1:ncol(verticesByLayers)) {
      dfMatrix[(row-1)*ncol(verticesByLayers) + col,1] = row
      dfMatrix[(row-1)*ncol(verticesByLayers) + col,2] = col
      dfMatrix[(row-1)*ncol(verticesByLayers) + col,3] = verticesByLayers[row,col]
    }
  }
  # Make matrix into dataframe
  df = as.data.frame(dfMatrix)
  names(df) <- c("vertex", "layer", "intensity")
  
  # Create a BM dataframe with 1 vertex ID col + 6 bound columns
  dfBounds = as.data.frame(as.numeric(seq(1,nrow(bounds))))
  names(dfBounds) <- c("vertex")
  dfBounds$obFirst = bounds[,1]
  dfBounds$obLast = bounds[,2]
  dfBounds$bmFirst = bounds[,3]
  dfBounds$bmLast = bounds[,4]
  dfBounds$ibFirst = bounds[,5]
  dfBounds$ibLast = bounds[,6]
  
  # Only retain the vertex range that want to plot
  dfToPlot=subset(df,firstVertexToPlot<= vertex & vertex<=lastVertexToPlot)
  dfBoundstoPlot=subset(dfBounds,firstVertexToPlot<= vertex & vertex<=lastVertexToPlot)
  
  # Create an extra column in each dataframe to determine the panel it will appear in
  dfToPlot$panel = floor(dfToPlot$vertex / nbVerticesInPanel)
  dfBoundstoPlot$panel = floor(dfBoundstoPlot$vertex / nbVerticesInPanel)
  
  lineSize=0.3
  # Plots intensities of layers at every vertex (with heatmap colours) and bounds of BM as lines
  ggplot() + geom_raster(data=dfToPlot, aes(x=vertex-panel*nbVerticesInPanel,y=layer, fill=intensity))  +  theme_bw() + scale_fill_gradientn(colours=c("dark blue","blue","cyan","red", "orange","gold","yellow"),limits=c(0,maxInt),breaks=c(seq(0,250, by=25)))  + scale_y_reverse(expand = c(0, 0)) + scale_x_continuous(expand = c(0, 0)) + geom_hline(yintercept=5, colour="white", size=lineSize) + geom_hline(yintercept=10, colour="white", size=lineSize) + geom_hline(yintercept=15, colour="white", size=lineSize) + geom_hline(yintercept=20, colour="white", size=lineSize) + geom_hline(yintercept=25, colour="white", size=lineSize)  + geom_hline(yintercept=30, colour="white", size=lineSize)  + geom_hline(yintercept=35, colour="white", size=lineSize)  + geom_hline(yintercept=40, colour="white", size=lineSize) + geom_hline(yintercept=45, colour="white", size=lineSize)  + 
    geom_line(data=dfBoundstoPlot, aes(x=vertex-panel*nbVerticesInPanel, y=obFirst), colour="green", size=lineSize) + geom_line(data=dfBoundstoPlot, aes(x=vertex-panel*nbVerticesInPanel, y=obLast), colour="green", size=lineSize) + geom_line(data=dfBoundstoPlot, aes(x=vertex-panel*nbVerticesInPanel, y=bmFirst), colour="white", size=lineSize) + 
    geom_line(data=dfBoundstoPlot, aes(x=vertex-panel*nbVerticesInPanel, y=bmLast), colour="white", size=lineSize) + geom_line(data=dfBoundstoPlot, aes(x=vertex-panel*nbVerticesInPanel, y=ibFirst), colour="green", size=lineSize) + geom_line(data=dfBoundstoPlot, aes(x=vertex-panel*nbVerticesInPanel, y=ibLast), colour="green", size=lineSize) +  facet_grid(rows = vars(panel)) + ggtitle(title)
  
  
} # END of function for plotting intensities


# Plot intensities and bounds (include some of above stats in plot)
cat("Plotting intensities and bounds for calvarial vertices.\n")

obIntMean = as.numeric(subset(summaryStatsAcrossCalvarialVertices, statName=="obIntMean")$mean)
bmIntMean = as.numeric(subset(summaryStatsAcrossCalvarialVertices, statName=="bmIntMean")$mean)
ibIntMean = as.numeric(subset(summaryStatsAcrossCalvarialVertices, statName=="ibIntMean")$mean)

qcMetrics = paste("NbVertices: ", nbVs, " - AvgIntOB: ", round(obIntMean,1), " - AvgIntBM: ", round(bmIntMean,1), " - AvgIntIB: ", round(ibIntMean,1))
title = paste(basename(intensitiesFilePath), " - ", qcMetrics)
# Call the plotting function and save to file
plotMatrix(intensityMatrixRedux, boundsRedux, 1, nrow(intensityMatrixRedux), title)
cat(">>>>> Outputing to file:", outputFilePath_5, "\n\n")
ggsave(file=outputFilePath_5, height=20, width=30, units="cm")

cat("Run complete\n\n")


