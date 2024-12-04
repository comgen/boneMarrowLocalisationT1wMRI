####################### Generated simulated data for non-suture vertices ##############################

# CODE FUNCTIONALITY:
# - This code simulates intensity vertices (from MRI imaging) of vectors perpendicular to the skull surface
# - Currently it simulates 50 layers where each later is 0.5 mm thick.
# - To be more generally usable these variables would need to be used to perform a normalisation / scaling of the data.

minInt = 0; maxInt = 255;# minimum and maximum intensities
nbLayers = 50 # this variable is used to ensure that not more than 50 layers are generated
layerThicknessMm = 0.5 # layer thickness in mm. This variable is currently not used, but would need to be to generalise the code to other layer thicknesses

# Simulation size
nbVerticesToSimul=5000 # most heads have 10000 vertices but only about 5000 of these are on the calvaria
vertexByLayer = matrix(nrow=nbVerticesToSimul, ncol=nbLayers) # stores the intensities for all layers at each vertex
bmFirstLast = matrix(nrow=nbVerticesToSimul, ncol=6) # stores the boundaries of the 3 RoIs: outerBone, boneMarrow, innerBone

# Output directory for the simulation files.
outputDir="/Users/timothyh/Google\ Drive/proj_bma_drive/code_v1.2_simulData50a"

###################################################################
# Intensity and thickness characteristics of anatomical structures
# (Set from inspection of large number of real brain scans)
###################################################################

# Air
airMean=35;airSd=7;
airThickMean=1.5;airThickSd=1; airThickMin=0
# Skin
skinMean=90;skinSd=10;skinMax=150;# use max to define maximum value for skin intensity
skinThickMean=3;skinThickSd=2; skinThickMin=1 # Skin can be anything from 1 to 8 layers thick
# Apo - when not present, mystery or thicker areolar are present
apoMean=45;apoSd=5;
apoThickMean=2;apoThickSd=2; apoThickMin=1 # Scaled by size # apoThickSd=1.3
# Areolar
areolarMean=180; areolarSd=35; areolarMin=100; # min to define a minimum areolar intensity to avoid overlap with skin intensity
areolarThickMean=6;areolarThickSd=2;areolarThickMin=1 # Scaled by adip
# Mystery - only present if aponeurosis is absent
mysteryMean=80; mysterySd=4; # has its own intensity that differs from areolar
mysteryThickMean=4;mysteryThickSd=3;mysteryThickMin=0 # Scaled by size. Absence of aponeurosis layer results in either thicker areolar layer or the mystery layer.
# Outer Bone
outerBoneMean=25; outerBoneSd=8;
outerBoneThickMean=4;outerBoneThickSd=1.5; outerBoneThickMin=2 # scaled by size
# Bone marrow
boneMarrowMean=110; boneMarrowSd=15; # scaled by BMA
boneMarrowThickMean=9;boneMarrowThickSd=3;boneMarrowThickMin=3 # scaled by size
# Inner Bone
#innerBoneMean=25; innerBoneSd=8; No longer needed because using same intensity characteristics for inner and outer bone intensity
innerBoneThickMean=3;innerBoneThickSd=1; innerBoneThickMin=1  # scaled by size # Inner bone shows a lot more variability than outer bone thus the higher sd
# Rest (Anything below the innerbone: arachnoid space and below)
restMean=60; restSd=15;
restThickMin=1 # rest thickness parameters not necessary as determined by thickness of previous layers and max nb layers


# Scaling of intensities:
# - BMA scaled by bmAdip

# Scaling of thicknesses: 
# - most layer thicknesses scale by sizeScale
# - areolar layer thickness scale by adipScale
# - do NOT apply scaling to air (as air presence is an error that should not scale with head size)

# Function that simulates the intensities of all layers, vertex by vertex.
# This function makes use of variables defined in this script (!!NB!!: they are not passed as parameters to the function)
simulate = function(){
  settings=paste("Simulating / ", " size: ",  sizeScale, " - adip: ", adipScale, " - bma: ", boneMarrowMean, sep="")
  print(settings)

  # Initialisation
  # The default is that aponeurosis is present, so we initialise with 0 aponeurosisAbsentVertices (this will randomly get set to >0 in code below).
  aponeurosisAbsentVertices = 0; # want to generate such vertices about 25% of time and them to be about 10 vertices long (this is what was observed in UKB)
  vertex=1;
  while (vertex<=nbVerticesToSimul){ # Iterate over all vertices to be simulated
    # Each of the following lines generates an array of intensities for a given anatomical structure
    # with the length of the array being the number of layers covered by the structure.
    # pmin and pmax are used to prevent the intensity of different layers from overlapping e.g. skin and areolar

    air = rnorm(n=pmax(floor(rnorm(n=1, mean = airThickMean, sd = airThickSd)),airThickMin), mean=rnorm(n=1, mean=airMean,sd=airSd), sd=3) # sometimes there is no air
    skin = pmin(rnorm(n=pmax(floor(rnorm(n=1, mean = skinThickMean, sd = skinThickSd)),skinThickMin), mean=rnorm(n=1, mean=skinMean,sd=skinSd), sd=2),skinMax)

    # Default is to assume that aponeurosis is present and mystery layer is absent
    aponeurosis = rnorm(n=pmax(floor(rnorm(n=1, mean = apoThickMean, sd = apoThickSd)*sizeScale),apoThickMin), mean=rnorm(n=1, mean=apoMean,sd=apoSd), sd=3)
    mystery = c(0); length(mystery)=0
    # Store the mean and standard deviation of areolar intensity, because need it further down for case where mystery has areolar intensities
    areolarMeanSet=rnorm(n=1, mean=areolarMean,sd=areolarSd) # pmax on mean is to ensure that areolar does not enter the range of skin
    areolarIntraVertexSd=pmax(rnorm(n=1, mean=13, sd=5),7) # This is allowing for very high variation in intensity within a vertex which is perhaps not necessary. Note that a lot of the variation in the areolar structure (yellow in core and red on edges) is due to the blurring effect (which is incorporated into the simulation).
    # In default case where mystery is not present, then areolar layer is just standard thickness
    areolar = pmax(rnorm(n=pmax(floor(rnorm(n=1, mean=areolarThickMean, sd=areolarThickSd)*adipScale),apoThickMin), mean=areolarMeanSet, sd=areolarIntraVertexSd),areolarMin) # areolarMin to ensure adip intensity does not overlap with skin intensity
    
    # Determine if we are going to start a new run of aponeurosis
    if(floor(runif(n=1,min=1,max=25))==1){ # 1/40 chance of run beginning
	    aponeurosisAbsentVertices = runif(n=1,min=10,max=20) # length of run
    }
    
    # aponeurosis absent >> then need to sort out intensity of the mystery layer
    if (aponeurosisAbsentVertices>0){
	    length(aponeurosis) = 0; # Obliterate aponeurosis layers
	    # The mystery layer can be up to 8 layers thick and can have the intensity features of either mystery or areolar
	    mysteryThickness = pmax(floor(rnorm(n=1, mean=mysteryThickMean, sd=mysteryThickSd)*sizeScale),mysteryThickMin)
	    if(runif(n=1, min = 0, max = 1) > 0.5){ # Randomly chose 50/50 which type of intensity to generate
		    mystery = pmin(rnorm(n=mysteryThickness, mean=rnorm(n=1, mean=mysteryMean,sd=mysterySd), sd=4),100) # mystery intensities not to exceed 100 in intensity
	    }else{
		    mystery = pmax(rnorm(n=mysteryThickness, mean=areolarMeanSet, sd=areolarIntraVertexSd),areolarMin) # areolar intensities
	    }
	    ## reduce number to generate
	    aponeurosisAbsentVertices = aponeurosisAbsentVertices - 1; 
    }

    # The 3 RoIs (each with minimum thickness requirements)
    obAndIbMeanInt = rnorm(n=1, mean=outerBoneMean,sd=outerBoneSd)
    outerBone = rnorm(n=pmax(floor(rnorm(n=1, mean=outerBoneThickMean, sd=outerBoneThickSd)*sizeScale),outerBoneThickMin), mean=obAndIbMeanInt, sd=4) # in some rare cases the outer bone can be extremely thin
    # !!! IMPORTANT !!! the sd of intensity within a vertex is very important - it is rare that this variance is high
    boneMarrow =  rnorm(n=pmax(floor(rnorm(n=1, mean=boneMarrowThickMean, sd= boneMarrowThickSd)*sizeScale),boneMarrowThickMin), mean=rnorm(n=1, mean=boneMarrowMean,sd=boneMarrowSd), sd=pmax(rnorm(n=1, mean=13, sd=4),1)) 
    innerBone = rnorm(n=pmax(floor(rnorm(n=1, mean=innerBoneThickMean, sd=innerBoneThickSd)*sizeScale),innerBoneThickMin), mean=obAndIbMeanInt, sd=4)

    # Extract the boundaries of the RoIs (the truth used in NN training)
    outerBoneFirstLayer = length(c(air, skin, aponeurosis, areolar, mystery)) +  1
    outerBoneLastLayer = outerBoneFirstLayer + length(outerBone) - 1
    boneMarrowFirstLayer = outerBoneLastLayer +  1 # add 1 as BM starts in next layer
    boneMarrowLastLayer = boneMarrowFirstLayer + length(boneMarrow) - 1
    innerBoneFirstLayer = boneMarrowLastLayer + 1
    innerBoneLastLayer = innerBoneFirstLayer + length(innerBone) - 1

    # aggregate all layers
    layers=c(air, skin, aponeurosis, areolar, mystery, outerBone, boneMarrow, innerBone)

    # Check that we have not generated too many layers >> will not use vertex if generated too many layers
    # >> then generate rest layer
    # >> then check no intensities out of range
    # >> then average layers
    # >> then use vertex
    if(length(layers) + restThickMin <= nbLayers){ # make sure that there are always a few rest layers
    	# need to fill with rest layers
    	# Cap Intensities at 100 (to ensure rest layers with low intensity)
    	rest = pmin(rnorm(nbLayers-length(layers), mean=rnorm(n=1, mean=restMean,sd=restSd), sd=pmax(rnorm(n=1, mean=5, sd=5),1)), 110) # mostly low intra vertex variability but can be higher
    	layers = c(layers, rest) ## add the rest to the layers
    	layers=pmax(pmin(layers, maxInt), minInt) # Make sure that none of the intensities are out of range
    	
    	# Average out the layers by taking a weighted mean of layer and layers above and below
    	averagedLayers=layers; # important to ensure first and last layer have value (don't get one from avg in next line)
    	for (i in 2:(nbLayers-1)){averagedLayers[i] = 0.3*layers[i-1]+0.4*layers[i]+0.3*layers[i+1] }
    	# Updated the two main matrices
    	vertexByLayer[vertex,] <<- averagedLayers;
    	bmFirstLast[vertex,] <<- c(outerBoneFirstLayer, outerBoneLastLayer, boneMarrowFirstLayer,boneMarrowLastLayer,innerBoneFirstLayer, innerBoneLastLayer)
    	vertex=vertex+1 
    }

  } # END of WHILE loop over vertices

} ## end of function


# Setup different different head types: different size, different bone marrow intensities, different body fat thickness

# Important to include very small person with very little fat as they will have the OB at a very shallow depth (and this can be hard to detect)
sizeScales = c(0.6,0.8,1.0,1.2,1.4,1.6) #c(0.8,1.0,1.2) # Largest person is 50% taller than smallest (roughly matches extremes of height distribution). Adding an extra very small head to boost very shallow FB numbers.
boneMarrowMeans = c(40, 50, 70, 90, 110, 140, 170)#c(30,50,70,100,180) # low adiposity BM is difficult to predict, so skew the simulated heads towards low adiposity. Bone marrow must not be allowed to be too intense because gets confused with adipose layer.
adipScales= c(0.3,0.5, 0.7, 1.0, 1.5, 2.0, 3.0, 4.0)# 0.7 and 1.4 # this is body adiposity (not bone marrow) - fatest person is many times fatter than thinnest (see BMI distribution). Have seen examples of scans that have 25 layers of fat

#sizeScales = c(0.8, 0.9, 1.0, 1.1, 1.2) # Largest person is 50% taller than smallest (roughly matches extremes of height distribution)
#boneMarrowMeans = c(30, 40, 50, 100, 180) # low adiposity BM is difficult to predict, so skew the simulated heads towards low adiposity
#adipScales=c(0.7, 1.0, 1.4) # this is body adiposity (not bone marrow) - fatest person is 100% fatter than thinnest (see BMI distribution)


for( sizeScale in sizeScales){
	for( boneMarrowMean in boneMarrowMeans){
		for( adipScale in adipScales){
		  
		simulate(); # simulate datapoints given the scales set in this for loop
		title=paste("Size:", sizeScale, " - BMA:" , boneMarrowMean, " -  Adip: ", adipScale)
		# update to screen which simulation working on
		# print(title) 
		# plotMatrix(vertexByLayer, bmFirstLast, 1, 2000, title)
		
		# Compute normalised data
		vertexByLayerNorm = vertexByLayer / maxInt
		bmFirstLastNorm = bmFirstLast / nbLayers
		
		# Setup output files
		base=paste("size",sizeScale, "-bma", boneMarrowMean, "-adip", adipScale,sep="")
		plotFilename = paste(base, "_plot.png",sep="")
		intensityFilename = paste(base,"_intensities.tab",sep="")
		boundsFilename = paste(base, "_boundsTrue.tab", sep="")
		intensityNormFilename = paste(base,"_intensitiesNorm.tab",sep="")
		boundsNormFilename = paste(base, "_boundsTrueNorm.tab", sep="")		
		
		# Write plot to disk
		#ggsave(file=paste(dir,"/",plotFilename,sep=""))
		# Write data to disk
		# write.table(vertexByLayer, file=paste(outputDir,"/",intensityFilename,sep=""), row.names=FALSE, col.names=FALSE, sep="\t")
		# write.table(bmFirstLast, file=paste(outputDir,"/",boundsFilename,sep=""), row.names=FALSE, col.names=FALSE, sep="\t")
		# Stopped outputitng the non-normalised files.
		write.table(vertexByLayerNorm, file=paste(outputDir,"/",intensityNormFilename,sep=""), row.names=FALSE, col.names=FALSE, sep="\t")
		write.table(bmFirstLastNorm, file=paste(outputDir,"/",boundsNormFilename,sep=""), row.names=FALSE, col.names=FALSE, sep="\t")
		}
	}
}



