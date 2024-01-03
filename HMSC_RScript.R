#This script gives the basic workflow used to run the HMSC model and generate output from this modelling framework.
#As an example, the AFBirds datasets are used here, but the same script should work for all other taxonomic groups (by changing AFBirds to e.g. Bats or Prim)

#Set your working directory
setwd("C:/...")

#Set a seed (for reproducability)
set.seed(1)

#Load required packages
library(Hmsc)
library(ape)
library(plyr)
library(dplyr)
library(Matrix)
library(ggplot2)

#Set directories
localDir = "."
dataDir = file.path(localDir, "data")
modelDir = file.path(localDir, "models")
if(!dir.exists(modelDir)) dir.create(modelDir)
if(!dir.exists(dataDir)) dir.create(dataDir)

#Reading in the data
Y <- read.csv(file.path(dataDir, "Y_AFBirds.csv"), row.names = 1)
envdata <- read.csv(file.path(dataDir, "X_AFBirds.csv"), row.names = 1, stringsAsFactors = TRUE)
envdata <- transform(envdata, siteID = as.integer(match(paste0(Latitude, Longitude), unique(paste0(Latitude, Longitude))))) #Creating unique site IDs
Xdata <- data.frame(ID = rownames(envdata), siteID = envdata$siteID, BSR = envdata$BSR, 
                    ForestCover = envdata$X5kmForestCover, FragmentAge = envdata$PropAge, 
                    MeanTemp = envdata$MeanTemperature, TempSeas = envdata$TemperatureSeasonality, 
                    AnnPrec = envdata$AnnualPrecipitation, PrecSeas = envdata$PrecipitationSeasonality, 
                    Elevation = envdata$Elevation, Slope = envdata$Slope)
rownames(Xdata) <- rownames(envdata)
Xdata$ID <- as.factor(Xdata$ID)
Xdata$siteID <- as.factor(Xdata$siteID)
envdata_distinct <- distinct(envdata, siteID, .keep_all = TRUE)
xy <- as.matrix(cbind(envdata_distinct$Longitude, envdata_distinct$Latitude))
rownames(xy) <- envdata_distinct$siteID
colnames(xy) <- c("x-coordinate","y-coordinate")
TrData <- read.csv(file.path(dataDir, "T_AFBirds.csv"), row.names = 1)
rownames(TrData) <- sub(" ", "_", rownames(TrData))
colnames(Y) <- rownames(TrData)
phyloTree <- read.tree(file.path(dataDir, "C_AFBirds.txt"))
colnames(Y) == rownames(TrData)

hist(TrData$bodymass)
TrData$bodymass <- log(TrData$bodymass) #Bodymass heavily left-skewed, fix using log transformation

#Set up the model
XFormula = ~ BSR + ForestCover + FragmentAge + MeanTemp + TempSeas + AnnPrec + PrecSeas + Elevation + Slope

TrFormula = ~ bodymass + specialisation + forstrat + dispersal

studyDesign <- data.frame(siteID = Xdata$siteID)
rownames(studyDesign) <- rownames(Y)

rL = HmscRandomLevel(sData = xy) #Random effect using the siteIDs, not the original IDs from AFBirds!!!

# Use the Hmsc model constructor to define a model
m = Hmsc(Y = Y, distr = "probit", 
         XData = Xdata,  
         XFormula = XFormula, 
         TrData = TrData, 
         TrFormula = TrFormula, 
         phyloTree = phyloTree,
         studyDesign = studyDesign, 
         ranLevels = list(siteID = rL))
m

#Saving the model
AFBirdmodels = list(m)
names(AFBirdmodels) = c("AFBird")
save(AFBirdmodels, file = file.path(modelDir, "unfitted_AFBirdmodels.RData"))

#Testing model for errors
for(i in 1:length(AFBirdmodels)){
  print(i)
  sampleMcmc(AFBirdmodels[[i]],samples=2)
}

### Fitting the model ###
nParallel = NULL

load(file=file.path(modelDir,"unfitted_AFBirdmodels.RData"))
nm = length(AFBirdmodels)
samples_list = c(5,250,250,250)
thin_list = c(1,1,10,100)
nChains = 4
if(is.null(nParallel)) nParallel = nChains
Lst = 1
while(Lst <= length(samples_list)){
  thin = thin_list[Lst]
  samples = samples_list[Lst]
  print(paste0("thin = ",as.character(thin),"; samples = ",as.character(samples)))
  filename = file.path(modelDir,paste("AFBirdmodels_thin_", as.character(thin),
                                      "_samples_", as.character(samples),
                                      "_chains_",as.character(nChains),
                                      ".Rdata",sep = ""))
  if(file.exists(filename)){
    print("model had been fitted already")
  } else {
    print(date())
    for (mi in 1:nm) {
      print(paste0("AFBirdmodel = ",names(AFBirdmodels)[mi]))
      m = AFBirdmodels[[mi]]
      m = sampleMcmc(m, samples = samples, thin=thin,
                     adaptNf=rep(ceiling(0.4*samples*thin),m$nr), 
                     transient = ceiling(0.5*samples*thin),
                     nChains = nChains,
                     nParallel = nParallel) 
      AFBirdmodels[[mi]] = m
    }
    save(AFBirdmodels,file=filename)
  }
  Lst = Lst + 1
}


### If wanted, MCMC convergence can be calculated after this step. This is not included in this script. ###
### The same holds for cross-validation. See HMSC website for these scripts. ###

#### Predicting species occurrences ###

filename = file.path(modelDir, "AFBirdmodels_thin_100_samples_250_chains_4.Rdata", sep = "")
load(filename)

grid <- read.csv(file.path(dataDir, "SaoPauloPlots_covariates.csv"), stringsAsFactors=TRUE, sep = ";")
levels(grid$ECO_NAME)
grid = droplevels(subset(grid,!(ECO_NAME=="Cerrado"))) #Removing ecoregion absent in the training data
xy.grid <- as.matrix(cbind(grid$xcoord, grid$ycoord))
colnames(xy.grid) <- c("x-coordinate","y-coordinate")
XData.grid = data.frame(ID = as.factor(grid$Plot_ID), BSR = grid$ECO_NAME, ForestCover = grid$X5kmfores_1, 
                        FragmentAge = grid$propage.20 , MeanTemp = grid$MeanTemp1, 
                        TempSeas = grid$TempSeas1, AnnPrec = grid$AnnPrec1, 
                        PrecSeas = grid$PrecSeas1, Elevation = grid$Elevation1, 
                        Slope = grid$Slope1, stringsAsFactors = TRUE)

levels(XData.grid$BSR)[1] <- "Alto Parana Atlantic forest" #Fixing an issue in the name of this ecoregion
levels(XData.grid$BSR)
XData.grid$Slope <- XData.grid$Slope/100 #Correcting slope values

AFBirdmodels
Gradient = prepareGradient(AFBirdmodels[[1]], XDataNew = XData.grid, sDataNew = list(siteID=xy.grid))

nParallel=2
predY = predict(AFBirdmodels[[1]], Gradient=Gradient, expected = TRUE, nParallel=nParallel)
class(predY)
length(predY)
dim(predY[[1]])
head(predY[[1]])
predY

EpredY=Reduce("+",predY)/length(predY)
dim(EpredY)
EpredY

CEpredY <- (EpredY > 0.5)*1
CEpredY

rownames(CEpredY) <- XData.grid$ID
nonzeros <- apply(CEpredY==1, 1, function(x) {
  if(any(x)) {
    c(names(which(x))) 
  }
  else NA
})
nonzeros

write.csv(CEpredY, "PredAFBirds.csv")