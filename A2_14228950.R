#===============================================================================
####Script to disentangle environmental, management and spatial drivers of beetle community composition across Scotland
# Student ID: 14228950
#===============================================================================

#Set working directory
setwd("W:/1UOM/71922 Spactial Ecology/assessment2/Assessment2_Data_GEOG71922/Beetles")

#community ordination and variation partitioning
library(vegan)

#variance inflation factors
library(usdm)

#===================== Data import and integrity checks ========================

#read in the community matrix
comm = read.csv("scot_beetle_community.csv", row.names = 1)

#read in the environment table
env = read.csv("scot_beetle_env.csv", row.names = 1)

#inspect
head(comm)
head(env)

#confirm dimensions
dim(comm)
dim(env)

#check site labels 
all(comm$Sites == env$Sites)

#check for missing values
sum(is.na(comm))
sum(is.na(env))

#check for duplicated sites
any(duplicated(env$Sites))

#drop the Sites column
spe = comm[, -1]

#check for negative abundances
min(spe)

#confirm the species matrix
dim(spe)

#======================== Community data processing ============================

#count site occurrences per species
colSums(spe > 0)

#Hellinger transform the abundances
spe.hel = decostand(spe, method = "hellinger")

#inspect
head(spe.hel)

#====================== Environmental data processing ==========================

#inspect distributions of the candidate predictors
summary(env)

#log transform strictly positive variables
logVars = c("Org", "AvailP", "AvailK", "Litter", "Plants_m2", "Elevation")
env[logVars] = log(env[logVars])

#square root transform the skewed cover variable containing zeros
env$Bryophyte = sqrt(env$Bryophyte)

#standardise all continuous predictors to zero mean and unit variance
env.z = decostand(env[, 2:15], method = "standardize")

#build the four predictor groups
abiotic = env.z[, c("Elevation", "pH", "Moist", "Org")]
manage  = env.z[, "Management", drop = FALSE]
habitat = env.z[, c("AvailP", "AvailK", "Litter", "Bryophyte", "Plants_m2",
                    "CanopyHeight", "Stem.density", "Biom0_5", "Repro.biom")]

#keep the coordinates (untransformed) for the spatial section
coords = env[, c("X", "Y")]

#========================= Collinearity screening ==============================

#inspect correlations within the habitat group
round(cor(habitat), 2)

#screen the abiotic group for collinearity
abiotic.vif = vifstep(abiotic, th = 10)
abiotic.vif
abiotic = exclude(abiotic, abiotic.vif)

#screen the habitat group for collinearity
habitat.vif = vifstep(habitat, th = 10)
habitat.vif
habitat = exclude(habitat, habitat.vif)

#============================= Spatial variable ================================

#build distance matrix
dist.xy = dist(coords)

#approximate spatial structure
pcnm.xy = pcnm(dist.xy)

#extract the PCNM axes as candidate spatial predictors
space = data.frame(scores(pcnm.xy))

#generat number of spatial predictors
ncol(space)