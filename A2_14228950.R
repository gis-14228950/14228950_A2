#===============================================================================
####Script to disentangle environmental, management and spatial drivers of beetle community composition across Scotland
# Student ID: 14228950
#===============================================================================

#Set working directory
setwd("W:/1UOM/71922 Spactial Ecology/assessment2/Assessment2_Data_GEOG71922/Beetles")

#community ordination and variation partitioning
library(vegan)

#spatial eigenvector maps and forward selection
library(adespatial)

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
