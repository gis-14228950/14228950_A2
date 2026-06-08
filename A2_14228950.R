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