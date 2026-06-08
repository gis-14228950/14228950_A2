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
abiotic.vif = vifstep(abiotic, th = 5)
abiotic.vif
abiotic = exclude(abiotic, abiotic.vif)

#screen the habitat group for collinearity
habitat.vif = vifstep(habitat, th = 5)
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

#===================== Forward selection within groups =========================

set.seed(11)

#function to forward select variables within a predictor group using the double stopping rule
selectVars = function(group){
  #full and null models for the group
  full = rda(spe.hel ~ ., data = group)
  null = rda(spe.hel ~ 1, data = group)
  #global test of the group
  print(anova(full, permutations = 999))
  #forward selection bounded by the global adjusted R squared
  sel = ordiR2step(null, scope = formula(full), R2scope = TRUE,
                   direction = "forward", permutations = 999, trace = FALSE)
  #return the group reduced to its selected variables
  group[, attr(terms(formula(sel)), "term.labels"), drop = FALSE]
}

#forward select each multi-variable group
abiotic.sel = selectVars(abiotic)
habitat.sel = selectVars(habitat)
space.sel   = selectVars(space)

#inspect the retained variables in each group
names(abiotic.sel)
names(habitat.sel)
names(space.sel)

#========================= Test the unique fractions ===========================

set.seed(11)

#pure abiotic fraction, conditioned on the other three groups
anova(rda(spe.hel, abiotic.sel, cbind(manage, habitat.sel, space.sel)), permutations = 999)

#pure management fraction
anova(rda(spe.hel, manage, cbind(abiotic.sel, habitat.sel, space.sel)), permutations = 999)

#pure habitat fraction
anova(rda(spe.hel, habitat.sel, cbind(abiotic.sel, manage, space.sel)), permutations = 999)

#pure spatial fraction
anova(rda(spe.hel, space.sel, cbind(abiotic.sel, manage, habitat.sel)), permutations = 999)