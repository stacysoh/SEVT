setwd("C:/Users/user/Desktop/EVT")
rm(list=ls())

#Load packages
library(dplyr)
library(tidyr)
library(SpatialExtremes)
library(extRemes)

#Load output from 01_SEVT_Max_Stable_Models first

##################################################################################
### To obtain the per unit change in the spatial covariate of interest:        ###
### First, impose a unit change in the spatial covariate of interest, then use ###
### the locCoeff for all spatial covariates and matrix multiply with the new   ### 
### spatial covariate of interest(+0.01) to obtain the new parameters          ###
### (location, scale, shape) for each location (for each X and Y).             ###
### Finally, inverse the qgev distribution to get the new return levels.       ###
##################################################################################

### Best model with lowest deviance: ms4.1t_brown ###

########################################################
#### Spatial covariates with positive coefficients #####
########################################################
## Freshwater surfaces
## impervious Surfaces
## Vegetation with human management (with tree canopy)
## Population size
## Median building age 

###################################
###### FRESHWATER SURFACES ########
###################################
loc.form.4.1 <- y ~ xcoord + ls_v1
scale.form.4.1 <- y ~  xcoord + ycoord + ls_v1
shape.form.4.1 <- y ~ 1

data <- ms4.1t_brown$coord 

loc.dsgnmat.ms4.1t_brown <- modeldef(data, loc.form.4.1)$dsgn.mat 
scale.dsgnmat.ms4.1t_brown <- modeldef(data, scale.form.4.1)$dsgn.mat 
shape.dsgnmat.ms4.1t_brown <- modeldef(data, shape.form.4.1)$dsgn.mat 


####to get loc.pred #### 
unit_change = 0.01 #impose a change of 1%
unit_change = 0.05 #impose a change of 5%
unit_change = 0.1 #impose a change of 10%
View(loc.dsgnmat.ms4.1t_brown)
loc.design.ms4.1t_brown = cbind(loc.dsgnmat.ms4.1t_brown[,c(1:2)], 
                                loc.dsgnmat.ms4.1t_brown[, 3] + unit_change, #add 0.01 to col 3 freshwater
                                loc.dsgnmat.ms4.1t_brown[,c(4:11)]) 
colnames(loc.design.ms4.1t_brown)[3] = "freshwater"
View(loc.design.ms4.1t_brown)

loc_coef.ms4.1t_brown = t(ms4.1t_brown$param) #11 coefficients
View(loc_coef.ms4.1t_brown)
loc_coef.ms4.1t_brown = loc_coef.ms4.1t_brown[,c(3:13)] #extract columns 4 to 12 
loc_coef_pred.ms4.1t_brown = as.matrix(loc.design.ms4.1t_brown) %*% loc_coef.ms4.1t_brown #matrix multiply (1618x1)

####to get scale.pred #### 
View(scale.dsgnmat.ms4.1t_brown)
scale.design.ms4.1t_brown = cbind(scale.dsgnmat.ms4.1t_brown[,c(1:3)], 
                     scale.dsgnmat.ms4.1t_brown[, 4] + unit_change, #add 0.01 to col 4 freshwater
                     scale.dsgnmat.ms4.1t_brown[,c(5:12)]) 
colnames(scale.design.ms4.1t_brown)[4] = "freshwater"
View(scale.design.ms4.1t_brown)
scale_coef.ms4.1t_brown = t(ms4.1t_brown$param) #11  coefficients
View(scale_coef.ms4.1t_brown)
scale_coef.ms4.1t_brown = scale_coef.ms4.1t_brown[,c(14:25)]
scale_coef_pred.ms4.1t_brown = as.matrix(scale.design.ms4.1t_brown) %*% scale_coef.ms4.1t_brown

####to get shape.pred #### constant
shape.value.ms4.1t_brown = ms4.1t_brown$param["shapeCoeff1"] 

#Invert qgev distribution to produce new 30-year return levels, requires location, scale and shape
ret.per=30
for (T in ret.per) ret.lev.ms4.1t_brown <-  qgev(1 - 1/T, loc_coef_pred.ms4.1t_brown, 
                                    scale_coef_pred.ms4.1t_brown, shape.value.ms4.1t_brown)
predicted_ms4.1t_brown_1 <- as.data.frame(predicted_ms4.1t_brown)

fw_change = ret.lev.ms4.1t_brown - predicted_ms4.1t_brown_1$Q30 
fw_percent_change = ((ret.lev.ms4.1t_brown - predicted_ms4.1t_brown_1$Q30)/predicted_ms4.1t_brown_1$Q30)*100
summary(fw_change) #0.3 increase in return levels for a 10% increase in freshwater surfaces
summary(fw_percent_change) #mean of 0.7% for a 10% increase in freshwater surfaces


### Plot per unit change effect of freshwater on ret levels ###
library(rgeos)
library(rgdal)
library(tmap)
library(grid)

#Load the hexagons shapefile into R
spat <- readOGR(dsn="C:/Users/user/Desktop/EVT/hexagonsData", layer="hexagonsData")

#Change in return levels (per unit change in spatial cov): freshwater of BEST MODEL 
ret.lev.ms4.1t_brown = as.data.frame(ret.lev.ms4.1t_brown)
predicted_ms4.1t_brown_df_rl_fw <- cbind(df_xy, ret.lev.ms4.1t_brown)
predicted_ms4.1t_brown_df_rl_fw["V1"][predicted_ms4.1t_brown_df_rl_fw["V1"]<0] <- 0 #force 0 for neg values
predicted_ms4.1t_brown_df_rl_fw$V1<- round(predicted_ms4.1t_brown_df_rl_fw$V1,
                                           0)
colnames(predicted_ms4.1t_brown_df_rl_fw)[4] = "Return Levels"
ms4.1t_brown_rl_fw = merge(spat, predicted_ms4.1t_brown_df_rl_fw, by.x="cell_id", by.y="cell_id")

#1 unit change
ms4.1t_brown_retlev_fw_1 = tm_shape(ms4.1t_brown_rl_fw) +
  tm_fill("Return Levels", palette = "Blues", style="fixed",title="", colorNA=grey(0.9),
          breaks = c(13, 35, 45, 55, 66, 92), 
          textNA = "None recorded")+
  tm_layout(legend.text.size = 1, legend.title.size = 1.5,
            title="A",title.size = 2, frame=F, legend.show=F) 
ms4.1t_brown_retlev_fw_1


#5 unit change
ms4.1t_brown_retlev_fw_5 = tm_shape(ms4.1t_brown_rl_fw) +
  tm_fill("Return Levels", palette = "Blues", style="fixed",title="", colorNA=grey(0.9),
          breaks = c(13, 35, 45, 55, 66, 92), 
          textNA = "None recorded")+
  tm_layout(legend.text.size = 1, legend.title.size = 1.5,
            title="B",title.size = 2, frame=F, legend.show=F) 
ms4.1t_brown_retlev_fw_5


#10 unit change
ms4.1t_brown_retlev_fw_10 = tm_shape(ms4.1t_brown_rl_fw) +
  tm_fill("Return Levels", palette = "Blues", style="fixed",title="", colorNA=grey(0.9),
          breaks = c(13, 35, 45, 55, 66, 92), #0, 12, 30, 55, 64, 92
          textNA = "None recorded")+
  tm_layout(legend.text.size = 1, legend.title.size = 1.5,
            title="C",title.size = 2, frame=F, legend.show=T) 
ms4.1t_brown_retlev_fw_10

###################################
###### IMPERVIOUS SURFACES ########
###################################

####to get loc.pred #### 
unit_change = 0.01
unit_change = 0.05
unit_change = 0.1
View(loc.dsgnmat.ms4.1t_brown)
loc.design.ms4.1t_brown = cbind(loc.dsgnmat.ms4.1t_brown[,c(1:3)], 
                                loc.dsgnmat.ms4.1t_brown[, 4] + unit_change, #add 0.01 to col 4
                                loc.dsgnmat.ms4.1t_brown[,c(5:11)]) 
colnames(loc.design.ms4.1t_brown)[4] = "impervious"
View(loc.design.ms4.1t_brown)

loc_coef.ms4.1t_brown = t(ms4.1t_brown$param) #9 coefficients
View(loc_coef.ms4.1t_brown)
loc_coef.ms4.1t_brown = loc_coef.ms4.1t_brown[,c(3:13)] #extract columns 3 to 11 
loc_coef_pred.ms4.1t_brown = as.matrix(loc.design.ms4.1t_brown) %*% loc_coef.ms4.1t_brown #matrix multiply (1618x1)

####to get scale.pred #### 
View(scale.dsgnmat.ms4.1t_brown)
scale.design.ms4.1t_brown = cbind(scale.dsgnmat.ms4.1t_brown[,c(1:4)], 
                                  scale.dsgnmat.ms4.1t_brown[, 5] + unit_change, #add 0.01 to col 4 freshwater
                                  scale.dsgnmat.ms4.1t_brown[,c(6:12)]) 
colnames(scale.design.ms4.1t_brown)[5] = "impervious"
View(scale.design.ms4.1t_brown)
scale_coef.ms4.1t_brown = t(ms4.1t_brown$param) #11  coefficients
View(scale_coef.ms4.1t_brown)
scale_coef.ms4.1t_brown = scale_coef.ms4.1t_brown[,c(14:25)]
scale_coef_pred.ms4.1t_brown = as.matrix(scale.design.ms4.1t_brown) %*% scale_coef.ms4.1t_brown

####to get shape.pred #### constant
shape.value.ms4.1t_brown = ms4.1t_brown$param["shapeCoeff1"] 

#Invert qgev distribution to produce new 30-year return levels, requires location, scale and shape
ret.per=30
for (T in ret.per) ret.lev.ms4.1t_brown <-  qgev(1 - 1/T, loc_coef_pred.ms4.1t_brown, 
                                                 scale_coef_pred.ms4.1t_brown, shape.value.ms4.1t_brown)
predicted_ms4.1t_brown_1 <- as.data.frame(predicted_ms4.1t_brown)

imp_change = ret.lev.ms4.1t_brown - predicted_ms4.1t_brown_1$Q30 
imp_percent_change = ((ret.lev.ms4.1t_brown - predicted_ms4.1t_brown_1$Q30)/predicted_ms4.1t_brown_1$Q30)*100

summary(imp_change) #1.6 increase in return levels for a 10% increase in impervious surfaces
summary(imp_percent_change) #3.3% mean change for a 10% increase

### Plot per unit change effect of impervious surfaces on return levels ###
#Change in return levels (per unit change in spatial cov: impervious surfaces of BEST MODEL 
ret.lev.ms4.1t_brown = as.data.frame(ret.lev.ms4.1t_brown)
predicted_ms4.1t_brown_df_rl_imp <- cbind(df_xy, ret.lev.ms4.1t_brown)
#View(predicted_ms4.1t_brown_df_rl_imp)
predicted_ms4.1t_brown_df_rl_imp["V1"][predicted_ms4.1t_brown_df_rl_imp["V1"]<0] <- 0 #force 0 for neg values
predicted_ms4.1t_brown_df_rl_imp$V1<- round(predicted_ms4.1t_brown_df_rl_imp$V1,
                                           0)
colnames(predicted_ms4.1t_brown_df_rl_imp)[4] = "Return Levels"
ms4.1t_brown_rl_imp = merge(spat, predicted_ms4.1t_brown_df_rl_imp, by.x="cell_id", by.y="cell_id")


#1 unit change
ms4.1t_brown_retlev_imp_1 = tm_shape(ms4.1t_brown_rl_imp) +
  tm_fill("Return Levels", palette = "BuPu", style="fixed",title="", colorNA=grey(0.9),
          breaks = c(13, 35, 45, 55, 66, 94), 
          textNA = "None recorded")+
  tm_layout(legend.text.size = 1, legend.title.size = 1.5,
            title="D",title.size = 2, frame=F, legend.show=F) 
ms4.1t_brown_retlev_imp_1



#5 unit change
ms4.1t_brown_retlev_imp_5 = tm_shape(ms4.1t_brown_rl_imp) +
  tm_fill("Return Levels", palette = "BuPu", style="fixed",title="", colorNA=grey(0.9),
          breaks = c(13, 35, 45, 55, 66, 94), 
          textNA = "None recorded")+
  tm_layout(legend.text.size = 1, legend.title.size = 1.5,
            title="E",title.size = 2, frame=F, legend.show=F) 
ms4.1t_brown_retlev_imp_5


#10 unit change
ms4.1t_brown_retlev_imp_10 = tm_shape(ms4.1t_brown_rl_imp) +
  tm_fill("Return Levels", palette = "BuPu", style="fixed",title="", colorNA=grey(0.9),
          breaks = c(13, 35, 45, 55, 66, 94), 
          textNA = "None recorded"
          )+
  tm_layout(legend.text.size = 1, legend.title.size = 1.5,
            title="F",title.size = 2, frame=F, legend.show=T) 
ms4.1t_brown_retlev_imp_10


##########################################################
### VEGETATION WITH HUMAN MANGEMENT (WITH TREE CANOPY) ###
##########################################################

####to get loc.pred #### 
unit_change = 0.01
unit_change = 0.05
unit_change = 0.1
View(loc.dsgnmat.ms4.1t_brown)
loc.design.ms4.1t_brown = cbind(loc.dsgnmat.ms4.1t_brown[,c(1:5)], 
                                loc.dsgnmat.ms4.1t_brown[, 6] + unit_change,
                                loc.dsgnmat.ms4.1t_brown[,c(7:11)])

colnames(loc.design.ms4.1t_brown)[6] = "veghmtc"
View(loc.design.ms4.1t_brown)

loc_coef.ms4.1t_brown = t(ms4.1t_brown$param) #9 coefficients
View(loc_coef.ms4.1t_brown)
loc_coef.ms4.1t_brown = loc_coef.ms4.1t_brown[,c(3:13)] #extract columns 3 to 11 
loc_coef_pred.ms4.1t_brown = as.matrix(loc.design.ms4.1t_brown) %*% loc_coef.ms4.1t_brown #matrix multiply (1618x1)

####to get scale.pred #### 
View(scale.dsgnmat.ms4.1t_brown)
scale.design.ms4.1t_brown = cbind(scale.dsgnmat.ms4.1t_brown[,c(1:6)], 
                                  scale.dsgnmat.ms4.1t_brown[, 7] 
                                  + unit_change, scale.dsgnmat.ms4.1t_brown[,c(8:12)] ) 
colnames(scale.design.ms4.1t_brown)[7] = "veghmtc"
View(scale.design.ms4.1t_brown)
scale_coef.ms4.1t_brown = t(ms4.1t_brown$param) #11  coefficients
View(scale_coef.ms4.1t_brown)
scale_coef.ms4.1t_brown = scale_coef.ms4.1t_brown[,c(14:25)]
scale_coef_pred.ms4.1t_brown = as.matrix(scale.design.ms4.1t_brown) %*% scale_coef.ms4.1t_brown

####to get shape.pred #### constant
shape.value.ms4.1t_brown = ms4.1t_brown$param["shapeCoeff1"] 

#Invert qgev distribution to produce new 30-year return levels, requires location, scale and shape
ret.per=30
for (T in ret.per) ret.lev.ms4.1t_brown <-  qgev(1 - 1/T, 
                                                 loc_coef_pred.ms4.1t_brown, 
                                                 scale_coef_pred.ms4.1t_brown, 
                                                 shape.value.ms4.1t_brown)
predicted_ms4.1t_brown_1 <- as.data.frame(predicted_ms4.1t_brown)

veghmtc_change = ret.lev.ms4.1t_brown - predicted_ms4.1t_brown_1$Q30 
veghmtc_percent_change = ((ret.lev.ms4.1t_brown - predicted_ms4.1t_brown_1$Q30)/predicted_ms4.1t_brown_1$Q30)*100

summary(veghmtc_change) #0.7 increase in return levels for a 10% increase in veghmtc surfaces
summary(veghmtc_percent_change) #mean of 1.4% for a 10% increase in veghmtc surfaces


### Plot per unit change effect of veghmtc on ret levels ###

#Change in return levels (per unit change in spatial cov: veghmtc of BEST MODEL 
ret.lev.ms4.1t_brown = as.data.frame(ret.lev.ms4.1t_brown)
predicted_ms4.1t_brown_df_rl_veghmtc <- cbind(df_xy, ret.lev.ms4.1t_brown)
predicted_ms4.1t_brown_df_rl_veghmtc["V1"][predicted_ms4.1t_brown_df_rl_veghmtc["V1"]<0] <- 0 #force 0 for neg values
predicted_ms4.1t_brown_df_rl_veghmtc$V1<- round(predicted_ms4.1t_brown_df_rl_veghmtc$V1,
                                                   0)
colnames(predicted_ms4.1t_brown_df_rl_veghmtc)[4] = "Return Levels"
ms4.1t_brown_rl_veghmtc = merge(spat, predicted_ms4.1t_brown_df_rl_veghmtc, by.x="cell_id", by.y="cell_id")


#1 unit change
ms4.1t_brown_retlev_veghmtc_1 = tm_shape(ms4.1t_brown_rl_veghmtc) +
  tm_fill("Return Levels", palette = "BuGn", style="fixed",title="", colorNA=grey(0.9),
          breaks = c(13, 35, 45, 55, 65, 93), 
          textNA = "None recorded")+
  tm_layout(legend.text.size = 1, legend.title.size = 1.5,
            title="G",title.size = 2, frame=F, legend.show=F) 
ms4.1t_brown_retlev_veghmtc_1


#5 unit change
ms4.1t_brown_retlev_veghmtc_5 = tm_shape(ms4.1t_brown_rl_veghmtc) +
  tm_fill("Return Levels", palette = "BuGn", style="fixed",title="", colorNA=grey(0.9),
          breaks = c(13, 35, 45, 55, 65, 93), 
          textNA = "None recorded")+
  tm_layout(legend.text.size = 1, legend.title.size = 1.5,
            title="H",title.size = 2, frame=F, legend.show=F) 
ms4.1t_brown_retlev_veghmtc_5


#10 unit change
ms4.1t_brown_retlev_veghmtc_10 = tm_shape(ms4.1t_brown_rl_veghmtc) +
  tm_fill("Return Levels", palette = "BuGn", style="fixed",title="", colorNA=grey(0.9),
          breaks = c(13, 35, 45, 55, 65, 93), 
          textNA = "None recorded"
  )+
  tm_layout(legend.text.size = 1, legend.title.size = 1.5,
            title="I",title.size = 2, frame=F, legend.show=T) 
ms4.1t_brown_retlev_veghmtc_10


########################
### TOTAL POPULATION ###
########################

####to get loc.pred #### 
unit_change = 0.01
unit_change = 0.05
unit_change = 0.1
View(loc.dsgnmat.ms4.1t_brown)
loc.design.ms4.1t_brown = cbind(loc.dsgnmat.ms4.1t_brown[,c(1:9)], 
                                loc.dsgnmat.ms4.1t_brown[, 10]*(1+unit_change),
                                loc.dsgnmat.ms4.1t_brown[,c(11)])

colnames(loc.design.ms4.1t_brown)[10] = "totalpop"
View(loc.design.ms4.1t_brown)

loc_coef.ms4.1t_brown = t(ms4.1t_brown$param) #9 coefficients
View(loc_coef.ms4.1t_brown)
loc_coef.ms4.1t_brown = loc_coef.ms4.1t_brown[,c(3:13)] #extract columns 3 to 11 
loc_coef_pred.ms4.1t_brown = as.matrix(loc.design.ms4.1t_brown) %*% loc_coef.ms4.1t_brown #matrix multiply (1618x1)

####to get scale.pred #### 
View(scale.dsgnmat.ms4.1t_brown)
scale.design.ms4.1t_brown = cbind(scale.dsgnmat.ms4.1t_brown[,c(1:10)], 
                                  scale.dsgnmat.ms4.1t_brown[, 11]*(1+unit_change), 
                                  scale.dsgnmat.ms4.1t_brown[,c(12)] ) 
colnames(scale.design.ms4.1t_brown)[11] = "totalpop"
View(scale.design.ms4.1t_brown)
scale_coef.ms4.1t_brown = t(ms4.1t_brown$param) #11  coefficients
View(scale_coef.ms4.1t_brown)
scale_coef.ms4.1t_brown = scale_coef.ms4.1t_brown[,c(14:25)]
scale_coef_pred.ms4.1t_brown = as.matrix(scale.design.ms4.1t_brown) %*% scale_coef.ms4.1t_brown

####to get shape.pred #### constant
shape.value.ms4.1t_brown = ms4.1t_brown$param["shapeCoeff1"] 

#Invert qgev distribution to produce new 30-year return levels, requires location, scale and shape
ret.per=30
for (T in ret.per) ret.lev.ms4.1t_brown <-  qgev(1 - 1/T, 
                                                 loc_coef_pred.ms4.1t_brown, 
                                                 scale_coef_pred.ms4.1t_brown, 
                                                 shape.value.ms4.1t_brown)
predicted_ms4.1t_brown_1 <- as.data.frame(predicted_ms4.1t_brown)

pop_change = ret.lev.ms4.1t_brown - predicted_ms4.1t_brown_1$Q30 
pop_percent_change = ((ret.lev.ms4.1t_brown - predicted_ms4.1t_brown_1$Q30)/predicted_ms4.1t_brown_1$Q30)*100

summary(pop_change) #0.1 increase in return levels for a 10% increase in population
summary(pop_percent_change) #abs value mean of 0.3% for a 10% increase in population

### Plot per unit change effect of population on ret levels ###

#Change in return levels (per unit change in spatial cov: population of BEST MODEL 
ret.lev.ms4.1t_brown = as.data.frame(ret.lev.ms4.1t_brown)
predicted_ms4.1t_brown_df_rl_totalpop <- cbind(df_xy, ret.lev.ms4.1t_brown)
predicted_ms4.1t_brown_df_rl_totalpop["V1"][predicted_ms4.1t_brown_df_rl_totalpop["V1"]<0] <- 0 #force 0 for neg values
predicted_ms4.1t_brown_df_rl_totalpop$V1<- round(predicted_ms4.1t_brown_df_rl_totalpop$V1,
                                                0)
colnames(predicted_ms4.1t_brown_df_rl_totalpop)[4] = "Return Levels"
ms4.1t_brown_rl_totalpop = merge(spat, predicted_ms4.1t_brown_df_rl_totalpop, by.x="cell_id", by.y="cell_id")


#1 unit change
ms4.1t_brown_retlev_totalpop_1 = tm_shape(ms4.1t_brown_rl_totalpop) +
  tm_fill("Return Levels", palette = "RdPu", style="fixed",title="", colorNA=grey(0.9),
          breaks = c(13, 35, 45, 55, 65, 92), 
          textNA = "None recorded")+
  tm_layout(legend.text.size = 1, legend.title.size = 1.5,
            title="J",title.size = 2, frame=F, legend.show=F) 
ms4.1t_brown_retlev_totalpop_1


#5 unit change
ms4.1t_brown_retlev_totalpop_5 = tm_shape(ms4.1t_brown_rl_totalpop) +
  tm_fill("Return Levels", palette = "RdPu", style="fixed",title="", colorNA=grey(0.9),
          breaks = c(13, 35, 45, 55, 65, 92), 
          textNA = "None recorded")+
  tm_layout(legend.text.size = 1, legend.title.size = 1.5,
            title="K",title.size = 2, frame=F, legend.show=F) 
ms4.1t_brown_retlev_totalpop_5


#10 unit change
ms4.1t_brown_retlev_totalpop_10 = tm_shape(ms4.1t_brown_rl_totalpop) +
  tm_fill("Return Levels", palette = "RdPu", style="fixed",title="", colorNA=grey(0.9),
          breaks = c(13, 35, 45, 55, 65, 92), 
          textNA = "None recorded"
  )+
  tm_layout(legend.text.size = 1, legend.title.size = 1.5,
            title="L",title.size = 2, frame=F, legend.show=T) 
ms4.1t_brown_retlev_totalpop_10


#######################################
### MEDIAN AGE OF PUBLIC APARTMENTS ###
#######################################

####to get loc.pred #### 
unit_change = 1 #Impose an increase of 1 year
unit_change = 2 #Impose an increase of 2 years
unit_change = 3 #Impose an increase of 3 years
View(loc.dsgnmat.ms4.1t_brown)
loc.design.ms4.1t_brown = cbind(loc.dsgnmat.ms4.1t_brown[,c(1:10)], 
                                loc.dsgnmat.ms4.1t_brown[, 11] + unit_change)

colnames(loc.design.ms4.1t_brown)[11] = "hdbmedage"
View(loc.design.ms4.1t_brown)

loc_coef.ms4.1t_brown = t(ms4.1t_brown$param) #9 coefficients
View(loc_coef.ms4.1t_brown)
loc_coef.ms4.1t_brown = loc_coef.ms4.1t_brown[,c(3:13)] #extract columns 3 to 11 
loc_coef_pred.ms4.1t_brown = as.matrix(loc.design.ms4.1t_brown) %*% loc_coef.ms4.1t_brown #matrix multiply (1618x1)

####to get scale.pred #### 
View(scale.dsgnmat.ms4.1t_brown)
scale.design.ms4.1t_brown = cbind(scale.dsgnmat.ms4.1t_brown[,c(1:11)], 
                                  scale.dsgnmat.ms4.1t_brown[, 12] 
                                  + unit_change) 
colnames(scale.design.ms4.1t_brown)[12] = "hdbmedage"
View(scale.design.ms4.1t_brown)
scale_coef.ms4.1t_brown = t(ms4.1t_brown$param) #11  coefficients
View(scale_coef.ms4.1t_brown)
scale_coef.ms4.1t_brown = scale_coef.ms4.1t_brown[,c(14:25)]
scale_coef_pred.ms4.1t_brown = as.matrix(scale.design.ms4.1t_brown) %*% scale_coef.ms4.1t_brown

####to get shape.pred #### constant,no need
shape.value.ms4.1t_brown = ms4.1t_brown$param["shapeCoeff1"] 

#Invert qgev distribution to produce new 30-year return levels, requires location, scale and shape
ret.per=30
for (T in ret.per) ret.lev.ms4.1t_brown <-  qgev(1 - 1/T, 
                                                 loc_coef_pred.ms4.1t_brown, 
                                                 scale_coef_pred.ms4.1t_brown, 
                                                 shape.value.ms4.1t_brown)
predicted_ms4.1t_brown_1 <- as.data.frame(predicted_ms4.1t_brown)

hdbage_change = ret.lev.ms4.1t_brown - predicted_ms4.1t_brown_1$Q30 
hdbage_percent_change = ((ret.lev.ms4.1t_brown - predicted_ms4.1t_brown_1$Q30)/predicted_ms4.1t_brown_1$Q30)*100

summary(hdbage_change) #1.8 increase in return levels for a 3 year increase in building age
summary(hdbage_percent_change) #3.8% mean for a 3 year increase in building age

### Plot per unit change effect of age of public apartment on ret levels ###

#Change in return levels (per unit change in spatial cov: building age of BEST MODEL 
ret.lev.ms4.1t_brown = as.data.frame(ret.lev.ms4.1t_brown)
predicted_ms4.1t_brown_df_rl_hdbmedage <- cbind(df_xy, ret.lev.ms4.1t_brown)
predicted_ms4.1t_brown_df_rl_hdbmedage["V1"][predicted_ms4.1t_brown_df_rl_hdbmedage["V1"]<0] <- 0 #force 0 for neg values
predicted_ms4.1t_brown_df_rl_hdbmedage$V1<- round(predicted_ms4.1t_brown_df_rl_hdbmedage$V1,
                                                 0)
colnames(predicted_ms4.1t_brown_df_rl_hdbmedage)[4] = "Return Levels"
ms4.1t_brown_rl_hdbmedage = merge(spat, predicted_ms4.1t_brown_df_rl_hdbmedage, by.x="cell_id", by.y="cell_id")


#1 year increase
ms4.1t_brown_retlev_hdbmedage_1 = tm_shape(ms4.1t_brown_rl_hdbmedage) +
  tm_fill("Return Levels", palette = "Oranges", style="fixed",title="", colorNA=grey(0.9),
          breaks = c(14, 36, 46, 56, 66, 94), 
          textNA = "None recorded")+
  tm_layout(legend.text.size = 1, legend.title.size = 1.5,
            title="M",title.size = 2, frame=F, legend.show=F) 
ms4.1t_brown_retlev_hdbmedage_1


#2 year increase
ms4.1t_brown_retlev_hdbmedage_2 = tm_shape(ms4.1t_brown_rl_hdbmedage) +
  tm_fill("Return Levels", palette = "Oranges", style="fixed",title="", colorNA=grey(0.9),
          breaks = c(14, 36, 46, 56, 66, 94), 
          textNA = "None recorded")+
  tm_layout(legend.text.size = 1, legend.title.size = 1.5,
            title="N",title.size = 2, frame=F, legend.show=F) 
ms4.1t_brown_retlev_hdbmedage_2


#3 year increase
ms4.1t_brown_retlev_hdbmedage_3 = tm_shape(ms4.1t_brown_rl_hdbmedage) +
  tm_fill("Return Levels", palette = "Oranges", style="fixed",title="", colorNA=grey(0.9),
          breaks = c(14, 36, 46, 56, 66, 94), 
          textNA = "None recorded"
  )+
  tm_layout(legend.text.size = 1, legend.title.size = 1.5,
            title="O",title.size = 2, frame=F, legend.show=T) 
ms4.1t_brown_retlev_hdbmedage_3


#For MANUSCRIPT overleaf; combined into 1 plot
png("Per_unitchange_spatcov NEW.png", units="cm", width=40, height=50, res=300)
tmap_arrange(ms4.1t_brown_retlev_fw_1, ms4.1t_brown_retlev_fw_5,
             ms4.1t_brown_retlev_fw_10, 
             ms4.1t_brown_retlev_imp_1, ms4.1t_brown_retlev_imp_5,
             ms4.1t_brown_retlev_imp_10,
             ms4.1t_brown_retlev_veghmtc_1, ms4.1t_brown_retlev_veghmtc_5,
             ms4.1t_brown_retlev_veghmtc_10,
             ms4.1t_brown_retlev_totalpop_1, ms4.1t_brown_retlev_totalpop_5, 
             ms4.1t_brown_retlev_totalpop_10, 
             ms4.1t_brown_retlev_hdbmedage_1, ms4.1t_brown_retlev_hdbmedage_2,
             ms4.1t_brown_retlev_hdbmedage_3,
             ncol=3)
dev.off()


