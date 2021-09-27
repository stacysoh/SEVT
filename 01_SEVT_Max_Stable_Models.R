getwd()
setwd("C:/Users/user/Desktop/EVT")

rm(list=ls())

#####################################
### Load datasets and format data ###
#####################################

dengue_1320 = read.csv("DEN-Hexagon XY DBF 2013_2020_complete.csv") #1563 unique hexagons 
dengue_0712 = read.csv("DEN-Hexagon XY DBF 2007_2012_complete.csv") #1368 unique hexagons
dengue_0712 = dengue_0712[,-c(16:17)] #remove columns with the wrong X and Y coords
colnames(dengue_0712)[17]="X.output"  #rename X coord column to match 2013-2020 file (dengue_1320)
colnames(dengue_0712)[18]="Y.output"  #rename X coord column to match 2013-2020 file (dengue_1320)

dengue_0712$POSTAL = as.integer(dengue_0712$POSTAL) #change class for postal code to interger

dengue = dplyr::bind_rows(dengue_0712, dengue_1320) #combined unique hexagons: 1618

#Remove duplicated rows 
dup = duplicated(dengue$CASE_ID)
dengue$CASE_ID[dup] #remove 53 duplicated entries
dengue = dengue[!dup,]

#Load packages
library(dplyr)
library(tidyr)
library(SpatialExtremes)
library(extRemes)


#########################################
###Obtain max yearly counts by polygons## 
#########################################

#Group by onset year and week, and cell_id
df_1 = dengue %>%
  group_by(Onset_EYear, Onset_EWeek, cell_id) %>%
  dplyr::summarize(count = n())


#Obtain the yearly max counts for each cell_id and remove duplicates for cell_id
df_2 =  df_1  %>%
  group_by(Onset_EYear, cell_id) %>%
  dplyr::summarise_each(funs(max)) 
#Need to sort by Onset year and cell_id
df_2 <-df_2[order(df_2$Onset_EYear, df_2$cell_id ),]
df_2 = df_2[,-3] #remove corresponding e-week column for yearly max

filter(df_2, cell_id == 3137)  #check if properly formatted

#Transpose the unique cell_id to rows 
df_counts <- df_2 %>% spread(key = cell_id, value = count)
#Remove incomplete 2021 data
df_counts<-df_counts[!(df_counts$Onset_EYear=="2006" | 
                         df_counts$Onset_EYear=="2021"),]
#Convert NAs to 0s
df_counts[is.na(df_counts)] <- 0

#Convert yearly max counts to matrix, remove year column
count_matrix = data.matrix(df_counts[,-1]) #SpatialExtremes package requirement

#Remove the intermediate dataframes
rm("df_1", "df_2", "dup")


####################
##Obtain XY coords##
####################

#Group by cell_id and XY
df_xy = dengue %>%
  group_by(cell_id, X.output, Y.output) %>%
  dplyr::distinct(cell_id) %>%
  relocate(cell_id, .before = X.output)

#Remove duplicated cell_id columns and reorder
df_xy = df_xy[!duplicated(df_xy$cell_id), ]
df_xy = df_xy[order(df_xy$cell_id, df_xy$X.output, df_xy$Y.output ),]

#Convert xy to matrix, remove cell_id column?
xy_matrix = cbind(df_xy$X.output, df_xy$Y.output)
colnames(xy_matrix) <- cbind("xcoord", "ycoord")
plot(df_xy$X.output, df_xy$Y.output) #visualize XY coords in plot


##############################
##Obtain temporal covariates##
##############################

###########Serotype data#######
sero = read.csv("Serotyped cases.csv")
sero[is.na(sero)] <- 0 #Convert NAs to 0s

#Obtain the yearly average
df_sero = aggregate(sero[, 3:6], list(sero$Year), mean)
colnames(df_sero)[1] <- "Year" #rename the year column

#Match the data timeframe
df_sero <- df_sero[!(df_sero$Year<="2006" | 
                       df_sero$Year=="2021"),]
#Convert serotype data to matrix, remove year column
sero_matrix = data.matrix(df_sero[,-1])
colnames(sero_matrix) <- c("d1", "d2", "d3", "d4")


##########Weather data##########
weather <- read.csv("Climate_2007_2020.csv")

temp_RH = aggregate(cbind(MeanT, MaxT, MinT, RH, AH)~Year, 
                    FUN=mean, data=weather, na.rm=TRUE) #obtain yearly mean
rain = aggregate(Daily.Rainfall~Year, 
                 FUN=sum, weather, na.rm=TRUE) #obtain yearly cumulative for rainfall
colnames(rain)[2] = "Rainfall"
weather_yearly = merge(temp_RH, rain, by=c("Year"))
weather_matrix = data.matrix(weather_yearly[,-1])


###############################
###Obtain spatial covariates###
#Obtain unique hexagon cell_ids#
################################
landsurface = read.csv("land surfaces.csv")

hexagon = distinct(dengue, cell_id) #obtain unique cell_ids: 1618
hexagon<-hexagon[order(hexagon$cell_id),] #sort in ascending order
hexagon = as.data.frame(hexagon)
colnames(hexagon)[1] = "cell_id"

ls_all = merge(hexagon, landsurface, by=c("cell_id") ) #1618 rows (n.sites); marg.cov


###Population data
popdata = read.csv("pop_data_2020.csv")
popdata$totalpop = popdata$scpr + popdata$nr #singapore citizens + non-residents
popdata_hex = (popdata[,-c(1,3:43)]) #remove all other columns

ls_all_1 = merge(ls_all, popdata_hex, by=c("cell_id") ) #1618 rows (n.sites); marg.cov

##Building age 
hdbage = read.csv("median building age.csv")
nrow(distinct(hdbage, cell_id)) #obtain unique cell_ids: 1483
hdbage = hdbage %>% 
  filter(Year==2019) %>% 
  distinct(cell_id, HDB_med_age) #only left with 1433 unique cell_ids
hdbage$HDB_med_age_2020 = hdbage$HDB_med_age 

hdbage$HDB_med_age_2020[
  hdbage$HDB_med_age_2020 >= 1] <- hdbage$HDB_med_age_2020[
    hdbage$HDB_med_age_2020 >= 1] + 1 #add 1 year to get the 2020 median building age
hdbage = (hdbage[,-c(2)]) #remove column with median building age as of 2019

ls_all_2 = merge(ls_all_1, hdbage, by=c("cell_id"), all.x=T) #merge spatial covariates together
ls_all_2[is.na(ls_all_2)] <- 0 #Convert NAs to 0s


#Convert spat cov data to matrix, remove cell_id column (col 1)
ls_v1 = data.matrix(ls_all_2[,-c(1,2,6,7)]) #remove irrelvant land cover spatial covariates
#remove marine (mosquitoes don't breed in seas) - col 2, 
#remove freshwater swamps and marshes - col 6 , 
#remove mangroves (only in southern islands and lck)- col 7
#retain freshwater surfaces, impervious surfaces
#retain non-vegetated pervious surfaces, vegetation with structure dominated by human management with tree canopy
#retain vegetation with structure dominated by human management without tree canopy, 
#retain vegetation with limited human management with tree canopy
#retain vegetation with limited human management with tree canopy



#######################################################################
### Max Stable Models (best 4 models with lowest pairwise deviance) ###
#######################################################################

### Modelling the trend surfaces (location, scale and shape parameters) ###

### ms8.1t ###
loc.form.8.1 <- y ~ xcoord + ls_v1 #location parameter
scale.form.8.1 <- y ~  xcoord*ycoord + ls_v1 #scale parameter
shape.form.8.1 <- y ~ 1 #shape parameter taken to be constant

#Schlather model with powered exponetial covariance function; takes approx. 17 hours
ms8.1t <- fitmaxstab(data=count_matrix, coord=xy_matrix, temp.cov=weather_matrix,
                     cov.mod="powexp", loc.form.8.1, scale.form.8.1, shape.form.8.1) 
predicted_ms8.1t <- predict(ms8.1t, ret.per=30) #to produce 30-year return levels

#Brown-resnick model: takes approx. 17 hours 
ms8.1t_brown <- fitmaxstab(data=count_matrix, coord=xy_matrix, temp.cov=weather_matrix,
                           cov.mod="brown", loc.form.8.1, scale.form.8.1, shape.form.8.1)
predicted_ms8.1t_brown <- predict(ms8.1t_brown, ret.per=30) #to produce 30-year return levels



### Modelling the trend surfaces (location, scale and shape parameters) ###

### ms 4.1 ###
loc.form.4.1 <- y ~ xcoord + ls_v1 #location parameter
scale.form.4.1 <- y ~  xcoord + ycoord + ls_v1 #scale parameter
shape.form.4.1 <- y ~ 1 #shape parameter taken to be constant

#Schlather model with powered exponetial covariance function; takes approx. 17 hours
ms4.1t <- fitmaxstab(data=count_matrix, coord=xy_matrix, temp.cov=weather_matrix,
                     cov.mod="powexp", loc.form.4.1, scale.form.4.1, shape.form.4.1) 
predicted_ms4.1t <- predict(ms4.1t, ret.per=30) #to produce 30-year return levels  

#Brown-resnick model: takes approx. 17 hours 
ms4.1t_brown <- fitmaxstab(data=count_matrix, coord=xy_matrix, temp.cov=weather_matrix,
                           cov.mod="brown", loc.form.4.1, scale.form.4.1, shape.form.4.1) 
predicted_ms4.1t_brown <- predict(ms4.1t_brown, ret.per=30) #to produce 30-year return levels



##########################################
### Max Stable Models (8 other models) ###
##########################################

### Modelling the trend surfaces (location, scale and shape parameters) ###

### ms10.1t ###
loc.form.10.1 <- y ~ xcoord*ycoord + ls_v1 #location parameter
scale.form.10.1 <- y ~  xcoord + ls_v1 #scale parameter
shape.form.10.1 <- y ~ 1 #shape parameter taken as constant

#Schlather model with powered exponetial covariance function
ms10.1t <- fitmaxstab(data=count_matrix, coord=xy_matrix, temp.cov=weather_matrix,
                      cov.mod="powexp", loc.form.10.1, scale.form.10.1, shape.form.10.1) 
predicted_ms10.1t <- predict(ms10.1t, ret.per=30) #to produce 30-year return levels

#Brown model
ms10.1t_brown <- fitmaxstab(data=count_matrix, coord=xy_matrix, temp.cov=weather_matrix,
                            cov.mod="brown", loc.form.10.1, scale.form.10.1, shape.form.10.1) 
predicted_ms10.1t_brown <- predict(ms10.1t_brown, ret.per=30) #to produce 30-year return levels

#Smith model 
ms10.1t_smith <- fitmaxstab(data=count_matrix, coord=xy_matrix, temp.cov=weather_matrix,
                            cov.mod="gauss", loc.form.10.1, scale.form.10.1, shape.form.10.1) 
predicted_ms10.1t_smith <- predict(ms10.1t_smith, ret.per=30) #to produce 30-year return levels


### Modelling the trend surfaces (location, scale and shape parameters) ###

### ms5.1t ###
loc.form.5.1 <- y ~ xcoord + ls_v1 #location parameter
scale.form.5.1 <- y ~  xcoord + ls_v1 #scale parameter
shape.form.5.1 <- y ~ 1 #shape parameter

#Schlather model 
ms5.1t <- fitmaxstab(data=count_matrix, coord=xy_matrix, temp.cov=weather_matrix,
                     cov.mod="powexp", loc.form.5.1, scale.form.5.1, shape.form.5.1) 
predicted_ms5.1t <- predict(ms5.1t, ret.per=30) #to produce 30-year return levels

#Brown-resnick model
ms5.1t_brown <- fitmaxstab(data=count_matrix, coord=xy_matrix, temp.cov=weather_matrix,
                           cov.mod="brown", loc.form.5.1, scale.form.5.1, shape.form.5.1) 
predicted_ms5.1t_brown <- predict(ms5.1t_brown, ret.per=30) #to produce 30-year return levels

#Smith model 
ms5.1t_smith <- fitmaxstab(data=count_matrix, coord=xy_matrix, temp.cov=weather_matrix,
                           cov.mod="gauss", loc.form.5.1, scale.form.5.1, shape.form.5.1) 
predicted_ms5.1t_smith <- predict(ms5.1t_smith, ret.per=30) #to produce 30-year return levels


#Smith model
ms8.1t_smith <- fitmaxstab(data=count_matrix, coord=xy_matrix, temp.cov=weather_matrix,
                           cov.mod="gauss", loc.form.8.1, scale.form.8.1, shape.form.8.1) 
predicted_ms8.1t_smith <- predict(ms8.1t_smith, ret.per=30) #to produce 30-year return levels

#Smith model 
ms4.1t_smith <- fitmaxstab(data=count_matrix, coord=xy_matrix, temp.cov=weather_matrix,
                           cov.mod="gauss", loc.form.4.1, scale.form.4.1, shape.form.4.1) 
predicted_ms4.1t_smith <- predict(ms4.1t_smith, ret.per=30) #to produce 30-year return levels 


