setwd("C:/Users/user/Desktop/EVT")
rm(list=ls())

#Load the [new] EVT_results inside
load("C:/Users/user/Desktop/EVT/R Scripts and R.data/[new] EVT_results NEW.RData")

################################################
### Plot spatial covariates and dengue cases ###
################################################

library(rgeos)
library(rgdal)
library(tmap)
library(grid)
library(dplyr)

#Load the hexagons shapefile
spat <- readOGR(dsn="C:/Users/user/Desktop/EVT/hexagonsData", layer="hexagonsData")
spat$v2 <- round(spat$v2, 2)
spat$v3 <- round(spat$v3, 2)
spat$v4 <- round(spat$v4, 2)
spat$v7 <- round(spat$v7, 2)
spat$v8 <- round(spat$v8, 2)
spat$v9 <- round(spat$v9, 2)
spat$v10 <- round(spat$v10, 2)


#Plots for the land covers

#Freshwater
plotv2 = tm_shape(spat) +
  tm_fill("v2", palette = "Blues", style="jenks", title="Proportion") +
  tm_layout(legend.text.size = 1, legend.title.size = 1.5,
            title="B",title.size = 2, frame=F)

#Impervious
plotv3 = tm_shape(spat) +
  tm_fill("v3", palette = "BuPu", style="jenks", title="Proportion") +
  tm_layout(legend.text.size = 1, legend.title.size = 1.5,
            title="C",title.size = 2, frame=F)

#Non-vegetated pervious surfaces
plotv4 = tm_shape(spat) +
  tm_fill("v4", palette = "RdPu", style="jenks", title="Proportion") +
  tm_layout(legend.text.size = 1.1, legend.title.size = 1.5,
            title="D",title.size = 2, frame=F)

#Vegetation with structure dominated by human management (w Tree Canopy)
plotv7 = tm_shape(spat) +
  tm_fill("v7", palette = "YlGn", style="jenks", title="Proportion") +
  tm_layout(legend.text.size = 1, legend.title.size = 1.5,
            title="E",title.size = 2, frame=F)

#Vegetation with structure dominated by human management (w/o Tree Canopy)
plotv8 = tm_shape(spat) +
  tm_fill("v8", palette = "YlGnBu", style="jenks", title="Proportion") +
  tm_layout(legend.text.size = 1, legend.title.size = 1.5,
            title="F",title.size = 2, frame=F)

#Vegetation with limited human management (w Tree Canopy)
plotv9 = tm_shape(spat) +
  tm_fill("v9", palette = "Greens", style="jenks", title="Proportion") +
  tm_layout(legend.text.size = 1, legend.title.size = 1.5,
            title="G",title.size = 2, frame=F)

#Vegetation with limited human management (w/o Tree Canopy)
plotv10 = tm_shape(spat) +
  tm_fill("v10", palette = "BuGn", style="jenks", title="Proportion") +
  tm_layout(legend.text.size = 1, legend.title.size = 1.5,
            title="H",title.size = 2, frame=F)

#Obtain the average maximum weekly dengue counts (2007-2020)
df_counts_t = t(df_counts)
df_counts_t_1 =df_counts_t[-1,]
df_counts_t_1 = rowMeans(df_counts_t_1)
df_counts_t_2= round(df_counts_t_1,digits=0)
df_counts_avemax  <- cbind(df_xy, df_counts_t_2)
colnames(df_counts_avemax)[4] = "ave_max_counts"
#Merge the average maximum weekly dengue counts (2007-2020) with the Hexagons shapefile
spat_merged = merge(spat, df_counts_avemax, by.x="cell_id", by.y="cell_id")

#Plot the average maximum weekly dengue counts (2007-2020)
plot_aveMaxDen = tm_shape(spat_merged) +
  tm_fill("ave_max_counts", style = "fixed", breaks = c(0, 1, 2, 3, 9), title="Count", 
          colorNA = grey(0.9), textNA = "None recorded") +
  tm_layout(legend.text.size = 1, legend.title.size = 1.5,
            title="A",title.size = 2, frame=F)


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


hdbage_plot = merge(spat, hdbage, by.x="cell_id", by.y="cell_id")
hdbage_plot_1 = tm_shape(hdbage_plot) +
  tm_fill("HDB_med_age_2020", palette = "Oranges", style="fixed", title="Age", 
          colorNA =grey(0.9), textNA = "None recorded",
          breaks = c(0, 20, 40, 60, 83)) +
  tm_layout(legend.text.size = 1, legend.title.size = 1.5,
            title="I",title.size = 2, frame=F)


###Population data
popdata = read.csv("pop_data_2020.csv")
popdata$totalpop = popdata$scpr + popdata$nr #singapore citizens + non-residents
popdata_hex = (popdata[,-c(1,3:43)]) #remove all other cols
popdata_hex_plot = merge(spat, popdata_hex, by.x="cell_id", by.y="cell_id")
popdata_hex_plot_1 = tm_shape(popdata_hex_plot) +
  tm_fill("totalpop", palette = "Purples", style="fixed", title="Count", 
          colorNA =grey(0.9), textNA = "None recorded",
          breaks = c(0, 1001, 4001, 7001, 11001, 17000)) +
  tm_layout(legend.text.size = 1, legend.title.size = 1.5,
            title="J",title.size = 2, frame=F)


png("Dengue_spatialcovariates.png", units="cm", width=50, height=20, res=300)
tmap_arrange(plot_aveMaxDen, plotv2, plotv3, plotv4, plotv7, 
             plotv8, plotv9, plotv10, hdbage_plot_1,popdata_hex_plot_1,
             ncol=5)
dev.off()


###########################################################################
### Coefficient plot of the paramters of the model with lowest deviance ###
############## plots the location and scale parameters ####################
###################### model ms4.1t_brown #################################

predicted_ms4.1t_brown = predict(ms4.1t_brown, ret.per=30)

#Location parameter
predicted_ms4.1t_brown_df <- as.data.frame(predicted_ms4.1t_brown)
predicted_ms4.1t_brown_df_loc <- cbind(df_xy, predicted_ms4.1t_brown_df$loc)
colnames(predicted_ms4.1t_brown_df_loc)[4] = "location"
predicted_ms4.1t_brown_df_loc$location<- round(predicted_ms4.1t_brown_df_loc$location
                                               , 2)
ms4.1t_brown_location = merge(spat, predicted_ms4.1t_brown_df_loc, by.x="cell_id", by.y="cell_id")

#Location plot
ms4.1t_brown_loc = tm_shape(ms4.1t_brown_location) +
  tm_fill("location", palette = "BuPu", style="jenks",
          breaks = c(0, 5, 14, 22, 29, 36, 51), title="", colorNA =grey(0.9),
          textNA = "None recorded") +
  tm_layout(legend.text.size = 1, legend.title.size = 1.5,
            title="A",title.size = 2, frame=F, legend.show=T) 

#Scale parameter
predicted_ms4.1t_brown_df <- as.data.frame(predicted_ms4.1t_brown)
predicted_ms4.1t_brown_df_scale <- cbind(df_xy, predicted_ms4.1t_brown_df$scale)
colnames(predicted_ms4.1t_brown_df_scale)[4] = "scale"
predicted_ms4.1t_brown_df_scale$scale<- round(predicted_ms4.1t_brown_df_scale$scale
                                               , 2)
ms4.1t_brown_scale = merge(spat, predicted_ms4.1t_brown_df_scale, by.x="cell_id", by.y="cell_id")

#Scale plot
ms4.1t_brown_sca = tm_shape(ms4.1t_brown_scale) +
  tm_fill("scale", palette = "BuGn", style="jenks",
          breaks = c(0, 5, 14, 22, 29, 36, 51), title="", colorNA =grey(0.9),
          textNA = "None recorded") +
  tm_layout(legend.text.size = 1, legend.title.size = 1.5,
            title="B",title.size = 2, frame=F, legend.show=T) 
ms4.1t_brown_sca


png("loc_scale_NEW.png", units="cm", width=20, height=10, res=300)
tmap_arrange(ms4.1t_brown_loc, ms4.1t_brown_sca,
             ncol=2)
dev.off()


#########################################################################
### Return level plots for the 4 best models (with lowest deviance) #####
#########################################################################

## Best model with lowest paired deviance: ms4.1t brown ###

#5-year return levels
predicted_ms4.1t_brown_5<- predict(ms4.1t_brown, ret.per=5) 
predicted_ms4.1t_brown_5 <- as.data.frame(predicted_ms4.1t_brown_5)
predicted_ms4.1t_brown_5 <- cbind(df_xy, predicted_ms4.1t_brown_5$Q5)
colnames(predicted_ms4.1t_brown_5)[4] = "Q5"

predicted_ms4.1t_brown_5$Q5 <- round(predicted_ms4.1t_brown_5$Q5, 0) #round to integer
predicted_ms4.1t_brown_5["Q5"][predicted_ms4.1t_brown_5["Q5"]<0] <- 0 #force 0 for neg values

predicted_ms4.1t_brown_5 = merge(spat, predicted_ms4.1t_brown_5, by.x="cell_id", by.y="cell_id")
predicted_ms4.1t_brown_5_plot = tm_shape(predicted_ms4.1t_brown_5) +
  tm_fill("Q5", palette = "Reds", style="fixed",
          breaks = c(0, 10, 24, 42, 64, 92), title="", colorNA =grey(0.9),
          textNA = "None recorded") +
  tm_layout(legend.text.size = 1, legend.title.size = 1.5,
            title="M1 [A]",title.size = 1.5, frame=F, legend.show=F) 
predicted_ms4.1t_brown_5_plot

#10-year return levels
predicted_ms4.1t_brown_10<- predict(ms4.1t_brown, ret.per=10) 
predicted_ms4.1t_brown_10 <- as.data.frame(predicted_ms4.1t_brown_10)
predicted_ms4.1t_brown_10 <- cbind(df_xy, predicted_ms4.1t_brown_10$Q10)
colnames(predicted_ms4.1t_brown_10)[4] = "Q10"

predicted_ms4.1t_brown_10$Q10 <- round(predicted_ms4.1t_brown_10$Q10, 0) #round to integer
predicted_ms4.1t_brown_10["Q10"][predicted_ms4.1t_brown_10["Q10"]<0] <- 0 #force 0 for neg values

predicted_ms4.1t_brown_10 = merge(spat, predicted_ms4.1t_brown_10, by.x="cell_id", by.y="cell_id")
predicted_ms4.1t_brown_10_plot = tm_shape(predicted_ms4.1t_brown_10) +
  tm_fill("Q10", palette = "Reds", style="fixed",
          breaks = c(0, 10, 24, 42, 64, 92), title="", colorNA =grey(0.9),
          textNA = "None recorded") +
  tm_layout(legend.text.size = 1, legend.title.size = 1.5,
            title="M1 [B]",title.size = 1.5, frame=F, legend.show=F)  
predicted_ms4.1t_brown_10_plot

#20-year return levels
predicted_ms4.1t_brown_20<- predict(ms4.1t_brown, ret.per=20) 
predicted_ms4.1t_brown_20 <- as.data.frame(predicted_ms4.1t_brown_20)
predicted_ms4.1t_brown_20 <- cbind(df_xy, predicted_ms4.1t_brown_20$Q20)
colnames(predicted_ms4.1t_brown_20)[4] = "Q20"

predicted_ms4.1t_brown_20$Q20 <- round(predicted_ms4.1t_brown_20$Q20, 0) #round to integer
predicted_ms4.1t_brown_20["Q20"][predicted_ms4.1t_brown_20["Q20"]<0] <- 0 #force 0 for neg values

predicted_ms4.1t_brown_20 = merge(spat, predicted_ms4.1t_brown_20, by.x="cell_id", by.y="cell_id")
predicted_ms4.1t_brown_20_plot = tm_shape(predicted_ms4.1t_brown_20) +
  tm_fill("Q20", palette = "Reds", style="fixed",
          breaks = c(0, 10, 24, 42, 64, 92), title="", colorNA =grey(0.9),
          textNA = "None recorded") +
  tm_layout(legend.text.size = 1, legend.title.size = 1.5,
            title="M1 [C]",title.size = 1.5, frame=F, legend.show=F) 
predicted_ms4.1t_brown_20_plot

#30-year return levels
predicted_ms4.1t_brown_30<- predict(ms4.1t_brown, ret.per=30) 
predicted_ms4.1t_brown_30 <- as.data.frame(predicted_ms4.1t_brown_30)
predicted_ms4.1t_brown_30 <- cbind(df_xy, predicted_ms4.1t_brown_30$Q30)
colnames(predicted_ms4.1t_brown_30)[4] = "Q30"

predicted_ms4.1t_brown_30$Q30 <- round(predicted_ms4.1t_brown_30$Q30, 0) #round to integer
predicted_ms4.1t_brown_30["Q30"][predicted_ms4.1t_brown_30["Q30"]<0] <- 0 #force 0 for neg values

predicted_ms4.1t_brown_30 = merge(spat, predicted_ms4.1t_brown_30, by.x="cell_id", by.y="cell_id")
predicted_ms4.1t_brown_30_plot = tm_shape(predicted_ms4.1t_brown_30) +
  tm_fill("Q30", palette = "Reds", style="fixed",
          breaks = c(0, 10, 24, 42, 64, 92), title="", colorNA =grey(0.9),
          textNA = "None recorded") +
  tm_layout(legend.text.size = 1.1, legend.title.size = 1.5,
            title="M1 [D]",title.size = 1.5, frame=F, legend.show=T) 
predicted_ms4.1t_brown_30_plot #0, 12, 30, 55, 64, 92

tiff("Map_ms4.1t_brown_ret_level_NEW.tiff", units="cm", width=40, height=10, res=300)
tmap_arrange(predicted_ms4.1t_brown_5_plot, predicted_ms4.1t_brown_10_plot, 
             predicted_ms4.1t_brown_20_plot, predicted_ms4.1t_brown_30_plot,
             ncol=4) #to fit dimensions for ppt slides
dev.off()


### Second best model: ms8.1t sclather ###

#5-year return levels
predicted_ms8.1t_5<- predict(ms8.1t, ret.per=5) 
predicted_ms8.1t_5 <- as.data.frame(predicted_ms8.1t_5)
predicted_ms8.1t_5 <- cbind(df_xy, predicted_ms8.1t_5$Q5)
colnames(predicted_ms8.1t_5)[4] = "Q5"

predicted_ms8.1t_5$Q5 <- round(predicted_ms8.1t_5$Q5, 0) #round to integer
predicted_ms8.1t_5["Q5"][predicted_ms8.1t_5["Q5"]<0] <- 0 #force 0 for neg values

predicted_ms8.1t_5 = merge(spat, predicted_ms8.1t_5, by.x="cell_id", by.y="cell_id")
predicted_ms8.1t_5_plot = tm_shape(predicted_ms8.1t_5) +
  tm_fill("Q5", palette = "Reds", style="fixed", title="", colorNA =grey(0.9),
          breaks = c(0, 5, 10, 25, 45, 75), 
          textNA = "None recorded") +
  tm_layout(legend.text.size = 1, legend.title.size = 1.5,
            title="M2 [E]",title.size = 1.5, frame=F, legend.show=F) 
predicted_ms8.1t_5_plot

#10-year return levels
predicted_ms8.1t_10<- predict(ms8.1t, ret.per=10) 
predicted_ms8.1t_10 <- as.data.frame(predicted_ms8.1t_10)
predicted_ms8.1t_10 <- cbind(df_xy, predicted_ms8.1t_10$Q10)
colnames(predicted_ms8.1t_10)[4] = "Q10"

predicted_ms8.1t_10$Q10 <- round(predicted_ms8.1t_10$Q10, 0) #round to integer
predicted_ms8.1t_10["Q10"][predicted_ms8.1t_10["Q10"]<0] <- 0 #force 0 for neg values

predicted_ms8.1t_10 = merge(spat, predicted_ms8.1t_10, by.x="cell_id", by.y="cell_id")
predicted_ms8.1t_10_plot = tm_shape(predicted_ms8.1t_10) +
  tm_fill("Q10", palette = "Reds", style="fixed",title="", colorNA =grey(0.9),
          breaks = c(0, 5, 10, 25, 45, 75), 
          textNA = "None recorded") +
  tm_layout(legend.text.size = 1, legend.title.size = 1.5,
            title="M2 [F]",title.size = 1.5, frame=F, legend.show=F)  
predicted_ms8.1t_10_plot

#20-year return levels
predicted_ms8.1t_20<- predict(ms8.1t, ret.per=20) 
predicted_ms8.1t_20 <- as.data.frame(predicted_ms8.1t_20)
predicted_ms8.1t_20 <- cbind(df_xy, predicted_ms8.1t_20$Q20)
colnames(predicted_ms8.1t_20)[4] = "Q20"

predicted_ms8.1t_20$Q20 <- round(predicted_ms8.1t_20$Q20, 0) #round to integer
predicted_ms8.1t_20["Q20"][predicted_ms8.1t_20["Q20"]<0] <- 0 #force 0 for neg values

predicted_ms8.1_20 = merge(spat, predicted_ms8.1t_20, by.x="cell_id", by.y="cell_id")
predicted_ms8.1t_20_plot = tm_shape(predicted_ms8.1_20) +
  tm_fill("Q20", palette = "Reds", style="fixed",title="", colorNA =grey(0.9),
          breaks = c(0, 5, 10, 25, 45, 75), 
          textNA = "None recorded") +
  tm_layout(legend.text.size = 1, legend.title.size = 1.5,
            title="M2 [G]",title.size = 1.5, frame=F, legend.show=F) 
predicted_ms8.1t_20_plot

#30-year return levels
predicted_ms8.1t_30<- predict(ms8.1t, ret.per=30) 
predicted_ms8.1t_30 <- as.data.frame(predicted_ms8.1t_30)
predicted_ms8.1t_30 <- cbind(df_xy, predicted_ms8.1t_30$Q30)
colnames(predicted_ms8.1t_30)[4] = "Q30"

predicted_ms8.1t_30$Q30 <- round(predicted_ms8.1t_30$Q30, 0) #round to integer
predicted_ms8.1t_30["Q30"][predicted_ms8.1t_30["Q30"]<0] <- 0 #force 0 for neg values

predicted_ms8.1t_30 = merge(spat, predicted_ms8.1t_30, by.x="cell_id", by.y="cell_id")
predicted_ms8.1t_30_plot = tm_shape(predicted_ms8.1t_30) +
  tm_fill("Q30", palette = "Reds", style="fixed",title="", colorNA =grey(0.9),
          breaks = c(0, 5, 10, 25, 45, 75), #0, 10, 24, 42, 64, 92
          textNA = "None recorded") +
  tm_layout(legend.text.size = 1.1, legend.title.size = 1.5,
            title="M2 [H]",title.size = 1.5, frame=F, legend.show=T) 
predicted_ms8.1t_30_plot



### Third best model:  ms4.1 schlather ###

#5-year return levels
predicted_ms4.1t_5<- predict(ms4.1t, ret.per=5) 
predicted_ms4.1t_5 <- as.data.frame(predicted_ms4.1t_5)
predicted_ms4.1t_5 <- cbind(df_xy, predicted_ms4.1t_5$Q5)
colnames(predicted_ms4.1t_5)[4] = "Q5"

predicted_ms4.1t_5$Q5 <- round(predicted_ms4.1t_5$Q5, 0) #round to integer
predicted_ms4.1t_5["Q5"][predicted_ms4.1t_5["Q5"]<0] <- 0 #force 0 for neg values

predicted_ms4.1t_5 = merge(spat, predicted_ms4.1t_5, by.x="cell_id", by.y="cell_id")
predicted_ms4.1t_5_plot = tm_shape(predicted_ms4.1t_5) +
  tm_fill("Q5", palette = "Reds", style="fixed",
          breaks = c(0, 10, 26, 44, 66, 96), title="", colorNA =grey(0.9),
          textNA = "None recorded") +
  tm_layout(legend.text.size = 1, legend.title.size = 1.5,
            title="M3 [I]",title.size = 1.5, frame=F, legend.show=F) 
predicted_ms4.1t_5_plot

#10-year return levels
predicted_ms4.1t_10<- predict(ms4.1t, ret.per=10) 
predicted_ms4.1t_10 <- as.data.frame(predicted_ms4.1t_10)
predicted_ms4.1t_10 <- cbind(df_xy, predicted_ms4.1t_10$Q10)
colnames(predicted_ms4.1t_10)[4] = "Q10"

predicted_ms4.1t_10$Q10 <- round(predicted_ms4.1t_10$Q10, 0) #round to integer
predicted_ms4.1t_10["Q10"][predicted_ms4.1t_10["Q10"]<0] <- 0 #force 0 for neg values

predicted_ms4.1t_10 = merge(spat, predicted_ms4.1t_10, by.x="cell_id", by.y="cell_id")
predicted_ms4.1t_10_plot = tm_shape(predicted_ms4.1t_10) +
  tm_fill("Q10", palette = "Reds", style="fixed",
          breaks = c(0, 10, 26, 44, 66, 96), title="", colorNA =grey(0.9),
          textNA = "None recorded") +
  tm_layout(legend.text.size = 1, legend.title.size = 1.5,
            title="M3 [J]",title.size = 1.5, frame=F, legend.show=F)  
predicted_ms4.1t_10_plot

#20-year return levels
predicted_ms4.1t_20<- predict(ms4.1t, ret.per=20) 
predicted_ms4.1t_20 <- as.data.frame(predicted_ms4.1t_20)
predicted_ms4.1t_20 <- cbind(df_xy, predicted_ms4.1t_20$Q20)
colnames(predicted_ms4.1t_20)[4] = "Q20"

predicted_ms4.1t_20$Q20 <- round(predicted_ms4.1t_20$Q20, 0) #round to integer
predicted_ms4.1t_20["Q20"][predicted_ms4.1t_20["Q20"]<0] <- 0 #force 0 for neg values

predicted_ms4.1t_20 = merge(spat, predicted_ms4.1t_20, by.x="cell_id", by.y="cell_id")
predicted_ms4.1t_20_plot = tm_shape(predicted_ms4.1t_20) +
  tm_fill("Q20", palette = "Reds", style="fixed",
          breaks = c(0, 10, 26, 44, 66, 96), title="", colorNA =grey(0.9),
          textNA = "None recorded") +
  tm_layout(legend.text.size = 1, legend.title.size = 1.5,
            title="M3 [K]",title.size = 1.5, frame=F, legend.show=F) 
predicted_ms4.1t_20_plot

#30-year return levels
predicted_ms4.1t_30<- predict(ms4.1t, ret.per=30) 
predicted_ms4.1t_30 <- as.data.frame(predicted_ms4.1t_30)
predicted_ms4.1t_30 <- cbind(df_xy, predicted_ms4.1t_30$Q30)
colnames(predicted_ms4.1t_30)[4] = "Q30"

predicted_ms4.1t_30$Q30 <- round(predicted_ms4.1t_30$Q30, 0) #round to integer
predicted_ms4.1t_30["Q30"][predicted_ms4.1t_30["Q30"]<0] <- 0 #force 0 for neg values

predicted_ms4.1t_30 = merge(spat, predicted_ms4.1t_30, by.x="cell_id", by.y="cell_id")
predicted_ms4.1t_30_plot = tm_shape(predicted_ms4.1t_30) +
  tm_fill("Q30", palette = "Reds", style="fixed",
          breaks = c(0, 10, 26, 44, 66, 96), title="", colorNA =grey(0.9),
          textNA = "None recorded") +
  tm_layout(legend.text.size = 1.1, legend.title.size = 1.5,
            title="M3 [L]",title.size = 1.5, frame=F, legend.show=T) 
predicted_ms4.1t_30_plot #0, 10, 24, 42, 64, 92


### Fourth best model: ms8.1t brown ###

#5-year return levels
predicted_ms8.1t_brown_5<- predict(ms8.1t_brown, ret.per=5) 
predicted_ms8.1t_brown_5 <- as.data.frame(predicted_ms8.1t_brown_5)
predicted_ms8.1t_brown_5 <- cbind(df_xy, predicted_ms8.1t_brown_5$Q5)
colnames(predicted_ms8.1t_brown_5)[4] = "Q5"

predicted_ms8.1t_brown_5$Q5 <- round(predicted_ms8.1t_brown_5$Q5, 0) #round to integer
predicted_ms8.1t_brown_5["Q5"][predicted_ms8.1t_brown_5["Q5"]<0] <- 0 #force 0 for neg values

predicted_ms8.1t_brown_5 = merge(spat, predicted_ms8.1t_brown_5, by.x="cell_id", by.y="cell_id")
predicted_ms8.1t_brown_5_plot = tm_shape(predicted_ms8.1t_brown_5) +
  tm_fill("Q5", palette = "Reds", style="fixed", title="", colorNA =grey(0.9),
          breaks = c(0, 10, 15, 25, 40, 65), 
          textNA = "None recorded") +
  tm_layout(legend.text.size = 1, legend.title.size = 1.5,
            title="M4 [M]",title.size = 1.5, frame=F, legend.show=F) 
predicted_ms8.1t_brown_5_plot

#10-year return levels
predicted_ms8.1t_brown_10<- predict(ms8.1t_brown, ret.per=10) 
predicted_ms8.1t_brown_10 <- as.data.frame(predicted_ms8.1t_brown_10)
predicted_ms8.1t_brown_10 <- cbind(df_xy, predicted_ms8.1t_brown_10$Q10)
colnames(predicted_ms8.1t_brown_10)[4] = "Q10"

predicted_ms8.1t_brown_10$Q10 <- round(predicted_ms8.1t_brown_10$Q10, 0) #round to integer
predicted_ms8.1t_brown_10["Q10"][predicted_ms8.1t_brown_10["Q10"]<0] <- 0 #force 0 for neg values

predicted_ms8.1t_brown_10 = merge(spat, predicted_ms8.1t_brown_10, by.x="cell_id", by.y="cell_id")
predicted_ms8.1t_brown_10_plot = tm_shape(predicted_ms8.1t_brown_10) +
  tm_fill("Q10", palette = "Reds", style="fixed",title="", colorNA =grey(0.9),
          breaks = c(0, 10, 15, 25, 40, 65), 
          textNA = "None recorded") +
  tm_layout(legend.text.size = 1, legend.title.size = 1.5,
            title="M4 [N]",title.size = 1.5, frame=F, legend.show=F)  
predicted_ms8.1t_brown_10_plot

#20-year return levels
predicted_ms8.1t_brown_20<- predict(ms8.1t_brown, ret.per=20) 
predicted_ms8.1t_brown_20 <- as.data.frame(predicted_ms8.1t_brown_20)
predicted_ms8.1t_brown_20 <- cbind(df_xy, predicted_ms8.1t_brown_20$Q20)
colnames(predicted_ms8.1t_brown_20)[4] = "Q20"

predicted_ms8.1t_brown_20$Q20 <- round(predicted_ms8.1t_brown_20$Q20, 0) #round to integer
predicted_ms8.1t_brown_20["Q20"][predicted_ms8.1t_brown_20["Q20"]<0] <- 0 #force 0 for neg values

predicted_ms8.1t_brown_20 = merge(spat, predicted_ms8.1t_brown_20, by.x="cell_id", by.y="cell_id")
predicted_ms8.1t_brown_20_plot = tm_shape(predicted_ms8.1t_brown_20) +
  tm_fill("Q20", palette = "Reds", style="fixed",title="", colorNA =grey(0.9),
          breaks = c(0, 10, 15, 25, 40, 65), 
          textNA = "None recorded") +
  tm_layout(legend.text.size = 1, legend.title.size = 1.5,
            title="M4 [O]",title.size = 1.5, frame=F, legend.show=F) 
predicted_ms8.1t_brown_20_plot

#30-year return levels
predicted_ms8.1t_brown_30<- predict(ms8.1t_brown, ret.per=30) 
predicted_ms8.1t_brown_30 <- as.data.frame(predicted_ms8.1t_brown_30)
predicted_ms8.1t_brown_30 <- cbind(df_xy, predicted_ms8.1t_brown_30$Q30)
colnames(predicted_ms8.1t_brown_30)[4] = "Q30"

predicted_ms8.1t_brown_30$Q30 <- round(predicted_ms8.1t_brown_30$Q30, 0) #round to integer
predicted_ms8.1t_brown_30["Q30"][predicted_ms8.1t_brown_30["Q30"]<0] <- 0 #force 0 for neg values

predicted_ms8.1t_brown_30 = merge(spat, predicted_ms8.1t_brown_30, by.x="cell_id", by.y="cell_id")
predicted_ms8.1t_brown_30_plot = tm_shape(predicted_ms8.1t_brown_30) +
  tm_fill("Q30", palette = "Reds", style="fixed",title="", colorNA =grey(0.9),
          breaks = c(0, 10, 15, 25, 40, 65), #0, 5, 10, 25, 45, 75
          textNA = "None recorded") +
  tm_layout(legend.text.size = 1.1, legend.title.size = 1.5,
            title="M4 [P]",title.size = 1.5, frame=F, legend.show=T) 
predicted_ms8.1t_brown_30_plot


#combine plots together to produce a png file
png("Return level plots (4 best models - 5,10,20,30 years) NEW_1.png", 
     units="cm", width=60, height=40, res=200) #origianl: res=300
tmap_arrange(predicted_ms4.1t_brown_5_plot, predicted_ms4.1t_brown_10_plot, 
             predicted_ms4.1t_brown_20_plot, predicted_ms4.1t_brown_30_plot,
             predicted_ms8.1t_5_plot, predicted_ms8.1t_10_plot, 
             predicted_ms8.1t_20_plot, predicted_ms8.1t_30_plot,
             predicted_ms4.1t_5_plot, predicted_ms4.1t_10_plot,
             predicted_ms4.1t_20_plot, predicted_ms4.1t_30_plot,
             predicted_ms8.1t_brown_5_plot, predicted_ms8.1t_brown_10_plot,
             predicted_ms8.1t_brown_20_plot, predicted_ms8.1t_brown_30_plot,
             ncol=4)
dev.off()
