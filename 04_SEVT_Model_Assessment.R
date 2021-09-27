rm(list=ls())
setwd("C:/Users/user/Desktop/EVT")

#Load packages
library(SpatialExtremes)

#Load output from 01_SEVT_Max_Stable_Models first

################################################################
######### Model assessment: QQ plots and F-madogram ############
################################################################
## Using the best model with lowest deviance: ms4.1t_brown 

#####################################################
### Individual QQ plots (from 9 random locations) ###
#####################################################
n.obs <- nrow(ms4.1t_brown$data[, sites])
model <- ms4.1t_brown$model

range <- ms4.1t_brown$par["range"]
smooth <- ms4.1t_brown$par["smooth"]
sim.maxstab <- rmaxstab(n.obs*1000, ms4.1t_brown$coord[sites,],
                        ms4.1t_brown$cov.mod, range=range, smooth=smooth)
#this produces simulations from the fitted model, takes approx. 16 hours

#Transforms the simulations to the fitted model into an array with 5 sheets
sim.maxstab_1 <- array(log(sim.maxstab), c(n.obs, 1000, 5)) 

#Transforms to unit Fretchet margins and log transforms it to the gumbel scale from the fitted model
gumb <- log(apply(ms4.1t_brown$data[, sites], 2, gev2frech, emp = TRUE))

#Margins for manuscript plot 
par(mar = c(3, 3, 0.3, 0.3))  
par(mfrow=c(3,3))

for (i in 1:4) {
  for (j in (i + 1):5) {
    pair.max <- sort(apply(gumb[, c(i, j)], 1, max)) #observed pairwise max
    sim.pair.max <- apply(pmax(sim.maxstab_1[, , i], 
                               sim.maxstab_1[, , j]), 2, sort)#sim pairwise max
    dummy <- rowMeans(sim.pair.max) #sim pairwise max
    ci <- apply(sim.pair.max, 1, quantile, c(0.025, 
                                             0.975))#sim pairwise max
    matplot(dummy, t(ci), pch = "-", col = 1, 
            xlab = "", ylab = "", las=1) #sim pairwise max
    points(dummy, pair.max, pch=20, cex=0.8) #observed pairwise max
    abline(0, 1)
    #h <- distance(ms4.1t_brown$coord[sites[c(i, j)], ])
    #legend("bottomright", paste("h =", 
    #                            round(h, 2)), bty = "n")
    title(ylab="Observed", line=2)
    title(xlab="Modelled", line=2) ##for pdf, 5x5.5 for pdf inches
    
  }
} 



###################################################
####Individual QQ plots for all 1618 locations#####
###################################################
sim.maxstab_2 <- array(log(sim.maxstab), c(n.obs, 1000, 1618)) #takes awhile
#simulations to produce an array with 1618 dimensions

for (i in 1:1617) { 
  for (j in (i + 1):1618) {
    pair.max <- sort(apply(gumb[, ], 1, max))
    sim.pair.max <- apply(pmax(sim.maxstab_2[, , i], 
                               sim.maxstab_2[, , j]), 2, sort)
    dummy <- rowMeans(sim.pair.max)
    ci <- apply(sim.pair.max, 1, quantile, c(0.025, 
                                             0.975))
    matplot(dummy, t(ci), pch = "-", col = 1, 
            xlab = "Model", ylab = "Observed")
    points(dummy, pair.max)
    abline(0, 1) 
  }
}

##########################################
## Extremal coefficient and F madogram ###
##########################################
par(mfrow=c(1,1))
n.site = ncol(ms4.1t_brown$data)
fmadogram(ms4.1t_brown$data, ms4.1t_brown$coord, col="lightgrey") 
fmadogram(fitted=ms4.1t_brown, add=TRUE, n.bins=n.site) 

fmadogram(ms4.1t_brown$data, ms4.1t_brown$coord, which="ext", col="lightgrey") 
fmadogram(fitted=ms4.1t_brown, which="ext", add=TRUE, n.bins=n.site) 
