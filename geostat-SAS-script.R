# Geostatistical assessment of the highly structured spatial pattern of salt-affected soils
# 
# Author:       Gabor Szatmari and Laszlo Pasztor
# E-mail:       szatmari@rissac.hu
# Affiliation:  Institute for Soil Sciences, Centre for Agricultural Research, Budapest, Hungary


# 1. Necessary packages ----
library(spcosa)
library(readxl)
library(sp)
library(gstat)
library(maptools)
library(sp)
library(rgdal)
library(spatstat)
library(pROC)
library(ggplot2)
library(gridExtra)



# 2. Load ground truth ----
gt <- read_xlsx("data/GT.xlsx") # read ground truth

coordinates(gt) <- ~x+y
gridded(gt)=TRUE

spplot(gt["SAS"], scales=list(draw=TRUE), main="Ground truth") # show ground truth



# 3. Spatial coverage sampling ----
size <- c(50,100,200,300,400,500,600,700,800,900,1000) # sample sizes


set.seed(1234)
samples <- lapply(1:length(size), function(i){
  spcosa::spsample(stratify(gt, nStrata=size[i]))
}) # generate sampling designs... this may take a while (approx. 1 min)

plot(samples[[8]], scales=list(draw=TRUE)) # show one of the sampling designs

data <- lapply(1:length(size), function(i){
  data.frame(as(samples[[i]], "data.frame"), data.frame(over(x=as(samples[[i]], "SpatialPoints"), y=gt)))
})



# 4. The indicator approach ----
coordinates(data[[1]]) <- ~x+y # set coordinates
coordinates(data[[2]]) <- ~x+y
coordinates(data[[3]]) <- ~x+y
coordinates(data[[4]]) <- ~x+y
coordinates(data[[5]]) <- ~x+y
coordinates(data[[6]]) <- ~x+y
coordinates(data[[7]]) <- ~x+y
coordinates(data[[8]]) <- ~x+y
coordinates(data[[9]]) <- ~x+y
coordinates(data[[10]]) <- ~x+y
coordinates(data[[11]]) <- ~x+y


## 4.1. Variography ====
temp <- as(gt, "SpatialPointsDataFrame")

plot(variogram(I(SAS)~1, temp, width=200, cutoff=1000, map=TRUE)) # show variogram map... this may take a while (< 1 min)

v <- variogram(I(SAS)~1, temp, width=100, cutoff=1500, alpha=c(0, 45, 90, 135)) # compute directional indicator variograms... this may also take a while (< 1 min)

vgm <- fit.variogram(v, vgm(psill=0.11,
                            model="Exp",
                            range=500,
                            nugget=0.01,
                            anis=c(135,0.5)))

plot(v,vgm, main="Directional indicator variograms"); vgm # show the fitted variogram model and its parameters
rm(temp)


## 4.2. Indicator kriging ====
ik <- lapply(1:length(data), function(i){
  krige(I(SAS)~1, data[[i]], gt, vgm, debug.level=-1)
}) # perform indicator kriging... this may take a while (2-3 min)

# Post-processing indicator kriging:
ik <- lapply(1:length(ik), function(i){
  ik[[i]]@data[,1] <- ifelse(ik[[i]]@data[,1] < 0, yes=0, no=ik[[i]]@data[,1]) # set lower limit to 0
  ik[[i]]@data[,1] <- ifelse(ik[[i]]@data[,1] > 1, yes=1, no=ik[[i]]@data[,1]) # set upper limit to 1
  return(ik[[i]])
})


spplot(ik[[8]][1], scales=list(draw=TRUE), main="Indicator kriging (n=700)") # show one of the probability maps given by indicator kriging



# 5. Validation ----

## 5.1. ROC curve and AUC for indicator kriging ====
ik.roc <- lapply(1:length(ik), function(i){
  roc(response=gt@data$SAS,
      predictor=ik[[i]]@data[,1])
})

ik.auc <- lapply(1:length(ik.roc), function(i){
  auc(ik.roc[[i]])
})
ik.auc <- do.call(rbind, ik.auc)


## 5.2. ROC curve and AUC for single normal equation simulation ====
load("snesim.RData") # read SGeMS results

spplot(snesim["X.700"], scales=list(draw=TRUE), main="Single normal equation simulation (n=700)") # show one of the probability maps given by single normal equation simulation

snesim.roc <- lapply(2:ncol(snesim), function(i){
  roc(response=gt@data$SAS,
      predictor=snesim@data[,i])
})

snesim.auc <- lapply(1:length(snesim.roc), function(i){
  auc(snesim.roc[[i]])
})
snesim.auc <- do.call(rbind, snesim.auc)


# Show ROC curves for both geostatistical methods:
names(ik.roc) <- as.character(size)
ROC.curves.ik <- ggroc(ik.roc, legacy.axes=TRUE) +
  labs(x = "False-positive rate", y = "True-positive rate", title = "ROC curves for IK") +
  scale_color_discrete(name="Sample size") +
  theme_bw()

names(snesim.roc) <- as.character(size)
ROC.curves.snesim <- ggroc(snesim.roc, legacy.axes=TRUE) +
  labs(x = "False-positive rate", y = "True-positive rate", title = "ROC curves for SNESIM") +
  scale_color_discrete(name="Sample size") +
  theme_bw()

grid.arrange(ROC.curves.ik, ROC.curves.snesim, ncol=2)


## 5.3. AUC and the nearest neighbor distance as function of the sample size ====
mean.dist <- lapply(1:length(data), function(i){
  mean(nndist(as.ppp(data[[i]])))
}) # compute the average nearest neighbour distance for each sampling design

mean.dist <- unlist(mean.dist)

par(mfrow=c(1,1))
par(mar=c(5,5,1,5))
plot(x=size, y=ik.auc, type="b", pch=1, ylim=c(0.5, 1), xlab="Sample size", ylab="AUC", col="forestgreen", lwd=1.5, main="")
points(x=size, y=snesim.auc, type="b", xlab="", ylab="", col="orange", lwd=1.5, pch=2)
par(new=TRUE)
plot(x=size, y=mean.dist, pch=0, type="b", lty=2, xlab="", ylab="", axes=FALSE, lwd=1.5, col="grey50")
axis(side=4, at=pretty(range(mean.dist)))
mtext("Average nearest neighbour distance (meter)", side=4, line=3)
legend("topright", legend=c("IK", "SNESIM", "NNDIST"), col=c("forestgreen", "orange", "grey50"), lty=c(1,1,2), lwd=1.5, box.lty=0, bty="n")
