#THE PURPOSE OF THIS SCRIPT IS TO OUTPUT THE FIGURES AND THE RESULTS OF THE STATISTICAL TESTS FOUND IN THE MAIN TEXT

library(nlme)
library(multcomp)
library(correlation)

setwd("D:/Meerkat_TurnSpeed_Influence")
source("scripts/functions.R")

discretizationStep <- 10

load(paste0("output/spatialMetrics_",discretizationStep,"m.RData"))
load(paste0("output/Influence_logistic_model_fits_",discretizationStep,"m.RData"))

sessions <- c("HM2017","HM2019","L2019","ZU2021","NQ2021")
status = c("DominantF","DominantM","Adult","Yearling","Sub-Adult")

spatialMetrics$session <- factor(spatialMetrics$session,levels=sessions)

groupCol = hcl.colors(5)
statusCol=rev(c("#000000","#900C3F","#C70039","#FF5733","#FFC300"))
saveImage = T
cex=2.5

plotCIs = T
format = "png"

dir.create(path="figures/",showWarnings = F)

#----DISTRIBUTION OF TIME SPENT IN THE FRONT----

allInd = modelParam_MovTurn$ind
allIndShort <- gsub("20","",allInd)

totPropFront <- tapply(spatialMetrics$inFrontHalf,factor(spatialMetrics$indUniqID,levels=allInd),mean,na.rm=T)

windowSize = 60*60
windowPropFront <- sapply(allInd,function(f){
  idx <- which(spatialMetrics$indUniqID==f)
  as.numeric(tapply(spatialMetrics$inFrontHalf[idx],rep(1:(length(idx)/windowSize),each=windowSize),mean,na.rm=T))
})

if(saveImage){
  if(format=="pdf")pdf(file = "figures/Figure_3.pdf", width = 18, height = 23.5)
  if(format=="svg")svg(file="figures/Figure_3.svg", width = 18, height = 23.5)
  if(format=="png")png(file="figures/Figure_3.png", width = 18, height = 23.5,units="in",res=300)
}

par(mar=c(5,8,4,2),mfrow=c(1,1))

propFrontTable <- data.frame(x=unlist(windowPropFront),y=factor(rep(allIndShort,sapply(windowPropFront,length)),levels=allIndShort))
propFrontTable$status <- factor(allIndInfo$status[match(propFrontTable$y,gsub("20","",allIndInfo$uniqueID))],status)
propFrontTable$session <- factor(rep(sapply(names(windowPropFront),function(f)strsplit(f,"_")[[1]][1]),sapply(windowPropFront,length)),levels=sessions)
propFrontTable <- propFrontTable[-which(is.na(propFrontTable$x)),]

colo = statusCol[match(modelParam_MovSpeed$status,status)]

myRidgeline(propFrontTable$x,propFrontTable$y,palette=colo,grouping=modelParam_MovSpeed$session,spacing=2,axes=T,xlab="Proportion of time in the front half",
            ylab="",meanNotMode = T,mode=T,modeCol="lightblue",yaxt="n",labCex=2,cex.lab=3,cex.axis=2)
box(lwd=10)
abline(v=0.5,lty="dashed",lwd=4)
legend("topright",legend=status,pch=19,col=statusCol,box.lwd=5,cex=2,bg="white")

if(saveImage)dev.off()

mod <- lme(x~status,random=~1|session,data=propFrontTable)
anova(mod)
summary(mod)
qqnorm(residuals(mod))
qqline(residuals(mod))

#----MOVEMENT TURN INFLUENCE SCORES----

if(saveImage){
  if(format=="pdf")pdf(file = "figures/Figure_2a.pdf", width = 20, height = 33.6)
  if(format=="svg")svg(file="figures/Figure_2a.svg", width = 20, height = 33.6)
  if(format=="png")png(file="figures/Figure_2a.png", width = 20, height = 33.6,units="in",res=300)
}

par(mfrow=c(1,1), mar=c(5,6.3,4,2),cex=cex)

spacing <- 1

yAx <- rev(1:nrow(modelParam_MovTurn))*spacing - match(modelParam_MovTurn$session,unique(modelParam_MovTurn$session))*spacing
yAx <- yAx + abs(min(yAx))
xAx <- modelParam_MovTurn$inflScore
pntSize = normalize(modelParam_MovTurn$N,1,3)
colo <- statusCol[match(modelParam_MovTurn$status,status)]

if(plotCIs){
  xMin = min(modelParam_MovTurn$lowerCI) - 0.01
  xMax = max(modelParam_MovTurn$upperCI) + 0.01
}else{
  xMin = min(xAx-0.05)
  xMax = max(xAx+0.1)
}

plot(NULL,xlim=c(xMin,xMax),ylim=c(1,max(yAx)-1),yaxt="n",ylab="",xlab="Turn influence score",
     main="",cex.axis=2,cex.lab=2.3)

segments(0,yAx,xAx,yAx,col="grey",lty="dotted",lwd=6)

if(plotCIs){
  arrows(xAx,yAx,modelParam_MovTurn$lowerCI,yAx,angle=90,length=0.1,col=colo,lwd=5)
  arrows(xAx,yAx,modelParam_MovTurn$upperCI,yAx,angle=90,length=0.1,col=colo,lwd=5)
}

points(xAx,yAx,bg=colo,pch=21,cex=pntSize)

box(lwd=10)
sessionLines = yAx[which(diff(yAx)==min(diff(yAx)))]-1
abline(h=sessionLines,lty="dashed",lwd=5)
axis(2,at=yAx,gsub("20","",modelParam_MovTurn$ind),tick=F,las=1,cex.axis=1.2)
#legend("topright",legend=status,pch=19,col=statusCol,pt.cex=cex,cex=1.3)
abline(v=0.5,lty="longdash",lwd=2)

if(saveImage)dev.off()

#GLMM
modelParam_MovTurn$propFront <- totPropFront
modelParam_MovTurn$status <- factor(modelParam_MovTurn$status,levels=status)

#NOT INCLUDING TIME IN FRONT

mod <- lme(inflScore~status,random=~1|session,data=modelParam_MovTurn)
anova(mod)
summary(mod)
qqnorm(residuals(mod))
qqline(residuals(mod))

#paiwise status comparison
pairWise <- glht(mod, linfct = mcp(status = "Tukey"))
summary(pairWise)

#INCLUDING TIME IN FRONT

mod <- lme(inflScore~(status+propFront),random=~1|session,data=modelParam_MovTurn)
anova(mod)
summary(mod)
qqnorm(residuals(mod))
qqline(residuals(mod))

#paiwise status comparison
pairWise <- glht(mod, linfct = mcp(status = "Tukey"))
summary(pairWise)


#----MOVEMENT SPEED INFLUENCE SCORES----

if(saveImage){
  if(format=="pdf")pdf(file = "figures/Figure_2b.pdf", width = 20, height = 33.6)
  if(format=="svg")svg(file="figures/Figure_2b.svg", width = 20, height = 33.6)
  if(format=="png")png(file="figures/Figure_2b.png", width = 20, height = 33.6,units="in",res=300)
}
par(mfrow=c(1,1), mar=c(5,6.3,4,2),cex=cex)

spacing <- 1

yAx <- rev(1:nrow(modelParam_MovSpeed))*spacing - match(modelParam_MovSpeed$session,unique(modelParam_MovSpeed$session))*spacing
yAx <- yAx + abs(min(yAx))
xAx <- modelParam_MovSpeed$inflScore
pntSize = normalize(modelParam_MovSpeed$N,1,3)
colo <- statusCol[match(modelParam_MovSpeed$status,status)]

if(plotCIs){
  xMin = min(modelParam_MovSpeed$lowerCI) - 0.01
  xMax = max(modelParam_MovSpeed$upperCI) + 0.01
}else{
  xMin = min(xAx-0.05)
  xMax = max(xAx+0.1)
}

plot(NULL,xlim=c(xMin,xMax),ylim=c(1,max(yAx)-1),yaxt="n",ylab="",xlab="Speed influence score",
     main="",cex.axis=2,cex.lab=2.3)

segments(0,yAx,xAx,yAx,col="grey",lty="dotted",lwd=6)

if(plotCIs){
  arrows(xAx,yAx,modelParam_MovSpeed$lowerCI,yAx,angle=90,length=0.1,col=colo,lwd=5)
  arrows(xAx,yAx,modelParam_MovSpeed$upperCI,yAx,angle=90,length=0.1,col=colo,lwd=5)
}

points(xAx,yAx,bg=colo,pch=21,cex=pntSize)

box(lwd=10)
sessionLines <- yAx[which(diff(yAx)==-2)]-1
abline(h=sessionLines,lty="dashed",lwd=5)
axis(2,at=yAx,gsub("20","",modelParam_MovSpeed$ind),tick=F,las=1,cex.axis=1.2)
legend("topright",legend=status,pch=19,col=statusCol,pt.cex=cex,cex=1.3,box.lwd=5)
segments(unique(modelParam_MovSpeed$gamma),c(100,sessionLines),unique(modelParam_MovSpeed$gamma),c(sessionLines,-5),lty="longdash",lwd=2)

if(saveImage)dev.off()

#GLMM
modelParam_MovSpeed$propFront <- totPropFront
modelParam_MovSpeed$status <- factor(modelParam_MovSpeed$status,levels=status)

#NOT INCLUDING TIME IN FRONT

mod <- lme(inflScore~status,random=~1|session,data=modelParam_MovSpeed)
anova(mod)
summary(mod)
qqnorm(residuals(mod))
qqline(residuals(mod))

#paiwise status comparison
pairWise <- glht(mod, linfct = mcp(status = "Tukey"))
summary(pairWise)

#INCLUDING TIME IN FRONT

mod <- lme(inflScore~(status+propFront),random=~1|session,data=modelParam_MovSpeed)
anova(mod)
summary(mod)
qqnorm(residuals(mod))
qqline(residuals(mod))

#paiwise status comparison
pairWise <- glht(mod, linfct = mcp(status = "Tukey"))
summary(pairWise)

#----MOVEMENT TURN VS MOVEMENT SPEED----

if(saveImage){
  if(format=="pdf")pdf(file = "figures/Figure_4a.pdf", width = 20, height = 20)
  if(format=="svg")svg(file="figures/Figure_4a.svg", width = 20, height = 20)
  if(format=="png")png(file="figures/Figure_4a.png", width = 20, height = 20,units="in",res=300)
}

par(mfrow=c(1,1), mar=c(5,5,4,2),cex=cex)

x = modelParam_MovSpeed$inflScore
y = modelParam_MovTurn$inflScore

borderCol <- statusCol
borderCol[c(3,4,5)] <- 1

plot(x,y,bg=groupCol[match(modelParam_MovTurn$session,sessions)],col=borderCol[match(modelParam_MovTurn$status,status)],pch=(21:25)[match(modelParam_MovTurn$session,sessions)],
     xlim=c(min(x)-0.05,max(x)+0.05),ylim=c(min(y)-0.05,max(y)+0.05),lwd=4,
     xlab="Speed influence score",ylab="Turn influence score",
     cex=2,cex.axis=2,cex.lab=2,cex.main=2)
box(lwd=10)

cor_test(data.frame(x,y,factor(modelParam_MovTurn$session,levels=sessions)),"x","y",include_factors = T,multilevel=T)

if(saveImage)dev.off()


#----MOVEMENT TURN INFLUENCE AS A FUNCTION OF TIME SPENT IN THE FRONT----

if(saveImage){
  if(format=="pdf")pdf(file = "figures/Figure_4b.pdf", width = 20, height = 20)
  if(format=="svg")svg(file="figures/Figure_4b.svg", width = 20, height = 20)
  if(format=="png")png(file="figures/Figure_4b.png", width = 20, height = 20,units="in",res=300)
}

par(mfrow=c(1,1), mar=c(5,5,4,2),cex=cex)

x = totPropFront
y = modelParam_MovTurn$inflScore

plot(x,y,bg=groupCol[match(modelParam_MovTurn$session,sessions)],col=borderCol[match(modelParam_MovTurn$status,status)],pch=(21:25)[match(modelParam_MovTurn$session,sessions)],
     xlim=c(min(x)-0.1,max(x)+0.1),ylim=c(min(y)-0.1,max(y)+0.1),lwd=4,
     xlab="Proportion of time in the front",ylab="Turn influence score",
     cex=2,cex.axis=2,cex.lab=2,cex.main=2)
box(lwd=10)

cor_test(data.frame(x,y,factor(modelParam_MovTurn$session,levels=sessions)),"x","y",include_factors = T,multilevel=T)

if(saveImage)dev.off()


#----MOVEMENT SPEED INFLUENCE AS A FUNCTION OF TIME SPENT IN THE FRONT----

if(saveImage){
  if(format=="pdf")pdf(file = "figures/Figure_4c.pdf", width = 20, height = 20)
  if(format=="svg")svg(file="figures/Figure_4c.svg", width = 20, height = 20)
  if(format=="png")png(file="figures/Figure_4c.png", width = 20, height = 20,units="in",res=300)
}
par(mfrow=c(1,1), mar=c(5,5,4,2),cex=cex)

x = totPropFront
y = modelParam_MovSpeed$inflScore

plot(x,y,bg=groupCol[match(modelParam_MovTurn$session,sessions)],col=borderCol[match(modelParam_MovTurn$status,status)],pch=(21:25)[match(modelParam_MovTurn$session,sessions)],
     xlim=c(min(x)-0.1,max(x)+0.1),ylim=c(min(y)-0.1,max(y)+0.1),lwd=4,
     xlab="Proportion of time in the front",ylab="Speed influence score",
     cex=2,cex.axis=2,cex.lab=2,cex.main=2)
box(lwd=10)
legend("topright",legend=sessions,pch=21:25,pt.lwd=4,pt.bg=groupCol,cex=1.5,bg="white",box.lwd=5)

cor_test(data.frame(x,y,factor(modelParam_MovTurn$session,levels=sessions)),"x","y",include_factors = T,multilevel=T)

if(saveImage)dev.off()






