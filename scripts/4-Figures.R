library(fields)
library(correlation)

setwd("C:/Users/baverly/Desktop/INFLUENCE_PAPER")
source("scripts/functions.R")

discretizationStep <- 10

load(paste0("output/spatialMetrics_",discretizationStep,"m.RData"))
load(paste0("output/Influence_logistic_model_fits_",discretizationStep,"m.RData"))

sessions <- c("HM2017","HM2019","L2019","ZU2021","NQ2021")
status = c("DominantF","DominantM","Adult","Yearling","Sub-Adult")

groupCol = hcl.colors(5)
statusCol=rev(c("#000000","#900C3F","#C70039","#FF5733","#FFC300"))
saveImage = T
cex=2.5

#----HEATMAP OF THE PROBABILITY OF THE GROUP TO TURN RIGHT AS A FUNCTION OF LEFT-RIGHT POSITION AND LEFT-RIGHT MOVEMENT----

if(saveImage)pdf(file = "figures/allGroups_turningHeatMap.pdf", width = 11, height = 10)
  
par(mfrow=c(1,1), mar=c(5,5,4,6))

movMax <- ceiling(max(abs(quantile(spatialMetrics$leftRightMovement,c(0.01,0.99),na.rm=T))))
mov1 <- seq(-movMax,movMax,0.1)
posMax <- ceiling(max(abs(quantile(spatialMetrics$leftRightPosition,c(0.01,0.99),na.rm=T))))
pos1 <- seq(-posMax,posMax,0.1)
grid1 <- outer(mov1,pos1,FUN=flatLogis2Variables,alpha=modelParam_Total$alpha[1],beta1=modelParam_Total$beta1[1],beta2=modelParam_Total$beta2[1],gamma=modelParam_Total$gamma[1])
image.plot(x=mov1,y=pos1,z=grid1,col=hcl.colors(100, "viridis", rev = TRUE),
           cex.axis=cex,cex.lab=cex,cex.main=cex,axis.args=(list(cex.axis=2)),legend.width=2,legend.mar=7,
           xlab="Individual left-right movement (m/min)",ylab="Individual left-right positon (m)",zlim=c(min(grid1)-0.05,max(grid1)+0.05),main="P(group turns right)")

if(saveImage)dev.off()

#----HEATMAP OF THE PROBILITY OF THE GROUP TO SPEED UP AS A FUNCTION OF FRONT-BACK POSITION AND FRONT-BACK MOVEMENT----

if(saveImage)pdf(file = "figures/allGroups_speedingHeatMap.pdf", width = 11, height = 10)
  
par(mfrow=c(1,1), mar=c(5,5,4,6))

movMax <- ceiling(max(abs(quantile(spatialMetrics$frontBackMovement,c(0.01,0.99),na.rm=T))))
mov2 <- seq(-movMax,movMax,0.1)
posMax <- ceiling(max(abs(quantile(spatialMetrics$frontBackPosition,c(0.01,0.99),na.rm=T))))
pos2 <- seq(-posMax,posMax,0.1)
grid2 <-  outer(mov2,pos2,FUN=flatLogis2Variables,alpha=modelParam_Total$alpha[2],beta1=modelParam_Total$beta1[2],beta2=modelParam_Total$beta2[2],gamma=modelParam_Total$gamma[2])
image.plot(x=mov2,y=pos2,z=grid2,col=hcl.colors(100, "viridis", rev = TRUE),
           cex.axis=cex,cex.lab=cex,legend.cex=cex,cex.main=cex,axis.args=(list(cex.axis=2)),legend.width=2,legend.mar=7,
           xlab="Individual front-back movement (m/min)",ylab="Individual front-back posiiton (m)",zlim=c(min(grid2)-0.05,max(grid2)+0.05),main="P(group speeds up)")

if(saveImage)dev.off()


#----MOVEMENT TURNING INFLUENCE SCORES----


if(saveImage)pdf(file = "figures/Movement_turning_influence_scores.pdf", width = 20, height = 35)

par(mfrow=c(1,1), mar=c(5,6.3,4,2),cex=cex)

spacing <- 1

yAx <- rev(1:nrow(modelParam_MovTurn))*spacing - match(modelParam_MovTurn$session,unique(modelParam_MovTurn$session))*spacing
yAx <- yAx + abs(min(yAx))
pntSize = normalize(modelParam_MovTurn$N,1,3)

plot(NULL,xlim=c(min(modelParam_MovTurn$inflScore-0.05),max(modelParam_MovTurn$inflScore+0.1)),ylim=c(1,max(yAx)-1),yaxt="n",ylab="",xlab="Turning influence score",
     main="",cex.axis=2,cex.lab=2.3)
for(i in 1:nrow(modelParam_MovTurn)){
  
  ind = modelParam_MovTurn$ind[i]
  colo = statusCol[match(modelParam_MovTurn$status[i],status)]
  
  xAx <- modelParam_MovTurn$inflScore[i]
  
  segments(0,yAx[i],xAx,yAx[i],col="grey",lty="dotted",lwd=6)
  #abline(h=yAx[i],col="grey",lty="dotted",lwd=4)
  points(xAx,yAx[i],bg=colo,pch=21,cex=pntSize[i])
  
}
abline(h=yAx[which(diff(yAx)==min(diff(yAx)))]-1,lty="dashed",lwd=4)
axis(2,at=yAx,gsub("20","",modelParam_MovTurn$ind),tick=F,las=1,cex.axis=1.5)
#legend("topright",legend=status,pch=19,col=statusCol,pt.cex=cex,cex=1.3)
abline(v=0.5,lty="longdash",lwd=2)

if(saveImage)dev.off()


#----MOVEMENT SPEEDING INFLUENCE SCORES----

if(saveImage)pdf(file = "figures/Movement_speeding_influence_scores.pdf", width = 20, height = 35)

par(mfrow=c(1,1), mar=c(5,6.3,4,2),cex=cex)

spacing <- 1

yAx <- rev(1:nrow(modelParam_MovSpeed))*spacing - match(modelParam_MovSpeed$session,unique(modelParam_MovSpeed$session))*spacing
yAx <- yAx + abs(min(yAx))
pntSize = normalize(modelParam_MovSpeed$N,1,3)

plot(NULL,xlim=c(min(modelParam_MovSpeed$inflScore-0.1),max(modelParam_MovSpeed$inflScore+0.05)),ylim=c(1,max(yAx)-1),yaxt="n",ylab="",xlab="Speeding influence score",
     main="",cex.axis=2,cex.lab=2.3)
for(i in 1:nrow(modelParam_MovSpeed)){
  
  ind = modelParam_MovSpeed$ind[i]
  colo = statusCol[match(modelParam_MovSpeed$status[i],status)]
  
  xAx <- modelParam_MovSpeed$inflScore[i]
  
  segments(0,yAx[i],xAx,yAx[i],col="grey",lty="dotted",lwd=6)
  #abline(h=yAx[i],col="grey",lty="dotted",lwd=4)
  points(xAx,yAx[i],bg=colo,pch=21,cex=pntSize[i])
  
}
sessionLines = yAx[which(diff(yAx)==-2)]-1
abline(h=sessionLines,lty="dashed",lwd=4)
axis(2,at=yAx,gsub("20","",modelParam_MovTurn$ind),tick=F,las=1,cex.axis=1.5)
legend("topright",legend=status,pch=19,col=statusCol,pt.cex=cex,cex=1.3)
segments(unique(modelParam_MovSpeed$gamma),c(100,sessionLines),unique(modelParam_MovSpeed$gamma),c(sessionLines,-5),lty="longdash",lwd=2)

if(saveImage)dev.off()


#----REGRESSION MOVEMENT TURNING VS MOVEMENT SPEEDING ACCROSS GROUPS

if(saveImage)pdf(file = "figures/Movement_turning_influence_VS_movement_speeding_influence_ACROSS_GROUPS.pdf", width = 20, height = 20)

par(mfrow=c(1,1), mar=c(5,5,4,2),cex=cex)

x = modelParam_MovSpeed$inflScore
y = modelParam_MovTurn$inflScore

plot(x,y,bg=statusCol[match(modelParam_MovTurn$status,status)],pch=(21:25)[match(modelParam_MovTurn$session,sessions)],
     xlim=c(min(x)-0.05,max(x)+0.05),ylim=c(min(y)-0.05,max(y)+0.05),
     xlab="Speeding influence score",ylab="Turning influence score",
     cex=2,cex.axis=2,cex.lab=2,cex.main=2)
#legend("topright",legend=status,fill=statusCol,bty="n")
legend("topright",legend=status,pch=19,col=statusCol,cex=1.5)
legend("bottomright",legend=sessions,pch=21:25,cex=1.5)

regLine <- lm(y~x)
abline(regLine,lwd=4,lty="dashed")
summary(regLine)

if(saveImage)dev.off()

#----REGRESSION MOVEMENT TURNING VS MOVEMENT SPEEDING PER GROUPS

if(saveImage)pdf(file = "figures/Movement_turning_influence_VS_movement_speeding_influence_PER_GROUP.pdf", width = 20, height = 20)

par(mfrow=c(1,1), mar=c(5,5,4,2),cex=cex)

x = modelParam_MovSpeed$inflScore
y = modelParam_MovTurn$inflScore

plot(x,y,bg=groupCol[match(modelParam_MovTurn$session,sessions)],pch=(21:25)[match(modelParam_MovTurn$session,sessions)],
     xlim=c(min(x)-0.05,max(x)+0.05),ylim=c(min(y)-0.05,max(y)+0.05),
     xlab="Speeding influence score",ylab="Turning influence score",
     cex=2,cex.axis=2,cex.lab=2,cex.main=2)

cor_test(data.frame(x,y,factor(modelParam_MovTurn$session,levels=sessions)),"x","y",include_factors = T,multilevel=T)

#regTurnSpeed <- data.frame("group"=sessions,"f"=NA,"DF"=NA,"pValue"=NA)
#for(session in sessions){
#  
#  idx <- which(modelParam_MovTurn$session==session)
#  
#  xSession <- x[idx]
#  ySession <- y[idx]
#  
#  regLine <- lm(ySession~xSession)
#  
#  minX <- min(x[idx]-0.05)
#  maxX <- max(x[idx]+0.05)
#  #segments(minX,predict(regLine,data.frame(xSession=minX)),maxX,predict(regLine,data.frame(xSession=maxX)),lwd=6,lty="dashed",col=groupCol[match(session,sessions)])
#  
#  #abline(regLine,lwd=4,lty="dashed",col=groupCol[match(session,sessions)])
#  summa <- anova(regLine)
#  regTurnSpeed$f[match(session,sessions)] <- summa$`F value`[1]
#  regTurnSpeed$DF[match(session,sessions)] <- summa$Df[2]
#  regTurnSpeed$pValue[match(session,sessions)] <- summa$`Pr(>F)`[1]
#}
#legend("topright",legend=sessions,pch=21:25,pt.bg=groupCol,cex=1.5)

if(saveImage)dev.off()

#----DISTRIBUTION OF TIME SPENT IN THE FRONT----

allInd = modelParam_MovTurn$ind
allIndShort <- gsub("20","",allInd)

totPropFront <- tapply(spatialMetrics$inFrontHalf,factor(spatialMetrics$indUniqID,levels=allInd),mean,na.rm=T)

windowSize = 60*60
windowPropFront <- sapply(allInd,function(f){
  idx <- which(spatialMetrics$indUniqID==f)
  as.numeric(tapply(spatialMetrics$inFrontHalf[idx],rep(1:(length(idx)/windowSize),each=windowSize),mean,na.rm=T))
})


if(saveImage)pdf(file = "figures/PropFront_ridgeLine.pdf", width = 15, height = 20)

par(mar=c(5,8,4,2),mfrow=c(1,1))

x <- unlist(windowPropFront)
y <- factor(rep(allIndShort,sapply(windowPropFront,length)),levels=allIndShort)
y<- y[which(!is.na(x))]
x<- x[which(!is.na(x))]

colo = statusCol[match(modelParam_MovSpeed$status,status)]

myRidgeline(x,y,palette=colo,grouping=modelParam_MovSpeed$session,spacing=2,axes=T,xlab="Proportion of time in the front half",
            ylab="",meanNotMode = T,mode=T,modeCol="lightblue",yaxt="n",labCex=2,cex.lab=2.5,cex.axis=2.5)

abline(v=0.5,lty="dashed",lwd=4)
legend("topright",legend=status,pch=19,col=statusCol,cex=2,bg="white")

if(saveImage)dev.off()


#----MOVEMENT TURNING INFLUENCE AS A FUNCTION OF TIME SPENT IN THE FRONT ACROSS GROUPS


if(saveImage)pdf(file = "figures/Movement_turning_influence_VS_position.pdf", width = 20, height = 20)

par(mfrow=c(1,1), mar=c(5,5,4,2),cex=cex)

x = totPropFront
y = modelParam_MovTurn$inflScore

plot(x,y,bg=statusCol[match(modelParam_MovTurn$status,status)],pch=(21:25)[match(modelParam_MovTurn$session,sessions)],
     xlim=c(min(x)-0.05,max(x)+0.05),ylim=c(min(y)-0.05,max(y)+0.05),xlab="Proportion of time in the front",ylab="Turning influence Score",
     cex=2,cex.axis=2,cex.lab=2,cex.main=2)
#legend("topright",legend=status,pch=19,col=statusCol,cex=1.5)
#legend("bottomright",legend=sessions,pch=21:25,cex=1.5)

regLine <- lm(y~x)
abline(regLine,lwd=4,lty="dashed")
summary(regLine)

if(saveImage)dev.off()

#----MOVEMENT TURNING INFLUENCE AS A FUNCTION OF TIME SPENT IN THE FRONT PER GROUP

if(saveImage)pdf(file = "figures/Movement_turning_influence_VS_position_PER_GROUP.pdf", width = 20, height = 20)

par(mfrow=c(1,1), mar=c(5,5,4,2),cex=cex)

x = totPropFront
y = modelParam_MovTurn$inflScore

plot(x,y,bg=groupCol[match(modelParam_MovTurn$session,sessions)],pch=(21:25)[match(modelParam_MovTurn$session,sessions)],
     xlim=c(min(x)-0.1,max(x)+0.1),ylim=c(min(y)-0.1,max(y)+0.1),
     xlab="Proportion of time in the front",ylab="Turning influence score",
     cex=2,cex.axis=2,cex.lab=2,cex.main=2)

cor_test(data.frame(x,y,factor(modelParam_MovTurn$session,levels=sessions)),"x","y",include_factors = T,multilevel=T)

#regTurnFront <- data.frame("group"=sessions,"f"=NA,"DF"=NA,"pValue"=NA)
#for(session in sessions){
#  
#  idx <- which(modelParam_MovTurn$session==session)
#  
#  xSession <- x[idx]
#  ySession <- y[idx]
#  
#  regLine <- lm(ySession~xSession)
#  
#  minX <- min(x[idx]-0.05)
#  maxX <- max(x[idx]+0.05)
#  #segments(minX,predict(regLine,data.frame(xSession=minX)),maxX,predict(regLine,data.frame(xSession=maxX)),lwd=6,lty="dashed",col=groupCol[match(session,sessions)])
#  
#  #abline(regLine,lwd=4,lty="dashed",col=groupCol[match(session,sessions)])
#  
#  summa <- anova(regLine)
#  regTurnFront$f[match(session,sessions)] <- summa$`F value`[1]
#  regTurnFront$DF[match(session,sessions)] <- summa$Df[2]
#  regTurnFront$pValue[match(session,sessions)] <- summa$`Pr(>F)`[1]
#}
#legend("topright",legend=sessions,pch=21:25,pt.bg=groupCol,cex=1.5)

if(saveImage)dev.off()

#----MOVEMENT SPEEDING INFLUENCE AS A FUNCTION OF TIME SPENT IN THE FRONT ACROSS GROUPS

if(saveImage)pdf(file = "figures/Movement_speeding_influence_VS_position.pdf", width = 20, height = 20)

par(mfrow=c(1,1), mar=c(5,5,4,2),cex=cex)

x = totPropFront
y = modelParam_MovSpeed$inflScore

plot(x,y,bg=statusCol[match(modelParam_MovSpeed$status,status)],pch=(21:25)[match(modelParam_MovSpeed$session,sessions)],
     xlim=c(min(x)-0.1,max(x)+0.1),ylim=c(min(y)-0.1,max(y)+0.1),xlab="Proportion of time in the front",ylab="Turning influence Score",
     cex=2,cex.axis=2,cex.lab=2,cex.main=2)
#legend("topright",legend=status,fill=statusCol,bty="n")
legend("topright",legend=status,pch=19,col=statusCol,cex=1.4)
legend("bottomright",legend=sessions,pch=21:25,cex=1.4)

regLine <- lm(y~x)
abline(regLine,lwd=4,lty="dashed")
summary(regLine)

if(saveImage)dev.off()


#----MOVEMENT SPEEDING INFLUENCE AS A FUNCTION OF TIME SPENT IN THE FRONT PER GROUPS

if(saveImage)pdf(file = "figures/Movement_speeding_influence_VS_position_PER_GROUP.pdf", width = 20, height = 20)

par(mfrow=c(1,1), mar=c(5,5,4,2),cex=cex)

x = totPropFront
y = modelParam_MovSpeed$inflScore

plot(x,y,bg=groupCol[match(modelParam_MovTurn$session,sessions)],pch=(21:25)[match(modelParam_MovTurn$session,sessions)],
     xlim=c(min(x)-0.1,max(x)+0.1),ylim=c(min(y)-0.1,max(y)+0.1),
     xlab="Proportion of time in the front",ylab="Speeding influence score",
     cex=2,cex.axis=2,cex.lab=2,cex.main=2)

cor_test(data.frame(x,y,factor(modelParam_MovTurn$session,levels=sessions)),"x","y",include_factors = T,multilevel=T)


#regSpeedFront <- data.frame("group"=sessions,"f"=NA,"DF"=NA,"pValue"=NA)
#for(session in sessions){
#  
#  idx <- which(modelParam_MovTurn$session==session)
#  
#  xSession <- x[idx]
#  ySession <- y[idx]
#  
#  regLine <- lm(ySession~xSession)
#  
#  minX <- min(x[idx]-0.05)
#  maxX <- max(x[idx]+0.05)
#  #segments(minX,predict(regLine,data.frame(xSession=minX)),maxX,predict(regLine,data.frame(xSession=maxX)),lwd=6,lty="dashed",col=groupCol[match(session,sessions)])
#  
#  #abline(regLine,lwd=4,lty="dashed",col=groupCol[match(session,sessions)])
#  
#  #summa <- anova(regLine)
#  #regSpeedFront$f[match(session,sessions)] <- summa$`F value`[1]
#  #regSpeedFront$DF[match(session,sessions)] <- summa$Df[2]
#  #regSpeedFront$pValue[match(session,sessions)] <- summa$`Pr(>F)`[1]
#  
#}
#legend("topright",legend=sessions,pch=21:25,pt.bg=groupCol,cex=1.5)

if(saveImage)dev.off()






