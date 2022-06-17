#THE PURPOSE OF THIS SCRIPT IS TO OUTPUT THE FIGURES FOUND IN THE SUPPLEMENTARY MATERIAL

library(fields)

setwd("C:/Users/baverly/Desktop/INFLUENCE_PAPER")
source("scripts/functions.R")

discretizationStep <- 10

load(paste0("output/spatialMetrics_",discretizationStep,"m.RData"))
load(paste0("output/Influence_logistic_model_fits_",discretizationStep,"m.RData"))

sessions <- c("HM2017","HM2019","L2019","ZU2021","NQ2021")
status = c("DominantF","DominantM","Adult","Yearling","Sub-Adult")

statusCol=rev(c("#000000","#900C3F","#C70039","#FF5733","#FFC300"))
saveImage = T
cex=2.5
lwd=4

allInd <- modelParam_MovTurn$ind

dir.create(path="figures/SuppMat/",showWarnings = F)

#----POSITIONAL TURN INFLUENCE SCORES----

dir.create(path="figures/SuppMat/positional_turn_and_speed_influence/",showWarnings = F)
if(saveImage)pdf(file = "figures/SuppMat/positional_turn_and_speed_influence/Positional_turn_influence_scores.pdf", width = 20, height = 35)

par(mfrow=c(1,1), mar=c(5,6.3,4,2),cex=cex)

spacing <- 1

yAx <- rev(1:nrow(modelParam_PosTurn))*spacing - match(modelParam_PosTurn$session,unique(modelParam_PosTurn$session))*spacing
yAx <- yAx + abs(min(yAx))
xAx <- modelParam_PosTurn$inflScore
pntSize = normalize(modelParam_PosTurn$N,1,3)
colo <- statusCol[match(modelParam_MovTurn$status,status)]

plot(NULL,xlim=c(min(xAx-0.02),max(xAx+0.02)),ylim=c(1,max(yAx)-1),yaxt="n",ylab="",xlab="Positional turn influence score",
     main="",cex.axis=2,cex.lab=2.3)

segments(0,yAx,xAx,yAx,col="grey",lty="dotted",lwd=6)

points(xAx,yAx,bg=colo,pch=21,cex=pntSize)

box(lwd=10)
sessionLines = yAx[which(diff(yAx)==min(diff(yAx)))]-1
abline(h=sessionLines,lty="dashed",lwd=4)
axis(2,at=yAx,gsub("20","",modelParam_PosTurn$ind),tick=F,las=1,cex.axis=1.5)
#legend("topright",legend=status,pch=19,col=statusCol,pt.cex=cex,cex=1.3)
abline(v=0.5,lty="longdash",lwd=2)

if(saveImage)dev.off()


#----POSITIONAL SPEED INFLUENCE SCORES----

if(saveImage)pdf(file = "figures/SuppMat/positional_turn_and_speed_influence/Positional_speed_influence_scores.pdf", width = 20, height = 35)

par(mfrow=c(1,1), mar=c(5,6.3,4,2),cex=cex)

yAx <- rev(1:nrow(modelParam_PosSpeed)) - match(modelParam_PosSpeed$session,unique(modelParam_PosSpeed$session)) + 4
xAx <- modelParam_PosSpeed$inflScore
pntSize = normalize(modelParam_PosSpeed$N,1,3)
colo <- statusCol[match(modelParam_PosSpeed$status,status)]

plot(NULL,xlim=c(min(modelParam_PosSpeed$gamma)-0.01,max(xAx)+0.02),ylim=c(1,max(yAx)-1),yaxt="n",ylab="",xlab="Positional speed influence score",
     main="",cex.axis=2,cex.lab=2.3)

segments(0,yAx,xAx,yAx,col="grey",lty="dotted",lwd=6)

points(xAx,yAx,bg=colo,pch=21,cex=pntSize)

box(lwd=10)
sessionLines <- yAx[which(diff(yAx)==-2)]-1
abline(h=sessionLines,lty="dashed",lwd=4)
axis(2,at=yAx,gsub("20","",modelParam_PosSpeed$ind),tick=F,las=1,cex.axis=1.5)
legend("topright",legend=status,pch=19,col=statusCol,pt.cex=cex,cex=1.3)
segments(unique(modelParam_PosSpeed$gamma),c(100,sessionLines),unique(modelParam_PosSpeed$gamma),c(sessionLines,-5),lty="longdash",lwd=2)

if(saveImage)dev.off()

#----HEATMAP OF THE PROBABILITY OF THE GROUP TO TURN RIGHT AS A FUNCTION OF LEFT-RIGHT POSITION AND LEFT-RIGHT MOVEMENT----

dir.create(path="figures/SuppMat/two_predictors_models/",showWarnings = F)
if(saveImage)pdf(file = "figures/SuppMat/two_predictors_models/allGroups_turningHeatMap.pdf", width = 11, height = 10)

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

if(saveImage)pdf(file = "figures/SuppMat/two_predictors_models/allGroups_speedingHeatMap.pdf", width = 11, height = 10)

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

#----DISTRIBUTION OF GROUP SPREAD----

nearNeigh <- list()
allDistances <- list()

for(session in sessions){
  
  dir.create(path="figures/SuppMat/distributions/",showWarnings = F)
  if(saveImage)png(filename = paste0("figures/SuppMat/distributions/groupSpread_",session,".png"), width = 1000, height = 1000)
  
  par(mfrow=c(1,1), mar=c(5,5,4,2))
  
  longNames <- c(load(paste0("data/movement/level1/",session,"_COORDINATES_SYNCHED_level1.RData")))
  shortNames <- simplifyNames(pat=paste(session,"_",sep=""))
  
  allDistances[[session]] <- mapply(allDist,as.data.frame(allX),as.data.frame(allY),symetry=T,SIMPLIFY="array")
  
  nearNeigh[[session]] <- apply(allDistances[[session]],c(1,3),min,na.rm=T)
  nearNeigh[[session]][which(is.infinite(nearNeigh[[session]]))] <- NA
  
  spread <- apply(allDistances[[session]],3,mean,na.rm=T)
  
  hist(spread,breaks=50,main=gsub("20","",session),xlab="Group spread (m)",xlim=c(0,70),cex.main=3,cex.lab=2.5,cex.axis=2.5)
  
  if(saveImage)dev.off()
}


#----DISTRIBUTION OF TIME SCALES FOR THE GIVEN DISCRETIZATION STEP----

for(session in sessions){
  
  if(saveImage)png(filename = paste0("figures/SuppMat/distributions/discretizationTimeSteps_",session,".png"), width = 1000, height = 1000)
  
  par(mfrow=c(1,1), mar=c(5,5,4,2))
  
  sessionTable <- spatialMetrics[which(spatialMetrics$session==session),]
  sessionTable <- sessionTable[match(unique(sessionTable$t),sessionTable$t),]
  
  hist(sessionTable$futurStepDuration,breaks=50,xlim=c(0,1500),
       main=gsub("20","",session),xlab="Duration of time steps for a spatial discretization of 10 m (s)",cex.main=3,cex.lab=2.5,cex.axis=2.5)
  
  if(saveImage)dev.off()
}

#----DISTRIBUTION OF GROUP TURNING ANGLES----
for(session in sessions){
  
  if(saveImage)png(filename = paste0("figures/SuppMat/distributions/turningAngle_",session,".png"), width = 1000, height = 1000)
  
  par(mfrow=c(1,1), mar=c(5,5,4,2))
  
  sessionTable <- spatialMetrics[which(spatialMetrics$session==session),]
  sessionTable <- sessionTable[match(unique(sessionTable$t),sessionTable$t),]
  
  metric <- sessionTable$groupTurnAngle
  
  hist(metric,breaks=50,
       main=gsub("20","",session),xlab="Turning angle of group centroid (rad)",cex.main=3,cex.lab=2.5,cex.axis=2.5)
  
  if(saveImage)dev.off()
}

#----INFLUENCE SCORES WHEN VARYING DISCRETIZATION STEP----

for(discretizationStep in c(5,10,15,20)){
  
  load(paste0("output/Influence_logistic_model_fits_",discretizationStep,"m.RData"))
  load(paste0("output/spatialMetrics_",discretizationStep,"m.RData"))
  
  #movement turn influence scores
  
  dir.create(path="figures/SuppMat/varying_discretization_step_length/",showWarnings = F)
  if(saveImage)pdf(file = paste0("figures/SuppMat/varying_discretization_step_length/Movement_turn_influence_scores_",discretizationStep,"m.pdf"),width = 20, height = 35)
  
  par(mfrow=c(1,1), mar=c(5,6.3,4,2),cex=cex)
  
  spacing <- 1
  
  yAx <- rev(1:nrow(modelParam_MovTurn))*spacing - match(modelParam_MovTurn$session,unique(modelParam_MovTurn$session))*spacing
  yAx <- yAx + abs(min(yAx))
  xAx <- modelParam_MovTurn$inflScore
  pntSize = normalize(modelParam_MovTurn$N,1,3)
  colo <- statusCol[match(modelParam_MovTurn$status,status)]
  
  plot(NULL,xlim=c(min(xAx-0.02),max(xAx+0.02)),ylim=c(1,max(yAx)-1),yaxt="n",ylab="",xlab="Turn influence score",
       main=paste0("Discretization step = ",discretizationStep,"m"),cex.main=2.5,cex.axis=2,cex.lab=2.3)
  
  segments(0,yAx,xAx,yAx,col="grey",lty="dotted",lwd=6)
  
  points(xAx,yAx,bg=colo,pch=21,cex=pntSize)
  
  box(lwd=10)
  sessionLines = yAx[which(diff(yAx)==min(diff(yAx)))]-1
  abline(h=sessionLines,lty="dashed",lwd=4)
  axis(2,at=yAx,gsub("20","",modelParam_MovTurn$ind),tick=F,las=1,cex.axis=1.5)
  #legend("topright",legend=status,pch=19,col=statusCol,pt.cex=cex,cex=1.3)
  abline(v=0.5,lty="longdash",lwd=2) 
  if(saveImage)dev.off()
  
  
  #movement speed influence scores
  
  if(saveImage)pdf(file = paste0("figures/SuppMat/varying_discretization_step_length/Movement_speed_influence_scores_",discretizationStep,"m.pdf"),width = 20, height = 35)
  
  par(mfrow=c(1,1), mar=c(5,6.3,4,2),cex=cex)
  
  spacing <- 1
  
  yAx <- rev(1:nrow(modelParam_MovSpeed))*spacing - match(modelParam_MovSpeed$session,unique(modelParam_MovSpeed$session))*spacing
  yAx <- yAx + abs(min(yAx))
  xAx <- modelParam_MovSpeed$inflScore
  pntSize = normalize(modelParam_MovSpeed$N,1,3)
  colo <- statusCol[match(modelParam_MovSpeed$status,status)]
  
  plot(NULL,xlim=c(min(modelParam_MovSpeed$gamma-0.01),max(xAx+0.02)),ylim=c(1,max(yAx)-1),yaxt="n",ylab="",xlab="Speed influence score",
       main=paste0("Discretization step = ",discretizationStep,"m"),cex.main=2.5,cex.axis=2,cex.lab=2.3)
  
  segments(0,yAx,xAx,yAx,col="grey",lty="dotted",lwd=6)
  
  points(xAx,yAx,bg=colo,pch=21,cex=pntSize)
  
  box(lwd=10)
  sessionLines= yAx[which(diff(yAx)==-2)]-1
  abline(h=sessionLines,lty="dashed",lwd=4)
  axis(2,at=yAx,gsub("20","",modelParam_MovSpeed$ind),tick=F,las=1,cex.axis=1.5)
  legend("topright",legend=status,pch=19,col=statusCol,pt.cex=cex,cex=1.3,box.lwd=5)
  segments(unique(modelParam_MovSpeed$gamma),c(100,sessionLines),unique(modelParam_MovSpeed$gamma),c(sessionLines,-5),lty="longdash",lwd=2)
  
  if(saveImage)dev.off()
  
  #distribution of time spent in the front
  allInd = modelParam_MovTurn$ind
  allIndShort <- gsub("20","",allInd)
  
  totPropFront <- tapply(spatialMetrics$inFrontHalf,factor(spatialMetrics$indUniqID,levels=allInd),mean,na.rm=T)
  
  windowSize = 60*60
  windowPropFront <- sapply(allInd,function(f){
    idx <- which(spatialMetrics$indUniqID==f)
    as.numeric(tapply(spatialMetrics$inFrontHalf[idx],rep(1:(length(idx)/windowSize),each=windowSize),mean,na.rm=T))
  })
  
  if(saveImage)pdf(file = paste0("figures/SuppMat/varying_discretization_step_length/PropFront_ridgeLine_",discretizationStep,"m.pdf"), width = 15, height = 20)
  
  par(mar=c(5,8,4,2),mfrow=c(1,1))
  
  x <- unlist(windowPropFront)
  y <- factor(rep(allIndShort,sapply(windowPropFront,length)),levels=allIndShort)
  y<- y[which(!is.na(x))]
  x<- x[which(!is.na(x))]
  
  colo = statusCol[match(modelParam_MovSpeed$status,status)]
  
  myRidgeline(x,y,palette=colo,grouping=modelParam_MovSpeed$session,spacing=2,axes=T,xlab="Proportion of time in the front half",
              ylab="",meanNotMode = T,mode=T,modeCol="lightblue",yaxt="n",labCex=2,cex.lab=2.3,cex.axis=2,cex.main=3,main=paste0("Discretization step = ",discretizationStep,"m"))
  
  box(lwd=10)
  abline(v=0.5,lty="dashed",lwd=4)
  if(discretizationStep==20)legend("topright",legend=status,pch=19,col=statusCol,box.lwd=5,cex=2,bg="white")
  
  if(saveImage)dev.off()
  
}

#----INDIVIDUAL LOGISTIC FITS FOR EACH 4 INFLUENCE TYPES----

discretizationStep <- 10

load(paste0("output/spatialMetrics_",discretizationStep,"m.RData"))
load(paste0("output/Influence_logistic_model_fits_",discretizationStep,"m.RData"))

for(ind in allInd){
  
  indIdx <- which(allInd==ind)
  session <- allIndInfo$session[which(allIndInfo$uniqueID==ind)]
  indData <-  spatialMetrics[which(spatialMetrics$indUniqID==ind),]
  
  dir.create(path="figures/SuppMat/individual_logistic_fits/",showWarnings = F)
  if(saveImage)png(filename = paste0("figures/SuppMat/individual_logistic_fits/",ind,"_logisticFit.png"), width = 1000, height = 1000)
  par(mfrow=c(2,2), mar=c(5,5,4,2))
  
  #position turning influence
  data <- data.frame(x=indData$leftRightPosition,y=as.numeric(indData$groupTurnsRight))
  data <- data[which(!is.na(data$x) & !is.na(data$y)),]

  gamma <- 0.5
  
  maxRange <- ceiling(max(abs(quantile(data$x,c(0.01,0.99)))))
  
  influencePlots(x=data$x,y=data$y,binSize=round(maxRange/6),totRange=c(-maxRange,maxRange),cex.main=3,cex.lab=2.5,lwd=lwd,
                 xlab="Individual left-right position (m)",ylab="P(group turns right)",main="Positional Turning Influence")
  segments(modelParam_PosTurn$quantile90[indIdx],-10,modelParam_PosTurn$quantile90[indIdx],modelParam_PosTurn$inflScore[indIdx],lty="dashed",col=4,lwd=3)
  segments(modelParam_PosTurn$quantile90[indIdx],modelParam_PosTurn$inflScore[indIdx],-100,modelParam_PosTurn$inflScore[indIdx],lty="dashed",col=4,lwd=3)
  lines(-100:100,flatLogis1Variable(-100:100,modelParam_PosTurn$alpha[indIdx],modelParam_PosTurn$beta[indIdx],gamma=gamma),col=2,lwd=lwd)
  legend("topright",legend=gsub("20","",ind),bty="n",cex=3)
  legend("topleft",cex=2,legend=c(paste0("alpha = ",round(modelParam_PosTurn$alpha[indIdx],3)),paste0("beta = ",round(modelParam_PosTurn$beta[indIdx],3)),
                                    paste0("Likelihood = ",round(modelParam_PosTurn$likelihood[indIdx],0)),paste("Influence score = ",round(modelParam_PosTurn$inflScore[indIdx],2))),bty="n")
  abline(v=0,h=gamma,lty="dashed")
  
  #movement turning influence
  data <- data.frame(x=indData$leftRightMovement,y=as.numeric(indData$groupTurnsRight))
  data <- data[which(!is.na(data$x) & !is.na(data$y)),]
  
  gamma <- 0.5
  
  maxRange <- ceiling(max(abs(quantile(data$x,c(0.01,0.99)))))
  
  influencePlots(x=data$x,y=data$y,binSize=round(maxRange/6),totRange=c(-maxRange,maxRange),cex.main=3,cex.lab=2.5,lwd=lwd,
                 xlab="Individual left-right movement (m/min)",ylab="P(group turns right)",main="Movement Turn Influence")
  segments(modelParam_MovTurn$quantile90[indIdx],-10,modelParam_MovTurn$quantile90[indIdx],modelParam_MovTurn$inflScore[indIdx],lty="dashed",col=4,lwd=3)
  segments(modelParam_MovTurn$quantile90[indIdx],modelParam_MovTurn$inflScore[indIdx],-100,modelParam_MovTurn$inflScore[indIdx],lty="dashed",col=4,lwd=3)
  lines(-100:100,flatLogis1Variable(-100:100,modelParam_MovTurn$alpha[indIdx],modelParam_MovTurn$beta[indIdx],gamma=gamma),col=2,lwd=lwd)
  legend("topright",legend=gsub("20","",ind),bty="n",cex=3)
  legend("topleft",cex=2,legend=c(paste0("alpha = ",round(modelParam_MovTurn$alpha[indIdx],3)),paste0("beta = ",round(modelParam_MovTurn$beta[indIdx],3)),
                                    paste0("Likelihood = ",round(modelParam_MovTurn$likelihood[indIdx],0)),paste("Influence score = ",round(modelParam_MovTurn$inflScore[indIdx],2))),bty="n")
  abline(v=0,h=gamma,lty="dashed")
  
  #position speeding influence
  data <- data.frame(x=indData$frontBackPosition,y=as.numeric(indData$groupSpeedsUp))
  data <- data[which(!is.na(data$x) & !is.na(data$y)),]
  
  gamma <- modelParam_PosSpeed$gamma[indIdx]
  
  maxRange <- ceiling(max(abs(quantile(data$x,c(0.01,0.99)))))
  
  influencePlots(x=data$x,y=data$y,binSize=round(maxRange/6),totRange=c(-maxRange,maxRange),cex.main=3,cex.lab=2.5,lwd=lwd,
                 xlab="Individual front-back position (m)",ylab="P(group speeds up)",main="Positional speed Influence")
  segments(modelParam_PosSpeed$quantile90[indIdx],-10,modelParam_PosSpeed$quantile90[indIdx],modelParam_PosSpeed$inflScore[indIdx],lty="dashed",col=4,lwd=3)
  segments(modelParam_PosSpeed$quantile90[indIdx],modelParam_PosSpeed$inflScore[indIdx],-100,modelParam_PosSpeed$inflScore[indIdx],lty="dashed",col=4,lwd=3)
  lines(-100:100,flatLogis1Variable(-100:100,modelParam_PosSpeed$alpha[indIdx],modelParam_PosSpeed$beta[indIdx],gamma=gamma),col=2,lwd=lwd)
  legend("topright",gsub("20","",ind),bty="n",cex=3)
  legend("topleft",cex=2,legend=c(paste0("alpha = ",round(modelParam_PosSpeed$alpha[indIdx],3)),paste0("beta = ",round(modelParam_PosSpeed$beta[indIdx],3)),
                                    paste0("Likelihood = ",round(modelParam_PosSpeed$likelihood[indIdx],0)),paste("Influence score = ",round(modelParam_PosSpeed$inflScore[indIdx],2))),bty="n")
  abline(v=0,h=gamma,lty="dashed")
  
  #movement speeding influence
  data <- data.frame(x=indData$frontBackMovement,y=as.numeric(indData$groupSpeedsUp))
  data <- data[which(!is.na(data$x) & !is.na(data$y)),]
  
  gamma <- modelParam_MovSpeed$gamma[indIdx]

  maxRange <- ceiling(max(abs(quantile(data$x,c(0.01,0.99)))))
  
  influencePlots(x=data$x,y=data$y,binSize=round(maxRange/6),totRange=c(-maxRange,maxRange),cex.main=3,cex.lab=2.5,lwd=lwd,
                 xlab="Individual front-back movement (m/min)",ylab="P(group speeds up)",main="Movement speed Influence")
  segments(modelParam_MovSpeed$quantile90[indIdx],-10,modelParam_MovSpeed$quantile90[indIdx],modelParam_MovSpeed$inflScore[indIdx],lty="dashed",col=4,lwd=3)
  segments(modelParam_MovSpeed$quantile90[indIdx],modelParam_MovSpeed$inflScore[indIdx],-100,modelParam_MovSpeed$inflScore[indIdx],lty="dashed",col=4,lwd=3)
  lines(-100:100,flatLogis1Variable(-100:100,modelParam_MovSpeed$alpha[indIdx],modelParam_MovSpeed$beta[indIdx],gamma=gamma),col=2,lwd=lwd)
  legend("topright",legend=gsub("20","",ind),bty="n",cex=3)
  legend("topleft",cex=2,legend=c(paste0("alpha = ",round(modelParam_MovSpeed$alpha[indIdx],3)),paste0("beta = ",round(modelParam_MovSpeed$beta[indIdx],3)),
                                    paste0("Likelihood = ",round(modelParam_MovSpeed$likelihood[indIdx],0)),paste("Influence score = ",round(modelParam_MovSpeed$inflScore[indIdx],2))),bty="n")
  abline(v=0,h=gamma,lty="dashed")
  
  if(saveImage)dev.off()
}








