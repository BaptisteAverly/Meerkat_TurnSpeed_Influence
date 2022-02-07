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

#----individual logistic fits for each 4 influence types----

for(ind in allInd){
  
  indIdx <- which(allInd==ind)
  session <- allIndInfo$session[which(allIndInfo$uniqueID==ind)]
  indData <-  spatialMetrics[which(spatialMetrics$indUniqID==ind),]
  
  if(saveImage)png(filename = paste0("figures/SuppMat/individual_logistic_fits/",ind,"_logisticFit.png"), width = 1000, height = 1000)
  par(mfrow=c(2,2), mar=c(5,5,4,2))
  
  #position turning influence
  data <- data.frame(x=indData$leftRightPosition,y=as.numeric(indData$groupTurnsRight))
  data <- data[which(!is.na(data$x) & !is.na(data$y)),]

  gamma <- 0.5
  
  maxRange <- ceiling(max(abs(quantile(data$x,c(0.01,0.99)))))
  
  influencePlots(x=data$x,y=data$y,binSize=round(maxRange/6),totRange=c(-maxRange,maxRange),cex.main=3,cex.lab=2.5,lwd=lwd,
                 xlab="Individual left-right position (m)",ylab="P(group turns right)",main="Position Turning Influence")
  segments(modelParam_PosTurn$quantile90[indIdx],-10,modelParam_PosTurn$quantile90[indIdx],modelParam_PosTurn$inflScore[indIdx],lty="dashed",col=4,lwd=3)
  segments(modelParam_PosTurn$quantile90[indIdx],modelParam_PosTurn$inflScore[indIdx],-100,modelParam_PosTurn$inflScore[indIdx],lty="dashed",col=4,lwd=3)
  lines(-100:100,flatLogis1Variable(-100:100,modelParam_PosTurn$alpha[indIdx],modelParam_PosTurn$beta[indIdx],gamma=gamma),col=2,lwd=lwd)
  legend("topright",legend=ind,bty="n",cex=3)
  legend("topleft",cex=2,legend=c(paste0("alpha = ",round(modelParam_PosTurn$alpha[indIdx],3)),paste0("beta = ",round(modelParam_PosTurn$beta[indIdx],3)),
                                    paste0("Likelihood = ",round(modelParam_PosTurn$likelihood[indIdx],0)),paste("Influence score = ",round(modelParam_PosTurn$inflScore[indIdx],2))),bty="n")
  abline(v=0,h=gamma,lty="dashed")
  
  #movement turning influence
  data <- data.frame(x=indData$leftRightMovement,y=as.numeric(indData$groupTurnsRight))
  data <- data[which(!is.na(data$x) & !is.na(data$y)),]
  
  gamma <- 0.5
  
  maxRange <- ceiling(max(abs(quantile(data$x,c(0.01,0.99)))))
  
  influencePlots(x=data$x,y=data$y,binSize=round(maxRange/6),totRange=c(-maxRange,maxRange),cex.main=3,cex.lab=2.5,lwd=lwd,
                 xlab="Individual left-right movement (m/min)",ylab="P(group turns right)",main="Movement Turning Influence")
  segments(modelParam_MovTurn$quantile90[indIdx],-10,modelParam_MovTurn$quantile90[indIdx],modelParam_MovTurn$inflScore[indIdx],lty="dashed",col=4,lwd=3)
  segments(modelParam_MovTurn$quantile90[indIdx],modelParam_MovTurn$inflScore[indIdx],-100,modelParam_MovTurn$inflScore[indIdx],lty="dashed",col=4,lwd=3)
  lines(-100:100,flatLogis1Variable(-100:100,modelParam_MovTurn$alpha[indIdx],modelParam_MovTurn$beta[indIdx],gamma=gamma),col=2,lwd=lwd)
  legend("topright",legend=ind,bty="n",cex=3)
  legend("topleft",cex=2,legend=c(paste0("alpha = ",round(modelParam_MovTurn$alpha[indIdx],3)),paste0("beta = ",round(modelParam_MovTurn$beta[indIdx],3)),
                                    paste0("Likelihood = ",round(modelParam_MovTurn$likelihood[indIdx],0)),paste("Influence score = ",round(modelParam_MovTurn$inflScore[indIdx],2))),bty="n")
  abline(v=0,h=gamma,lty="dashed")
  
  #position speeding influence
  data <- data.frame(x=indData$frontBackPosition,y=as.numeric(indData$groupSpeedsUp))
  data <- data[which(!is.na(data$x) & !is.na(data$y)),]
  
  gamma <- modelParam_PosSpeed$gamma[indIdx]
  
  maxRange <- ceiling(max(abs(quantile(data$x,c(0.01,0.99)))))
  
  influencePlots(x=data$x,y=data$y,binSize=round(maxRange/6),totRange=c(-maxRange,maxRange),cex.main=3,cex.lab=2.5,lwd=lwd,
                 xlab="Individual front-back position (m)",ylab="P(group speeds up)",main="Position speeding Influence")
  segments(modelParam_PosSpeed$quantile90[indIdx],-10,modelParam_PosSpeed$quantile90[indIdx],modelParam_PosSpeed$inflScore[indIdx],lty="dashed",col=4,lwd=3)
  segments(modelParam_PosSpeed$quantile90[indIdx],modelParam_PosSpeed$inflScore[indIdx],-100,modelParam_PosSpeed$inflScore[indIdx],lty="dashed",col=4,lwd=3)
  lines(-100:100,flatLogis1Variable(-100:100,modelParam_PosSpeed$alpha[indIdx],modelParam_PosSpeed$beta[indIdx],gamma=gamma),col=2,lwd=lwd)
  legend("topright",legend=ind,bty="n",cex=3)
  legend("topleft",cex=2,legend=c(paste0("alpha = ",round(modelParam_PosSpeed$alpha[indIdx],3)),paste0("beta = ",round(modelParam_PosSpeed$beta[indIdx],3)),
                                    paste0("Likelihood = ",round(modelParam_PosSpeed$likelihood[indIdx],0)),paste("Influence score = ",round(modelParam_PosSpeed$inflScore[indIdx],2))),bty="n")
  abline(v=0,h=gamma,lty="dashed")
  
  #movement speeding influence
  data <- data.frame(x=indData$frontBackMovement,y=as.numeric(indData$groupSpeedsUp))
  data <- data[which(!is.na(data$x) & !is.na(data$y)),]
  
  gamma <- modelParam_MovSpeed$gamma[indIdx]

  maxRange <- ceiling(max(abs(quantile(data$x,c(0.01,0.99)))))
  
  influencePlots(x=data$x,y=data$y,binSize=round(maxRange/6),totRange=c(-maxRange,maxRange),cex.main=3,cex.lab=2.5,lwd=lwd,
                 xlab="Individual front-back movement (m/min)",ylab="P(group speeds up)",main="Movement speeding Influence")
  segments(modelParam_MovSpeed$quantile90[indIdx],-10,modelParam_MovSpeed$quantile90[indIdx],modelParam_MovSpeed$inflScore[indIdx],lty="dashed",col=4,lwd=3)
  segments(modelParam_MovSpeed$quantile90[indIdx],modelParam_MovSpeed$inflScore[indIdx],-100,modelParam_MovSpeed$inflScore[indIdx],lty="dashed",col=4,lwd=3)
  lines(-100:100,flatLogis1Variable(-100:100,modelParam_MovSpeed$alpha[indIdx],modelParam_MovSpeed$beta[indIdx],gamma=gamma),col=2,lwd=lwd)
  legend("topright",legend=ind,bty="n",cex=3)
  legend("topleft",cex=2,legend=c(paste0("alpha = ",round(modelParam_MovSpeed$alpha[indIdx],3)),paste0("beta = ",round(modelParam_MovSpeed$beta[indIdx],3)),
                                    paste0("Likelihood = ",round(modelParam_MovSpeed$likelihood[indIdx],0)),paste("Influence score = ",round(modelParam_MovSpeed$inflScore[indIdx],2))),bty="n")
  abline(v=0,h=gamma,lty="dashed")
  
  if(saveImage)dev.off()
}

#----INFLUENCE SCORES WHEN VARYING DISCRETIZATION STEP----

for(discretizationStep in c(5,10,15,20)){
  
  load(paste0("output/Influence_logistic_model_fits_",discretizationStep,"m.RData"))
  load(paste0("output/spatialMetrics_",discretizationStep,"m.RData"))
  
  #movement turning influence scores
  
  if(saveImage)png(file = paste0("figures/SuppMat/varying_discretization_step_length/Movement_turning_influence_scores_",discretizationStep,"m.png"),units="in",res=100, width = 20, height = 30)
  
  par(mfrow=c(1,1), mar=c(5,7,4,2),cex=cex)
  
  yAx <- rev(1:nrow(modelParam_MovTurn)) - match(modelParam_MovTurn$session,unique(modelParam_MovTurn$session)) + 4
  pntSize = normalize(modelParam_MovTurn$N,1,3)
  
  plot(NULL,xlim=c(min(modelParam_MovTurn$inflScore-0.05),max(modelParam_MovTurn$inflScore+0.1)),ylim=c(1,max(yAx)-1),yaxt="n",ylab="",xlab="Movement turning influence score",
      cex.axis=2,cex.main=2,cex.lab=cex,main=paste0("Discretization step = ",discretizationStep,"m"))
  for(i in 1:nrow(modelParam_MovTurn)){
    
    ind = modelParam_MovTurn$ind[i]
    colo = statusCol[match(modelParam_MovTurn$status[i],status)]
    
    xAx <- modelParam_MovTurn$inflScore[i]
    
    segments(0,yAx[i],xAx,yAx[i],col="grey",lty="dotted",lwd=4)
    #abline(h=yAx[i],col="grey",lty="dotted",lwd=4)
    points(xAx,yAx[i],bg=colo,pch=21,cex=pntSize[i])
    
  }
  abline(h=yAx[which(diff(yAx)==-2)]-1,lty="dashed",lwd=4)
  axis(2,at=yAx,modelParam_MovTurn$ind,tick=F,las=1,cex.axis=1.1)
  #legend("topright",legend=status,pch=19,col=statusCol,pt.cex=cex)
  abline(v=0.5,lty="longdash",lwd=2)
  
  if(saveImage)dev.off()
  
  
  #movement speeding influence scores
  
  if(saveImage)png(file = paste0("figures/SuppMat/varying_discretization_step_length/Movement_speeding_influence_scores_",discretizationStep,"m.png"),units="in",res=100, width = 20, height = 30)
  
  par(mfrow=c(1,1), mar=c(5,7,4,2),cex=cex)
  
  yAx <- rev(1:nrow(modelParam_MovSpeed)) - match(modelParam_MovSpeed$session,unique(modelParam_MovSpeed$session)) + 4
  pntSize = normalize(modelParam_MovSpeed$N,1,3)
  
  plot(NULL,xlim=c(min(modelParam_MovSpeed$inflScore-0.05),max(modelParam_MovSpeed$inflScore+0.1)),ylim=c(1,max(yAx)-1),yaxt="n",ylab="",xlab="Movement speeding influence score",
       cex.axis=2,cex.lab=2,cex.main=2,cex.lab=cex,main=paste0("Discretization step = ",discretizationStep,"m"))
  for(i in 1:nrow(modelParam_MovSpeed)){
    
    ind = modelParam_MovSpeed$ind[i]
    colo = statusCol[match(modelParam_MovSpeed$status[i],status)]
    
    xAx <- modelParam_MovSpeed$inflScore[i]
    
    segments(0,yAx[i],xAx,yAx[i],col="grey",lty="dotted",lwd=4)
    #abline(h=yAx[i],col="grey",lty="dotted",lwd=4)
    points(xAx,yAx[i],bg=colo,pch=21,cex=pntSize[i])
    
  }
  sessionLines = yAx[which(diff(yAx)==-2)]-1
  abline(h=sessionLines,lty="dashed",lwd=4)
  axis(2,at=yAx,modelParam_MovSpeed$ind,tick=F,las=1,cex.axis=1.1)
  legend("topright",legend=status,pch=19,col=statusCol,pt.cex=cex)
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
  
  if(saveImage)png(file = paste0("figures/SuppMat/varying_discretization_step_length/PropFront_ridgeLine_",discretizationStep,"m.png"), width = 15, height = 20,units = "in",res=100)
  
  par(mar=c(5,8,4,2),mfrow=c(1,1))
  
  x <- unlist(windowPropFront)
  y <- factor(rep(allIndShort,sapply(windowPropFront,length)),levels=allIndShort)
  y<- y[which(!is.na(x))]
  x<- x[which(!is.na(x))]
  
  colo = statusCol[match(modelParam_MovSpeed$status,status)]
  
  myRidgeline(x,y,palette=colo,grouping=modelParam_MovSpeed$session,spacing=2,axes=T,xlab="Proportion of time in the front half",
              ylab="",meanNotMode = T,mode=T,modeCol="lightblue",yaxt="n",labCex=2,cex.lab=3,cex.axis=3,cex.main=3,main=paste0("Discretization step = ",discretizationStep,"m"))
  
  abline(v=0.5,lty="dashed",lwd=4)
  legend("topright",legend=status,pch=19,col=statusCol,cex=2,bg="white")
  
  if(saveImage)dev.off()
  
}

discretizationStep <- 10

load(paste0("output/spatialMetrics_",discretizationStep,"m.RData"))
load(paste0("output/Influence_logistic_model_fits_",discretizationStep,"m.RData"))

#----POSITION TURNING INFLUENCE SCORES----

if(saveImage)pdf(file = "figures/SuppMat/position_turning_and_speeding_influence/Position_turning_influence_scores.pdf", width = 20, height = 35)

par(mfrow=c(1,1), mar=c(5,6.3,4,2),cex=cex)

yAx <- rev(1:nrow(modelParam_PosTurn)) - match(modelParam_PosTurn$session,unique(modelParam_PosTurn$session)) + 4
pntSize = normalize(modelParam_PosTurn$N,1,3)

plot(NULL,xlim=c(min(modelParam_PosTurn$inflScore-0.05),max(modelParam_PosTurn$inflScore+0.1)),ylim=c(1,max(yAx)-1),yaxt="n",ylab="",xlab="Position turning influence score",
     main="",cex.axis=2,cex.lab=2.3)
for(i in 1:nrow(modelParam_PosTurn)){
  
  ind = modelParam_PosTurn$ind[i]
  colo = statusCol[match(modelParam_PosTurn$status[i],status)]
  
  xAx <- modelParam_PosTurn$inflScore[i]
  
  segments(0,yAx[i],xAx,yAx[i],col="grey",lty="dotted",lwd=4)
  #abline(h=yAx[i],col="grey",lty="dotted",lwd=4)
  points(xAx,yAx[i],bg=colo,pch=21,cex=pntSize[i])
  
}
abline(h=yAx[which(diff(yAx)==min(diff(yAx)))]-1,lty="dashed",lwd=4)
axis(2,at=yAx,gsub("20","",modelParam_MovTurn$ind),tick=F,las=1,cex.axis=1.5)
#legend("topright",legend=status,pch=19,col=statusCol,pt.cex=cex)
abline(v=0.5,lty="longdash",lwd=2)

if(saveImage)dev.off()


#----POSITION SPEEDING INFLUENCE SCORES----

if(saveImage)pdf(file = "figures/SuppMat/position_turning_and_speeding_influence/Position_speeding_influence_scores.pdf", width = 20, height = 35)

par(mfrow=c(1,1), mar=c(5,6.3,4,2),cex=cex)

yAx <- rev(1:nrow(modelParam_PosSpeed)) - match(modelParam_PosSpeed$session,unique(modelParam_PosSpeed$session)) + 4
pntSize = normalize(modelParam_PosSpeed$N,1,3)

plot(NULL,xlim=c(min(modelParam_PosSpeed$inflScore-0.1),max(modelParam_PosSpeed$inflScore+0.05)),ylim=c(1,max(yAx)-1),yaxt="n",ylab="",xlab="Position speeding influence score",
     main="",cex.axis=2,cex.lab=2.3)
for(i in 1:nrow(modelParam_PosSpeed)){
  
  ind = modelParam_PosSpeed$ind[i]
  colo = statusCol[match(modelParam_PosSpeed$status[i],status)]
  
  xAx <- modelParam_PosSpeed$inflScore[i]
  
  segments(0,yAx[i],xAx,yAx[i],col="grey",lty="dotted",lwd=4)
  #abline(h=yAx[i],col="grey",lty="dotted",lwd=4)
  points(xAx,yAx[i],bg=colo,pch=21,cex=pntSize[i])
  
}
sessionLines = yAx[which(diff(yAx)==min(diff(yAx)))]-1
abline(h=sessionLines,lty="dashed",lwd=4)
axis(2,at=yAx,gsub("20","",modelParam_MovTurn$ind),tick=F,las=1,cex.axis=1.5)
legend("topright",legend=status,pch=19,col=statusCol,pt.cex=cex,cex=1.3)
segments(unique(modelParam_PosSpeed$gamma),c(100,sessionLines),unique(modelParam_PosSpeed$gamma),c(sessionLines,-5),lty="longdash",lwd=2)

if(saveImage)dev.off()


#----distribution of group spreads----

nearNeigh <- list()
allDistances <- list()

for(session in sessions){
  
  if(saveImage)png(filename = paste0("figures/SuppMat/distributions/groupSpread_",session,".png"), width = 1000, height = 1000)
  
  par(mfrow=c(1,1), mar=c(5,5,4,2))
  
  longNames <- c(load(paste0("data/movement/level1/",session,"_COORDINATES_SYNCHED_level1.RData")))
  shortNames <- simplifyNames(pat=paste(session,"_",sep=""))
  
  allDistances[[session]] <- mapply(allDist,as.data.frame(allX),as.data.frame(allY),symetry=T,SIMPLIFY="array")
  
  nearNeigh[[session]] <- apply(allDistances[[session]],c(1,3),min,na.rm=T)
  nearNeigh[[session]][which(is.infinite(nearNeigh[[session]]))] <- NA
  
  spread <- apply(allDistances[[session]],3,mean,na.rm=T)
  
  hist(spread,breaks=50,main=session,xlab="Group spread (mean inter-individual distance in meters)",xlim=c(0,70),cex.main=3,cex.lab=2.5,cex.axis=2.5)
  
  if(saveImage)dev.off()
}


#----distribution of time scales for the given discretization step----

for(session in sessions){
  
  if(saveImage)png(filename = paste0("figures/SuppMat/distributions/discretizationTimeSteps_",session,".png"), width = 1000, height = 1000)
  
  par(mfrow=c(1,1), mar=c(5,5,4,2))
  
  sessionTable <- spatialMetrics[which(spatialMetrics$session==session),]
  sessionTable <- sessionTable[match(unique(sessionTable$t),sessionTable$t),]
  
  hist(sessionTable$futurStepDuration,breaks=50,xlim=c(0,1500),
       main=session,xlab="Duration of time steps for a spatial discretization of 10m (s)",cex.main=3,cex.lab=2.5,cex.axis=2.5)
  
  if(saveImage)dev.off()
}


#distribution of group turning angles

for(session in sessions){
  
  if(saveImage)png(filename = paste0("figures/SuppMat/distributions/turningAngle_",session,".png"), width = 1000, height = 1000)
  
  par(mfrow=c(1,1), mar=c(5,5,4,2))
  
  sessionTable <- spatialMetrics[which(spatialMetrics$session==session),]
  sessionTable <- sessionTable[match(unique(sessionTable$t),sessionTable$t),]
  
  metric <- sessionTable$groupTurnAngle
  
  hist(metric,breaks=50,
       main=session,xlab="Turning angle of group centroid (rad)",cex.main=3,cex.lab=2.5,cex.axis=2.5)
  
  if(saveImage)dev.off()
}





