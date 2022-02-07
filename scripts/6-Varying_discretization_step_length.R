

setwd("C:/Users/baverly/Desktop/INFLUENCE_PAPER")
source('scripts/functions.R')

sessions <- c('HM2017','HM2019','L2019','ZU2021',"NQ2021")
status = c("DominantF","DominantM","Adult","Yearling","Sub-Adult")

#combining individual info tables for all sessions
print('creating info table for all sessions')
allIndInfo <- data.frame()

for(session in sessions){
  
  longNames <- c(load(paste0("data/movement/level1/",session,"_COORDINATES_SYNCHED_level1.RData")))
  shortNames <- simplifyNames(pat=paste(session,"_",sep=""))
  indInfo$session <- session
  indInfo$idx <- 1:nrow(indInfo)
  indInfo$uniqueID <- paste(indInfo$session,indInfo$idx,sep="_")
  allIndInfo <- rbind(allIndInfo, indInfo)
}

discretizationStep = 10 #lag for the spatial discretization (in meters) - default 10

for(discretizationStep  in c(5,15,20)){
  
  print('starting pre-processing')
  spatialMetrics <- data.frame()
  for(session in sessions){
    
    print(paste('session:',session))
    timestamp()
    
    #loading spatial data and making names easier
    longNames <- c(load(paste0("data/movement/level1/",session,"_COORDINATES_SYNCHED_level1.RData")))
    shortNames <- simplifyNames(pat=paste(session,"_",sep=""))
    
    #get number of individuals
    nInds <- nrow(indInfo)
    
    #create data frame for that session and put in features that don't depend on whether the focal is left out
    sessionTable <- data.frame(session = rep(session, length(timeLine)*nInds),date=rep(as.Date(timeLine),each=nInds), t = rep(timeLine, each = nInds), tIdx = rep(1:length(timeLine), each = nInds), indIdx = rep(1:nInds, length(timeLine)))
    sessionTable$status = rep(indInfo$status, length(timeLine))
    sessionTable$indUniqID <- paste(sessionTable$session, sessionTable$indIdx, sep = '_')
    sessionTable$statusUniqID <- paste(sessionTable$session, sessionTable$status, sep = '_')
    
    #computing a bunch of individual and group spatial measures from the past and the future
    spatialPast <- relativePos(allX,allY,step=discretizationStep,discrSpatial=T,futur=F,timeline=timeLine,discrByCentroid=F,centroidSpeed = T, removeInd = T) #metrics based on past
    spatialFutur <- relativePos(allX,allY,step=discretizationStep,discrSpatial=T,futur=T,timeline=timeLine,discrByCentroid=F,centroidSpeed = T, removeInd = T) #metrics based in future
    
    spatialPastNoRemove <- relativePos(allX,allY,step=discretizationStep,discrSpatial=T,futur=F,timeline=timeLine,discrByCentroid=T,centroidSpeed = T, removeInd = F) #metrics based on past
    rankAlongAxis <- apply(spatialPastNoRemove$relative_ind_X,2,normalize)
    
    #add the individual-level metrics to the table (which do not depend on whether focal is left out of group-level measures or not)
    sessionTable$pastStepDuration <- spatialPast$ind_stepDuration[cbind(sessionTable$indIdx, sessionTable$tIdx)]
    sessionTable$futurStepDuration <- spatialFutur$ind_stepDuration[cbind(sessionTable$indIdx, sessionTable$tIdx)]
    sessionTable$indSpeedPast <- spatialPast$ind_speed[cbind(sessionTable$indIdx, sessionTable$tIdx)]
    sessionTable$indSpeedFutur <- spatialFutur$ind_speed[cbind(sessionTable$indIdx, sessionTable$tIdx)]
    sessionTable$indHeadPast <- spatialPast$ind_heading[cbind(sessionTable$indIdx, sessionTable$tIdx)]
    sessionTable$indHeadFutur <- spatialFutur$ind_heading[cbind(sessionTable$indIdx, sessionTable$tIdx)]
    sessionTable$IndRankAlongMvmtAxis <- rankAlongAxis[cbind(sessionTable$indIdx, sessionTable$tIdx)]
    
    #add the group-level metrics to the table
    sessionTable$groupSpeedPast <- spatialPast$group_speed[cbind(sessionTable$indIdx, sessionTable$tIdx)]
    sessionTable$groupSpeedFutur <- spatialFutur$group_speed[cbind(sessionTable$indIdx, sessionTable$tIdx)]
    sessionTable$groupHeadPast <- spatialPast$group_heading[cbind(sessionTable$indIdx, sessionTable$tIdx)]
    sessionTable$groupHeadFutur <- spatialFutur$group_heading[cbind(sessionTable$indIdx, sessionTable$tIdx)]
    sessionTable$relX <- spatialPast$relative_ind_X[cbind(sessionTable$indIdx, sessionTable$tIdx)]
    sessionTable$relY <- spatialPast$relative_ind_Y[cbind(sessionTable$indIdx, sessionTable$tIdx)]
    sessionTable$relativeAngle <- spatialPast$relative_angle[cbind(sessionTable$indIdx, sessionTable$tIdx)]
    
    spatialMetrics <- rbind(spatialMetrics, sessionTable)
  }
  
  #COMPUTING LEADERSHIP METRICS
  
  #Get components of past individual motion parallel and perpendicular to past group heading
  dxIndPast <- cos(spatialMetrics$indHeadPast) * spatialMetrics$indSpeedPast #full x component (m / s)
  dyIndPast <- sin(spatialMetrics$indHeadPast) * spatialMetrics$indSpeedPast #full y component (m / s)
  dxGroupPast <- cos(spatialMetrics$groupHeadPast) #unit x component (unitless)
  dyGroupPast <- sin(spatialMetrics$groupHeadPast) #unit y component (unitless)
  spatialMetrics$indPastVecParallelToGroupAxis <- dxIndPast * dxGroupPast + dyIndPast * dyGroupPast #project ind past heading onto group past heading (m / s)
  spatialMetrics$indPastVecPerpToGroupAxis <- dyIndPast * dxGroupPast - dxIndPast * dyGroupPast #'"vector rejection" of ind past heading onto group past heading (projection onto the perpendicular) (m / s)
  
  #Get components of group motion parallel and perpendicular to past group heading
  dxGroupFutur <- cos(spatialMetrics$groupHeadFutur) * spatialMetrics$groupSpeedFutur #full x component (m / s)
  dyGroupFutur <- sin(spatialMetrics$groupHeadFutur) * spatialMetrics$groupSpeedFutur #full y component (m / s)
  spatialMetrics$groupFuturVecParallelToGroupAxis <- dxGroupFutur * dxGroupPast + dyGroupFutur * dyGroupPast #projection of group future heading into group past heading (m / s)
  spatialMetrics$groupFuturVecPerpToGroupAxis <- dyGroupFutur * dxGroupPast - dxGroupFutur * dyGroupPast #"vector rejection" of group future heading onto group past heading (projection onto the perpendicular) (m / s)
  
  #Get how much the group speed up or slowed down along its direction of travel
  spatialMetrics$groupSpeedUpAlongGroupAxis <- (spatialMetrics$groupFuturVecParallelToGroupAxis - spatialMetrics$groupSpeedPast) # (m / s)
  
  #Get the difference between individual's speed along group past axis and group speed along its axis
  spatialMetrics$indGroupSpeedDiffAlongGroupAxis <- (spatialMetrics$indPastVecParallelToGroupAxis - spatialMetrics$groupSpeedPast) # (m / s)
  
  #group turn angle
  spatialMetrics$groupTurnAngle <- atan2(spatialMetrics$groupFuturVecPerpToGroupAxis, spatialMetrics$groupFuturVecParallelToGroupAxis)
  
  #angle between individual position and group direction
  spatialMetrics$indAngleFromGroupMovement <- atan2(spatialMetrics$relY,spatialMetrics$relX)
  
  #binary, wether individual is in the front half
  spatialMetrics$inFrontHalf <- spatialMetrics$relX > 0
  
  #binary, wether individual is in the front quarter
  spatialMetrics$inFrontQuarter <- abs(spatialMetrics$indAngleFromGroupMovement) < pi/4
  
  #binary, whether group sped up or slowed down
  spatialMetrics$groupSpeedsUp <- as.numeric(spatialMetrics$groupSpeedUpAlongGroupAxis > 0)
  
  #binary, whether group turned right or not
  spatialMetrics$groupTurnsRight <- spatialMetrics$groupFuturVecPerpToGroupAxis < 0
  
  #Individul left-right position
  spatialMetrics$leftRightPosition <- -spatialMetrics$relY
  
  #Individual left-right movement
  spatialMetrics$leftRightMovement <- -spatialMetrics$indPastVecPerpToGroupAxis*60
  
  #Individul front-back position
  spatialMetrics$frontBackPosition <- spatialMetrics$relX
  
  #Individual front-back movement
  spatialMetrics$frontBackMovement <- spatialMetrics$indGroupSpeedDiffAlongGroupAxis*60
  
  spatialMetrics$t <- as.POSIXct(spatialMetrics$t,tz="UTC")
  spatialMetrics <- spatialMetrics[order(spatialMetrics$t),]
  
  save(file = paste0('output/spatialMetrics_',discretizationStep,'m.RData'), list = c('spatialMetrics','allIndInfo'))
  
  
  #----MODELLING OF GROUP OUTCOME (TURNING OR SPEEDING) AS A FUNCTION OF BOTH POSITION AND MOVEMENT (ALL INDIVDUALS COMBINED)----
  
  modelParam_Total <-  data.frame(model=c("Turning_position+movement","Speeding_position+movement"),session="all",N=NA,alpha=NA,beta1=NA,beta2=NA,gamma=NA,likelihood=NA)
  
  tlIdx <- match(unique(spatialMetrics$t),spatialMetrics$t)
  pSpeedUpTotal <- mean(spatialMetrics$groupSpeedsUp[tlIdx],na.rm=T)
  
  #position & movement turning influence
  data1 <- data.frame(x1=spatialMetrics$leftRightMovement,x2=spatialMetrics$leftRightPosition,y=as.numeric(spatialMetrics$groupTurnsRight))
  data1 <- data1[which(!is.na(data1$x1) & !is.na(data1$x2) & !is.na(data1$y)),]
  
  gamma1 <- 0.5
  
  optim1 <- optim(par=c(0.5,0.5,0.5),fn=logLikelihood2Variables,data=data1,gamma=gamma1,fun=flatLogis2Variables)
  
  modelParam_Total$N[1] <- nrow(data1)
  modelParam_Total$alpha[1] <- optim1$par[1]
  modelParam_Total$beta1[1] <- optim1$par[2]
  modelParam_Total$beta2[1] <- optim1$par[3]
  modelParam_Total$gamma[1] <- gamma1
  modelParam_Total$likelihood[1] <- optim1$value
  
  #position & movement speeding influence
  data2 <- data.frame(x1=spatialMetrics$frontBackMovement,x2=spatialMetrics$frontBackPosition,y=as.numeric(spatialMetrics$groupSpeedsUp))
  data2 <- data2[which(!is.na(data2$x1) & !is.na(data2$x2) & !is.na(data2$y)),]
  
  gamma2 <- pSpeedUpTotal
  
  optim2 <- optim(par=c(0.5,0.5,0.5),fn=logLikelihood2Variables,data=data2,fun=flatLogis2Variables,gamma=gamma2)
  
  modelParam_Total$N[2] <- nrow(data2)
  modelParam_Total$alpha[2] <- optim2$par[1]
  modelParam_Total$beta1[2] <- optim2$par[2]
  modelParam_Total$beta2[2] <- optim2$par[3]
  modelParam_Total$gamma[2] <- gamma2
  modelParam_Total$likelihood[2] <- optim2$value
  
  
  #----MODELLING OF POSITION TURNING INFLUENCE, MOVEMENT TURNING INFLUENCE, POSITION SPEEDING INFLUENCE, MOVEMENT SPEEDING INFLUENCE(SEPARETELLY FOR EACH INDIVIDUAL)----
  
  allInd <- allIndInfo$uniqueID
  
  #overall probability to speed up for each group
  pGroupSpeedUp <- tapply(spatialMetrics$groupSpeedsUp[tlIdx], spatialMetrics$session[tlIdx],mean,na.rm=T)
  
  modelParam_PosTurn <- modelParam_MovTurn <- modelParam_PosSpeed <- modelParam_MovSpeed <- data.frame(ind=allIndInfo$uniqueID,session=allIndInfo$session,status=allIndInfo$status,N=NA,alpha=NA,beta=NA,gamma=NA,likelihood=NA)
  for(ind in allInd){
    
    print(ind)
    
    indIdx <- which(allInd==ind)
    session <- allIndInfo$session[which(allIndInfo$uniqueID==ind)]
    indData <-  spatialMetrics[which(spatialMetrics$indUniqID==ind),]
    
    #position turning influence
    data <- data.frame(x=indData$leftRightPosition,y=as.numeric(indData$groupTurnsRight))
    data <- data[which(!is.na(data$x) & !is.na(data$y)),]
    
    if(nrow(data)==0){
      next()
    }
    
    gamma <- 0.5
    optim <- optim(par=c(0.5,0.5),fn=logLikelihood1Variable,data=data,fun=flatLogis1Variable,gamma=gamma)
    
    modelParam_PosTurn$N[indIdx] <- nrow(data)
    modelParam_PosTurn$alpha[indIdx] <- optim$par[1]
    modelParam_PosTurn$beta[indIdx] <- optim$par[2]
    modelParam_PosTurn$gamma[indIdx] <- gamma
    modelParam_PosTurn$likelihood[indIdx] <- optim$value
    
    #movement turning influence
    data <- data.frame(x=indData$leftRightMovement,y=as.numeric(indData$groupTurnsRight))
    data <- data[which(!is.na(data$x) & !is.na(data$y)),]
    
    gamma <- 0.5
    optim <- optim(par=c(0.5,0.5),fn=logLikelihood1Variable,data=data,fun=flatLogis1Variable,gamma=gamma)
    
    modelParam_MovTurn$N[indIdx] <- nrow(data)
    modelParam_MovTurn$alpha[indIdx] <- optim$par[1]
    modelParam_MovTurn$beta[indIdx] <- optim$par[2]
    modelParam_MovTurn$gamma[indIdx] <- gamma
    modelParam_MovTurn$likelihood[indIdx] <- optim$value
    
    #position speeding influence
    data <- data.frame(x=indData$frontBackPosition,y=as.numeric(indData$groupSpeedsUp))
    data <- data[which(!is.na(data$x) & !is.na(data$y)),]
    
    gamma <- pGroupSpeedUp[which(names(pGroupSpeedUp)==session)]
    optim <- optim(par=c(0.5,0.5),fn=logLikelihood1Variable,data=data,fun=flatLogis1Variable,gamma=gamma)
    
    modelParam_PosSpeed$N[indIdx] <- nrow(data)
    modelParam_PosSpeed$alpha[indIdx] <- optim$par[1]
    modelParam_PosSpeed$beta[indIdx] <- optim$par[2]
    modelParam_PosSpeed$gamma[indIdx] <- gamma
    modelParam_PosSpeed$likelihood[indIdx] <- optim$value
    
    #movement speeding influence
    data <- data.frame(x=indData$frontBackMovement,y=as.numeric(indData$groupSpeedsUp))
    data <- data[which(!is.na(data$x) & !is.na(data$y)),]
    
    gamma <- pGroupSpeedUp[which(names(pGroupSpeedUp)==session)]
    optim <- optim(par=c(0.5,0.5),fn=logLikelihood1Variable,data=data,fun=flatLogis1Variable,gamma=gamma)
    
    modelParam_MovSpeed$N[indIdx] <- nrow(data)
    modelParam_MovSpeed$alpha[indIdx] <- optim$par[1]
    modelParam_MovSpeed$beta[indIdx] <- optim$par[2]
    modelParam_MovSpeed$gamma[indIdx] <- gamma
    modelParam_MovSpeed$likelihood[indIdx] <- optim$value
    
  }
  
  modelParam_PosTurn <- modelParam_PosTurn[which(!is.na(modelParam_PosTurn$alpha)),]
  modelParam_MovTurn <- modelParam_MovTurn[which(!is.na(modelParam_MovTurn$alpha)),]
  modelParam_PosSpeed <- modelParam_PosSpeed[which(!is.na(modelParam_PosSpeed$alpha)),]
  modelParam_MovSpeed <- modelParam_MovSpeed[which(!is.na(modelParam_MovSpeed$alpha)),]
  
  
  quantileValue = 0.9
  
  quantileLeftRightPos <- tapply(abs(spatialMetrics$leftRightPosition),spatialMetrics$session,quantile,probs=quantileValue,na.rm=T)
  quantileLeftRightMov <- tapply(abs(spatialMetrics$leftRightMovement),spatialMetrics$session,quantile,probs=quantileValue,na.rm=T)
  quantileFrontBackPos <- tapply(abs(spatialMetrics$frontBackPosition),spatialMetrics$session,quantile,probs=quantileValue,na.rm=T)
  quantileFrontBackMov <- tapply(abs(spatialMetrics$frontBackMovement),spatialMetrics$session,quantile,probs=quantileValue,na.rm=T)
  
  modelParam_PosTurn$quantile90 <- rep(quantileLeftRightPos,as.numeric(table(modelParam_PosTurn$session)))
  modelParam_MovTurn$quantile90 <- rep(quantileLeftRightMov,as.numeric(table(modelParam_MovTurn$session)))
  modelParam_PosSpeed$quantile90 <- rep(quantileFrontBackPos,as.numeric(table(modelParam_PosSpeed$session)))
  modelParam_MovSpeed$quantile90 <- rep(quantileFrontBackMov,as.numeric(table(modelParam_MovSpeed$session)))
  
  modelParam_PosTurn$inflScore <- flatLogis1Variable(modelParam_PosTurn$quantile90,modelParam_PosTurn$alpha,modelParam_PosTurn$beta,modelParam_PosTurn$gamma)
  modelParam_MovTurn$inflScore <- flatLogis1Variable(modelParam_MovTurn$quantile90,modelParam_MovTurn$alpha,modelParam_MovTurn$beta,modelParam_MovTurn$gamma)
  modelParam_PosSpeed$inflScore <- flatLogis1Variable(modelParam_PosSpeed$quantile90,modelParam_PosSpeed$alpha,modelParam_PosSpeed$beta,modelParam_PosSpeed$gamma)
  modelParam_MovSpeed$inflScore <- flatLogis1Variable(modelParam_MovSpeed$quantile90,modelParam_MovSpeed$alpha,modelParam_MovSpeed$beta,modelParam_MovSpeed$gamma)
  
  
  save(modelParam_Total,modelParam_PosTurn,modelParam_MovTurn,modelParam_PosSpeed,modelParam_MovSpeed,file=paste0("output/Influence_logistic_modelFits_",discretizationStep,"m.RData"))

}
