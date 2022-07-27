#THE PURPOSE OF THIS SCRIPT IS TO COMPUTE GROUP AND INDIVIDUAL SPATIAL METRICS (MAINLY VELOCITY VECTORS) BY SPATIALY DISCRETIZING THE TRAJECTORIES

setwd("D:/Meerkat_TurnSpeed_Influence")
source('scripts/functions.R')

sessions <- c('HM2017','HM2019','L2019','ZU2021',"NQ2021") 
discretizationStep = 10 #lag for the spatial discretization (in meters) - default 10

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
  
  #computing a bunch of individual and group spatial measures from the past and the future using spatial discretization
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
spatialMetrics$groupTurnsRight <- as.numeric(spatialMetrics$groupFuturVecPerpToGroupAxis < 0)

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

spatialMetrics$session <- factor(spatialMetrics$session,levels=sessions)

dir.create(path="output/",showWarnings = F)
save(file = paste0('output/spatialMetrics_',discretizationStep,'m.RData'), list = c('spatialMetrics','allIndInfo'))


