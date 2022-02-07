
setwd("C:/Users/baverly/Desktop/INFLUENCE_PAPER")
source("scripts/functions.R")

discretizationStep <- 10

load(paste0("output/spatialMetrics_",discretizationStep,"m.RData"))

sessions <- c("HM2017","HM2019","L2019","ZU2021","NQ2021")
status = c("DominantF","DominantM","Adult","Yearling","Sub-Adult")
tlIdx <- match(unique(spatialMetrics$t),spatialMetrics$t)

#----MODELLING OF GROUP OUTCOME (TURNING OR SPEEDING) AS A FUNCTION OF BOTH POSITION AND MOVEMENT (ALL INDIVDUALS COMBINED)----

modelParam_Total <-  data.frame(model=c("Turning_position+movement","Speeding_position+movement"),session="all",N=NA,alpha=NA,beta1=NA,beta2=NA,gamma=NA,likelihood=NA)

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


save(modelParam_Total,modelParam_PosTurn,modelParam_MovTurn,modelParam_PosSpeed,modelParam_MovSpeed,file=paste0("output/Influence_logistic_model_fits_",discretizationStep,"m.RData"))





