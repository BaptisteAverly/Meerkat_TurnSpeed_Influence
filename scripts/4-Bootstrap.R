
setwd("C:/Users/baverly/Desktop/INFLUENCE_PAPER")
source("scripts/functions.R")

discretizationStep <- 10

load(paste0("output/spatialMetrics_",discretizationStep,"m.RData"))
load(paste0("output/Influence_logistic_model_fits_",discretizationStep,"m.RData"))

sessions <- c("HM2017","HM2019","L2019","ZU2021","NQ2021")
status = c("DominantF","DominantM","Adult","Yearling","Sub-Adult")
tlIdx <- match(unique(spatialMetrics$t),spatialMetrics$t)


#spatialMetrics$groupTurnsRight <- as.numeric(spatialMetrics$groupTurnsRight)
#spatialMetrics$groupTurnsRight[which(spatialMetrics$groupTurnsRight==0)] <- -1
#spatialMetrics$groupSpeedsUp[which(spatialMetrics$groupSpeedsUp==0)] <- -1
#
#spatialMetrics$movTurn <- sign(spatialMetrics$leftRightMovement) * sign(spatialMetrics$groupTurnsRight)
#spatialMetrics$movSpeed <- sign(spatialMetrics$frontBackMovement) * sign(spatialMetrics$groupSpeedsUp)
#
#autoCorr <- data.frame(matrix(ncol=3,nrow=0))
#colnames(autoCorr) <- c("ind","date","lag")
#  
#for(ind in unique(spatialMetrics$indUniqID)){
#  
#  indData <- spatialMetrics[which(spatialMetrics$indUniqID==ind),]
#  
#  indData <- indData[order(as.POSIXct(spatialMetrics$t,tz="UTC")),]
#  
#  for(d in unique(indData$date)){
#    
#    subDate <- indData[which(indData$date==d),]
#    
#    if(length(which(!is.na(subDate$movTurn)))>0){
#      
#      series <- subDate$movTurn
#      
#      corr <- acf(series,na.action=na.pass,lag.max=600,main=ind)
#      
#      significance <- qnorm((1 + 0.95)/2)/sqrt(sum(!is.na(series)))
#      
#      idx <- nrow(autoCorr)+1
#      autoCorr[idx,] <- NA
#      autoCorr$ind[idx] <- ind
#      autoCorr$date[idx] <- d
#      lag <- head(which(corr$acf<significance),1)
#      if(length(lag>0)){
#        autoCorr$lag[idx] <- lag
#      }else{
#        autoCorr$lag[idx] <- NA
#      }
#    }
#  }
#}


chunk <- 240
iter <- 1000

movTurn_list <- movSpeed_list <- data.frame(ind=allIndInfo$uniqueID,session=allIndInfo$session,status=allIndInfo$status,N=NA,alpha=NA,beta=NA,gamma=NA,likelihood=NA)
movTurn_list <- lapply(1:iter, function(x) movTurn_list)
movSpeed_list <- lapply(1:iter, function(x) movSpeed_list)

#overall probability to speed up for each group
pGroupSpeedUp <- tapply(spatialMetrics$groupSpeedsUp[tlIdx], spatialMetrics$session[tlIdx],mean,na.rm=T)
quantileValue = 0.9
quantileLeftRightMov <- tapply(abs(spatialMetrics$leftRightMovement),spatialMetrics$session,quantile,probs=quantileValue,na.rm=T)
quantileFrontBackMov <- tapply(abs(spatialMetrics$frontBackMovement),spatialMetrics$session,quantile,probs=quantileValue,na.rm=T)

allInd <- unique(spatialMetrics$indUniqID)
for(ind in allInd){
  
  print(ind)
  indIdx <- which(allInd==ind)
  session <- allIndInfo$session[which(allIndInfo$uniqueID==ind)]
  sessionIdx <- which(sessions == session)
  indData <-  spatialMetrics[which(spatialMetrics$indUniqID==ind),]
  indData <- indData[order(as.POSIXct(indData$t)),]
  
  totT <- nrow(indData)
  
  if(all(is.na(indData$leftRightMovement)))next()
  
  fact <- factor(rep(1:(totT/chunk),each=chunk))
  
  splittedX <- split(indData$leftRightMovement,fact)
  splittedIdx <- split(1:totT,fact)
  
  drawFrom <- as.vector(which(sapply(splittedX,function(f)!all(is.na(f)))))

  for(i in 1:iter){
    
    print(i)
    boot <- as.vector(unlist(splittedIdx[sample(drawFrom,length(drawFrom),replace=T)]))
    
    #movement turning influence
    data <- data.frame(x=indData$leftRightMovement[boot],y=as.numeric(indData$groupTurnsRight[boot]))
    data <- data[which(!is.na(data$x) & !is.na(data$y)),]
    
    gamma <- 0.5
    optim <- optim(par=c(0.5,0.5),fn=logLikelihood1Variable,data=data,fun=flatLogis1Variable,gamma=gamma)
    
    movTurn_list[[i]]$N[indIdx] <- nrow(data)
    movTurn_list[[i]]$alpha[indIdx] <- optim$par[1]
    movTurn_list[[i]]$beta[indIdx] <- optim$par[2]
    movTurn_list[[i]]$gamma[indIdx] <- gamma
    movTurn_list[[i]]$likelihood[indIdx] <- optim$value
    movTurn_list[[i]]$quantile90[indIdx] <- quantileLeftRightMov[session]
    movTurn_list[[i]]$inflScore[indIdx] <- flatLogis1Variable(quantileLeftRightMov[session],optim$par[1],optim$par[2],gamma)
    
    #movement speeding influence
    data <- data.frame(x=indData$frontBackMovement[boot],y=as.numeric(indData$groupSpeedsUp[boot]))
    data <- data[which(!is.na(data$x) & !is.na(data$y)),]
    
    gamma <- pGroupSpeedUp[which(names(pGroupSpeedUp)==session)]
    optim <- optim(par=c(0.5,0.5),fn=logLikelihood1Variable,data=data,fun=flatLogis1Variable,gamma=gamma)
    
    movSpeed_list[[i]]$N[indIdx] <- nrow(data)
    movSpeed_list[[i]]$alpha[indIdx] <- optim$par[1]
    movSpeed_list[[i]]$beta[indIdx] <- optim$par[2]
    movSpeed_list[[i]]$gamma[indIdx] <- gamma
    movSpeed_list[[i]]$likelihood[indIdx] <- optim$value
    movSpeed_list[[i]]$quantile90[indIdx] <- quantileFrontBackMov[session]
    movSpeed_list[[i]]$inflScore[indIdx] <- flatLogis1Variable(quantileFrontBackMov[session],optim$par[1],optim$par[2],gamma)
    
  }
}

for(i in 1:iter){
  
  movTurn_list[[i]] <- movTurn_list[[i]][which(!is.na(movTurn_list[[i]]$N)),]
  movSpeed_list[[i]] <- movSpeed_list[[i]][which(!is.na(movSpeed_list[[i]]$N)),]
  
}

modelParam_MovTurn$lowerCI <- sapply(1:nrow(modelParam_MovTurn),function(ind)quantile(sapply(movTurn_list,function(f)f[ind,"inflScore"]),0.05))
modelParam_MovTurn$upperCI <- sapply(1:nrow(modelParam_MovTurn),function(ind)quantile(sapply(movTurn_list,function(f)f[ind,"inflScore"]),0.95))
modelParam_MovSpeed$lowerCI <- sapply(1:nrow(modelParam_MovSpeed),function(ind)quantile(sapply(movSpeed_list,function(f)f[ind,"inflScore"]),0.05))
modelParam_MovSpeed$upperCI <- sapply(1:nrow(modelParam_MovSpeed),function(ind)quantile(sapply(movSpeed_list,function(f)f[ind,"inflScore"]),0.95))

save(modelParam_Total,modelParam_PosTurn,modelParam_MovTurn,modelParam_PosSpeed,modelParam_MovSpeed,file=paste0("output/Influence_logistic_model_fits_",discretizationStep,"m.RData"))
save(movTurn_list,movSpeed_list,file=paste0("output/bootstraps_",discretizationStep,"m.RData"))




