#THE PURPOSE OF THIS SCRIPT IS TO GET CONFIDENCE INTERVALS ON THE INFLUENCE SCORES BY BOOTSTRAPPING THE DATA
#IT FIRST LOOKS AT THE DISTRIBUTION OF AUTOCORELATION LAG FOR THE TURN AND SPEED INFLUENCE ACROSSE DATES AND INDIVIDUALS
#THEN USES THE MEAN OF THESE DISTRIBUTIONS AS THE SIZE OF DATA CHUNKS THAT ARE DRAWN DURING THE BOOTSTRAP

setwd("C:/Users/baverly/Desktop/INFLUENCE_PAPER")
source("scripts/functions.R")

discretizationStep <- 10

load(paste0("output/spatialMetrics_",discretizationStep,"m.RData"))
load(paste0("output/Influence_logistic_model_fits_",discretizationStep,"m.RData"))

sessions <- c("HM2017","HM2019","L2019","ZU2021","NQ2021")
status = c("DominantF","DominantM","Adult","Yearling","Sub-Adult")
tlIdx <- match(unique(spatialMetrics$t),spatialMetrics$t)


# ----Selecting an appropriate chunk size for the bootstrap by looking at the autocorrelation lag of the turn and speed influence----

spatialMetrics$movTurn <- as.numeric((spatialMetrics$leftRightMovement>0) == spatialMetrics$groupTurnsRight)
spatialMetrics$movSpeed <- as.numeric((spatialMetrics$frontBackMovement>0) == spatialMetrics$groupSpeedsUp)

autoCorr <- data.frame(matrix(ncol=4,nrow=0))
colnames(autoCorr) <- c("ind","date","lagTurn","lagSpeed")
  
for(ind in unique(spatialMetrics$indUniqID)){
  
  indData <- spatialMetrics[which(spatialMetrics$indUniqID==ind),]
  
  indData <- indData[order(as.POSIXct(indData$t,tz="UTC")),]
  
  for(d in unique(indData$date)){
    
    subDate <- indData[which(indData$date==d),]
    
    if(length(which(!is.na(subDate$movTurn)))>0){
      
      seriesTurn <- subDate$movTurn
      seriesSpeed <- subDate$movSpeed
      
      corrTurn <- acf(seriesTurn,na.action=na.pass,lag.max=1500,main=ind)
      corrSpeed <- acf(seriesSpeed,na.action=na.pass,lag.max=1500,main=ind)
      
      significanceTurn <- qnorm((1 + 0.95)/2)/sqrt(sum(!is.na(seriesTurn)))
      significanceSpeed <- qnorm((1 + 0.95)/2)/sqrt(sum(!is.na(seriesSpeed)))
      
      idx <- nrow(autoCorr)+1
      autoCorr[idx,] <- NA
      autoCorr$ind[idx] <- ind
      autoCorr$date[idx] <- d
      lagTurn <- head(which(corrTurn$acf<significanceTurn),1)
      lagSpeed <- head(which(corrSpeed$acf<significanceSpeed),1)
      if(length(lagTurn>0)){
        autoCorr$lagTurn[idx] <- lagTurn
      }else{
        autoCorr$lagTurn[idx] <- NA
      }
      if(length(lagSpeed>0)){
        autoCorr$lagSpeed[idx] <- lagSpeed
      }else{
        autoCorr$lagSpeed[idx] <- NA
      }
    }
  }
}

#mean autocorrelation lag, that we are going to use as chunck size for the bootstrapping
print(mean(autoCorr$lagTurn,na.rm=T))
print(mean(autoCorr$lagSpeed,na.rm=T))


#----Performing the bootstrap----

chunk <- 240 #rounding it to a full minute to make it easier for the bootstrapping
iter <- 1000 #number of iterations

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
  
  #splitting the data into chunks of size chunk
  fact <- factor(rep(1:(totT/chunk),each=chunk))
  splittedX <- split(indData$leftRightMovement,fact)
  splittedIdx <- split(1:totT,fact)
  
  #only drawing from chunks that don't have only NAs
  drawFrom <- as.vector(which(sapply(splittedX,function(f)!all(is.na(f)))))

  for(i in 1:iter){
    
    print(i)
    
    #for each iteration, drawing N chunks with replacement, then recomputing the influence scores
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

#getting the 0.05-0.95% quantile of the influence score distribution from the bootstrap
modelParam_MovTurn$lowerCI <- sapply(1:nrow(modelParam_MovTurn),function(ind)quantile(sapply(movTurn_list,function(f)f[ind,"inflScore"]),0.05))
modelParam_MovTurn$upperCI <- sapply(1:nrow(modelParam_MovTurn),function(ind)quantile(sapply(movTurn_list,function(f)f[ind,"inflScore"]),0.95))
modelParam_MovSpeed$lowerCI <- sapply(1:nrow(modelParam_MovSpeed),function(ind)quantile(sapply(movSpeed_list,function(f)f[ind,"inflScore"]),0.05))
modelParam_MovSpeed$upperCI <- sapply(1:nrow(modelParam_MovSpeed),function(ind)quantile(sapply(movSpeed_list,function(f)f[ind,"inflScore"]),0.95))

save(modelParam_Total,modelParam_PosTurn,modelParam_MovTurn,modelParam_PosSpeed,modelParam_MovSpeed,file=paste0("output/Influence_logistic_model_fits_",discretizationStep,"m.RData"))
save(movTurn_list,movSpeed_list,file=paste0("output/bootstraps_",discretizationStep,"m.RData"))




