library(stringr)

#function to simplify the name of certain objects in the environment (removing "pat" argument from their names) if reverse=F
#OR adding "pat" argument at the beggining of all object provided in argument "names" if reverse=T
#returns a vector of the new names
simplifyNames <- function(pat,reverse=F,items=ls(name=.GlobalEnv, pattern=pat)){
  
  for(item in items){
    if(reverse){
      assign(paste(pat,item,sep=""),get(item,envir=.GlobalEnv),envir=.GlobalEnv)
    }else{
      assign(str_remove(item,pat),get(item,envir=.GlobalEnv),envir=.GlobalEnv)
    }
    rm(list=item,envir=.GlobalEnv)
  }
  if(reverse){
    return(paste(pat,items,sep=""))
  }else{
    return(str_remove(items,pat))
  }
}


#quick day index
dayIndex <- function(day){
  if(!is.numeric(day) | length(day)==0){
    return(NA)
  }else{
    dayIdx[day]:(dayIdx[day+1]-1)
  }
}

#calculating various group- and individual- level spatial metrics 
#input: two matrices of UTM coordinates (one for x, one for y), with rows= inds (N), col=time frames (T)
#       discrSpatial: boolean, whether to calculate speed etc with time intervals or using spatial discretization
#       step:  the step used to calculate the velocity vectors, 
#              either in seconds for temporal discretization (if discrSpatial is True) 
#              or in meters for spatial discretization (if discrSpatial is True)
#       accelerationStep: time step in seconds used for temporal discretization to calculate acceleration
#       futur: boolean, whether to calculate speed, direction etc... further or backward
#       normalize: boolean, indicates if the direction vectors should be normalize to unit vectors
#       timeLine: vector of timeStamps to be used if there are several days
#       removeInd: numeric vector of one or more individual to remove from the calculation of the centroid OR boolean specifying remove no one (F) or remove each individual sequentially for calculating group metrics (T)
#       centroidSpeed: boolean, wether to calculate the speed of the group by calculating the speed of the centroid (T) or the average speed of all individuals (F)
#       discrIdx: vector of size N giving the indices for discretization if these were computed outside the function
#       discrByCentroid: boolean, wether to use the centroid track for spatial discretization (T) or calculate separate indices for each individual (F)
#output: a list containing several matrices of group and individual movement metrics all with the same dimensions N*T. 
#For group metrics, rows will be different only if removeInd = T :
#     group_discrIdx: indices used for the discretization of individual tracks
#     ind_discrIdx: indices used for the discretization of individual tracks
#     group_stepLength: distance in meters travelled by the centroid between each discretized step (uniform and equal to step argument if discrSpatial = T)
#     ind_stepLength: distance in meters travelled by individuals between each discretized step (uniform and equal to step argument only if discrSpatial = T AND discrByCentroid = F)
#     group_stepDuration: time in seconds between each discretized step (uniform and equal to step argument if discrSpatial = F)
#     ind_stepDuration: time in seconds between each discretized step (uniform and equal to step argument if discrSpatial = F AND discrByCentroid = F)
#     group_speed: speed of the centroid for each fix
#     ind_speed: individual speed
#     group_acc: acceleration of the centroid for each fix
#     ind_acc: individual acceleration
#     distance_to_centroid: distance of every individual from centroid at times t
#     group_heading: angles between direction vector of centroid and x-axis for each fix
#     group_direction_x: component of the heading of the group
#     group_direction_y: component of the heading of the group
#     group_polarization: polarization of the group
#     ind_heading: angles between direction vector of individuals and x-axis for each fix
#     ind_direction_x: x component of the heading of the individuals
#     ind_direction_y: y component of the heading of the individuals
#     relative_ind_x: x coordinates of individuals in the rotated frame (x-axis=direction of centroid)
#     relative_ind_Y: y coordinates of individuals in the rotated frame (x-axis=direction of centroid)
#     relative_angle: relative angle of each ind to direction of movement 
#     northX: x component of the unit vector pointing north in the rotated frame 
#     northY: y component of the unit vector pointing north in the rotated frame 
relativePos <- function(x,y,discrSpatial=F,step,accelerationStep = 10,futur=T,normalize=T,timeline=rep("1970-01-01",ncol(x)),removeInd=F,centroidSpeed=T,discrIdx=NULL,discrByCentroid=T,...){
  
  n<-nrow(x)
  times<-ncol(x)
  multiCentroid = F
  #position of centroid
  if(!removeInd){
    xc <- rbind(colMeans(x,na.rm=T))
    xc <- xc[rep(1,n),]
    yc <- rbind(colMeans(y,na.rm=T))
    yc <- yc[rep(1,n),]
  }else if(is.numeric(removeInd)){
    xc<-rbind(colMeans(x[-removeInd,],na.rm=T))
    xc <- xc[rep(1,n),]
    yc<-rbind(colMeans(y[-removeInd,],na.rm=T))
    yc <- yc[rep(1,n),]
  }else{
    xc <- t(sapply(1:n,function(ind)colMeans(x[-ind,],na.rm=T)))
    yc <- t(sapply(1:n,function(ind)colMeans(y[-ind,],na.rm=T)))
    multiCentroid = T
  }
  
  #getting the discretization index
  #either spatial or temporal
  #for spatial, either based only on centroid trajectory (same for all individuals), or calculated independantly for each individual
  if(is.null(discrIdx)){
    
    #discretization of centroid track
    if(multiCentroid){
      centroidDiscr <- t(sapply(1:n,function(ind)discretization(xc[ind,],yc[ind,],step,spatial=discrSpatial,futur=futur,timeline=timeline,...)))
      centroidSteps <- t(apply(centroidDiscr,1,function(f)abs(f - (1:length(f)))))
    }else{
      centroidDiscr <- rbind(discretization(xc[1,],yc[1,],step,spatial=discrSpatial,futur=futur,timeline=timeline,...))
      centroidSteps <- abs(centroidDiscr - (1:length(centroidDiscr)))
      centroidDiscr <- centroidDiscr[rep(1,n),]
      centroidSteps <- centroidSteps[rep(1,n),]
    }
    
    #discretization of individual tracks
    if(!discrSpatial | discrByCentroid){
      discrIdx <- centroidDiscr
      steps <- centroidSteps
    }else{
      discrIdx <- t(sapply(1:n,function(ind)discretization(x[ind,],y[ind,],step,spatial=discrSpatial,futur=futur,timeline=timeline,...)))
      steps <- t(abs(t(discrIdx) - (1:ncol(discrIdx))))
    }
  }else{
    centroidDiscr <- cbind(discrIdx)
    centroidSteps <- abs(centroidDiscr - (1:length(centroidDiscr)))
    centroidDiscr <- centroidDiscr[rep(1,n),]
    centroidSteps <- centroidSteps[rep(1,n),]
    
    discrIdx <- centroidDiscr
    steps <- centroidSteps 
  }
  
  if(futur)futurMod <- 1 else futurMod <- -1
  #group and ind direction
  groupDirecX <- t(sapply(1:n,function(ind)(xc[ind,centroidDiscr[ind,]]-xc[ind,])*futurMod))
  groupDirecY <- t(sapply(1:n,function(ind)(yc[ind,centroidDiscr[ind,]]-yc[ind,])*futurMod))
  
  indDirecX <- t(sapply(1:n,function(ind)(x[ind,discrIdx[ind,]]-x[ind,])*futurMod))
  indDirecY <- t(sapply(1:n,function(ind)(y[ind,discrIdx[ind,]]-y[ind,])*futurMod))
  
  #distance to centroid
  dist2centr <- sqrt((x-xc)^2 + (y-yc)^2)
  
  #linear distance travelled + speed between each step
  indDist <- sqrt(indDirecX^2+indDirecY^2)
  
  groupDist <- sqrt(groupDirecX^2+groupDirecY^2) 
  
  indSpeed <- indDist/steps
  if(centroidSpeed){
    groupSpeed<-groupDist/centroidSteps
  }else{
    groupSpeed <- colMeans(indSpeed,na.rm=T)
  }
  
  #group and individual speed
  discrAcc <- discretization(xc[1,],yc[1,],accelerationStep,spatial=F,futur=futur,timeline=timeline)
  
  groupAcc <- ((groupSpeed[,discrAcc] - groupSpeed) * (futurMod))/accelerationStep
  indAcc <- ((indSpeed[,discrAcc] - indSpeed) * (futurMod))/accelerationStep
  
  if(normalize){
    #we divide by the distance to obtain unit vectors
    groupDirecX<-groupDirecX/groupDist
    groupDirecY<-groupDirecY/groupDist
    indDirecX <- indDirecX / indDist
    indDirecY <- indDirecY / indDist
  }
  indAngles <- atan2(indDirecY,indDirecX)
  
  groupPolar <- sqrt(colSums(indDirecX,na.rm = T)^2 + colSums(indDirecY,na.rm = T)^2) / apply(x,2,function(f)length(which(!is.na(f))))
  
  #position of individuals relative to centroid
  relativeX <- x-xc#if positive --> ind East of centroid
  relativeY <- y-yc #if positive --> ind North
  
  #ROTATING THE GROUP REFERENCE FRAME
  #here we want that at each time, the x-axis correspond to the direction of travel of the group
  
  #getting the angle of group direction relative to x-axis
  angles<-atan2(groupDirecY,groupDirecX) #function gives angle between x-axis and a vector
  
  #then do the rotation using matrix rotation formula:
  # = we calculate the coordinates of each individual in the new system
  
  northX <- rep(0,times)
  northY <- rep(1,times)
  northX<-cos(angles)*northX+sin(angles)*northY
  northY<-(-sin(angles)*northX)+cos(angles)*northY
  
  #we use the relative distance of ind to centroid so center of system at all times = centroid
  rotatX <- cos(angles)*relativeX+sin(angles)*relativeY
  rotatY <- (-sin(angles)*relativeX)+cos(angles)*relativeY
  relAngle <- atan2(rotatY,rotatX)
  #--> if ind at the origin = at the centroid
  #--> if ind on the right = at the front of the group
  #--> if ind on the top = at the left of the group
  
  list(group_discrIdx=centroidDiscr,ind_discrIdx=discrIdx,group_stepLength=groupDist,ind_stepLength=indDist,group_stepDuration=centroidSteps,ind_stepDuration=steps,
       group_speed=groupSpeed,ind_speed=indSpeed,group_acc=groupAcc,ind_acc=indAcc,distance_to_centroid=dist2centr,
       group_heading=angles,group_direction_X=groupDirecX,group_direction_Y=groupDirecY,group_polarization=groupPolar,
       ind_heading=indAngles,ind_direction_X=indDirecX,ind_direction_Y=indDirecY,
       relative_ind_X=rotatX,relative_ind_Y=rotatY,relative_angle=relAngle,
       northX=northX,northY=northY)
}


#computing spatial discretization of a trajectory
#inputs: Two vector of coordinates (x and y) of the same length
#        step: the distance used as threshold for discretization
#       futur: boolean, whether to calculate speed, direction etc... further or backward
#       timeLine: vector of timeStamps to be used if there are several days  
#       chunkSize: number of distances which are calculated at the same time. Fiddling woth this can make the function run faster or slwoer, normaly default works quite well
#       referenceX and referenceY: optional trajectory to use as reference
#       buffer: within which range should the distance be acceptable (mainly relevant if there are NAs in the data). Range is from (step-buffer*step) to (step+buffer*step)
spatialDiscr <- function(x,y,step,futur=T,timeline=rep("1970-01-01",length(x)),chunkSize=200,referenceX=x,referenceY=y,buffer=0.25,useBuffer=T){

  totalTime <- length(x)
  dayFactor <- factor(substr(timeline,1,10))
  
  result <- tapply(1:length(x),dayFactor,function(idx){
    
    vec <- rep(NA,length(idx))
    
    end = tail(idx,1)
    
    if(futur){
      for(i in idx[-length(idx)]){
        
        #if current timepoint is NA, move on to the next one
        if(is.na(referenceX[i])){
          vec[match(i,idx)] <- NA
          next()
        }
        
        t = i + 1
        findDist = F # whether the step distance has been reached
        reachEnd = F # whether the end of the vector has been reached
        overShoot = F # whether the first distance bigger than step is too big (as decided by the buffer argument)
        
        #calculate distances until dist threshold or end of vector is reached
        while(!findDist & !reachEnd & !overShoot){
          
          chunk <- chunkSize
          if(t+chunk >= end){ #if the end has been reached
            chunk <- end - t
            reachEnd = T
          }
          
          #calculating distances with current timepoint in batches of size chunk
          dists <- sqrt((x[t:(t+chunk)]-referenceX[i])^2 + (y[t:(t+chunk)]-referenceY[i])^2)
          
          if(any(dists >= step,na.rm=T)){ #if there is at least one distance that is bigger than step
            
            reachDist =  head(which(dists >= step),1) #index of that distance in the current chunk
            
            if(useBuffer & dists[reachDist] > (step + buffer*step)){ #if the distance is too big
              
              if(any(dists < step & dists >= (step - buffer*step),na.rm=T)){ #we check if there is a smaller distance within the buffer range
                vec[match(i,idx)] =  t + head(which(dists >= (step - buffer*step)),1) - 1
                findDist = T
              }else{ #otherwise we store NA
                vec[match(i,idx)] <- NA
                overShoot <- T
              }
              
            }else if(!findDist){ #otherwise we store the index of the distance
              vec[match(i,idx)] <- t + reachDist - 1
              findDist = T
            }
         
          }else if(reachEnd){ #if we reached the end
            
            if(any(dists >= (step - buffer*step),na.rm=T)){
              vec[match(i,idx)] =  t + head(which(dists >= (step - buffer*step)),1) - 1
            }else{
              vec[match(i,idx)] <- NA
            }
          }
          
          t = t+chunk+1
        }
      }
    }else{
      
      for(i in rev(idx[-1])){
        
        if(is.na(referenceX[i])){
          vec[match(i,idx)] <- NA
          next()
        }
        
        t = i - 1
        findDist = F # whether the step distance has been reached
        reachStart = F # whether the end of the vector has been reached
        overShoot = F # whether the first distance bigger than step is too big (as decided by the buffer argument)
        
        while(!reachStart & !findDist & !overShoot){
          
          chunk <- chunkSize
          if(t-chunk <= idx[1]){
            chunk <- t - idx[1]
            reachStart = T
          }
          
          dists <- sqrt((x[t:(t-chunk)]-referenceX[i])^2 + (y[t:(t-chunk)]-referenceY[i])^2)
          
          if(any(dists >= step,na.rm=T)){
            
            reachDist =  head(which(dists >= step),1) 
            
            if(useBuffer & dists[reachDist] > (step + buffer*step)){ #if the distance is too big
              
              if(any(dists < step & dists >= (step - buffer*step),na.rm=T)){ #we check if there is a smaller distance within the buffer range
                vec[match(i,idx)] =  t - head(which(dists >= (step - buffer*step)),1) + 1
                findDist = T
              }else{ #otherwise we store NA
                vec[match(i,idx)] <- NA
                overShoot <- T
              }
              
            }else if(!findDist){ #otherwise we store the index of the distance
              vec[match(i,idx)] <- t - reachDist + 1
              findDist = T
            }
            
          }else if(reachStart){ #if we reached the start
            
            if(any(dists >= (step - buffer*step),na.rm=T)){
              vec[match(i,idx)] =  t - head(which(dists >= (step - buffer*step)),1) + 1
            }else{
              vec[match(i,idx)] <- NA
            }
          }
          t = t-chunk-1
        }
      }
    }
    vec
  })
  as.vector(unlist(result))
}

#temporal discretization
temporalDiscr<-function(x,y,step,futur=T,timeline=rep("1970-01-01",length(x))){
  
  dayFactor <- factor(substr(timeline,1,10))
  
  result <- tapply(1:length(x),dayFactor,function(idx){
    
    if(futur){
      vec <- idx + step
      vec[which(vec>tail(idx,1))] <- NA
    }else{
      vec <- idx-step
      vec[which(vec<head(idx,1))] <- NA
    }
    vec
  })
  as.vector(unlist(result))
}

#wrapper for the two above discretizing functions
discretization<-function(x,y,step,spatial=T,futur=T,timeline=rep("1970-01-01",length(x)),...){
  if(spatial){
    return(spatialDiscr(x,y,step,futur,timeline,...))
  }else{
    return(temporalDiscr(x,y,step,futur,timeline))
  }
}

#normalize a vector x in [a,b]
normalize <- function(x,a=-1,b=1,na.rm=T,minX=min(x,na.rm=na.rm),maxX=max(x,na.rm=na.rm)){
  if(length(which(is.na(x)))>=(length(x)-1)){
    return(rep(NA,length(x)))
  }else{
    return((b-a)* ((x-minX)/(maxX-minX)) + a )
  }
}

#distance between all points in a given set
allDist <- function(x,y,symetry=T,diago=F){
  distances <- outer(1:length(x),1:length(x),function(ind1,ind2)sqrt((x[ind1]-x[ind2])^2 + (y[ind1]-y[ind2])^2 ))
  if(!diago)diag(distances) <- NA
  if(!symetry)distances[lower.tri(distances)] <- NA
  distances
}

#distances between one point (argument 'focal') and a set of others
focalDistances<-function(x,y,focal=1,keepFocal=T){
  dist <- rep(NA,length(x))
  dist[-focal] <- mapply(distance,x[-focal],y[-focal],MoreArgs = list(x[focal],y[focal]))
  if(keepFocal){
    return(dist)
  }else{
    return(dist[-focal])
  }
}

# "flattened logistic" functions = our model:

flatLogis2Variables <- function(x1,x2=0,alpha,beta1,beta2=0,gamma=0.5){
  a <- alpha * (1/(1+exp(-(beta1*x1 + beta2*x2)))-0.5+gamma) + (1-alpha)*gamma
  if(any(a<=0))a[which(a<=0)] <- .Machine$double.eps
  if(any(a>=1))a[which(a>=1)] <-  1 - .Machine$double.eps
  a
}

flatLogis1Variable <- function(x,alpha,beta,gamma=0.5){
  a <- alpha * (1/(1+exp(-beta*x))-0.5+gamma) + (1-alpha)*gamma
  if(any(a<=0))a[which(a<=0)] <- .Machine$double.eps
  if(any(a>=1))a[which(a>=1)] <-  1 - .Machine$double.eps
  a
}

#log likelihood = optimizers:

logLikelihood2Variables <- function(par,data,fun,...){
  -sum( log(fun(data$x1,data$x2,par[1],par[2],par[3],...) * data$y + (1-fun(data$x1,data$x2,par[1],par[2],par[3],...)) * (1-data$y)))
}

logLikelihood1Variable <- function(par,data,fun,...){
  -sum( log(fun(data$x,par[1],par[2],...) * data$y + (1-fun(data$x,par[1],par[2],...)) * (1-data$y)))
}

#Produces ridgeLine plots
#adapted from ridgeline package by R-CoderDotCom on github
myRidgeline <- function(x, y, bw = "nrd0", mode = FALSE,col = "gray", border,meanNotMode = T,
                         palette, grouping=NULL,spacing=1,rangeX=NULL,labCex=1,lwd=2,lty=1,xlim=c(-0.25,1.25),labels=NULL,modeCol=1,...) {
  
  dens <- tapply(x, y, density, bw = bw)
  
  xs <- Map(getElement, dens, "x")
  ys <- Map(getElement, dens, "y")
  ys <- Map(function(x) (x - min(x)) / max(x - min(x)) * 1.5, ys)
  ys <- Map(`+`, ys, length(ys):1) 
  
  if(!is.null(grouping)){
    yModif <- (match(grouping,rev(unique(grouping))) -1)*spacing
    for(i in 1:length(ys)){
      ys[[i]]<-ys[[i]]+yModif[i]
    }
  }else{
    yModif <- rep(0,length(ys))
  }
  
  yMax <- max(sapply(ys,max))
  
  if(is.null(rangeX)){
    rangeX = range(xs)
  }
  
  op <- par(no.readonly = TRUE)
  par(mar = op$mar)
  plot(NULL,xlim = xlim, ylim = c(1, yMax + 1.5),...)
  abline(h = (length(ys):1)+yModif, col = col,lty = "dotted")
  
  if(missing(palette)) {
    cols <- hcl.colors(length(ys), "Zissou", alpha = 0.8)
  } else {
    if(length(palette) == 1) {
      cols <- rep(palette, length(ys))
    } else {
      cols <- palette
    }
  }
  
  if(missing(border)) {
    border <- rep(1, length(ys))
  } else {
    if(length(border) == 1) {
      border <- rep(border, length(ys))
    } else {
      border <- border
    }
  }
  
  if(isTRUE(mode)) {
    
    if(meanNotMode){
      modes <- tapply(x,y,mean)
    }else{
      modes <-  unlist(Map(which.max, ys))
    }
    
    sapply(1:length(dens), function(i) {
      polygon(xs[[i]], ys[[i]], col = cols[i], border = border[i], lty = lty, lwd = lwd)
      if(!meanNotMode){
        segments(x0 = xs[[i]][modes[i]], x1 = xs[[i]][modes[i]], y0 = (length(ys):1)[i]+yModif[i], y1 = as.numeric(Map(max, ys))[i], col = modeCol, lty = 2, lwd = 3)
      }else{
        segments(x0 = modes[i], x1 = modes[i], y0 = (length(ys):1)[i]+yModif[i], y1 = (length(ys):1)[i]+yModif[i]+1, col = modeCol, lty = 1, lwd = 3)
      }
    })
  } else {
    Map(polygon, xs, ys, col = cols, border = border, lty = lty, lwd = lwd)
  }
  
  if(!is.null(labels)) {
    if(length(labels) != length(names(dens))){
      stop("The number of labels must be equal to the number of groups")
    } else {
      mtext(labels, 2, at = (length(ys):1)+yModif, las = 2, padj = 0, line = 0.5,cex=labCex)
    }
  } else {
    mtext(names(dens), 2, at = (length(ys):1)+yModif, las = 2, padj = 0, line = 0.5,cex=labCex)
  }
  
  par(op)
}

#binning a numeric variable (y) based on categories from another numeric variable (x)
bining <- function(y=NULL,x,binSize,fun=function(f)length(which(!is.na(f))),overlap=F,count=T,range=NULL,binMean=T,quantileLimits=c(0,1),...){
  
  if(is.null(range)){
    range = quantile(x,quantileLimits,na.rm=T)
  }
  
  inRange <- which(x>=range[1] & x<=range[2])
  x <- x[inRange]
  if(length(dim(y))==2){
    y <- y[,inRange]
  }else{
    y <- y[inRange]
  }
  
  if(overlap){
    a <- numeric(0)
    bins <- (range%/%binSize)[1] : (range%/%binSize)[2]
    if(binMean){binNames=(bins*binSize + binSize/2) - i}else{binNames=bins*binSize - i}
    for(i in 0:(binSize-1)){
      
      binFactor <- factor((x+1)%/%binSize,levels=bins)
      
      if(length(dim(y))==2){
        temp <- t(apply(y,1,function(f)tapply(f,binFactor,fun,...)))
        colnames(temp) <- binNames
        a <- cbind(a,temp)
        
      }else{
        temp <- tapply(y,binFactor,fun,...)
        names(temp) <- binNames
        a <- c(a,temp)
      }
    }
    if(length(dim(y))==2){
      a <- a[,order(as.numeric(colnames(a)))]
    }else{
      a <- a[order(as.numeric(names(a)))]
    }
    
  }else{
    
    bins <- (range%/%binSize)[1] : (range%/%binSize)[2]
    binFactor <- factor(x%/%binSize,levels=bins)
    if(binMean){binNames=(bins*binSize + binSize/2)}else{binNames=bins*binSize}
    
    if(is.null(y)){
      if(count){
        a <- tapply(binFactor,binFactor,function(f)length(which(!is.na(f))),...)
      }else{
        a <- table(binFactor)/length(which(!is.na(binFactor)))
      }
      names(a) <- binNames
    }else if(length(dim(y))==2){
      a <- t(apply(y,1,function(f)tapply(f,binFactor,fun,...)))
      colnames(a) <- binNames
    }else{
      a <- tapply(y,binFactor,fun,...)
      names(a) <- binNames
    }
  }
  return(a)
}

#creating plots of one response variable as a function of a binned numerical response variable, but binned
influencePlots <- function(x,y,binSize=5,totRange = NULL,showDataSize=T,v=0,h=0.5,ylim=c(0,1),sizeLegendPos="bottomright",
                           quantileLimits = c(0.025,0,975),binFunction=mean,col1=1,lwd=2,legendPntSize = c(0.1,1,2,3),binMean=T,...){
  
  if(showDataSize){type="b"}else{type="l"}
    
  if(is.null(totRange)){
    totRange <- quantile(x, c(0.025,0.975),na.rm=T)
  }
  
  N <- bining(y,x,binSize=binSize,range=totRange,count=T)
  bins <- bining(y,x,binSize=binSize,binFunction,range=totRange,count=F,binMean=binMean,na.rm=T)
  
  pointSlope = (max(N,na.rm=T)- min(N,na.rm=T) ) / (max(legendPntSize)-min(legendPntSize))
  pointOo <- max(N,na.rm=T) - pointSlope * max(legendPntSize)
  
  plot(names(bins),bins,ylim=ylim,cex=(N-pointOo)/pointSlope,cex.lab=2,xlim=c(totRange[1]-1,totRange[2]+1),cex.axis=2,lwd=lwd,pch=16,type=type,col=col1,...)
  if(showDataSize)legend(sizeLegendPos,legend=round(legendPntSize*pointSlope + pointOo),pch=16,pt.cex=legendPntSize,bty="n",cex=2)
 
}

