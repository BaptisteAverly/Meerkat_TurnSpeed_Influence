setwd("C:/Users/baverly/Desktop/Meerkat_leadership_hackathon/INFLUENCE_PAPER")

source("scripts/functions.R")
sessions <- c("HM2017","HM2019","L2019","ZU2021","NQ2021")

#scan data
scanFiles <- list.files('data/scans',pattern = c("SCAN_WHOLE"),recursive=T)
synchInfo <- read.csv("data/metadata/total_synch_info.csv")

for(session in sessions){
  
  #loading movement data
  longNames <- c(load(paste0("data/movement/level0/",session,"_COORDINATES_SYNCHED_level0.RData")))
  shortNames <- simplifyNames(pat=paste(session,"_",sep=""))

  row.names(allX) <- row.names(allY) <- indInfo$code
  
  #----KEEPING ONLY MOMENTS WHEN ENOUGH INDIVIDUALS ARE RECORDED----
  
  nonJuv <- which(indInfo$status!="Juvenile")
  keepThresh <- 2/3
  
  # removing data if (number of non-juveniles recorded / number of non-juveniles present that day) < keepThresh
  lowRecord <- as.numeric(unlist(tapply(1:length(timeLine),substr(timeLine,1,10),function(idx){
    day <- gsub("-","",substr(head(timeLine[idx],1),1,10))
    nonJuvPresent = length(which(movSummary[nonJuv,match(day,colnames(movSummary))]!="Absent"))
    
    nbrInd <- apply(allX[nonJuv,idx],2,function(f)length(which(complete.cases(f))))
    idx[which(nbrInd/nonJuvPresent < keepThresh)]
  })))
  
  allX[,lowRecord] <- NA
  allY[,lowRecord] <- NA
  
  #removing days when all data were removed
  goodDays <- which(tapply(1:length(timeLine),substr(timeLine,1,10),function(idx){
    !all(is.na(allX[,idx]))
  }))
  
  keep <- as.vector(sapply(goodDays,function(f)dayIdx[f]:(dayIdx[f+1]-1)))
  allX <- allX[,keep]
  allY <- allY[,keep]
  timeLine <- timeLine[keep]
  dayIdx <- match(unique(substr(timeLine,1,10)),substr(timeLine,1,10))
  dayIdx <- c(dayIdx,length(timeLine)+1)
  dates <- dates[goodDays]
  
  
  #----REMOVING SPECIAL EVENTS----
  
  if(session=="HM2017"){
    
    #focal had not started
    allX[1,dayIndex(which(dates=="20170903"))][1:900] <- NA
    allY[1,dayIndex(which(dates=="20170903"))][1:900] <- NA
    
  }else if(session=="HM2019"){
    
    #rovers
    allX[c(3,6),(dayIndex(which(dates=="20190625")))[5750:8000]] <- NA
    allY[c(3,6),(dayIndex(which(dates=="20190625")))[5750:8000]] <- NA
    
    absentDante <- c((dayIndex(which(dates=="20190627")))[8520:10260],
                     dayIndex(which(dates=="20190629")),
                     (dayIndex(which(dates=="20190718")))[1:1600])
    allX[3,absentDante] <- NA
    allY[3,absentDante] <- NA
    
    #removing the day where the DM and 2 other males went roving
    noKeep <- dayIndex(which(dates=="20190712"))
    allX <- allX[,-noKeep]
    allY <- allY[,-noKeep]
    timeLine <- timeLine[-noKeep]
    dayIdx <- match(unique(substr(timeLine,1,10)),substr(timeLine,1,10))
    dayIdx <- c(dayIdx,length(timeLine)+1)
    dates <- dates[-which(dates=="20190712")]
    
  }else if(session=="L2019"){
    
    #Day they encountered HM
    allX[3,dayIndex(which(dates=="20190810"))][1:1080] <- NA
    allY[,dayIndex(which(dates=="20190810"))][,4200:10800] <- NA
    allX[,dayIndex(which(dates=="20190810"))][,4200:10800] <- NA
    allY[,dayIndex(which(dates=="20190810"))][,4200:10800] <- NA
    
  }else if(session=="ZU2021"){
    
    #crashing times
    allX[,dayIndex(which(dates=="20210519"))][,8100:10800] <- NA
    allY[,dayIndex(which(dates=="20210519"))][,8100:10800] <- NA
    allX[,dayIndex(which(dates=="20210520"))][,7800:10800] <- NA
    allY[,dayIndex(which(dates=="20210520"))][,7800:10800] <- NA
    
  }else if(session=="NQ2021"){
    
    #crashing times
    allX[,dayIndex(which(dates=="20210811"))][,6000:8000] <- NA
    allY[,dayIndex(which(dates=="20210811"))][,6000:8000] <- NA
  }
  
  #----REMOVING ALARM EVENTS FROM SCAN DATA----
  
  alarmEvents <- rep(NA,ncol(allX))
  
  for (file in scanFiles){
    date <- substr(tail(strsplit(file,"_")[[1]],1),1,8)
    if(date %in% dates){
      
      scans <- read.csv(paste0("data/scans/",file))
      if(substr(as.character(scans$Time[1]),2,2)=="/"){
        scans$Time <- paste("0",scans$Time,sep="")
        scans$Time <- as.POSIXlt(as.character(scans$Time),format="%m/%d/%y %T",tz="")
      }else{
        scans$Time <- as.POSIXlt(as.character(scans$Time),format="%d/%m/%Y %T",tz="")
      }
      
      #computing the difference between tablet time and gps time
      synchDiff<- as.POSIXct(synchInfo$Device.Time.SAST[match(date,synchInfo$Date)],tz="UTC") - as.POSIXct(synchInfo$GPS.Time.UTC[match(date,synchInfo$Date)],tz="UTC")
      
      scans$GPS_time_UTC <- scans$Time - synchDiff
      scans$GPS_time_UTC <- as.POSIXct(as.character(scans$GPS_time_UTC),tz="UTC")
      
      alarms <- scans[grep("Alarm running",scans$Event.type),]
      
      if(nrow(alarms)!=0 & !all(is.na(alarms$Time))){
        
        for(l in grep("start",alarms$Timing,ignore.case = T)){
          if(grepl("stop",alarms$Timing[l+1])){
            alarmTL <- as.POSIXlt((alarms$Time[l] - 20) : (alarms$Time[l+1]+40),origin="1970-01-01 00:00.00",tz="UTC")
          }else{
            alarmTL <- as.POSIXlt((alarms$Time[l] - 20) : (alarms$Time[l]+40),origin="1970-01-01 00:00.00",tz="UTC")
          }
          alarmEvents[match(alarmTL,timeLine)] <- "alarm"
        }
        
        for(l in grep("ponctual",alarms$Timing,ignore.case = T)){
          
          alarmTL <- as.POSIXlt((alarms$Time[l] - 20) : (alarms$Time[l]+40),origin="1970-01-01 00:00.00",tz="UTC")
          alarmEvents[match(alarmTL,timeLine)] <- "alarm"
        }
      }
    }
  }
  
  allX[,which(!is.na(alarmEvents))] <- NA
  allY[,which(!is.na(alarmEvents))] <- NA
  
  obj <- simplifyNames(paste0(session,"_"),reverse= T,items = shortNames)
  
  save(list=obj,file=paste0("data/movement/level1/",session,"_COORDINATES_SYNCHED_level1.RData"))
}



