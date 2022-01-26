library(nlme)
library(multcomp)

setwd("C:/Users/baverly/Desktop/Meerkat_leadership_hackathon/INFLUENCE_PAPER")
source("scripts/functions.R")

discretizationStep <- 10

load(paste0("output/spatialMetrics_",discretizationStep,"m.RData"))
load(paste0("output/Influence_logistic_model_fits_",discretizationStep,"m.RData"))

sessions <- c("HM2017","HM2019","L2019","ZU2021","NQ2021")
status = c("DominantF","DominantM","Adult","Yearling","Sub-Adult")

allInd = modelParam_MovTurn$ind

totPropFront <- tapply(spatialMetrics$inFrontHalf,factor(spatialMetrics$indUniqID,levels=allInd),mean,na.rm=T)


#----MOVEMENT TURNING INFLUENCE----

data <- modelParam_MovTurn

data$status <- factor(data$status,levels=status)

mod <- lme(inflScore~status,random=~1|session,data=data)
anova(mod)
qqnorm(residuals(mod))
qqline(residuals(mod))

#comparaison des status par paires
pairWise <- glht(mod, linfct = mcp(status = "Tukey"))
summary(pairWise)


#----MOVEMENT SPEEDING INFLUENCE----

data <- modelParam_MovSpeed

data$status <- factor(data$status,levels=status)

mod <- lme(inflScore~status,random=~1|session,data=data)
anova(mod)
qqnorm(residuals(mod))
qqline(residuals(mod))

#comparaison des status par paires
pairWise <- glht(mod, linfct = mcp(status = "Tukey"))
summary(pairWise)
