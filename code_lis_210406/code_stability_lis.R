# modified from the script accompanies the paper on "Ecological Stability Properties of
# Microbial Communities Assessed by Flow Cytometry" by Liu et al., 2018
# http://msphere.asm.org/content/3/1/e00564-17

# #set working directory
setwd(".../code_lis_210406")

# 1. functions
# Function for computing resilience
computeRL <- function(dx0xt, dx0x1) {
  return(2.0*dx0x1/(dx0x1+dx0xt)-1.0)
}

# ***the following code gives the result of L1 in phase1, for other reactors and phases, codes are run repeatedly from line to with changes on data and timepoints selection 
# 2. input data (relative abundance)
abundance<-read.table("RA.txt",header=TRUE,na.strings = "NA",blank.lines.skip = FALSE)

# for phase1-Insular I
data <- subset(abundance,Reactor%in%c('L1','PRE'))[,c(3,5:84)]#or
#data <- subset(abundance,Reactor%in%c('L2','PRE'))[,c(3,5:84)]#or
#data <- subset(abundance,Reactor%in%c('L3','PRE'))[,c(3,5:84)]#or
#data <- subset(abundance,Reactor%in%c('L4','PRE'))[,c(3,5:84)]#or
#data <- subset(abundance,Reactor%in%c('L5','PRE'))[,c(3,5:84)]#or
#data <- subset(abundance,Reactor=='R')[,c(3,5:84)]

#for other phases
#data <- subset(abundance,Reactor%in%c('L1'))[,c(3,5:84)]#or
#data <- subset(abundance,Reactor%in%c('L2'))[,c(3,5:84)]#or
#data <- subset(abundance,Reactor%in%c('L3'))[,c(3,5:84)]#or
#data <- subset(abundance,Reactor%in%c('L4'))[,c(3,5:84)]#or
#data <- subset(abundance,Reactor%in%c('L5'))[,c(3,5:84)]#or
#data <- subset(abundance,Reactor=='R')[,c(3,5:84)]

#3. specify timepoints
#   the start of the experiment (experimentStart),
#   the end of the experiment (experimentEnd),
#   maximal deviation (d_max) using Euclidean distance (overrideMaxEuclidean),
#   maximal deviation (d_max) using Canberra distance (overrideMaxCanberra), and
#   the disturbance event (tref).
#
#disturbance0, for phase1-Insular I
experimentStart <- 0
experimentEnd <- 23
overrideMaxCanberra <- -1
tref <- 0

#disturbance0 for R
#experimentStart <- 9
#experimentEnd <- 23
#overrideMaxCanberra <- -1
#tref <- 9

#disturbance1, for phase2-RC10
#experimentStart <- 26
#experimentEnd <- 44
#tref <- 26

#disturbance2, for phase3-RC50
#experimentStart <- 47
#experimentEnd <- 63
#tref <- 47

#disturbance3, for phase4-RC80
#experimentStart <- 64
#experimentEnd <- 86
#tref <- 64

#disturbance4, for phase5-Insular II
#experimentStart <- 89
#experimentEnd <- 110
#tref <- 89

# To manually select gates to be plotted instead of automatically selecting the
# two most dominant gates at the start and end of the experiment, uncomment the
# following line(s) (remove the '#' character at the beginning of the line) and
# provide gate names as they appear in the column headings in the input file.
#
#plotGatesStart <- c("G1", "G2")
#plotGatesEnd <- c("G3", "G4")
#
# USER INPUT ENDS HERE ========================================================

# Characterizing the reference time point
numberOfGates <- ncol(data) - 1
data$referencePhase <- FALSE
data[data[1] >= experimentStart & data[1] <= tref, "referencePhase"] <- TRUE

referenceState <- colMeans(data[data$referencePhase==TRUE, 2:(numberOfGates + 1)])

# restricting data set to one experiment
data <- data[data[1]>=experimentStart & data[1]<=experimentEnd,]

# compute deviation from reference state
data$canberra <- apply(data[2:(numberOfGates+1)], 1, function(x) dist(rbind(referenceState, x), method = "canberra"))
data$canberra <- data$canberra / numberOfGates

#extract time reaching max deviation and corresponding state
if (overrideMaxCanberra == -1) {
  maxCanberraId <- order(data$canberra, decreasing=TRUE)[1]
} else {
  maxCanberraId <- match(overrideMaxCanberra, data[[1]])
}

maxDevStateCanberra <- data[maxCanberraId, 2:(numberOfGates+1)]

maxDevStateCanberraTime <- data[maxCanberraId, 1]

# compute online version of resilience
data$maxCanberraOnline[1] <- data$canberra[1]
for (i in 2:nrow(data)) {
  if (data$canberra[i]< data$maxCanberraOnline[i-1]) {
    data$maxCanberraOnline[i] <- data$maxCanberraOnline[i-1]
  }else {
    data$maxCanberraOnline[i] <- data$canberra[i]
  }
}

# compute Resilience RL
data$RLcanberra <- apply(data, 1, function(x) computeRL(x["canberra"], data[maxCanberraId, "canberra"]))
data$RLcanberraOnline <- apply(data, 1, function(x) computeRL(x["canberra"], x["maxCanberraOnline"]))

#conbine data for latter plot
#for the first time code running, save online data as stabilitydata
stabilitydata<-data
#or repeated code running except for the first time, add online data to stabilitydata
#stabilitydata<-rbind(stabilitydata,data)

# computing stability properties
#resistance
RScanberra <- 1.0 - data[maxCanberraId, "canberra"]
#displacement speed
DisSpeedcanberra <- data[maxCanberraId, "canberra"] / (data[maxCanberraId, 1] - tref)
#elasticity
ElasticityCanberra <- (data[maxCanberraId, "canberra"] - data[nrow(data), "canberra"]) / (data[nrow(data), 1] - data[maxCanberraId, 1])

# ensure output data as plotData 
plotData <- data[data[1] >= tref,]
plotData[1,2:(numberOfGates+1)] <- referenceState
plotData$canberra[1] <- 0
plotData$RLcanberraOnline[1] <- 0

#save all stability properties in stabilitytbl
stabilitytbl<-data.frame(Reactor=rep(c('L1','L2','L3','L4','L5','R'),each=4),Disturbance=rep(c(1:4),6),RS=rep(NA,24),DS=rep(NA,24),RL=rep(NA,24),E=rep(NA,24))
stabilitytbl[1,3:6]<-c(RScanberra,DisSpeedcanberra,plotData$RLcanberra[nrow(plotData)],ElasticityCanberra)#or
#stabilitytbl[2,3:6]<-c(RScanberra,DisSpeedcanberra,plotData$RLcanberra[nrow(plotData)],ElasticityCanberra)#or
#stabilitytbl[3,3:6]<-c(RScanberra,DisSpeedcanberra,plotData$RLcanberra[nrow(plotData)],ElasticityCanberra)#or
#stabilitytbl[4,3:6]<-c(RScanberra,DisSpeedcanberra,plotData$RLcanberra[nrow(plotData)],ElasticityCanberra)#or
#stabilitytbl[5,3:6]<-c(RScanberra,DisSpeedcanberra,plotData$RLcanberra[nrow(plotData)],ElasticityCanberra)#or
#stabilitytbl[6,3:6]<-c(RScanberra,DisSpeedcanberra,plotData$RLcanberra[nrow(plotData)],ElasticityCanberra)#or
#stabilitytbl[7,3:6]<-c(RScanberra,DisSpeedcanberra,plotData$RLcanberra[nrow(plotData)],ElasticityCanberra)#or
#stabilitytbl[8,3:6]<-c(RScanberra,DisSpeedcanberra,plotData$RLcanberra[nrow(plotData)],ElasticityCanberra)#or
#stabilitytbl[9,3:6]<-c(RScanberra,DisSpeedcanberra,plotData$RLcanberra[nrow(plotData)],ElasticityCanberra)#or
#stabilitytbl[10,3:6]<-c(RScanberra,DisSpeedcanberra,plotData$RLcanberra[nrow(plotData)],ElasticityCanberra)#or
#stabilitytbl[11,3:6]<-c(RScanberra,DisSpeedcanberra,plotData$RLcanberra[nrow(plotData)],ElasticityCanberra)#or
#stabilitytbl[12,3:6]<-c(RScanberra,DisSpeedcanberra,plotData$RLcanberra[nrow(plotData)],ElasticityCanberra)#or
#stabilitytbl[13,3:6]<-c(RScanberra,DisSpeedcanberra,plotData$RLcanberra[nrow(plotData)],ElasticityCanberra)#or
#stabilitytbl[14,3:6]<-c(RScanberra,DisSpeedcanberra,plotData$RLcanberra[nrow(plotData)],ElasticityCanberra)#or
#stabilitytbl[15,3:6]<-c(RScanberra,DisSpeedcanberra,plotData$RLcanberra[nrow(plotData)],ElasticityCanberra)#or
#stabilitytbl[16,3:6]<-c(RScanberra,DisSpeedcanberra,plotData$RLcanberra[nrow(plotData)],ElasticityCanberra)#or
#stabilitytbl[17,3:6]<-c(RScanberra,DisSpeedcanberra,plotData$RLcanberra[nrow(plotData)],ElasticityCanberra)#or
#stabilitytbl[18,3:6]<-c(RScanberra,DisSpeedcanberra,plotData$RLcanberra[nrow(plotData)],ElasticityCanberra)#or
#stabilitytbl[19,3:6]<-c(RScanberra,DisSpeedcanberra,plotData$RLcanberra[nrow(plotData)],ElasticityCanberra)#or
#stabilitytbl[20,3:6]<-c(RScanberra,DisSpeedcanberra,plotData$RLcanberra[nrow(plotData)],ElasticityCanberra)#or
#stabilitytbl[21,3:6]<-c(RScanberra,DisSpeedcanberra,plotData$RLcanberra[nrow(plotData)],ElasticityCanberra)#or
#stabilitytbl[22,3:6]<-c(RScanberra,DisSpeedcanberra,plotData$RLcanberra[nrow(plotData)],ElasticityCanberra)#or
#stabilitytbl[23,3:6]<-c(RScanberra,DisSpeedcanberra,plotData$RLcanberra[nrow(plotData)],ElasticityCanberra)#or
#stabilitytbl[24,3:6]<-c(RScanberra,DisSpeedcanberra,plotData$RLcanberra[nrow(plotData)],ElasticityCanberra)

#organise stabilitydata
P<-c('Closing-I','Recycle-10%','Recycle-50%','Recycle-80%','Closing-II')
s<-c(P[1],as.character(subset(abundance,Time%in%c(1:110)&Reactor=='L1')[,4]))
s_r<-as.character(subset(abundance,Reactor=='R')[,4])
stabilitydata<-data.frame(Reactor=c(rep(reactorname[1:5],each=77),rep('R',68)),Phase=c(rep(s,5),s_r),stabilitydata)
stabilitydata<-stabilitydata[,-c(5:84)]

#output 
write.csv(stabilitytbl,'stabilitytbl.csv',row.names = FALSE)
write.csv(stabilitydata,'stabilitydata.csv',row.names = FALSE)