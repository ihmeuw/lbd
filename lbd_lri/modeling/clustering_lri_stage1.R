setwd('<<<< FILEPATH REDACTED >>>>')
rundate <- '<<<< RUN DATE REDACTED >>>>'
data1bydraw <- readRDS('<<<< FILEPATH REDACTED >>>>')
data1 <- read.csv('<<<< FILEPATH REDACTED >>>>')
NumDraws <- 250


SetToZerobydraw <- which(is.na(data1bydraw$mean_rate))
data1bydraw$mean_rate[SetToZerobydraw] <- 0
data1bydraw[SetToZerobydraw,1:(2*NumDraws)+12] <- 0
keep1bydraw <- which(data1bydraw$measure == "mortality")
data1bydraw <- data1bydraw[keep1bydraw,]



SetToZero <- which(is.na(data1$mean_rate))
data1$mean_rate[SetToZero] <- 0
keep1 <- which(data1$measure == "mortality")
data1 <- data1[keep1,]

data1 <- data1[-c(13295:13304),]
data1bydraw <- data1bydraw[-c(13295:13304),]

if (!require("pacman")) install.packages("pacman")
pacman::p_load(maptools,
               RColorBrewer,
               raster,
               animation,
               spdep)


mask <- raster('<<<< FILEPATH REDACTED >>>>')
A1 <- readRDS('<<<< FILEPATH REDACTED >>>>')
A0 <- readRDS('<<<< FILEPATH REDACTED >>>>')

UA1 <- unique(data1$ADM1_CODE)
UA0 <- unique(data1$ADM0_CODE)

Keep0 <- sapply(UA0, function(x) which(A0$ADM0_CODE == x))

A0 <- A0[Keep0,]

Keep <- sapply(UA1, function(x) which(A1$ADM1_CODE == x))

A1 <- A1[Keep,]

A1nb <- poly2nb(A1,snap=.007)
listA1 <- nb2listw(A1nb, style="B",zero.policy = TRUE)



Big <- 0
for (Year in c(2000,2005,2010,2016)){
  datasub <- data1[which(data1$year == Year),]
  NewOrder <- sapply(A1$ADM1_CODE, function(x) which(datasub$ADM1_CODE == x)[1])
  datasub <- datasub[NewOrder,]
  test2 <- localG(datasub$mean_rate, listA1,zero.policy = TRUE)
  Big <- max(Big,round(max(abs(range(test2,na.rm=TRUE)))))
}
COLS <- c(rev(brewer.pal(9,"Blues")[-c(6:7)]),"white",brewer.pal(9,"Reds")[-c(6:7)])

Big <- 5.25
BINS <- seq(-Big,Big,length=length(COLS)+1)

jpeg(file='<<<< FILEPATH REDACTED >>>>', width=1600,height=1600)
layout(matrix(c(1,2,5,3,4,5),2,3,byrow=TRUE),width=c(5,5,1))
par(mar=c(0,0,2.1,0),oma=c(0,0,0,3.1))
for (Year in c(2000,2005,2010,2016)){
  datasub <- data1[which(data1$year == Year),]
  NewOrder <- sapply(A1$ADM1_CODE, function(x) which(datasub$ADM1_CODE == x)[1])
  datasub <- datasub[NewOrder,]
  test <- globalG.test(datasub$mean_rate, listA1,zero.policy = TRUE)
  test2 <- localG(datasub$mean_rate, listA1,zero.policy = TRUE)
  ColLoc <- sapply(test2,function(x) findInterval(x,BINS,all.inside=TRUE))
  plot(A1,col=COLS[ColLoc],ylim=c(-33,36),border="grey70")
  plot(A0,add=TRUE,lwd=2.5)
  mtext(Year,3,line=-5,cex=2)
}

plot(1,type="n",ann=FALSE,axes=FALSE,ylim=c(0,length(COLS)),xlim=c(-.5,1.5))
for (i in 1:length(COLS)){
  rect(0,i-1,.5,i,col=COLS[i])
  text(.5,i-1,round(BINS[i],1),pos=4,cex=1.75)
}
text(.5,length(COLS),round(BINS[length(COLS)+1]),pos=4,cex=1.75)
mtext("Local G statistic",4,line=-1,cex=2)
dev.off()


Cut <- abs(qnorm(1-(1-0.05)^(1/781)))

tvec <- 2000:2016
GStat <- matrix(0,17,3)
LocalG <- vector("list",17)
for (i in 1:17){
  Year <- i+1999
  datasub <- data1[which(data1$year == Year),]
  NewOrder <- sapply(A1$ADM1_CODE, function(x) which(datasub$ADM1_CODE == x)[1])
  datasub <- datasub[NewOrder,]
  test <- globalG.test(datasub$mean_rate, listA1,zero.policy = TRUE)
  LocalG[[i]] <- localG(datasub$mean_rate, listA1,zero.policy = TRUE)
  GStat[i,] <- test$estimate
}

Year <- 2000
datasub <- data1[which(data1$year == Year),]
NewOrder <- sapply(A1$ADM1_CODE, function(x) which(datasub$ADM1_CODE == x))
datasub <- datasub[NewOrder,]
Hot2000 <- which(LocalG[[1]] > Cut)
Y2000 <- table(datasub$ADM0_NAME[Hot2000])
Y2000 <- Y2000[Y2000>0]
Tot2000 <- sapply(names(Y2000),function(x)length(which(datasub$ADM0_NAME == x)))

TotDeath2000 <- sum(datasub$mean_count)
TotPop2000 <- sum(datasub$pop)

HotDeath2000 <- sum(datasub$mean_count[Hot2000])
HotPop2000 <- sum(datasub$pop[Hot2000])


Year <- 2016
datasub <- data1[which(data1$year == Year),]
NewOrder <- sapply(A1$ADM1_CODE, function(x) which(datasub$ADM1_CODE == x))
datasub <- datasub[NewOrder,]
Hot2016 <- which(LocalG[[16]] > Cut)
Y2016 <- table(datasub$ADM0_NAME[Hot2016])
Y2016 <- Y2016[Y2016>0]
Tot2016 <- sapply(names(Y2016),function(x)length(which(datasub$ADM0_NAME == x)))


TotDeath2016 <- sum(datasub$mean_count)
TotPop2016 <- sum(datasub$pop)

HotDeath2016 <- sum(datasub$mean_count[Hot2016])
HotPop2016 <- sum(datasub$pop[Hot2016])

Upper <- GStat[,1] + 1.96 * sqrt(GStat[,3])
Lower <- GStat[,1] - 1.96 * sqrt(GStat[,3])

pdf(file='<<<< FILEPATH REDACTED >>>>', width=8,height=8)

plot(tvec,GStat[,1],ylim=range(c(Lower,Upper)),pch=19,ylab="Global G test statistic",xlab="Year")
for (i in 1:17){
  Year <- i+1999
  arrows(Year,GStat[i,1],Year,Lower[i],angle=90,length=.1)
  
  arrows(Year,GStat[i,1],Year,Upper[i],angle=90,length=.1)
}

dev.off()

pdf(file='<<<< FILEPATH REDACTED >>>>', width=8,height=8)

mod <- lm(GStat[,1]~ tvec, weights = 1/GStat[,3])
plot(tvec,GStat[,1],ylim=range(c(Lower,Upper)),pch=19,ylab="Global G test statistic",xlab="Year")
for (i in 1:17){
  Year <- i+1999
  arrows(Year,GStat[i,1],Year,Lower[i],angle=90,length=.1)
  arrows(Year,GStat[i,1],Year,Upper[i],angle=90,length=.1)
}

abline(mod,lwd=2,lty=2)
dev.off()

#### By draw

Cut <- abs(qnorm(1-(1-0.05)^(1/781)))


Years <- 2000:2016
NLoc1 <- length(which(data1$year == 2000))


Out1 <- matrix(0, length(Years), NLoc1)


for (i in 1:length(Years)){
  datathisyear <- data1bydraw[which(data1bydraw$year == 1999+i),]
  for (j in 1:NumDraws){
    datasub <- datathisyear[,c(2,16+j)]
    NewOrder <- sapply(A1$ADM1_CODE, function(x) which(datasub$ADM1_CODE == x))
    datasub <- datasub[NewOrder,]
    tmplocal <- localG(datasub[,2], listA1,zero.policy = TRUE)
    sig <- which(tmplocal > Cut)
    Out1[i,sig] <- Out1[i,sig] + 1
  }
}


colfunc <- colorRampPalette(c("white","red"))
COLS <- colfunc(NumDraws+1)
COLS[floor(NumDraws*.95):(NumDraws+1)] <- "purple"

par(mfrow=c(2,2),mar=rep(0,4))
for (i in c(1,6,11,17)){
  plot(A1,col=COLS[Out1[i,]+1])
}


jpeg(file='<<<< FILEPATH REDACTED >>>>', width=1600,height=1600)
layout(matrix(c(1,2,5,3,4,5),2,3,byrow=TRUE),width=c(5,5,1))
par(mar=c(0,0,2.1,0),oma=c(0,0,0,3.1))
for (i in c(1,6,11,17)){
  Year <- 1999+i
  plot(A1,col=COLS[Out1[i,]+1],ylim=c(-33,36),border="grey70")
  plot(A0,add=TRUE,lwd=2.5)
  mtext(Year,3,line=-5,cex=2)
}

plot(1,type="n",ann=FALSE,axes=FALSE,ylim=c(0,length(COLS)),xlim=c(-.5,1.5))

for (i in 1:length(COLS)){
  rect(0,i-1,.5,i,col=COLS[i],border=NA)
}
PlotVec <- pretty(c(0,NumDraws))
for (i in PlotVec[-length(PlotVec)]){
  text(.5,i,100*i/NumDraws,pos=4,cex=1.75)
}
rect(0,floor(NumDraws*.95),.5,NumDraws+1,col=COLS[NumDraws+1],border=NA)
text(.5,length(COLS)*.975,">95",pos=4,cex=1.75)

mtext("Percent of draws identified
as hotspot",4,line=-1,cex=2)

dev.off()


