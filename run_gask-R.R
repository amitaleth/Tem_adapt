# arguments
alst <- commandArgs(trailingOnly=TRUE)
fdat <- as.character(alst[1])
fout <- as.character(alst[2])
nboo <- as.numeric(alst[3])

# DEFINE gask
gask <- function(x,m,lo,hi,sup,sdn,b)
{
 if(lo<0 | hi<=lo | sup<=0 | sdn<=0 | b<=0){ return(-1e-100)}
 ss <- sup+(sdn-sup)/(1+exp(-(x-m)/b))
 xx <- (x-m)/ss
 lo + (hi-lo)*exp(-xx*xx)
}

# DEFINE err4
err4 <- function(par,mm,ll,xxx,t,ntem,nrep)
{
 yyy <- t
 for(i in 1:ntem){ yyy[i]<-gask(t[i],mm,par[1],par[2],par[3],par[4],ll)}
 ddd <- xxx
 for(i in 1:nrep){ ddd[,i] <- xxx[,i] - yyy}
 sum(ddd^2)
}

# DEFINE err5
err5 <- function(par,ll,xxx,t,ntem,nrep)
{
 yyy <- t
 for(i in 1:ntem){ yyy[i]<-gask(t[i],par[5],par[1],par[2],par[3],par[4],ll)}
 ddd <- xxx
 for(i in 1:nrep){ ddd[,i] <- xxx[,i] - yyy}
 sum(ddd^2)
}

# data
raw <- read.table(fdat,sep="\t",header=FALSE)
raw.dim <- dim(raw)
ntem <- raw.dim[1]
nrep <- raw.dim[2]-1

t <- as.numeric(raw[,1])
dat <- raw[,2:(nrep+1)]

#dat <- read.table(fdat,sep="\t",header=FALSE)
#t<-seq(from=0,length.out=ntem,by=10)

# scale
# just copy, no need to scale
dax <- dat
#for(i in 1:nrep){
# lf<-lsfit(y=dat[,2],x=dat[,i])
# dax[,i] = dat[,i]*lf$coef[2]+lf$coef[1]
#}

# medians
me <- t
for(i in 1:ntem){ me[i]<-median(as.numeric(dax[i,]))}

# estimates
fmin <- min(me)
fmax <- max(me)
mi <- which(me==fmax)
mt <- t[mi]
fthr <- fmin + (fmax-fmin)/5

sup <- 5
i <- mi - 1;
while(i>=1){
 sup <- mt - t[i]
 if(me[i]<=fthr){
  sup <- (mt - t[i])/2
  break
 }
 i <- i - 1
}

sdn <- 5
i <- mi + 1;
while(i<=ntem){
 sdn <- t[i] - mt
 if(me[i]<=fthr){
  sdn <- (t[i] - mt)/2
  break
 }
 i <- i + 1
}

# optimize 4 parameters
para = c(fmin,fmax,sup,sdn)
m4 <- optim(par=para,err4,mm=mt,ll=5,xxx=dax,t=t,ntem=ntem,nrep=nrep,control=list(maxit=10000))

# optimize 5 parameters
parx <- c(m4$par,mt)
m5 <- optim(par=parx,err5,ll=5,xxx=dax,t=t,ntem=ntem,nrep=nrep,control=list(maxit=10000))

# output
nout <- length(m5$par)
write(t(m5$par),file=fout,sep="\t",ncol=nout)

# resample
day <- dax
for(ii in 1:nboo){
 for(i in 1:ntem){ day[i,] <- sample(dax[i,],replace=TRUE)}
 for(i in 1:ntem){ me[i]<-median(as.numeric(dax[i,]))}
 mf <- max(me)
 mi <- which(me==mf)
 mt <- t[mi]
 m4 <- optim(par=para,err4,mm=mt,ll=5,xxx=day,t=t,ntem=ntem,nrep=nrep,control=list(maxit=10000))
 parx <- c(m4$par,mt)
 m5 <- optim(par=parx,err5,ll=5,xxx=day,t=t,ntem=ntem,nrep=nrep,control=list(maxit=10000))
 nout <- length(m5$par)
 write(t(m5$par),file=fout,sep="\t",ncol=nout,append=TRUE)
}
