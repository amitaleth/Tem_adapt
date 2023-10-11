t <- NULL

for(i in 1:71){ t[i]<-i-1 }
gask <- function(x,m,lo,hi,sup,sdn,b)
  
{
  
  if(lo<0 | hi<=lo | sup<=0 | sdn<=0 | b<=0){ return(-1e-100)}
  
  ss <- sup+(sdn-sup)/(1+exp(-(x-m)/b))
  
  xx <- (x-m)/ss
  
  lo + (hi-lo)*exp(-xx*xx)
  
}



raw <- read.table("C:/Users/17340516030/Desktop/out.168-hermo.tab",sep="\t",header=FALSE)



ymax <- max(raw[,2])



y <- gask(t,raw[2,5],raw[2,1],raw[2,2],raw[2,3],raw[2,4],5)



plot(t,y,type="l",xlim=c(0,70),ylim=c(0,ymax),col="grey",xaxt="n")
axis(1,at=c(0,10,20,30,40,50,60,70),las=2)



for(i in 3:1001){
  
  y <- gask(t,raw[i,5],raw[i,1],raw[i,2],raw[i,3],raw[i,4],5)
  
  lines(t,y,type="l",col="grey")
  
}

y <- gask(t,raw[1,5],raw[1,1],raw[1,2],raw[1,3],raw[1,4],5)

lines(t,y,type="l",col="red",lwd=3)
