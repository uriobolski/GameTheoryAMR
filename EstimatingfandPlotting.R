library(tidyverse)
library(ggplot2)
library(truncdist)
library(ReIns)
library(nimble)
library(grid)
library(ggthemes)
#some ggplot commands for aesthetic plots
source('./theme_Publication.R')
#### cauculate and plot T, E, PH, PL for the theoretical distributions ####

#functions for numeric integration of densities
integratefrom=function(lowlim,myfunc){
  return(integrate(myfunc,lower = lowlim,upper=0.9999)[[1]])
}

integrateto=function(uplim,myfunc){
  return(integrate(myfunc,lower = 0.0001,upper=uplim,stop.on.error = F)[[1]])
}



## Main text examples (figure 2)
#run a loop on combinations of beta parameters
combs=matrix(ncol=2,c(0.5,0.5,20,20,2,5,1,1),byrow = TRUE)
for (i in 1:dim(combs)[1]){
  alpha=combs[i,1]
  beta=combs[i,2]
  xvals=seq(0.0001,0.9999,length.out = 1000)
  fpbeta=dbeta(xvals,shape1 =alpha ,shape2 =beta )
  betafun=data.frame(cbind(xvals,fpbeta))
#plot top panels
    ggplot(betafun,aes(x=xvals,y=fpbeta))+
    geom_ribbon(data=betafun,aes(ymax=fpbeta),fill="red", ymin=0)+
    theme_Publication()
  ggsave(paste0('alpha=',alpha,', beta=',beta,'dist.png'),width=7,height=3.5)
  
  #approximate functions for convenience in integration later on
  fphat=approxfun(x=xvals,y=fpbeta)
  fptimesphat=approxfun(xvals,fpbeta*xvals)
  
  #calculate expected posterior
  E=integrate(fptimesphat,lower = 0.0001,upper=0.9999)[[1]]
  #T values to be checked
  Thresh=seq(0.0001,0.9999,0.0001)
  
  #high and low posterior probabilities 
  pH=sapply(Thresh,integratefrom,fptimesphat)/sapply(Thresh,integratefrom,fphat)
  pL=sapply(Thresh,integrateto,fptimesphat)/sapply(Thresh,integrateto,fphat)
  
  Finaldat=data.frame(cbind(pH,pL,Thresh))
  #find the first T value satisfying the condition
  intersection=Thresh[which(pH-pL-E>0)[1]]
  #plot bottom panels 
  ggplot(Finaldat,aes(x=Thresh))+geom_line(aes(y=pH),size=1.5,linetype = "dashed")+
    geom_line(aes(y=pL),size=1.5,linetype = "dashed")+ geom_line(aes(y=pH-pL),size=2,col='red')+
    geom_hline(yintercept=E,size=2)+
    geom_segment(x=intersection,y=0,xend=intersection,yend=E,size=2,color='blue')+
    ylab('') +theme_Publication()
  ggsave(paste0('alpha=',alpha,', beta=',beta,'f.png'),width=7,height=3.5)

}

##supplementary examples

#we keep the 'fpbeta' name vec for convenience and choose
#from the commented distributions 
# fpbeta=dtexp(xvals,rate =1/10 ,endpoint = 1)
# fpbeta=dtrunc(xvals,'norm',a=0,b=1,mean=0.5,sd=1)
# fpbeta=dtrunc(xvals,'norm',a=0,b=1,mean=0.2,sd=0.05)
# fpbeta=dtrunc(xvals,'dexp',a=0,b=1,location=0.4,scale=0.4)
fpbeta=dtrunc(xvals,'exp',a=0,b=1,rate=5)

#the rest of the code below is as in the code for
#the beta distributions above
fphat=approxfun(x=xvals,y=fpbeta)
fptimesphat=approxfun(xvals,fpbeta*xvals)

betafun=data.frame(cbind(xvals,fpbeta))

#top panels
ggplot(betafun,aes(x=xvals,y=fpbeta))+
  geom_ribbon(data=betafun,aes(ymax=fpbeta),fill="red", ymin=0)+
  # ylim(c(0,max(fpbeta)))+
  theme_Publication()
ggsave("SIplottop.png",width=7,height=3.5)


E=integrate(fptimesphat,lower = 0.0001,upper=0.9999)[[1]]

Thresh=seq(0.0001,0.9999,0.0001)

pH=sapply(Thresh,integratefrom,fptimesphat)/sapply(Thresh,integratefrom,fphat)
pL=sapply(Thresh,integrateto,fptimesphat)/sapply(Thresh,integrateto,fphat)

Finaldat=data.frame(cbind(pH,pL,Thresh))
intersection=Thresh[which(pH-pL-E>0)[1]]
#bottom panels
ggplot(Finaldat,aes(x=Thresh))+geom_line(aes(y=pH),size=1.5,linetype = "dashed")+
  geom_line(aes(y=pL),size=1.5,linetype = "dashed")+ geom_line(aes(y=pH-pL),size=2,col='red')+
  geom_hline(yintercept=E,size=2)+
  geom_segment(x=intersection,y=0,xend=intersection,yend=E,size=2,color='blue')+
  ylab('') +theme_Publication()
ggsave("SIplotbottom.png",width=7,height=3.5)



#### plot an empirical f ####
##plotting the empirical example

##load estimated probabilities of a bacterial infection
## the probabilities are obtained as explained in the Methods
pdataxgbost=read.csv("fhat.csv",header=FALSE)
pdataxgbost=pdataxgbost[[1]]

## create a smooth density function
fpxgb=density(as.numeric(na.omit(pdataxgbost)),
              kernel ='gaussian',
              from = 0,to = 1)

#fix fp to have a 0-1 support using ~1/(2n)
fpxgb$y=fpxgb$y+1/(2500)
#renormalize fp
fpxgb$y=fpxgb$y/sum(fpxgb$y*(fpxgb$x[2]-fpxgb$x[1]))
RSVfpxgb=data.frame(cbind(fpxgb$x,fpxgb$y))
labelpdat=data.frame(cbind(pdataxgbost,as.character(datrsv$xray_diag)))
names(labelpdat)=c('f','diagnosis')
labelpdat$f=as.numeric(as.character(labelpdat$f))
names(RSVfpxgb)=c('x','y')

#create a density plot+rug plt for predictions (figure 3)
ggplot(RSVfpxgb,aes(x=x,y=y))+
  geom_ribbon(data=RSVfpxgb,aes(ymax=y),fill="red", ymin=0)+
  geom_rug(data=labelpdat,aes(y=-0.5,x=f,color=diagnosis),alpha=0.5,size=1.5,length = unit(0.1, "npc"))+
  theme_Publication()
ggsave("empiricdist.png")

#### cauculate and plot T, E, PH, PL for the empiric example ####

#create functions from the smooothed density
fphat=approxfun(fpxgb)
fptimesphat=approxfun(fpxgb$x,fpxgb$y*fpxgb$x)

#the rest as in the theoretical examples
E=integrate(fptimesphat,lower = 0.001,upper=0.999)[[1]]

pH=sapply(Thresh,integratefrom,fptimesphat)/sapply(Thresh,integratefrom,fphat)
pL=sapply(Thresh,integrateto,fptimesphat)/sapply(Thresh,integrateto,fphat)

Finaldat=data.frame(cbind(pH,pL,Thresh))
intersection=Thresh[which(pH-pL-E>0)[1]]

ggplot(Finaldat,aes(x=Thresh))+geom_line(aes(y=pH),size=1.5,linetype = "dashed")+
  geom_line(aes(y=pL),size=1.5,linetype = "dashed")+ geom_line(aes(y=pH-pL),size=2,col='red')+
  geom_hline(yintercept=E,size=2)+
  geom_segment(x=intersection,y=0,xend=intersection,yend=E,size=2,color='blue')+
  ylab('') +theme_Publication()
ggsave("empiricf.png",width=7,height=3.5)
