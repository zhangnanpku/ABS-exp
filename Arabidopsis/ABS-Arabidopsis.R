### Arabidopsis gene, two generations, sub-region ###
rm(list=ls())

dat=read.table('data_generation.txt', header=T)
n=dim(dat)[1]
attach(dat)

library(gss)

#########################################
my.x <- position
my.z <- generation

my.N <- h
my.y <- mc
my.p <- my.y/my.N

## slicing 
nslice <- 10

my.p.slice <- seq(range(my.p)[1], range(my.p)[2],length.out=nslice+1)
my.p.slice[nslice+1]=my.p.slice[nslice+1]+1
my.sample <- c(1, length(my.p))
# nbasis <- max(30, ceiling(10 * n^(2/9)))
nobs.slice <- 10

set.seed(3458)
for(i in 1:(nslice)){  
  my.which <- which( my.p.slice[i]<= my.p   &  my.p  < my.p.slice[i+1])
  if( length(my.which) >= nobs.slice){
    my.sample <- union(my.sample,sample(my.which, nobs.slice ))
  } else {    
    my.sample <- sort(union(my.sample, my.which))
  }  
}
length(my.sample)

t0=proc.time()
my.abs <- gssanova(cbind(my.y, my.N-my.y)~my.x*my.z, id.basis=my.sample,
                   family="binomial",  skip.iter=F)
my.est.abs <- predict(my.abs, data.frame(my.x=my.x, my.z=my.z), se.fit=T)
t.abs=proc.time()-t0

t.abs

Fit.abs=1-1/(1+exp(my.est.abs$fit))


sum.abs = summary(my.abs)
rsq.abs = round(1-sum.abs$deviance/sum.abs$dev.null, digits=2)
rsq.abs

## Plot ##

gen1=1:max(which(dat$generation=='119_1'))
gen2=min(which(dat$generation!='119_1')):n

gen11=which(dat$group=='119_1_r1')
gen12=which(dat$group=='119_1_r2')
gen21=which(dat$group=='12_55_r1')
gen22=which(dat$group=='12_55_r2')


par(mfrow=c(2,2), mar=c(5,4,4,2)+.1, mex=.75)
my.range=range(dat$mc,Fit.abs*dat$h)
plot(dat$position[gen11]-13700000, dat$mc[gen11], type='l',col='gray', ylim=my.range,
     xlab='', ylab='Observed count', main='Generation 1')
plot(dat$position[gen21]-13700000, dat$mc[gen21], type='l',col='gray', ylim=my.range,
     xlab='', ylab='', main='Generation 2')

plot(dat$position[gen11]-13700000, (Fit.abs*dat$h)[gen11], type='l', ylim=my.range,
     xlab='', ylab='Expected count', main='')
plot(dat$position[gen22]-13700000, (Fit.abs*dat$h)[gen22], type='l', ylim=my.range,
     xlab='', ylab='', main='')
