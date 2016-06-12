# Hsc70-4 V5~V16
rm(list=ls())

library(gss)

dat <- read.csv("Hsc70-4.csv")
attach(dat)

gene.names = unique(V1)
gene.length = rep(NULL, length(gene.names))

for(i in 1:length(gene.names))
{
  gene.length[i] = sum(V1==gene.names[i])
}

my.dist = 76

Fit1 = matrix(NA, 12, dim(dat)[1])
Fit2 = matrix(NA, 12, dim(dat)[1])
time1 = time2 = NA

dev1 = dev.null1 = Rsquare1 = NA
dev2 = dev.null2 = Rsquare2 = NA


for(ii in 5:16)
{
  print(ii)
  
  my.y <- dat[,ii]
  n <- length(my.y)
  
  # Appearances of GC in far mid and near neighbourhood
  G.far = G.mid = G.near = rep(0, length(my.y))
  C.far = C.mid = C.near = rep(0, length(my.y))
  # Total counts of GC in far mid and near neighbourhood
  G.far.count = G.mid.count = G.near.count = rep(0, length(my.y))
  C.far.count = C.mid.count = C.near.count = rep(0, length(my.y))
  
  for(k in 1:length(gene.names))
  {
    start = min(which(V1==gene.names[k]))
    end   = max(which(V1==gene.names[k]))
    
    
    for(i in start:end)
    {
      far  = c(max(start, i-3*my.dist):max(start, i-2*my.dist-1), 
               min(end, i+2*my.dist+1):min(end, i+3*my.dist))
      mid  = c(max(start, i-2*my.dist):max(start, i-1*my.dist-1), 
               min(end, i+1*my.dist+1):min(end, i+2*my.dist))
      near = c(max(start, i-1*my.dist):max(start, i-0*my.dist-1), 
               min(end, i+0*my.dist+1):min(end, i+1*my.dist))
      
      G.far.index <- which(V4[far]=="G")
      G.mid.index <- which(V4[mid]=="G")
      G.near.index <- which(V4[near]=="G")
      C.far.index <- which(V4[far]=="C")
      C.mid.index <- which(V4[mid]=="C")
      C.near.index <- which(V4[near]=="C")
      
      G.far[i] <- length(G.far.index)
      G.mid[i] <- length(G.mid.index)
      G.near[i] <- length(G.near.index)
      C.far[i] <- length(C.far.index)
      C.mid[i] <- length(C.mid.index)
      C.near[i] <- length(C.near.index)
      
      G.far.count[i] <- sum(my.y[G.far.index])
      G.mid.count[i] <- sum(my.y[G.mid.index])
      G.near.count[i] <- sum(my.y[G.near.index])
      C.far.count[i] <- sum(my.y[C.far.index])
      C.mid.count[i] <- sum(my.y[C.mid.index])
      C.near.count[i] <- sum(my.y[C.near.index])
      
    }
  }
  
  GC.far = G.far + C.far
  GC.mid = G.mid + C.mid
  GC.near= G.near+ C.near
  
  GC.far.count = G.far.count + C.far.count
  GC.mid.count = G.mid.count + C.mid.count
  GC.near.count= G.near.count+ C.near.count
  
  
  
  ## slicing
  my.y.slice <- hist(my.y,breaks="Scott", plot=FALSE)$breaks
  
  my.sample <- c(1, length(my.y))
  
  nslice <- length(my.y.slice)-1
  nbasis <- max(30, ceiling(10* n^(2/9)))
  nobs.slice <- ceiling(nbasis/nslice)
  
  set.seed(3458)
  for(i in 1:(nslice)){  
    my.which <- which( my.y.slice[i]<= my.y   &  my.y  < my.y.slice[i+1])
    if( length(my.which) >= nobs.slice){
      my.sample <- union(my.sample,sample(my.which, nobs.slice ))
    } else {    
      my.sample <- sort(union(my.sample, my.which))
    }  
  }
  print(length(my.sample))
  
  
  # poisson 1
  t0=proc.time()
  my.ss1 <- gssanova(my.y ~ GC.far * GC.mid * GC.near, 
                     id.basis=my.sample, type="linear", family="poisson", skip.iter=F)
  my.est1 <- predict(my.ss1, data.frame(GC.far=GC.far, GC.mid=GC.mid, GC.near=GC.near), se.fit=T)
  t1=proc.time()-t0
  
  my.ss1.sum <- summary(my.ss1)
  Fit1[ii-4,]=exp(my.est1$fit)
  
  # poisson 2
  t0=proc.time()
  my.ss2 <- gssanova(my.y ~ GC.far.count * GC.mid.count * GC.near.count, 
                     id.basis=my.sample, type="linear", family="poisson", skip.iter=F)
  my.est2 <- predict(my.ss2, data.frame(GC.far.count=GC.far.count, GC.mid.count=GC.mid.count, 
                                        GC.near.count=GC.near.count), se.fit=T)
  t2=proc.time()-t0
  
  my.ss2.sum <- summary(my.ss2)
  Fit2[ii-4,]=exp(my.est2$fit)
  
  time1[ii-4] = round(t1, digits=0)[3]
  time2[ii-4] = round(t2, digits=0)[3]
  
  dev1[ii-4] = my.ss1.sum$deviance
  dev.null1[ii-4] = my.ss1.sum$dev.null
  Rsquare1[ii-4] = round(1-my.ss1.sum$deviance/my.ss1.sum$dev.null, digits=2)
  
  dev2[ii-4] = my.ss2.sum$deviance
  dev.null2[ii-4] = my.ss2.sum$dev.null
  Rsquare2[ii-4] = round(1-my.ss2.sum$deviance/my.ss2.sum$dev.null, digits=2)
  
}

Rsquare=data.frame(poisson_app=Rsquare1, poisson_count=Rsquare2)

## Plot ##
par(mfrow=c(1,2), mar=c(5,4,4,2)+.1, mex=.75)

load('Hsc70-4_each_pois.RData')
ii=3
plot(1:dim(dat)[1],dat[,ii+4], type='l',col='gray', ylim=c(0,800), lwd=2,
     xlab='', ylab='Count', main='Time course 3')
lines(1:dim(dat)[1], Fit1[ii,])
ii=6
plot(1:dim(dat)[1],dat[,ii+4], type='l',col='gray', ylim=c(0,800), lwd=2,
     xlab='', ylab='', main='Time course 6')
lines(1:dim(dat)[1], Fit1[ii,])

