source('R/core.R')
run_limiting_power <- function(batches=20,per=5000,grids=c(256,1024,4096),outdir='results/R') {
 dir.create(outdir,recursive=TRUE,showWarnings=FALSE);set.seed(129122026);H<-c(-.5,0,.1,.25,.5,.75,1,1.5,2,2.5,3,4);rows<-list();G<-max(grids)
 for(b in seq_len(batches)){
  stages<-vector('list',2)
  for(stage in 1:2){vals<-lapply(grids,function(g)matrix(NA_real_,per,2))
   for(i in seq_len(per)){w<-c(0,cumsum(rnorm(G)/sqrt(G)));bridge<-w-seq(0,1,length.out=G+1)*tail(w,1)
    for(j in seq_along(grids)){g<-grids[j];x<-bridge[seq(1,G+1,by=G/g)];vals[[j]][i,]<-c(diff(range(x)),sqrt(trapz(x*x,seq(0,1,length.out=g+1))))}
   };stages[[stage]]<-vals
  }
  for(j in seq_along(grids))for(m in 1:2){tr<-stages[[1]][[j]][,m];ev<-stages[[2]][[j]][,m];q<-uniroot(function(q)mean(pnorm(-q*tr))-.05,c(.01,30),tol=1e-10)$root
   for(h in H)rows[[length(rows)+1]]<-data.frame(batch=b,grid=grids[j],method=c('Range-SN','Q-SN')[m],h=h,q=q,power=mean(pnorm(h-q*ev)),slope=mean(dnorm(q*ev-h)))
  };cat('R Brownian batch',b,'of',batches,'\n')
 };a<-do.call(rbind,rows);write.csv(a,file.path(outdir,'power_batches.csv'),row.names=FALSE)
 f<-function(x)c(mean=mean(x),sem=sd(x)/sqrt(length(x)))
 s<-aggregate(power~grid+method+h,a,f);write.csv(data.frame(s[1:3],power=s$power[,1],power_se=s$power[,2]),file.path(outdir,'limiting_power.csv'),row.names=FALSE)
 z<-reshape(a[c('batch','grid','method','h','power')],idvar=c('batch','grid','h'),timevar='method',direction='wide');z$difference<-z[['power.Range-SN']]-z[['power.Q-SN']];s<-aggregate(difference~grid+h,z,f);write.csv(data.frame(s[1:2],mean=s$difference[,1],sem=s$difference[,2]),file.path(outdir,'power_differences.csv'),row.names=FALSE)
}
