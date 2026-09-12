source('R/primary_simulations.R')
run_variance_diagnostics <- function() {
rows<-list()
for(f in list.files(results_dir,pattern='^cp_.*rds$',full.names=TRUE)){
 run<-readRDS(f);c<-run$config;if(c$state!='boundary')next
 n<-c$n;k0<-floor(n/2);a2<-n*diff(profile(seq(0,1,length.out=n+1),c$profile));rho<-if(c$dep=='switch')ifelse(seq_len(n)<=floor(.6*n),.2,.8) else rep(as.numeric(sub('ar','',c$dep)),n)
 vv<-a2*if(c$dep=='switch')(1+rho)/(1-rho)/4.5 else 1
 target<-4*(sum(vv[1:k0])/k0^2+sum(vv[(k0+1):n])/(n-k0)^2)
 # A known-split lower bound's distance differs from the estimated-split distance;
 # reconstruct the exact known-split distance from the original deterministic data seed.
 set.seed(config_seed(c,'revision_data'));seeds<-sample.int(.Machine$integer.max,nrow(run$lower_bounds));ratios<-numeric(length(seeds))
 for(i in seq_along(seeds)){set.seed(seeds[i]);x<-generate_noise(n,c$dep,c$tail);x[,1]<-sqrt(a2)*x[,1];x[1:k0,1]<-x[1:k0,1]+1;s<-segment_stats(x,k0,max(2,round(n^.4)));ratios[i]<-s$se2/target}
 rows[[length(rows)+1]]<-data.frame(c,median_variance_ratio=median(ratios),q10=quantile(ratios,.1),q90=quantile(ratios,.9))
}
ans<-do.call(rbind,rows)
write.csv(ans,file.path(results_dir,'variance_ratios.csv'),row.names=FALSE)
invisible(ans)
}
