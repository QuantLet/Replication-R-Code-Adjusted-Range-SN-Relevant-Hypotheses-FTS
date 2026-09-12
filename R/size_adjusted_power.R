source('R/primary_simulations.R')
run_size_adjusted_power <- function(boot=400){files<-list.files(results_dir,pattern='^cp_.*rds$',full.names=TRUE);runs<-lapply(files,readRDS);rows<-list()
 for(tr in runs){c<-tr$config;if(c$state!='boundary')next;ca<-c;ca$state<-'local_power';ev<-readRDS(file.path(results_dir,paste0('cp_',digest::digest(ca,algo='sha256'),'.rds')))
  labels<-c('fitted_AR','fitted_Q','fitted_Wald');score<-function(obj){M<-obj$lower_bounds; (M[,'distance']-1)/(M[,'distance']-M[,labels])}
  train<-score(tr);test<-score(ev);q<-apply(train,2,quantile,.95,type=1);power<-colMeans(sweep(test,2,q,`>`))
  set.seed(config_seed(c,'revision_paired_bootstrap'));diffs<-replicate(boot,{u<-train[sample.int(nrow(train),replace=TRUE),];w<-test[sample.int(nrow(test),replace=TRUE),];qq<-apply(u,2,quantile,.95,type=1);pp<-colMeans(sweep(w,2,qq,`>`));c(pp[1]-pp[2],pp[1]-pp[3])});ci<-t(apply(diffs,1,quantile,c(.025,.975)))
  rows[[length(rows)+1]]<-data.frame(c,AR=power[1],Q=power[2],Wald=power[3],AR_minus_Q_low=ci[1,1],AR_minus_Q_high=ci[1,2],AR_minus_Wald_low=ci[2,1],AR_minus_Wald_high=ci[2,2],train_n=nrow(train),evaluation_n=nrow(test),bootstrap=boot)
 };write.csv(do.call(rbind,rows),file.path(results_dir,'adjusted_power.csv'),row.names=FALSE)
}
