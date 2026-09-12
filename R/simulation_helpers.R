source('R/core.R')
QREF<-c(AR=1.3973901267,Q=5.32267990053)
`%or%`<-function(x,y)if(is.null(x))y else x
config_seed<-function(c,stream='data')strtoi(substr(digest::digest(list(c,stream),algo='sha256'),1,7),16L)
# A deterministic configuration seed and separate reference stream support resumption.
simulate_design <- function(c,reps=c$reps) {
 n<-c$n;r<-seq(0,1,length.out=(c$grid%or%40)+1);v<-profile(seq(0,1,length.out=n+1),c$profile);a<-sqrt(n*pmax(diff(v),0))
 model<-c$model%or%if(c$suite=='cp')'C' else 'L';dep<-c$dep%or%paste0('ar',c$rho%or%.4)
 d<-c[["d"]]%or%c$distance%or%switch(c$state%or%'boundary',boundary=1,local_power=1+if(model=='C')8/sqrt(n) else 4/sqrt(n),interior_null=.8,no_break=0,fixed_power=1.5)
 K<-c$K%or%3;out<-matrix(NA_real_,reps,if(model=='C')11 else 3)
 set.seed(config_seed(c));seeds<-sample.int(.Machine$integer.max,reps)
 for(i in seq_len(reps)){
  set.seed(seeds[i]);e<-errors(n,dep,c$tail,K)
  if(K>3){j<-seq_len(K-1);weights<-.6*j^(-.85)/sum(j^(-1.275))^(2/3);e[,-1]<-sweep(e[,-1,drop=FALSE],2,weights,'*')}
  if(model=='L')out[i,]<-loading_stats(e,a,d)$piv else{
   k0<-round((c$theta%or%.4)*n);e[,1]<-a*e[,1];e[1:k0,1]<-e[1:k0,1]+sqrt(d)
   k<-estimate_split(e,c$trim%or%.15);b<-max(2,round((c$block_multiplier%or%1)*n^.4));s<-cp_stats(e,k,b,r)
   set.seed(config_seed(c,'reference')+i);ref<-brownian_refs(s$v,r,c$reference_draws%or%999)
   known<-cp_stats(e,k0,b,r);rawref<-brownian_refs(s$raw_v,r,z=matrix(rnorm((c$reference_draws%or%999)*(length(r)-1)),ncol=length(r)-1))
   # Known-split and fitted-split references share Gaussian increments.
   set.seed(config_seed(c,'reference')+i);knownref<-brownian_refs(known$v,r,c$reference_draws%or%999)
   out[i,]<-c(rank_p(ref,s$piv),k,s$piv,rank_p(knownref,known$piv),rank_p(rawref,s$piv),s$d,abs(k-k0))
  }
 };colnames(out)<-if(model=='C')c('AR_p','Q_p','split','AR_pivot','Q_pivot','known_AR_p','known_Q_p','raw_AR_p','raw_Q_p','distance','split_error') else c('AR','Q','calendarQ');out
}
run_plan <- function(plan_path,reps_override=NULL,outdir='results/R/designs') {
 dir.create(outdir,recursive=TRUE,showWarnings=FALSE);plan<-jsonlite::fromJSON(plan_path,simplifyVector=FALSE);if(!is.null(plan$plan))plan<-plan$plan
 for(j in seq_along(plan)){c<-plan[[j]];if(!is.null(reps_override))c$reps<-reps_override
  path<-file.path(outdir,paste0(digest::digest(c,algo='sha256'),'.rds'));if(file.exists(path))next
  # High-dimensional experiments preserve paired dimensions using shared raw innovations.
  if(!is.null(c$model))res<-paired_highdim(c) else res<-simulate_design(c)
  path<-file.path(outdir,paste0(digest::digest(c,algo='sha256'),'.rds'))
  saveRDS(list(config=c,values=res),path);cat('R design',j,'of',length(plan),'completed\n')
 }
}
paired_highdim <- function(c) {
 n<-c$n;v<-profile(seq(0,1,length.out=n+1),c$profile);a<-sqrt(n*pmax(diff(v),0));r<-seq(0,1,length.out=41)
 out<-list('50'=matrix(NA_real_,c$reps,if(c$model=='L')2 else 3),'100'=matrix(NA_real_,c$reps,if(c$model=='L')2 else 3))
 set.seed(config_seed(c));seeds<-sample.int(2^30,c$reps)
 for(i in seq_len(c$reps)){
  set.seed(seeds[i]);raw<-matrix(rnorm((n+300)*100),n+300,100);if(c$tail=='heavy')raw[,-1]<-matrix(rt((n+300)*99,1.6),n+300,99)
  e100<-errors(n,paste0('ar',c$rho),c$tail,100,raw)
  set.seed(config_seed(c,'reference')+i);zref<-matrix(rnorm(999*40),999,40)
  for(K in c(50,100)){
   e<-e100[,1:K];j<-seq_len(K-1);e[,-1]<-sweep(e[,-1],2,.6*j^(-.85)/sum(j^(-1.275))^(2/3),'*')
   if(c$model=='L')out[[as.character(K)]][i,]<-loading_stats(e,a,c[["d"]])$piv[1:2] else{
    e[,1]<-e[,1]*a;e[1:floor(.4*n),1]<-e[1:floor(.4*n),1]+sqrt(c[["d"]]);k<-estimate_split(e,.15);s<-cp_stats(e,k,round(n^.4),r);out[[as.character(K)]][i,]<-c(rank_p(brownian_refs(s$v,r,z=zref),s$piv),k)
   }
  }
 };out
}
confirmation <- function(reps=5000,boot=500,outdir='results/R') {
 dir.create(outdir,recursive=TRUE,showWarnings=FALSE);rows<-list()
 for(p in c('uniform','r5','step','plateau'))for(dep in c('ar0.4','ar0.8','ma0.7'))for(tail in c('gaussian','orthogonal_t1.6')){
  c<-list(n=1000,profile=p,dep=dep,tail=tail,suite='matched',reps=reps,d=1)
  train<-simulate_design(c)[,1:2];c[["d"]]<-1+4/sqrt(1000);test<-simulate_design(c)[,1:2]
  q<-apply(train,2,quantile,.95);size<-colMeans(sweep(train,2,QREF,`>`));power<-colMeans(sweep(test,2,q,`>`))
  set.seed(config_seed(c,'bootstrap'));diffs<-replicate(boot,{tr<-train[sample.int(reps,reps,TRUE),];ev<-test[sample.int(reps,reps,TRUE),];pr<-colMeans(sweep(ev,2,apply(tr,2,quantile,.95),`>`));pr[1]-pr[2]})
  ci<-quantile(diffs,c(.025,.975));rows[[length(rows)+1]]<-data.frame(profile=p,dep=dep,tail=tail,n=1000,reps=reps,ar_size=size[1],q_size=size[2],ar_adjusted_power=power[1],q_adjusted_power=power[2],adjusted_difference=power[1]-power[2],difference_ci_low=ci[1],difference_ci_high=ci[2])
 };write.csv(do.call(rbind,rows),file.path(outdir,'paired_results.csv'),row.names=FALSE)
}
ordinary_study <- function(reps=1500,B=2048,outdir='results/R') {
 dir.create(outdir,recursive=TRUE,showWarnings=FALSE);n<-1000;r<-seq(0,1,length.out=41);ix<-floor(n*r);rows<-list();set.seed(7192026);zref<-qnorm(qrng::sobol(B,40,randomize='Owen',seed=7192026))
 oldref<-function(v){w<-cbind(0,t(apply(sweep(zref,2,sqrt(pmax(diff(v),0)),'*'),1,cumsum)));g<-(w-outer(w[,41],r))*rep(r,each=B);w[,41]/sqrt(rowMeans(g[,2:40]^2))}
 qstationary<-quantile(oldref(r),.95)
 for(p in c('uniform','r5','step','plateau'))for(state in c('boundary','local_power')){
  d<-if(state=='boundary')1 else 1+4/sqrt(n);a<-sqrt(n*diff(profile(seq(0,1,length.out=n+1),p)));hits<-matrix(FALSE,reps,4)
  for(i in seq_len(reps)){x<-matrix(rnorm(n*3),n,3);x<-sweep(x,2,c(1,.3,.2),'*');x[,1]<-sqrt(d)+a*x[,1];mu<-colMeans(x);sc<-as.vector(2*sweep(x,2,mu)%*%mu);path<-c(0,cumsum(sc))[ix+1]/n;path[c(1,41)]<-0;v<-c(0,cumsum(sc^2))[ix+1]/sum(sc^2)
   G<-rowSums((csrows(x)[ix+1,,drop=FALSE]/n)^2)-r^2*sum(mu^2);num<-sum(mu^2)-1;orig<-num/sqrt(mean(G[2:40]^2));piv<-num/c(diff(range(path)),sqrt(mean(path[2:40]^2)));crit<-apply(brownian_refs(v,r,z=zref),2,quantile,.95)
   hits[i,]<-c(orig>qstationary,orig>quantile(oldref(v),.95),piv>crit)
  };rows[[length(rows)+1]]<-data.frame(profile=p,state=state,method=c('original_stationary_Q','updated_Q','projected_AR','projected_Q'),rate=colMeans(hits),reps=reps)
 };write.csv(do.call(rbind,rows),file.path(outdir,'ordinary_mean.csv'),row.names=FALSE)
}
run_wild_size <- function(reps=1000,B=499,outdir='results/R') {
 rows<-list();r<-seq(0,1,length.out=41)
 for(p in c('uniform','r5'))for(dep in c('ar0.4','ar0.8'))for(tail in c('gaussian','orthogonal_t1.6')){
  set.seed(config_seed(list(p,dep,tail)));hits<-matrix(FALSE,reps,2);a<-sqrt(70*diff(profile(seq(0,1,length.out=71),p)))
  for(i in seq_len(reps)){x<-errors(70,dep,tail);x[,1]<-x[,1]*a;x[1:35,1]<-x[1:35,1]+1;k<-estimate_split(x,.2);s<-cp_stats(x,k,5);ref<-wild_refs(x,k,5,B);hits[i,]<-rank_p(ref,s$piv)<=.05}
  rows[[length(rows)+1]]<-data.frame(profile=p,dep=dep,tail=tail,method=c('AR','Q'),rejection_rate=colMeans(hits),reps=reps,reference_draws=B)
 };dir.create(outdir,recursive=TRUE,showWarnings=FALSE);write.csv(do.call(rbind,rows),file.path(outdir,'wild_size.csv'),row.names=FALSE)
}
