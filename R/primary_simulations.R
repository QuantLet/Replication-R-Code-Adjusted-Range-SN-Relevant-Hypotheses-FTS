# Primary simulation experiments for the paper.
source('R/core.R')
source('R/simulation_helpers.R')
results_dir <- 'results'
dir.create(results_dir,recursive=TRUE,showWarnings=FALSE)
ref_quantile <- function(x) { chk<-reference_check(x);if(!chk$valid)stop(chk$reason);chk$q }
# Exact observation indices and a shared deterministic b across candidate splits.
segment_stats <- function(x,k,b,r=seq(0,1,length.out=41),direction=NULL) {
 n<-nrow(x);m<-c(k,n-k);mu<-list(colMeans(x[1:k,,drop=FALSE]),colMeans(x[(k+1):n,,drop=FALSE]));D<-mu[[1]]-mu[[2]]
 u<-if(is.null(direction))D else direction
 sc<-list(as.vector(2*sweep(x[1:k,,drop=FALSE],2,mu[[1]])%*%u),as.vector(2*sweep(x[(k+1):n,,drop=FALSE],2,mu[[2]])%*%u))
 rx<-lapply(m,function(mm)floor(mm*r)/mm)
 paths<-lapply(1:2,function(j)c(0,cumsum(sc[[j]]))[floor(m[j]*r)+1]/m[j]);g<-paths[[1]]-paths[[2]];g[c(1,length(r))]<-0
 qb<-lapply(1:2,function(j)block_cumulative(sc[[j]],rx[[j]],b)/m[j]^2);tot<-vapply(qb,tail,numeric(1),1);se2<-sum(tot)
 valid<-is.finite(se2)&&se2>0
 v<-lapply(1:2,function(j)if(tot[j]>0)qb[[j]]/tot[j] else rx[[j]])
 list(d=sum(D^2),D=D,k=k,m=m,b=b,r=r,rx=rx,v=v,weights=if(valid)tot/se2 else c(NA,NA),se2=se2,
      path=g,den=c(AR=diff(range(g)),Q=sqrt(mean(g[-c(1,length(g))]^2))),valid=valid)
}
# Two independent segment Brownian paths; endpoints and bridges share increments.
segment_refs <- function(s,z) {
 if(!s$valid)stop('Zero or invalid total block variance; inference unavailable')
 g<-matrix(0,nrow(z),length(s$r));endpoint<-numeric(nrow(z));M<-length(s$r)-1
 for(j in 1:2){w<-cbind(0,matrixStats::rowCumsums(sweep(z[,((j-1)*M+1):(j*M),drop=FALSE],2,sqrt(pmax(diff(s$v[[j]]),0)),'*')))
  end<-w[,M+1];sgn<-if(j==1)1 else -1;fac<-sgn*sqrt(s$weights[j]);endpoint<-endpoint+fac*end;g<-g+fac*(w-outer(end,s$rx[[j]]))
 }
 den<-cbind(AR=matrixStats::rowMaxs(g)-matrixStats::rowMins(g),Q=sqrt(rowMeans(g[,-c(1,ncol(g)),drop=FALSE]^2)))
 if(any(!is.finite(den))||any(den<=0))stop('Degenerate Gaussian reference; no nominal inference returned')
 endpoint/den
}
generate_noise <- function(n,dep,tail) {
 raw<-matrix(rnorm((n+300)*3),n+300,3)
 if(tail%in%c('orthogonal_t1.6','directional_t2.5','orthogonal_t1.4'))raw[,2]<-rt(n+300,if(tail=='orthogonal_t1.4')1.4 else 1.6)
 if(tail=='directional_t2.5')raw[,1]<-rt(n+300,2.5)/sqrt(5)
 errors(n,dep,'gaussian',3,raw)
}
oracle_profile <- function(s,variance) {
 ix<-list(seq_len(s$k),(s$k+1):sum(s$m));total<-vapply(1:2,function(j)sum(variance[ix[[j]]])/s$m[j]^2,numeric(1))
 s$weights<-total/sum(total)
 s$v<-lapply(1:2,function(j){a<-c(0,cumsum(variance[ix[[j]]]));if(tail(a,1)>0)a[floor(s$m[j]*s$r)+1]/tail(a,1) else s$rx[[j]]})
 s
}
run_change_point_study <- function(c,reps=2000,B=999) {
 key<-digest::digest(c,algo='sha256');dest<-file.path(results_dir,paste0('cp_',key,'.rds'));if(file.exists(dest)){cached<-readRDS(dest);if(nrow(cached$lower_bounds)==reps&&cached$summary$B[1]==B)return(cached$summary)}
 set.seed(config_seed(c,'revision_data'));seeds<-sample.int(.Machine$integer.max,reps)
 reference_seed<-config_seed(c,'revision_reference')
 n<-c$n;k0<-floor(.5*n);d<-if(c$state=='boundary')1 else 1+8/sqrt(n);v<-profile(seq(0,1,length.out=n+1),c$profile);a<-sqrt(n*diff(v))
 rho<-if(c$dep=='switch')ifelse(seq_len(n)<=floor(.6*n),.2,.8) else rep(as.numeric(sub('ar','',c$dep)),n)
 varprof<-a^2*if(c$dep=='switch')(1+rho)/(1-rho)/4.5 else 1
 # Extra b choices diagnose sensitivity in a prespecified subset.
 multipliers<-if(c$profile=='uniform'&&c$dep=='ar0.4'&&c$tail=='gaussian'&&n%in%c(70,1000))c(.5,1,2) else 1
 variant<-c('fitted','known_split','oracle_profile','oracle_direction','calendar_reference')
 method<-c('AR','Q','Wald');labels<-as.vector(outer(variant,method,paste,sep='_'));labels<-c(labels,unlist(lapply(setdiff(multipliers,1),function(bm)paste0('block',bm,'_',method))))
 out<-matrix(NA_real_,reps,length(labels)+5,dimnames=list(NULL,c(labels,'split_error','distance','remainder_scaled','direction_error','invalid')))
 for(i in seq_len(reps)){
  set.seed(seeds[i]);noise<-generate_noise(n,c$dep,c$tail);noise[,1]<-a*noise[,1];x<-noise;x[1:k0,1]<-x[1:k0,1]+sqrt(d)
  set.seed(reference_seed+i);z<-matrix(rnorm(B*80),B,80)
  k<-estimate_split(x,.2);b<-max(2,round(n^.4));s<-segment_stats(x,k,b);known<-segment_stats(x,k0,b);od<-segment_stats(x,k,b,direction=c(sqrt(d),0,0));sp<-oracle_profile(s,varprof)
  liststats<-list(fitted=s,known_split=known,oracle_profile=sp,oracle_direction=od,calendar_reference=s)
  for(nm in names(liststats)){
   st<-liststats[[nm]];if(!st$valid)next
   ref<-if(nm=='calendar_reference'){old<-cp_stats(x,k,b);brownian_refs(old$v,z=z[,1:40])} else segment_refs(st,z)
   q<-ref_quantile(ref);lower<-c(st$d-q*st$den,st$d-qnorm(.95)*sqrt(st$se2))
   out[i,paste(nm,method,sep='_')]<-lower
  }
  for(bm in setdiff(multipliers,1)){st<-segment_stats(x,k,max(2,round(bm*n^.4)));q<-ref_quantile(segment_refs(st,z));out[i,paste0('block',bm,'_',method)]<-c(st$d-q*st$den,st$d-qnorm(.95)*sqrt(st$se2))}
  rem<-colMeans(noise[1:k0,])-colMeans(noise[(k0+1):n,]);out[i,c('split_error','distance','remainder_scaled','direction_error','invalid')]<-c(abs(k-k0),s$d,sqrt(n)*sum(rem^2),sqrt(sum((s$D-c(sqrt(d),0,0))^2)),!s$valid)
 }
 rows<-lapply(labels,function(label){low<-out[,label];ok<-is.finite(low);rate<-mean(low[ok]>1);data.frame(c,comparison=label,rate=rate,mcse=sqrt(rate*(1-rate)/sum(ok)),coverage=mean(low[ok]<=d),reps=reps,valid=sum(ok),B=B,median_split_error=median(out[,'split_error']),q90_split_error=quantile(out[,'split_error'],.9),median_remainder=median(out[,'remainder_scaled']),q90_remainder=quantile(out[,'remainder_scaled'],.9),median_direction_error=median(out[,'direction_error']))})
 ans<-do.call(rbind,rows);saveRDS(list(config=c,lower_bounds=out,summary=ans),dest);cat('Finished CP',n,c$profile,c$dep,c$tail,c$state,'\n');ans
}
run_primary_simulations <- function(reps=2000,B=999,cores=6) {
 configs<-expand.grid(n=c(70,200,1000),profile=c('uniform','r5'),dep=c('ar0.4','ar0.8'),tail='gaussian',state=c('boundary','local_power'),stringsAsFactors=FALSE)
 configs<-rbind(configs,expand.grid(n=c(200,1000),profile='uniform',dep='ar0.4',tail=c('directional_t2.5','orthogonal_t1.4'),state=c('boundary','local_power'),stringsAsFactors=FALSE),expand.grid(n=c(200,1000),profile='step',dep='switch',tail='orthogonal_t1.6',state=c('boundary','local_power'),stringsAsFactors=FALSE))
 write.csv(configs,file.path(results_dir,'cp_designs.csv'),row.names=FALSE)
 ans<-parallel::mclapply(seq_len(nrow(configs)),function(i)run_change_point_study(as.list(configs[i,]),reps,B),mc.cores=cores,mc.preschedule=FALSE)
 if(any(vapply(ans,inherits,logical(1),'try-error')))stop('At least one simulation design failed')
 write.csv(do.call(rbind,ans),file.path(results_dir,'cp_results.csv'),row.names=FALSE)
}
loading_refs <- function(v,B=19999,seed=9192026) {
 set.seed(seed);n<-length(v)-1;ans<-matrix(NA_real_,B,2,dimnames=list(NULL,c('AR','Q')))
 for(lo in seq(1,B,by=1000)){ix<-lo:min(B,lo+999);z<-matrix(rnorm(length(ix)*n),length(ix));w<-cbind(0,matrixStats::rowCumsums(sweep(z,2,sqrt(diff(v)),'*')));end<-w[,n+1];g<-w-outer(end,v)
  den<-cbind(matrixStats::rowMaxs(g)-matrixStats::rowMins(g),sqrt(as.vector((g[,-1,drop=FALSE]^2+g[,-(n+1),drop=FALSE]^2)%*%diff(v)/2)))
  ans[ix,]<-end/den
 };ans
}
run_loading_simulations <- function(reps=2000,B=19999) {
 rows<-list()
 for(n in c(70,200,1000))for(p in c('uniform','r5')){
  v<-profile(seq(0,1,length.out=n+1),p);a<-sqrt(n*diff(v));q<-ref_quantile(loading_refs(v,B,config_seed(list(n,p),'loading_reference')))
  for(dep in c('iid','ar0.4'))for(state in c('boundary','local_power')){
   c<-list(n=n,profile=p,dep=dep,state=state);set.seed(config_seed(c,'loading_data'));d<-if(state=='boundary')1 else 1+4/sqrt(n);hit<-matrix(FALSE,reps,4)
   for(i in seq_len(reps)){st<-loading_stats(generate_noise(n,dep,'gaussian'),a,d);hit[i,]<-c(st$piv[1:2]>QREF,st$piv[1:2]>q)}
   rows[[length(rows)+1]]<-data.frame(c,method=c('continuous_AR','continuous_Q','grid_AR','grid_Q'),rate=colMeans(hit),mcse=sqrt(colMeans(hit)*(1-colMeans(hit))/reps),reps=reps,B=B)
  };cat('Loading grids',n,p,'completed\n')
 };write.csv(do.call(rbind,rows),file.path(results_dir,'loading_results.csv'),row.names=FALSE)
}
run_mean_simulations <- function(reps=2000,B=8192) {
 n<-1000;r<-seq(0,1,length.out=41);ix<-floor(n*r);b<-round(n^.4);rows<-list()
 z<-qnorm(qrng::sobol(B,40,randomize='Owen',seed=7199226))
 make_ref<-function(v){w<-cbind(0,matrixStats::rowCumsums(sweep(z,2,sqrt(pmax(diff(v),0)),'*')));end<-w[,41];g<-w-outer(end,r);orig<-sweep(g,2,r,'*');cbind(original_AR=end/(matrixStats::rowMaxs(orig)-matrixStats::rowMins(orig)),original_Q=end/sqrt(rowMeans(orig[,2:40]^2)),projected_AR=end/(matrixStats::rowMaxs(g)-matrixStats::rowMins(g)),projected_Q=end/sqrt(rowMeans(g[,2:40]^2)))}
 qstat<-ref_quantile(make_ref(r))[2]
 for(p in c('uniform','r5','step','plateau'))for(state in c('boundary','local_power')){
  c<-list(n=n,profile=p,state=state);set.seed(config_seed(c,'ordinary_revision'));d<-if(state=='boundary')1 else 1+4/sqrt(n);a<-sqrt(n*diff(profile(seq(0,1,length.out=n+1),p)));hit<-matrix(FALSE,reps,6)
  for(i in seq_len(reps)){x<-generate_noise(n,'iid','gaussian');x[,1]<-sqrt(d)+a*x[,1];mu<-colMeans(x);sc<-as.vector(2*sweep(x,2,mu)%*%mu);g<-c(0,cumsum(sc))[ix+1]/n;g[c(1,41)]<-0;qb<-block_cumulative(sc,r,b);v<-qb/tail(qb,1);q<-ref_quantile(make_ref(v));G<-rowSums((csrows(x)[ix+1,,drop=FALSE]/n)^2)-r^2*sum(mu^2);num<-sum(mu^2)-1
   den<-c(diff(range(G)),sqrt(mean(G[2:40]^2)),diff(range(g)),sqrt(mean(g[2:40]^2)))
   hit[i,]<-c(num>q*den,num>qnorm(.95)*sqrt(tail(qb,1))/n,num>qstat*den[2])
  };rate<-colMeans(hit);rows[[length(rows)+1]]<-data.frame(c,method=c('original_AR','original_Q','projected_AR','projected_Q','Wald','stationary_Q'),rate=rate,mcse=sqrt(rate*(1-rate)/reps),reps=reps,B=B,block=b);cat('Ordinary',p,state,'completed\n')
 };write.csv(do.call(rbind,rows),file.path(results_dir,'ordinary_results.csv'),row.names=FALSE)
}
# Floating-point tie screen, not a proof of analytic quantile regularity.
reference_check <- function(ref,alpha=.05,tol=1e-10) {
 if(any(!is.finite(ref)))return(list(valid=FALSE,reason='nonfinite reference',q=c(NA,NA)))
 q<-apply(ref,2,quantile,1-alpha,type=1,names=FALSE)
 mass<-colMeans(abs(sweep(ref,2,q,'-'))<=tol*pmax(1,rep(abs(q),each=nrow(ref))))
 list(valid=all(mass<alpha),reason=if(any(mass>=alpha))'atom at selected quantile' else 'numerical screen passed',q=q,estimated_atom=mass)
}
