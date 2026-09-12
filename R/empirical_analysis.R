source('R/primary_simulations.R')
source('R/data_preparation.R')
run_empirical_analysis <- function(B=49999) {
 download_sources();dir.create(results_dir,recursive=TRUE,showWarnings=FALSE)
 a<-audit_data();write.csv(a$audit,file.path(results_dir,'data_audit.csv'),row.names=FALSE)
 rows<-list();curvesout<-list();summary<-list()
 dates0<-seq(as.Date('2020-02-04'),as.Date('2020-05-11'),by='day');dates0<-dates0[!as.POSIXlt(dates0)$wday%in%c(0,6)]
 # England/Wales bank holidays in the date ranges being matched. GOV.UK list.
 holidays<-as.Date(c('2019-04-19','2019-04-22','2019-05-06','2020-04-10','2020-04-13','2020-05-08'))
 cases<-c('primary','two_year_reference','placebo_2019','placebo_2018','window_earlier','window_later','exclude_holidays')
 writeLines(c('Sensitivity rules fixed before this recomputation:', 'Shift both window endpoints by -7 or +7 calendar days, retaining weekdays.', 'Exclude a pair if either 2020 date or 364-day reference date is an England/Wales bank holiday.', 'These are exploratory sensitivity checks, not a newly preregistered confirmatory analysis.','Holiday source: https://www.gov.uk/bank-holidays'),file.path(results_dir,'empirical_sensitivity_rules.txt'))
 for(cc in cases){off<-switch(cc,placebo_2019=-364,placebo_2018=-728,window_earlier=-7,window_later=7,0);dates<-dates0+off
  if(cc=='exclude_holidays')dates<-dates[!(dates%in%holidays|(dates-364)%in%holidays)]
  xx<-t(vapply(as.character(dates),function(day){d<-as.Date(day);refs<-as.character(d-c(364,if(cc=='two_year_reference')728));stopifnot(all(c(day,refs)%in%names(a$curves)));a$curves[[day]]-Reduce('+',a$curves[refs])/length(refs)},numeric(48)))
  x<-xx/sqrt(48);k<-estimate_split(x,.2);D<-colMeans(xx[1:k,])-colMeans(xx[(k+1):nrow(xx),]);shape<-D-mean(D)
  summary[[cc]]<-data.frame(case=cc,n=nrow(x),split=k,date=as.character(dates[k+1]),RMS=sqrt(mean(D^2)),level=mean(D),shape_RMS=sqrt(mean(shape^2)),level_share=mean(D)^2/mean(D^2))
  curvesout[[cc]]<-data.frame(case=cc,period=1:48,contrast=D,shape=shape)
  for(b in c(3,5,7,10)){
   st<-segment_stats(x,k,b);set.seed(config_seed(list(cc,b),'empirical_revision_reference'));z<-matrix(rnorm(B*80),B,80);ref<-segment_refs(st,z);q<-ref_quantile(ref)
   lows<-c(st$d-q*st$den,st$d-qnorm(.95)*sqrt(st$se2))
   for(delta in c(.5,1,1.5,2,2.5,3,4,5)){
    piv<-(st$d-delta^2)/st$den;pv<-c(rank_p(ref,piv),pnorm((st$d-delta^2)/sqrt(st$se2),lower.tail=FALSE))
    rows[[length(rows)+1]]<-data.frame(case=cc,n=nrow(x),split=k,date=as.character(dates[k+1]),block=b,method=c('AR','Q','Wald'),threshold=delta,RMS=sqrt(st$d),lower_RMS=sqrt(pmax(lows,0)),pvalue=pv,B=B)
   }
  };cat('Empirical case',cc,'completed\n')
 };write.csv(do.call(rbind,rows),file.path(results_dir,'empirical_tests.csv'),row.names=FALSE);write.csv(do.call(rbind,summary),file.path(results_dir,'empirical_summary.csv'),row.names=FALSE);write.csv(do.call(rbind,curvesout),file.path(results_dir,'empirical_shapes.csv'),row.names=FALSE)
}
