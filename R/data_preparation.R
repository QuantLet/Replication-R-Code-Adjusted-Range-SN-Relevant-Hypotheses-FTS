source('R/core.R')
# Public NESO files used by the empirical analysis. Checksums fix the sample.
neso_sources <- data.frame(
 year=2017:2021,
 file=paste0('demanddata_',2017:2021,'.csv'),
 url=c(
  'https://api.neso.energy/dataset/8f2fe0af-871c-488d-8bad-960426f24601/resource/2f0f75b8-39c5-46ff-a914-ae38088ed022/download/demanddata_2017.csv',
  'https://api.neso.energy/dataset/8f2fe0af-871c-488d-8bad-960426f24601/resource/fcb12133-0db0-4f27-a4a5-1669fd9f6d33/download/demanddata_2018.csv',
  'https://api.neso.energy/dataset/8f2fe0af-871c-488d-8bad-960426f24601/resource/dd9de980-d724-415a-b344-d8ae11321432/download/demanddata_2019.csv',
  'https://api.neso.energy/dataset/8f2fe0af-871c-488d-8bad-960426f24601/resource/33ba6857-2a55-479f-9308-e5c4c53d4381/download/demanddata_2020.csv',
  'https://api.neso.energy/dataset/8f2fe0af-871c-488d-8bad-960426f24601/resource/18c69c42-f20d-46f0-84e9-e279045befc6/download/demanddata_2021.csv'),
 sha256=c(
  'ad6c47de0264aa1183336c7ca438631f90056c41f438c86f34da147a51351de4',
  'a2e19b4eb4157a190d1e7f2c946ee92c4be2ae196f72d55545c870399d3efb2f',
  '0084d546dcdba0af21b7414f1a1b0ab308a9de8460884aa8e550fd0e858b38da',
  '7121b6fbf317862febe364f6b702184836b174327bf503ac143de1191f61f141',
  'b9f7285fecec48b05244c5bd3206c14d2409ab6da9cfc5c3b976104d5d5734ef'),
 stringsAsFactors=FALSE)

download_sources <- function() {
 dir.create('data/raw',recursive=TRUE,showWarnings=FALSE)
 for(i in seq_len(nrow(neso_sources))){z<-neso_sources[i,];dest<-file.path('data/raw',z$file)
  if(!file.exists(dest))download.file(z$url,dest,mode='wb',quiet=TRUE)
  stopifnot(identical(digest::digest(file=dest,algo='sha256'),z$sha256))
 };invisible(TRUE)
}
audit_data <- function() {
 oldlocale<-Sys.getlocale('LC_TIME');on.exit(Sys.setlocale('LC_TIME',oldlocale));Sys.setlocale('LC_TIME','C')
 curves<-list();audit<-list()
 for(f in list.files('data/raw',pattern='^demanddata_.*csv$',full.names=TRUE)){
  a<-read.csv(f);year<-as.integer(sub('.*([0-9]{4})\\.csv$','\\1',f));a$date<-as.Date(a$SETTLEMENT_DATE,'%d-%b-%Y')
  expected<-seq(as.Date(paste0(year,'-01-01')),as.Date(paste0(year,'-12-31')),by='day')
  stopifnot(!anyNA(a$date),setequal(expected,a$date),!anyDuplicated(a[c('date','SETTLEMENT_PERIOD')]),all(is.finite(a$ND)&a$ND>0))
  days<-split(a,as.character(a$date));atypical<-0
  for(day in names(days)){z<-days[[day]];z<-z[order(z$SETTLEMENT_PERIOD),];stopifnot(identical(as.integer(z$SETTLEMENT_PERIOD),seq_len(nrow(z))))
   if(nrow(z)==48)curves[[day]]<-z$ND/1000 else{stopifnot(as.POSIXlt(as.Date(day))$wday==0,nrow(z)%in%c(46,50));atypical<-atypical+1}
  };audit[[length(audit)+1]]<-data.frame(year=year,observations=nrow(a),dates=length(days),clock_change_days=atypical)
 };list(curves=curves,audit=do.call(rbind,audit))
}
construct_series <- function(curves,offset=0,two_year=FALSE) {
 dates<-seq(as.Date('2020-02-04'),as.Date('2020-05-11'),by='day');dates<-dates[!as.POSIXlt(dates)$wday%in%c(0,6)]-offset
 stopifnot(length(dates)==70)
 x<-t(vapply(as.character(dates),function(day){d<-as.Date(day);refs<-as.character(d-c(364,if(two_year)728));stopifnot(all(c(day,refs)%in%names(curves)));curves[[day]]-Reduce(`+`,curves[refs])/length(refs)},numeric(48)))
 list(dates=dates,x=x)
}
