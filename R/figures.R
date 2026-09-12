suppressPackageStartupMessages({library(ggplot2);library(patchwork)})
source('R/primary_simulations.R')
source('R/data_preparation.R')

paper_colours <- c(AR='#0072B2',Q='#D55E00',Wald='#009E73')
paper_lines <- c(AR='solid',Q='dashed',Wald='dotdash')
paper_shapes <- c(AR=16,Q=0,Wald=17)
paper_theme <- theme_bw(base_size=11)+theme(
 panel.grid.minor=element_blank(),panel.grid.major=element_line(colour='grey92',linewidth=.25),
 legend.position='bottom',legend.title=element_blank(),legend.key.width=grid::unit(1.5,'cm'),
 strip.background=element_rect(fill='grey96'),plot.title=element_text(size=11))

save_figure <- function(plot,name,width=6.5,height=4.8) {
 dir.create('figures',showWarnings=FALSE)
 ggsave(file.path('figures',paste0(name,'.pdf')),plot,width=width,height=height,device='pdf')
 ggsave(file.path('figures',paste0(name,'.png')),plot,width=width,height=height,dpi=180,bg='white')
}

make_figures <- function() {
 required<-file.path(results_dir,c('cp_results.csv','adjusted_power.csv','empirical_tests.csv','empirical_shapes.csv'))
 if(any(!file.exists(required)))stop('Run run_simulations.R and run_empirical.R before making figures.')

 cp<-read.csv(file.path(results_dir,'cp_results.csv'))
 fitted<-subset(cp,comparison%in%c('fitted_AR','fitted_Q','fitted_Wald'))
 fitted$method<-sub('fitted_','',fitted$comparison)
 adj<-read.csv(file.path(results_dir,'adjusted_power.csv'))
 power<-do.call(rbind,lapply(names(paper_colours),function(method)
  data.frame(adj[,c('n','profile','dep','tail')],method=method,rate=adj[[method]],panel='Size-adjusted power')))
 size<-transform(subset(fitted,tail=='gaussian'&dep=='ar0.4'&state=='boundary',
  select=c(n,profile,dep,tail,method,rate)),panel='Boundary rejection')
 x<-rbind(size,subset(power,tail=='gaussian'&dep=='ar0.4'))
 x$panel<-factor(x$panel,c('Boundary rejection','Size-adjusted power'))
 x$profile<-factor(x$profile,c('uniform','r5'))
 nominal<-data.frame(panel=factor('Boundary rejection',levels=levels(x$panel)),
  profile=factor(c('uniform','r5'),levels=levels(x$profile)))
 p1<-ggplot(x,aes(n,100*rate,colour=method,linetype=method,shape=method))+
  geom_line()+geom_point(size=2)+geom_hline(data=nominal,aes(yintercept=5),linetype='dotted',colour='grey45')+
  facet_grid(panel~profile,scales='free_y',labeller=labeller(
   profile=as_labeller(c(uniform='Uniform',r5='r^5'),label_parsed),panel=label_value))+
  scale_x_log10(breaks=c(70,200,1000),labels=c('70','200','1000'))+
  scale_colour_manual(values=paper_colours)+scale_linetype_manual(values=paper_lines)+
  scale_shape_manual(values=paper_shapes)+labs(x='Number of curves N',y='Percent')+paper_theme
 save_figure(p1,'unknown_break_comparison')

 download_sources();audited<-audit_data();series<-construct_series(audited$curves)
 dates<-series$dates;xx<-series$x;k<-estimate_split(xx/sqrt(48),.2)
 daily<-data.frame(date=dates,y=rowMeans(xx))
 p2a<-ggplot(daily,aes(date,y))+geom_hline(yintercept=0,colour='grey70')+
  geom_line(colour='#0072B2')+geom_point(colour='#0072B2',size=.7)+
  geom_vline(xintercept=dates[k+1],linetype='dashed')+
  scale_x_date(date_breaks='1 month',date_labels='%b')+
  labs(x='2020 weekday',y='Daily average\ndifference (GW)',title='(A) Same-weekday annual differences')+paper_theme
 means<-rbind(data.frame(hour=((1:48)-.5)/2,y=colMeans(xx[1:k,]),series='Before 24 March'),
  data.frame(hour=((1:48)-.5)/2,y=colMeans(xx[(k+1):nrow(xx),]),series='From 24 March'))
 p2b<-ggplot(means,aes(hour,y,colour=series,linetype=series))+geom_line()+
  scale_colour_manual(values=c('#0072B2','#D55E00'))+scale_linetype_manual(values=c('solid','dashed'))+
  labs(x='Local hour',y='Mean difference (GW)',title='(B) Segment mean curves')+paper_theme
 shapes<-read.csv(file.path(results_dir,'empirical_shapes.csv'))
 shapes<-subset(shapes,case%in%c('primary','two_year_reference'))
 shapes$case<-factor(shapes$case,c('primary','two_year_reference'),c('One-year reference','Two-year reference'))
 p2c<-ggplot(shapes,aes((period-.5)/2,shape,colour=case,linetype=case))+geom_hline(yintercept=0,colour='grey65')+
  geom_line()+scale_colour_manual(values=c('#0072B2','#D55E00'))+scale_linetype_manual(values=c('solid','dashed'))+
  labs(x='Local hour',y='Centered change (GW)',title='(C) Shape component')+paper_theme
 p2<-(p2a/p2b/p2c)&theme(legend.margin=margin(0,0,0,0),legend.box.spacing=grid::unit(1,'pt'),plot.margin=margin(3,3,3,3))
 save_figure(p2,'demand_decomposition',height=5.7)

 tests<-subset(read.csv(file.path(results_dir,'empirical_tests.csv')),
  case%in%c('primary','two_year_reference')&block==5)
 tests$Reference<-factor(tests$case,c('primary','two_year_reference'),c('One-year reference','Two-year reference'))
 p3<-ggplot(tests,aes(threshold,pvalue,colour=method,linetype=method,shape=method))+
  geom_hline(yintercept=.05,linetype='dotted',colour='grey40')+geom_line()+geom_point(size=1.8)+
  facet_wrap(~Reference,nrow=1)+scale_y_log10()+scale_colour_manual(values=paper_colours)+
  scale_linetype_manual(values=paper_lines)+scale_shape_manual(values=paper_shapes)+
  labs(x='RMS relevance threshold (GW)',y='Reference p-value')+paper_theme
 save_figure(p3,'empirical_thresholds')
 invisible(c('unknown_break_comparison','demand_decomposition','empirical_thresholds'))
}
