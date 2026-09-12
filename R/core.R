# Functional tests in coefficient space; rows are observations.
profile <- function(r, type) switch(type, uniform=r, r5=r^5,
 step=ifelse(r<=.5,r/4,.125+1.75*(r-.5)),
 plateau=ifelse(r<=.3,5*r/3,ifelse(r<=.7,.5,.5+5*(r-.7)/3)),stop('Unknown profile'))
trapz <- function(y,x) sum(diff(x)*(head(y,-1)+tail(y,-1))/2)
csrows <- function(x) rbind(0,apply(as.matrix(x),2,cumsum))
estimate_split <- function(x,trim=.15) {
 n<-nrow(x);ks<-seq(max(1,ceiling(trim*n)),min(n-1,floor((1-trim)*n)))
 cs<-csrows(x);r<-ks/n;g<-cs[ks+1,,drop=FALSE]-outer(r,cs[n+1,])
 ks[which.max(rowSums(g*g)/(n*r*(1-r)))]
}
block_cumulative <- function(scores,r,block) {
 n<-length(scores);edges<-unique(c(seq(0,n,by=block),n));cs<-c(0,cumsum(scores))
 approx(edges/n,c(0,cumsum(diff(cs[edges+1])^2)),xout=r,rule=2)$y
}
cp_stats <- function(x,k,b,r=seq(0,1,length.out=41),delta=1) {
 n<-nrow(x);xx<-x[seq_len(k),,drop=FALSE];yy<-x[(k+1):n,,drop=FALSE]
 mx<-colMeans(xx);my<-colMeans(yy);D<-mx-my
 sx<-as.vector(2*sweep(xx,2,mx)%*%D);sy<-as.vector(2*sweep(yy,2,my)%*%D)
 path<-c(0,cumsum(sx))[floor(k*r)+1]/k-c(0,cumsum(sy))[floor((n-k)*r)+1]/(n-k)
 path[c(1,length(r))]<-0;den<-c(AR=diff(range(path)),Q=sqrt(mean(path[-c(1,length(r))]^2)))
 qb<-block_cumulative(sx,r,b)/k^2+block_cumulative(sy,r,b)/(n-k)^2
 qr<-block_cumulative(sx,r,1)/k^2+block_cumulative(sy,r,1)/(n-k)^2
 list(d=sum(D^2),den=den,piv=(sum(D^2)-delta^2)/pmax(den,1e-15),v=qb/tail(qb,1),raw_v=qr/tail(qr,1),path=path,direction=D,scores=c(sx,sy))
}
brownian_refs <- function(v,r=seq(0,1,length.out=length(v)),B=999,z=NULL) {
 if(is.null(z))z<-matrix(rnorm(B*(length(r)-1)),B)
 B<-nrow(z);w<-cbind(0,t(apply(sweep(z,2,sqrt(pmax(diff(v),0)),'*'),1,cumsum)))
 end<-w[,ncol(w)];g<-w-outer(end,r)
 cbind(AR=end/apply(g,1,function(a)diff(range(a))),Q=end/sqrt(rowMeans(g[,-c(1,ncol(g)),drop=FALSE]^2)))
}
# Explicit comparison avoids partial matching of primitive operator names.
rank_p <- function(ref,piv) (1+colSums(sweep(ref,2,piv,`>=`)))/(nrow(ref)+1)
errors <- function(n,dep='ar0.4',tail='gaussian',K=3,raw=NULL) {
 burn<-300
 if(is.null(raw)) {
  raw<-matrix(rnorm((n+burn)*K),n+burn,K)
  if(tail=='orthogonal_t1.6')raw[,2]<-rt(n+burn,1.6)
  if(tail=='heavy')raw[,-1]<-matrix(rt((n+burn)*(K-1),1.6),n+burn,K-1)
 }
 if(K==3)raw<-sweep(raw,2,c(1,.3,.2),'*')
 if(startsWith(dep,'ar')) {
  rho<-as.numeric(sub('ar','',dep));e<-apply(raw,2,function(z)as.numeric(filter((1-rho)*z,rho,method='recursive',init=0)))
 } else if(dep=='ma0.7') e<-(raw+.7*rbind(0,head(raw,-1)))/1.7
 else if(dep=='iid')e<-raw
 else if(dep=='switch') {
  e<-raw;prev<-rep(0,K)
  for(i in seq_len(n+burn)){rho<-if(i<=burn+floor(.6*n)).2 else .8;prev<-rho*prev+sqrt(1-rho^2)*raw[i,];e[i,]<-prev}
  e<-e/sqrt(4.5)
 }else stop('Unknown filter')
 e[(burn+1):(n+burn),,drop=FALSE]
}
loading_stats <- function(e,a,d) {
 n<-nrow(e);w<-e*a;m<-colMeans(w);mu<-m;mu[1]<-mu[1]+sqrt(d)
 sc<-as.vector(2*w%*%mu)-2*a^2*sum(mu*m);g<-c(0,cumsum(sc))/n;g[n+1]<-0
 v<-c(0,cumsum(a^2))/n;den<-c(AR=diff(range(g)),Q=sqrt(trapz(g*g,v)),calendarQ=sqrt(trapz(g*g,seq(0,1,length.out=n+1))))
 list(piv=(sum(mu^2)-1)/pmax(den,1e-15),d=sum(mu^2),path=g,den=den)
}
wild_refs <- function(x,k,b,B=4999,delta=1,trim=.2,r=seq(0,1,length.out=41),weights=NULL) {
 n<-nrow(x);ma<-colMeans(x[1:k,,drop=FALSE]);mz<-colMeans(x[(k+1):n,,drop=FALSE]);D<-ma-mz;u<-D/sqrt(sum(D^2))
 res<-x-rbind(matrix(ma,k,ncol(x),byrow=TRUE),matrix(mz,n-k,ncol(x),byrow=TRUE))
 null<-outer(c(rep(.5*delta,k),rep(-.5*delta,n-k)),u)
 out<-matrix(NA_real_,B,2,dimnames=list(NULL,c('AR','Q')))
 for(j in seq_len(B)) {
  if(is.null(weights)){z<-rnorm(n+b-1);cs<-c(0,cumsum(z));w<-(cs[(b+1):(n+b)]-cs[1:n])/sqrt(b)}else w<-weights[j,]
  xb<-null+res*w;kb<-estimate_split(xb,trim);out[j,]<-cp_stats(xb,kb,b,r,delta)$piv
 };out
}
