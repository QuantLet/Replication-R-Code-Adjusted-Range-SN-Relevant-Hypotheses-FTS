source('R/primary_simulations.R')
run_validation <- function() {
dir.create('validation',showWarnings=FALSE)
set.seed(4711);x<-matrix(rnorm(73*3),73,3);st<-segment_stats(x,31,6)
# Independent explicit sums verify total variance and rounded score paths.
mx<-colMeans(x[1:31,]);my<-colMeans(x[32:73,]);D<-mx-my
sx<-as.vector(2*sweep(x[1:31,],2,mx)%*%D);sy<-as.vector(2*sweep(x[32:73,],2,my)%*%D)
blocks<-function(u,b){groups<-split(u,ceiling(seq_along(u)/b));sum(vapply(groups,sum,numeric(1))^2)}
se<-blocks(sx,6)/31^2+blocks(sy,6)/42^2
manualpath<-vapply(st$r,function(r)sum(sx[seq_len(floor(31*r))])/31-sum(sy[seq_len(floor(42*r))])/42,numeric(1));manualpath[c(1,41)]<-0
stopifnot(abs(se-st$se2)<1e-12,max(abs(manualpath-st$path))<1e-12)
# Independent scalar construction tests two-segment covariance and endpoint sharing.
z<-matrix(rnorm(80*27),27,80);ref<-segment_refs(st,z);brute<-matrix(NA_real_,27,2)
for(i in 1:27){g<-numeric(41);E<-0;for(j in 1:2){w<-c(0,cumsum(z[i,((j-1)*40+1):(j*40)]*sqrt(diff(st$v[[j]]))));fac<-c(1,-1)[j]*sqrt(st$weights[j]);E<-E+fac*w[41];g<-g+fac*(w-st$rx[[j]]*w[41])};brute[i,]<-E/c(diff(range(g)),sqrt(mean(g[2:40]^2)))}
stopifnot(max(abs(brute-ref))<1e-12,anyDuplicated(st$rx[[1]])>0)
zero<-segment_stats(matrix(1,70,3),35,5);stopifnot(!zero$valid)
atom<-matrix(rep(c(-2,2),100),ncol=2);stopifnot(!reference_check(atom)$valid)
# Data scale changes multiply squared distance, variances, and denominators appropriately.
st2<-segment_stats(3*x,31,6);stopifnot(abs(st2$d-9*st$d)<1e-10,abs(st2$se2-81*st$se2)<1e-10,max(abs(st2$den-9*st$den))<1e-10)
# Deliberately nonfinite total variance must not return a nominal inference.
y<-matrix(0,70,3);y[1,1]<-Inf;invalid<-try(segment_stats(y,35,5),silent=TRUE);stopifnot(inherits(invalid,'try-error')||!invalid$valid)
checks<-data.frame(check=c('segment variance','rounded path','joint segment Gaussian reference','zero total variance','quantile atom detection','scale consistency'),passed=TRUE)
write.csv(checks,file.path('validation','checks.csv'),row.names=FALSE)
cat('All segment, reference, degeneracy and scaling checks passed.\n')
invisible(checks)
}
