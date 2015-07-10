"mq" <- function(x,lag){
# Compute multivariate Ljung-Box test statistics
#
LB <- matrix(nrow=lag,ncol=3)
if(!is.matrix(x))x=as.matrix(x)
nr=nrow(x)
nc=ncol(x)
g0=var(x)
ginv=solve(g0)
qm=0.0
#print("m,       Q(m) and  p-value:")
df = 0
for (i in 1:lag){
  x1=x[(i+1):nr,]
  x2=x[1:(nr-i),]
  g = cov(x1,x2)
  g = g*(nr-i-1)/(nr-1)
  h=t(g)%*%ginv%*%g%*%ginv
  qm=qm+nr*nr*sum(diag(h))/(nr-i)
  df=df+nc*nc
  pv=1-pchisq(qm,df)
   LB[i,]<- c(i,qm,pv)
  
 # print(c(i,qm,pv))
}
return(LB)
}