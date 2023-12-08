inde_four_stats<-function(X,Y)
{
  n<-dim(X)[1]
  p<-dim(X)[2]
  q<-dim(Y)[2]
  v=n*(n-3)/2
  
  numeratorZZYS=0
  denominatorZZYS1=0
  denominatorZZYS2=0
  
  numeratorHoeffdingD=0
  denominatorHoeffdingD1=0
  denominatorHoeffdingD2=0
  
  numeratorBKR=0
  denominatorBKR1=0
  denominatorBKR2=0
  
  numeratorTauStar=0
  denominatorTauStar1=0
  denominatorTauStar2=0
  
  for (i in 1:p)
  {
    for (j in 1:q)
    {
      numeratorZZYS=numeratorZZYS+dcov2d(X[,i],Y[,j],type="U")
      numeratorHoeffdingD=numeratorHoeffdingD+HoeffdingD(X[,i],Y[,j])
      numeratorBKR=numeratorBKR+(tStar(X[,i],Y[,j])-12*HoeffdingD(X[,i],Y[,j]))/24
      numeratorTauStar=numeratorTauStar+tStar(X[,i],Y[,j])
    }
  }
  
  for (i in 1:p)
  {
    for (j in 1:p)
    {
      denominatorZZYS1=denominatorZZYS1+dcov2d(X[,i],X[,j],type="U")
      denominatorBKR1=denominatorBKR1+(tStar(X[,i],X[,j])-12*HoeffdingD(X[,i],X[,j]))/24
    }
  }
  
  for (i in 1:q)
  {
    for (j in 1:q)
    {
      denominatorZZYS2=denominatorZZYS2+dcov2d(Y[,i],Y[,j],type="U")
      denominatorBKR2=denominatorBKR2+(tStar(Y[,i],Y[,j])-12*HoeffdingD(Y[,i],Y[,j]))/24
    }
  }
  
  uDcor2=numeratorZZYS/sqrt(denominatorZZYS1*denominatorZZYS2)
  TuDcor2=sqrt(v-1)*uDcor2/sqrt(1-uDcor2^2)
  
  HoeffdingD=numeratorHoeffdingD/sqrt(denominatorBKR1*denominatorBKR2)
  THoeffdingD=sqrt(n*(n-1)/2)*HoeffdingD
  
  BKR=numeratorBKR/sqrt(denominatorBKR1*denominatorBKR2)
  TBKR=sqrt(n*(n-1)/2)*BKR
  
  TauStar=(1/36)*numeratorTauStar/sqrt(denominatorBKR1*denominatorBKR2)
  TTauStar=sqrt(n*(n-1)/2)*TauStar
  
  return(list(TuDcor2=TuDcor2,THoeffdingD=THoeffdingD,TBKR=TBKR,TTauStar=TTauStar,v=v))
}

HoeffdingD=function(x,y)
{
  n=length(y)
  posit=order(y,decreasing=FALSE)
  x=x[posit]
  xM=x%*%matrix(1,1,n)
  IM=(xM<t(xM))
  iV=apply(IM*upper.tri(matrix(0,n,n),diag=T),2,sum)
  iVall=apply(IM,2,sum)
  H1=(sum(iV^2)-sum(iV))/(n*(n-1)*(n-2))
  H2=(sum(iV*iVall*(1:n-1))-sum(iV^2)-sum(iV*(iVall+1:n-3)))/(n*(n-1)*(n-2)*(n-3))
  H3=(sum(iVall*(iVall-1)*(1:n-1)*(1:n-2))-2*n*(n-1)*(n-2)*H1-4*n*(n-1)*(n-2)*(n-3)*H2)/(n*(n-1)*(n-2)*(n-3)*(n-4))
  return(H1-2*H2+H3)
}

###R package
library(Rcpp)
library(MASS)
library(energy)
library(TauStar)



###sample size
n=100
###dimensions 
p=50
q=50
signal=floor(p/3)

###the degree of heteroscedasticity
DD<-c(0,0.25,0.5,0.75,1)

###The empirical power
Mpower=matrix(0,5,5)

for (k in 1:5){
set.seed(2021)

###The repeated time
reptime=500
pvaluesZZYS=rep(0,reptime)
pvaluesSR=rep(0,reptime)
pvaluesHoeffdingD=rep(0,reptime)
pvaluesBKR=rep(0,reptime)
pvaluesTauStar=rep(0,reptime)

####AR structure
delta1=DD[k]
delta2=DD[k]
rho=0.5
Sig11=diag((1:p)^(delta1))
for (i in 1:p){
  for (j in 1:p){
    Sig11[i,j]=rho^{abs(i-j)}*sqrt(Sig11[i,i]*Sig11[j,j])
  }
}
Sig22=diag((1:(q-signal))^(delta2))
for (i in 1:(q-signal)){
  for (j in 1:(q-signal)){
    Sig22[i,j]=rho^{abs(i-j)}*sqrt(Sig22[i,i]*Sig22[j,j])
  }
}


for (i in 1:reptime)
{
  print(i)
  set.seed(10*i)
  ###generate x and y
  x=mvrnorm(n, rep(0, p), Sig11)
  y=matrix(0,n,q)
  ###Linear
  y[,1:signal]=0.5*x[,1:signal]+mvrnorm(n, rep(0, signal), diag(1:signal))
  y[,(signal+1):q]=mvrnorm(n, rep(0,(q-signal)), Sig22)
  ###Nonlinear 
  #y[,1:signal]=x[,1:signal]^2+mvrnorm(n, rep(0, signal), diag(1:signal))
  #y[,(signal+1):q]=mvrnorm(n, rep(0,(q-signal)), Sig22)^2
  
  ####test 
  results=inde_four_stats(x,y)
  ###Hoeffding  D
  pvaluesHoeffdingD[i]=1-pnorm(results$THoeffdingD,0,1)
  ###Blum–Kiefer–Rosenblatt R
  pvaluesBKR[i]=1-pnorm(results$TBKR,0,1)
  ###tauStar
  pvaluesTauStar[i]=1-pnorm(results$TTauStar,0,1)
  ###aggregated distance covariance (Zhu et al. (2020)).
  pvaluesZZYS[i]=1-pt(results$TuDcor2,results$v)
  ###bias-corrected distance correlation (Székely and Rizzo (2013))
  pvaluesSR[i]=dcorT.test(x,y)$p.value
  
  print(c(mean(pvaluesHoeffdingD[1:i]<=0.05),mean(pvaluesBKR[1:i]<=0.05),mean(pvaluesTauStar[1:i]<=0.05),
    mean(pvaluesZZYS[1:i]<=0.05),mean(pvaluesSR[1:i]<=0.05)))
}
Mpower[k,]=c(mean(pvaluesHoeffdingD<=0.05),mean(pvaluesBKR<=0.05),mean(pvaluesTauStar<=0.05),
             mean(pvaluesZZYS<=0.05),mean(pvaluesSR<=0.05))
print(Mpower)
}






