Wildbootfctn<-function(d){
  n<-dim(d)[1]
  kep<-dim(d)[2]-1
  INFL2<-d[,1]
  INFL.ep<-d[,2:(1+kep)]
  INFL.ep<-as.matrix(INFL.ep)
  COV.ep<-rep(0,kep)
  V.ep<-matrix(0,kep,kep)

  for(jj in 1:kep){
    COV.ep[jj]<-mean(INFL.ep[,jj]*INFL2)/n
  }
  for(jj in 1:kep){
    for(kk in jj:kep){
      V.ep[jj,kk]<-V.ep[kk,jj]<-mean(INFL.ep[,jj]*INFL.ep[,kk])/n
    }
  }
  VE2<-mean(INFL2^2)/n

  return(list(VE2=VE2,V.ep=V.ep,COV.ep=COV.ep))
}


