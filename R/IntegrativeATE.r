#' Integrative Average Treatment Effect (IntegrativeATE)
#'
#' Implements integrative analyses for the average treatment effect combining main data with unmeasured confounder
#' and validation data with supplementary information on these confounders.
#' @param I the vector of the binary indicator of the validation sample membership; i.e., 1 if the unit belongs to the validation sample, and 0 otherwise  (n x 1)
#' @param x the matrix of confounders fully observed  (n x dim(x))
#' @param u the matrix of confounders observed only for the validation sample; and NA otherwise (n x dim(u))
#' @param y the vector of outcome (n x 1)
#' @param A the vector of binary treatment (n x 1)
#' @param method_val the estimation method for ATE on the validation sample
#'
#' select 1 method from c("reg","ipw","aipw","matching")
#'
#' \code{"reg"}: a linear regression imputation estimator of the ATE
#'
#' \code{"ipw"}: the inverse probability of treatment weighting estimator of the ATE, where the propensity score follows a logistic regression model
#'
#' \code{"aipw"}: the augmented inverse probability of treatment weighting estimator (Lunceford Davidian, 2004) of the ATE
#'
#' \code{"matching"}: the matching estimator (Abadie and Imbens, 2006) of the ATE, using matching based on (X,U) with replacement with the number of matches fixed at M.
#'
#' @param method_ep  the error prone estimation method applied to the validation sample and the main sample
#'
#' select >=1 methods from c("none","reg","ipw","aipw","matching")
#'
#' \code{"none"}: return the initial estimator based soly on the validation sample
#'
#' @param nboot the number of bootstrap samples; if nboot=0, then return only the variance estimator based on the asymptotic result
#'
#' @return
#' \itemize{
#'
#' \item \code{est}: estimate of the ATE
#'
#' \item \code{ve}:  variance estimate for \code{est} based on the asymptotic result
#'
#' \item \code{ve_boot}: variance estimate for \code{est} based on wild bootstrap
#' }
#'
#' @details
#' Under the unconfoundedness assumption
#' with completely observed confounders, the smaller validation data allow for constructing consistent estimators
#' for causal effects, but the big main data can only give error-prone estimators in general.
#'
#' The integrative estimator leverages the information in the big main data to improve the estimation
#' efficiencies yet preserve the consistencies of the initial estimators based solely on the validation data.
#' Specifically, it uses the difference of the error-prone estimators applied to the main data and the validation data,
#' which is consistent for zero assuming that the main data and the validation data are representative of the same target population.
#'
#' The framework applies to asymptotically normal estimators, including the commonly-used regression imputation,
#' weighting, and matching estimators.
#'
#' @import MASS ks Matching stats
#'
#'@references
#'
#'Yang, S. and Ding, P. (2018). Combining multiple observational data sources to estimate causal effects.
#' \url{https://arxiv.org/abs/1801.00802}
#'
#'
#' @examples
#'
#' n<-1000  # the (combined) main sample size (Samples 1 and 2)
#' n2<-500  # the validation sample size (Sample 2)
#'
#' ## generate covariates (x,u)
#'
#' x<-runif(n,0,2)
#' u<-2*sin(x)+runif(n,0,1)
#'
#' ## generate treatment A
#'
#' lps<-cbind(1,x,u)%*%c(0,-1,1)
#' ps<-exp(lps)/(1+exp(lps))
#' A<-rbinom(n,1,ps)
#'
#' ## generate potential and observed outcomes
#'
#' loc.a1<-which(A==1)
#' loc.a0<-which(A==0)
#' y1<- -x+2*u+rnorm(n,0,1)
#' y0<- -x+1*u+rnorm(n,0,1)
#' y<-y1*A+y0*(1-A)
#' true<-mean(y1-y0)
#'
#' ## generate indicator of membership of the validation sample (Sample 2)
#'
#' I<-rep(1,n)
#' I[((1+(n-n2)):n)]<-0
#'
#' ## u is not observed for Sample 1
#' loc1<-which(I==1)
#' u[loc1]<-NA
#'
#'
#' method_val <-c("reg")
#' method_ep<-c("reg","ipw")
#'
#' true
#' out<-IntegrativeCI::IntegrativeATE(I,x,u,y,A,method_val,method_ep,nboot=50)
#' out$est
#' out$ve
#' out$ve_boot
#'
#'
#' @export

IntegrativeATE<-function(I,x,u,y,A,method_val,method_ep,nboot){

  ## check argument

  method_val<-method_val[which(method_val%in%c("reg","ipw","aipw","matching")==1)]
  method_ep <-method_ep[which(method_ep%in%c("none","reg","ipw","aipw","matching")==1)]
  if(length(method_val)==0){return(cat("Error: specify a valid argument for method_val"))}
  if(length(method_val)>1) {return(cat("Error: allow only for one estimation method based on the validation sample; check method_val"))}
  if(length(method_ep)==0) {return(cat("Error: specify a valid argument for method_ep"))}

  loc1<-which(I==1)
  loc2<-which(I==0)
  n<-length(y)     # the (combined) main sample size (Samples 1 and 2)
  n2<-length(loc2) # the validation sample size (Sample 2)

  ## extract Sample 2 data

  x<-as.matrix(x)
  u<-as.matrix(u)

  AA2<-A[loc2]
  yy2<-y[loc2]
  xx2<-x[loc2,]
  uu2<-u[loc2,]
  xx2<-as.matrix(xx2)
  uu2<-as.matrix(uu2)

  kep<-length(method_ep)

  if( ("none" %in%method_ep)&(kep==1) ){
    kep=0
    EST.ep<-0
    COV.ep<-0
    V.ep<-1
    INFL.ep<-rep(0,n)
  }else{
    EST.ep<-rep(0,kep)
    COV.ep<-rep(0,kep)
    V.ep<-matrix(0,kep,kep)
    INFL.ep<-matrix(0,kep,n)
  }


  if(method_val=="reg"){
    lm.y<- yy2
    lm.x<- cbind(xx2,uu2)
    lm.x<- as.matrix(lm.x)
    dimxu2<-dim(lm.x)[2]
    lm.out1<-lm(lm.y[which(AA2==1)]~lm.x[which(AA2==1),])
    reg2.mu1<-cbind(1,lm.x)%*%lm.out1$coefficients
    lm.out0<-lm(lm.y[which(AA2==0)]~lm.x[which(AA2==0),])
    reg2.mu0<-cbind(1,lm.x)%*%lm.out0$coefficients
    reg2<-mean(reg2.mu1-reg2.mu0)
    B<-apply(cbind(1,lm.x),2,mean)
    B1mat<-cbind(AA2,lm.x*AA2)
    B1<-matrix(0,1+dimxu2,1+dimxu2)
    for(jj in 1:(1+dimxu2)){
      for(kk in jj:(1+dimxu2)){
        B1[jj,kk]<-B1[kk,jj]<-mean(B1mat[,jj]*B1mat[,kk])
      }
    }
    B0mat<-cbind(1-AA2,lm.x*(1-AA2))
    B0<-matrix(0,1+dimxu2,1+dimxu2)
    for(jj in 1:(1+dimxu2)){
      for(kk in jj:(1+dimxu2)){
        B0[jj,kk]<-B0[kk,jj]<-mean(B0mat[,jj]*B0mat[,kk])
      }
    }
    b1<-B%*%ginv(B1)
    b0<-B%*%ginv(B0)
    reg2.resid1<-reg2.resid0<-rep(0,n2)
    reg2.resid1[AA2==1]<-(lm.out1$resid)
    reg2.resid0[AA2==0]<-(lm.out0$resid)
    reg2.adj<-apply( B1mat*matrix(b1,n2,1+dimxu2,byrow=TRUE),1,sum )*AA2*reg2.resid1-
      apply( B0mat*matrix(b0,n2,1+dimxu2,byrow=TRUE),1,sum )*(1-AA2)*reg2.resid0
    infl2.reg<-rep(0,n)
    infl2.reg[loc2]<-(( reg2.mu1-reg2.mu0)+reg2.adj-mean(( reg2.mu1-reg2.mu0)))/n2*n
    ve_reg2<-mean(((reg2.mu1-reg2.mu0+reg2.adj)- mean((reg2.mu1-reg2.mu0+reg2.adj)))*
                    (( reg2.mu1-reg2.mu0)+reg2.adj -mean((reg2.mu1-reg2.mu0+reg2.adj))))/n2
    EST2<-reg2
    VE2<-ve_reg2
    INFL2<-infl2.reg
  }


  if(method_val=="ipw"){
    Rxu<-cbind(xx2,uu2)
    Rxu<-as.matrix(Rxu)
    dimxu2<-dim(Rxu)[2]
    glm.out<-glm(AA2~Rxu,family="binomial")
    pshat2<-glm.out$fitted.values
    psw2<-mean(yy2*AA2/pshat2)-mean(yy2*(1-AA2)/(1-pshat2))
    ## using Lunceford Davidian (2004)'s method for variance estimation of psw2
    H<-c(mean( AA2*yy2*(1-pshat2)/pshat2+(1-AA2)*yy2*(pshat2)/(1-pshat2) ),
         apply( (AA2*yy2*(1-pshat2)/pshat2+(1-AA2)*yy2*(pshat2)/(1-pshat2))*Rxu , 2 , mean ) )

    E<-matrix(0,1+dimxu2,1+dimxu2)
    forE<-matrix(1,1+dimxu2,n2)
    forE[2:(1+dimxu2),]<-t(Rxu)
    for(jj in 1:(1+dimxu2)){
      for(kk in jj:(1+dimxu2)){
        E[jj,kk]<-E[kk,jj]<-mean(pshat2*(1-pshat2)*forE[jj,]*forE[kk,])
      }
    }
    alpha2<- -(AA2-pshat2) * H%*%ginv(E)%*%forE
    infl2.ipw<-rep(0,n);
    psi2_psw<-yy2*AA2/pshat2-yy2*(1-AA2)/(1-pshat2)
    infl2.ipw[loc2]<- ( psi2_psw+alpha2-psw2 )/n2*n
    ve_psw2<-mean((psi2_psw+alpha2-psw2)*(psi2_psw+alpha2-psw2))/n2
    EST2<-psw2
    VE2<-ve_psw2
    INFL2<-infl2.ipw
  }


  if(method_val=="aipw"){
    lm.y<- yy2
    lm.x<- cbind(xx2,uu2)
    lm.out1<-lm(lm.y[which(AA2==1)]~lm.x[which(AA2==1),])
    mu1<-cbind(1,lm.x)%*%lm.out1$coefficients
    lm.out0<-lm(lm.y[which(AA2==0)]~lm.x[which(AA2==0),])
    mu0<-cbind(1,lm.x)%*%lm.out0$coefficients
    Rxu<-cbind(xx2,uu2)
    glm.out<-glm(AA2~Rxu,family="binomial")
    pshat2<-glm.out$fitted.values
    psi2.aug<-yy2*AA2/pshat2-yy2*(1-AA2)/(1-pshat2)+( mu1*(1-AA2/pshat2)-mu0*(1-(1-AA2)/(1-pshat2)) )
    aipw2<-mean(psi2.aug)
    infl2.aug<-rep(0,n);
    infl2.aug[loc2]<- (psi2.aug-aipw2)/n2*n
    ve_aipw2<- mean((psi2.aug-aipw2)*(psi2.aug-aipw2))/n2
    EST2<-aipw2
    VE2<-ve_aipw2
    INFL2<-infl2.aug
  }


  if(method_val=="matching"){
    lm.y<- yy2
    lm.x<- cbind(xx2,uu2)
    lm.out1<-lm(lm.y[which(AA2==1)]~lm.x[which(AA2==1),])
    mu1<-cbind(1,lm.x)%*%lm.out1$coefficients
    lm.out0<-lm(lm.y[which(AA2==0)]~lm.x[which(AA2==0),])
    mu0<-cbind(1,lm.x)%*%lm.out0$coefficients
    ## NNI for Y(1)
    thisps<-cbind(xx2,uu2)
    loc.o<-which(AA2==1)
    Di.nni<-rep(0,length=n2)
    out1<-Match(Y=yy2,Tr=1-AA2,X=thisps,distance.tolerance=0,ties=FALSE, Weight=2,M=1)
    mdata1<-out1$mdata
    fromto<-loc.o
    Di.nni[fromto]<-table(factor(out1$index.control,levels=fromto))
    Di.nni[loc.o]<-Di.nni[loc.o]+1
    est.nni<-sum(yy2*Di.nni,na.rm=TRUE)/n2
    Ey<-mu1
    bias.nni<-sum(Ey[out1$index.control]-Ey[out1$index.treated])/n2
    est1.nni<-est.nni-bias.nni
    truemboot1<-(yy2-Ey)*Di.nni+Ey
    ## NNI for Y(0)
    loc.o<-which(AA2==0)
    Di.nni<-rep(0,length=n2)
    out1<-Matching::Match(Y=yy2,Tr=AA2,X=thisps,distance.tolerance=0,ties=FALSE,Weight=2,M=1)
    mdata1<-out1$mdata
    fromto<-loc.o
    Di.nni[fromto]<-table(factor(out1$index.control,levels=fromto))
    Di.nni[loc.o]<-Di.nni[loc.o]+1
    est.nni<-sum(yy2*Di.nni,na.rm=TRUE)/n2
    Ey<-mu0
    bias.nni<-sum(Ey[out1$index.control]-Ey[out1$index.treated])/n2
    est0.nni<-est.nni-bias.nni
    truemboot0<-(yy2-Ey)*Di.nni+(Ey)
    truemat<-est1.nni-est0.nni
    truemboot<-truemboot1-truemboot0
    mat2<-truemat
    infl2.mat<-rep(0,n);
    infl2.mat[loc2]<- ( truemboot1-truemboot0 -mean(truemboot1-truemboot0))/n2*n
    infl2.matforboot<-infl2.mat
    ve_mat2<-mean((truemboot1-truemboot0-mat2)*(truemboot1-truemboot0-mat2))/n2
    EST2<-mat2
    VE2<-ve_mat2
    INFL2<-infl2.mat
  }

  jjep<-0

  if("reg"%in%method_ep){
    jjep<-jjep+1
    ## error-prone reg -main
    Rx<-cbind(x,x^2)
    Rx<-as.matrix(Rx)
    glm.out<-glm(A~Rx,family="binomial")
    pshatx<-glm.out$fitted.values
    pswx<-mean(y*A/pshatx)-mean(y*(1-A)/(1-pshatx))
    psix<-y*A/pshatx-y*(1-A)/(1-pshatx)
    lm.y<- y
    lm.x<- cbind(x,x^2)
    lm.x<- as.matrix(lm.x)
    dimx<- dim(lm.x)[2]
    lm.out1<-lm(lm.y[which(A==1)]~lm.x[which(A==1),])
    mu1x<-cbind(1,lm.x)%*%lm.out1$coefficients
    lm.out0<-lm(lm.y[which(A==0)]~lm.x[which(A==0),])
    mu0x<-cbind(1,lm.x)%*%lm.out0$coefficients
    B<-apply(cbind(1,lm.x),2,mean)
    B1mat<-cbind(A,lm.x*A)
    B1<-matrix(0,1+dimx,1+dimx)
    for(jj in 1:(1+dimx)){
      for(kk in jj:(1+dimx)){
        B1[jj,kk]<-B1[kk,jj]<-mean(B1mat[,jj]*B1mat[,kk])
      }
    }
    B0mat<-cbind(1-A,lm.x*(1-A))
    B0<-matrix(0,1+dimx,1+dimx)
    for(jj in 1:(1+dimx)){
      for(kk in jj:(1+dimx)){
        B0[jj,kk]<-B0[kk,jj]<-mean(B0mat[,jj]*B0mat[,kk])
      }
    }
    b1<-B%*%ginv(B1)
    b0<-B%*%ginv(B0)
    regx.resid1<-regx.resid0<-rep(0,n)
    regx.resid1[A==1]<-(lm.out1$resid)
    regx.resid0[A==0]<-(lm.out0$resid)
    regx.adj<-(apply( B1mat*matrix(b1,n,1+dimx,byrow=TRUE),1,sum))*regx.resid1*A-
      (apply( B0mat*matrix(b0,n,1+dimx,byrow=TRUE),1,sum))*regx.resid0*(1-A)
    H<-c(mean( A*y*(1-pshatx)/pshatx+(1-A)*y*(pshatx)/(1-pshatx) ),
         apply( (A*y*(1-pshatx)/pshatx+(1-A)*y*(pshatx)/(1-pshatx))*lm.x,2,mean ))
    E<-matrix(0,1+dimx,1+dimx)
    forE<-matrix(1,1+dimx,n);forE[2:(1+dimx),]<-t(lm.x);
    for(jj in 1:(1+dimx)){
      for(kk in jj:(1+dimx)){
        E[jj,kk]<-E[kk,jj]<-mean(pshatx*(1-pshatx)*forE[jj,]*forE[kk,])
      }
    }
    alphax<- -(A-pshatx) * H%*%ginv(E)%*%forE
    infl1.reg<- ( psix+alphax- pswx )/n*n
    mean(infl1.reg)

    ## error prone reg val
    Rxu<-cbind(xx2,xx2^2)
    glm.out<-glm(AA2~Rxu,family="binomial")
    pshatx2<-glm.out$fitted.values
    pswx2<-mean(yy2*AA2/pshatx2)-mean(yy2*(1-AA2)/(1-pshatx2))
    pswx2
    psix2<-yy2*AA2/pshatx2-yy2*(1-AA2)/(1-pshatx2)
    lm.y<- yy2
    lm.x<- cbind(xx2,xx2^2)
    lm.x<-as.matrix(lm.x)
    dimx2<-dim(lm.x)[2]
    lm.out1<-lm(lm.y[which(AA2==1)]~lm.x[which(AA2==1),])
    mu1x2<-cbind(1,lm.x)%*%lm.out1$coefficients
    lm.out0<-lm(lm.y[which(AA2==0)]~lm.x[which(AA2==0),])
    mu0x2<-cbind(1,lm.x)%*%lm.out0$coefficients
    B<-apply(cbind(1,lm.x),2,mean)
    B1mat<-cbind(AA2,lm.x*AA2)
    B1<-matrix(0,1+dimx2,1+dimx2)
    for(jj in 1:(1+dimx2)){
      for(kk in jj:(1+dimx2)){
        B1[jj,kk]<-B1[kk,jj]<-mean(B1mat[,jj]*B1mat[,kk])
      }
    }
    B0mat<-cbind(1-AA2,lm.x*(1-AA2))
    B0<-matrix(0,1+dimx2,1+dimx2)
    for(jj in 1:(1+dimx2)){
      for(kk in jj:(1+dimx2)){
        B0[jj,kk]<-B0[kk,jj]<-mean(B0mat[,jj]*B0mat[,kk])
      }
    }
    b1<-B%*%ginv(B1)
    b0<-B%*%ginv(B0)
    regx2.resid1<-regx2.resid0<-rep(0,n2)
    regx2.resid1[AA2==1]<-(lm.out1$resid)
    regx2.resid0[AA2==0]<-(lm.out0$resid)
    regx2.adj<-apply(B1mat*matrix(b1,n2,1+dimx2,byrow=TRUE),1,sum)*regx2.resid1*AA2-
      apply(B0mat*matrix(b0,n2,1+dimx2,byrow=TRUE),1,sum)*regx2.resid0*(1-AA2)
    H<-c(mean( AA2*yy2*(1-pshatx2)/pshatx2+(1-AA2)*yy2*(pshatx2)/(1-pshatx2) ),
         apply( (AA2*yy2*(1-pshatx2)/pshatx2+(1-AA2)*yy2*(pshatx2)/(1-pshatx2))*lm.x,2,mean ))
    E<-matrix(0,(1+dimx2),(1+dimx2))
    forE<-matrix(1,(1+dimx2),n2);forE[2:(1+dimx2),]<-t(lm.x);
    for(jj in 1:(1+dimx2)){
      for(kk in jj:(1+dimx2)){
        E[jj,kk]<-E[kk,jj]<-mean(pshatx2*(1-pshatx2)*forE[jj,]*forE[kk,])
      }
    }
    alphax2<- -(AA2-pshatx2) * H%*%ginv(E)%*%forE
    infl2x.reg<-( psix2+alphax2- pswx2 )
    infl1.reg[loc2]<- infl1.reg[loc2]-( psix2+alphax2- pswx2 )/n2*n
    EST.ep[jjep]<-mean( mu1x-mu0x )-mean( mu1x2-mu0x2)
    INFL.ep[jjep,]<-infl1.reg

  }


  if("ipw"%in%method_ep){
    jjep<-jjep+1
    ## error-prone psw -main
    Rx<-cbind(x,x^2)
    Rx<-as.matrix(Rx)
    dimx2<-dim(Rx)[2]
    glm.out<-glm(A~Rx,family="binomial")
    pshatx<-glm.out$fitted.values
    pswx<-mean(y*A/pshatx)-mean(y*(1-A)/(1-pshatx))
    psix<-y*A/pshatx-y*(1-A)/(1-pshatx)
    H<-c(mean( A*y*(1-pshatx)/pshatx+(1-A)*y*(pshatx)/(1-pshatx) ),
         apply( (A*y*(1-pshatx)/pshatx+(1-A)*y*(pshatx)/(1-pshatx))*Rx,2,mean ))
    E<-matrix(0,1+dimx2,1+dimx2)
    forE<-matrix(1,(1+dimx2),n);forE[2:(1+dimx2),]<-t(Rx);
    for(jj in 1:(1+dimx2)){
      for(kk in jj:(1+dimx2)){
        E[jj,kk]<-E[kk,jj]<-mean(pshatx*(1-pshatx)*forE[jj,]*forE[kk,])
      }
    }
    alphax<- -(A-pshatx) * H%*%ginv(E)%*%forE
    mean(alphax)
    infl1.ipw<- ( psix+alphax- pswx )/n*n

    ## error-prone psw -val
    Rxu<-cbind(xx2,xx2^2)
    Rxu<-as.matrix(Rxu)
    dimx2<-dim(Rxu)[2]
    glm.out<-glm(AA2~Rxu,family="binomial")
    pshatx2<-glm.out$fitted.values
    pswx2<-mean(yy2*AA2/pshatx2)-mean(yy2*(1-AA2)/(1-pshatx2))
    pswx2
    psix2<-yy2*AA2/pshatx2-yy2*(1-AA2)/(1-pshatx2)
    H<-c(mean( AA2*yy2*(1-pshatx2)/pshatx2+(1-AA2)*yy2*(pshatx2)/(1-pshatx2) ),
         apply( (AA2*yy2*(1-pshatx2)/pshatx2+(1-AA2)*yy2*(pshatx2)/(1-pshatx2))*Rxu,2,mean ))
    E<-matrix(0,1+ dimx2,1+ dimx2)
    forE<-matrix(1,(1+dimx2),n2);forE[2:(1+dimx2),]<-t(Rxu);
    for(jj in 1:(1+dimx2)){
      for(kk in jj:(1+dimx2)){
        E[jj,kk]<-E[kk,jj]<-mean(pshatx2*(1-pshatx2)*forE[jj,]*forE[kk,])
      }
    }
    alphax2<- -(AA2-pshatx2) * H%*%ginv(E)%*%forE
    infl2x.ipw<-( psix2+alphax2- pswx2 )
    infl1.ipw[loc2]<- infl1.ipw[loc2]-( psix2+alphax2- pswx2 )/n2*n

    EST.ep[jjep]<-pswx-pswx2
    INFL.ep[jjep,]<-infl1.ipw

  }


  if("aipw"%in%method_ep){
    jjep<-jjep+1

    ## error-prone reg -main
    Rx<-cbind(x,x^2)
    Rx<-as.matrix(Rx)
    glm.out<-glm(A~Rx,family="binomial")
    pshatx<-glm.out$fitted.values
    pswx<-mean(y*A/pshatx)-mean(y*(1-A)/(1-pshatx))
    psix<-y*A/pshatx-y*(1-A)/(1-pshatx)
    lm.y<- y
    lm.x<- cbind(x,x^2)
    lm.x<- as.matrix(lm.x)
    dimx<- dim(lm.x)[2]
    lm.out1<-lm(lm.y[which(A==1)]~lm.x[which(A==1),])
    mu1x<-cbind(1,lm.x)%*%lm.out1$coefficients
    lm.out0<-lm(lm.y[which(A==0)]~lm.x[which(A==0),])
    mu0x<-cbind(1,lm.x)%*%lm.out0$coefficients

    ## error prone reg val
    Rxu<-cbind(xx2,xx2^2)
    glm.out<-glm(AA2~Rxu,family="binomial")
    pshatx2<-glm.out$fitted.values
    pswx2<-mean(yy2*AA2/pshatx2)-mean(yy2*(1-AA2)/(1-pshatx2))
    psix2<-yy2*AA2/pshatx2-yy2*(1-AA2)/(1-pshatx2)
    lm.y<- yy2
    lm.x<- cbind(xx2,xx2^2)
    lm.x<-as.matrix(lm.x)
    dimx2<-dim(lm.x)[2]
    lm.out1<-lm(lm.y[which(AA2==1)]~lm.x[which(AA2==1),])
    mu1x2<-cbind(1,lm.x)%*%lm.out1$coefficients
    lm.out0<-lm(lm.y[which(AA2==0)]~lm.x[which(AA2==0),])
    mu0x2<-cbind(1,lm.x)%*%lm.out0$coefficients

    ## error-prone psw -main
    Rx<-cbind(x,x^2)
    Rx<-as.matrix(Rx)
    dimx2<-dim(Rx)[2]
    glm.out<-glm(A~Rx,family="binomial")
    pshatx<-glm.out$fitted.values
    pswx<-mean(y*A/pshatx)-mean(y*(1-A)/(1-pshatx))
    psix<-y*A/pshatx-y*(1-A)/(1-pshatx)

    ## error-prone psw -val
    Rxu<-cbind(xx2,xx2^2)
    Rxu<-as.matrix(Rxu)
    dimx2<-dim(Rxu)[2]
    glm.out<-glm(AA2~Rxu,family="binomial")
    pshatx2<-glm.out$fitted.values
    pswx2<-mean(yy2*AA2/pshatx2)-mean(yy2*(1-AA2)/(1-pshatx2))
    psix2<-yy2*AA2/pshatx2-yy2*(1-AA2)/(1-pshatx2)
    infl1.aug<-temp1<-y*A/pshatx-y*(1-A)/(1-pshatx)+( mu1x*(1-A/pshatx)-mu0x*(1-(1-A)/(1-pshatx)) )
    infl1.aug<-infl1.aug-mean(infl1.aug)
    adj.aug<-yy2*AA2/pshatx2-yy2*(1-AA2)/(1-pshatx2)+( mu1x2*(1-AA2/pshatx2)-mu0x2*(1-(1-AA2)/(1-pshatx2)) )
    infl1.aug[loc2]<- infl1.aug[loc2]-( adj.aug -mean(adj.aug))/n2*n

    EST.ep[jjep]<- s1.aug<-mean(temp1)-mean(adj.aug)
    INFL.ep[jjep,]<-infl1.aug
  }

  if("matching"%in%method_ep){
    jjep<-jjep+1

    Rx<-cbind(x,x^2)
    Rx<-as.matrix(Rx)
    glm.out<-glm(A~Rx,family="binomial")
    pshatx<-glm.out$fitted.values
    pswx<-mean(y*A/pshatx)-mean(y*(1-A)/(1-pshatx))
    psix<-y*A/pshatx-y*(1-A)/(1-pshatx)
    lm.y<- y
    lm.x<- cbind(x,x^2)
    lm.x<- as.matrix(lm.x)
    dimx<- dim(lm.x)[2]
    lm.out1<-lm(lm.y[which(A==1)]~lm.x[which(A==1),])
    mu1x<-cbind(1,lm.x)%*%lm.out1$coefficients
    lm.out0<-lm(lm.y[which(A==0)]~lm.x[which(A==0),])
    mu0x<-cbind(1,lm.x)%*%lm.out0$coefficients
    Rxu<-cbind(xx2,xx2^2)
    glm.out<-glm(AA2~Rxu,family="binomial")
    pshatx2<-glm.out$fitted.values
    pswx2<-mean(yy2*AA2/pshatx2)-mean(yy2*(1-AA2)/(1-pshatx2))
    psix2<-yy2*AA2/pshatx2-yy2*(1-AA2)/(1-pshatx2)
    lm.y<- yy2
    lm.x<- cbind(xx2,xx2^2)
    lm.x<-as.matrix(lm.x)
    dimx2<-dim(lm.x)[2]
    lm.out1<-lm(lm.y[which(AA2==1)]~lm.x[which(AA2==1),])
    mu1x2<-cbind(1,lm.x)%*%lm.out1$coefficients
    lm.out0<-lm(lm.y[which(AA2==0)]~lm.x[which(AA2==0),])
    mu0x2<-cbind(1,lm.x)%*%lm.out0$coefficients

    ## error-prone mat -val
    thisps<-xx2
    loc.o<-which(AA2==1)
    Di.nni<-rep(0,length=n2)
    out1<-Matching::Match(Y=yy2,Tr=1-AA2,X=thisps,distance.tolerance=0,ties=FALSE, Weight=2,M=1)
    mdata1<-out1$mdata
    fromto<-loc.o
    Di.nni[fromto]<-table(factor(out1$index.control,levels=fromto))
    Di.nni[loc.o]<-Di.nni[loc.o]+1
    est.nni<-sum(yy2*Di.nni,na.rm=TRUE)/n2
    Ey<-mu1x2
    bias.nni<-sum(Ey[out1$index.control]-Ey[out1$index.treated])/n2
    est1.nni<-est.nni-bias.nni
    truemboot1<-(yy2-Ey)*Di.nni+Ey
    loc.o<-which(AA2==0)
    Di.nni<-rep(0,length=n2)
    out1<-Matching::Match(Y=yy2,Tr=AA2,X=thisps,distance.tolerance=0,ties=FALSE,Weight=2,M=1)
    mdata1<-out1$mdata
    fromto<-loc.o
    Di.nni[fromto]<-table(factor(out1$index.control,levels=fromto))
    Di.nni[loc.o]<-Di.nni[loc.o]+1
    est.nni<-sum(yy2*Di.nni,na.rm=TRUE)/n2
    Ey<-mu1x2
    bias.nni<-sum(Ey[out1$index.control]-Ey[out1$index.treated])/n2
    est0.nni<-est.nni-bias.nni
    truemboot0<-(yy2-Ey)*Di.nni+(Ey)
    truemat<-est1.nni-est0.nni
    truemboot<-truemboot1-truemboot0
    matx2<-truemat
    infl1.mat<-rep(0,n);infl1.mat[loc2]<- -( truemboot1-truemboot0-mean(truemboot1-truemboot0 ) )/n2*n

    infl1.mat2forboot<-infl1.mat
    infl1.mat2forboot[loc2]<- ( truemboot1-truemboot0-mean(truemboot1-truemboot0 ) )/n2*n

    ## error-prone mat -main
    thisps<-x
    loc.o<-which(A==1)
    Di.nni<-rep(0,length=n)
    out1<-Matching::Match(Y=y,Tr=1-A,X=thisps,distance.tolerance=0,ties=FALSE, Weight=2,M=1)
    mdata1<-out1$mdata
    fromto<-loc.o
    Di.nni[fromto]<-table(factor(out1$index.control,levels=fromto))
    Di.nni[loc.o]<-Di.nni[loc.o]+1
    est.nni<-sum(y*Di.nni,na.rm=TRUE)/n
    Ey<-mu1x
    bias.nni<-sum(Ey[out1$index.control]-Ey[out1$index.treated])/n
    est1.nni<-est.nni-bias.nni
    est1.nni
    truemboot1<-(y-Ey)*Di.nni+Ey
    loc.o<-which(A==0)
    Di.nni<-rep(0,length=n)
    out1<-Matching::Match(Y=y,Tr=A,X=thisps,distance.tolerance=0,ties=FALSE,Weight=2,M=1)
    mdata1<-out1$mdata
    fromto<-loc.o
    Di.nni[fromto]<-table(factor(out1$index.control,levels=fromto))
    Di.nni[loc.o]<-Di.nni[loc.o]+1
    est.nni<-sum(y*Di.nni,na.rm=TRUE)/n
    Ey<-mu0x
    bias.nni<-sum(Ey[out1$index.control]-Ey[out1$index.treated])/n
    est0.nni<-est.nni-bias.nni
    truemboot0<-(y-Ey)*Di.nni+(Ey)
    truemat<-est1.nni-est0.nni
    truemboot<-truemboot1-truemboot0
    matx<-truemat
    infl1.mat<- infl1.mat+( truemboot1-truemboot0 -mean(truemboot1-truemboot0) )/n*n
    infl1.mat1forboot<-( truemboot1-truemboot0 -mean(truemboot1-truemboot0) )

    EST.ep[jjep]<- matx-matx2
    INFL.ep[jjep,]<-infl1.mat

  }

  if(kep==0){
    EST<-EST2
    VE<-VE2
  }else{

    for(jj in 1:kep){
      for(kk in jj:kep){
        V.ep[jj,kk]<-V.ep[kk,jj]<-mean(INFL.ep[jj,]*INFL.ep[kk,])/n
      }
    }
    for(jj in 1:kep){
      COV.ep[jj]<-mean(INFL.ep[jj,]*INFL2)/n
    }
    EST<-EST2-COV.ep%*%ginv(V.ep)%*%EST.ep
    VE <-VE2-COV.ep%*%ginv(V.ep)%*%COV.ep
  }

  if(nboot==0){
    VE_boot<-NA
  }

  if(nboot>0){

    if(is.vector(INFL.ep)){
      bootdata<-cbind(INFL2,(INFL.ep))
      kkep<-1
    }else{
      bootdata<-cbind(INFL2,t(INFL.ep))
      kkep<-dim(INFL.ep)[1]
    }

    bootVE2    <-rep(0,1)
    bootCOV.ep <-rep(0,kkep)
    bootV.ep   <-matrix(0,kkep,kkep)
    for(jjboot in 1:nboot){
      ## bootstrap option 1
      bootid<-sample(1:n,n,replace = TRUE)
      ## bootstrap option 2
      bootid<-c(sample(which(I==0),n-n2,replace = TRUE),sample( which(I==1) ,n2,replace = TRUE))
      d<-bootdata[bootid,]
      boot.out<-Wildbootfctn(d)
      bootVE2<-bootVE2+boot.out$VE2
      bootCOV.ep<-bootCOV.ep+boot.out$COV.ep
      bootV.ep<-bootV.ep+boot.out$V.ep
    }
    bootVE2<-bootVE2/nboot
    bootCOV.ep<-bootCOV.ep/nboot
    bootV.ep<-bootV.ep/nboot
    VE_boot<- bootVE2-bootCOV.ep%*% ginv(bootV.ep) %*%bootCOV.ep
  }
  return(list(est=EST,ve=VE,ve_boot=VE_boot))
}



