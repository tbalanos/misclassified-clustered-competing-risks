############################################################
############################################################
##############################
##############################   bssmle
##############################
############################################################
############################################################

#Elements of the package intccr

library(alabama)
library(splines)
library(survival)

naive_b<-function(data,x,q){
  t<-x
  nk=floor(length(t)^(1/3))
  max<-nk+1
  knots<-quantile(t,seq(0,1,by=1/(nk+1)))[2:max]
  B<-bs(t,knots=knots,degree=3,intercept=TRUE)
  
  ## Generate beta
  b<-seq(from=0.0001,to=0.875,by=((0.875-0.0001)/(dim(B)[2]-1)))
  b<-log(b^3)
  b<-c(b,b,rep(0,times=2*q)) 
  b
}

bssmle_clustered <- function(data, covariates = c("z1", "z2", "z3")) {
  # extract misclassification probabilities
  p21 <- data$p21
  p12 <- data$p12
  
  x <- data$x
  d <- (data$c_obs > 0)
  
  # Create design matrix dynamically based on variable names
  Z <- as.matrix(data[, covariates, drop = FALSE])
  
  ## B-spline basis matrix 
  t<-x
  nk=floor(length(t)^(1/3))
  max<-nk+1
  knots<-quantile(t,seq(0,1,by=1/(nk+1)))[2:max]
  Bt<-bs(t,knots=knots,degree=3,intercept=TRUE)
  
  ## B-spline basis matrix for predictions
  times<-seq(from=min(t),to=max(t),
             by=((max(t)-min(t))/99))
  Bpred<-predict(Bt,times)
  
  ## Derivative of the B-spline matrix
  B<-bs(t,knots=knots,degree=2,intercept=TRUE)
  nbs=dim(B)[2]
  dBt=B
  augkn<-c(rep(attr(B,"Boundary.knots")[1],times=3)
           ,attr(B,"knots"),rep(attr(B,"Boundary.knots")[2],times=3))
  for(i in 1:nbs){
    dBt[,i]<-(3/(augkn[i+3]-augkn[i]))*B[,i]
  }
  rm(B)
  
  ## Generate auxiliary variables
  n<-dim(Bt)[2]
  q<-dim(Z)[2]
  d1<-(data$c_obs==1)
  d2<-(data$c_obs==2)
  
  
  ## Define the starting values for the optimizer
  b0<-naive_b(data,x,q)
  
  ## Define the function for the negative log-likelihood
  nLL<-function(b) {
    phi1<-b[1:n]
    phi2<-b[(n+1):(2*n)]
    b1<-b[(2*n+1):(2*n+q)]
    b2<-b[(2*n+q+1):(2*n+2*q)]
    dphi1<-phi1
    for(i in 2:n){
      dphi1[i]<-phi1[i]-phi1[i-1]
    }
    dphi1<-dphi1[2:n]
    dphi2<-phi2
    for(i in 2:n){
      dphi2[i]<-phi2[i]-phi2[i-1]
    }
    dphi2<-dphi2[2:n]
    
    #Create cumulative hazards
    BS1<-Bt%*%phi1
    dBS1<-dBt%*%dphi1
    BS2<-Bt%*%phi2
    dBS2<-dBt%*%dphi2
    
    bz_1<-Z%*%b1
    bz_2<-Z%*%b2
    H1<-exp(BS1+bz_1)
    h1<-exp(BS1+bz_1)*dBS1
    H2<-exp(BS2+bz_2)
    h2<-exp(BS2+bz_2)*dBS2
    
    #Calculate the loglikelihood	
    #ill<-d1*log(h1)+d2*log(h2)-(H1+H2)
    h1<-as.vector(h1)
    h2<-as.vector(h2)
    H1<-as.vector(H1)
    H2<-as.vector(H2)
    ill<-d1*log((1-p21)*h1+p12*h2)+
      d2*log(p21*h1+(1-p12)*h2)-(H1+H2)
    
    # ICS correction
    # nll<- -sum(ill)
    nll <- -sum(ill/data$Mi)
    nll
  }
  
  
  
  ## Define the function for the score
  Grad<-function(b) {
    phi1<-b[1:n]
    phi2<-b[(n+1):(2*n)]
    b1<-b[(2*n+1):(2*n+q)]
    b2<-b[(2*n+q+1):(2*n+2*q)]
    dphi1<-phi1
    for(i in 2:n){
      dphi1[i]<-phi1[i]-phi1[i-1]
    }
    dphi1<-dphi1[2:n]
    dphi2<-phi2
    for(i in 2:n){
      dphi2[i]<-phi2[i]-phi2[i-1]
    }
    dphi2<-dphi2[2:n]
    
    #Create cumulative hazards
    BS1<-Bt%*%phi1
    dBS1<-dBt%*%dphi1
    BS2<-Bt%*%phi2
    dBS2<-dBt%*%dphi2
    
    bz_1<-Z%*%b1
    bz_2<-Z%*%b2
    H1<-exp(BS1+bz_1)
    h1<-exp(BS1+bz_1)*dBS1
    H2<-exp(BS2+bz_2)
    h2<-exp(BS2+bz_2)*dBS2
    
    h1<-as.vector(h1)
    h2<-as.vector(h2)
    H1<-as.vector(H1)
    H2<-as.vector(H2)
    #Create derivative matrix for B(t) w.r.t. theta
    der1B1<-Bt
    der1B1<-cbind(der1B1,matrix(rep(0,times=(n*length(data$x))),ncol=n))
    der1B2<-Bt
    der1B2<-cbind(matrix(rep(0,times=(n*length(data$x))),ncol=n),der1B2)
    
    #Create derivative matrix for dB(t)/dt w.r.t. theta
    dB1plus<-cbind(dBt,rep(0,times=dim(data)[1]))
    dB2plus<-cbind(dBt,rep(0,times=dim(data)[1]))
    
    der2B1<- cbind(rep(0,times=length(data$x)),dBt)-dB1plus
    der2B1<-cbind(der2B1,matrix(rep(0,times=(n*length(data$x))),ncol=n))
    der2B2<- cbind(rep(0,times=length(data$x)),dBt)-dB2plus
    der2B2<-cbind(matrix(rep(0,times=(n*length(data$x))),ncol=n),der2B2)
    
    #Gradient for the spline parameters
    deriv1<-as.data.frame((der2B1+as.vector(dBS1)*der1B1))
    deriv2<-as.data.frame((der2B2+as.vector(dBS2)*der1B2))
    
    iG1<-(d1*((1-p21)*H1*deriv1)/((1-p21)*h1+p12*h2)+
            d2*(p21*H1*deriv1)/(p21*h1+(1-p12)*h2)-
            H1*as.data.frame(der1B1))/data$Mi
    iG2<-(d1*(p12*H2*deriv2)/((1-p21)*h1+p12*h2)+
            d2*((1-p12)*H2*deriv2)/(p21*h1+(1-p12)*h2)-
            H2*as.data.frame(der1B2))/data$Mi
    iG<-iG1+iG2
    G<- -colSums(iG)
    
    #Gradient for the covariate effects
    
    G1<--t(Z)%*%((d1*((1-p21)*h1)/((1-p21)*h1+p12*h2)+d2*(p21*h1)/(p21*h1+(1-p12)*h2)-H1)/data$Mi)
    G2<--t(Z)%*%((d1*(p12*h2)/((1-p21)*h1+p12*h2)+d2*((1-p12)*h2)/(p21*h1+(1-p12)*h2)-H2)/data$Mi)
    
    G<- c(G,G1,G2)
    G
  }
  
  ## Monotonicity constraints
  p=2*n
  ui<-matrix(rep(0,times=(p*(p-2))),ncol=p,nrow=(p-2),byrow=TRUE)
  for(i in 1:(n-1)){
    for(j in 1:n){
      ui[i,j]<-(i==j-1)-(i==j)
    }
  }
  for(i in n:(n+n-2)){
    for(j in (n+1):(2*n)){
      ui[i,j]<-(i==j-2)-(i==j-1)
    }
  }
  zeros<-matrix(rep(0,times=dim(ui)[1]*2*q),nrow=dim(ui)[1])
  ui<-cbind(ui,zeros)
  ci<-rep(0.0000001,times=dim(ui)[1])
  
  est<-try(constrOptim(b0, f=nLL, grad=Grad, ui=ui,
                       ci=ci,control=list(maxit=2000)),
           silent=TRUE)
  if(class(est)!="try-error"){
    if(est$convergence==0){
      beta<-est$par
      #print(-est$value)
    } else {
      beta<-NA
    }
  } else {
    beta<-NA
  }
  return(beta)
}