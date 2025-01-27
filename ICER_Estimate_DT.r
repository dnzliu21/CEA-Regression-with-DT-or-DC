####################################################################
## Estimate regression coefficient 						## 
####################################################################

###### estimate K-M at each value of follow-up time ######
KM<-function(X,delta)
  # X=follow-up time, delta=indicator of event
  # X=X.F, delta=delta.F
{
  n=length(X)
  K=numeric(n)
  k<-1
  ord<-order(X)
  u=1
  while(u<=n)
  {
    s<-u+1      
    while((s<=n)&(X[ord[u]]==X[ord[s]])) s=s+1  # X[ord[u]] and X[ord[s]] are neighbors in ordered X
    k<-k*(1-sum(delta[ord[u:(s-1)]]==1)/(n-u+1)) 
    K[ord[u:(s-1)]]<-k
    u=s  # to next time point in X
  }
  return(K)
}



###### calculate accumulated cost/QAL for subject j at C_i ###### 
CalAccum<-function(n,X,delta,total,diag,obs,term_true,T_true)
{
  Accum=matrix(0,n,n)
  N.largerX=numeric(n)
  
  for (i in 1:n) {
    N.largerX[i]=sum(X>=X[i])   # Y(C_i)
    # calculating cost history
    if (delta[i]==0) {
      for (j in 1:n) {
        if (X[j]>=X[i])  {
          if (X[i]<=floor(X[j]))  {
            (Accum[i,j]=diag[j]  # diag cost at time 0
             +(apply(as.matrix(obs[j,0:floor(X[i])]),1,sum)) # accumulated cost from 0 to floor(X[i])
             +obs[j,ceiling(X[i])]*(X[i]-floor(X[i]))
             +term_true[j]*((max(0,X[i]-max(0,T_true[j]-1)))/(T_true[j]-max(0,T_true[j]-1))) )  # proportional term cost happened at last year
          }
          else {
            (Accum[i,j]=diag[j]  # diag cost at time 0
             +(apply(as.matrix(obs[j,0:floor(X[i])]),1,sum)) # accumulated cost from 0 to floor(X[i])
             +obs[j,ceiling(X[i])]*((X[i]-floor(X[i]))/(X[j]-floor(X[j])))
             +term_true[j]*((max(0,X[i]-max(0,T_true[j]-1)))/(T_true[j]-max(0,T_true[j]-1)))  )# proportional term cost happened at last year
          }
        }
      }
    }
  }
  return(Accum)
}




###### estimate simple weighted covariance matrix and se ###### 
CalCovComp_SW<-function(X,delta,Dih,w0,S,K,n,p)
{
  tmp1=matrix(0,p,p)
  for (i in 1:n) {
    # attention: sometime matrix multiplication sometime dot multiplication
    if (delta[i]==0) {
      tmp1.1=t(Dih)%*%(Dih*w0*(X>=X[i]))/(n*S[i])
      tmp1.2=apply(w0*Dih*(X>=X[i]),2,sum)/(n*S[i])
      tmp1=tmp1+((1-delta[i])/(K[i]^2))*(tmp1.1-(tmp1.2)%*%(t(tmp1.2)))
    }
  }
  
  return(list(tmp1=tmp1))
}

###### estimate IMP covariance matrix and se ###### 
CalCovComp_IMP<-function(Accum,est.imp,covar,X,delta,Dih,w0,S,K,n,p)
{
  tmp1=tmp2=tmp3=tmp4=matrix(0,p,p)
  N.largerX=numeric(n)
  for (i in 1:n) {
    N.largerX[i]=sum(X>=X[i])
    # attention: sometime matrix multiplication sometime dot multiplication
    if (delta[i]==0) {
      
      tmp1.1=t(Dih)%*%(Dih*w0*(X>=X[i]))/(n*S[i])
      tmp1.2=apply(w0*Dih*(X>=X[i]),2,sum)/(n*S[i])
      tmp1=tmp1+((1-delta[i])/(K[i]^2))*(tmp1.1-tmp1.2%*%t(tmp1.2))
      
      # calculate D_j^h(C_i)
      DjhCi=matrix(0,n,p)
      DjhCi=covar*as.vector(Accum[i,]-covar%*%est.imp)
      
      tmp2.1=t(DjhCi)%*%(DjhCi*(X>=X[i]))/N.largerX[i]
      tmp2.2=apply(DjhCi*(X>=X[i]),2,sum)/N.largerX[i]
      tmp2=tmp2+((1-delta[i])/(K[i]^2))*(tmp2.1-tmp2.2%*%t(tmp2.2))
      
      tmp3.1=t(Dih)%*%(DjhCi*w0*(X>=X[i]))/(n*S[i])
      tmp3=tmp3+((1-delta[i])/(K[i]^2))*(tmp3.1-tmp1.2%*%t(tmp2.2))
      
      tmp4.1=t(DjhCi)%*%(Dih*w0*(X>=X[i]))/(n*S[i])
      tmp4=tmp4+((1-delta[i])/(K[i]^2))*(tmp4.1-tmp2.2%*%t(tmp1.2))
    }
  }
  
  return(list(tmp1=tmp1, tmp2=tmp2, tmp3=tmp3, tmp4=tmp4))
}


###### estimate beta_M and beta_Q covariance matrix and se ###### 
CalCovComp_ICER_SW<-function(Dih_SW_M,Dih_SW_Q,covar,n,p,w0,delta_F,S_F,K_F,X_F,X)
{
  tmp1=matrix(0,p,p)
  N.largerX=numeric(n)
  for (i in 1:n) {
    N.largerX[i]=sum(X_F>=X[i])
    # attention: sometime matrix multiplication sometime dot multiplication
    if (delta_F[i]==0) {
       
      tmp1.1=t(Dih_SW_M)%*%(Dih_SW_Q*w0*(X_F>=X[i]))/(n*S_F[i])
      tmp1.2=apply(w0*Dih_SW_M*(X_F>=X[i]),2,sum)/(n*S_F[i])
      tmp1.3=apply(w0*Dih_SW_Q*(X_F>=X[i]),2,sum)/(n*S_F[i])
      tmp1=tmp1+((1-delta_F[i])/(K_F[i]^2))*(tmp1.1-tmp1.2%*%t(tmp1.3))
    }
  }
  return(list(tmp1=tmp1))
}



CalCovComp_ICER_IMP<-function(M_Accum,Q_Accum,Dih_IMP_M,Dih_IMP_Q,est_IMP_M,est_IMP_Q,covar,n,p,w0,delta_F,S_F,K_F,X_F,X)
{
  tmp1=tmp2=tmp3=tmp4=matrix(0,p,p)
  N.largerX=numeric(n)
  for (i in 1:n) {
    N.largerX[i]=sum(X_F>=X[i])
    # attention: sometime matrix multiplication sometime dot multiplication
    if (delta_F[i]==0) {
      
      tmp1.1=t(Dih_IMP_M)%*%(Dih_IMP_Q*w0*(X_F>=X[i]))/(n*S_F[i])
      tmp1.2=apply(w0*Dih_IMP_M*(X_F>=X[i]),2,sum)/(n*S_F[i])
      tmp1.3=apply(w0*Dih_IMP_Q*(X_F>=X[i]),2,sum)/(n*S_F[i])
      tmp1=tmp1+((1-delta_F[i])/(K_F[i]^2))*(tmp1.1-tmp1.2%*%t(tmp1.3))
      
      # calculate D_j^h(C_i)
      DjhCi_M=DjhCi_Q=matrix(0,n,p)
      DjhCi_M=covar*as.vector(M_Accum[i,]-covar%*%est_IMP_M)
      DjhCi_Q=covar*as.vector(Q_Accum[i,]-covar%*%est_IMP_Q)
      
      tmp2.1=t(DjhCi_M)%*%(DjhCi_Q*(X_F>=X[i]))/N.largerX[i]
      tmp2.2=apply(DjhCi_M*(X_F>=X[i]),2,sum)/N.largerX[i]
      tmp2.3=apply(DjhCi_Q*(X_F>=X[i]),2,sum)/N.largerX[i]
      tmp2=tmp2+((1-delta_F[i])/(K_F[i]^2))*(tmp2.1-tmp2.2%*%t(tmp2.3))
      
      tmp3.1=t(Dih_IMP_M)%*%(DjhCi_Q*w0*(X_F>=X[i]))/(n*S_F[i])
      tmp3=tmp3+((1-delta_F[i])/(K_F[i]^2))*(tmp3.1-tmp1.2%*%t(tmp2.3))
      
      tmp4.1=t(DjhCi_M)%*%(Dih_IMP_Q*w0*(X_F>=X[i]))/(n*S_F[i])
      tmp4=tmp4+((1-delta_F[i])/(K_F[i]^2))*(tmp4.1-tmp2.2%*%t(tmp1.3))
    }
  }
  
  return(list(tmp1=tmp1, tmp2=tmp2, tmp3=tmp3, tmp4=tmp4))
}




CalCovComp_ICER_SWIMP <-function(Dih_SW_X_F,M_accum,Dih_IMP_M,est_IMP_M,covar,n,p,w0,delta_F,S_F,K_F,X_F,X)
{
  tmp1=tmp2=matrix(0,p,p)
  N.largerX=numeric(n)
  for (i in 1:n) {
    N.largerX[i]=sum(X_F>=X[i])
    # attention: sometime matrix multiplication sometime dot multiplication
    if (delta_F[i]==0) {
      
      tmp1.1=t(Dih_IMP_M)%*%(Dih_SW_X_F*w0*(X_F>=X[i]))/(n*S_F[i])
      tmp1.2=apply(w0*Dih_IMP_M*(X_F>=X[i]),2,sum)/(n*S_F[i])
      tmp1.3=apply(w0*Dih_SW_X_F*(X_F>=X[i]),2,sum)/(n*S_F[i])
      tmp1=tmp1+((1-delta_F[i])/(K_F[i]^2))*(tmp1.1-tmp1.2%*%t(tmp1.3))
      
      
      # calculate D_j^h(C_i)
      DjhCi_M=matrix(0,n,p)
      DjhCi_M=covar*as.vector(M_accum[i,]-covar%*%est_IMP_M)
      
      tmp2.1=t(DjhCi_M)%*%(Dih_SW_X_F*w0*(X_F>=X[i]))/(n*S_F[i])
      tmp2.2=apply(DjhCi_M*(X_F>=X[i]),2,sum)/N.largerX[i]
      tmp2=tmp2+((1-delta_F[i])/(K_F[i]^2))*(tmp2.1-tmp2.2%*%t(tmp1.3))
      
    }
  }
  
  return(list(tmp1=tmp1, tmp2=tmp2))
}


###### estimate ICER confidence interval using Fieller's theorem###### 
CalCI_ICER_Fieller<-function(Cov_mat_M,Cov_mat_Q,Cov_mat_MandQ,est_M,est_Q,alpha,Z_group)
{
  z_alpha=qnorm(1-alpha/2)
  
  if (Z_group==0) {
    x=est_M[3]/1000
    y=est_Q[3]
    s_xx=Cov_mat_M[3,3]/(1000^2)
    s_yy=Cov_mat_Q[3,3]
    s_xy=Cov_mat_MandQ[3,3]/1000
  }
  else {
    x=(est_M[3]+est_M[4])/1000
    y=est_Q[3]+est_Q[4]
    s_xx=(Cov_mat_M[3,3]+Cov_mat_M[4,4]+2*Cov_mat_M[3,4])/(1000^2)
    s_yy=Cov_mat_Q[3,3]+Cov_mat_Q[4,4]+2*Cov_mat_Q[3,4]
    s_xy=(Cov_mat_MandQ[3,3]+Cov_mat_MandQ[3,4]+Cov_mat_MandQ[4,3]+Cov_mat_MandQ[4,4])/1000
  }
  
  revert=0
  if ((y^2-(z_alpha^2)*s_yy)<0) {
    revert=1
  }
  
  inbnd=(x*y-(z_alpha^2)*s_xy)^2-(x^2-(z_alpha^2)*s_xx)*(y^2-(z_alpha^2)*s_yy)
  if (inbnd>=0) {
    bnd=sqrt(inbnd)
    upperbnd=((x*y-(z_alpha^2)*s_xy)+bnd)/(y^2-(z_alpha^2)*s_yy)
    lowerbnd=((x*y-(z_alpha^2)*s_xy)-bnd)/(y^2-(z_alpha^2)*s_yy)
  }
  else if (inbnd<0) {
    upperbnd=Inf
    lowerbnd=-Inf
    revert=1
  }
  
  return(list(upperbnd=upperbnd,lowerbnd=lowerbnd,revert=revert))
}



###### estimate ICER confidence interval using delta method ###### 
CalCI_ICER_Delta<-function(Cov_mat_M,Cov_mat_Q,Cov_mat_MandQ,est_M,est_Q,alpha,Z_group)
{
  z_alpha=qnorm(1-alpha/2)
  
  if (Z_group==0) {
    x=est_M[3]/1000
    y=est_Q[3]
    s_xx=Cov_mat_M[3,3]/(1000^2)
    s_yy=Cov_mat_Q[3,3]
    s_xy=Cov_mat_MandQ[3,3]/1000
  }
  else {
    x=(est_M[3]+est_M[4])/1000
    y=est_Q[3]+est_Q[4]
    s_xx=(Cov_mat_M[3,3]+Cov_mat_M[4,4]+2*Cov_mat_M[3,4])/(1000^2)
    s_yy=Cov_mat_Q[3,3]+Cov_mat_Q[4,4]+2*Cov_mat_Q[3,4]
    s_xy=(Cov_mat_MandQ[3,3]+Cov_mat_MandQ[3,4]+Cov_mat_MandQ[4,3]+Cov_mat_MandQ[4,4])/1000
  }
  
  gamma=x/y
  sigma=s_xx/(y^2)-2*x*s_xy/(y^3)+(x^2)*s_yy/(y^4)
  
  
  upperbnd=gamma+z_alpha*sqrt(sigma)
  lowerbnd=gamma-z_alpha*sqrt(sigma)
  
  
  return(list(upperbnd=upperbnd,lowerbnd=lowerbnd))
}




#######  calculate ICER estimators and its confidence interval
#######  calculate ICER by calculating SW and IMP estimators of coefficient beta
#######  calculate s.e of betas and covariance of betas
ICER <- function(X,delta,M_total,M_diag,M_obs,M_term_true,  # cost
                 X_F,delta_F,Q_total,Q_obs,  # QAL and survival time
                 T_true,censor_true,
                 Trt,Z,max.time,type="interaction",covariate=TRUE,interaction=TRUE,Sep.K=FALSE)
  
{
  n=length(X)
  A=Trt
  # set regression model
  ZN=1-Z
  covar=cbind(1,Z,A*ZN,A*Z)  # no interaction
  if(type=="interaction"){covar=cbind(1,Z,A,A*Z)} # with interaction
  if(type=="test"){covar=cbind(1,ZN,A,A*ZN)}
  if(interaction==FALSE){covar=cbind(1,Z,A)}
  if(covariate==FALSE){covar=cbind(1,A)}
  
  p=dim(covar)[2]
  
  # estimate K by K-M
  if(!Sep.K) 
    {K=KM(X,1-delta)
     K_F=KM(X_F,1-delta_F)}	#estimate overall K 
  if(Sep.K){	#estimate K separately within each treatment
    K=numeric(n)
    K[A==0]=KM(X[A==0], 1-delta[A==0])
    K[A==1]=KM(X[A==1], 1-delta[A==1])
    K_F=numeric(n)
    K_F[A==0]=KM(X_F[A==0], 1-delta_F[A==0])
    K_F[A==1]=KM(X_F[A==1], 1-delta_F[A==1])}
  
  # estimate S by K-M
  if(!Sep.K) 
    {S=KM(X,delta)
     S_F=KM(X_F,delta_F)}	#estimate overall S
  if(Sep.K){	#estimate S separately within each treatment
    S=numeric(n)
    S[A==0]=KM(X[A==0], delta[A==0])
    S[A==1]=KM(X[A==1], delta[A==1])
    S_F=numeric(n)
    S_F[A==0]=KM(X_F[A==0], delta_F[A==0])
    S_F[A==1]=KM(X_F[A==1], delta_F[A==1])
  }
  
  w0=delta/K				#weight for event K(T_i)
  w1=(1-delta)/K    #weight for censored K(C_i)
  w0_F=delta_F/K_F
  w1_F=(1-delta_F)/K_F
  
  ######## SW estimator for cost
  est_SW_M=solve(t(covar)%*%(covar*w0))%*%(apply(covar*w0*M_total,2,sum)) 
  ######## SW estimator for QAL
  est_SW_Q=solve(t(covar)%*%(covar*w0_F))%*%(apply(covar*w0_F*Q_total,2,sum)) 
  ######## SW estimator for survival time
  est_SW_X_F=solve(t(covar)%*%(covar*w0_F))%*%(apply(covar*w0_F*X_F,2,sum)) 
  ###############################################################
  
  ######## SW: calculate D_i^h(beta) for M, Q_F, and X_F
  Dih_SW_M=sweep(covar,1,M_total-covar%*%est_SW_M,'*')
  Dih_SW_Q=sweep(covar,1,Q_total-covar%*%est_SW_Q,'*')
  Dih_SW_X_F=sweep(covar,1,X_F-covar%*%est_SW_X_F,'*')
  
  A_mat=t(covar)%*%covar  # \hat I_0*n: pxp matrix
  
  # SW: components for calculating covariance
  Comp_SW_M<-CalCovComp_SW(X,delta,Dih_SW_M,w0,S,K,n,p)
  Comp_SW_Q<-CalCovComp_SW(X_F,delta_F,Dih_SW_Q,w0_F,S_F,K_F,n,p)
  Comp_SW_X_F<-CalCovComp_SW(X_F,delta_F,Dih_SW_X_F,w0_F,S_F,K_F,n,p)
  
  ######## Obtain the covariance of SW estimator
  S_mat_SW_M=(t(Dih_SW_M)%*%(Dih_SW_M*w0))+Comp_SW_M$tmp1
  Cov_SW_M=solve(A_mat)%*%S_mat_SW_M%*%solve(A_mat)
  se_SW_M=sqrt(diag(Cov_SW_M))
  
  S_mat_SW_Q=(t(Dih_SW_Q)%*%(Dih_SW_Q*w0_F))+Comp_SW_Q$tmp1
  Cov_SW_Q=solve(A_mat)%*%S_mat_SW_Q%*%solve(A_mat)
  se_SW_Q=sqrt(diag(Cov_SW_Q))
  
  S_mat_SW_X_F=(t(Dih_SW_X_F)%*%(Dih_SW_X_F*w0_F))+Comp_SW_X_F$tmp1
  Cov_SW_X_F=solve(A_mat)%*%S_mat_SW_X_F%*%solve(A_mat)
  se_SW_X_F=sqrt(diag(Cov_SW_X_F))
  ###############################################################
  
  
  ######## the accumulated cost and QAL
  M_accum<-CalAccum(n,X,delta,M_total,M_diag,M_obs,M_term_true,T_true)
  Q_accum<-CalAccum(n,X_F,delta_F,Q_total,numeric(n),Q_obs,numeric(n),T_true)
  ###############################################################

  
  ######## IMP estimator for cost and QAL
  cterm.M.1=matrix(0,p,p)
  cterm.M.2=matrix(0,p,1)
  cterm.Q.1=matrix(0,p,p)
  cterm.Q.2=matrix(0,p,1)
  for (i in 1:n) {
    if (delta[i]==0) {
      Y_M_C_i=sum(X>=X[i])
      cterm.M.1=cterm.M.1+(w1[i]/Y_M_C_i)*(t(covar)%*%(covar*(X>=X[i]))) 
      cterm.M.2=cterm.M.2+(w1[i])*(apply(covar*(X>=X[i])*M_accum[i,],2,sum)/Y_M_C_i) 
    }
    
    if (delta_F[i]==0) {
      Y_Q_C_i=sum(X_F>=X_F[i])
      cterm.Q.1=cterm.Q.1+(w1_F[i]/Y_Q_C_i)*(t(covar)%*%(covar*(X_F>=X_F[i]))) 
      cterm.Q.2=cterm.Q.2+(w1_F[i])*(apply(covar*(X_F>=X_F[i])*Q_accum[i,],2,sum)/Y_Q_C_i) 
    }
  }
  
  C_M_1=t(covar)%*%(covar*w0)+t(covar)%*%(covar*w1)-cterm.M.1
  C_M_2=apply(covar*w0*M_total,2,sum)+apply(covar*w1*M_total,2,sum)-cterm.M.2
  
  est_IMP_M=solve(C_M_1)%*%(C_M_2)
  
  C_Q_1=t(covar)%*%(covar*w0_F)+t(covar)%*%(covar*w1_F)-cterm.Q.1
  C_Q_2=apply(covar*w0_F*Q_total,2,sum)+apply(covar*w1_F*Q_total,2,sum)-cterm.Q.2
  
  est_IMP_Q=solve(C_Q_1)%*%(C_Q_2)
  ###############################################################
  
  ######## IMP: calculate D_i^h(beta) for M, Q, and X_F
  Dih_IMP_M=sweep(covar,1,M_total-covar%*%est_IMP_M,'*')
  Dih_IMP_Q=sweep(covar,1,Q_total-covar%*%est_IMP_Q,'*')
  
  # IMP: components for calculating covariance
  Comp_IMP_M=CalCovComp_IMP(M_accum,est_IMP_M,covar,X,delta,Dih_IMP_M,w0,S,K,n,p)
  Comp_IMP_Q=CalCovComp_IMP(Q_accum,est_IMP_Q,covar,X_F,delta_F,Dih_IMP_Q,w0_F,S_F,K_F,n,p)
  
  ######## Estimate covariance of IMP estimator
  S_mat_IMP_M=(t(Dih_IMP_M)%*%(Dih_IMP_M*w0))+Comp_IMP_M$tmp1+Comp_IMP_M$tmp2-Comp_IMP_M$tmp3-Comp_IMP_M$tmp4
  Cov_IMP_M=solve(A_mat)%*%S_mat_IMP_M%*%solve(A_mat)
  se_IMP_M=sqrt(diag(Cov_IMP_M))
  
  S_mat_IMP_Q=(t(Dih_IMP_Q)%*%(Dih_IMP_Q*w0_F))+Comp_IMP_Q$tmp1+Comp_IMP_Q$tmp2-Comp_IMP_Q$tmp3-Comp_IMP_Q$tmp4
  Cov_IMP_Q=solve(A_mat)%*%S_mat_IMP_Q%*%solve(A_mat)
  se_IMP_Q=sqrt(diag(Cov_IMP_Q))
  ###############################################################
  
  
  ######## ICER point estimate: 
  ICER_SW_MandQ_Z0=(est_SW_M[3]/1000)/(est_SW_Q[3])
  ICER_SW_MandQ_Z1=(est_SW_M[3]/1000+est_SW_M[4]/1000)/(est_SW_Q[3]+est_SW_Q[4])
  ICER_SW_MandXF_Z0=(est_SW_M[3]/1000)/(est_SW_X_F[3])
  ICER_SW_MandXF_Z1=(est_SW_M[3]/1000+est_SW_M[4]/1000)/(est_SW_X_F[3]+est_SW_X_F[4])
  
  ICER_IMP_MandQ_Z0=(est_IMP_M[3]/1000)/(est_IMP_Q[3])
  ICER_IMP_MandQ_Z1=(est_IMP_M[3]/1000+est_IMP_M[4]/1000)/(est_IMP_Q[3]+est_IMP_Q[4])
  ICER_SWIMP_MandXF_Z0=(est_IMP_M[3]/1000)/(est_SW_X_F[3])
  ICER_SWIMP_MandXF_Z1=(est_IMP_M[3]/1000+est_IMP_M[4]/1000)/(est_SW_X_F[3]+est_SW_X_F[4])
  
  
  
  ######## Covariance of beta's (regression coefficients)
  ## SW (M and Q_F)
  Comp_SW_ICER_MandQ=CalCovComp_ICER_SW(Dih_SW_M,Dih_SW_Q,covar,n,p,w0,delta_F,S_F,K_F,X_F,X)
  S_mat_SW_MandQ=( (t(Dih_SW_M)%*%(Dih_SW_Q*w0))-(apply(Dih_SW_M*w0,2,sum))%*%t((apply(Dih_SW_Q*w0_F,2,sum)))/n + Comp_SW_ICER_MandQ$tmp1)
  Cov_SW_MandQ=solve(A_mat)%*%S_mat_SW_MandQ%*%solve(A_mat)
  
  ## SW (M and X_F)
  Comp_SW_ICER_MandXF=CalCovComp_ICER_SW(Dih_SW_M,Dih_SW_X_F,covar,n,p,w0,delta_F,S_F,K_F,X_F,X)
  S_mat_SW_MandXF=( (t(Dih_SW_M)%*%(Dih_SW_X_F*w0))-(apply(Dih_SW_M*w0,2,sum))%*%t((apply(Dih_SW_X_F*w0_F,2,sum)))/n + Comp_SW_ICER_MandXF$tmp1)
  Cov_SW_MandXF=solve(A_mat)%*%S_mat_SW_MandXF%*%solve(A_mat)
  
  
  
  ## IMP (M and Q_F)
  Comp_IMP_ICER_MandQ=CalCovComp_ICER_IMP(M_accum,Q_accum,Dih_IMP_M,Dih_IMP_Q,est_IMP_M,est_IMP_Q,covar,n,p,w0,delta_F,S_F,K_F,X_F,X)
  S_mat_IMP_MandQ=( (t(Dih_IMP_M)%*%(Dih_IMP_Q*w0))-(apply(Dih_IMP_M*w0,2,sum))%*%t((apply(Dih_IMP_Q*w0_F,2,sum)))/n
                      +Comp_IMP_ICER_MandQ$tmp1+Comp_IMP_ICER_MandQ$tmp2-Comp_IMP_ICER_MandQ$tmp3-Comp_IMP_ICER_MandQ$tmp4 )
  Cov_IMP_MandQ=solve(A_mat)%*%S_mat_IMP_MandQ%*%solve(A_mat)
 
  ## IMP_M and SW_X_F
  Comp_SWIMP_ICER_MandXF=CalCovComp_ICER_SWIMP(Dih_SW_X_F,M_accum,Dih_IMP_M,est_IMP_M,covar,n,p,w0,delta_F,S_F,K_F,X_F,X)
  S_mat_SWIMP_MandXF=( (t(Dih_IMP_M)%*%(Dih_SW_X_F*w0))-(apply(Dih_IMP_M*w0,2,sum))%*%t((apply(Dih_SW_X_F*w0_F,2,sum)))/n
                    +Comp_SWIMP_ICER_MandXF$tmp1-Comp_SWIMP_ICER_MandXF$tmp2 )
  Cov_SWIMP_MandXF=solve(A_mat)%*%S_mat_SWIMP_MandXF%*%solve(A_mat)
  
  
  
  # 95%
  ######## CI by Fieller for ICER
  ICER_CI_SW_MandQ_Z095=CalCI_ICER_Fieller(Cov_SW_M,Cov_SW_Q,Cov_SW_MandQ,est_SW_M,est_SW_Q,0.05,Z_group=0)
  ICER_CI_SW_MandQ_Z195=CalCI_ICER_Fieller(Cov_SW_M,Cov_SW_Q,Cov_SW_MandQ,est_SW_M,est_SW_Q,0.05,Z_group=1)
  
  ICER_CI_SW_MandXF_Z095=CalCI_ICER_Fieller(Cov_SW_M,Cov_SW_X_F,Cov_SW_MandXF,est_SW_M,est_SW_X_F,0.05,Z_group=0)
  ICER_CI_SW_MandXF_Z195=CalCI_ICER_Fieller(Cov_SW_M,Cov_SW_X_F,Cov_SW_MandXF,est_SW_M,est_SW_X_F,0.05,Z_group=1)
  
  ICER_CI_IMP_MandQ_Z095=CalCI_ICER_Fieller(Cov_IMP_M,Cov_IMP_Q,Cov_IMP_MandQ,est_IMP_M,est_IMP_Q,0.05,Z_group=0)
  ICER_CI_IMP_MandQ_Z195=CalCI_ICER_Fieller(Cov_IMP_M,Cov_IMP_Q,Cov_IMP_MandQ,est_IMP_M,est_IMP_Q,0.05,Z_group=1)
  
  ICER_CI_SWIMP_MandXF_Z095=CalCI_ICER_Fieller(Cov_IMP_M,Cov_SW_X_F,Cov_SWIMP_MandXF,est_IMP_M,est_SW_X_F,0.05,Z_group=0)
  ICER_CI_SWIMP_MandXF_Z195=CalCI_ICER_Fieller(Cov_IMP_M,Cov_SW_X_F,Cov_SWIMP_MandXF,est_IMP_M,est_SW_X_F,0.05,Z_group=1)
  
  
  ######## CI by Delta method for ICER
  ICER_CI_SW_MandQ_Z095_Delta=CalCI_ICER_Delta(Cov_SW_M,Cov_SW_Q,Cov_SW_MandQ,est_SW_M,est_SW_Q,0.05,Z_group=0)
  ICER_CI_SW_MandQ_Z195_Delta=CalCI_ICER_Delta(Cov_SW_M,Cov_SW_Q,Cov_SW_MandQ,est_SW_M,est_SW_Q,0.05,Z_group=1)
  
  ICER_CI_SW_MandXF_Z095_Delta=CalCI_ICER_Delta(Cov_SW_M,Cov_SW_X_F,Cov_SW_MandXF,est_SW_M,est_SW_X_F,0.05,Z_group=0)
  ICER_CI_SW_MandXF_Z195_Delta=CalCI_ICER_Delta(Cov_SW_M,Cov_SW_X_F,Cov_SW_MandXF,est_SW_M,est_SW_X_F,0.05,Z_group=1)
  
  ICER_CI_IMP_MandQ_Z095_Delta=CalCI_ICER_Delta(Cov_IMP_M,Cov_IMP_Q,Cov_IMP_MandQ,est_IMP_M,est_IMP_Q,0.05,Z_group=0)
  ICER_CI_IMP_MandQ_Z195_Delta=CalCI_ICER_Delta(Cov_IMP_M,Cov_IMP_Q,Cov_IMP_MandQ,est_IMP_M,est_IMP_Q,0.05,Z_group=1)
  
  ICER_CI_SWIMP_MandXF_Z095_Delta=CalCI_ICER_Delta(Cov_IMP_M,Cov_SW_X_F,Cov_SWIMP_MandXF,est_IMP_M,est_SW_X_F,0.05,Z_group=0)
  ICER_CI_SWIMP_MandXF_Z195_Delta=CalCI_ICER_Delta(Cov_IMP_M,Cov_SW_X_F,Cov_SWIMP_MandXF,est_IMP_M,est_SW_X_F,0.05,Z_group=1)
  
  
  # 90%
  ######## CI by Fieller for ICER
  ICER_CI_SW_MandQ_Z090=CalCI_ICER_Fieller(Cov_SW_M,Cov_SW_Q,Cov_SW_MandQ,est_SW_M,est_SW_Q,0.1,Z_group=0)
  ICER_CI_SW_MandQ_Z190=CalCI_ICER_Fieller(Cov_SW_M,Cov_SW_Q,Cov_SW_MandQ,est_SW_M,est_SW_Q,0.1,Z_group=1)
  
  ICER_CI_SW_MandXF_Z090=CalCI_ICER_Fieller(Cov_SW_M,Cov_SW_X_F,Cov_SW_MandXF,est_SW_M,est_SW_X_F,0.1,Z_group=0)
  ICER_CI_SW_MandXF_Z190=CalCI_ICER_Fieller(Cov_SW_M,Cov_SW_X_F,Cov_SW_MandXF,est_SW_M,est_SW_X_F,0.1,Z_group=1)
  
  ICER_CI_IMP_MandQ_Z090=CalCI_ICER_Fieller(Cov_IMP_M,Cov_IMP_Q,Cov_IMP_MandQ,est_IMP_M,est_IMP_Q,0.1,Z_group=0)
  ICER_CI_IMP_MandQ_Z190=CalCI_ICER_Fieller(Cov_IMP_M,Cov_IMP_Q,Cov_IMP_MandQ,est_IMP_M,est_IMP_Q,0.1,Z_group=1)
  
  ICER_CI_SWIMP_MandXF_Z090=CalCI_ICER_Fieller(Cov_IMP_M,Cov_SW_X_F,Cov_SWIMP_MandXF,est_IMP_M,est_SW_X_F,0.1,Z_group=0)
  ICER_CI_SWIMP_MandXF_Z190=CalCI_ICER_Fieller(Cov_IMP_M,Cov_SW_X_F,Cov_SWIMP_MandXF,est_IMP_M,est_SW_X_F,0.1,Z_group=1)
  
  
  ######## CI by Delta for ICER
  ICER_CI_SW_MandQ_Z090_Delta=CalCI_ICER_Delta(Cov_SW_M,Cov_SW_Q,Cov_SW_MandQ,est_SW_M,est_SW_Q,0.1,Z_group=0)
  ICER_CI_SW_MandQ_Z190_Delta=CalCI_ICER_Delta(Cov_SW_M,Cov_SW_Q,Cov_SW_MandQ,est_SW_M,est_SW_Q,0.1,Z_group=1)
  
  ICER_CI_SW_MandXF_Z090_Delta=CalCI_ICER_Delta(Cov_SW_M,Cov_SW_X_F,Cov_SW_MandXF,est_SW_M,est_SW_X_F,0.1,Z_group=0)
  ICER_CI_SW_MandXF_Z190_Delta=CalCI_ICER_Delta(Cov_SW_M,Cov_SW_X_F,Cov_SW_MandXF,est_SW_M,est_SW_X_F,0.1,Z_group=1)
  
  ICER_CI_IMP_MandQ_Z090_Delta=CalCI_ICER_Delta(Cov_IMP_M,Cov_IMP_Q,Cov_IMP_MandQ,est_IMP_M,est_IMP_Q,0.1,Z_group=0)
  ICER_CI_IMP_MandQ_Z190_Delta=CalCI_ICER_Delta(Cov_IMP_M,Cov_IMP_Q,Cov_IMP_MandQ,est_IMP_M,est_IMP_Q,0.1,Z_group=1)
  
  ICER_CI_SWIMP_MandXF_Z090_Delta=CalCI_ICER_Delta(Cov_IMP_M,Cov_SW_X_F,Cov_SWIMP_MandXF,est_IMP_M,est_SW_X_F,0.1,Z_group=0)
  ICER_CI_SWIMP_MandXF_Z190_Delta=CalCI_ICER_Delta(Cov_IMP_M,Cov_SW_X_F,Cov_SWIMP_MandXF,est_IMP_M,est_SW_X_F,0.1,Z_group=1)
  
  
  # 80%
  ######## CI by Fieller for ICER
  ICER_CI_SW_MandQ_Z080=CalCI_ICER_Fieller(Cov_SW_M,Cov_SW_Q,Cov_SW_MandQ,est_SW_M,est_SW_Q,0.2,Z_group=0)
  ICER_CI_SW_MandQ_Z180=CalCI_ICER_Fieller(Cov_SW_M,Cov_SW_Q,Cov_SW_MandQ,est_SW_M,est_SW_Q,0.2,Z_group=1)
  
  ICER_CI_SW_MandXF_Z080=CalCI_ICER_Fieller(Cov_SW_M,Cov_SW_X_F,Cov_SW_MandXF,est_SW_M,est_SW_X_F,0.2,Z_group=0)
  ICER_CI_SW_MandXF_Z180=CalCI_ICER_Fieller(Cov_SW_M,Cov_SW_X_F,Cov_SW_MandXF,est_SW_M,est_SW_X_F,0.2,Z_group=1)
  
  ICER_CI_IMP_MandQ_Z080=CalCI_ICER_Fieller(Cov_IMP_M,Cov_IMP_Q,Cov_IMP_MandQ,est_IMP_M,est_IMP_Q,0.2,Z_group=0)
  ICER_CI_IMP_MandQ_Z180=CalCI_ICER_Fieller(Cov_IMP_M,Cov_IMP_Q,Cov_IMP_MandQ,est_IMP_M,est_IMP_Q,0.2,Z_group=1)
  
  ICER_CI_SWIMP_MandXF_Z080=CalCI_ICER_Fieller(Cov_IMP_M,Cov_SW_X_F,Cov_SWIMP_MandXF,est_IMP_M,est_SW_X_F,0.2,Z_group=0)
  ICER_CI_SWIMP_MandXF_Z180=CalCI_ICER_Fieller(Cov_IMP_M,Cov_SW_X_F,Cov_SWIMP_MandXF,est_IMP_M,est_SW_X_F,0.2,Z_group=1)
  
  
  ######## CI by Delta for ICER
  ICER_CI_SW_MandQ_Z080_Delta=CalCI_ICER_Delta(Cov_SW_M,Cov_SW_Q,Cov_SW_MandQ,est_SW_M,est_SW_Q,0.2,Z_group=0)
  ICER_CI_SW_MandQ_Z180_Delta=CalCI_ICER_Delta(Cov_SW_M,Cov_SW_Q,Cov_SW_MandQ,est_SW_M,est_SW_Q,0.2,Z_group=1)
  
  ICER_CI_SW_MandXF_Z080_Delta=CalCI_ICER_Delta(Cov_SW_M,Cov_SW_X_F,Cov_SW_MandXF,est_SW_M,est_SW_X_F,0.2,Z_group=0)
  ICER_CI_SW_MandXF_Z180_Delta=CalCI_ICER_Delta(Cov_SW_M,Cov_SW_X_F,Cov_SW_MandXF,est_SW_M,est_SW_X_F,0.2,Z_group=1)
  
  ICER_CI_IMP_MandQ_Z080_Delta=CalCI_ICER_Delta(Cov_IMP_M,Cov_IMP_Q,Cov_IMP_MandQ,est_IMP_M,est_IMP_Q,0.2,Z_group=0)
  ICER_CI_IMP_MandQ_Z180_Delta=CalCI_ICER_Delta(Cov_IMP_M,Cov_IMP_Q,Cov_IMP_MandQ,est_IMP_M,est_IMP_Q,0.2,Z_group=1)
  
  ICER_CI_SWIMP_MandXF_Z080_Delta=CalCI_ICER_Delta(Cov_IMP_M,Cov_SW_X_F,Cov_SWIMP_MandXF,est_IMP_M,est_SW_X_F,0.2,Z_group=0)
  ICER_CI_SWIMP_MandXF_Z180_Delta=CalCI_ICER_Delta(Cov_IMP_M,Cov_SW_X_F,Cov_SWIMP_MandXF,est_IMP_M,est_SW_X_F,0.2,Z_group=1)
  
  
  
  return(list(est_SW_M=est_SW_M, est_SW_Q=est_SW_Q, est_SW_X_F=est_SW_X_F, 
              est_IMP_M=est_IMP_M, est_IMP_Q=est_IMP_Q,
              Cov_SW_M=Cov_SW_M, Cov_SW_Q=Cov_SW_Q, Cov_SW_X_F=Cov_SW_X_F, 
              Cov_IMP_M=Cov_IMP_M, Cov_IMP_Q=Cov_IMP_Q,
              
              ICER_SW_MandQ_Z0=ICER_SW_MandQ_Z0, ICER_SW_MandQ_Z1=ICER_SW_MandQ_Z1, 
              ICER_IMP_MandQ_Z0=ICER_IMP_MandQ_Z0, ICER_IMP_MandQ_Z1=ICER_IMP_MandQ_Z1,
              ICER_SW_MandXF_Z0=ICER_SW_MandXF_Z0, ICER_SW_MandXF_Z1=ICER_SW_MandXF_Z1,
              ICER_SWIMP_MandXF_Z0=ICER_SWIMP_MandXF_Z0, ICER_SWIMP_MandXF_Z1=ICER_SWIMP_MandXF_Z1,
              
              Cov_SW_MandQ=Cov_SW_MandQ, Cov_SW_MandXF=Cov_SW_MandXF, 
              Cov_IMP_MandQ=Cov_IMP_MandQ, Cov_SWIMP_MandXF=Cov_SWIMP_MandXF,
              
              
              # 95%
              ICER_CI_SW_MandQ_Z095_up=ICER_CI_SW_MandQ_Z095$upperbnd, ICER_CI_SW_MandQ_Z095_low=ICER_CI_SW_MandQ_Z095$lowerbnd,
              ICER_CI_SW_MandQ_Z095_revert=ICER_CI_SW_MandQ_Z095$revert, 
              ICER_CI_SW_MandQ_Z195_up=ICER_CI_SW_MandQ_Z195$upperbnd, ICER_CI_SW_MandQ_Z195_low=ICER_CI_SW_MandQ_Z195$lowerbnd,
              ICER_CI_SW_MandQ_Z195_revert=ICER_CI_SW_MandQ_Z195$revert, 
              
              ICER_CI_SW_MandXF_Z095_up=ICER_CI_SW_MandXF_Z095$upperbnd, ICER_CI_SW_MandXF_Z095_low=ICER_CI_SW_MandXF_Z095$lowerbnd,
              ICER_CI_SW_MandXF_Z095_revert=ICER_CI_SW_MandXF_Z095$revert, 
              ICER_CI_SW_MandXF_Z195_up=ICER_CI_SW_MandXF_Z195$upperbnd, ICER_CI_SW_MandXF_Z195_low=ICER_CI_SW_MandXF_Z195$lowerbnd,
              ICER_CI_SW_MandXF_Z195_revert=ICER_CI_SW_MandXF_Z195$revert, 
              
              ICER_CI_IMP_MandQ_Z095_up=ICER_CI_IMP_MandQ_Z095$upperbnd, ICER_CI_IMP_MandQ_Z095_low=ICER_CI_IMP_MandQ_Z095$lowerbnd,
              ICER_CI_IMP_MandQ_Z095_revert=ICER_CI_IMP_MandQ_Z095$revert,
              ICER_CI_IMP_MandQ_Z195_up=ICER_CI_IMP_MandQ_Z195$upperbnd, ICER_CI_IMP_MandQ_Z195_low=ICER_CI_IMP_MandQ_Z195$lowerbnd,
              ICER_CI_IMP_MandQ_Z195_revert=ICER_CI_IMP_MandQ_Z195$revert,
              
              ICER_CI_SWIMP_MandXF_Z095_up=ICER_CI_SWIMP_MandXF_Z095$upperbnd, ICER_CI_SWIMP_MandXF_Z095_low=ICER_CI_SWIMP_MandXF_Z095$lowerbnd,
              ICER_CI_SWIMP_MandXF_Z095_revert=ICER_CI_SWIMP_MandXF_Z095$revert, 
              ICER_CI_SWIMP_MandXF_Z195_up=ICER_CI_SWIMP_MandXF_Z195$upperbnd, ICER_CI_SWIMP_MandXF_Z195_low=ICER_CI_SWIMP_MandXF_Z195$lowerbnd,
              ICER_CI_SWIMP_MandXF_Z195_revert=ICER_CI_SWIMP_MandXF_Z195$revert,
              
              
              ICER_CI_SW_MandQ_Z095_Delta_up=ICER_CI_SW_MandQ_Z095_Delta$upperbnd, ICER_CI_SW_MandQ_Z095_Delta_low=ICER_CI_SW_MandQ_Z095_Delta$lowerbnd,
              ICER_CI_SW_MandQ_Z195_Delta_up=ICER_CI_SW_MandQ_Z195_Delta$upperbnd, ICER_CI_SW_MandQ_Z195_Delta_low=ICER_CI_SW_MandQ_Z195_Delta$lowerbnd,
              
              ICER_CI_SW_MandXF_Z095_Delta_up=ICER_CI_SW_MandXF_Z095_Delta$upperbnd, ICER_CI_SW_MandXF_Z095_Delta_low=ICER_CI_SW_MandXF_Z095_Delta$lowerbnd,
              ICER_CI_SW_MandXF_Z195_Delta_up=ICER_CI_SW_MandXF_Z195_Delta$upperbnd, ICER_CI_SW_MandXF_Z195_Delta_low=ICER_CI_SW_MandXF_Z195_Delta$lowerbnd,
              
              ICER_CI_IMP_MandQ_Z095_Delta_up=ICER_CI_IMP_MandQ_Z095_Delta$upperbnd, ICER_CI_IMP_MandQ_Z095_Delta_low=ICER_CI_IMP_MandQ_Z095_Delta$lowerbnd,
              ICER_CI_IMP_MandQ_Z195_Delta_up=ICER_CI_IMP_MandQ_Z195_Delta$upperbnd, ICER_CI_IMP_MandQ_Z195_Delta_low=ICER_CI_IMP_MandQ_Z195_Delta$lowerbnd,
              
              ICER_CI_SWIMP_MandXF_Z095_Delta_up=ICER_CI_SWIMP_MandXF_Z095_Delta$upperbnd, ICER_CI_SWIMP_MandXF_Z095_Delta_low=ICER_CI_SWIMP_MandXF_Z095_Delta$lowerbnd,
              ICER_CI_SWIMP_MandXF_Z195_Delta_up=ICER_CI_SWIMP_MandXF_Z195_Delta$upperbnd, ICER_CI_SWIMP_MandXF_Z195_Delta_low=ICER_CI_SWIMP_MandXF_Z195_Delta$lowerbnd,
              
              
              
              # 90%
              ICER_CI_SW_MandQ_Z090_up=ICER_CI_SW_MandQ_Z090$upperbnd, ICER_CI_SW_MandQ_Z090_low=ICER_CI_SW_MandQ_Z090$lowerbnd,
              ICER_CI_SW_MandQ_Z090_revert=ICER_CI_SW_MandQ_Z090$revert, 
              ICER_CI_SW_MandQ_Z190_up=ICER_CI_SW_MandQ_Z190$upperbnd, ICER_CI_SW_MandQ_Z190_low=ICER_CI_SW_MandQ_Z190$lowerbnd,
              ICER_CI_SW_MandQ_Z190_revert=ICER_CI_SW_MandQ_Z190$revert, 
              
              ICER_CI_SW_MandXF_Z090_up=ICER_CI_SW_MandXF_Z090$upperbnd, ICER_CI_SW_MandXF_Z090_low=ICER_CI_SW_MandXF_Z090$lowerbnd,
              ICER_CI_SW_MandXF_Z090_revert=ICER_CI_SW_MandXF_Z090$revert, 
              ICER_CI_SW_MandXF_Z190_up=ICER_CI_SW_MandXF_Z190$upperbnd, ICER_CI_SW_MandXF_Z190_low=ICER_CI_SW_MandXF_Z190$lowerbnd,
              ICER_CI_SW_MandXF_Z190_revert=ICER_CI_SW_MandXF_Z190$revert, 
              
              ICER_CI_IMP_MandQ_Z090_up=ICER_CI_IMP_MandQ_Z090$upperbnd, ICER_CI_IMP_MandQ_Z090_low=ICER_CI_IMP_MandQ_Z090$lowerbnd,
              ICER_CI_IMP_MandQ_Z090_revert=ICER_CI_IMP_MandQ_Z090$revert,
              ICER_CI_IMP_MandQ_Z190_up=ICER_CI_IMP_MandQ_Z190$upperbnd, ICER_CI_IMP_MandQ_Z190_low=ICER_CI_IMP_MandQ_Z190$lowerbnd,
              ICER_CI_IMP_MandQ_Z190_revert=ICER_CI_IMP_MandQ_Z190$revert,
              
              ICER_CI_SWIMP_MandXF_Z090_up=ICER_CI_SWIMP_MandXF_Z090$upperbnd, ICER_CI_SWIMP_MandXF_Z090_low=ICER_CI_SWIMP_MandXF_Z090$lowerbnd,
              ICER_CI_SWIMP_MandXF_Z090_revert=ICER_CI_SWIMP_MandXF_Z090$revert, 
              ICER_CI_SWIMP_MandXF_Z190_up=ICER_CI_SWIMP_MandXF_Z190$upperbnd, ICER_CI_SWIMP_MandXF_Z190_low=ICER_CI_SWIMP_MandXF_Z190$lowerbnd,
              ICER_CI_SWIMP_MandXF_Z190_revert=ICER_CI_SWIMP_MandXF_Z190$revert,
              
              
              ICER_CI_SW_MandQ_Z090_Delta_up=ICER_CI_SW_MandQ_Z090_Delta$upperbnd, ICER_CI_SW_MandQ_Z090_Delta_low=ICER_CI_SW_MandQ_Z090_Delta$lowerbnd,
              ICER_CI_SW_MandQ_Z190_Delta_up=ICER_CI_SW_MandQ_Z190_Delta$upperbnd, ICER_CI_SW_MandQ_Z190_Delta_low=ICER_CI_SW_MandQ_Z190_Delta$lowerbnd,
              
              ICER_CI_SW_MandXF_Z090_Delta_up=ICER_CI_SW_MandXF_Z090_Delta$upperbnd, ICER_CI_SW_MandXF_Z090_Delta_low=ICER_CI_SW_MandXF_Z090_Delta$lowerbnd,
              ICER_CI_SW_MandXF_Z190_Delta_up=ICER_CI_SW_MandXF_Z190_Delta$upperbnd, ICER_CI_SW_MandXF_Z190_Delta_low=ICER_CI_SW_MandXF_Z190_Delta$lowerbnd,
              
              ICER_CI_IMP_MandQ_Z090_Delta_up=ICER_CI_IMP_MandQ_Z090_Delta$upperbnd, ICER_CI_IMP_MandQ_Z090_Delta_low=ICER_CI_IMP_MandQ_Z090_Delta$lowerbnd,
              ICER_CI_IMP_MandQ_Z190_Delta_up=ICER_CI_IMP_MandQ_Z190_Delta$upperbnd, ICER_CI_IMP_MandQ_Z190_Delta_low=ICER_CI_IMP_MandQ_Z190_Delta$lowerbnd,
              
              ICER_CI_SWIMP_MandXF_Z090_Delta_up=ICER_CI_SWIMP_MandXF_Z090_Delta$upperbnd, ICER_CI_SWIMP_MandXF_Z090_Delta_low=ICER_CI_SWIMP_MandXF_Z090_Delta$lowerbnd,
              ICER_CI_SWIMP_MandXF_Z190_Delta_up=ICER_CI_SWIMP_MandXF_Z190_Delta$upperbnd, ICER_CI_SWIMP_MandXF_Z190_Delta_low=ICER_CI_SWIMP_MandXF_Z190_Delta$lowerbnd,
              
              
              # 80%
              ICER_CI_SW_MandQ_Z080_up=ICER_CI_SW_MandQ_Z080$upperbnd, ICER_CI_SW_MandQ_Z080_low=ICER_CI_SW_MandQ_Z080$lowerbnd,
              ICER_CI_SW_MandQ_Z080_revert=ICER_CI_SW_MandQ_Z080$revert, 
              ICER_CI_SW_MandQ_Z180_up=ICER_CI_SW_MandQ_Z180$upperbnd, ICER_CI_SW_MandQ_Z180_low=ICER_CI_SW_MandQ_Z180$lowerbnd,
              ICER_CI_SW_MandQ_Z180_revert=ICER_CI_SW_MandQ_Z180$revert, 
              
              ICER_CI_SW_MandXF_Z080_up=ICER_CI_SW_MandXF_Z080$upperbnd, ICER_CI_SW_MandXF_Z080_low=ICER_CI_SW_MandXF_Z080$lowerbnd,
              ICER_CI_SW_MandXF_Z080_revert=ICER_CI_SW_MandXF_Z080$revert, 
              ICER_CI_SW_MandXF_Z180_up=ICER_CI_SW_MandXF_Z180$upperbnd, ICER_CI_SW_MandXF_Z180_low=ICER_CI_SW_MandXF_Z180$lowerbnd,
              ICER_CI_SW_MandXF_Z180_revert=ICER_CI_SW_MandXF_Z180$revert, 
              
              ICER_CI_IMP_MandQ_Z080_up=ICER_CI_IMP_MandQ_Z080$upperbnd, ICER_CI_IMP_MandQ_Z080_low=ICER_CI_IMP_MandQ_Z080$lowerbnd,
              ICER_CI_IMP_MandQ_Z080_revert=ICER_CI_IMP_MandQ_Z080$revert,
              ICER_CI_IMP_MandQ_Z180_up=ICER_CI_IMP_MandQ_Z180$upperbnd, ICER_CI_IMP_MandQ_Z180_low=ICER_CI_IMP_MandQ_Z180$lowerbnd,
              ICER_CI_IMP_MandQ_Z180_revert=ICER_CI_IMP_MandQ_Z180$revert,
              
              ICER_CI_SWIMP_MandXF_Z080_up=ICER_CI_SWIMP_MandXF_Z080$upperbnd, ICER_CI_SWIMP_MandXF_Z080_low=ICER_CI_SWIMP_MandXF_Z080$lowerbnd,
              ICER_CI_SWIMP_MandXF_Z080_revert=ICER_CI_SWIMP_MandXF_Z080$revert, 
              ICER_CI_SWIMP_MandXF_Z180_up=ICER_CI_SWIMP_MandXF_Z180$upperbnd, ICER_CI_SWIMP_MandXF_Z180_low=ICER_CI_SWIMP_MandXF_Z180$lowerbnd,
              ICER_CI_SWIMP_MandXF_Z180_revert=ICER_CI_SWIMP_MandXF_Z180$revert,
              
              
              ICER_CI_SW_MandQ_Z080_Delta_up=ICER_CI_SW_MandQ_Z080_Delta$upperbnd, ICER_CI_SW_MandQ_Z080_Delta_low=ICER_CI_SW_MandQ_Z080_Delta$lowerbnd,
              ICER_CI_SW_MandQ_Z180_Delta_up=ICER_CI_SW_MandQ_Z180_Delta$upperbnd, ICER_CI_SW_MandQ_Z180_Delta_low=ICER_CI_SW_MandQ_Z180_Delta$lowerbnd,
              
              ICER_CI_SW_MandXF_Z080_Delta_up=ICER_CI_SW_MandXF_Z080_Delta$upperbnd, ICER_CI_SW_MandXF_Z080_Delta_low=ICER_CI_SW_MandXF_Z080_Delta$lowerbnd,
              ICER_CI_SW_MandXF_Z180_Delta_up=ICER_CI_SW_MandXF_Z180_Delta$upperbnd, ICER_CI_SW_MandXF_Z180_Delta_low=ICER_CI_SW_MandXF_Z180_Delta$lowerbnd,
              
              ICER_CI_IMP_MandQ_Z080_Delta_up=ICER_CI_IMP_MandQ_Z080_Delta$upperbnd, ICER_CI_IMP_MandQ_Z080_Delta_low=ICER_CI_IMP_MandQ_Z080_Delta$lowerbnd,
              ICER_CI_IMP_MandQ_Z180_Delta_up=ICER_CI_IMP_MandQ_Z180_Delta$upperbnd, ICER_CI_IMP_MandQ_Z180_Delta_low=ICER_CI_IMP_MandQ_Z180_Delta$lowerbnd,
              
              ICER_CI_SWIMP_MandXF_Z080_Delta_up=ICER_CI_SWIMP_MandXF_Z080_Delta$upperbnd, ICER_CI_SWIMP_MandXF_Z080_Delta_low=ICER_CI_SWIMP_MandXF_Z080_Delta$lowerbnd,
              ICER_CI_SWIMP_MandXF_Z180_Delta_up=ICER_CI_SWIMP_MandXF_Z180_Delta$upperbnd, ICER_CI_SWIMP_MandXF_Z180_Delta_low=ICER_CI_SWIMP_MandXF_Z180_Delta$lowerbnd
              
              
              ))
  
}



#######  calculate if true ICER is covered in the confidence interval
# Fieller's method
CalCoverage_Fieller <- function(lowerbnd, upperbnd, revert, truevalue) {
  if (revert==0) {(lowerbnd<=truevalue)&(upperbnd>=truevalue)}
  else {(lowerbnd<=truevalue)|(upperbnd>=truevalue)}
}
# Delta's method
CalCoverage_Delta <- function(lowerbnd, upperbnd, truevalue) {
  {(lowerbnd<=truevalue)&(upperbnd>=truevalue)}
}





