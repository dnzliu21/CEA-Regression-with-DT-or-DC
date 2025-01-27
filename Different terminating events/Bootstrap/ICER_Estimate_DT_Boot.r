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
  
  
  
  ######## ICER point estimate: 
  ICER_SW_MandQ_Z0=(est_SW_M[3]/1000)/(est_SW_Q[3])
  ICER_SW_MandQ_Z1=(est_SW_M[3]/1000+est_SW_M[4]/1000)/(est_SW_Q[3]+est_SW_Q[4])
  ICER_SW_MandXF_Z0=(est_SW_M[3]/1000)/(est_SW_X_F[3])
  ICER_SW_MandXF_Z1=(est_SW_M[3]/1000+est_SW_M[4]/1000)/(est_SW_X_F[3]+est_SW_X_F[4])
  
  ICER_IMP_MandQ_Z0=(est_IMP_M[3]/1000)/(est_IMP_Q[3])
  ICER_IMP_MandQ_Z1=(est_IMP_M[3]/1000+est_IMP_M[4]/1000)/(est_IMP_Q[3]+est_IMP_Q[4])
  ICER_SWIMP_MandXF_Z0=(est_IMP_M[3]/1000)/(est_SW_X_F[3])
  ICER_SWIMP_MandXF_Z1=(est_IMP_M[3]/1000+est_IMP_M[4]/1000)/(est_SW_X_F[3]+est_SW_X_F[4])
  
  
  
  return(list(est_SW_M=est_SW_M, est_SW_Q=est_SW_Q, est_SW_X_F=est_SW_X_F, 
              est_IMP_M=est_IMP_M, est_IMP_Q=est_IMP_Q,
              
              ICER_SW_MandQ_Z0=ICER_SW_MandQ_Z0, ICER_SW_MandQ_Z1=ICER_SW_MandQ_Z1, 
              ICER_IMP_MandQ_Z0=ICER_IMP_MandQ_Z0, ICER_IMP_MandQ_Z1=ICER_IMP_MandQ_Z1,
              ICER_SW_MandXF_Z0=ICER_SW_MandXF_Z0, ICER_SW_MandXF_Z1=ICER_SW_MandXF_Z1,
              ICER_SWIMP_MandXF_Z0=ICER_SWIMP_MandXF_Z0, ICER_SWIMP_MandXF_Z1=ICER_SWIMP_MandXF_Z1
              ))
  
}




