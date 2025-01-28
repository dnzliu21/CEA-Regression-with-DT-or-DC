####################################################################
## Estimate regression coefficient 						## 
####################################################################

###### estimate K-M at each value of follow-up time ######
KM<-function(X,delta)
  # X=follow-up time, delta=indicator of event
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
  
  for (i in 1:n) {

    # calculating cost history
    if (delta[i]==0) {
      for (j in 1:n) {
        if (X[j]>=X[i])  {
          if (X[i]<=floor(X[j]))  {
            (Accum[i,j]=diag[j]  # diag cost happened at time 0
             +(apply(as.matrix(obs[j,0:floor(X[i])]),1,sum)) # accumulated cost from 0 to floor(X[i])
             +obs[j,ceiling(X[i])]*(X[i]-floor(X[i]))
             +term_true[j]*((max(0,X[i]-max(0,T_true[j]-1)))/(T_true[j]-max(0,T_true[j]-1))) )  # proportional term cost happened at last year
          }
          else {
            (Accum[i,j]=diag[j]  # diag cost happened at time 0
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
ICER <- function(X0,delta0,M_total,M_diag,M_obs,M_term_true,  # cost
                 X,delta,Q_total,Q_obs,T_true,Q_total0,Q_obs0,  # QAL and survival time
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
     K_0=KM(X0,1-delta0)}	#estimate overall K 
  if(Sep.K){	#estimate K separately within each treatment
    K=numeric(n)
    K[A==0]=KM(X[A==0], 1-delta[A==0])
    K[A==1]=KM(X[A==1], 1-delta[A==1])
    K_0=numeric(n)
    K_0[A==0]=KM(X0[A==0], 1-delta0[A==0])
    K_0[A==1]=KM(X0[A==1], 1-delta0[A==1])}
  
  # estimate S by K-M
  if(!Sep.K) 
    {S=KM(X,delta)
     S_0=KM(X0,delta0)}	#estimate overall S
  if(Sep.K){	#estimate S separately within each treatment
    S=numeric(n)
    S[A==0]=KM(X[A==0], delta[A==0])
    S[A==1]=KM(X[A==1], delta[A==1])
    S_0=numeric(n)
    S_0[A==0]=KM(X0[A==0], delta0[A==0])
    S_0[A==1]=KM(X0[A==1], delta0[A==1])
  }
  
  w0=delta/K				#weight for event K(T_i)
  w1=(1-delta)/K    #weight for censored K(C_i)
  w0_0=delta0/K_0
  w1_0=(1-delta0)/K_0
  
  ######## SW estimator for cost
  est_SW_M=solve(t(covar)%*%(covar*w0_0))%*%(apply(covar*w0_0*M_total,2,sum)) 
  ######## SW estimator for QAL
  est_SW_Q=solve(t(covar)%*%(covar*w0))%*%(apply(covar*w0*Q_total,2,sum)) 
  ######## SW estimator for survival time
  est_SW_X=solve(t(covar)%*%(covar*w0))%*%(apply(covar*w0*X,2,sum)) 
  ###############################################################

  
  ######## the accumulated cost and QAL
  M_accum<-CalAccum(n,X0,delta0,M_total,M_diag,M_obs,M_term_true,T_true)
  Q_accum<-CalAccum(n,X,delta,Q_total,numeric(n),Q_obs,numeric(n),T_true)
  ###############################################################
  
  
  ######## IMP estimator for cost and QAL
  cterm.M.1=matrix(0,p,p)
  cterm.M.2=matrix(0,p,1)
  cterm.Q.1=matrix(0,p,p)
  cterm.Q.2=matrix(0,p,1)
  for (i in 1:n) {
    if (delta0[i]==0) {
      Y_M_C0_i=sum(X0>=X0[i])
      cterm.M.1=cterm.M.1+(w1_0[i]/Y_M_C0_i)*(t(covar)%*%(covar*(X0>=X0[i]))) 
      cterm.M.2=cterm.M.2+(w1_0[i])*(apply(covar*(X0>=X0[i])*M_accum[i,],2,sum)/Y_M_C0_i) 
    }
    
    if (delta[i]==0) {
      Y_Q_C_i=sum(X>=X[i])
      cterm.Q.1=cterm.Q.1+(w1[i]/Y_Q_C_i)*(t(covar)%*%(covar*(X>=X[i]))) 
      cterm.Q.2=cterm.Q.2+(w1[i])*(apply(covar*(X>=X[i])*Q_accum[i,],2,sum)/Y_Q_C_i) 
    }
  }
  
  C_M_1=t(covar)%*%(covar*w0_0)+t(covar)%*%(covar*w1_0)-cterm.M.1
  C_M_2=apply(covar*w0_0*M_total,2,sum)+apply(covar*w1_0*M_total,2,sum)-cterm.M.2
  
  est_IMP_M=solve(C_M_1)%*%(C_M_2)
  
  C_Q_1=t(covar)%*%(covar*w0)+t(covar)%*%(covar*w1)-cterm.Q.1
  C_Q_2=apply(covar*w0*Q_total,2,sum)+apply(covar*w1*Q_total,2,sum)-cterm.Q.2
  
  est_IMP_Q=solve(C_Q_1)%*%(C_Q_2)
  ###############################################################
  
  
  
  ######## ICER point estimate: 
  ICER_SW_MandQ_Z0=(est_SW_M[3]/1000)/(est_SW_Q[3])
  ICER_SW_MandQ_Z1=(est_SW_M[3]/1000+est_SW_M[4]/1000)/(est_SW_Q[3]+est_SW_Q[4])
  ICER_SW_MandX_Z0=(est_SW_M[3]/1000)/(est_SW_X[3])
  ICER_SW_MandX_Z1=(est_SW_M[3]/1000+est_SW_M[4]/1000)/(est_SW_X[3]+est_SW_X[4])
  
  ICER_IMP_MandQ_Z0=(est_IMP_M[3]/1000)/(est_IMP_Q[3])
  ICER_IMP_MandQ_Z1=(est_IMP_M[3]/1000+est_IMP_M[4]/1000)/(est_IMP_Q[3]+est_IMP_Q[4])
  ICER_SWIMP_MandX_Z0=(est_IMP_M[3]/1000)/(est_SW_X[3])
  ICER_SWIMP_MandX_Z1=(est_IMP_M[3]/1000+est_IMP_M[4]/1000)/(est_SW_X[3]+est_SW_X[4])

  
  
  return(list(est_SW_M=est_SW_M, est_SW_Q=est_SW_Q, est_SW_X=est_SW_X, 
              est_IMP_M=est_IMP_M, est_IMP_Q=est_IMP_Q,
              
              ICER_SW_MandQ_Z0=ICER_SW_MandQ_Z0, ICER_SW_MandQ_Z1=ICER_SW_MandQ_Z1, 
              ICER_IMP_MandQ_Z0=ICER_IMP_MandQ_Z0, ICER_IMP_MandQ_Z1=ICER_IMP_MandQ_Z1,
              ICER_SW_MandX_Z0=ICER_SW_MandX_Z0, ICER_SW_MandX_Z1=ICER_SW_MandX_Z1,
              ICER_SWIMP_MandX_Z0=ICER_SWIMP_MandX_Z0, ICER_SWIMP_MandX_Z1=ICER_SWIMP_MandX_Z1
              ))
  
}






