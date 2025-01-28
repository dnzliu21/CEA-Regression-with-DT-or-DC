sim.CostEff.Obs <- function(n, Coef.Trt, Coef.Surv.MainEff, Coef.Surv.InterEff, Cost.diag, Cost.fix, Cost.rand, Cost.term, Par.Censor, Par.Censor.0, max.time, k)  {
#n: sample size
#Coef.Trt: coefficients for treatment assignment model using logistic regression
#Coef.Surv.MainEff: coefficients of main effects for generating survival time 
#Coef.Surv.InterEff: coefficients of interaction effects for generating survival time
#Par.Censor: parameter for max of C=censoring time (following uniform distribution [1,Par.Censor]) (Light censoring=25, heavy censoring=14)
#Par.Censor.0: parameter for max of C_0=censoring time for early-censored cost
#max.time: time limit horizon 
#k interval to group outcomes (ie, costs within Day(or Year) 1~k, costs within Day(or Year) (k+1)~2k, etc.)
  
  
  #generate independent binary covariate
  Z=matrix(rbinom(n,1,prob=0.5),n,1)
  
  # generate treatment assignment Trt from logistic model
  PS.lin = cbind(1,Z) %*% Coef.Trt 
  Trt.Prob = exp(PS.lin)/(1+exp(PS.lin))  		
  Trt = rbinom(n,1,prob=Trt.Prob)		#Trt=1 or 0
  
  # generate survival time T, also get potential uncensored outcomes (counterfactuals) used to evaluate methods
  T = rexp(n, 1/exp(cbind(1,Z)%*%Coef.Surv.MainEff+Trt*(cbind(1,Z)%*%Coef.Surv.InterEff)))  # survival time T for treatment group
  
  # generate censoring time C
  C=runif(n)*Par.Censor  #independent censored time for survival time 
  C0=runif(n)*Par.Censor.0 #independent censored time for cost
  C0[C0>C]=C[C0>C] #If C0 > C then we set C0 = C
  
  L=max.time
  X=delta=X0=delta0=TL=numeric(n)
  
  for(i in 1:n)
  {
    TL[i]=min(T[i],L)		#L-truncated survival
    X[i]=min(TL[i],C[i])	#L-truncated follow-up survival time
    delta[i]=(TL[i]<=C[i])	#L-truncated death indicator
    X0[i]=min(TL[i],C0[i])	#L-truncated follow-up cost accum
    delta0[i]=(TL[i]<=C0[i])	#L-truncated cost indicator
  }
  delta[which(max(X)==X)]=1 # force the last observation to be dead to prevent possibly error due to inversed probability of 1/0 
  delta0[which(max(X0)==X0)]=1
  #this error can be prevented by reducing L in model fitting
  

  #### generate grouped costs, each within a k-time interval, also generate potential costs #####
  t.int=(1:10)/10*L                        #10 time points partition time 0-L
  
  # parameters for fixed costs in a time interval depending on treatment and Z
  par_fix=rep(Cost.fix[1],n)
  par_fix[(Z==1)&(Trt==0)]=Cost.fix[2]
  par_fix[(Z==0)&(Trt==1)]=Cost.fix[3]	# If in treatment group, has lower fixed cost
  par_fix[(Z==1)&(Trt==1)]=Cost.fix[4]
  # fixed annual cost
  M_ann_fixed = matrix(rep(rlnorm(n,par_fix,0.2),length(t.int)),n,length(t.int))  
  
  # random annual cost
  M_ann_rand = matrix(rlnorm(length(t.int)*n,Cost.rand,0.2),n,length(t.int))
  
  # M_ann[i,j] is sum of fixed and random cost during time t[j-1] to t[j] for patient i          
  M_ann=M_ann_fixed+M_ann_rand
  
  ## diagnosis cost at time 0
  t.int=c(0,t.int)	#add time 0
  par_diag=rep(Cost.diag[1],n)	#diagnosis cost at time 0
  par_diag[(Trt==1)]=Cost.diag[3]	# If in treatment group, has more diagnosis cost
  M_ann=cbind(rlnorm(n,par_diag,0.2),M_ann)	# add diagnosis cost at time 0
  
  # terminal cost at time T
  M_term=rlnorm(n,Cost.term,0.4)	
  
  # calculate observed and true costs for each year
  # M_obs[i,j]=observed cost of ith people accumulated in jth interval
  M_obs=M_true=M_ann
  M_term_true=M_term

  
  for(i in 1:n)
  {
    for(j in 2:length(t.int))
    {
      # truncate cost by follow-up time X to obtain observed cost until X
      if(X0[i]<=t.int[j]) 
      {
        if (X0[i]<=t.int[j-1]) M_obs[i,j]=0 
        else 
        {
          M_obs[i,j]=M_ann[i,j]*(X0[i]-t.int[j-1])/(t.int[j]-t.int[j-1])	# prorate cost if X is in middle of an interval 
        }
      }
      
      # truncate cost by truncated survival time TL to obtain true cost until TL
      if(TL[i]<=t.int[j])  
      {
        if (TL[i]<=t.int[j-1]) M_true[i,j]=0 
        else 
        {
          M_true[i,j]=M_ann[i,j]*(TL[i]-t.int[j-1])/(t.int[j]-t.int[j-1])	# prorate cost if TL is in middle of an interval 
        }
      }
    }  
  }
  
  # calculate the true total cost
  M_total=numeric(n)
  M_total=apply(M_obs,1,sum)
  for (i in 1:n) {
    if (delta0[i]==1 & X0[i]!=10) {
      M_total[i]=M_total[i]+M_term_true[i]
    }
    if (delta0[i]==1 & X0[i]==10 & (T[i]<max.time+1)) {
      M_total[i]=M_total[i]+M_term_true[i]*(max.time-(T[i]-1))
    }
    if (delta0[i]==0 & (T[i]<C0[i]+1)) {
      M_total[i]=M_total[i]+M_term_true[i]*((C0[i]-max(0,T[i]-1))/(T[i]-max(0,T[i]-1)) )  
    }
  }
  

  
  #### generate grouped quality-adjusted life time, each within a k interval #####
  # generate heart failure time T.HF, using similar coef for survival  
  q_par_fix=rep(0.6,n)
  q_par_fix[(Trt==1)]=0.9 	# If in treatment group, has higher QAL
  
  q_fix=matrix(rep(q_par_fix,length(t.int)-1),n,length(t.int)-1)	#fixed quality for a person
  q_ran=matrix(runif((length(t.int)-1)*n)*0.1,n,length(t.int)-1)   #random quality for a person
  
  
  QALY_obs=QALY_obs0=QALY_true=matrix(0,n,(length(t.int)-1))
  for(i in 1:n)
  {
    for(j in 2:length(t.int))
    {
      # true uncensored grouped quality-adjusted life time 
      if(TL[i]>=t.int[j])
      {
        QALY_true[i,j-1]=(t.int[j]-t.int[j-1])*(q_fix[i,j-1]+q_ran[i,j-1])
       }
      if(TL[i]<t.int[j])
      {
        if (TL[i]<t.int[j-1]) QALY_true[i,j-1]=0 
        else 
        {
          QALY_true[i,j-1]=(TL[i]-t.int[j-1])*(q_fix[i,j-1]+q_ran[i,j-1])
       }
      }
      
   
      # observed grouped quality-adjusted life time until follow-up time X
      # should be truncated at X, the QAL after X doesn't matter anymore
      if(X[i]>=t.int[j])
      {
        QALY_obs[i,j-1]=(t.int[j]-t.int[j-1])*(q_fix[i,j-1]+q_ran[i,j-1])
       }
      if(X[i]<t.int[j])
      {
        if (X[i]<t.int[j-1]) QALY_obs[i,j-1]=0 
        else 
        {
          QALY_obs[i,j-1]=(X[i]-t.int[j-1])*(q_fix[i,j-1]+q_ran[i,j-1])
         }
      }
      
      
      # observed grouped quality-adjusted life time until follow-up time X0/Censor survival time early
      if(X0[i]>=t.int[j])
      {
        QALY_obs0[i,j-1]=(t.int[j]-t.int[j-1])*(q_fix[i,j-1]+q_ran[i,j-1])
      }
      if(X0[i]<t.int[j])
      {
        if (X0[i]<t.int[j-1]) QALY_obs0[i,j-1]=0 
        else 
        {
          QALY_obs0[i,j-1]=(X0[i]-t.int[j-1])*(q_fix[i,j-1]+q_ran[i,j-1])
        }
      }
      
      
    }  
  }  
  
  
  return(list(Patition.times=t.int, Follow.up=X, Follow.up.0=X0, Covariate=Z, Death=delta, Death.0=delta0,
              Trt=Trt, Cost.grp.obs=M_obs, M_total=M_total, QASurv.grp.obs=QALY_obs,
              true_T=T, TL=TL, M_true=M_true, M_term_true=M_term_true, true_C=C, true_C0=C0,
              QALY_true=QALY_true, Trt.Prob=Trt.Prob, QASurv.grp.obs0=QALY_obs0		
              ))	# potential outcomes, used for evaluation 

}

  
  
