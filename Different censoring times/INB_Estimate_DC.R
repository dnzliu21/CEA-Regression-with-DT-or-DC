GetEstimateResult=readRDS("/home/dnliu21/code/CEA/simulation/DC/Results/Results_ICER_Est_Heavy_N1200_nsim2000_logn_xxxx_xx_xx.RDS")
EstimateResult = data.frame(do.call(what = rbind, args = GetEstimateResult))

n=1200
WTP15=10
nsim=2000
alpha=0.05
z_alpha=qnorm(1-alpha/2)
rundate="xxxx_xx_xx"

#coef for survival model
Coef.Surv.MainEff=c(2.2, -0.5) #coef for survival model, 1st for intercept, others for main effects of covariates
Coef.Surv.InterEff=c(-0.3, 1.2)  #coef for survival model, 1st for main effect of treatment, others for treatment-covariate interactions

# Parameter for max censoring time=10 years
Par.Censor=30	 
Censorstatus="Heavy"
Par.Censor.0=15


#coef for treatment assignment model
Coef.Trt=c(0,0)		# using logistic regression
k=1					#k=1-year interval to group outcomes
max.time=10				

####################   Calculate true values of coefficients, only need to run once ############## 
true.mean.M=matrix(0,2,2) 	#true.mean.M[i,j]=mean cost with ith Z and jth A (Z=0,1; A=0,1) 
true.mean.TL=matrix(0,2,2) 	#true.mean.TL[i,j]=mean survival (trucated by L) with ith Z and jth A (Z=0,1; A=0,1) 
true.mean.Q=matrix(0,2,2)   #true.mean.Q[i,j]=mean QAL (based TL) with ith Z and jth A (Z=0,1; A=0,1) 
Z0=A0=0:1 

Q_fix_q_mean=matrix(c(0.6,0.6,0.9,0.9),2,2)
Q_rand_q_mean=0.1/2

Cost.diag=c(7.5,7.5,8.5,8.5)
Cost.fix=c(6.8,6.5,7.2,6)
Cost.rand=4
Cost.term=7.2

# log-normal distribution for cost
M_diag_mean=matrix(exp(Cost.diag+(0.2^2)/2),2,2)
M_term_mean=exp(Cost.term+(0.4^2)/2)
M_fix_mean=matrix(exp(Cost.fix+(0.2^2)/2),2,2)
M_rand_mean=exp(Cost.rand+(0.2^2)/2)


for(i in 1:2) {
  for(j in 1:2)
  { 
    lambda=1/exp(cbind(1,Z0[i])%*%Coef.Surv.MainEff+(cbind(A0[j],A0[j]*Z0[i])%*%Coef.Surv.InterEff))	#parameter in exp dist of T
    
    true.mean.TL[i,j]=(1-exp(-lambda*max.time))/lambda
    true.mean.M[i,j]= (M_diag_mean[i,j]+M_term_mean*(1+(exp(-lambda*(max.time+1))-exp(-lambda*max.time))/lambda)+
                         true.mean.TL[i,j]*(M_fix_mean[i,j]+M_rand_mean))
    true.mean.Q[i,j]=(Q_fix_q_mean[i,j]+Q_rand_q_mean)*true.mean.TL[i,j]
    
  }
}

true.beta.M=c(true.mean.M[1,1], true.mean.M[2,1]-true.mean.M[1,1], true.mean.M[1,2]-true.mean.M[1,1], 
              true.mean.M[2,2]-true.mean.M[1,2]-true.mean.M[2,1]+true.mean.M[1,1])
true.beta.TL=c(true.mean.TL[1,1], true.mean.TL[2,1]-true.mean.TL[1,1], true.mean.TL[1,2]-true.mean.TL[1,1], 
               true.mean.TL[2,2]-true.mean.TL[1,2]-true.mean.TL[2,1]+true.mean.TL[1,1])
true.beta.Q=c(true.mean.Q[1,1], true.mean.Q[2,1]-true.mean.Q[1,1], true.mean.Q[1,2]-true.mean.Q[1,1], 
              true.mean.Q[2,2]-true.mean.Q[1,2]-true.mean.Q[2,1]+true.mean.Q[1,1])


true.INB.M.Q.Z0=WTP15*true.beta.Q[3]-true.beta.M[3]/1000
true.INB.M.Q.Z1=WTP15*(true.beta.Q[3]+true.beta.Q[4])-(true.beta.M[3]/1000+true.beta.M[4]/1000)
true.INB.M.X.Z0=WTP15*true.beta.TL[3]-true.beta.M[3]/1000
true.INB.M.X.Z1=WTP15*(true.beta.TL[3]+true.beta.TL[4])-(true.beta.M[3]/1000+true.beta.M[4]/1000)



# Estimated covariance
CovSW_INB_MQ_Z0=matrix(unlist(EstimateResult$ICER_SW_MandQ_Cov_Z0),nsim,1,byrow=TRUE)
CovSW_INB_MQ_Z1=matrix(unlist(EstimateResult$ICER_SW_MandQ_Cov_Z1),nsim,1,byrow=TRUE)
CovIMP_INB_MQ_Z0=matrix(unlist(EstimateResult$ICER_IMP_MandQ_Cov_Z0),nsim,1,byrow=TRUE)
CovIMP_INB_MQ_Z1=matrix(unlist(EstimateResult$ICER_IMP_MandQ_Cov_Z1),nsim,1,byrow=TRUE)

CovSW_INB_MX_Z0=matrix(unlist(EstimateResult$ICER_SW_MandX_Cov_Z0),nsim,1,byrow=TRUE)
CovSW_INB_MX_Z1=matrix(unlist(EstimateResult$ICER_SW_MandX_Cov_Z1),nsim,1,byrow=TRUE)
CovSWIMP_INB_MX_Z0=matrix(unlist(EstimateResult$ICER_SWIMP_MandX_Cov_Z0),nsim,1,byrow=TRUE)
CovSWIMP_INB_MX_Z1=matrix(unlist(EstimateResult$ICER_SWIMP_MandX_Cov_Z1),nsim,1,byrow=TRUE)


# Sample cov
estDiff_SW_M_Z0=matrix(unlist(EstimateResult$estDiff_SW_M_Z0),nsim,1,byrow=TRUE)
estDiff_SW_M_Z1=matrix(unlist(EstimateResult$estDiff_SW_M_Z1),nsim,1,byrow=TRUE)
estDiff_SW_Q_Z0=matrix(unlist(EstimateResult$estDiff_SW_Q_Z0),nsim,1,byrow=TRUE)
estDiff_SW_Q_Z1=matrix(unlist(EstimateResult$estDiff_SW_Q_Z1),nsim,1,byrow=TRUE)
estDiff_SW_X_Z0=matrix(unlist(EstimateResult$estDiff_SW_X_Z0),nsim,1,byrow=TRUE)
estDiff_SW_X_Z1=matrix(unlist(EstimateResult$estDiff_SW_X_Z1),nsim,1,byrow=TRUE)
estDiff_IMP_M_Z0=matrix(unlist(EstimateResult$estDiff_IMP_M_Z0),nsim,1,byrow=TRUE)
estDiff_IMP_M_Z1=matrix(unlist(EstimateResult$estDiff_IMP_M_Z1),nsim,1,byrow=TRUE)
estDiff_IMP_Q_Z0=matrix(unlist(EstimateResult$estDiff_IMP_Q_Z0),nsim,1,byrow=TRUE)
estDiff_IMP_Q_Z1=matrix(unlist(EstimateResult$estDiff_IMP_Q_Z1),nsim,1,byrow=TRUE)


# Estimated Variance
Sxx_SW_M_Z0=matrix(unlist(EstimateResult$Sxx_SW_M_Z0),nsim,1,byrow=TRUE)
Syy_SW_Q_Z0=matrix(unlist(EstimateResult$Syy_SW_Q_Z0),nsim,1,byrow=TRUE)
Syy_SW_X_Z0=matrix(unlist(EstimateResult$Syy_SW_X_Z0),nsim,1,byrow=TRUE)
Sxx_SW_M_Z1=matrix(unlist(EstimateResult$Sxx_SW_M_Z1),nsim,1,byrow=TRUE)
Syy_SW_Q_Z1=matrix(unlist(EstimateResult$Syy_SW_Q_Z1),nsim,1,byrow=TRUE)
Syy_SW_X_Z1=matrix(unlist(EstimateResult$Syy_SW_X_Z1),nsim,1,byrow=TRUE)
Sxx_IMP_M_Z0=matrix(unlist(EstimateResult$Sxx_IMP_M_Z0),nsim,1,byrow=TRUE)
Syy_IMP_Q_Z0=matrix(unlist(EstimateResult$Syy_IMP_Q_Z0),nsim,1,byrow=TRUE)
Sxx_IMP_M_Z1=matrix(unlist(EstimateResult$Sxx_IMP_M_Z1),nsim,1,byrow=TRUE)
Syy_IMP_Q_Z1=matrix(unlist(EstimateResult$Syy_IMP_Q_Z1),nsim,1,byrow=TRUE)


GetResult <- function(estDiff_Eff_Z,estDiff_M_Z,true.INB.Z,nsim,WTP,Syy_Eff_Z,Sxx_M_Z,Cov_INB_Z,z_alpha) {
  INB_Z=Var_INB_Z_nsim=cover_INB_Z_nsim=numeric(nsim)
  for (i in 1:nsim) {
    INB_Z[i]=WTP*estDiff_Eff_Z[i]-estDiff_M_Z[i]/1000
    Var_INB_Z_nsim[i]=(WTP^2)*Syy_Eff_Z[i]+Sxx_M_Z[i]/(1000^2)-2*WTP*Cov_INB_Z[i]/1000
    upperbnd=INB_Z[i]+z_alpha*sqrt(Var_INB_Z_nsim[i])
    lowerbnd=INB_Z[i]-z_alpha*sqrt(Var_INB_Z_nsim[i])
    cover_INB_Z_nsim[i]=ifelse((lowerbnd<=true.INB.Z)&(upperbnd>=true.INB.Z),1,0)
  }
  # Mean
  mean_INB_Z=mean(INB_Z,na.rm=TRUE)
  # bias
  bias_INB_Z=mean_INB_Z-true.INB.Z
  # ESE
  Var_INB_Z=mean(Var_INB_Z_nsim,na.rm=TRUE)
  ESE_INB_Z=sqrt(Var_INB_Z)
  # SSE
  SSE_INB_Z=sd(INB_Z,na.rm=TRUE)
  # CP
  CP_INB_Z=sum(cover_INB_Z_nsim)/nsim
  
  return(list(mean_INB_Z=mean_INB_Z,bias_INB_Z=bias_INB_Z, ESE_INB_Z=ESE_INB_Z, SSE_INB_Z=SSE_INB_Z, CP_INB_Z=CP_INB_Z))
}


ResSW_MQ_Z0=GetResult(estDiff_SW_Q_Z0,estDiff_SW_M_Z0,true.INB.M.Q.Z0,nsim,WTP15,Syy_SW_Q_Z0,Sxx_SW_M_Z0,CovSW_INB_MQ_Z0,z_alpha)
ResIMP_MQ_Z0=GetResult(estDiff_IMP_Q_Z0,estDiff_IMP_M_Z0,true.INB.M.Q.Z0,nsim,WTP15,Syy_IMP_Q_Z0,Sxx_IMP_M_Z0,CovIMP_INB_MQ_Z0,z_alpha)

ResSW_MQ_Z1=GetResult(estDiff_SW_Q_Z1,estDiff_SW_M_Z1,true.INB.M.Q.Z1,nsim,WTP15,Syy_SW_Q_Z1,Sxx_SW_M_Z1,CovSW_INB_MQ_Z1,z_alpha)
ResIMP_MQ_Z1=GetResult(estDiff_IMP_Q_Z1,estDiff_IMP_M_Z1,true.INB.M.Q.Z1,nsim,WTP15,Syy_IMP_Q_Z1,Sxx_IMP_M_Z1,CovIMP_INB_MQ_Z1,z_alpha)

ResSW_MX_Z0=GetResult(estDiff_SW_X_Z0,estDiff_SW_M_Z0,true.INB.M.X.Z0,nsim,WTP15,Syy_SW_X_Z0,Sxx_SW_M_Z0,CovSW_INB_MX_Z0,z_alpha)
ResSWIMP_MX_Z0=GetResult(estDiff_SW_X_Z0,estDiff_IMP_M_Z0,true.INB.M.X.Z0,nsim,WTP15,Syy_SW_X_Z0,Sxx_IMP_M_Z0,CovSWIMP_INB_MX_Z0,z_alpha)

ResSW_MX_Z1=GetResult(estDiff_SW_X_Z1,estDiff_SW_M_Z1,true.INB.M.X.Z1,nsim,WTP15,Syy_SW_X_Z1,Sxx_SW_M_Z1,CovSW_INB_MX_Z1,z_alpha)
ResSWIMP_MX_Z1=GetResult(estDiff_SW_X_Z1,estDiff_IMP_M_Z1,true.INB.M.X.Z1,nsim,WTP15,Syy_SW_X_Z1,Sxx_IMP_M_Z1,CovSWIMP_INB_MX_Z1,z_alpha)

results15 = data.frame(
  rep(c("INB"),1,each=4),
  rep(c(10),1,each=4),
  rep(c("SW","IMP"),2,each=1),
  rep(c("Z=0","Z=1"),1,each=2),
  round(t(cbind(ResSW_MQ_Z0$mean_INB_Z, ResIMP_MQ_Z0$mean_INB_Z, ResSW_MQ_Z1$mean_INB_Z, ResIMP_MQ_Z1$mean_INB_Z)),2),
  round(t(cbind(ResSW_MQ_Z0$bias_INB_Z, ResIMP_MQ_Z0$bias_INB_Z, ResSW_MQ_Z1$bias_INB_Z, ResIMP_MQ_Z1$bias_INB_Z)),2),
  round(t(cbind(ResSW_MQ_Z0$SSE_INB_Z, ResIMP_MQ_Z0$SSE_INB_Z, ResSW_MQ_Z1$SSE_INB_Z, ResIMP_MQ_Z1$SSE_INB_Z)),2),
  round(t(cbind(ResSW_MQ_Z0$ESE_INB_Z, ResIMP_MQ_Z0$ESE_INB_Z, ResSW_MQ_Z1$ESE_INB_Z, ResIMP_MQ_Z1$ESE_INB_Z)),2),
  round(t(cbind(ResSW_MQ_Z0$CP_INB_Z, ResIMP_MQ_Z0$CP_INB_Z, ResSW_MQ_Z1$CP_INB_Z, ResIMP_MQ_Z1$CP_INB_Z)),3),
  
  round(t(cbind(ResSW_MX_Z0$mean_INB_Z, ResSWIMP_MX_Z0$mean_INB_Z, ResSW_MX_Z1$mean_INB_Z, ResSWIMP_MX_Z1$mean_INB_Z)),2),
  round(t(cbind(ResSW_MX_Z0$bias_INB_Z, ResSWIMP_MX_Z0$bias_INB_Z, ResSW_MX_Z1$bias_INB_Z, ResSWIMP_MX_Z1$bias_INB_Z)),2),
  round(t(cbind(ResSW_MX_Z0$SSE_INB_Z, ResSWIMP_MX_Z0$SSE_INB_Z, ResSW_MX_Z1$SSE_INB_Z, ResSWIMP_MX_Z1$SSE_INB_Z)),2),
  round(t(cbind(ResSW_MX_Z0$ESE_INB_Z, ResSWIMP_MX_Z0$ESE_INB_Z, ResSW_MX_Z1$ESE_INB_Z, ResSWIMP_MX_Z1$ESE_INB_Z)),2),
  round(t(cbind(ResSW_MX_Z0$CP_INB_Z, ResSWIMP_MX_Z0$CP_INB_Z, ResSW_MX_Z1$CP_INB_Z, ResSWIMP_MX_Z1$CP_INB_Z)),3)
  )


colnames(results15)=c("Variable","WTP","Method","Subgroup",
                    "EMean INB M&Q","Bias M&Q","SSE M&Q","ESE M&Q","CP M&Q",
                    "EMean INB M&X","Bias M&X","SSE M&X","ESE M&X","CP M&X")




WTP30=20
true.INB.M.Q.Z0=WTP30*true.beta.Q[3]-true.beta.M[3]/1000
true.INB.M.Q.Z1=WTP30*(true.beta.Q[3]+true.beta.Q[4])-(true.beta.M[3]/1000+true.beta.M[4]/1000)
true.INB.M.X.Z0=WTP30*true.beta.TL[3]-true.beta.M[3]/1000
true.INB.M.X.Z1=WTP30*(true.beta.TL[3]+true.beta.TL[4])-(true.beta.M[3]/1000+true.beta.M[4]/1000)

ResSW_MQ_Z0=GetResult(estDiff_SW_Q_Z0,estDiff_SW_M_Z0,true.INB.M.Q.Z0,nsim,WTP30,Syy_SW_Q_Z0,Sxx_SW_M_Z0,CovSW_INB_MQ_Z0,z_alpha)
ResIMP_MQ_Z0=GetResult(estDiff_IMP_Q_Z0,estDiff_IMP_M_Z0,true.INB.M.Q.Z0,nsim,WTP30,Syy_IMP_Q_Z0,Sxx_IMP_M_Z0,CovIMP_INB_MQ_Z0,z_alpha)

ResSW_MQ_Z1=GetResult(estDiff_SW_Q_Z1,estDiff_SW_M_Z1,true.INB.M.Q.Z1,nsim,WTP30,Syy_SW_Q_Z1,Sxx_SW_M_Z1,CovSW_INB_MQ_Z1,z_alpha)
ResIMP_MQ_Z1=GetResult(estDiff_IMP_Q_Z1,estDiff_IMP_M_Z1,true.INB.M.Q.Z1,nsim,WTP30,Syy_IMP_Q_Z1,Sxx_IMP_M_Z1,CovIMP_INB_MQ_Z1,z_alpha)

ResSW_MX_Z0=GetResult(estDiff_SW_X_Z0,estDiff_SW_M_Z0,true.INB.M.X.Z0,nsim,WTP30,Syy_SW_X_Z0,Sxx_SW_M_Z0,CovSW_INB_MX_Z0,z_alpha)
ResSWIMP_MX_Z0=GetResult(estDiff_SW_X_Z0,estDiff_IMP_M_Z0,true.INB.M.X.Z0,nsim,WTP30,Syy_SW_X_Z0,Sxx_IMP_M_Z0,CovSWIMP_INB_MX_Z0,z_alpha)

ResSW_MX_Z1=GetResult(estDiff_SW_X_Z1,estDiff_SW_M_Z1,true.INB.M.X.Z1,nsim,WTP30,Syy_SW_X_Z1,Sxx_SW_M_Z1,CovSW_INB_MX_Z1,z_alpha)
ResSWIMP_MX_Z1=GetResult(estDiff_SW_X_Z1,estDiff_IMP_M_Z1,true.INB.M.X.Z1,nsim,WTP30,Syy_SW_X_Z1,Sxx_IMP_M_Z1,CovSWIMP_INB_MX_Z1,z_alpha)

results30 = data.frame(
  rep(c("INB"),1,each=4),
  rep(c(20),1,each=4),
  rep(c("SW","IMP"),2,each=1),
  rep(c("Z=0","Z=1"),1,each=2),
  round(t(cbind(ResSW_MQ_Z0$mean_INB_Z, ResIMP_MQ_Z0$mean_INB_Z, ResSW_MQ_Z1$mean_INB_Z, ResIMP_MQ_Z1$mean_INB_Z)),2),
  round(t(cbind(ResSW_MQ_Z0$bias_INB_Z, ResIMP_MQ_Z0$bias_INB_Z, ResSW_MQ_Z1$bias_INB_Z, ResIMP_MQ_Z1$bias_INB_Z)),2),
  round(t(cbind(ResSW_MQ_Z0$SSE_INB_Z, ResIMP_MQ_Z0$SSE_INB_Z, ResSW_MQ_Z1$SSE_INB_Z, ResIMP_MQ_Z1$SSE_INB_Z)),2),
  round(t(cbind(ResSW_MQ_Z0$ESE_INB_Z, ResIMP_MQ_Z0$ESE_INB_Z, ResSW_MQ_Z1$ESE_INB_Z, ResIMP_MQ_Z1$ESE_INB_Z)),2),
  round(t(cbind(ResSW_MQ_Z0$CP_INB_Z, ResIMP_MQ_Z0$CP_INB_Z, ResSW_MQ_Z1$CP_INB_Z, ResIMP_MQ_Z1$CP_INB_Z)),3),
  
  round(t(cbind(ResSW_MX_Z0$mean_INB_Z, ResSWIMP_MX_Z0$mean_INB_Z, ResSW_MX_Z1$mean_INB_Z, ResSWIMP_MX_Z1$mean_INB_Z)),2),
  round(t(cbind(ResSW_MX_Z0$bias_INB_Z, ResSWIMP_MX_Z0$bias_INB_Z, ResSW_MX_Z1$bias_INB_Z, ResSWIMP_MX_Z1$bias_INB_Z)),2),
  round(t(cbind(ResSW_MX_Z0$SSE_INB_Z, ResSWIMP_MX_Z0$SSE_INB_Z, ResSW_MX_Z1$SSE_INB_Z, ResSWIMP_MX_Z1$SSE_INB_Z)),2),
  round(t(cbind(ResSW_MX_Z0$ESE_INB_Z, ResSWIMP_MX_Z0$ESE_INB_Z, ResSW_MX_Z1$ESE_INB_Z, ResSWIMP_MX_Z1$ESE_INB_Z)),2),
  round(t(cbind(ResSW_MX_Z0$CP_INB_Z, ResSWIMP_MX_Z0$CP_INB_Z, ResSW_MX_Z1$CP_INB_Z, ResSWIMP_MX_Z1$CP_INB_Z)),3)
)


colnames(results30)=c("Variable","WTP","Method","Subgroup",
                      "EMean INB M&Q","Bias M&Q","SSE M&Q","ESE M&Q","CP M&Q",
                      "EMean INB M&X","Bias M&X","SSE M&X","ESE M&X","CP M&X")

results=rbind(results15,results30)

write.csv(results, paste("/home/dnliu21/code/CEA/simulation/DC/Results/",
                         "Results_INB_",Censorstatus,"_N", n, "_nsim", nsim, "_logn_", rundate,".csv",sep=""), row.names = F)




