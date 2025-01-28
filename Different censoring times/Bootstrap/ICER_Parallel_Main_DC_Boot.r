
##########################
## 		main 	  	##
##########################

library(parallel)

setwd("/home/dnliu21/Code/DCcode/Bootstrap")

source("DataGen_DC_Boot.r")	
source("ICER_Estimate_DC_Boot.r")

rundate = "xxxx_xx_xx"


GetBootEstimate = function(c) {
  print(c)
  # generate data
  data<-sim.CostEff.Obs.V1(n=n,	Coef.Trt=Coef.Trt,
                           Coef.Surv.MainEff=Coef.Surv.MainEff, Coef.Surv.InterEff=Coef.Surv.InterEff,
                           Cost.diag=Cost.diag, Cost.fix=Cost.fix, Cost.rand=Cost.rand, Cost.term=Cost.term,
                           Par.Censor=Par.Censor,	Par.Censor.0=Par.Censor.0,
                           k=k, max.time=max.time)
  
  
  # save the observed censored data
  datacsv=data.frame(id=1:n, survival=data$Follow.up, dead=data$Death, censor=1-data$Death,
                     survival.cost=data$Follow.up.0, dead.cost=data$Death.0, censor.cost=1-data$Death.0,
                     Trt=data$Trt, covariate=data$Covariate,
                     cost=data$Cost.grp.obs, cost_true.term=data$M_term_true,
                     tot.cost=data$M_total,
                     QALY=data$QASurv.grp.obs, tot.QALY=apply(data$QASurv.grp.obs,1,sum),
                     survival_true=data$true_T, censor_time_true=data$true_C, censor_time_cost_true=data$true_C0,
                     cost_true=data$M_true,
                     QALY0=data$QASurv.grp.obs0, tot.QALY0=apply(data$QASurv.grp.obs0,1,sum))

  
  
  for (b in 1:B) {
    resample_order=sample(n,n,replace=TRUE)
    datacsv_Boot <- datacsv[resample_order,]
    
    ICER_All_Boot <- ICER(X0=datacsv_Boot$survival.cost ,delta0=datacsv_Boot$dead.cost, M_total=datacsv_Boot$tot.cost, 
                       M_diag=datacsv_Boot[,10], M_obs=datacsv_Boot[,11:20], M_term_true=datacsv_Boot$cost_true.term,  # cost
                       X=datacsv_Boot$survival, delta=datacsv_Boot$dead, Q_total=datacsv_Boot$tot.QALY, Q_obs=datacsv_Boot[,23:32], 
                       T_true=datacsv_Boot$survival_true, Q_total0=datacsv_Boot$tot.QALY0, Q_obs0=datacsv_Boot[,48:57],  # QAL and survival time
                       Trt=datacsv_Boot$Trt, Z=datacsv_Boot$covariate, max.time=10, type="interaction", covariate=TRUE, interaction=TRUE, Sep.K=FALSE)


    est_SW_M_Bootlist[b,]=ICER_All_Boot$est_SW_M
    est_SW_Q_Bootlist[b,]=ICER_All_Boot$est_SW_Q
    est_SW_X_Bootlist[b,]=ICER_All_Boot$est_SW_X
    est_IMP_M_Bootlist[b,]=ICER_All_Boot$est_IMP_M
    est_IMP_Q_Bootlist[b,]=ICER_All_Boot$est_IMP_Q
    
    ICER_SW_MandQ_Z0_Bootlist[b]=ICER_All_Boot$ICER_SW_MandQ_Z0
    ICER_SW_MandQ_Z1_Bootlist[b]=ICER_All_Boot$ICER_SW_MandQ_Z1
    ICER_IMP_MandQ_Z0_Bootlist[b]=ICER_All_Boot$ICER_IMP_MandQ_Z0
    ICER_IMP_MandQ_Z1_Bootlist[b]=ICER_All_Boot$ICER_IMP_MandQ_Z1
    ICER_SW_MandX_Z0_Bootlist[b]=ICER_All_Boot$ICER_SW_MandX_Z0
    ICER_SW_MandX_Z1_Bootlist[b]=ICER_All_Boot$ICER_SW_MandX_Z1
    ICER_SWIMP_MandX_Z0_Bootlist[b]=ICER_All_Boot$ICER_SWIMP_MandX_Z0
    ICER_SWIMP_MandX_Z1_Bootlist[b]=ICER_All_Boot$ICER_SWIMP_MandX_Z1
    
  }
  
  ICER_All <- ICER(X0=datacsv$survival.cost ,delta0=datacsv$dead.cost, M_total=datacsv$tot.cost, 
                   M_diag=datacsv[,10], M_obs=datacsv[,11:20], M_term_true=datacsv$cost_true.term,  # cost
                   X=datacsv$survival, delta=datacsv$dead, Q_total=datacsv$tot.QALY, Q_obs=datacsv[,23:32], 
                   T_true=datacsv$survival_true, Q_total0=datacsv$tot.QALY0, Q_obs0=datacsv[,48:57],  # QAL and survival time
                   Trt=datacsv$Trt, Z=datacsv$covariate, max.time=10, type="interaction", covariate=TRUE, interaction=TRUE, Sep.K=FALSE)
  
  
  return(list(ICER_SW_MandQ_Z0_Boot=ICER_SW_MandQ_Z0_Bootlist, ICER_SW_MandQ_Z1_Boot=ICER_SW_MandQ_Z1_Bootlist,
              ICER_IMP_MandQ_Z0_Boot=ICER_IMP_MandQ_Z0_Bootlist, ICER_IMP_MandQ_Z1_Boot=ICER_IMP_MandQ_Z1_Bootlist,
              ICER_SW_MandX_Z0_Boot=ICER_SW_MandX_Z0_Bootlist, ICER_SW_MandX_Z1_Boot=ICER_SW_MandX_Z1_Bootlist,
              ICER_SWIMP_MandX_Z0_Boot=ICER_SWIMP_MandX_Z0_Bootlist, ICER_SWIMP_MandX_Z1_Boot=ICER_SWIMP_MandX_Z1_Bootlist,

              est_SW_M=est_SW_M_Bootlist,
              est_SW_Q=est_SW_Q_Bootlist,
              est_SW_X=est_SW_X_Bootlist,
              est_IMP_M=est_IMP_M_Bootlist,
              est_IMP_Q=est_IMP_Q_Bootlist,

              ICER_SW_MandQ_Z0=ICER_All$ICER_SW_MandQ_Z0, ICER_SW_MandQ_Z1=ICER_All$ICER_SW_MandQ_Z1,
              ICER_IMP_MandQ_Z0=ICER_All$ICER_IMP_MandQ_Z0, ICER_IMP_MandQ_Z1=ICER_All$ICER_IMP_MandQ_Z1,
              ICER_SW_MandX_Z0=ICER_All$ICER_SW_MandX_Z0, ICER_SW_MandX_Z1=ICER_All$ICER_SW_MandX_Z1,
              ICER_SWIMP_MandX_Z0=ICER_All$ICER_SWIMP_MandX_Z0, ICER_SWIMP_MandX_Z1=ICER_All$ICER_SWIMP_MandX_Z1
         ))
}


nsim=2000 	# number of simulation replications
B=1000
n=400		#sample size

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



true.ICER.M.Q.Z0=(true.beta.M[3]/1000)/(true.beta.Q[3])
true.ICER.M.Q.Z1=(true.beta.M[3]/1000+true.beta.M[4]/1000)/(true.beta.Q[3]+true.beta.Q[4])

true.ICER.M.X.Z0=(true.beta.M[3]/1000)/(true.beta.TL[3])
true.ICER.M.X.Z1=(true.beta.M[3]/1000+true.beta.M[4]/1000)/(true.beta.TL[3]+true.beta.TL[4])


est_SW_M_Bootlist=est_SW_Q_Bootlist=est_SW_X_Bootlist=est_IMP_M_Bootlist=est_IMP_Q_Bootlist=matrix(NA,B,4)
ICER_SW_MandQ_Z0_Bootlist=ICER_SW_MandQ_Z1_Bootlist=ICER_IMP_MandQ_Z0_Bootlist=ICER_IMP_MandQ_Z1_Bootlist=numeric(B)
ICER_SW_MandX_Z0_Bootlist=ICER_SW_MandX_Z1_Bootlist=ICER_SWIMP_MandX_Z0_Bootlist=ICER_SWIMP_MandX_Z1_Bootlist=numeric(B)


cl <- makeCluster(30)

clusterEvalQ(cl = cl, expr = library(mvtnorm))
clusterEvalQ(cl = cl, expr = library(survival))
clusterEvalQ(cl = cl, expr = library(geepack))

# Will need to export all of these to the cluster
clusterExport(cl = cl, varlist = c("GetBootEstimate", "sim.CostEff.Obs.V1","ICER"))
clusterExport(cl = cl, varlist = c("KM","CalAccum"))
clusterExport(cl = cl, varlist = c("Coef.Surv.MainEff", "Coef.Surv.InterEff","Par.Censor",'Par.Censor.0'))
clusterExport(cl = cl, varlist = c("Coef.Trt", "max.time", "k", "n","nsim","B"))
clusterExport(cl = cl, varlist = c("true.ICER.M.Q.Z0", "true.ICER.M.Q.Z1","true.ICER.M.X.Z0","true.ICER.M.X.Z1"))
clusterExport(cl = cl, varlist = c("Cost.diag","Cost.fix","Cost.rand","Cost.term"))
clusterExport(cl = cl, varlist = c("ICER_SW_MandQ_Z0_Bootlist","ICER_SW_MandQ_Z1_Bootlist","ICER_IMP_MandQ_Z0_Bootlist","ICER_IMP_MandQ_Z1_Bootlist"))
clusterExport(cl = cl, varlist = c("ICER_SW_MandX_Z0_Bootlist","ICER_SW_MandX_Z1_Bootlist","ICER_SWIMP_MandX_Z0_Bootlist","ICER_SWIMP_MandX_Z1_Bootlist"))
clusterExport(cl = cl, varlist = c("est_SW_M_Bootlist","est_SW_Q_Bootlist","est_SW_X_Bootlist","est_IMP_M_Bootlist","est_IMP_Q_Bootlist"))


clusterSetRNGStream(cl = cl, iseed = 12345)
paste0("Started on: ", Sys.time())
start = proc.time()
GetBootEstimateResult = parLapply(cl = cl, X = 1:nsim, fun = GetBootEstimate)
stop = proc.time()
paste0("Ended on: ", Sys.time())
runtime = stop - start
runtime

# Save a copy of the result since it took ~4 hours 
saveRDS(object = GetBootEstimateResult, 
   file = paste("/home/dnliu21/Code/DTcode/Bootstrap/Results/",
    paste0("Results_Bootstrap_ICER_CP_",Censorstatus,"_N", n, "_nsim", nsim, "_logn_", rundate),".RDS",sep=""))



stopCluster(cl)


BootEstimateResult = data.frame(do.call(what = rbind, args = GetBootEstimateResult))
# Bootstrap samples of simulations
ICER_SW_MandQ_Z0_Boot=matrix(unlist(BootEstimateResult$ICER_SW_MandQ_Z0_Boot),nsim,B,byrow=TRUE)
ICER_SW_MandQ_Z1_Boot=matrix(unlist(BootEstimateResult$ICER_SW_MandQ_Z1_Boot),nsim,B,byrow=TRUE)
ICER_IMP_MandQ_Z0_Boot=matrix(unlist(BootEstimateResult$ICER_IMP_MandQ_Z0_Boot),nsim,B,byrow=TRUE)
ICER_IMP_MandQ_Z1_Boot=matrix(unlist(BootEstimateResult$ICER_IMP_MandQ_Z1_Boot),nsim,B,byrow=TRUE)
ICER_SW_MandX_Z0_Boot=matrix(unlist(BootEstimateResult$ICER_SW_MandX_Z0_Boot),nsim,B,byrow=TRUE)
ICER_SW_MandX_Z1_Boot=matrix(unlist(BootEstimateResult$ICER_SW_MandX_Z1_Boot),nsim,B,byrow=TRUE)
ICER_SWIMP_MandX_Z0_Boot=matrix(unlist(BootEstimateResult$ICER_SWIMP_MandX_Z0_Boot),nsim,B,byrow=TRUE)
ICER_SWIMP_MandX_Z1_Boot=matrix(unlist(BootEstimateResult$ICER_SWIMP_MandX_Z1_Boot),nsim,B,byrow=TRUE)
# simulation values
ICER_SW_MandQ_Z0=matrix(unlist(BootEstimateResult$ICER_SW_MandQ_Z0),nsim,1,byrow=TRUE)
ICER_SW_MandQ_Z1=matrix(unlist(BootEstimateResult$ICER_SW_MandQ_Z1),nsim,1,byrow=TRUE)
ICER_IMP_MandQ_Z0=matrix(unlist(BootEstimateResult$ICER_IMP_MandQ_Z0),nsim,1,byrow=TRUE)
ICER_IMP_MandQ_Z1=matrix(unlist(BootEstimateResult$ICER_IMP_MandQ_Z1),nsim,1,byrow=TRUE)
ICER_SW_MandX_Z0=matrix(unlist(BootEstimateResult$ICER_SW_MandX_Z0),nsim,1,byrow=TRUE)
ICER_SW_MandX_Z1=matrix(unlist(BootEstimateResult$ICER_SW_MandX_Z1),nsim,1,byrow=TRUE)
ICER_SWIMP_MandX_Z0=matrix(unlist(BootEstimateResult$ICER_SWIMP_MandX_Z0),nsim,1,byrow=TRUE)
ICER_SWIMP_MandX_Z1=matrix(unlist(BootEstimateResult$ICER_SWIMP_MandX_Z1),nsim,1,byrow=TRUE)


# obtain confidence interval based on Bootstrap
# if the simulation result is within the CI and calculate the CP

# percentile bootstrap CI
CI_percentile_cover<-function(ICER_Boot, true_value, low, up) {
  ICER_Boot_sorted<-sort(ICER_Boot)
  BootCI<-quantile(ICER_Boot_sorted, c(low, up))
  {(true_value<BootCI[2]) & (true_value>BootCI[1])} 
}

# 95%
cover_SW_MandQ_Z0_Boot_percentile95=cover_SW_MandQ_Z1_Boot_percentile95=cover_IMP_MandQ_Z0_Boot_percentile95=cover_IMP_MandQ_Z1_Boot_percentile95=numeric(nsim)
cover_SW_MandX_Z0_Boot_percentile95=cover_SW_MandX_Z1_Boot_percentile95=cover_SWIMP_MandX_Z0_Boot_percentile95=cover_SWIMP_MandX_Z1_Boot_percentile95=numeric(nsim)

for (i in 1:nsim) {
  cover_SW_MandQ_Z0_Boot_percentile95[i]=CI_percentile_cover(ICER_SW_MandQ_Z0_Boot[i,], true.ICER.M.Q.Z0,0.025,0.975)
  cover_SW_MandQ_Z1_Boot_percentile95[i]=CI_percentile_cover(ICER_SW_MandQ_Z1_Boot[i,], true.ICER.M.Q.Z1,0.025,0.975)
  cover_IMP_MandQ_Z0_Boot_percentile95[i]=CI_percentile_cover(ICER_IMP_MandQ_Z0_Boot[i,], true.ICER.M.Q.Z0,0.025,0.975)
  cover_IMP_MandQ_Z1_Boot_percentile95[i]=CI_percentile_cover(ICER_IMP_MandQ_Z1_Boot[i,], true.ICER.M.Q.Z1,0.025,0.975)
  
  cover_SW_MandX_Z0_Boot_percentile95[i]=CI_percentile_cover(ICER_SW_MandX_Z0_Boot[i,],true.ICER.M.X.Z0,0.025,0.975)
  cover_SW_MandX_Z1_Boot_percentile95[i]=CI_percentile_cover(ICER_SW_MandX_Z1_Boot[i,],true.ICER.M.X.Z1,0.025,0.975)
  cover_SWIMP_MandX_Z0_Boot_percentile95[i]=CI_percentile_cover(ICER_SWIMP_MandX_Z0_Boot[i,],true.ICER.M.X.Z0,0.025,0.975)
  cover_SWIMP_MandX_Z1_Boot_percentile95[i]=CI_percentile_cover(ICER_SWIMP_MandX_Z1_Boot[i,],true.ICER.M.X.Z1,0.025,0.975)
}


CP_SW_MandQ_Z0_Boot_percentile95=sum(cover_SW_MandQ_Z0_Boot_percentile95)/nsim
CP_SW_MandQ_Z1_Boot_percentile95=sum(cover_SW_MandQ_Z1_Boot_percentile95)/nsim
CP_IMP_MandQ_Z0_Boot_percentile95=sum(cover_IMP_MandQ_Z0_Boot_percentile95)/nsim
CP_IMP_MandQ_Z1_Boot_percentile95=sum(cover_IMP_MandQ_Z1_Boot_percentile95)/nsim

CP_SW_MandX_Z0_Boot_percentile95=sum(cover_SW_MandX_Z0_Boot_percentile95)/nsim
CP_SW_MandX_Z1_Boot_percentile95=sum(cover_SW_MandX_Z1_Boot_percentile95)/nsim
CP_SWIMP_MandX_Z0_Boot_percentile95=sum(cover_SWIMP_MandX_Z0_Boot_percentile95)/nsim
CP_SWIMP_MandX_Z1_Boot_percentile95=sum(cover_SWIMP_MandX_Z1_Boot_percentile95)/nsim

# 90%
cover_SW_MandQ_Z0_Boot_percentile90=cover_SW_MandQ_Z1_Boot_percentile90=cover_IMP_MandQ_Z0_Boot_percentile90=cover_IMP_MandQ_Z1_Boot_percentile90=numeric(nsim)
cover_SW_MandX_Z0_Boot_percentile90=cover_SW_MandX_Z1_Boot_percentile90=cover_SWIMP_MandX_Z0_Boot_percentile90=cover_SWIMP_MandX_Z1_Boot_percentile90=numeric(nsim)

for (i in 1:nsim) {
  cover_SW_MandQ_Z0_Boot_percentile90[i]=CI_percentile_cover(ICER_SW_MandQ_Z0_Boot[i,], true.ICER.M.Q.Z0,0.05,0.95)
  cover_SW_MandQ_Z1_Boot_percentile90[i]=CI_percentile_cover(ICER_SW_MandQ_Z1_Boot[i,], true.ICER.M.Q.Z1,0.05,0.95)
  cover_IMP_MandQ_Z0_Boot_percentile90[i]=CI_percentile_cover(ICER_IMP_MandQ_Z0_Boot[i,], true.ICER.M.Q.Z0,0.05,0.95)
  cover_IMP_MandQ_Z1_Boot_percentile90[i]=CI_percentile_cover(ICER_IMP_MandQ_Z1_Boot[i,], true.ICER.M.Q.Z1,0.05,0.95)
  
  cover_SW_MandX_Z0_Boot_percentile90[i]=CI_percentile_cover(ICER_SW_MandX_Z0_Boot[i,],true.ICER.M.X.Z0,0.05,0.95)
  cover_SW_MandX_Z1_Boot_percentile90[i]=CI_percentile_cover(ICER_SW_MandX_Z1_Boot[i,],true.ICER.M.X.Z1,0.05,0.95)
  cover_SWIMP_MandX_Z0_Boot_percentile90[i]=CI_percentile_cover(ICER_SWIMP_MandX_Z0_Boot[i,],true.ICER.M.X.Z0,0.05,0.95)
  cover_SWIMP_MandX_Z1_Boot_percentile90[i]=CI_percentile_cover(ICER_SWIMP_MandX_Z1_Boot[i,],true.ICER.M.X.Z1,0.05,0.95)
}

CP_SW_MandQ_Z0_Boot_percentile90=sum(cover_SW_MandQ_Z0_Boot_percentile90)/nsim
CP_SW_MandQ_Z1_Boot_percentile90=sum(cover_SW_MandQ_Z1_Boot_percentile90)/nsim
CP_IMP_MandQ_Z0_Boot_percentile90=sum(cover_IMP_MandQ_Z0_Boot_percentile90)/nsim
CP_IMP_MandQ_Z1_Boot_percentile90=sum(cover_IMP_MandQ_Z1_Boot_percentile90)/nsim

CP_SW_MandX_Z0_Boot_percentile90=sum(cover_SW_MandX_Z0_Boot_percentile90)/nsim
CP_SW_MandX_Z1_Boot_percentile90=sum(cover_SW_MandX_Z1_Boot_percentile90)/nsim
CP_SWIMP_MandX_Z0_Boot_percentile90=sum(cover_SWIMP_MandX_Z0_Boot_percentile90)/nsim
CP_SWIMP_MandX_Z1_Boot_percentile90=sum(cover_SWIMP_MandX_Z1_Boot_percentile90)/nsim



# 80%
cover_SW_MandQ_Z0_Boot_percentile80=cover_SW_MandQ_Z1_Boot_percentile80=cover_IMP_MandQ_Z0_Boot_percentile80=cover_IMP_MandQ_Z1_Boot_percentile80=numeric(nsim)
cover_SW_MandX_Z0_Boot_percentile80=cover_SW_MandX_Z1_Boot_percentile80=cover_SWIMP_MandX_Z0_Boot_percentile80=cover_SWIMP_MandX_Z1_Boot_percentile80=numeric(nsim)

for (i in 1:nsim) {
  cover_SW_MandQ_Z0_Boot_percentile80[i]=CI_percentile_cover(ICER_SW_MandQ_Z0_Boot[i,], true.ICER.M.Q.Z0,0.1,0.9)
  cover_SW_MandQ_Z1_Boot_percentile80[i]=CI_percentile_cover(ICER_SW_MandQ_Z1_Boot[i,], true.ICER.M.Q.Z1,0.1,0.9)
  cover_IMP_MandQ_Z0_Boot_percentile80[i]=CI_percentile_cover(ICER_IMP_MandQ_Z0_Boot[i,], true.ICER.M.Q.Z0,0.1,0.9)
  cover_IMP_MandQ_Z1_Boot_percentile80[i]=CI_percentile_cover(ICER_IMP_MandQ_Z1_Boot[i,], true.ICER.M.Q.Z1,0.1,0.9)
  
  cover_SW_MandX_Z0_Boot_percentile80[i]=CI_percentile_cover(ICER_SW_MandX_Z0_Boot[i,],true.ICER.M.X.Z0,0.1,0.9)
  cover_SW_MandX_Z1_Boot_percentile80[i]=CI_percentile_cover(ICER_SW_MandX_Z1_Boot[i,],true.ICER.M.X.Z1,0.1,0.9)
  cover_SWIMP_MandX_Z0_Boot_percentile80[i]=CI_percentile_cover(ICER_SWIMP_MandX_Z0_Boot[i,],true.ICER.M.X.Z0,0.1,0.9)
  cover_SWIMP_MandX_Z1_Boot_percentile80[i]=CI_percentile_cover(ICER_SWIMP_MandX_Z1_Boot[i,],true.ICER.M.X.Z1,0.1,0.9)
}

CP_SW_MandQ_Z0_Boot_percentile80=sum(cover_SW_MandQ_Z0_Boot_percentile80)/nsim
CP_SW_MandQ_Z1_Boot_percentile80=sum(cover_SW_MandQ_Z1_Boot_percentile80)/nsim
CP_IMP_MandQ_Z0_Boot_percentile80=sum(cover_IMP_MandQ_Z0_Boot_percentile80)/nsim
CP_IMP_MandQ_Z1_Boot_percentile80=sum(cover_IMP_MandQ_Z1_Boot_percentile80)/nsim

CP_SW_MandX_Z0_Boot_percentile80=sum(cover_SW_MandX_Z0_Boot_percentile80)/nsim
CP_SW_MandX_Z1_Boot_percentile80=sum(cover_SW_MandX_Z1_Boot_percentile80)/nsim
CP_SWIMP_MandX_Z0_Boot_percentile80=sum(cover_SWIMP_MandX_Z0_Boot_percentile80)/nsim
CP_SWIMP_MandX_Z1_Boot_percentile80=sum(cover_SWIMP_MandX_Z1_Boot_percentile80)/nsim


# reordered bootstrap CI
CI_reorder_cover<-function(ICER_Boot, true_value, low, up) {
  
  ICER_Boot_sorted<-sort(ICER_Boot)
  
  if (sum(ICER_Boot_sorted<0)==0 | sum(ICER_Boot_sorted>0)==0) {
    BootCI<-quantile(ICER_Boot_sorted, c(low, up))
    output=({(true_value<BootCI[2]) & (true_value>BootCI[1])})
  }
  else {
    pos_idx <- which(ICER_Boot_sorted>0)[1]
    neg_idx <- which(ICER_Boot_sorted<0)[length(which(ICER_Boot_sorted<0))]
    Boot_sorted_reordered=c(sort(ICER_Boot_sorted[1:neg_idx],decreasing=TRUE),sort(ICER_Boot_sorted[pos_idx:length(ICER_Boot_sorted)],decreasing=TRUE))
    
    
    #quant=floor(quantile(1:B,c(low, up)))
    quant=floor(quantile(1:length(Boot_sorted_reordered),c(low, up)))
    BootCI_low=Boot_sorted_reordered[quant[1]]
    BootCI_up=Boot_sorted_reordered[quant[2]]

    # including infinity
    if (BootCI_low>0) {
      output={(true_value>BootCI_up) & (true_value<BootCI_low)}
    }
    if (BootCI_up<0) {
      output={(true_value>BootCI_up) & (true_value<BootCI_low)}
    }
    if (BootCI_low<0 & BootCI_up>0){
      output={(true_value<BootCI_low) | (true_value>BootCI_up)}
    }
    
  }
  output
}


# 95%
cover_SW_MandQ_Z0_Boot_reorder95=cover_SW_MandQ_Z1_Boot_reorder95=cover_IMP_MandQ_Z0_Boot_reorder95=cover_IMP_MandQ_Z1_Boot_reorder95=numeric(nsim)
cover_SW_MandX_Z0_Boot_reorder95=cover_SW_MandX_Z1_Boot_reorder95=cover_SWIMP_MandX_Z0_Boot_reorder95=cover_SWIMP_MandX_Z1_Boot_reorder95=numeric(nsim)

for (i in 1:nsim) {
  cover_SW_MandQ_Z0_Boot_reorder95[i]=CI_reorder_cover(ICER_SW_MandQ_Z0_Boot[i,], true.ICER.M.Q.Z0,0.025,0.975)
  cover_SW_MandQ_Z1_Boot_reorder95[i]=CI_reorder_cover(ICER_SW_MandQ_Z1_Boot[i,], true.ICER.M.Q.Z1,0.025,0.975)
  cover_IMP_MandQ_Z0_Boot_reorder95[i]=CI_reorder_cover(ICER_IMP_MandQ_Z0_Boot[i,], true.ICER.M.Q.Z0,0.025,0.975)
  cover_IMP_MandQ_Z1_Boot_reorder95[i]=CI_reorder_cover(ICER_IMP_MandQ_Z1_Boot[i,], true.ICER.M.Q.Z1,0.025,0.975)
  
  cover_SW_MandX_Z0_Boot_reorder95[i]=CI_reorder_cover(ICER_SW_MandX_Z0_Boot[i,],true.ICER.M.X.Z0,0.025,0.975)
  cover_SW_MandX_Z1_Boot_reorder95[i]=CI_reorder_cover(ICER_SW_MandX_Z1_Boot[i,],true.ICER.M.X.Z1,0.025,0.975)
  cover_SWIMP_MandX_Z0_Boot_reorder95[i]=CI_reorder_cover(ICER_SWIMP_MandX_Z0_Boot[i,],true.ICER.M.X.Z0,0.025,0.975)
  cover_SWIMP_MandX_Z1_Boot_reorder95[i]=CI_reorder_cover(ICER_SWIMP_MandX_Z1_Boot[i,],true.ICER.M.X.Z1,0.025,0.975)
}

CP_SW_MandQ_Z0_Boot_reorder95=sum(cover_SW_MandQ_Z0_Boot_reorder95)/nsim
CP_SW_MandQ_Z1_Boot_reorder95=sum(cover_SW_MandQ_Z1_Boot_reorder95)/nsim
CP_IMP_MandQ_Z0_Boot_reorder95=sum(cover_IMP_MandQ_Z0_Boot_reorder95)/nsim
CP_IMP_MandQ_Z1_Boot_reorder95=sum(cover_IMP_MandQ_Z1_Boot_reorder95)/nsim

CP_SW_MandX_Z0_Boot_reorder95=sum(cover_SW_MandX_Z0_Boot_reorder95)/nsim
CP_SW_MandX_Z1_Boot_reorder95=sum(cover_SW_MandX_Z1_Boot_reorder95)/nsim
CP_SWIMP_MandX_Z0_Boot_reorder95=sum(cover_SWIMP_MandX_Z0_Boot_reorder95)/nsim
CP_SWIMP_MandX_Z1_Boot_reorder95=sum(cover_SWIMP_MandX_Z1_Boot_reorder95)/nsim

# 90%
cover_SW_MandQ_Z0_Boot_reorder90=cover_SW_MandQ_Z1_Boot_reorder90=cover_IMP_MandQ_Z0_Boot_reorder90=cover_IMP_MandQ_Z1_Boot_reorder90=numeric(nsim)
cover_SW_MandX_Z0_Boot_reorder90=cover_SW_MandX_Z1_Boot_reorder90=cover_SWIMP_MandX_Z0_Boot_reorder90=cover_SWIMP_MandX_Z1_Boot_reorder90=numeric(nsim)

for (i in 1:nsim) {
  cover_SW_MandQ_Z0_Boot_reorder90[i]=CI_reorder_cover(ICER_SW_MandQ_Z0_Boot[i,], true.ICER.M.Q.Z0,0.05,0.95)
  cover_SW_MandQ_Z1_Boot_reorder90[i]=CI_reorder_cover(ICER_SW_MandQ_Z1_Boot[i,], true.ICER.M.Q.Z1,0.05,0.95)
  cover_IMP_MandQ_Z0_Boot_reorder90[i]=CI_reorder_cover(ICER_IMP_MandQ_Z0_Boot[i,], true.ICER.M.Q.Z0,0.05,0.95)
  cover_IMP_MandQ_Z1_Boot_reorder90[i]=CI_reorder_cover(ICER_IMP_MandQ_Z1_Boot[i,], true.ICER.M.Q.Z1,0.05,0.95)
  
  cover_SW_MandX_Z0_Boot_reorder90[i]=CI_reorder_cover(ICER_SW_MandX_Z0_Boot[i,],true.ICER.M.X.Z0,0.05,0.95)
  cover_SW_MandX_Z1_Boot_reorder90[i]=CI_reorder_cover(ICER_SW_MandX_Z1_Boot[i,],true.ICER.M.X.Z1,0.05,0.95)
  cover_SWIMP_MandX_Z0_Boot_reorder90[i]=CI_reorder_cover(ICER_SWIMP_MandX_Z0_Boot[i,],true.ICER.M.X.Z0,0.05,0.95)
  cover_SWIMP_MandX_Z1_Boot_reorder90[i]=CI_reorder_cover(ICER_SWIMP_MandX_Z1_Boot[i,],true.ICER.M.X.Z1,0.05,0.95)
}

CP_SW_MandQ_Z0_Boot_reorder90=sum(cover_SW_MandQ_Z0_Boot_reorder90)/nsim
CP_SW_MandQ_Z1_Boot_reorder90=sum(cover_SW_MandQ_Z1_Boot_reorder90)/nsim
CP_IMP_MandQ_Z0_Boot_reorder90=sum(cover_IMP_MandQ_Z0_Boot_reorder90)/nsim
CP_IMP_MandQ_Z1_Boot_reorder90=sum(cover_IMP_MandQ_Z1_Boot_reorder90)/nsim

CP_SW_MandX_Z0_Boot_reorder90=sum(cover_SW_MandX_Z0_Boot_reorder90)/nsim
CP_SW_MandX_Z1_Boot_reorder90=sum(cover_SW_MandX_Z1_Boot_reorder90)/nsim
CP_SWIMP_MandX_Z0_Boot_reorder90=sum(cover_SWIMP_MandX_Z0_Boot_reorder90)/nsim
CP_SWIMP_MandX_Z1_Boot_reorder90=sum(cover_SWIMP_MandX_Z1_Boot_reorder90)/nsim



# 80%
cover_SW_MandQ_Z0_Boot_reorder80=cover_SW_MandQ_Z1_Boot_reorder80=cover_IMP_MandQ_Z0_Boot_reorder80=cover_IMP_MandQ_Z1_Boot_reorder80=numeric(nsim)
cover_SW_MandX_Z0_Boot_reorder80=cover_SW_MandX_Z1_Boot_reorder80=cover_SWIMP_MandX_Z0_Boot_reorder80=cover_SWIMP_MandX_Z1_Boot_reorder80=numeric(nsim)

for (i in 1:nsim) {
  cover_SW_MandQ_Z0_Boot_reorder80[i]=CI_reorder_cover(ICER_SW_MandQ_Z0_Boot[i,], true.ICER.M.Q.Z0,0.1,0.9)
  cover_SW_MandQ_Z1_Boot_reorder80[i]=CI_reorder_cover(ICER_SW_MandQ_Z1_Boot[i,], true.ICER.M.Q.Z1,0.1,0.9)
  cover_IMP_MandQ_Z0_Boot_reorder80[i]=CI_reorder_cover(ICER_IMP_MandQ_Z0_Boot[i,], true.ICER.M.Q.Z0,0.1,0.9)
  cover_IMP_MandQ_Z1_Boot_reorder80[i]=CI_reorder_cover(ICER_IMP_MandQ_Z1_Boot[i,], true.ICER.M.Q.Z1,0.1,0.9)
  
  cover_SW_MandX_Z0_Boot_reorder80[i]=CI_reorder_cover(ICER_SW_MandX_Z0_Boot[i,],true.ICER.M.X.Z0,0.1,0.9)
  cover_SW_MandX_Z1_Boot_reorder80[i]=CI_reorder_cover(ICER_SW_MandX_Z1_Boot[i,],true.ICER.M.X.Z1,0.1,0.9)
  cover_SWIMP_MandX_Z0_Boot_reorder80[i]=CI_reorder_cover(ICER_SWIMP_MandX_Z0_Boot[i,],true.ICER.M.X.Z0,0.1,0.9)
  cover_SWIMP_MandX_Z1_Boot_reorder80[i]=CI_reorder_cover(ICER_SWIMP_MandX_Z1_Boot[i,],true.ICER.M.X.Z1,0.1,0.9)
}

CP_SW_MandQ_Z0_Boot_reorder80=sum(cover_SW_MandQ_Z0_Boot_reorder80)/nsim
CP_SW_MandQ_Z1_Boot_reorder80=sum(cover_SW_MandQ_Z1_Boot_reorder80)/nsim
CP_IMP_MandQ_Z0_Boot_reorder80=sum(cover_IMP_MandQ_Z0_Boot_reorder80)/nsim
CP_IMP_MandQ_Z1_Boot_reorder80=sum(cover_IMP_MandQ_Z1_Boot_reorder80)/nsim

CP_SW_MandX_Z0_Boot_reorder80=sum(cover_SW_MandX_Z0_Boot_reorder80)/nsim
CP_SW_MandX_Z1_Boot_reorder80=sum(cover_SW_MandX_Z1_Boot_reorder80)/nsim
CP_SWIMP_MandX_Z0_Boot_reorder80=sum(cover_SWIMP_MandX_Z0_Boot_reorder80)/nsim
CP_SWIMP_MandX_Z1_Boot_reorder80=sum(cover_SWIMP_MandX_Z1_Boot_reorder80)/nsim



results = data.frame(
  rep(c("ICER"),1,each=4),
  rep(c("SW","IMP"),2,each=1),
  rep(c("Z=0","Z=1"),1,each=2),
  round(t(cbind(CP_SW_MandQ_Z0_Boot_percentile95, CP_IMP_MandQ_Z0_Boot_percentile95,CP_SW_MandQ_Z1_Boot_percentile95, CP_IMP_MandQ_Z1_Boot_percentile95)),4),
  round(t(cbind(CP_SW_MandQ_Z0_Boot_reorder95, CP_IMP_MandQ_Z0_Boot_reorder95,CP_SW_MandQ_Z1_Boot_reorder95, CP_IMP_MandQ_Z1_Boot_reorder95)),4),
  round(t(cbind(CP_SW_MandQ_Z0_Boot_percentile90, CP_IMP_MandQ_Z0_Boot_percentile90,CP_SW_MandQ_Z1_Boot_percentile90, CP_IMP_MandQ_Z1_Boot_percentile90)),4),
  round(t(cbind(CP_SW_MandQ_Z0_Boot_reorder90, CP_IMP_MandQ_Z0_Boot_reorder90,CP_SW_MandQ_Z1_Boot_reorder90, CP_IMP_MandQ_Z1_Boot_reorder90)),4),
  round(t(cbind(CP_SW_MandQ_Z0_Boot_percentile80, CP_IMP_MandQ_Z0_Boot_percentile80,CP_SW_MandQ_Z1_Boot_percentile80, CP_IMP_MandQ_Z1_Boot_percentile80)),4),
  round(t(cbind(CP_SW_MandQ_Z0_Boot_reorder80, CP_IMP_MandQ_Z0_Boot_reorder80,CP_SW_MandQ_Z1_Boot_reorder80, CP_IMP_MandQ_Z1_Boot_reorder80)),4),
  round(t(cbind(CP_SW_MandX_Z0_Boot_percentile95, CP_SWIMP_MandX_Z0_Boot_percentile95,CP_SW_MandX_Z1_Boot_percentile95, CP_SWIMP_MandX_Z1_Boot_percentile95)),4),
  round(t(cbind(CP_SW_MandX_Z0_Boot_reorder95, CP_SWIMP_MandX_Z0_Boot_reorder95,CP_SW_MandX_Z1_Boot_reorder95, CP_SWIMP_MandX_Z1_Boot_reorder95)),4),
  round(t(cbind(CP_SW_MandX_Z0_Boot_percentile90, CP_SWIMP_MandX_Z0_Boot_percentile90,CP_SW_MandX_Z1_Boot_percentile90, CP_SWIMP_MandX_Z1_Boot_percentile90)),4),
  round(t(cbind(CP_SW_MandX_Z0_Boot_reorder90, CP_SWIMP_MandX_Z0_Boot_reorder90,CP_SW_MandX_Z1_Boot_reorder90, CP_SWIMP_MandX_Z1_Boot_reorder90)),4),
  round(t(cbind(CP_SW_MandX_Z0_Boot_percentile80, CP_SWIMP_MandX_Z0_Boot_percentile80,CP_SW_MandX_Z1_Boot_percentile80, CP_SWIMP_MandX_Z1_Boot_percentile80)),4),
  round(t(cbind(CP_SW_MandX_Z0_Boot_reorder80, CP_SWIMP_MandX_Z0_Boot_reorder80,CP_SW_MandX_Z1_Boot_reorder80, CP_SWIMP_MandX_Z1_Boot_reorder80)),4)
)




colnames(results)=c("Variable","Method","Subgroup",
                    "M & Q - Percentile 95%", "M & Q - Reordered 95%", "M & Q - Percentile 90%", "M & Q - Reordered 90%","M & Q - Percentile 80%", "M & Q - Reordered 80%",
                    "M & X - Percentile 95%", "M & X - Reordered 95%", "M & X - Percentile 90%", "M & X - Reordered 90%", "M & X - Percentile 80%", "M & X - Reordered 80%")



write.csv(results, paste("/home/dnliu21/Code/DTcode/Bootstrap/Results/",
            "Results_Bootstrap_ICER_CP_",Censorstatus,"_N", n, "_nsim", nsim, "_logn_", rundate,".csv",sep=""), row.names = F)



