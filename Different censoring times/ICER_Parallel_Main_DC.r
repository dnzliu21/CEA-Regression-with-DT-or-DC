
##########################
## 		main 	  	##
##########################

library(parallel)

setwd("/home/dnliu21/code/CEA/simulation/DC")

source("DataGen_DC.r")	
source("ICER_Estimate_DC.r")
#source("ICER_Estimate_DC_censorearly.r")

rundate = "xxxx_xx_xx"

GetEstimate = function(b) {
  print(b)
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
  

  
  ICER_All <- ICER(X0=datacsv$survival.cost ,delta0=datacsv$dead.cost, M_total=datacsv$tot.cost, 
    M_diag=datacsv[,10], M_obs=datacsv[,11:20], M_term_true=datacsv$cost_true.term,  # cost
    X=datacsv$survival, delta=datacsv$dead, Q_total=datacsv$tot.QALY, Q_obs=datacsv[,23:32], 
    T_true=datacsv$survival_true, Q_total0=datacsv$tot.QALY0, Q_obs0=datacsv[,48:57],  # QAL and survival time
    Trt=datacsv$Trt, Z=datacsv$covariate, max.time=10, type="interaction", covariate=TRUE, interaction=TRUE, Sep.K=FALSE)
  
  
  # 95%
  coverSW_MandQ_Z095_Fieller=CalCoverage_Fieller(ICER_All$ICER_CI_SW_MandQ_Z095_low,ICER_All$ICER_CI_SW_MandQ_Z095_up,ICER_All$ICER_CI_SW_MandQ_Z095_revert,true.ICER.M.Q.Z0)
  coverSW_MandQ_Z195_Fieller=CalCoverage_Fieller(ICER_All$ICER_CI_SW_MandQ_Z195_low,ICER_All$ICER_CI_SW_MandQ_Z195_up,ICER_All$ICER_CI_SW_MandQ_Z195_revert,true.ICER.M.Q.Z1)
  coverIMP_MandQ_Z095_Fieller=CalCoverage_Fieller(ICER_All$ICER_CI_IMP_MandQ_Z095_low,ICER_All$ICER_CI_IMP_MandQ_Z095_up,ICER_All$ICER_CI_IMP_MandQ_Z095_revert,true.ICER.M.Q.Z0)
  coverIMP_MandQ_Z195_Fieller=CalCoverage_Fieller(ICER_All$ICER_CI_IMP_MandQ_Z195_low,ICER_All$ICER_CI_IMP_MandQ_Z195_up,ICER_All$ICER_CI_IMP_MandQ_Z195_revert,true.ICER.M.Q.Z1)
  
  coverSW_MandX_Z095_Fieller=CalCoverage_Fieller(ICER_All$ICER_CI_SW_MandX_Z095_low,ICER_All$ICER_CI_SW_MandX_Z095_up,ICER_All$ICER_CI_SW_MandX_Z095_revert,true.ICER.M.X.Z0)
  coverSW_MandX_Z195_Fieller=CalCoverage_Fieller(ICER_All$ICER_CI_SW_MandX_Z195_low,ICER_All$ICER_CI_SW_MandX_Z195_up,ICER_All$ICER_CI_SW_MandX_Z195_revert,true.ICER.M.X.Z1)
  coverSWIMP_MandX_Z095_Fieller=CalCoverage_Fieller(ICER_All$ICER_CI_SWIMP_MandX_Z095_low,ICER_All$ICER_CI_SWIMP_MandX_Z095_up,ICER_All$ICER_CI_SWIMP_MandX_Z095_revert,true.ICER.M.X.Z0)
  coverSWIMP_MandX_Z195_Fieller=CalCoverage_Fieller(ICER_All$ICER_CI_SWIMP_MandX_Z195_low,ICER_All$ICER_CI_SWIMP_MandX_Z195_up,ICER_All$ICER_CI_SWIMP_MandX_Z195_revert,true.ICER.M.X.Z1)
  
  
  coverSW_MandQ_Z095_Delta=CalCoverage_Delta(ICER_All$ICER_CI_SW_MandQ_Z095_Delta_low,ICER_All$ICER_CI_SW_MandQ_Z095_Delta_up,true.ICER.M.Q.Z0)
  coverSW_MandQ_Z195_Delta=CalCoverage_Delta(ICER_All$ICER_CI_SW_MandQ_Z195_Delta_low,ICER_All$ICER_CI_SW_MandQ_Z195_Delta_up,true.ICER.M.Q.Z1)
  coverIMP_MandQ_Z095_Delta=CalCoverage_Delta(ICER_All$ICER_CI_IMP_MandQ_Z095_Delta_low,ICER_All$ICER_CI_IMP_MandQ_Z095_Delta_up,true.ICER.M.Q.Z0)
  coverIMP_MandQ_Z195_Delta=CalCoverage_Delta(ICER_All$ICER_CI_IMP_MandQ_Z195_Delta_low,ICER_All$ICER_CI_IMP_MandQ_Z195_Delta_up,true.ICER.M.Q.Z1)
  
  coverSW_MandX_Z095_Delta=CalCoverage_Delta(ICER_All$ICER_CI_SW_MandX_Z095_Delta_low,ICER_All$ICER_CI_SW_MandX_Z095_Delta_up,true.ICER.M.X.Z0)
  coverSW_MandX_Z195_Delta=CalCoverage_Delta(ICER_All$ICER_CI_SW_MandX_Z195_Delta_low,ICER_All$ICER_CI_SW_MandX_Z195_Delta_up,true.ICER.M.X.Z1)
  coverSWIMP_MandX_Z095_Delta=CalCoverage_Delta(ICER_All$ICER_CI_SWIMP_MandX_Z095_Delta_low,ICER_All$ICER_CI_SWIMP_MandX_Z095_Delta_up,true.ICER.M.X.Z0)
  coverSWIMP_MandX_Z195_Delta=CalCoverage_Delta(ICER_All$ICER_CI_SWIMP_MandX_Z195_Delta_low,ICER_All$ICER_CI_SWIMP_MandX_Z195_Delta_up,true.ICER.M.X.Z1)
  
  
  # 90%
  coverSW_MandQ_Z090_Fieller=CalCoverage_Fieller(ICER_All$ICER_CI_SW_MandQ_Z090_low,ICER_All$ICER_CI_SW_MandQ_Z090_up,ICER_All$ICER_CI_SW_MandQ_Z090_revert,true.ICER.M.Q.Z0)
  coverSW_MandQ_Z190_Fieller=CalCoverage_Fieller(ICER_All$ICER_CI_SW_MandQ_Z190_low,ICER_All$ICER_CI_SW_MandQ_Z190_up,ICER_All$ICER_CI_SW_MandQ_Z190_revert,true.ICER.M.Q.Z1)
  coverIMP_MandQ_Z090_Fieller=CalCoverage_Fieller(ICER_All$ICER_CI_IMP_MandQ_Z090_low,ICER_All$ICER_CI_IMP_MandQ_Z090_up,ICER_All$ICER_CI_IMP_MandQ_Z090_revert,true.ICER.M.Q.Z0)
  coverIMP_MandQ_Z190_Fieller=CalCoverage_Fieller(ICER_All$ICER_CI_IMP_MandQ_Z190_low,ICER_All$ICER_CI_IMP_MandQ_Z190_up,ICER_All$ICER_CI_IMP_MandQ_Z190_revert,true.ICER.M.Q.Z1)
  
  coverSW_MandX_Z090_Fieller=CalCoverage_Fieller(ICER_All$ICER_CI_SW_MandX_Z090_low,ICER_All$ICER_CI_SW_MandX_Z090_up,ICER_All$ICER_CI_SW_MandX_Z090_revert,true.ICER.M.X.Z0)
  coverSW_MandX_Z190_Fieller=CalCoverage_Fieller(ICER_All$ICER_CI_SW_MandX_Z190_low,ICER_All$ICER_CI_SW_MandX_Z190_up,ICER_All$ICER_CI_SW_MandX_Z190_revert,true.ICER.M.X.Z1)
  coverSWIMP_MandX_Z090_Fieller=CalCoverage_Fieller(ICER_All$ICER_CI_SWIMP_MandX_Z090_low,ICER_All$ICER_CI_SWIMP_MandX_Z090_up,ICER_All$ICER_CI_SWIMP_MandX_Z090_revert,true.ICER.M.X.Z0)
  coverSWIMP_MandX_Z190_Fieller=CalCoverage_Fieller(ICER_All$ICER_CI_SWIMP_MandX_Z190_low,ICER_All$ICER_CI_SWIMP_MandX_Z190_up,ICER_All$ICER_CI_SWIMP_MandX_Z190_revert,true.ICER.M.X.Z1)
  
  
  coverSW_MandQ_Z090_Delta=CalCoverage_Delta(ICER_All$ICER_CI_SW_MandQ_Z090_Delta_low,ICER_All$ICER_CI_SW_MandQ_Z090_Delta_up,true.ICER.M.Q.Z0)
  coverSW_MandQ_Z190_Delta=CalCoverage_Delta(ICER_All$ICER_CI_SW_MandQ_Z190_Delta_low,ICER_All$ICER_CI_SW_MandQ_Z190_Delta_up,true.ICER.M.Q.Z1)
  coverIMP_MandQ_Z090_Delta=CalCoverage_Delta(ICER_All$ICER_CI_IMP_MandQ_Z090_Delta_low,ICER_All$ICER_CI_IMP_MandQ_Z090_Delta_up,true.ICER.M.Q.Z0)
  coverIMP_MandQ_Z190_Delta=CalCoverage_Delta(ICER_All$ICER_CI_IMP_MandQ_Z190_Delta_low,ICER_All$ICER_CI_IMP_MandQ_Z190_Delta_up,true.ICER.M.Q.Z1)
  
  coverSW_MandX_Z090_Delta=CalCoverage_Delta(ICER_All$ICER_CI_SW_MandX_Z090_Delta_low,ICER_All$ICER_CI_SW_MandX_Z090_Delta_up,true.ICER.M.X.Z0)
  coverSW_MandX_Z190_Delta=CalCoverage_Delta(ICER_All$ICER_CI_SW_MandX_Z190_Delta_low,ICER_All$ICER_CI_SW_MandX_Z190_Delta_up,true.ICER.M.X.Z1)
  coverSWIMP_MandX_Z090_Delta=CalCoverage_Delta(ICER_All$ICER_CI_SWIMP_MandX_Z090_Delta_low,ICER_All$ICER_CI_SWIMP_MandX_Z090_Delta_up,true.ICER.M.X.Z0)
  coverSWIMP_MandX_Z190_Delta=CalCoverage_Delta(ICER_All$ICER_CI_SWIMP_MandX_Z190_Delta_low,ICER_All$ICER_CI_SWIMP_MandX_Z190_Delta_up,true.ICER.M.X.Z1)
  
  
  # 80%
  coverSW_MandQ_Z080_Fieller=CalCoverage_Fieller(ICER_All$ICER_CI_SW_MandQ_Z080_low,ICER_All$ICER_CI_SW_MandQ_Z080_up,ICER_All$ICER_CI_SW_MandQ_Z080_revert,true.ICER.M.Q.Z0)
  coverSW_MandQ_Z180_Fieller=CalCoverage_Fieller(ICER_All$ICER_CI_SW_MandQ_Z180_low,ICER_All$ICER_CI_SW_MandQ_Z180_up,ICER_All$ICER_CI_SW_MandQ_Z180_revert,true.ICER.M.Q.Z1)
  coverIMP_MandQ_Z080_Fieller=CalCoverage_Fieller(ICER_All$ICER_CI_IMP_MandQ_Z080_low,ICER_All$ICER_CI_IMP_MandQ_Z080_up,ICER_All$ICER_CI_IMP_MandQ_Z080_revert,true.ICER.M.Q.Z0)
  coverIMP_MandQ_Z180_Fieller=CalCoverage_Fieller(ICER_All$ICER_CI_IMP_MandQ_Z180_low,ICER_All$ICER_CI_IMP_MandQ_Z180_up,ICER_All$ICER_CI_IMP_MandQ_Z180_revert,true.ICER.M.Q.Z1)
  
  coverSW_MandX_Z080_Fieller=CalCoverage_Fieller(ICER_All$ICER_CI_SW_MandX_Z080_low,ICER_All$ICER_CI_SW_MandX_Z080_up,ICER_All$ICER_CI_SW_MandX_Z080_revert,true.ICER.M.X.Z0)
  coverSW_MandX_Z180_Fieller=CalCoverage_Fieller(ICER_All$ICER_CI_SW_MandX_Z180_low,ICER_All$ICER_CI_SW_MandX_Z180_up,ICER_All$ICER_CI_SW_MandX_Z180_revert,true.ICER.M.X.Z1)
  coverSWIMP_MandX_Z080_Fieller=CalCoverage_Fieller(ICER_All$ICER_CI_SWIMP_MandX_Z080_low,ICER_All$ICER_CI_SWIMP_MandX_Z080_up,ICER_All$ICER_CI_SWIMP_MandX_Z080_revert,true.ICER.M.X.Z0)
  coverSWIMP_MandX_Z180_Fieller=CalCoverage_Fieller(ICER_All$ICER_CI_SWIMP_MandX_Z180_low,ICER_All$ICER_CI_SWIMP_MandX_Z180_up,ICER_All$ICER_CI_SWIMP_MandX_Z180_revert,true.ICER.M.X.Z1)
  
  
  coverSW_MandQ_Z080_Delta=CalCoverage_Delta(ICER_All$ICER_CI_SW_MandQ_Z080_Delta_low,ICER_All$ICER_CI_SW_MandQ_Z080_Delta_up,true.ICER.M.Q.Z0)
  coverSW_MandQ_Z180_Delta=CalCoverage_Delta(ICER_All$ICER_CI_SW_MandQ_Z180_Delta_low,ICER_All$ICER_CI_SW_MandQ_Z180_Delta_up,true.ICER.M.Q.Z1)
  coverIMP_MandQ_Z080_Delta=CalCoverage_Delta(ICER_All$ICER_CI_IMP_MandQ_Z080_Delta_low,ICER_All$ICER_CI_IMP_MandQ_Z080_Delta_up,true.ICER.M.Q.Z0)
  coverIMP_MandQ_Z180_Delta=CalCoverage_Delta(ICER_All$ICER_CI_IMP_MandQ_Z180_Delta_low,ICER_All$ICER_CI_IMP_MandQ_Z180_Delta_up,true.ICER.M.Q.Z1)
  
  coverSW_MandX_Z080_Delta=CalCoverage_Delta(ICER_All$ICER_CI_SW_MandX_Z080_Delta_low,ICER_All$ICER_CI_SW_MandX_Z080_Delta_up,true.ICER.M.X.Z0)
  coverSW_MandX_Z180_Delta=CalCoverage_Delta(ICER_All$ICER_CI_SW_MandX_Z180_Delta_low,ICER_All$ICER_CI_SW_MandX_Z180_Delta_up,true.ICER.M.X.Z1)
  coverSWIMP_MandX_Z080_Delta=CalCoverage_Delta(ICER_All$ICER_CI_SWIMP_MandX_Z080_Delta_low,ICER_All$ICER_CI_SWIMP_MandX_Z080_Delta_up,true.ICER.M.X.Z0)
  coverSWIMP_MandX_Z180_Delta=CalCoverage_Delta(ICER_All$ICER_CI_SWIMP_MandX_Z180_Delta_low,ICER_All$ICER_CI_SWIMP_MandX_Z180_Delta_up,true.ICER.M.X.Z1)
  
  
  
  return(list(ICER_SW_MandQ_Z0=ICER_All$ICER_SW_MandQ_Z0, ICER_SW_MandQ_Z1=ICER_All$ICER_SW_MandQ_Z1,
              ICER_IMP_MandQ_Z0=ICER_All$ICER_IMP_MandQ_Z0, ICER_IMP_MandQ_Z1=ICER_All$ICER_IMP_MandQ_Z1,
              ICER_SW_MandX_Z0=ICER_All$ICER_SW_MandX_Z0, ICER_SW_MandX_Z1=ICER_All$ICER_SW_MandX_Z1,
              ICER_SWIMP_MandX_Z0=ICER_All$ICER_SWIMP_MandX_Z0, ICER_SWIMP_MandX_Z1=ICER_All$ICER_SWIMP_MandX_Z1,
              
              # estimated covariance between cost and QAL for group Z=0 and Z=1
              ICER_SW_MandQ_Cov_Z0=ICER_All$Cov_SW_MandQ[3,3], 
              ICER_SW_MandQ_Cov_Z1=ICER_All$Cov_SW_MandQ[3,3]+ICER_All$Cov_SW_MandQ[3,4]+ICER_All$Cov_SW_MandQ[4,3]+ICER_All$Cov_SW_MandQ[4,4],
              ICER_IMP_MandQ_Cov_Z0=ICER_All$Cov_IMP_MandQ[3,3], 
              ICER_IMP_MandQ_Cov_Z1=ICER_All$Cov_IMP_MandQ[3,3]+ICER_All$Cov_IMP_MandQ[3,4]+ICER_All$Cov_IMP_MandQ[4,3]+ICER_All$Cov_IMP_MandQ[4,4],

              ICER_SW_MandX_Cov_Z0=ICER_All$Cov_SW_MandX[3,3], 
              ICER_SW_MandX_Cov_Z1=ICER_All$Cov_SW_MandX[3,3]+ICER_All$Cov_SW_MandX[3,4]+ICER_All$Cov_SW_MandX[4,3]+ICER_All$Cov_SW_MandX[4,4],
              ICER_SWIMP_MandX_Cov_Z0=ICER_All$Cov_SWIMP_MandX[3,3], 
              ICER_SWIMP_MandX_Cov_Z1=ICER_All$Cov_SWIMP_MandX[3,3]+ICER_All$Cov_SWIMP_MandX[3,4]+ICER_All$Cov_SWIMP_MandX[4,3]+ICER_All$Cov_SWIMP_MandX[4,4],
              
              
              Sxx_SW_M_Z0=ICER_All$Cov_SW_M[3,3], 
              Syy_SW_Q_Z0=ICER_All$Cov_SW_Q[3,3], 
              Syy_SW_X_Z0=ICER_All$Cov_SW_X[3,3],
              Sxx_SW_M_Z1=ICER_All$Cov_SW_M[3,3]+ICER_All$Cov_SW_M[4,4]+2*ICER_All$Cov_SW_M[3,4], 
              Syy_SW_Q_Z1=ICER_All$Cov_SW_Q[3,3]+ICER_All$Cov_SW_Q[4,4]+2*ICER_All$Cov_SW_Q[3,4], 
              Syy_SW_X_Z1=ICER_All$Cov_SW_X[3,3]+ICER_All$Cov_SW_X[4,4]+2*ICER_All$Cov_SW_X[3,4],
              
              Sxx_IMP_M_Z0=ICER_All$Cov_IMP_M[3,3], 
              Syy_IMP_Q_Z0=ICER_All$Cov_IMP_Q[3,3], 
              Sxx_IMP_M_Z1=ICER_All$Cov_IMP_M[3,3]+ICER_All$Cov_IMP_M[4,4]+2*ICER_All$Cov_IMP_M[3,4], 
              Syy_IMP_Q_Z1=ICER_All$Cov_IMP_Q[3,3]+ICER_All$Cov_IMP_Q[4,4]+2*ICER_All$Cov_IMP_Q[3,4], 
              
              # obtain estimated beta for calculating sample variance
              estDiff_SW_M_Z0=ICER_All$est_SW_M[3], 
              estDiff_SW_Q_Z0=ICER_All$est_SW_Q[3], 
              estDiff_SW_X_Z0=ICER_All$est_SW_X[3],
              estDiff_SW_M_Z1=ICER_All$est_SW_M[3]+ICER_All$est_SW_M[4],
              estDiff_SW_Q_Z1=ICER_All$est_SW_Q[3]+ICER_All$est_SW_Q[4],
              estDiff_SW_X_Z1=ICER_All$est_SW_X[3]+ICER_All$est_SW_X[4],
              
              estDiff_IMP_M_Z0=ICER_All$est_IMP_M[3], 
              estDiff_IMP_Q_Z0=ICER_All$est_IMP_Q[3], 
              estDiff_IMP_M_Z1=ICER_All$est_IMP_M[3]+ICER_All$est_IMP_M[4],
              estDiff_IMP_Q_Z1=ICER_All$est_IMP_Q[3]+ICER_All$est_IMP_Q[4],
             
              # confidence limits
              ICER_CI_SW_MandQ_Z095_low=ICER_All$ICER_CI_SW_MandQ_Z095_low, ICER_CI_SW_MandQ_Z095_up=ICER_All$ICER_CI_SW_MandQ_Z095_up,
              ICER_CI_SW_MandQ_Z095_revert=ICER_All$ICER_CI_SW_MandQ_Z095_revert,
              
              ICER_CI_SW_MandQ_Z195_low=ICER_All$ICER_CI_SW_MandQ_Z195_low, ICER_CI_SW_MandQ_Z195_up=ICER_All$ICER_CI_SW_MandQ_Z195_up,
              ICER_CI_SW_MandQ_Z195_revert=ICER_All$ICER_CI_SW_MandQ_Z195_revert,
              
              ICER_CI_IMP_MandQ_Z095_low=ICER_All$ICER_CI_IMP_MandQ_Z095_low, ICER_CI_IMP_MandQ_Z095_up=ICER_All$ICER_CI_IMP_MandQ_Z095_up,
              ICER_CI_IMP_MandQ_Z095_revert=ICER_All$ICER_CI_IMP_MandQ_Z095_revert,
              
              ICER_CI_IMP_MandQ_Z195_low=ICER_All$ICER_CI_IMP_MandQ_Z195_low, ICER_CI_IMP_MandQ_Z195_up=ICER_All$ICER_CI_IMP_MandQ_Z195_up,
              ICER_CI_IMP_MandQ_Z195_revert=ICER_All$ICER_CI_IMP_MandQ_Z195_revert,
              
              ICER_CI_SW_MandX_Z095_low=ICER_All$ICER_CI_SW_MandX_Z095_low, ICER_CI_SW_MandX_Z095_up=ICER_All$ICER_CI_SW_MandX_Z095_up,
              ICER_CI_SW_MandX_Z095_revert=ICER_All$ICER_CI_SW_MandX_Z095_revert,
              
              ICER_CI_SW_MandX_Z195_low=ICER_All$ICER_CI_SW_MandX_Z195_low,ICER_CI_SW_MandX_Z195_up=ICER_All$ICER_CI_SW_MandX_Z195_up,
              ICER_CI_SW_MandX_Z195_revert=ICER_All$ICER_CI_SW_MandX_Z195_revert,
              
              ICER_CI_SWIMP_MandX_Z095_low=ICER_All$ICER_CI_SWIMP_MandX_Z095_low, ICER_CI_SWIMP_MandX_Z095_up=ICER_All$ICER_CI_SWIMP_MandX_Z095_up,
              ICER_CI_SWIMP_MandX_Z095_revert=ICER_All$ICER_CI_SWIMP_MandX_Z095_revert,
              
              ICER_CI_SWIMP_MandX_Z195_low=ICER_All$ICER_CI_SWIMP_MandX_Z195_low, ICER_CI_SWIMP_MandX_Z195_up=ICER_All$ICER_CI_SWIMP_MandX_Z195_up,
              ICER_CI_SWIMP_MandX_Z195_revert=ICER_All$ICER_CI_SWIMP_MandX_Z195_revert,
              # coverage probabilities
              # 95%
              coverSW_MandQ_Z095_Fieller=coverSW_MandQ_Z095_Fieller, coverSW_MandQ_Z195_Fieller=coverSW_MandQ_Z195_Fieller, 
              coverIMP_MandQ_Z095_Fieller=coverIMP_MandQ_Z095_Fieller, coverIMP_MandQ_Z195_Fieller=coverIMP_MandQ_Z195_Fieller,
              coverSW_MandX_Z095_Fieller=coverSW_MandX_Z095_Fieller, coverSW_MandX_Z195_Fieller=coverSW_MandX_Z195_Fieller, 
              coverSWIMP_MandX_Z095_Fieller=coverSWIMP_MandX_Z095_Fieller, coverSWIMP_MandX_Z195_Fieller=coverSWIMP_MandX_Z195_Fieller,
              
              coverSW_MandQ_Z095_Delta=coverSW_MandQ_Z095_Delta, coverSW_MandQ_Z195_Delta=coverSW_MandQ_Z195_Delta, 
              coverIMP_MandQ_Z095_Delta=coverIMP_MandQ_Z095_Delta, coverIMP_MandQ_Z195_Delta=coverIMP_MandQ_Z195_Delta,
              coverSW_MandX_Z095_Delta=coverSW_MandX_Z095_Delta, coverSW_MandX_Z195_Delta=coverSW_MandX_Z195_Delta, 
              coverSWIMP_MandX_Z095_Delta=coverSWIMP_MandX_Z095_Delta, coverSWIMP_MandX_Z195_Delta=coverSWIMP_MandX_Z195_Delta,
              
              # 90%
              coverSW_MandQ_Z090_Fieller=coverSW_MandQ_Z090_Fieller, coverSW_MandQ_Z190_Fieller=coverSW_MandQ_Z190_Fieller, 
              coverIMP_MandQ_Z090_Fieller=coverIMP_MandQ_Z090_Fieller, coverIMP_MandQ_Z190_Fieller=coverIMP_MandQ_Z190_Fieller,
              coverSW_MandX_Z090_Fieller=coverSW_MandX_Z090_Fieller, coverSW_MandX_Z190_Fieller=coverSW_MandX_Z190_Fieller, 
              coverSWIMP_MandX_Z090_Fieller=coverSWIMP_MandX_Z090_Fieller, coverSWIMP_MandX_Z190_Fieller=coverSWIMP_MandX_Z190_Fieller,
              
              coverSW_MandQ_Z090_Delta=coverSW_MandQ_Z090_Delta, coverSW_MandQ_Z190_Delta=coverSW_MandQ_Z190_Delta, 
              coverIMP_MandQ_Z090_Delta=coverIMP_MandQ_Z090_Delta, coverIMP_MandQ_Z190_Delta=coverIMP_MandQ_Z190_Delta,
              coverSW_MandX_Z090_Delta=coverSW_MandX_Z090_Delta, coverSW_MandX_Z190_Delta=coverSW_MandX_Z190_Delta, 
              coverSWIMP_MandX_Z090_Delta=coverSWIMP_MandX_Z090_Delta, coverSWIMP_MandX_Z190_Delta=coverSWIMP_MandX_Z190_Delta,
              
              # 80%
              coverSW_MandQ_Z080_Fieller=coverSW_MandQ_Z080_Fieller, coverSW_MandQ_Z180_Fieller=coverSW_MandQ_Z180_Fieller, 
              coverIMP_MandQ_Z080_Fieller=coverIMP_MandQ_Z080_Fieller, coverIMP_MandQ_Z180_Fieller=coverIMP_MandQ_Z180_Fieller,
              coverSW_MandX_Z080_Fieller=coverSW_MandX_Z080_Fieller, coverSW_MandX_Z180_Fieller=coverSW_MandX_Z180_Fieller, 
              coverSWIMP_MandX_Z080_Fieller=coverSWIMP_MandX_Z080_Fieller, coverSWIMP_MandX_Z180_Fieller=coverSWIMP_MandX_Z180_Fieller,
              
              coverSW_MandQ_Z080_Delta=coverSW_MandQ_Z080_Delta, coverSW_MandQ_Z180_Delta=coverSW_MandQ_Z180_Delta, 
              coverIMP_MandQ_Z080_Delta=coverIMP_MandQ_Z080_Delta, coverIMP_MandQ_Z180_Delta=coverIMP_MandQ_Z180_Delta,
              coverSW_MandX_Z080_Delta=coverSW_MandX_Z080_Delta, coverSW_MandX_Z180_Delta=coverSW_MandX_Z180_Delta, 
              coverSWIMP_MandX_Z080_Delta=coverSWIMP_MandX_Z080_Delta, coverSWIMP_MandX_Z180_Delta=coverSWIMP_MandX_Z180_Delta
              
              
              
              
         ))
}


nsim=2000 	# number of simulation replications
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


cl <- makeCluster(30)

# Will need to export all of these to the cluster
clusterExport(cl = cl, varlist = c("GetEstimate", "sim.CostEff.Obs.V1","ICER","CalCoverage_Fieller","CalCoverage_Delta"))
clusterExport(cl = cl, varlist = c("KM","CalCovComp_SW","CalCovComp_IMP","CalAccum"))
clusterExport(cl = cl, varlist = c("CalCovComp_ICER_SW","CalCovComp_ICER_IMP","CalCovComp_ICER_SWIMP","CalCI_ICER_Fieller","CalCI_ICER_Delta"))
clusterExport(cl = cl, varlist = c("Coef.Surv.MainEff", "Coef.Surv.InterEff","Par.Censor","Par.Censor.0"))
clusterExport(cl = cl, varlist = c("Cost.diag","Cost.fix","Cost.rand","Cost.term"))
clusterExport(cl = cl, varlist = c("Coef.Trt", "max.time", "k", "n","nsim"))
clusterExport(cl = cl, varlist = c("true.ICER.M.Q.Z0", "true.ICER.M.Q.Z1","true.ICER.M.X.Z0","true.ICER.M.X.Z1"))


clusterSetRNGStream(cl = cl, iseed = 12345)
paste0("Started on: ", Sys.time())
start = proc.time()
GetEstimateResult = parLapply(cl = cl, X = 1:nsim, fun = GetEstimate)
stop = proc.time()
paste0("Ended on: ", Sys.time())
runtime = stop - start
runtime

# Save a copy of the result since it took ~4 hours 
saveRDS(object = GetEstimateResult, 
   file = paste("/home/dnliu21/code/CEA/simulation/DC/Results/",
    paste0("Results_ICER_Est_",Censorstatus,"_N", n, "_nsim", nsim, "_logn_", rundate),".RDS",sep=""))



stopCluster(cl)


EstimateResult = data.frame(do.call(what = rbind, args = GetEstimateResult))

estSW_ICER_MQ_Z0=matrix(unlist(EstimateResult$ICER_SW_MandQ_Z0),nsim,1,byrow=TRUE)
estSW_ICER_MQ_Z1=matrix(unlist(EstimateResult$ICER_SW_MandQ_Z1),nsim,1,byrow=TRUE)
estIMP_ICER_MQ_Z0=matrix(unlist(EstimateResult$ICER_IMP_MandQ_Z0),nsim,1,byrow=TRUE)
estIMP_ICER_MQ_Z1=matrix(unlist(EstimateResult$ICER_IMP_MandQ_Z1),nsim,1,byrow=TRUE)

estSW_ICER_MX_Z0=matrix(unlist(EstimateResult$ICER_SW_MandX_Z0),nsim,1,byrow=TRUE)
estSW_ICER_MX_Z1=matrix(unlist(EstimateResult$ICER_SW_MandX_Z1),nsim,1,byrow=TRUE)
estSWIMP_ICER_MX_Z0=matrix(unlist(EstimateResult$ICER_SWIMP_MandX_Z0),nsim,1,byrow=TRUE)
estSWIMP_ICER_MX_Z1=matrix(unlist(EstimateResult$ICER_SWIMP_MandX_Z1),nsim,1,byrow=TRUE)


# Estimated covariance
CovSW_ICER_MQ_Z0=matrix(unlist(EstimateResult$ICER_SW_MandQ_Cov_Z0),nsim,1,byrow=TRUE)
CovSW_ICER_MQ_Z1=matrix(unlist(EstimateResult$ICER_SW_MandQ_Cov_Z1),nsim,1,byrow=TRUE)
CovIMP_ICER_MQ_Z0=matrix(unlist(EstimateResult$ICER_IMP_MandQ_Cov_Z0),nsim,1,byrow=TRUE)
CovIMP_ICER_MQ_Z1=matrix(unlist(EstimateResult$ICER_IMP_MandQ_Cov_Z1),nsim,1,byrow=TRUE)

CovSW_ICER_MX_Z0=matrix(unlist(EstimateResult$ICER_SW_MandX_Cov_Z0),nsim,1,byrow=TRUE)
CovSW_ICER_MX_Z1=matrix(unlist(EstimateResult$ICER_SW_MandX_Cov_Z1),nsim,1,byrow=TRUE)
CovSWIMP_ICER_MX_Z0=matrix(unlist(EstimateResult$ICER_SWIMP_MandX_Cov_Z0),nsim,1,byrow=TRUE)
CovSWIMP_ICER_MX_Z1=matrix(unlist(EstimateResult$ICER_SWIMP_MandX_Cov_Z1),nsim,1,byrow=TRUE)

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


# 95%
coverSW_MandQ_Z095_Fieller=matrix(unlist(EstimateResult$coverSW_MandQ_Z095_Fieller),nsim,1,byrow=TRUE)
coverSW_MandQ_Z195_Fieller=matrix(unlist(EstimateResult$coverSW_MandQ_Z195_Fieller),nsim,1,byrow=TRUE)
coverIMP_MandQ_Z095_Fieller=matrix(unlist(EstimateResult$coverIMP_MandQ_Z095_Fieller),nsim,1,byrow=TRUE)
coverIMP_MandQ_Z195_Fieller=matrix(unlist(EstimateResult$coverIMP_MandQ_Z195_Fieller),nsim,1,byrow=TRUE)
coverSW_MandX_Z095_Fieller=matrix(unlist(EstimateResult$coverSW_MandX_Z095_Fieller),nsim,1,byrow=TRUE)
coverSW_MandX_Z195_Fieller=matrix(unlist(EstimateResult$coverSW_MandX_Z195_Fieller),nsim,1,byrow=TRUE)
coverSWIMP_MandX_Z095_Fieller=matrix(unlist(EstimateResult$coverSWIMP_MandX_Z095_Fieller),nsim,1,byrow=TRUE)
coverSWIMP_MandX_Z195_Fieller=matrix(unlist(EstimateResult$coverSWIMP_MandX_Z195_Fieller),nsim,1,byrow=TRUE)

coverSW_MandQ_Z095_Delta=matrix(unlist(EstimateResult$coverSW_MandQ_Z095_Delta),nsim,1,byrow=TRUE)
coverSW_MandQ_Z195_Delta=matrix(unlist(EstimateResult$coverSW_MandQ_Z195_Delta),nsim,1,byrow=TRUE)
coverIMP_MandQ_Z095_Delta=matrix(unlist(EstimateResult$coverIMP_MandQ_Z095_Delta),nsim,1,byrow=TRUE)
coverIMP_MandQ_Z195_Delta=matrix(unlist(EstimateResult$coverIMP_MandQ_Z195_Delta),nsim,1,byrow=TRUE)
coverSW_MandX_Z095_Delta=matrix(unlist(EstimateResult$coverSW_MandX_Z095_Delta),nsim,1,byrow=TRUE)
coverSW_MandX_Z195_Delta=matrix(unlist(EstimateResult$coverSW_MandX_Z195_Delta),nsim,1,byrow=TRUE)
coverSWIMP_MandX_Z095_Delta=matrix(unlist(EstimateResult$coverSWIMP_MandX_Z095_Delta),nsim,1,byrow=TRUE)
coverSWIMP_MandX_Z195_Delta=matrix(unlist(EstimateResult$coverSWIMP_MandX_Z195_Delta),nsim,1,byrow=TRUE)


# 90%
coverSW_MandQ_Z090_Fieller=matrix(unlist(EstimateResult$coverSW_MandQ_Z090_Fieller),nsim,1,byrow=TRUE)
coverSW_MandQ_Z190_Fieller=matrix(unlist(EstimateResult$coverSW_MandQ_Z190_Fieller),nsim,1,byrow=TRUE)
coverIMP_MandQ_Z090_Fieller=matrix(unlist(EstimateResult$coverIMP_MandQ_Z090_Fieller),nsim,1,byrow=TRUE)
coverIMP_MandQ_Z190_Fieller=matrix(unlist(EstimateResult$coverIMP_MandQ_Z190_Fieller),nsim,1,byrow=TRUE)
coverSW_MandX_Z090_Fieller=matrix(unlist(EstimateResult$coverSW_MandX_Z090_Fieller),nsim,1,byrow=TRUE)
coverSW_MandX_Z190_Fieller=matrix(unlist(EstimateResult$coverSW_MandX_Z190_Fieller),nsim,1,byrow=TRUE)
coverSWIMP_MandX_Z090_Fieller=matrix(unlist(EstimateResult$coverSWIMP_MandX_Z090_Fieller),nsim,1,byrow=TRUE)
coverSWIMP_MandX_Z190_Fieller=matrix(unlist(EstimateResult$coverSWIMP_MandX_Z190_Fieller),nsim,1,byrow=TRUE)

coverSW_MandQ_Z090_Delta=matrix(unlist(EstimateResult$coverSW_MandQ_Z090_Delta),nsim,1,byrow=TRUE)
coverSW_MandQ_Z190_Delta=matrix(unlist(EstimateResult$coverSW_MandQ_Z190_Delta),nsim,1,byrow=TRUE)
coverIMP_MandQ_Z090_Delta=matrix(unlist(EstimateResult$coverIMP_MandQ_Z090_Delta),nsim,1,byrow=TRUE)
coverIMP_MandQ_Z190_Delta=matrix(unlist(EstimateResult$coverIMP_MandQ_Z190_Delta),nsim,1,byrow=TRUE)
coverSW_MandX_Z090_Delta=matrix(unlist(EstimateResult$coverSW_MandX_Z090_Delta),nsim,1,byrow=TRUE)
coverSW_MandX_Z190_Delta=matrix(unlist(EstimateResult$coverSW_MandX_Z190_Delta),nsim,1,byrow=TRUE)
coverSWIMP_MandX_Z090_Delta=matrix(unlist(EstimateResult$coverSWIMP_MandX_Z090_Delta),nsim,1,byrow=TRUE)
coverSWIMP_MandX_Z190_Delta=matrix(unlist(EstimateResult$coverSWIMP_MandX_Z190_Delta),nsim,1,byrow=TRUE)


# 80%
coverSW_MandQ_Z080_Fieller=matrix(unlist(EstimateResult$coverSW_MandQ_Z080_Fieller),nsim,1,byrow=TRUE)
coverSW_MandQ_Z180_Fieller=matrix(unlist(EstimateResult$coverSW_MandQ_Z180_Fieller),nsim,1,byrow=TRUE)
coverIMP_MandQ_Z080_Fieller=matrix(unlist(EstimateResult$coverIMP_MandQ_Z080_Fieller),nsim,1,byrow=TRUE)
coverIMP_MandQ_Z180_Fieller=matrix(unlist(EstimateResult$coverIMP_MandQ_Z180_Fieller),nsim,1,byrow=TRUE)
coverSW_MandX_Z080_Fieller=matrix(unlist(EstimateResult$coverSW_MandX_Z080_Fieller),nsim,1,byrow=TRUE)
coverSW_MandX_Z180_Fieller=matrix(unlist(EstimateResult$coverSW_MandX_Z180_Fieller),nsim,1,byrow=TRUE)
coverSWIMP_MandX_Z080_Fieller=matrix(unlist(EstimateResult$coverSWIMP_MandX_Z080_Fieller),nsim,1,byrow=TRUE)
coverSWIMP_MandX_Z180_Fieller=matrix(unlist(EstimateResult$coverSWIMP_MandX_Z180_Fieller),nsim,1,byrow=TRUE)

coverSW_MandQ_Z080_Delta=matrix(unlist(EstimateResult$coverSW_MandQ_Z080_Delta),nsim,1,byrow=TRUE)
coverSW_MandQ_Z180_Delta=matrix(unlist(EstimateResult$coverSW_MandQ_Z180_Delta),nsim,1,byrow=TRUE)
coverIMP_MandQ_Z080_Delta=matrix(unlist(EstimateResult$coverIMP_MandQ_Z080_Delta),nsim,1,byrow=TRUE)
coverIMP_MandQ_Z180_Delta=matrix(unlist(EstimateResult$coverIMP_MandQ_Z180_Delta),nsim,1,byrow=TRUE)
coverSW_MandX_Z080_Delta=matrix(unlist(EstimateResult$coverSW_MandX_Z080_Delta),nsim,1,byrow=TRUE)
coverSW_MandX_Z180_Delta=matrix(unlist(EstimateResult$coverSW_MandX_Z180_Delta),nsim,1,byrow=TRUE)
coverSWIMP_MandX_Z080_Delta=matrix(unlist(EstimateResult$coverSWIMP_MandX_Z080_Delta),nsim,1,byrow=TRUE)
coverSWIMP_MandX_Z180_Delta=matrix(unlist(EstimateResult$coverSWIMP_MandX_Z180_Delta),nsim,1,byrow=TRUE)



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



meanSW_ICER_MQ_Z0=apply(estSW_ICER_MQ_Z0,2,mean,na.rm=TRUE)
biasSW_ICER_MQ_Z0=apply(estSW_ICER_MQ_Z0,2,mean,na.rm=TRUE)-true.ICER.M.Q.Z0 
ECovSW_ICER_MQ_Z0=apply(CovSW_ICER_MQ_Z0,2,mean,na.rm=TRUE)
SCovSW_ICER_MQ_Z0=cov(estDiff_SW_M_Z0, estDiff_SW_Q_Z0)
# 95%
CovPSW_ICER_MQ_Z095_Fieller=apply(coverSW_MandQ_Z095_Fieller,2,mean,na.rm=TRUE)
CovPSW_ICER_MQ_Z095_Delta=apply(coverSW_MandQ_Z095_Delta,2,mean,na.rm=TRUE)
# 90%
CovPSW_ICER_MQ_Z090_Fieller=apply(coverSW_MandQ_Z090_Fieller,2,mean,na.rm=TRUE)
CovPSW_ICER_MQ_Z090_Delta=apply(coverSW_MandQ_Z090_Delta,2,mean,na.rm=TRUE)
# 80%
CovPSW_ICER_MQ_Z080_Fieller=apply(coverSW_MandQ_Z080_Fieller,2,mean,na.rm=TRUE)
CovPSW_ICER_MQ_Z080_Delta=apply(coverSW_MandQ_Z080_Delta,2,mean,na.rm=TRUE)
meanDiffSW_M_x_Z0=apply(estDiff_SW_M_Z0,2,mean,na.rm=TRUE)
biasSW_M_x_Z0=apply(estDiff_SW_M_Z0,2,mean,na.rm=TRUE)-true.beta.M[3]
meanDiffSW_Q_y_Z0=apply(estDiff_SW_Q_Z0,2,mean,na.rm=TRUE)
biasSW_Q_y_Z0=apply(estDiff_SW_Q_Z0,2,mean,na.rm=TRUE)-true.beta.Q[3]
ECovSW_M_Sxx_Z0=apply(Sxx_SW_M_Z0,2,mean,na.rm=TRUE)
SCovSW_M_Sxx_Z0=var(estDiff_SW_M_Z0)
ECovSW_Q_Syy_Z0=apply(Syy_SW_Q_Z0,2,mean,na.rm=TRUE)
SCovSW_Q_Syy_Z0=var(estDiff_SW_Q_Z0)

meanSW_ICER_MQ_Z1=apply(estSW_ICER_MQ_Z1,2,mean,na.rm=TRUE)
biasSW_ICER_MQ_Z1=apply(estSW_ICER_MQ_Z1,2,mean,na.rm=TRUE)-true.ICER.M.Q.Z1
ECovSW_ICER_MQ_Z1=apply(CovSW_ICER_MQ_Z1,2,mean,na.rm=TRUE)
SCovSW_ICER_MQ_Z1=cov(estDiff_SW_M_Z1, estDiff_SW_Q_Z1)
# 95%
CovPSW_ICER_MQ_Z195_Fieller=apply(coverSW_MandQ_Z195_Fieller,2,mean,na.rm=TRUE)
CovPSW_ICER_MQ_Z195_Delta=apply(coverSW_MandQ_Z195_Delta,2,mean,na.rm=TRUE)
# 90%
CovPSW_ICER_MQ_Z190_Fieller=apply(coverSW_MandQ_Z190_Fieller,2,mean,na.rm=TRUE)
CovPSW_ICER_MQ_Z190_Delta=apply(coverSW_MandQ_Z190_Delta,2,mean,na.rm=TRUE)
# 80%
CovPSW_ICER_MQ_Z180_Fieller=apply(coverSW_MandQ_Z180_Fieller,2,mean,na.rm=TRUE)
CovPSW_ICER_MQ_Z180_Delta=apply(coverSW_MandQ_Z180_Delta,2,mean,na.rm=TRUE)
meanDiffSW_M_x_Z1=apply(estDiff_SW_M_Z1,2,mean,na.rm=TRUE)
biasSW_M_x_Z1=apply(estDiff_SW_M_Z1,2,mean,na.rm=TRUE)-(true.beta.M[3]+true.beta.M[4])
meanDiffSW_Q_y_Z1=apply(estDiff_SW_Q_Z1,2,mean,na.rm=TRUE)
biasSW_Q_y_Z1=apply(estDiff_SW_Q_Z1,2,mean,na.rm=TRUE)-(true.beta.Q[3]+true.beta.Q[4])
ECovSW_M_Sxx_Z1=apply(Sxx_SW_M_Z1,2,mean,na.rm=TRUE)
SCovSW_M_Sxx_Z1=var(estDiff_SW_M_Z1)
ECovSW_Q_Syy_Z1=apply(Syy_SW_Q_Z1,2,mean,na.rm=TRUE)
SCovSW_Q_Syy_Z1=var(estDiff_SW_Q_Z1)

meanIMP_ICER_MQ_Z0=apply(estIMP_ICER_MQ_Z0,2,mean,na.rm=TRUE)
biasIMP_ICER_MQ_Z0=apply(estIMP_ICER_MQ_Z0,2,mean,na.rm=TRUE)-true.ICER.M.Q.Z0 
ECovIMP_ICER_MQ_Z0=apply(CovIMP_ICER_MQ_Z0,2,mean,na.rm=TRUE)
SCovIMP_ICER_MQ_Z0=cov(estDiff_IMP_M_Z0, estDiff_IMP_Q_Z0)
# 95%
CovPIMP_ICER_MQ_Z095_Fieller=apply(coverIMP_MandQ_Z095_Fieller,2,mean,na.rm=TRUE)
CovPIMP_ICER_MQ_Z095_Delta=apply(coverIMP_MandQ_Z095_Delta,2,mean,na.rm=TRUE)
# 90%
CovPIMP_ICER_MQ_Z090_Fieller=apply(coverIMP_MandQ_Z090_Fieller,2,mean,na.rm=TRUE)
CovPIMP_ICER_MQ_Z090_Delta=apply(coverIMP_MandQ_Z090_Delta,2,mean,na.rm=TRUE)
# 80%
CovPIMP_ICER_MQ_Z080_Fieller=apply(coverIMP_MandQ_Z080_Fieller,2,mean,na.rm=TRUE)
CovPIMP_ICER_MQ_Z080_Delta=apply(coverIMP_MandQ_Z080_Delta,2,mean,na.rm=TRUE)
meanDiffIMP_M_x_Z0=apply(estDiff_IMP_M_Z0,2,mean,na.rm=TRUE)
biasIMP_M_x_Z0=apply(estDiff_IMP_M_Z0,2,mean,na.rm=TRUE)-true.beta.M[3]
meanDiffIMP_Q_y_Z0=apply(estDiff_IMP_Q_Z0,2,mean,na.rm=TRUE)
biasIMP_Q_y_Z0=apply(estDiff_IMP_Q_Z0,2,mean,na.rm=TRUE)-true.beta.Q[3]
ECovIMP_M_Sxx_Z0=apply(Sxx_IMP_M_Z0,2,mean,na.rm=TRUE)
SCovIMP_M_Sxx_Z0=var(estDiff_IMP_M_Z0)
ECovIMP_Q_Syy_Z0=apply(Syy_IMP_Q_Z0,2,mean,na.rm=TRUE)
SCovIMP_Q_Syy_Z0=var(estDiff_IMP_Q_Z0)

meanIMP_ICER_MQ_Z1=apply(estIMP_ICER_MQ_Z1,2,mean,na.rm=TRUE)
biasIMP_ICER_MQ_Z1=apply(estIMP_ICER_MQ_Z1,2,mean,na.rm=TRUE)-true.ICER.M.Q.Z1
ECovIMP_ICER_MQ_Z1=apply(CovIMP_ICER_MQ_Z1,2,mean,na.rm=TRUE)
SCovIMP_ICER_MQ_Z1=cov(estDiff_IMP_M_Z1, estDiff_IMP_Q_Z1)
# 95%
CovPIMP_ICER_MQ_Z195_Fieller=apply(coverIMP_MandQ_Z195_Fieller,2,mean,na.rm=TRUE)
CovPIMP_ICER_MQ_Z195_Delta=apply(coverIMP_MandQ_Z195_Delta,2,mean,na.rm=TRUE)
# 90%
CovPIMP_ICER_MQ_Z190_Fieller=apply(coverIMP_MandQ_Z190_Fieller,2,mean,na.rm=TRUE)
CovPIMP_ICER_MQ_Z190_Delta=apply(coverIMP_MandQ_Z190_Delta,2,mean,na.rm=TRUE)
# 80%
CovPIMP_ICER_MQ_Z180_Fieller=apply(coverIMP_MandQ_Z180_Fieller,2,mean,na.rm=TRUE)
CovPIMP_ICER_MQ_Z180_Delta=apply(coverIMP_MandQ_Z180_Delta,2,mean,na.rm=TRUE)
meanDiffIMP_M_x_Z1=apply(estDiff_IMP_M_Z1,2,mean,na.rm=TRUE)
biasIMP_M_x_Z1=apply(estDiff_IMP_M_Z1,2,mean,na.rm=TRUE)-(true.beta.M[3]+true.beta.M[4])
meanDiffIMP_Q_y_Z1=apply(estDiff_IMP_Q_Z1,2,mean,na.rm=TRUE)
biasIMP_Q_y_Z1=apply(estDiff_IMP_Q_Z1,2,mean,na.rm=TRUE)-(true.beta.Q[3]+true.beta.Q[4])
ECovIMP_M_Sxx_Z1=apply(Sxx_IMP_M_Z1,2,mean,na.rm=TRUE)
SCovIMP_M_Sxx_Z1=var(estDiff_IMP_M_Z1)
ECovIMP_Q_Syy_Z1=apply(Syy_IMP_Q_Z1,2,mean,na.rm=TRUE)
SCovIMP_Q_Syy_Z1=var(estDiff_IMP_Q_Z1)

meanSW_ICER_MX_Z0=apply(estSW_ICER_MX_Z0,2,mean,na.rm=TRUE)
biasSW_ICER_MX_Z0=apply(estSW_ICER_MX_Z0,2,mean,na.rm=TRUE)-true.ICER.M.X.Z0 
ECovSW_ICER_MX_Z0=apply(CovSW_ICER_MX_Z0,2,mean,na.rm=TRUE)
SCovSW_ICER_MX_Z0=cov(estDiff_SW_M_Z0, estDiff_SW_X_Z0)
# 95%
CovPSW_ICER_MX_Z095_Fieller=apply(coverSW_MandX_Z095_Fieller,2,mean,na.rm=TRUE)
CovPSW_ICER_MX_Z095_Delta=apply(coverSW_MandX_Z095_Delta,2,mean,na.rm=TRUE)
# 90%
CovPSW_ICER_MX_Z090_Fieller=apply(coverSW_MandX_Z090_Fieller,2,mean,na.rm=TRUE)
CovPSW_ICER_MX_Z090_Delta=apply(coverSW_MandX_Z090_Delta,2,mean,na.rm=TRUE)
# 80%
CovPSW_ICER_MX_Z080_Fieller=apply(coverSW_MandX_Z080_Fieller,2,mean,na.rm=TRUE)
CovPSW_ICER_MX_Z080_Delta=apply(coverSW_MandX_Z080_Delta,2,mean,na.rm=TRUE)
meanDiffSW_M_x_Z0=apply(estDiff_SW_M_Z0,2,mean,na.rm=TRUE)
biasSW_M_x_Z0=apply(estDiff_SW_M_Z0,2,mean,na.rm=TRUE)-true.beta.M[3]
meanDiffSW_X_y_Z0=apply(estDiff_SW_X_Z0,2,mean,na.rm=TRUE)
biasSW_X_y_Z0=apply(estDiff_SW_X_Z0,2,mean,na.rm=TRUE)-true.beta.TL[3]
ECovSW_M_Sxx_Z0=apply(Sxx_SW_M_Z0,2,mean,na.rm=TRUE)
SCovSW_M_Sxx_Z0=var(estDiff_SW_M_Z0)
ECovSW_X_Syy_Z0=apply(Syy_SW_X_Z0,2,mean,na.rm=TRUE)
SCovSW_X_Syy_Z0=var(estDiff_SW_X_Z0)

meanSW_ICER_MX_Z1=apply(estSW_ICER_MX_Z1,2,mean,na.rm=TRUE)
biasSW_ICER_MX_Z1=apply(estSW_ICER_MX_Z1,2,mean,na.rm=TRUE)-true.ICER.M.X.Z1 
ECovSW_ICER_MX_Z1=apply(CovSW_ICER_MX_Z1,2,mean,na.rm=TRUE)
SCovSW_ICER_MX_Z1=cov(estDiff_SW_M_Z1, estDiff_SW_X_Z1)
# 95%
CovPSW_ICER_MX_Z195_Fieller=apply(coverSW_MandX_Z195_Fieller,2,mean,na.rm=TRUE)
CovPSW_ICER_MX_Z195_Delta=apply(coverSW_MandX_Z195_Delta,2,mean,na.rm=TRUE)
# 90%
CovPSW_ICER_MX_Z190_Fieller=apply(coverSW_MandX_Z190_Fieller,2,mean,na.rm=TRUE)
CovPSW_ICER_MX_Z190_Delta=apply(coverSW_MandX_Z190_Delta,2,mean,na.rm=TRUE)
# 80%
CovPSW_ICER_MX_Z180_Fieller=apply(coverSW_MandX_Z180_Fieller,2,mean,na.rm=TRUE)
CovPSW_ICER_MX_Z180_Delta=apply(coverSW_MandX_Z180_Delta,2,mean,na.rm=TRUE)
meanDiffSW_M_x_Z1=apply(estDiff_SW_M_Z1,2,mean,na.rm=TRUE)
biasSW_M_x_Z1=apply(estDiff_SW_M_Z1,2,mean,na.rm=TRUE)-(true.beta.M[3]+true.beta.M[4])
meanDiffSW_X_y_Z1=apply(estDiff_SW_X_Z1,2,mean,na.rm=TRUE)
biasSW_X_y_Z1=apply(estDiff_SW_X_Z1,2,mean,na.rm=TRUE)-(true.beta.TL[3]+true.beta.TL[4])
ECovSW_M_Sxx_Z1=apply(Sxx_SW_M_Z1,2,mean,na.rm=TRUE)
SCovSW_M_Sxx_Z1=var(estDiff_SW_M_Z1)
ECovSW_X_Syy_Z1=apply(Syy_SW_X_Z1,2,mean,na.rm=TRUE)
SCovSW_X_Syy_Z1=var(estDiff_SW_X_Z1)

meanSWIMP_ICER_MX_Z0=apply(estSWIMP_ICER_MX_Z0,2,mean,na.rm=TRUE)
biasSWIMP_ICER_MX_Z0=apply(estSWIMP_ICER_MX_Z0,2,mean,na.rm=TRUE)-true.ICER.M.X.Z0 
ECovSWIMP_ICER_MX_Z0=apply(CovSWIMP_ICER_MX_Z0,2,mean,na.rm=TRUE)
SCovSWIMP_ICER_MX_Z0=cov(estDiff_IMP_M_Z0, estDiff_SW_X_Z0)
# 95%
CovPSWIMP_ICER_MX_Z095_Fieller=apply(coverSWIMP_MandX_Z095_Fieller,2,mean,na.rm=TRUE)
CovPSWIMP_ICER_MX_Z095_Delta=apply(coverSWIMP_MandX_Z095_Delta,2,mean,na.rm=TRUE)
# 90%
CovPSWIMP_ICER_MX_Z090_Fieller=apply(coverSWIMP_MandX_Z090_Fieller,2,mean,na.rm=TRUE)
CovPSWIMP_ICER_MX_Z090_Delta=apply(coverSWIMP_MandX_Z090_Delta,2,mean,na.rm=TRUE)
# 80%
CovPSWIMP_ICER_MX_Z080_Fieller=apply(coverSWIMP_MandX_Z080_Fieller,2,mean,na.rm=TRUE)
CovPSWIMP_ICER_MX_Z080_Delta=apply(coverSWIMP_MandX_Z080_Delta,2,mean,na.rm=TRUE)
meanDiffIMP_M_x_Z0=apply(estDiff_IMP_M_Z0,2,mean,na.rm=TRUE)
biasIMP_M_x_Z0=apply(estDiff_IMP_M_Z0,2,mean,na.rm=TRUE)-true.beta.M[3]
meanDiffSW_X_y_Z0=apply(estDiff_SW_X_Z0,2,mean,na.rm=TRUE)
biasSW_X_y_Z0=apply(estDiff_SW_X_Z0,2,mean,na.rm=TRUE)-true.beta.TL[3]
ECovIMP_M_Sxx_Z0=apply(Sxx_IMP_M_Z0,2,mean,na.rm=TRUE)
SCovIMP_M_Sxx_Z0=var(estDiff_IMP_M_Z0)
ECovSW_X_Syy_Z0=apply(Syy_SW_X_Z0,2,mean,na.rm=TRUE)
SCovSW_X_Syy_Z0=var(estDiff_SW_X_Z0)

meanSWIMP_ICER_MX_Z1=apply(estSWIMP_ICER_MX_Z1,2,mean,na.rm=TRUE)
biasSWIMP_ICER_MX_Z1=apply(estSWIMP_ICER_MX_Z1,2,mean,na.rm=TRUE)-true.ICER.M.X.Z1 
ECovSWIMP_ICER_MX_Z1=apply(CovSWIMP_ICER_MX_Z1,2,mean,na.rm=TRUE)
SCovSWIMP_ICER_MX_Z1=cov(estDiff_IMP_M_Z1, estDiff_SW_X_Z1)
# 95%
CovPSWIMP_ICER_MX_Z195_Fieller=apply(coverSWIMP_MandX_Z195_Fieller,2,mean,na.rm=TRUE)
CovPSWIMP_ICER_MX_Z195_Delta=apply(coverSWIMP_MandX_Z195_Delta,2,mean,na.rm=TRUE)
# 90%
CovPSWIMP_ICER_MX_Z190_Fieller=apply(coverSWIMP_MandX_Z190_Fieller,2,mean,na.rm=TRUE)
CovPSWIMP_ICER_MX_Z190_Delta=apply(coverSWIMP_MandX_Z190_Delta,2,mean,na.rm=TRUE)
# 80%
CovPSWIMP_ICER_MX_Z180_Fieller=apply(coverSWIMP_MandX_Z180_Fieller,2,mean,na.rm=TRUE)
CovPSWIMP_ICER_MX_Z180_Delta=apply(coverSWIMP_MandX_Z180_Delta,2,mean,na.rm=TRUE)
meanDiffIMP_M_x_Z1=apply(estDiff_IMP_M_Z1,2,mean,na.rm=TRUE)
biasIMP_M_x_Z1=apply(estDiff_IMP_M_Z1,2,mean,na.rm=TRUE)-(true.beta.M[3]+true.beta.M[4])
meanDiffSW_X_y_Z1=apply(estDiff_SW_X_Z1,2,mean,na.rm=TRUE)
biasSW_X_y_Z1=apply(estDiff_SW_X_Z1,2,mean,na.rm=TRUE)-(true.beta.TL[3]+true.beta.TL[4])
ECovIMP_M_Sxx_Z1=apply(Sxx_IMP_M_Z1,2,mean,na.rm=TRUE)
SCovIMP_M_Sxx_Z1=var(estDiff_IMP_M_Z1)
ECovSW_X_Syy_Z1=apply(Syy_SW_X_Z1,2,mean,na.rm=TRUE)
SCovSW_X_Syy_Z1=var(estDiff_SW_X_Z1)




results = data.frame(
  rep(c("ICER"),1,each=4),
  rep(c("SW","IMP"),2,each=1),
  rep(c("Z=0","Z=1"),1,each=2),
  round(t(cbind(meanSW_ICER_MQ_Z0, meanIMP_ICER_MQ_Z0, meanSW_ICER_MQ_Z1, meanIMP_ICER_MQ_Z1)),4),
  round(t(cbind(biasSW_ICER_MQ_Z0, biasIMP_ICER_MQ_Z0, biasSW_ICER_MQ_Z1, biasIMP_ICER_MQ_Z1)),4),
  round(t(cbind(ECovSW_ICER_MQ_Z0, ECovIMP_ICER_MQ_Z0, ECovSW_ICER_MQ_Z1, ECovIMP_ICER_MQ_Z1)),4),
  round(t(cbind(SCovSW_ICER_MQ_Z0, SCovIMP_ICER_MQ_Z0, SCovSW_ICER_MQ_Z1, SCovIMP_ICER_MQ_Z1)),4),
  # 95%
  round(t(cbind(CovPSW_ICER_MQ_Z095_Fieller, CovPIMP_ICER_MQ_Z095_Fieller, CovPSW_ICER_MQ_Z195_Fieller, CovPIMP_ICER_MQ_Z195_Fieller)),4),
  round(t(cbind(CovPSW_ICER_MQ_Z095_Delta, CovPIMP_ICER_MQ_Z095_Delta, CovPSW_ICER_MQ_Z195_Delta, CovPIMP_ICER_MQ_Z195_Delta)),4),
  # 90%
  round(t(cbind(CovPSW_ICER_MQ_Z090_Fieller, CovPIMP_ICER_MQ_Z090_Fieller, CovPSW_ICER_MQ_Z190_Fieller, CovPIMP_ICER_MQ_Z190_Fieller)),4),
  round(t(cbind(CovPSW_ICER_MQ_Z090_Delta, CovPIMP_ICER_MQ_Z090_Delta, CovPSW_ICER_MQ_Z190_Delta, CovPIMP_ICER_MQ_Z190_Delta)),4),
  # 80%
  round(t(cbind(CovPSW_ICER_MQ_Z080_Fieller, CovPIMP_ICER_MQ_Z080_Fieller, CovPSW_ICER_MQ_Z180_Fieller, CovPIMP_ICER_MQ_Z180_Fieller)),4),
  round(t(cbind(CovPSW_ICER_MQ_Z080_Delta, CovPIMP_ICER_MQ_Z080_Delta, CovPSW_ICER_MQ_Z180_Delta, CovPIMP_ICER_MQ_Z180_Delta)),4),
  round(t(cbind(meanDiffSW_M_x_Z0, meanDiffIMP_M_x_Z0, meanDiffSW_M_x_Z1, meanDiffIMP_M_x_Z1)),4),
  round(t(cbind(biasSW_M_x_Z0, biasIMP_M_x_Z0, biasSW_M_x_Z1, biasIMP_M_x_Z1)),4),
  round(t(cbind(meanDiffSW_Q_y_Z0, meanDiffIMP_Q_y_Z0, meanDiffSW_Q_y_Z1, meanDiffIMP_Q_y_Z1)),4),
  round(t(cbind(biasSW_Q_y_Z0, biasIMP_Q_y_Z0, biasSW_Q_y_Z1, biasIMP_Q_y_Z1)),4),
  round(t(cbind(ECovSW_M_Sxx_Z0, ECovIMP_M_Sxx_Z0, ECovSW_M_Sxx_Z1, ECovIMP_M_Sxx_Z1)),4),
  round(t(cbind(SCovSW_M_Sxx_Z0, SCovIMP_M_Sxx_Z0, SCovSW_M_Sxx_Z1, SCovIMP_M_Sxx_Z1)),4),
  round(t(cbind(ECovSW_Q_Syy_Z0, ECovIMP_Q_Syy_Z0, ECovSW_Q_Syy_Z1, ECovIMP_Q_Syy_Z1)),4),
  round(t(cbind(SCovSW_Q_Syy_Z0, SCovIMP_Q_Syy_Z0, SCovSW_Q_Syy_Z1, SCovIMP_Q_Syy_Z1)),4),
  
  round(t(cbind(meanSW_ICER_MX_Z0, meanSWIMP_ICER_MX_Z0, meanSW_ICER_MX_Z1, meanSWIMP_ICER_MX_Z1)),4),
  round(t(cbind(biasSW_ICER_MX_Z0, biasSWIMP_ICER_MX_Z0, biasSW_ICER_MX_Z1, biasSWIMP_ICER_MX_Z1)),4),
  round(t(cbind(ECovSW_ICER_MX_Z0, ECovSWIMP_ICER_MX_Z0, ECovSW_ICER_MX_Z1, ECovSWIMP_ICER_MX_Z1)),4),
  round(t(cbind(SCovSW_ICER_MX_Z0, SCovSWIMP_ICER_MX_Z0, SCovSW_ICER_MX_Z1, SCovSWIMP_ICER_MX_Z1)),4),
  # 95%
  round(t(cbind(CovPSW_ICER_MX_Z095_Fieller, CovPSWIMP_ICER_MX_Z095_Fieller, CovPSW_ICER_MX_Z195_Fieller, CovPSWIMP_ICER_MX_Z195_Fieller)),4),
  round(t(cbind(CovPSW_ICER_MX_Z095_Delta, CovPSWIMP_ICER_MX_Z095_Delta, CovPSW_ICER_MX_Z195_Delta, CovPSWIMP_ICER_MX_Z195_Delta)),4),
  # 90%
  round(t(cbind(CovPSW_ICER_MX_Z090_Fieller, CovPSWIMP_ICER_MX_Z090_Fieller, CovPSW_ICER_MX_Z190_Fieller, CovPSWIMP_ICER_MX_Z190_Fieller)),4),
  round(t(cbind(CovPSW_ICER_MX_Z090_Delta, CovPSWIMP_ICER_MX_Z090_Delta, CovPSW_ICER_MX_Z190_Delta, CovPSWIMP_ICER_MX_Z190_Delta)),4),
  # 80%
  round(t(cbind(CovPSW_ICER_MX_Z080_Fieller, CovPSWIMP_ICER_MX_Z080_Fieller, CovPSW_ICER_MX_Z180_Fieller, CovPSWIMP_ICER_MX_Z180_Fieller)),4),
  round(t(cbind(CovPSW_ICER_MX_Z080_Delta, CovPSWIMP_ICER_MX_Z080_Delta, CovPSW_ICER_MX_Z180_Delta, CovPSWIMP_ICER_MX_Z180_Delta)),4),
  round(t(cbind(meanDiffSW_M_x_Z0, meanDiffIMP_M_x_Z0, meanDiffSW_M_x_Z1, meanDiffIMP_M_x_Z1)),4),
  round(t(cbind(biasSW_M_x_Z0, biasIMP_M_x_Z0, biasSW_M_x_Z1, biasIMP_M_x_Z1)),4),
  round(t(cbind(meanDiffSW_X_y_Z0, NA, meanDiffSW_X_y_Z1, NA)),4),
  round(t(cbind(biasSW_X_y_Z0, NA, biasSW_X_y_Z1, NA)),4),
  round(t(cbind(ECovSW_M_Sxx_Z0, ECovIMP_M_Sxx_Z0, ECovSW_M_Sxx_Z1, ECovIMP_M_Sxx_Z1)),4),
  round(t(cbind(SCovSW_M_Sxx_Z0, SCovIMP_M_Sxx_Z0, SCovSW_M_Sxx_Z1, SCovIMP_M_Sxx_Z1)),4),
  round(t(cbind(ECovSW_X_Syy_Z0, NA, ECovSW_X_Syy_Z1 ,NA)),4),
  round(t(cbind(SCovSW_X_Syy_Z0, NA, SCovSW_X_Syy_Z1 ,NA)),4)
  )


colnames(results)=c("Variable","Method","Sep.K",
                    "EMean ICER M&Q","Bias M&Q","ECov M&Q (S_xy)","SCov M&Q","CP M&Q Fieller 95%","CP M&Q Delta 95%","CP M&Q Fieller 90%","CP M&Q Delta 90%","CP M&Q Fieller 80%","CP M&Q Delta 80%",
                    "Ex_M","Bias M","Ey_Q","Bias Q","E_S_M_xx","S_S_M_xx","E_S_Q_yy","S_S_Q_yy",
                    "EMean ICER M&X","Bias M&X","ECov M&X (S_xy)","SCov M&X","CP M&X Fieller 95%","CP M&X Delta 95%","CP M&X Fieller 90%","CP M&X Delta 90%","CP M&X Fieller 80%","CP M&X Delta 80%",
                    "Ex_M","Bias M","Ey_X","Bias X","E_S_M_xx","S_S_M_xx","E_S_X_yy","S_S_X_yy"
                    )


write.csv(results, paste("/home/dnliu21/code/CEA/simulation/DC/Results/",
            "Results_ICER_Bias_",Censorstatus,"_N", n, "_nsim", nsim, "_logn_", rundate,".csv",sep=""), row.names = F)



