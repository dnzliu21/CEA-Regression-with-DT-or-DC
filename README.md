This repository contains codes necessary to replicate the simulation results of Liu and Chen (2025+): "Regression methods for cost-effectiveness analysis with different censoring times or terminating events for survival time and costs".

**Folder and files**

| Folder | Files | Details |
| --- | --- | --- |
| Different terminating events | Bootstrap | R code for generating results from Bootstrap methods |
|  | DataGen_DT.r | Data generation for simulation in DT case |
|  | ICER_Estimate_DT.r | Functions for obtaining estimators, variance, covariances, and CI in DT case |
|  | ICER_Parallel_Main_DT.r | Scripts for simulation process in DT case |
|  | INB_Estimate_DT.r | Scripts for obtaining INB resultsfor simulation in DT case |
| Different censoring times | Bootstrap |
|  | DataGen_DC.r |
|  | ICER_Estimate_DC.r |
|  | ICER_Parallel_Main_DC.r |
|  | INB_Estimate_DC.r |
