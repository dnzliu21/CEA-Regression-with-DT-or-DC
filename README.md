This repository contains codes necessary to replicate the simulation results of Liu and Chen (2025+): "Regression methods for cost-effectiveness analysis with different censoring times or terminating events for survival time and costs".

**Folder and files**

| Folder | Files | Details |
| --- | --- | --- |
| Different terminating events | Bootstrap | Script for Bootstrap methods in DT case |
|  | DataGen_DT.r | Data generation for simulation in DT case |
|  | ICER_Estimate_DT.r | Functions for estimators, variance, covariances, and CI in DT case |
|  | ICER_Parallel_Main_DT.r | Scripts for simulation process in DT case |
|  | INB_Estimate_DT.r | Scripts for INB results for simulation in DT case |
| Different censoring times | Bootstrap | Script for Bootstrap methods in DC case |
|  | DataGen_DC.r | Data generation for simulation in DC case |
|  | ICER_Estimate_DC.r | Functions for estimators, variance, covariances, and CI in DC case |
|  | ICER_Parallel_Main_DC.r | Scripts for obtaining INB results for simulation in DC case |
|  | INB_Estimate_DC.r | Scripts for INB results for simulation in DC case |


Note: DT:Different terminating events, DC: Different censoring times.
