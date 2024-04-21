# Mendelian-randomization-on-diet-and-T2D

Mendelian randomization approach to assess causality in diet and type 2 diabetes mellitus.
Dietary preferences as exposures 
T2D and 21 related cardiometabolic risk factors as outcomes

The MR approach consists of univariable and multivariable including two-step mediation analyses. 
  For univariable MR, TwoSampleMR R package was used.
  For multivariable MR, MVMR along with TwoSampleMR packages was used.
  For two-step MR mediation, RMediation R package was used to estimate indirect (mediated) effect and calculate confidence intervals.
  
To conduct these analyses the following codes/scripts were used.

  1. process_gwas_summary_stats.py (python script for rsID standardization using dbSNP build 155 to GRCh37, written by Mitch - cite)

  2. UVMR_diet_T2D_TwoSampleMR.R

  3. MVMR_diet_T2D_using_both_2samplemr_mvmr_with_proxy.R

  4. UVMR_MR-RAPS_diet_T2D_TwoSampleMR.R

  5. MVMR_MR-RAPS_diet_T2D_using_both_2samplemr_mvmr.R
